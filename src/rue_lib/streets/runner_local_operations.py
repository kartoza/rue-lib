# src/rue_lib/streets/runner_local_operations.py
"""Operations specific to the runner_local workflow."""

from pathlib import Path
from typing import Optional

import geopandas as gpd
from shapely import LineString, unary_union
from shapely.ops import polygonize

from rue_lib.core.geometry_sampling import extend_line

from .runner_utils import subtract_layer


def classify_blocks_by_road_type(
    output_path: Path,
    base_polygons_layer: str,
    roads_buffer_layer: str,
    dead_end_lines_layer: str,
    site_layer: str = "00_site",
    output_layer_name: str = "08_classified_blocks",
    road_type_column: str = "road_type",
) -> str:
    """Classify blocks based on which road type they are adjacent to.

    Classification logic:
    - If polygon is adjacent to arterial road (road_type='road_art') → block_type='on_grid_art'
    - If polygon is adjacent to secondary road (road_type='road_sec') → block_type='on_grid_sec'
    - Otherwise → block_type='off_grid'

    Adds additional attributes:
    - block_id: Unique identifier for each block
    - area: Area of the block in square meters
    - site_id: ID of the site polygon the block intersects with

    Args:
        output_path: Path to the GeoPackage
        base_polygons_layer: Name of layer containing base polygons to classify
        roads_buffer_layer: Name of layer containing buffered roads with road_type attribute
        dead_end_lines_layer: Name of layer containing dead-end lines
        site_layer: Name of layer containing site polygons (default: '00_site')
        output_layer_name: Name for the output classified blocks layer
        road_type_column: Column name containing road type (default: 'road_type')

    Returns:
        Name of the created classified blocks layer
    """
    output_path = str(output_path)

    print(f"  Loading base polygons from layer: {base_polygons_layer}")
    gdf_blocks = gpd.read_file(output_path, layer=base_polygons_layer)

    # Load site layer to get site_id information
    print(f"  Loading site polygons from layer: {site_layer}")
    try:
        gdf_sites = gpd.read_file(output_path, layer=site_layer)
    except Exception as e:
        print(f"  Warning: Could not load site layer '{site_layer}': {e}")
        gdf_sites = None

    if gdf_blocks.empty:
        raise RuntimeError(f"Layer '{base_polygons_layer}' is empty")

    # Add block attributes
    print("  Adding block attributes...")

    # Add unique block ID
    gdf_blocks["block_id"] = range(len(gdf_blocks))
    gdf_blocks["id"] = range(len(gdf_blocks))

    # Calculate area in square meters
    gdf_blocks["area"] = gdf_blocks.geometry.area

    # Determine site_id by spatial join with site layer
    if gdf_sites is not None and not gdf_sites.empty:
        print("  Determining site_id for each block...")

        # Ensure sites have an ID field
        if "id" not in gdf_sites.columns and "site_id" not in gdf_sites.columns:
            gdf_sites["site_id"] = range(len(gdf_sites))

        site_id_col = "site_id" if "site_id" in gdf_sites.columns else "id"

        # Initialize site_id column
        gdf_blocks["site_id"] = None

        # Find which site each block belongs to
        for block_idx, block in gdf_blocks.iterrows():
            if block.geometry is None or block.geometry.is_empty:
                continue

            # Find the site with maximum overlap
            max_overlap = 0
            assigned_site_id = None

            for _site_idx, site in gdf_sites.iterrows():
                if site.geometry is None or site.geometry.is_empty:
                    continue

                try:
                    if block.geometry.intersects(site.geometry):
                        overlap = block.geometry.intersection(site.geometry).area
                        if overlap > max_overlap:
                            max_overlap = overlap
                            assigned_site_id = site[site_id_col]
                except Exception as e:
                    print(f"    Warning: Error checking intersection: {e}")
                    continue

            if assigned_site_id is not None:
                gdf_blocks.at[block_idx, "site_id"] = assigned_site_id
    else:
        print("  Warning: Site layer not available, site_id will be null")
        gdf_blocks["site_id"] = None

    print(f"  Loading roads buffer from layer: {roads_buffer_layer}")
    gdf_roads = gpd.read_file(output_path, layer=roads_buffer_layer)

    if gdf_roads.empty:
        print("  Warning: Roads buffer layer is empty, classifying all blocks as 'off_grid'")
        gdf_blocks["block_type"] = "off_grid"
        gdf_blocks.to_file(output_path, layer=output_layer_name, driver="GPKG")
        return output_layer_name

    # Separate arterial and secondary roads
    gdf_arterial = gdf_roads[gdf_roads[road_type_column] == "road_art"].copy()
    gdf_secondary = gdf_roads[gdf_roads[road_type_column] == "road_sec"].copy()
    gdf_dead_ends = gpd.read_file(output_path, layer=dead_end_lines_layer)

    print(f"  Found {len(gdf_arterial)} arterial road buffer(s)")
    print(f"  Found {len(gdf_secondary)} secondary road buffer(s)")

    # Initialize all blocks as off_grid
    gdf_blocks["block_type"] = "off_grid"

    if not gdf_dead_ends.empty:
        print("  Classifying blocks adjacent to dead-end lines...")
        for idx, block in gdf_blocks.iterrows():
            if block.geometry is None or block.geometry.is_empty:
                continue
            for _, dead_end in gdf_dead_ends.iterrows():
                if dead_end.geometry is None or dead_end.geometry.is_empty:
                    continue
                if block.geometry.touches(
                    dead_end.geometry.buffer(0.1)
                ) or block.geometry.intersects(dead_end.geometry.buffer(0.1)):
                    gdf_blocks.at[idx, "block_type"] = "cold_boundary"
                    break

    # Check adjacency to arterial roads (highest priority)
    if not gdf_arterial.empty:
        print("  Classifying blocks adjacent to arterial roads...")
        for idx, block in gdf_blocks.iterrows():
            if block.geometry is None or block.geometry.is_empty:
                continue
            if gdf_blocks.at[idx, "block_type"] == "cold_boundary":
                continue
            for _, road in gdf_arterial.iterrows():
                if road.geometry is None or road.geometry.is_empty:
                    continue

                if block.geometry.touches(road.geometry.buffer(0.1)) or block.geometry.intersects(
                    road.geometry.buffer(0.1)
                ):
                    gdf_blocks.at[idx, "block_type"] = "on_grid_art"
                    break

    if not gdf_secondary.empty:
        print("  Classifying blocks adjacent to secondary roads...")
        for idx, block in gdf_blocks.iterrows():
            if block.geometry is None or block.geometry.is_empty:
                continue
            if gdf_blocks.at[idx, "block_type"] == "cold_boundary":
                continue
            for _, road in gdf_secondary.iterrows():
                if road.geometry is None or road.geometry.is_empty:
                    continue
                if block.geometry.touches(road.geometry.buffer(0.1)) or block.geometry.intersects(
                    road.geometry.buffer(0.1)
                ):
                    if gdf_blocks.at[idx, "block_type"] == "on_grid_art":
                        gdf_blocks.at[idx, "block_type"] = "on_grid_art_sec"
                    else:
                        gdf_blocks.at[idx, "block_type"] = "on_grid_sec"
                    break

    # Count classifications
    type_counts = gdf_blocks["block_type"].value_counts()
    print("  Classification results:")
    for block_type, count in type_counts.items():
        print(f"    {block_type}: {count} block(s)")

    # Save classified blocks
    gdf_blocks.to_file(output_path, layer=output_layer_name, driver="GPKG")
    print(f"  Saved classified blocks to layer: {output_layer_name}")

    return output_layer_name


def create_dead_end_lines(
    output_path: Path,
    site_minus_roads_layer: str,
    buffer_roads_layer: str,
    output_layer_name: str = "07_dead_end_lines",
    buffer_distance: float = 0.1,
) -> str:
    """Create dead-end lines by extracting site boundaries and subtracting buffered roads.

    This function:
    1. Extracts the envelope ring (boundary) of each site polygon
    2. Subtracts the buffered roads from these boundaries
    3. Returns the remaining lines (dead-end lines)

    Args:
        output_path: Path to the GeoPackage
        site_minus_roads_layer: Name of layer containing site polygons with roads removed
        buffer_roads_layer: Name of layer containing buffered roads
        output_layer_name: Name for the output dead-end lines layer
        buffer_distance: Buffer distance to apply when subtracting roads (default: 0.1m)

    Returns:
        Name of the created dead-end lines layer
    """
    output_path = str(output_path)

    print(f"  Loading site polygons from layer: {site_minus_roads_layer}")
    gdf_site = gpd.read_file(output_path, layer=site_minus_roads_layer)

    if gdf_site.empty:
        raise RuntimeError(f"Layer '{site_minus_roads_layer}' is empty")

    # Extract envelope ring (boundary) from each site polygon
    print(f"  Extracting envelope rings from {len(gdf_site)} site polygon(s)")
    line_geoms = []
    line_attrs = []

    for idx, row in gdf_site.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        if geom.geom_type == "Polygon":
            line_geoms.append(LineString(geom.exterior.coords))
            line_attrs.append({"site_id": idx, "source": "site_boundary"})
        elif geom.geom_type == "MultiPolygon":
            for poly_idx, poly in enumerate(geom.geoms):
                line_geoms.append(LineString(poly.exterior.coords))
                line_attrs.append(
                    {"site_id": idx, "poly_part": poly_idx, "source": "site_boundary"}
                )

    if not line_geoms:
        raise RuntimeError("No boundary lines extracted from site polygons")

    print(f"  Extracted {len(line_geoms)} boundary line(s)")

    # Create GeoDataFrame with boundary lines
    site_boundary_lines_layer = f"{output_layer_name}_site_boundaries"
    gdf_boundaries = gpd.GeoDataFrame(line_attrs, geometry=line_geoms, crs=gdf_site.crs)
    gdf_boundaries.to_file(output_path, layer=site_boundary_lines_layer, driver="GPKG")
    print(f"  Saved site boundaries to layer: {site_boundary_lines_layer}")

    # Subtract buffered roads from boundary lines
    print("  Subtracting buffered roads from boundaries...")
    result_layer = subtract_layer(
        output_path,
        site_boundary_lines_layer,
        buffer_roads_layer,
        output_path,
        output_layer_name,
        buffer_distance=buffer_distance,
    )

    print(f"  Created dead-end lines layer: {result_layer}")

    return result_layer


def build_base_polygons(
    output_gpkg: Path,
    site_minus_roads_layer: str,
    roads_buffer_layer: str,
    local_roads_layer: str,
    output_layer_name: str = "09_base_polygons",
    line_extension: float = 2.0,
) -> Optional[str]:
    """Create polygons from site-minus-roads boundaries merged with provided local roads.

    Optionally extend local road lines slightly to ensure intersections meet.
    """
    gdf_site_minus = gpd.read_file(output_gpkg, layer=site_minus_roads_layer)
    gdf_local_lines = gpd.read_file(output_gpkg, layer=local_roads_layer)
    gdf_roads_buffer = gpd.read_file(output_gpkg, layer=roads_buffer_layer)
    roads_buffer = unary_union(gdf_roads_buffer.geometry)

    line_geoms = []
    base_crs = gdf_site_minus.crs or gdf_local_lines.crs

    for geom in gdf_site_minus.geometry:
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "Polygon":
            line_geoms.append(LineString(geom.exterior))
            line_geoms.extend(LineString(ring) for ring in geom.interiors if not ring.is_empty)
        elif geom.geom_type == "MultiPolygon":
            for poly in geom.geoms:
                line_geoms.append(LineString(poly.exterior))
                line_geoms.extend(LineString(ring) for ring in poly.interiors if not ring.is_empty)

    for geom in gdf_local_lines.geometry:
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "LineString":
            should_extend = False
            geom_start = geom.interpolate(0)
            geom_end = geom.interpolate(geom.length)
            if geom_start.buffer(0.1).intersects(roads_buffer) or geom_end.buffer(0.1).intersects(
                roads_buffer
            ):
                should_extend = True
            line_geoms.append(extend_line(geom, line_extension) if should_extend else geom)
        elif geom.geom_type == "MultiLineString":
            line_geoms.extend([extend_line(g, line_extension) for g in geom.geoms])

    if not line_geoms:
        print("  Warning: no line geometries available to polygonize")
        return None

    union_lines = unary_union(line_geoms)
    polys = list(polygonize(union_lines))
    gdf_base = gpd.GeoDataFrame(geometry=polys, crs=base_crs)
    gdf_base.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
    print(f"  Created {len(gdf_base)} polygons in layer: {output_layer_name}")
    return output_layer_name
