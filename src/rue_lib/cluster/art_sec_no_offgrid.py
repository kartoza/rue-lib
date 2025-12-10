# src/rue_lib/cluster/art_sec_no_offgrid.py
"""Generate arterial and secondary parts without off-grid subdivision.

This module creates block parts (corner and side parts) for arterial and
secondary road blocks by offsetting roads inward and performing boolean
operations to create the different part types.
"""

from pathlib import Path
from typing import Optional

import geopandas as gpd
import numpy as np
import pyogrio.errors
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon

from rue_lib.cluster.block_edges import extract_block_edges
from rue_lib.cluster.classification import classify_part_type
from rue_lib.cluster.off_grid import extend_line
from rue_lib.core.definitions import RoadTypes, ColorTypes


def get_perimeter_lines_by_road_type(
        block: Polygon,
        block_edges_gdf: gpd.GeoDataFrame,
        road_type: str
) -> list[LineString]:
    """
    Get perimeter lines from block edges matching a specific road type.

    Args:
        block: Block polygon
        block_edges_gdf: GeoDataFrame with block edges and 'road_type' attribute
        road_type: Type of road to filter ('road_art', 'road_sec', 'road_loc')

    Returns:
        List of LineStrings representing continuous road segments
    """
    # Filter edges by road type
    matching_edges = block_edges_gdf[
        block_edges_gdf.get('road_type', '') == road_type
        ]
    if matching_edges.empty:
        return []

    # Group consecutive edges into polylines
    plines = []
    for _, edge in matching_edges.iterrows():
        if isinstance(edge.geometry, LineString):
            # Extend the line by 100 units (as in JS code)
            extended = extend_line(edge.geometry, 100.0)
            plines.append(extended)

    return plines


def get_perpendicular_local_lines(
        block: Polygon,
        block_edges_gdf: gpd.GeoDataFrame,
        block_type: str,
        ortho_direction: np.ndarray
) -> list[LineString]:
    """
    Get local road edges that are perpendicular to the primary direction.

    Args:
        block: Block polygon
        block_edges_gdf: GeoDataFrame with block edges
        block_type: Block type ('art' or 'sec')
        ortho_direction: Orthogonal direction vector

    Returns:
        List of perpendicular local road LineStrings
    """
    perp_dir = ortho_direction.copy()

    loc_perp_edges = []

    for _, edge in block_edges_gdf.iterrows():
        if edge.get('road_type') != 'road_loc':
            continue

        # Get edge direction vector
        coords = list(edge.geometry.coords)
        if len(coords) < 2:
            continue

        edge_vec = np.array([
            coords[-1][0] - coords[0][0],
            coords[-1][1] - coords[0][1],
            0.0
        ])

        edge_vec_norm = edge_vec / (np.linalg.norm(edge_vec) + 1e-10)

        # Check if edge is perpendicular (dot product close to 0.8 or higher)
        if abs(np.dot(edge_vec_norm, perp_dir)) > 0.8:
            loc_perp_edges.append(edge.geometry)

    # Extend all lines
    extended = [extend_line(line, 100.0) for line in loc_perp_edges]

    return extended


def trim_parts_with_local_roads(
        block: Polygon,
        parts: list[Polygon],
        local_depth: float,
        block_edges_gdf: gpd.GeoDataFrame,
        ortho_direction: np.ndarray
) -> list[Polygon]:
    """
    Trim parts using local roads offset.

    Args:
        block: Block polygon
        parts: List of part polygons to trim
        local_depth: Offset distance for local roads
        block_edges_gdf: GeoDataFrame with block edges
        ortho_direction: Orthogonal direction for perpendicular check

    Returns:
        Trimmed list of part polygons
    """
    # Get perpendicular local road lines
    plines_loc = get_perpendicular_local_lines(
        block, block_edges_gdf, '', ortho_direction
    )

    if not plines_loc:
        return parts

    # Create offset buffer around local roads
    loc_lines = MultiLineString(plines_loc) if len(plines_loc) > 1 else \
        plines_loc[0]
    loc_off = loc_lines.buffer(local_depth, cap_style='flat')

    # Trim parts by subtracting local road buffer
    parts_trim = []

    for part in parts:
        try:
            trimmed = part.difference(loc_off)

            if isinstance(trimmed, Polygon):
                parts_trim.append(trimmed)
            elif isinstance(trimmed, MultiPolygon):
                parts_trim.extend(list(trimmed.geoms))
        except Exception:
            parts_trim.append(part)

        # Also create intersection parts (local road parts)
        try:
            loc_part = loc_off.intersection(part)

            if isinstance(loc_part, Polygon) and loc_part.area > 1.0:
                parts_trim.append(loc_part)
            elif isinstance(loc_part, MultiPolygon):
                parts_trim.extend([p for p in loc_part.geoms if p.area > 1.0])
        except Exception:
            pass

    return parts_trim


def create_parts_from_block(
        block: Polygon,
        block_edges_gdf: gpd.GeoDataFrame,
        depths_dict: dict[str, float],
        ortho_direction: np.ndarray
) -> list[Polygon]:
    """
    Create parts from a block based on arterial and secondary roads.

    Args:
        block: Block polygon
        block_edges_gdf: GeoDataFrame with block edges and road_type attributes
        depths_dict: Dictionary mapping road types to offset depths
        ortho_direction: Orthogonal direction vector

    Returns:
        List of part polygons
    """
    # Get perimeter lines for arterial and secondary roads
    plines_art = get_perimeter_lines_by_road_type(
        block, block_edges_gdf, RoadTypes.Artery
    )
    plines_sec = get_perimeter_lines_by_road_type(
        block, block_edges_gdf, RoadTypes.Secondary
    )

    part_art_d = depths_dict.get(RoadTypes.Artery, 40)
    part_sec_d = depths_dict.get(RoadTypes.Secondary, 30)
    part_loc_d = depths_dict.get(RoadTypes.Local, 20)

    # Case 1: Only arterial roads
    if len(plines_art) > 0 and len(plines_sec) == 0:
        if not plines_art:
            return []

        # Create offset buffer
        art_lines = MultiLineString(plines_art) if len(plines_art) > 1 else \
            plines_art[0]
        art_off = art_lines.buffer(part_art_d, cap_style='flat')

        # Boolean operations
        try:
            art_bool1 = art_off.intersection(block)
            art_bool2 = block.difference(art_off)

            parts = []
            if isinstance(art_bool1, Polygon):
                parts.append(art_bool1)
            elif isinstance(art_bool1, MultiPolygon):
                parts.extend(list(art_bool1.geoms))

            if isinstance(art_bool2, Polygon):
                parts.append(art_bool2)
            elif isinstance(art_bool2, MultiPolygon):
                parts.extend(list(art_bool2.geoms))

            # Trim with local roads if needed
            result = trim_parts_with_local_roads(
                block, parts, part_loc_d, block_edges_gdf, ortho_direction
            )

            return [p for p in result if p.area > 1.0]

        except Exception:
            return []

    # Case 2: Only secondary roads
    elif len(plines_art) == 0 and len(plines_sec) > 0:
        if not plines_sec:
            return []

        # Create offset buffer
        sec_lines = MultiLineString(plines_sec) if len(plines_sec) > 1 else \
            plines_sec[0]
        sec_off = sec_lines.buffer(part_sec_d, cap_style='flat')

        # Boolean operations
        try:
            sec_bool1 = sec_off.intersection(block)
            sec_bool2 = block.difference(sec_off)

            parts = []
            if isinstance(sec_bool1, Polygon):
                parts.append(sec_bool1)
            elif isinstance(sec_bool1, MultiPolygon):
                parts.extend(list(sec_bool1.geoms))

            if isinstance(sec_bool2, Polygon):
                parts.append(sec_bool2)
            elif isinstance(sec_bool2, MultiPolygon):
                parts.extend(list(sec_bool2.geoms))

            # Trim with local roads
            result = trim_parts_with_local_roads(
                block, parts, part_loc_d, block_edges_gdf, ortho_direction
            )

            return [p for p in result if p.area > 1.0]

        except Exception:
            return []

    # Case 3: Both arterial and secondary roads
    elif len(plines_art) > 0 and len(plines_sec) > 0:
        try:
            # Create offset buffers
            art_lines = MultiLineString(plines_art) if len(plines_art) > 1 else \
                plines_art[0]
            art_off = art_lines.buffer(part_art_d, cap_style='flat')

            sec_lines = MultiLineString(plines_sec) if len(plines_sec) > 1 else \
                plines_sec[0]
            sec_off = sec_lines.buffer(part_sec_d, cap_style='flat')

            # Boolean operations
            art_bool1 = art_off.intersection(block)
            art_bool2 = block.difference(art_off)

            # Create corner and other parts
            art_corners = art_bool1.intersection(sec_off) if hasattr(
                art_bool1, 'intersection') else Polygon()
            other_corners = art_bool2.intersection(sec_off) if hasattr(
                art_bool2, 'intersection') else Polygon()

            art_bool1_trim = art_bool1.difference(sec_off) if hasattr(
                art_bool1, 'difference') else art_bool1
            art_bool2_trim = art_bool2.difference(sec_off) if hasattr(
                art_bool2, 'difference') else art_bool2

            # Check if art_bool2_trim area is significant
            art_bool2_trim_area = art_bool2_trim.area if hasattr(
                art_bool2_trim, 'area') else 0

            # Collect all parts
            all_parts = []

            for geom in [
                art_corners, other_corners, art_bool1_trim,
                art_bool2_trim
            ]:
                if isinstance(geom, Polygon) and geom.area > 1.0:
                    all_parts.append(geom)
                elif isinstance(geom, MultiPolygon):
                    all_parts.extend([p for p in geom.geoms if p.area > 1.0])

            # Only trim with local roads if art_bool2_trim area is large enough
            if art_bool2_trim_area > (part_loc_d * part_loc_d * 4):
                result = trim_parts_with_local_roads(
                    block, all_parts, part_loc_d, block_edges_gdf,
                    ortho_direction
                )
                return [p for p in result if p.area > 1.0]
            else:
                return all_parts
        except Exception:
            return []
    return []


def generate_art_sec_parts_no_offgrid(
        output_path: Path,
        blocks_layer_name: str,
        roads_layer_name: str,
        part_art_d: float,
        part_sec_d: float,
        part_loc_d: float,
        output_layer_name,
        ortho_direction: Optional[np.ndarray] = None,
):
    """
    Generate arterial and secondary block parts without off-grid subdivision.

    This function processes blocks and creates parts based on arterial and
    secondary road adjacency. It performs geometric operations to create
    corner parts, side parts, and other part types by offsetting roads inward
    and performing boolean operations. Results are saved to a GeoPackage file.

    Args:
        output_path: Path to the GeoPackage file containing input layers and
            where output will be saved.
        blocks_layer_name: Name of the layer containing block polygons.
        roads_layer_name: Name of the layer containing road geometries.
        road_local_width_m: Width of local roads in meters, used for extracting
            block edges.
        part_art_d: Offset depth for arterial roads in meters (typically 40m).
        part_sec_d: Offset depth for secondary roads in meters (typically 30m).
        part_loc_d: Offset depth for local roads in meters (typically 20m).
        output_layer_name: Name for the output layer containing generated parts.
        ortho_direction: Optional orthogonal direction vector for perpendicular
            checks. Defaults to [0, 1, 0] if not provided.
    """
    try:
        blocks_layer = gpd.read_file(
            output_path, layer=blocks_layer_name
        )
    except pyogrio.errors.DataLayerError:
        print(f"No off grid on on grid.")
        return

    roads_layer = gpd.read_file(
        output_path, layer=roads_layer_name
    )
    block_edges_gdf = extract_block_edges(blocks_layer, roads_layer)

    if ortho_direction is None:
        ortho_direction = np.array([0, 1, 0])

    if blocks_layer.empty:
        return gpd.GeoDataFrame(
            columns=['geometry', 'class', 'type', 'block_id'])

    depths_dict = {
        RoadTypes.Artery: part_art_d,
        RoadTypes.Secondary: part_sec_d,
        RoadTypes.Local: part_loc_d,
        'cold': 0
    }

    all_parts = []
    edges_for_blocks = []
    edges_for_parts = []

    for idx, block_row in blocks_layer.iterrows():
        block = block_row.geometry

        # Get edges for this block
        block_id = block_row.get('block_id', idx)
        edges_for_block = block_edges_gdf[
            block_edges_gdf.get('block_id', -1) == block_id
            ]

        # Create parts
        parts = create_parts_from_block(
            block, edges_for_block, depths_dict, ortho_direction
        )

        # TODO:
        #  Remove: Save it to layer
        for idx, edge in edges_for_block.iterrows():
            edges_for_blocks.append(edge)

        for part in parts:
            # Convert part polygon to GeoDataFrame
            part_gdf = gpd.GeoDataFrame(
                [{'geometry': part, 'block_id': block_id}],
                crs=blocks_layer.crs
            )
            part_edges = extract_block_edges(part_gdf, roads_layer)
            part_type = classify_part_type(
                part, part_edges,
                part_art_d=part_art_d,
                part_sec_d=part_sec_d,
                part_loc_d=part_loc_d
            )

            # TODO:
            #  Remove: Save it to layer
            for idx, edge in part_edges.iterrows():
                edges_for_parts.append(edge)

            try:
                color = ColorTypes[part_type]
            except KeyError:
                color = "rgb(255,255,255)"

            all_parts.append({
                'geometry': part,
                'class': 'part',
                'block_id': block_id,
                'block_type': block_row.get('block_type', ''),
                'site': block_row.get('site', ''),
                'type': part_type,
                'color': color,
            })
    # Create GeoDataFrame
    gdf_out = gpd.GeoDataFrame(edges_for_blocks, crs=blocks_layer.crs)
    gdf_out.to_file(
        output_path, layer=output_layer_name + "-edges", driver="GPKG"
    )
    gdf_out = gpd.GeoDataFrame(edges_for_parts, crs=blocks_layer.crs)
    gdf_out.to_file(
        output_path, layer=output_layer_name + "-edges-parts", driver="GPKG"
    )

    # Create GeoDataFrame
    gdf_out = gpd.GeoDataFrame(all_parts, crs=blocks_layer.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")
    return output_layer_name
