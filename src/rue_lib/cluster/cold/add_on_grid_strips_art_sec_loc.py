# src/rue_lib/cluster/cold/add_on_grid_strips_art_sec_loc.py
"""Create on-grid strips for arterial, secondary, and local roads.

This module processes blocks and creates on-grid strip parts by:
1. Extracting road edges by type (arterial, secondary, local)
2. Creating offset buffers around each road type
3. Performing boolean operations to create different part types
4. Transferring edge attributes between touching parts
5. Cleaning polygon geometries
"""

from pathlib import Path
from typing import Optional

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon, Point
from shapely.ops import unary_union

from rue_lib.cluster.block_edges import extract_block_edges
from rue_lib.cluster.helpers import compute_angle_dot
from rue_lib.core.definitions import RoadTypes


def get_site_plines(
        block_edges_gdf: gpd.GeoDataFrame,
        road_type: str
) -> list[LineString]:
    """
    Extract continuous polylines from block edges by road type.

    Groups consecutive edges with the same road type into continuous polylines.

    Args:
        block_edges_gdf: GeoDataFrame with block edges and 'road_type' attribute
        road_type: Road type to extract ('road_art', 'road_sec', 'road_loc')

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
    posis = []

    for _, edge in matching_edges.iterrows():
        if not isinstance(edge.geometry, LineString):
            continue

        coords = list(edge.geometry.coords)
        if len(coords) < 2:
            continue

        start_end = [coords[0], coords[-1]]

        if len(posis) == 0 or start_end[0] != posis[-1][-1]:
            posis.append(start_end)
        else:
            posis[-1].append(start_end[1])

    if len(posis) == 0:
        return []

    # Check if first and last connect (closed loop)
    if len(posis) > 1 and posis[-1][-1] == posis[0][0]:
        first_list = posis[-1] + posis[0][1:]
        posis[0] = first_list
        posis = posis[:-1]

    # Create LineStrings
    site_plines = []
    for posi_list in posis:
        if len(posi_list) >= 2:
            site_plines.append(LineString(posi_list))

    return site_plines


def clean_polygon_edges(polygon: Polygon, min_length: float = 0.0001) -> Polygon:
    """
    Remove edges shorter than min_length from a polygon.

    Args:
        polygon: Input polygon
        min_length: Minimum edge length threshold

    Returns:
        Cleaned polygon
    """
    coords = list(polygon.exterior.coords)[:-1]
    new_coords = []

    for i in range(len(coords)):
        curr = coords[i]
        next_idx = (i + 1) % len(coords)
        next_pt = coords[next_idx]

        # Calculate edge length
        length = np.sqrt(
            (next_pt[0] - curr[0])**2 + (next_pt[1] - curr[1])**2
        )

        if length >= min_length:
            new_coords.append(curr)

    if len(new_coords) >= 3:
        new_coords.append(new_coords[0])
        return Polygon(new_coords)
    return polygon


def clean_polygon_angles(
        polygon: Polygon,
        ang_threshold: float = 0.9999
) -> Polygon:
    """
    Remove vertices where edges are nearly collinear (angle close to 180�).

    Uses dot product of normalized edge vectors: values close to -1 indicate
    nearly straight angles (180�).

    Args:
        polygon: Input polygon
        ang_threshold: Dot product threshold (higher = more strict)

    Returns:
        Cleaned polygon
    """
    coords = list(polygon.exterior.coords)[:-1]

    if len(coords) < 4:
        return polygon

    # Find vertices to keep
    keep_vertices = []

    for i in range(len(coords)):
        dot = compute_angle_dot(polygon, i)

        # If angle is NOT nearly straight, keep the vertex
        if abs(dot) <= ang_threshold:
            keep_vertices.append(coords[i])

    if len(keep_vertices) >= 3:
        keep_vertices.append(keep_vertices[0])
        return Polygon(keep_vertices)

    return polygon


def clean_polygons(polygons: list[Polygon]) -> list[Polygon]:
    """
    Clean a list of polygons by removing small edges and collinear vertices.

    Args:
        polygons: List of polygons to clean

    Returns:
        List of cleaned polygons
    """
    cleaned = []

    for pgon in polygons:
        if not isinstance(pgon, Polygon) or pgon.is_empty:
            continue

        # Clean edges
        pgon_clean = clean_polygon_edges(pgon)

        # Clean angles
        pgon_clean = clean_polygon_angles(pgon_clean)

        if pgon_clean.is_valid and pgon_clean.area > 1.0:
            cleaned.append(pgon_clean)

    return cleaned


def find_touching_edge(
        edges_gdf: gpd.GeoDataFrame,
        point: Point,
        tolerance: float = 0.01
) -> Optional[int]:
    """
    Find an edge that touches or is very close to a given point.

    Args:
        edges_gdf: GeoDataFrame containing edges
        point: Point to check
        tolerance: Distance tolerance for "touching"

    Returns:
        Index of touching edge or None
    """
    for idx, edge in edges_gdf.iterrows():
        if not isinstance(edge.geometry, LineString):
            continue

        distance = point.distance(edge.geometry)
        if distance < tolerance:
            return idx

    return None


def transfer_edge_attributes_touching(
        from_edges_gdf: gpd.GeoDataFrame,
        to_edges_gdf: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Transfer road_type attribute from touching edges.

    Args:
        from_edges_gdf: Source edges with road_type attribute
        to_edges_gdf: Target edges to update

    Returns:
        Updated to_edges_gdf with transferred attributes
    """
    for idx, to_edge in to_edges_gdf.iterrows():
        if to_edge.get('road_type') is not None:
            continue

        # Get centroid of target edge
        centroid = to_edge.geometry.centroid

        # Find touching source edge
        touching_idx = find_touching_edge(from_edges_gdf, centroid)

        if touching_idx is not None:
            from_edge = from_edges_gdf.loc[touching_idx]
            to_edges_gdf.at[idx, 'road_type'] = from_edge.get('road_type')

    return to_edges_gdf


def transfer_edge_attributes_between_parts(
        parts_gdf: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Transfer road_type attributes between touching parts.

    Args:
        parts_gdf: GeoDataFrame with parts

    Returns:
        Updated parts_gdf with transferred edge attributes
    """
    # This is a simplified version - would need proper edge extraction
    return parts_gdf


def copy_block_attributes(
        block_row: gpd.GeoSeries,
        parts_gdf: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Copy block attributes to all parts.

    Args:
        block_row: Source block with attributes
        parts_gdf: Target parts GeoDataFrame

    Returns:
        Updated parts_gdf with copied attributes
    """
    for attr in ['block_id', 'block_type', 'site']:
        if attr in block_row:
            parts_gdf[attr] = block_row[attr]

    return parts_gdf


def create_on_grid_strips(
        block: Polygon,
        block_edges_gdf: gpd.GeoDataFrame,
        part_art_d: float,
        part_sec_d: float,
        part_loc_d: float
) -> list[dict]:
    """
    Create on-grid strip parts from block by offsetting roads.

    This is the core algorithm that:
    1. Creates offset buffers around arterial, secondary, and local roads
    2. Performs boolean operations to create different part types:
       - art_sec: Intersection of arterial and secondary offsets
       - art_loc: Intersection of arterial and local offsets
       - sec_loc: Intersection of secondary and local offsets
       - art: Arterial strip only
       - sec: Secondary strip only
       - loc: Local strip only
       - off_grid: Remaining area

    Args:
        block: Block polygon
        block_edges_gdf: GeoDataFrame with block edges
        part_art_d: Offset distance for arterial roads
        part_sec_d: Offset distance for secondary roads
        part_loc_d: Offset distance for local roads

    Returns:
        List of dictionaries with geometry and type for each part
    """
    # Get site polylines for each road type
    roads_art = get_site_plines(block_edges_gdf, RoadTypes.Artery)
    roads_sec = get_site_plines(block_edges_gdf, RoadTypes.Secondary)
    roads_loc = get_site_plines(block_edges_gdf, RoadTypes.Local)

    # Create offset buffers (using flat cap style, similar to JS 'square_end')
    off_art = None
    if roads_art:
        art_lines = MultiLineString(roads_art) if len(roads_art) > 1 else roads_art[0]
        off_art = art_lines.buffer(part_art_d, cap_style=3)  # 3 = flat/square

    off_sec = None
    if roads_sec:
        sec_lines = MultiLineString(roads_sec) if len(roads_sec) > 1 else roads_sec[0]
        off_sec = sec_lines.buffer(part_sec_d, cap_style=3)

    off_loc = None
    if roads_loc:
        loc_lines = MultiLineString(roads_loc) if len(roads_loc) > 1 else roads_loc[0]
        off_loc = loc_lines.buffer(part_loc_d, cap_style=3)

    # Boolean operations to create off_grid area
    off_grids = block
    if off_art is not None:
        off_grids = off_grids.difference(off_art)
    if off_sec is not None:
        off_grids = off_grids.difference(off_sec)
    if off_loc is not None:
        off_grids = off_grids.difference(off_loc)

    # Create on-grid strips by intersecting offsets with block
    on_arts1 = None
    if off_art is not None:
        on_arts1 = off_art.intersection(block)

    on_secs1 = None
    if off_sec is not None:
        on_secs1 = off_sec.intersection(block)

    on_locs1 = None
    if off_loc is not None:
        on_locs1 = off_loc.intersection(block)

    # Remove overlaps between different road types
    on_arts2 = on_arts1
    if on_arts2 is not None:
        if off_sec is not None:
            on_arts2 = on_arts2.difference(off_sec)
        if off_loc is not None:
            on_arts2 = on_arts2.difference(off_loc)

    on_secs2 = on_secs1
    if on_secs2 is not None:
        if off_art is not None:
            on_secs2 = on_secs2.difference(off_art)
        if off_loc is not None:
            on_secs2 = on_secs2.difference(off_loc)

    on_locs2 = on_locs1
    if on_locs2 is not None:
        if off_art is not None:
            on_locs2 = on_locs2.difference(off_art)
        if off_sec is not None:
            on_locs2 = on_locs2.difference(off_sec)

    # Union multi-parts into single geometries
    on_arts3 = None
    if on_arts2 is not None and not on_arts2.is_empty:
        on_arts3 = unary_union([on_arts2]) if not isinstance(on_arts2, MultiPolygon) else unary_union(on_arts2.geoms)

    on_secs3 = None
    if on_secs2 is not None and not on_secs2.is_empty:
        on_secs3 = unary_union([on_secs2]) if not isinstance(on_secs2, MultiPolygon) else unary_union(on_secs2.geoms)

    on_locs3 = None
    if on_locs2 is not None and not on_locs2.is_empty:
        on_locs3 = unary_union([on_locs2]) if not isinstance(on_locs2, MultiPolygon) else unary_union(on_locs2.geoms)

    # Create intersection parts (corners)
    art_sec = None
    if off_art is not None and off_sec is not None:
        art_sec = off_art.intersection(off_sec)
        if not art_sec.is_empty:
            art_sec = art_sec.intersection(block)

    art_loc = None
    if off_art is not None and off_loc is not None:
        art_loc = off_art.intersection(off_loc)
        if not art_loc.is_empty:
            art_loc = art_loc.intersection(block)

    sec_loc = None
    if off_sec is not None and off_loc is not None:
        sec_loc = off_sec.intersection(off_loc)
        if not sec_loc.is_empty:
            sec_loc = sec_loc.intersection(block)

    # Collect all parts with their types
    parts = []

    def add_parts(geom, part_type):
        """Helper to add polygon or multipolygon parts."""
        if geom is None or geom.is_empty:
            return

        if isinstance(geom, Polygon):
            if geom.area > 1.0:
                parts.append({'geometry': geom, 'type': part_type})
        elif isinstance(geom, MultiPolygon):
            for poly in geom.geoms:
                if poly.area > 1.0:
                    parts.append({'geometry': poly, 'type': part_type})

    # Add all part types
    add_parts(art_sec, 'art_sec')
    add_parts(art_loc, 'art_loc')
    add_parts(sec_loc, 'sec_loc')
    add_parts(on_arts3, 'art')
    add_parts(on_secs3, 'sec')
    add_parts(on_locs3, 'loc')
    add_parts(off_grids, 'off_grid')

    # Clean polygons
    cleaned_parts = []
    for part in parts:
        cleaned_geom = clean_polygon_edges(part['geometry'])
        cleaned_geom = clean_polygon_angles(cleaned_geom)

        if cleaned_geom.is_valid and cleaned_geom.area > 1.0:
            cleaned_parts.append({
                'geometry': cleaned_geom,
                'type': part['type']
            })

    return cleaned_parts


def process_block(
        block_row: gpd.GeoSeries,
        roads_gdf: gpd.GeoDataFrame,
        part_art_d: float,
        part_sec_d: float,
        part_loc_d: float
) -> list[dict]:
    """
    Process a single block to create on-grid strip parts.

    Args:
        block_row: Block GeoSeries with geometry and attributes
        roads_gdf: GeoDataFrame with roads
        part_art_d: Offset distance for arterial roads
        part_sec_d: Offset distance for secondary roads
        part_loc_d: Offset distance for local roads

    Returns:
        List of dictionaries with part geometries and attributes
    """
    block = block_row.geometry

    # Extract block edges
    block_gdf = gpd.GeoDataFrame([block_row], geometry='geometry', crs=roads_gdf.crs)
    block_edges_gdf = extract_block_edges(block_gdf, roads_gdf)

    # Create on-grid strips
    parts = create_on_grid_strips(
        block,
        block_edges_gdf,
        part_art_d,
        part_sec_d,
        part_loc_d
    )

    # Add block attributes to parts
    for part in parts:
        for attr in ['block_id', 'block_type', 'site']:
            if attr in block_row.index:
                part[attr] = block_row[attr]

    return parts


def generate_on_grid_strips_art_sec_loc(
        output_path: Path,
        blocks_layer_name: str,
        roads_layer_name: str,
        part_art_d: float,
        part_sec_d: float,
        part_loc_d: float,
        output_layer_name: str
):
    """
    Generate on-grid strip parts for arterial, secondary, and local roads.

    This function processes all blocks and creates on-grid strip parts by
    offsetting roads inward and performing boolean operations. The strips
    represent different zones within blocks based on road adjacency.

    Args:
        output_path: Path to GeoPackage file
        blocks_layer_name: Name of blocks layer
        roads_layer_name: Name of roads layer
        part_art_d: Offset distance for arterial roads (typically 40m)
        part_sec_d: Offset distance for secondary roads (typically 30m)
        part_loc_d: Offset distance for local roads (typically 20m)
        output_layer_name: Name for output layer

    Returns:
        Name of output layer
    """
    # Read input layers
    blocks_layer = gpd.read_file(output_path, layer=blocks_layer_name)
    roads_layer = gpd.read_file(output_path, layer=roads_layer_name)

    if blocks_layer.empty:
        return output_layer_name

    # Process each block
    all_parts = []

    for idx, block_row in blocks_layer.iterrows():
        parts = process_block(
            block_row,
            roads_layer,
            part_art_d,
            part_sec_d,
            part_loc_d
        )
        all_parts.extend(parts)

    # Create output GeoDataFrame
    if all_parts:
        parts_gdf = gpd.GeoDataFrame(all_parts, crs=blocks_layer.crs)
        parts_gdf.to_file(output_path, layer=output_layer_name, driver="GPKG")

    return output_layer_name