# src/rue_lib/cluster/art_sec_no_offgrid.py
"""Generate arterial and secondary parts without off-grid subdivision.

This module creates block parts (corner and side parts) for arterial and
secondary road blocks by offsetting roads inward and performing boolean
operations to create the different part types.
"""

from typing import Optional

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon

from rue_lib.cluster.off_grid import extend_line


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
    loc_lines = MultiLineString(plines_loc) if len(plines_loc) > 1 else plines_loc[0]
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
        block, block_edges_gdf, 'road_art'
    )
    plines_sec = get_perimeter_lines_by_road_type(
        block, block_edges_gdf, 'road_sec'
    )

    part_art_d = depths_dict.get('road_art', 40)
    part_sec_d = depths_dict.get('road_sec', 30)
    part_loc_d = depths_dict.get('road_loc', 20)

    # Case 1: Only arterial roads
    if len(plines_art) > 0 and len(plines_sec) == 0:
        if not plines_art:
            return []

        # Create offset buffer
        art_lines = MultiLineString(plines_art) if len(plines_art) > 1 else plines_art[0]
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
        sec_lines = MultiLineString(plines_sec) if len(plines_sec) > 1 else plines_sec[0]
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
            art_lines = MultiLineString(plines_art) if len(plines_art) > 1 else plines_art[0]
            art_off = art_lines.buffer(part_art_d, cap_style='flat')

            sec_lines = MultiLineString(plines_sec) if len(plines_sec) > 1 else plines_sec[0]
            sec_off = sec_lines.buffer(part_sec_d, cap_style='flat')

            # Boolean operations
            art_bool1 = art_off.intersection(block)
            art_bool2 = block.difference(art_off)

            # Create corner and other parts
            art_corners = art_bool1.intersection(sec_off) if hasattr(art_bool1, 'intersection') else Polygon()
            other_corners = art_bool2.intersection(sec_off) if hasattr(art_bool2, 'intersection') else Polygon()

            art_bool1_trim = art_bool1.difference(sec_off) if hasattr(art_bool1, 'difference') else art_bool1
            art_bool2_trim = art_bool2.difference(sec_off) if hasattr(art_bool2, 'difference') else art_bool2

            # Check if art_bool2_trim area is significant
            art_bool2_trim_area = art_bool2_trim.area if hasattr(art_bool2_trim, 'area') else 0

            # Collect all parts
            all_parts = []

            for geom in [art_corners, other_corners, art_bool1_trim, art_bool2_trim]:
                if isinstance(geom, Polygon) and geom.area > 1.0:
                    all_parts.append(geom)
                elif isinstance(geom, MultiPolygon):
                    all_parts.extend([p for p in geom.geoms if p.area > 1.0])

            # Only trim with local roads if art_bool2_trim area is large enough
            if art_bool2_trim_area > (part_loc_d * part_loc_d * 4):
                result = trim_parts_with_local_roads(
                    block, all_parts, part_loc_d, block_edges_gdf, ortho_direction
                )
                return [p for p in result if p.area > 1.0]
            else:
                return all_parts

        except Exception:
            return []

    return []


def get_edge_angle(coords: list, vertex_idx: int) -> float:
    """
    Calculate the angle at a vertex.

    Args:
        coords: List of polygon coordinates
        vertex_idx: Index of vertex

    Returns:
        Angle in degrees
    """
    n = len(coords) - 1

    p_prev = np.array(coords[(vertex_idx - 1) % n])
    p_curr = np.array(coords[vertex_idx % n])
    p_next = np.array(coords[(vertex_idx + 1) % n])

    vec0 = p_prev - p_curr
    vec1 = p_next - p_curr

    vec0_norm = vec0 / (np.linalg.norm(vec0) + 1e-10)
    vec1_norm = vec1 / (np.linalg.norm(vec1) + 1e-10)

    vec0_rev = -vec0_norm

    # Calculate angle
    dot = np.clip(np.dot(vec1_norm, vec0_rev), -1.0, 1.0)
    angle_rad = np.arccos(dot)

    # Convert to degrees
    angle_deg = np.degrees(angle_rad)

    return angle_deg


def classify_part_type(
    part: Polygon,
    part_edges_gdf: gpd.GeoDataFrame
) -> str:
    """
    Classify part type based on curved road edges.

    Args:
        part: Part polygon
        part_edges_gdf: GeoDataFrame with part edges and road_type attributes

    Returns:
        Part type string ('art', 'sec', 'loc', 'art_sec', etc.)
    """
    road_types = []

    # Check each road type
    for road_type_check in ['road_art', 'road_sec', 'road_loc']:
        # Get edges matching this road type
        matching_edges = part_edges_gdf[
            part_edges_gdf.get('road_type', '') == road_type_check
        ]

        if matching_edges.empty:
            continue

        # Group consecutive edges
        edge_groups = [[]]
        prev_found = False

        for _, edge in part_edges_gdf.iterrows():
            if edge.get('road_type') == road_type_check:
                edge_groups[-1].append(edge)
                prev_found = True
            else:
                if prev_found and len(edge_groups[-1]) > 0:
                    edge_groups.append([])
                prev_found = False

        # Remove empty groups
        edge_groups = [g for g in edge_groups if len(g) > 0]

        if not edge_groups:
            continue

        # Merge first and last if they connect
        if len(edge_groups) == 2:
            edges = edge_groups[1] + edge_groups[0]
        else:
            edges = edge_groups[0] if edge_groups else []

        if len(edges) == 0:
            continue

        if len(edges) == 1:
            road_types.append(road_type_check)
            continue

        # Check angles at vertices to count distinct road segments
        coords = list(part.exterior.coords)

        # For multiple edges, check angles to see if they're separate segments
        segment_count = 1
        for i in range(1, len(edges)):
            # Get vertex angle
            try:
                angle = get_edge_angle(coords, i)
                if abs(angle - 180) > 45:
                    segment_count += 1
            except Exception:
                pass

        road_types.extend([road_type_check] * segment_count)

    # Map to block type
    block_types_dict = {
        'road_art': 'art',
        'road_sec': 'sec',
        'road_ter': 'ter',
        'road_loc': 'loc'
    }

    if len(road_types) == 0:
        return 'off_grid'
    elif len(road_types) == 1:
        return block_types_dict.get(road_types[0], 'off_grid')
    else:
        # Multiple road types
        type0 = block_types_dict.get(road_types[0], '')
        type1 = block_types_dict.get(road_types[1], '')
        return f"{type0}_{type1}" if type0 and type1 else type0


def generate_art_sec_parts_no_offgrid(
    blocks_gdf: gpd.GeoDataFrame,
    block_edges_gdf: gpd.GeoDataFrame,
    part_art_d: float = 40.0,
    part_sec_d: float = 30.0,
    part_loc_d: float = 20.0,
    ortho_direction: Optional[np.ndarray] = None
) -> gpd.GeoDataFrame:
    """
    Generate arterial and secondary block parts without off-grid subdivision.

    This function processes blocks and creates parts based on arterial and
    secondary road adjacency. It performs geometric operations to create
    corner parts, side parts, and other part types.

    Args:
        blocks_gdf: GeoDataFrame with blocks (must have 'block_type' column)
        block_edges_gdf: GeoDataFrame with block edges (must have 'road_type' column)
        part_art_d: Offset depth for arterial roads (default: 40m)
        part_sec_d: Offset depth for secondary roads (default: 30m)
        part_loc_d: Offset depth for local roads (default: 20m)
        ortho_direction: Orthogonal direction vector for perpendicular checks

    Returns:
        GeoDataFrame with generated parts
    """
    if ortho_direction is None:
        ortho_direction = np.array([0, 1, 0])

    # Filter to only arterial and secondary blocks
    art_blocks = blocks_gdf[blocks_gdf['block_type'] == 'art']
    sec_blocks = blocks_gdf[blocks_gdf['block_type'] == 'sec']

    blocks_to_process = gpd.GeoDataFrame(
        pd.concat([art_blocks, sec_blocks], ignore_index=True)
    ) if not art_blocks.empty or not sec_blocks.empty else gpd.GeoDataFrame()

    if blocks_to_process.empty:
        return gpd.GeoDataFrame(columns=['geometry', 'class', 'type', 'block_id'])

    depths_dict = {
        'road_art': part_art_d,
        'road_sec': part_sec_d,
        'road_loc': part_loc_d,
        'cold': 0
    }

    areas_dict = {
        'art_art': part_art_d * part_art_d,
        'art_sec': part_art_d * part_sec_d,
        'art_loc': part_art_d * part_loc_d,
        'sec_sec': part_sec_d * part_sec_d,
        'sec_loc': part_sec_d * part_loc_d,
    }

    all_parts = []

    for idx, block_row in blocks_to_process.iterrows():
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

        for part in parts:
            # Get edges for this part (simplified - in real implementation would need spatial join)
            part_edges = edges_for_block  # Simplified

            # Classify part type
            part_type = classify_part_type(part, part_edges)

            # Check if corner part should be simplified to side part
            corner_expected_area = areas_dict.get(part_type)
            if corner_expected_area is not None:
                if part.area > (corner_expected_area * 3):
                    # It's too large for a corner, make it a side part
                    part_type = part_type[:3]  # Take first 3 chars (e.g., 'art_sec' -> 'art')

            # Check for local road length
            if len(part_type) == 7:  # e.g., 'art_loc'
                # Check local road edge lengths
                loc_edges = part_edges[part_edges.get('road_type') == 'road_loc']
                if not loc_edges.empty:
                    max_length = loc_edges.geometry.length.max()
                    if max_length > (part_loc_d * 2):
                        part_type = part_type[:3]

                # Check secondary road edge lengths
                sec_edges = part_edges[part_edges.get('road_type') == 'road_sec']
                if not sec_edges.empty:
                    max_length = sec_edges.geometry.length.max()
                    if max_length > (part_sec_d * 2):
                        part_type = part_type[:3]

            all_parts.append({
                'geometry': part,
                'class': 'part',
                'type': part_type,
                'block_id': block_id,
                'block_type': block_row.get('block_type', ''),
                'site': block_row.get('site', '')
            })

    # Create GeoDataFrame
    if all_parts:
        parts_gdf = gpd.GeoDataFrame(all_parts, crs=blocks_gdf.crs)
        return parts_gdf
    else:
        return gpd.GeoDataFrame(columns=['geometry', 'class', 'type', 'block_id'])