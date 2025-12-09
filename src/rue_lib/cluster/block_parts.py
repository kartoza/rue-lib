# src/rue_lib/cluster/block_parts.py
"""Create corner and side block parts from off-grid polygons."""

from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.geometry import MultiPolygon, Polygon

from rue_lib.core.geometry import remove_vertices_by_angle


def vector_length(v: np.ndarray) -> float:
    """
    Calculate the Euclidean length of a vector.

    Args:
        v: Input vector (any dimension)

    Returns:
        Length (magnitude) of the vector
    """
    return np.linalg.norm(v)


def normalize_vector(v: np.ndarray) -> np.ndarray:
    """
    Normalize a vector to unit length.

    Args:
        v: Input vector to normalize

    Returns:
        Unit vector in the same direction, or original vector if length is 0
    """
    length = vector_length(v)
    if length == 0:
        return v
    return v / length


def set_vector_length(v: np.ndarray, length: float) -> np.ndarray:
    """
    Scale a vector to a specific length.

    Args:
        v: Input vector
        length: Target length for the vector

    Returns:
        Vector with same direction but specified length
    """
    return normalize_vector(v) * length


def compute_internal_angle(
        coords: list[tuple[float, float]], idx: int
) -> float:
    """
    Compute internal angle at vertex idx in a polygon.

    Args:
        coords: list of (x, y) coordinates
        idx: Index of vertex to compute angle at

    Returns:
        Angle in degrees (0-360)
    """
    n = len(coords) - 1

    p_prev = np.array(coords[(idx - 1) % n])
    p_curr = np.array(coords[idx % n])
    p_next = np.array(coords[(idx + 1) % n])

    v1 = p_prev - p_curr
    v2 = p_next - p_curr

    v1_norm = normalize_vector(v1)
    v2_norm = normalize_vector(v2)

    dot = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
    angle_rad = np.arccos(dot)
    angle_deg = np.degrees(angle_rad)

    cross = np.cross(v1, v2)
    if cross < 0:
        angle_deg = 360 - angle_deg

    return angle_deg


def get_edge_direction(
        coords: list[tuple[float, float]], idx: int
) -> np.ndarray:
    """
    Get direction vector of edge starting at vertex idx.

    Args:
        coords: list of (x, y) coordinates
        idx: Index of starting vertex

    Returns:
        Normalized direction vector [dx, dy, 0]
    """
    n = len(coords) - 1
    p_curr = np.array(coords[idx % n])
    p_next = np.array(coords[(idx + 1) % n])

    direction = p_next - p_curr
    direction_3d = np.array([direction[0], direction[1], 0.0])

    return normalize_vector(direction_3d)


def create_corner_fan(
        origin: np.ndarray, dir0: np.ndarray, dir1: np.ndarray,
        distance: float = 1000.0,
        rorate_90: bool = True, rotate_180: bool = False
) -> Polygon:
    """
    Create a corner fan polygon at a sharp corner.

    Args:
        origin: Corner vertex position [x, y, z]
        dir0: Direction of first edge
        dir1: Direction of second edge
        distance: Distance to extend perpendiculars

    Returns:
        Corner fan polygon
    """
    if rorate_90:
        # Compute perpendicular vectors (rotate 90° counterclockwise)
        v0_perp = -np.array([dir0[1], -dir0[0], 0.0])
        v1_perp = -np.array([dir1[1], -dir1[0], 0.0])
    elif rotate_180:
        v0_perp = np.array([dir0[0], dir0[1], 0.0])
        v1_perp = -np.array([dir1[0], dir1[1], 0.0])

    v0_perp = set_vector_length(v0_perp, distance)
    v1_perp = set_vector_length(v1_perp, distance)

    v_mid_perp = v0_perp + v1_perp
    v_mid_perp = set_vector_length(v_mid_perp, distance)

    xyz0 = origin[:2] + v0_perp[:2]
    xyz1 = origin[:2] + v1_perp[:2]
    xyz_mid = origin[:2] + v_mid_perp[:2]

    fan_coords = [
        tuple(origin[:2]), tuple(xyz0), tuple(xyz_mid), tuple(xyz1),
        tuple(origin[:2])
    ]

    return Polygon(fan_coords)


def get_positions_between_edges(
        coords: list[tuple[float, float]], start_idx: int, end_idx: int,
        include_end: bool = True
) -> list[tuple[float, float]]:
    """
    Get all vertex positions along the perimeter between two edge indices.

    Args:
        coords: Polygon coordinates
        start_idx: Starting edge index
        end_idx: Ending edge index
        include_end: Whether to include end vertex

    Returns:
        list of coordinates along the perimeter
    """
    n = len(coords) - 1  # Exclude duplicate last point
    positions = []

    # Start from the vertex after start_idx
    current = (start_idx + 1) % n

    while current != end_idx:
        positions.append(coords[current])
        current = (current + 1) % n

    if include_end:
        positions.append(coords[end_idx])

    return positions


def create_side_polygon(
        a_origin: tuple[float, float],
        a_outer: tuple[float, float],
        b_outer: tuple[float, float],
        b_origin: tuple[float, float],
        perimeter_positions: list[tuple[float, float]],
) -> Polygon:
    """
    Create a side polygon between two corner fans.

    The polygon goes: A-origin -> A-outer -> B-outer -> B-origin -> perimeter (reversed)

    Args:
        a_origin: Origin point of first corner
        a_outer: Outer point of first corner fan
        b_outer: Outer point of second corner fan
        b_origin: Origin point of second corner
        perimeter_positions: Positions along the perimeter between corners

    Returns:
        Side polygon
    """
    # Build coordinate list
    side_coords = [
        a_origin,
        a_outer,
        b_outer,
        b_origin,
    ]

    # Add reversed perimeter positions
    side_coords.extend(reversed(perimeter_positions))

    # Close the polygon
    side_coords.append(a_origin)

    return Polygon(side_coords)


def get_side_parts(
        block_id: int,
        block: Polygon,
        off_grid: Polygon,
        corner_parts: list[Polygon],
        simplify_tolerance: float = 0.1,
) -> list[Polygon]:
    """
    Get side block parts by subtracting off-grid and corners from block.

    Computes: side_parts = block - off_grid - corners
    Then removes spikes by simplifying the geometry.

    Args:
        block: Original block polygon (outer boundary)
        off_grid: Off-grid polygon (inner part of block)
        corner_parts: List of corner polygons to subtract
        simplify_tolerance: Tolerance for simplifying geometry to remove spikes

    Returns:
        list of side strip polygons with spikes removed

    Example:
        >>> block = Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
        >>> off_grid = Polygon([(20, 20), (80, 20), (80, 80), (20, 80)])
        >>> corners = [...]
        >>> sides = get_side_parts(block, off_grid, corners)
    """
    try:
        result = block.difference(off_grid)

        # Subtract each corner part
        for corner in corner_parts:
            result = result.difference(corner)

        # Clean the polygon
        if result.geom_type == "Polygon":
            result = remove_vertices_by_angle(
                result, min_angle_threshold=10
            )
        elif result.geom_type == "MultiPolygon":
            cleaned_geoms = []
            for geom in result.geoms:
                clean_polygon = remove_vertices_by_angle(
                    geom, min_angle_threshold=10
                )
                cleaned_geoms.append(clean_polygon)
            result = MultiPolygon(cleaned_geoms)

        # Simplify to remove unnecessary vertices
        if not result.is_empty and result.area > 0:
            result = result.simplify(
                simplify_tolerance, preserve_topology=True
            )

        # Handle result which could be Polygon, MultiPolygon, or empty
        side_parts = []
        if not result.is_empty and result.area > 0:
            if result.geom_type == "Polygon":
                side_parts.append(result)
            elif result.geom_type == "MultiPolygon":
                side_parts.extend(list(result.geoms))
        return side_parts
    except Exception as e:
        print(f"Warning: Failed to compute side parts: {e}")
        return []


def create_block_parts_from_off_grid(
        block_id: int,
        block: Polygon,
        off_grid: Polygon,
        angle_threshold: float = 155.0,
        corner_distance: float = 1000.0,
) -> tuple[list[Polygon], list[Polygon]]:
    """
    Create corner and side block parts using the off-grid polygon.

    This function identifies sharp corners in the off-grid polygon and creates:
    - Corner parts: Small fan-shaped polygons at each sharp corner
    - Side parts: Strip polygons along the edges between corners
    - Off-grid part: The original off-grid polygon (center part)

    Args:
        block: Original block polygon (outer boundary)
        off_grid: Off-grid polygon (inner part of block)
        angle_threshold: Minimum angle (degrees) to consider a corner "sharp"
        corner_distance: Distance to extend corner fan perpendiculars

    Returns:
        tuple of (corner_parts, side_parts, off_grid_part)
        - corner_parts: list of corner fan polygons
        - side_parts: list of side strip polygons
        - off_grid_part: The original off-grid polygon

    Example:
        >>> block = Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
        >>> off_grid = Polygon([(20, 20), (80, 20), (80, 80), (20, 80)])
        >>> corners, sides, center = create_block_parts_from_off_grid(block, off_grid)
    """
    coords = list(off_grid.exterior.coords)
    n_vertices = len(coords) - 1

    corner_candidates = []
    corner_180_candidates = []
    side_anchor_data = []
    for i in range(n_vertices):
        angle = compute_internal_angle(coords, i)

        if angle >= angle_threshold:
            continue

        dir0 = get_edge_direction(coords, i - 1)
        dir1 = get_edge_direction(coords, i)

        origin = np.array([coords[i][0], coords[i][1], 0.0])

        corner_fan = create_corner_fan(origin, dir0, dir1, corner_distance)
        corner_fan_180 = create_corner_fan(
            origin, dir0, dir1, corner_distance, rotate_180=True,
            rorate_90=False
        )
        corner_candidates.append(corner_fan)
        corner_180_candidates.append(corner_fan_180)

        v0_perp = np.array([dir0[1], -dir0[0], 0.0])
        v1_perp = np.array([dir1[1], -dir1[0], 0.0])
        v0_perp = set_vector_length(v0_perp, corner_distance)
        v1_perp = set_vector_length(v1_perp, corner_distance)
        _v_mid_perp = set_vector_length(v0_perp + v1_perp, corner_distance)

        xyz0 = origin[:2] + v0_perp[:2]
        xyz1 = origin[:2] + v1_perp[:2]

        side_anchor_data.append(
            (
                tuple(origin[:2]),
                tuple(xyz0),
                i - 1,
            )
        )

        side_anchor_data.append(
            (
                tuple(origin[:2]),
                tuple(xyz1),
                i,
            )
        )

    side_candidate_polys = []
    n_anchors = len(side_anchor_data)

    if n_anchors > 0:
        for i in range(1, n_anchors, 2):
            a = side_anchor_data[i % n_anchors]
            b = side_anchor_data[(i + 1) % n_anchors]

            a_origin, a_outer, a_edge_idx = a
            b_origin, b_outer, b_edge_idx = b

            start_idx = (a_edge_idx + 1) % n_vertices
            end_idx = b_edge_idx % n_vertices

            perimeter_positions = get_positions_between_edges(
                coords, start_idx, end_idx, include_end=True
            )
            side_poly = create_side_polygon(
                a_origin, a_outer, b_outer, b_origin, perimeter_positions
            )
            side_candidate_polys.append(side_poly)
    else:
        try:
            frame = block.difference(off_grid)
            if not frame.is_empty and frame.area > 0:
                if frame.geom_type == "Polygon":
                    side_candidate_polys.append(frame)
                elif frame.geom_type == "MultiPolygon":
                    side_candidate_polys.extend(list(frame.geoms))
        except Exception as e:
            print(f"Warning: Failed to create frame: {e}")

    corner_parts = []
    for corner in corner_candidates:
        try:
            intersected = block.intersection(corner)
            if intersected.is_empty or intersected.area == 0:
                continue

            frame_part = intersected.difference(off_grid)
            if not frame_part.is_empty and frame_part.area > 0:
                if frame_part.geom_type == "Polygon":
                    corner_parts.append(frame_part)
                elif frame_part.geom_type == "MultiPolygon":
                    corner_parts.extend(list(frame_part.geoms))
        except Exception as e:
            print(f"Warning: Failed to intersect corner: {e}")
            continue

    side_parts = get_side_parts(block_id, block, off_grid, corner_parts)
    return corner_parts, side_parts, off_grid, corner_candidates, corner_180_candidates


def extract_block_parts_from_off_grid(
        output_path: Path,
        warm_grid_layer_name: str,
        off_grids_inner_layer_name: str,
        off_grid_frame_layer_name: str,
        angle_threshold: float,
        corner_distance: float,
        output_off_grid_layer_name: str,
        output_corner_layer_name: str,
        output_sides_layer_name: str,
        output_not_generated_block_layer_name: str,
):
    """Extract and categorize block parts (corners, sides, off-grids) from blocks with off-grid areas.

    Processes blocks containing off-grid areas to decompose them into distinct
    geometric parts: corner sections, side sections, and the central off-grid area.
    This subdivision is useful for urban planning and block subdivision analysis.

    Args:
        output_path: Path to the GeoPackage file containing input layers and
            where outputs will be saved.
        warm_grid_layer_name: Name of the layer containing the original block
            geometries.
        off_grids_inner_layer_name: Name of the layer containing inner off-grid
            geometries with their associated block_id values.
        off_grid_frame_layer_name: Name of the layer containing frame geometries
            (the perimeter area around off-grids).
        output_off_grid_layer_name: Name of the layer containing frame geometries
            (the perimeter area around off-grids).
        angle_threshold: Maximum angle in degrees to consider a vertex as a corner.
            Vertices with angles below this threshold are identified as corners.
        corner_distance: Distance in meters used to define the depth of corner
            sections extending from corner vertices.
        output_corner_layer_name: Name for the output layer containing corner
            part geometries.
        output_sides_layer_name: Name for the output layer containing side
            part geometries.

    Notes:
        - Each output feature includes metadata: block_id, part_type, area,
          and part_index (for corners and sides).
        - Corner parts are created at sharp vertices of the block boundary.
        - Side parts fill the space between corners along the block edges.
        - Off-grid parts represent the central void within each block.
        - Errors during part creation for individual blocks are caught and logged.
    """
    warm_grid_layer = gpd.read_file(
        output_path, layer=warm_grid_layer_name
    )
    off_grids_inner_layer = gpd.read_file(
        output_path, layer=off_grids_inner_layer_name
    )
    off_grid_frame_layer = gpd.read_file(
        output_path, layer=off_grid_frame_layer_name
    )
    all_corner_parts = []
    all_side_parts = []
    all_off_grid_parts = []
    all_corner_candidates = []
    all_corner_180_candidates = []

    used_block_ids = []
    for idx, frame_row in off_grid_frame_layer.iterrows():
        block_id = frame_row.get("block_id")
        # Find corresponding off-grid
        off_grid_row = off_grids_inner_layer[
            off_grids_inner_layer["block_id"] == block_id].iloc[0]

        # Get original block (without hole)
        original_block = warm_grid_layer.loc[block_id]

        try:
            # Pass original block (not frame) - the function creates parts from block and off_grid
            corner_parts, side_parts, off_grid_final, corner_candidates, corner_180_candidates = create_block_parts_from_off_grid(
                block_id,
                original_block.geometry,  # Original block boundary
                off_grid_row.geometry,  # Off-grid center polygon
                angle_threshold=angle_threshold,
                corner_distance=corner_distance,
            )

            print(
                f"  ✓ Created {len(corner_parts)} corners, {len(side_parts)} sides")

            # Collect corner for candidates with 90 degree rotation
            # Currently we are using this
            for i, corner in enumerate(corner_candidates):
                all_corner_candidates.append(
                    {
                        "geometry": corner,
                        "block_id": block_id,
                        "part_type": "corner",
                        "part_index": i,
                        "area": corner.area,
                    }
                )

            # Collect corner for candidates with 180 degree rotation
            for i, corner in enumerate(corner_180_candidates):
                all_corner_180_candidates.append(
                    {
                        "geometry": corner,
                        "block_id": block_id,
                        "part_type": "corner",
                        "part_index": i,
                        "area": corner.area,
                    }
                )

            # Collect corner parts
            for i, corner in enumerate(corner_parts):
                all_corner_parts.append(
                    {
                        "geometry": corner,
                        "block_id": block_id,
                        "part_type": "corner",
                        "part_index": i,
                        "area": corner.area,
                    }
                )

            # Collect side parts
            for i, side in enumerate(side_parts):
                all_side_parts.append(
                    {
                        "geometry": side,
                        "block_id": block_id,
                        "part_type": "side",
                        "part_index": i,
                        "area": side.area,
                    }
                )

            # Collect off-grid part
            all_off_grid_parts.append(
                {
                    "geometry": off_grid_final,
                    "block_id": block_id,
                    "part_type": "off_grid",
                    "area": off_grid_final.area,
                }
            )
            used_block_ids.append(block_id)
        except Exception as e:
            print(f"  ✗ Error creating parts: {e}")

    if all_corner_candidates:
        gdf_out = gpd.GeoDataFrame(
            all_corner_candidates, crs=warm_grid_layer.crs
        )
        gdf_out.to_file(
            output_path, layer="06_all_corner_candidates", driver="GPKG"
        )

    if all_corner_180_candidates:
        gdf_out = gpd.GeoDataFrame(
            all_corner_180_candidates, crs=warm_grid_layer.crs
        )
        gdf_out.to_file(
            output_path, layer="06_all_corner_180_candidates", driver="GPKG"
        )

    if all_corner_parts:
        gdf_out = gpd.GeoDataFrame(all_corner_parts, crs=warm_grid_layer.crs)
        gdf_out.to_file(
            output_path, layer=output_corner_layer_name, driver="GPKG"
        )
    if all_side_parts:
        gdf_out = gpd.GeoDataFrame(all_side_parts, crs=warm_grid_layer.crs)
        gdf_out.to_file(
            output_path, layer=output_sides_layer_name, driver="GPKG"
        )
    if all_off_grid_parts:
        gdf_out = gpd.GeoDataFrame(all_off_grid_parts, crs=warm_grid_layer.crs)
        gdf_out.to_file(
            output_path, layer=output_off_grid_layer_name, driver="GPKG"
        )

    all_not_used_block = []
    for idx, frame_row in warm_grid_layer.iterrows():
        block_id = idx
        if block_id not in used_block_ids:
            all_not_used_block.append(frame_row)

    if all_not_used_block:
        gdf_out = gpd.GeoDataFrame(all_not_used_block, crs=warm_grid_layer.crs)
        gdf_out.to_file(
            output_path, layer=output_not_generated_block_layer_name,
            driver="GPKG"
        )
