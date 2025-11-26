# src/rue_lib/cluster/block_parts.py
"""Create corner and side block parts from off-grid polygons."""

import numpy as np
from shapely.geometry import Polygon


def vector_length(v: np.ndarray) -> float:
    """Calculate length of a vector."""
    return np.linalg.norm(v)


def normalize_vector(v: np.ndarray) -> np.ndarray:
    """Normalize a vector to unit length."""
    length = vector_length(v)
    if length == 0:
        return v
    return v / length


def set_vector_length(v: np.ndarray, length: float) -> np.ndarray:
    """Set vector to a specific length."""
    return normalize_vector(v) * length


def compute_internal_angle(coords: list[tuple[float, float]], idx: int) -> float:
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


def get_edge_direction(coords: list[tuple[float, float]], idx: int) -> np.ndarray:
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
    origin: np.ndarray, dir0: np.ndarray, dir1: np.ndarray, distance: float = 1000.0
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
    # Compute perpendicular vectors (rotate 90° counterclockwise)
    v0_perp = -np.array([dir0[1], -dir0[0], 0.0])
    v1_perp = -np.array([dir1[1], -dir1[0], 0.0])

    v0_perp = set_vector_length(v0_perp, distance)
    v1_perp = set_vector_length(v1_perp, distance)

    v_mid_perp = v0_perp + v1_perp
    v_mid_perp = set_vector_length(v_mid_perp, distance)

    xyz0 = origin[:2] + v0_perp[:2]
    xyz1 = origin[:2] + v1_perp[:2]
    xyz_mid = origin[:2] + v_mid_perp[:2]

    fan_coords = [tuple(origin[:2]), tuple(xyz0), tuple(xyz_mid), tuple(xyz1), tuple(origin[:2])]

    return Polygon(fan_coords)


def get_positions_between_edges(
    coords: list[tuple[float, float]], start_idx: int, end_idx: int, include_end: bool = True
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


def create_block_parts_from_off_grid(
    block: Polygon,
    off_grid: Polygon,
    angle_threshold: float = 155.0,
    corner_distance: float = 1000.0,
) -> tuple[list[Polygon], list[Polygon], Polygon]:
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
    side_anchor_data = []

    for i in range(n_vertices):
        angle = compute_internal_angle(coords, i)

        if angle >= angle_threshold:
            continue

        dir0 = get_edge_direction(coords, i - 1)
        dir1 = get_edge_direction(coords, i)

        origin = np.array([coords[i][0], coords[i][1], 0.0])

        corner_fan = create_corner_fan(origin, dir0, dir1, corner_distance)
        corner_candidates.append(corner_fan)

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

    side_parts = []
    for side in side_candidate_polys:
        try:
            intersected = block.intersection(side)
            if intersected.is_empty or intersected.area == 0:
                continue
            frame_part = intersected.difference(off_grid)
            if not frame_part.is_empty and frame_part.area > 0:
                if frame_part.geom_type == "Polygon":
                    side_parts.append(frame_part)
                elif frame_part.geom_type == "MultiPolygon":
                    side_parts.extend(list(frame_part.geoms))
        except Exception as e:
            print(f"Warning: Failed to intersect side: {e}")
            continue

    return corner_parts, side_parts, off_grid
