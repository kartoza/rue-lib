# src/rue_lib/cluster/subdiv_at_concave_corner.py
"""
Subdivide blocks at concave corners.

This module identifies concave corners (angles > 250 degrees) in blocks and
subdivides them by creating cut lines perpendicular to the block edges at those corners.
"""

from pathlib import Path
from typing import Optional

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Polygon

from rue_lib.core.definitions import RoadTypes
from rue_lib.core.geometry import remove_vertices_by_angle


def get_vertex_angle(coords: list, vertex_idx: int) -> float:
    """
    Calculate the angle at a vertex in a polygon.

    Uses the angle between the incoming and outgoing edge vectors.

    Args:
        coords: List of polygon coordinates
        vertex_idx: Index of vertex to calculate angle at

    Returns:
        Angle in degrees (0-360)
    """
    n = len(coords) - 1  # Exclude duplicate last point

    # Get previous, current, and next points
    p_prev = np.array(coords[(vertex_idx - 1) % n])
    p_curr = np.array(coords[vertex_idx % n])
    p_next = np.array(coords[(vertex_idx + 1) % n])

    # Get edge vectors
    vec0 = p_prev - p_curr  # Incoming edge
    vec1 = p_next - p_curr  # Outgoing edge

    # Normalize vectors
    vec0_norm = vec0 / (np.linalg.norm(vec0) + 1e-10)
    vec1_norm = vec1 / (np.linalg.norm(vec1) + 1e-10)

    # Reverse vec0 to point away from vertex
    vec0_rev = -vec0_norm

    # Calculate angle using dot product
    dot = np.clip(np.dot(vec1_norm, vec0_rev), -1.0, 1.0)
    angle_rad = np.arccos(dot)

    # Use cross product to determine sign
    cross = np.cross(
        np.append(vec0_rev, 0),
        np.append(vec1_norm, 0)
    )

    # Convert to degrees
    angle_deg = np.degrees(angle_rad)

    # Adjust angle based on cross product sign
    if cross[2] < 0:
        angle_deg = 360 - angle_deg

    return angle_deg


def get_edges_on_roads(
        block: Polygon,
        roads: gpd.GeoDataFrame
) -> tuple[list[int], list[int]]:
    """
    Separate block edges into road edges and other edges.

    Groups consecutive edges that are on roads (have road_type attribute)
    and separates them from edges that are not on roads.

    Args:
        block: Block polygon

    Returns:
        Tuple of (road_edge_indices, other_edge_indices)
    """
    road_edge_lists = [[]]
    other_edges = []

    coords = list(block.exterior.coords)[:-1]
    n_edges = len(coords)

    for i in range(n_edges):
        # TODO
        #  Fix this
        edge_road = RoadTypes.Artery
        if edge_road in [RoadTypes.Artery, RoadTypes.Local, RoadTypes.Local]:
            road_edge_lists[-1].append(i)
        else:
            other_edges.append(i)
            if len(road_edge_lists[-1]) != 0:
                road_edge_lists.append([])

    # Flatten road edges
    if len(road_edge_lists) == 0:
        road_edges = []
    elif len(road_edge_lists) == 1:
        road_edges = road_edge_lists[0]
    else:
        # Join last and first groups if they're consecutive
        road_edges = road_edge_lists[1:][0] + road_edge_lists[0]
        for group in road_edge_lists[2:]:
            road_edges.extend(group)

    return road_edges, other_edges


def ray_intersect_edges(
        ray_origin: np.ndarray,
        ray_direction: np.ndarray,
        edges: list[LineString],
        exclude_edges: list[LineString] = None
) -> tuple[Optional[np.ndarray], Optional[int]]:
    """
    Find the nearest intersection between a ray and a list of edges.

    Args:
        ray_origin: Ray origin point [x, y]
        ray_direction: Ray direction vector [x, y]
        edges: List of LineString edges to test
        exclude_edges: List of edges to exclude from intersection test

    Returns:
        Tuple of (intersection_point, edge_index) or (None, None) if no intersection
    """
    if exclude_edges is None:
        exclude_edges = []

    # Check if ray intersects any excluded edges
    for _exclude_edge in exclude_edges:
        # Simple check - in a full implementation would do proper ray-line intersection
        pass

    isect_min = None
    dist_min = np.inf
    edge_min_idx = None

    # Extend ray to a very long line
    ray_origin + ray_direction * 10000

    for edge_idx, edge in enumerate(edges):
        # Get edge coordinates
        edge_coords = list(edge.coords)
        if len(edge_coords) < 2:
            continue

        p1 = np.array(edge_coords[0])
        p2 = np.array(edge_coords[1])

        # Calculate intersection using parametric line equations
        # Ray: P = ray_origin + t * ray_direction (t >= 0)
        # Edge: Q = p1 + s * (p2 - p1) (0 <= s <= 1)

        edge_vec = p2 - p1
        edge_len = np.linalg.norm(edge_vec)

        if edge_len < 1e-10:
            continue

        # Solve for intersection
        # ray_origin + t * ray_direction = p1 + s * edge_vec
        # This is a 2x2 linear system

        denom = ray_direction[0] * edge_vec[1] - ray_direction[1] * edge_vec[0]

        if abs(denom) < 1e-10:
            # Lines are parallel
            continue

        diff = p1 - ray_origin
        t = (diff[0] * edge_vec[1] - diff[1] * edge_vec[0]) / denom
        s = (diff[0] * ray_direction[1] - diff[1] * ray_direction[0]) / denom

        # Check if intersection is valid
        if t >= 0 and 0 <= s <= 1:
            isect = ray_origin + t * ray_direction
            dist = np.linalg.norm(isect - ray_origin)

            if dist < dist_min:
                isect_min = isect
                dist_min = dist
                edge_min_idx = edge_idx

    return isect_min, edge_min_idx


def get_positions_from_ring(
        positions: list,
        idx_start: int,
        idx_end: int
) -> list:
    """
    Get positions along a ring between two indices.

    Args:
        positions: List of positions in the ring
        idx_start: Starting index (exclusive)
        idx_end: Ending index (exclusive)

    Returns:
        List of positions between the indices
    """
    num_posis = len(positions)

    idx0 = (idx_start + 1) % num_posis
    idx1 = (idx_end + 1) % num_posis

    if idx0 == idx1:
        return []

    ring = []

    if idx0 < idx1:
        ring.extend(positions[idx0:idx1])
    else:
        ring.extend(positions[idx0:])
        ring.extend(positions[:idx1])

    return ring


def make_subdivision(
        block_coords: list,
        cut0: Optional[tuple[int, int]],
        cut1: Optional[tuple[int, int]],
        posis0: Optional[list],
        posis1: Optional[list]
) -> Optional[Polygon]:
    """
    Create a subdivision polygon from block coordinates and cut positions.

    Args:
        block_coords: Original block coordinates
        cut0: First cut edge indices (start_edge_idx, end_edge_idx)
        cut1: Second cut edge indices (start_edge_idx, end_edge_idx)
        posis0: Positions for first cut line
        posis1: Positions for second cut line

    Returns:
        Subdivision polygon or None if invalid
    """
    ring = []

    if cut0 is None and cut1 is None:
        # Return the whole block
        return Polygon(block_coords)

    elif cut0 is None:
        # Only cut1 exists
        ring = list(posis1)

        # Get positions from ring
        posis_list = get_positions_from_ring(
            block_coords[:-1],  # Exclude duplicate last point
            cut1[1],
            cut1[0]
        )
        ring.extend(posis_list)

    elif cut1 is None:
        # Only cut0 exists
        ring = list(reversed(posis0))

        # Get positions from ring
        posis_list = get_positions_from_ring(
            block_coords[:-1],
            cut0[0],
            cut0[1]
        )
        ring.extend(posis_list)

    else:
        # Both cuts exist
        ring = list(posis1)

        # Get positions between cuts
        posis_list = get_positions_from_ring(
            block_coords[:-1],
            cut1[1],
            cut0[1]
        )
        ring.extend(posis_list)

        # Add reversed posis0
        ring.extend(reversed(posis0))

        # Get remaining positions
        posis_list = get_positions_from_ring(
            block_coords[:-1],
            cut0[0],
            cut1[0]
        )
        ring.extend(posis_list)

    if len(ring) < 3:
        return None

    # Close the ring
    ring.append(ring[0])

    return Polygon(ring)


def subdivide_block_at_concave_corners(
        block: Polygon,
        roads: gpd.GeoDataFrame,
        plot_loc_w: float = 20.0,
        part_loc_d: float = 20.0,
        angle_threshold: float = 250.0,
        max_cut_distance: float = 300.0
) -> list[Polygon]:
    """
    Subdivide a block at concave corners.

    This function identifies concave corners (angles > threshold) in a block
    and creates subdivisions by cutting perpendicular to the edges at those corners.

    Args:
        block: Block polygon to subdivide
        plot_loc_w: Plot width on local roads (meters)
        part_loc_d: Partition depth on local roads (meters)
        angle_threshold: Minimum angle to consider a corner concave (degrees)
        max_cut_distance: Maximum distance for cut lines (meters)

    Returns:
        List of subdivision polygons

    Example:
        >>> block = Polygon([(0, 0), (100, 0), (100, 100), (50, 50), (0, 100)])
        >>> edge_types = ["road_loc", "road_loc", None, None, "road_loc"]
        >>> subdivs = subdivide_block_at_concave_corners(block)
    """
    coords = list(block.exterior.coords)
    block_coords = coords[:-1]  # Exclude duplicate last point

    # Get edges on roads
    road_edges, other_edges = get_edges_on_roads(block, roads=roads)
    if len(road_edges) == 0:
        return []

    # Calculate edge lengths
    edge_lengths = []
    for i in road_edges:
        p0 = np.array(block_coords[i])
        p1 = np.array(block_coords[(i + 1) % len(block_coords)])
        edge_lengths.append(np.linalg.norm(p1 - p0))

    # Find longest road edge
    longest_idx = np.argmax(edge_lengths)
    edge_long_idx = road_edges[longest_idx]

    # Get longest edge vector
    p0 = np.array(block_coords[edge_long_idx])
    p1 = np.array(block_coords[(edge_long_idx + 1) % len(block_coords)])
    vec_long = p1 - p0
    vec_long = vec_long / np.linalg.norm(vec_long)

    # Perpendicular vector
    np.array([-vec_long[1], vec_long[0]])

    # Create edges as LineStrings
    other_edge_lines = []
    for i in other_edges:
        p0 = block_coords[i]
        p1 = block_coords[(i + 1) % len(block_coords)]
        other_edge_lines.append(LineString([p0, p1]))

    subdivisions = []
    cut0 = None
    cut1 = None
    posis0 = None
    posis1 = None

    # Iterate through road edges looking for concave corners
    for i in range(1, len(road_edges)):
        cut0 = cut1
        posis0 = posis1

        # Get vertex at this road edge
        vert_idx = road_edges[i]
        xyz = np.array(block_coords[vert_idx])

        # Get edges at this vertex
        edge_idx0 = (vert_idx - 1) % len(block_coords)
        edge_idx1 = vert_idx

        # Calculate vertex angle
        angle = get_vertex_angle(coords, vert_idx)

        # Check if concave corner (angle > threshold)
        if angle <= angle_threshold:
            continue

        # Get edge vectors at vertex
        p_prev = np.array(block_coords[(vert_idx - 1) % len(block_coords)])
        p_next = np.array(block_coords[(vert_idx + 1) % len(block_coords)])

        vec0 = xyz - p_prev
        vec1 = p_next - xyz

        vec0 = vec0 / (np.linalg.norm(vec0) + 1e-10)
        vec1 = vec1 / (np.linalg.norm(vec1) + 1e-10)

        # Perpendicular vectors (rotate 90 degrees)
        perp_vec0 = np.array([vec0[1], -vec0[0]])
        perp_vec1 = np.array([vec1[1], -vec1[0]])

        # Try cutting from both sides
        xyz_a = None
        isect_a = None
        edge_a_idx = None

        if edge_lengths[i - 1] > part_loc_d:
            # Create ray from edge 0
            -perp_vec0 * 1000
            xyz_a = xyz - vec0 * (plot_loc_w / 2)

            isect_a, edge_a_idx = ray_intersect_edges(
                xyz_a, -perp_vec0, other_edge_lines
            )

            # Check distance
            if isect_a is not None and np.linalg.norm(
                    isect_a - xyz_a) > max_cut_distance:
                isect_a = None
                edge_a_idx = None

        xyz_b = None
        isect_b = None
        edge_b_idx = None

        if i < len(edge_lengths) and edge_lengths[i] > part_loc_d:
            # Create ray from edge 1
            -perp_vec1 * 1000
            xyz_b = xyz + vec1 * (plot_loc_w / 2)

            isect_b, edge_b_idx = ray_intersect_edges(
                xyz_b, -perp_vec1, other_edge_lines
            )

            # Check distance
            if isect_b is not None and np.linalg.norm(
                    isect_b - xyz_b) > max_cut_distance:
                isect_b = None
                edge_b_idx = None

        # Adjust cuts if only one side worked
        if isect_a is not None and isect_b is None:
            xyz_a = xyz
            isect_a, edge_a_idx = ray_intersect_edges(
                xyz_a, -perp_vec0, other_edge_lines
            )

        if isect_a is None and isect_b is not None:
            xyz_b = xyz
            isect_b, edge_b_idx = ray_intersect_edges(
                xyz_b, -perp_vec1, other_edge_lines
            )

        # Create subdivisions based on intersections
        if isect_a is None and isect_b is None:
            continue

        if isect_a is not None and isect_b is None:
            # Only cut from edge 0
            posis1 = [tuple(xyz_a), tuple(isect_a)]
            idx_a = other_edges[edge_a_idx] if edge_a_idx is not None else 0
            cut1 = (edge_idx0, idx_a)

            subdiv = make_subdivision(coords, cut0, cut1, posis0, posis1)
            if subdiv is not None:
                subdivisions.append(subdiv)

        elif isect_a is None and isect_b is not None:
            # Only cut from edge 1
            posis1 = [tuple(xyz_b), tuple(isect_b)]
            idx_b = other_edges[edge_b_idx] if edge_b_idx is not None else 0
            cut1 = (edge_idx1, idx_b)

            subdiv = make_subdivision(coords, cut0, cut1, posis0, posis1)
            if subdiv is not None:
                subdivisions.append(subdiv)

        else:
            # Both cuts work - create corner subdivision
            # First subdivision
            posis1 = [tuple(xyz_a), tuple(isect_a)]
            idx_a = other_edges[edge_a_idx] if edge_a_idx is not None else 0
            cut1 = (edge_idx0, idx_a)

            subdiv = make_subdivision(coords, cut0, cut1, posis0, posis1)
            if subdiv is not None:
                subdivisions.append(subdiv)

            # Update for corner
            posis0 = posis1
            cut0 = cut1

            # Second subdivision (corner)
            posis1 = [tuple(xyz_b), tuple(isect_b)]
            idx_b = other_edges[edge_b_idx] if edge_b_idx is not None else 0
            cut1 = (edge_idx1, idx_b)

            subdiv = make_subdivision(coords, cut0, cut1, posis0, posis1)
            if subdiv is not None:
                subdivisions.append(subdiv)

    # Create final subdivision
    posis0 = posis1
    cut0 = cut1
    cut1 = None

    subdiv = make_subdivision(coords, cut0, cut1, posis0, posis1)
    if subdiv is not None:
        subdivisions.append(subdiv)

    return subdivisions


def subdivide_blocks_at_concave_corners(
        output_path: Path,
        input_layer_name: str,
        roads_layer_name: str,
        output_layer_name: str,
        plot_loc_w: float = 20.0,
        part_loc_d: float = 20.0,
        angle_threshold: float = 250.0,
        max_cut_distance: float = 300.0
) -> str:
    """
    Subdivide blocks at concave corners and save to GeoPackage.

    Reads blocks from a GeoPackage layer, identifies concave corners, creates
    subdivisions, and saves the results to a new layer.

    Args:
        output_path: Path to GeoPackage file
        input_layer_name: Name of layer containing input blocks
        output_layer_name: Name of layer for output subdivisions
        plot_loc_w: Plot width on local roads (meters)
        part_loc_d: Partition depth on local roads (meters)
        angle_threshold: Minimum angle to consider a corner concave (degrees)
        max_cut_distance: Maximum distance for cut lines (meters)

    Returns:
        Name of the output layer created
    """
    # Read input blocks
    blocks_layer = gpd.read_file(output_path, layer=input_layer_name)
    road_layer = gpd.read_file(output_path, layer=roads_layer_name)

    all_subdivisions = []
    count = 0

    for idx, block_row in blocks_layer.iterrows():
        # Convert PolygonZ to Polygon (remove Z dimension)
        geom = block_row.geometry
        if geom.has_z:
            coords_2d = [(x, y) for x, y, *_ in geom.exterior.coords]
            geom = Polygon(coords_2d)

        block = remove_vertices_by_angle(
            geom, min_angle_threshold=5
        )
        block_id = idx

        # Subdivide block
        subdivs = subdivide_block_at_concave_corners(
            block,
            roads=road_layer,
            plot_loc_w=plot_loc_w,
            part_loc_d=part_loc_d
        )

        # Add subdivisions to results
        for i, subdiv in enumerate(subdivs):
            all_subdivisions.append({
                'geometry': subdiv,
                'block_id': f"{block_id}_{i}",
                'original_block_id': block_id,
                'subdivision_index': i,
                'area': subdiv.area,
                'type': 'block_corner' if i > 0 and i < len(
                    subdivs) - 1 else 'block'
            })
            count += 1

        print(f"  Block {block_id}: Created {len(subdivs)} subdivisions")

    # Save results
    if all_subdivisions:
        gdf_out = gpd.GeoDataFrame(all_subdivisions, crs=blocks_layer.crs)
        gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")
        print(f"  ✓ Saved to layer: {output_layer_name}")
    else:
        print("  ✗ No subdivisions created")

    return output_layer_name
