# src/rue_lib/cluster/off_grid_subdivision.py
"""Subdivide off-grid areas into smaller plot clusters."""

from typing import Optional

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Polygon

from rue_lib.core.definitions import RoadTypes


def compute_angle_dot(polygon: Polygon, vertex_idx: int) -> float:
    """
    Compute the dot product of normalized edge vectors at a vertex.

    This helps identify sharp corners (low dot product = sharp angle).

    Args:
        polygon: Input polygon
        vertex_idx: Index of vertex to check

    Returns:
        Dot product of normalized edge vectors (-1 to 1)
    """
    coords = list(polygon.exterior.coords)[:-1]
    n = len(coords)

    p_prev = np.array(coords[(vertex_idx - 1) % n])
    p_curr = np.array(coords[vertex_idx % n])
    p_next = np.array(coords[(vertex_idx + 1) % n])

    vec0 = p_curr - p_prev
    vec1 = p_next - p_curr

    vec0_norm = vec0 / np.linalg.norm(vec0) if np.linalg.norm(
        vec0) > 0 else vec0
    vec1_norm = vec1 / np.linalg.norm(vec1) if np.linalg.norm(
        vec1) > 0 else vec1

    return np.dot(vec0_norm, vec1_norm)


def convert_to_quadrilateral(
        polygon: Polygon, min_area: float = 100.0
) -> Optional[Polygon]:
    """
    Convert an off-grid polygon to a quadrilateral by removing vertices with sharpest angles.

    This simplifies complex polygons to rectangles/quads suitable for grid subdivision.

    Args:
        polygon: Input off-grid polygon
        min_area: Minimum area threshold (m²)

    Returns:
        Quadrilateral polygon or None if area too small or too few vertices
    """
    if polygon.area < min_area:
        return None

    coords = list(polygon.exterior.coords)[:-1]
    n_vertices = len(coords)

    if n_vertices < 4:
        return None

    # Compute angle dot products for all vertices
    vertex_dots = []
    for i in range(n_vertices):
        dot = compute_angle_dot(polygon, i)
        vertex_dots.append((i, dot))

    sorted_vertices = sorted(vertex_dots, key=lambda x: x[1], reverse=True)

    vertices_to_keep = sorted([v[0] for v in sorted_vertices[-4:]])

    quad_coords = [coords[i] for i in vertices_to_keep]
    quad_coords.append(quad_coords[0])
    quad = Polygon(quad_coords)

    worst_dot = sorted_vertices[-1][1]

    if worst_dot > -0.7:
        return quad

    return quad


def get_roads_near_block(
        block: Polygon, roads: gpd.GeoDataFrame, road_type: str,
        max_distance: float = 10.0
) -> list[LineString]:
    """
    Get roads of a specific type that are near (within max_distance) of the block.

    Args:
        block: Block polygon
        roads: GeoDataFrame with roads (must have 'road_type' column)
        road_type: Type of road to filter ('road_art', 'road_sec', 'road_loc')
        max_distance: Maximum distance from block to consider

    Returns:
        list of LineStrings near the block
    """
    # Filter roads by type
    roads_filtered = (
        roads[roads["road_type"] == road_type]
        if "road_type" in roads.columns
        else gpd.GeoDataFrame()
    )

    if roads_filtered.empty:
        return []

    # Find roads that are near the block
    nearby_roads = []
    block_buffered = block.buffer(max_distance)

    for _, road in roads_filtered.iterrows():
        road_geom = road.geometry

        # Check if road intersects the buffered block
        if road_geom.intersects(block_buffered):
            nearby_roads.append(road_geom)

    return nearby_roads


def find_closest_road_type(
        edge: LineString, roads: gpd.GeoDataFrame, max_distance: float = 20.0
) -> Optional[str]:
    """
    Find the road type that is closest to the center point of a given edge.

    Args:
        edge: Edge line of the block
        roads: GeoDataFrame with roads
        max_distance: Maximum distance to consider

    Returns:
        Road type string or None if no road within max_distance
    """
    if roads.empty or "road_type" not in roads.columns:
        return None

    # Use the center point of the edge
    edge_center = edge.interpolate(0.5, normalized=True)

    min_distance = float("inf")
    closest_type = RoadTypes.Local

    for _, road in roads.iterrows():
        dist = edge_center.distance(road.geometry)
        if dist < min_distance:
            min_distance = dist
            closest_type = road["road_type"]

    # Return None if the closest road is too far
    if min_distance > max_distance:
        return None

    if closest_type == "road_art":
        closest_type = RoadTypes.Artery
    elif closest_type == "road_sec":
        closest_type = RoadTypes.Secondary
    return closest_type
