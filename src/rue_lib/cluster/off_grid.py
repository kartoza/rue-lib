# src/rue_lib/cluster/off_grid.py
"""Create off-grid area by offsetting roads inward from block perimeter."""

from typing import Optional

import geopandas as gpd
from shapely.geometry import LineString, Polygon


def get_roads_near_block(
    block: Polygon, roads: gpd.GeoDataFrame, road_type: str, max_distance: float = 10.0
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


def extend_line(line: LineString, extension: float) -> LineString:
    """
    Extend a LineString at both ends.

    Args:
        line: LineString to extend
        extension: Distance to extend at each end

    Returns:
        Extended LineString
    """
    if line.is_empty or line.length == 0:
        return line

    coords = list(line.coords)
    if len(coords) < 2:
        return line

    # Extend start
    start = coords[0]
    second = coords[1]
    dx = start[0] - second[0]
    dy = start[1] - second[1]
    length = (dx**2 + dy**2) ** 0.5
    if length > 0:
        dx /= length
        dy /= length
        new_start = (start[0] + dx * extension, start[1] + dy * extension)
    else:
        new_start = start

    # Extend end
    end = coords[-1]
    second_last = coords[-2]
    dx = end[0] - second_last[0]
    dy = end[1] - second_last[1]
    length = (dx**2 + dy**2) ** 0.5
    if length > 0:
        dx /= length
        dy /= length
        new_end = (end[0] + dx * extension, end[1] + dy * extension)
    else:
        new_end = end

    # Create extended line
    extended_coords = [new_start] + coords[1:-1] + [new_end]
    return LineString(extended_coords)


def create_extended_roads(
    block: Polygon, roads: gpd.GeoDataFrame, extension: float = 100.0, max_distance: float = 10.0
) -> tuple[list[LineString], list[LineString], list[LineString]]:
    """
    Get roads near the block for each road type.

    Args:
        block: Block polygon
        roads: GeoDataFrame with roads
        extension: Distance to extend lines at ends (meters) - not used but kept for compatibility
        max_distance: Maximum distance from block to consider roads

    Returns:
        tuple of (arterial_roads, secondary_roads, local_roads)
    """
    # Get roads near block for each type
    roads_art = get_roads_near_block(block, roads, "road_art", max_distance)
    roads_sec = get_roads_near_block(block, roads, "road_sec", max_distance)
    roads_loc = get_roads_near_block(block, roads, "road_loc", max_distance)

    return roads_art, roads_sec, roads_loc


def clean_small_polygons(polygons: list[Polygon], min_area: float = 1.0) -> list[Polygon]:
    """
    Remove very small polygons that are likely artifacts.

    Args:
        polygons: list of polygons to clean
        min_area: Minimum area to keep

    Returns:
        Cleaned list of polygons
    """
    return [p for p in polygons if p.area >= min_area]


def get_block_edges(block: Polygon) -> list[LineString]:
    """
    Get the four edges of a rectangular block.

    Args:
        block: Block polygon

    Returns:
        list of LineStrings representing edges [bottom, right, top, left]
    """
    coords = list(block.exterior.coords)
    if len(coords) < 5:  # Not a proper closed polygon
        return []

    # For a rectangle, we have 4 edges
    edges = []
    for i in range(len(coords) - 1):
        edge = LineString([coords[i], coords[i + 1]])
        edges.append(edge)

    return edges


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
    closest_type = None

    for _, road in roads.iterrows():
        dist = edge_center.distance(road.geometry)
        if dist < min_distance:
            min_distance = dist
            closest_type = road["road_type"]

    # Return None if the closest road is too far
    if min_distance > max_distance:
        return None

    return closest_type


def create_off_grid_area(
    block: Polygon,
    roads: gpd.GeoDataFrame,
    part_art_d: float = 40.0,
    part_sec_d: float = 30.0,
    part_loc_d: float = 20.0,
    extension: float = 100.0,
) -> Optional[Polygon]:
    """
    Create the off-grid polygon by buffering each edge of the block inward.

    This function works edge-by-edge:
    1. For each edge of the block, find the closest road
    2. Buffer that edge inward by the appropriate distance (art/sec/loc)
    3. Progressively shrink the block by subtracting each buffered edge
    4. Return the final off-grid polygon

    Args:
        block: Block polygon
        roads: GeoDataFrame with roads (must have 'road_type' column)
        part_art_d: Depth of partition along arterial roads
        part_sec_d: Depth of partition along secondary roads
        part_loc_d: Depth of partition along local roads
        extension: Not used, kept for compatibility

    Returns:
        Off-grid polygon or None if no valid area created

    Example:
        >>> block = Polygon([(0, 0), (100, 0), (100, 100), (0, 100)])
        >>> roads = gpd.GeoDataFrame(...)
        >>> off_grid = create_off_grid_area(block, roads)
    """
    # Get edges of the block
    edges = get_block_edges(block)

    if not edges:
        print("Warning: Could not extract block edges")
        return None

    # Start with the original block
    current_polygon = block

    # Process each edge
    for i, edge in enumerate(edges):
        # Find closest road type
        road_type = find_closest_road_type(edge, roads, max_distance=20.0)

        # Determine buffer distance
        if road_type == "road_art":
            buffer_dist = part_art_d
            road_label = "arterial"
        elif road_type == "road_sec":
            buffer_dist = part_sec_d
            road_label = "secondary"
        elif road_type == "road_loc":
            buffer_dist = part_loc_d
            road_label = "local"
        else:
            # No road nearby, use default (local)
            buffer_dist = part_loc_d
            road_label = "default(local)"

        # Create a buffer strip from this edge
        try:
            # Buffer the edge line
            edge_buffer = edge.buffer(buffer_dist, cap_style=2)  # flat caps

            # Subtract from current polygon
            new_polygon = current_polygon.difference(edge_buffer)

            # Validate result
            if new_polygon.is_empty:
                print(f"  Edge {i} ({road_label}): Buffer consumed entire polygon")
                return None

            # Handle MultiPolygon - take largest
            if new_polygon.geom_type == "MultiPolygon":
                polygons = list(new_polygon.geoms)
                polygons = clean_small_polygons(polygons, min_area=1.0)
                if not polygons:
                    return None
                new_polygon = max(polygons, key=lambda p: p.area)
            elif new_polygon.geom_type != "Polygon":
                print(f"  Edge {i} ({road_label}): Invalid geometry type {new_polygon.geom_type}")
                return None

            current_polygon = new_polygon

        except Exception as e:
            print(f"  Warning: Failed to process edge {i} ({road_label}): {e}")
            continue

    # Final validation
    if current_polygon.area < 1.0:
        return None

    return current_polygon
