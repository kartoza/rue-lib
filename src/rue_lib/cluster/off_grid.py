# src/rue_lib/cluster/off_grid.py
"""Create off-grid area by offsetting roads inward from block perimeter."""

from pathlib import Path
from typing import Optional

import geopandas as gpd
from shapely.geometry import LineString, Polygon

from rue_lib.core.geometry import clean_angles


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
    length = (dx ** 2 + dy ** 2) ** 0.5
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
    length = (dx ** 2 + dy ** 2) ** 0.5
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
        block: Polygon, roads: gpd.GeoDataFrame, extension: float = 100.0,
        max_distance: float = 10.0
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


def clean_small_polygons(polygons: list[Polygon], min_area: float = 1.0) -> \
        list[Polygon]:
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


def create_off_grid_inner_layer(
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
    edges = get_block_edges(block)

    if not edges:
        print("Warning: Could not extract block edges")
        return None

    current_polygon = block

    for i, edge in enumerate(edges):
        road_type = find_closest_road_type(edge, roads, max_distance=20.0)

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
            buffer_dist = part_loc_d
            road_label = "default(local)"

        try:
            edge_buffer = edge.buffer(buffer_dist, cap_style=2)

            new_polygon = current_polygon.difference(edge_buffer)

            if new_polygon.is_empty:
                print(
                    f"  Edge {i} ({road_label}): Buffer consumed entire polygon")
                return None

            if new_polygon.geom_type == "MultiPolygon":
                polygons = list(new_polygon.geoms)
                polygons = clean_small_polygons(polygons, min_area=1.0)
                if not polygons:
                    return None
                new_polygon = max(polygons, key=lambda p: p.area)
            elif new_polygon.geom_type != "Polygon":
                print(
                    f"  Edge {i} ({road_label}): Invalid geometry type {new_polygon.geom_type}")
                return None

            current_polygon = new_polygon

        except Exception as e:
            print(f"  Warning: Failed to process edge {i} ({road_label}): {e}")
            continue

    # Final validation
    if current_polygon.area < 1.0:
        return None

    return clean_angles(current_polygon, 0.01)


def create_off_grid_inner_layers(
        blocks: gpd.GeoDataFrame,
        roads: gpd.GeoDataFrame,
        part_art_d: float = 40.0,
        part_sec_d: float = 30.0,
        part_loc_d: float = 20.0,
        extension: float = 100.0,
) -> list[dict]:
    """
    Create off-grid areas for multiple blocks by buffering edges inward based on adjacent road types.

    This function processes each block in the input GeoDataFrame and creates an off-grid area
    by progressively buffering each edge of the block inward. The buffer distance depends on
    the type of road adjacent to each edge (arterial, secondary, or local).

    :param blocks: GeoDataFrame containing block polygons to process. Must have a geometry column.
    :type blocks: gpd.GeoDataFrame
    :param roads: GeoDataFrame containing road linestrings with road type information.
                  Must have 'road_type' column with values like 'road_art', 'road_sec', 'road_loc'.
    :type roads: gpd.GeoDataFrame
    :param part_art_d: Buffer distance (in meters) for edges adjacent to arterial roads.
    :type part_art_d: float
    :param part_sec_d: Buffer distance (in meters) for edges adjacent to secondary roads.
    :type part_sec_d: float
    :param part_loc_d: Buffer distance (in meters) for edges adjacent to local roads.
    :type part_loc_d: float
    :param extension: Extension distance for roads (not currently used, kept for compatibility).
    :type extension: float
    :return: List of dictionaries, one per input block, each containing:
             - geometry: The resulting off-grid polygon (or original block if no valid off-grid created)
             - block_id: Identifier of the block (from 'id' column or row index)
             - original_area: Area of the original block
             - off_grid_area: Area of the off-grid polygon (0 if no valid off-grid created)
             - reduction_pct: Percentage reduction from original block area
    :rtype: list[dict]

    :Example:

    >>> import geopandas as gpd
    >>> from shapely.geometry import Polygon, LineString
    >>> blocks = gpd.GeoDataFrame({
    ...     'id': [1, 2],
    ...     'geometry': [Polygon([(0, 0), (100, 0), (100, 100), (0, 100)]),
    ...                  Polygon([(200, 0), (300, 0), (300, 100), (200, 100)])]
    ... })
    >>> roads = gpd.GeoDataFrame({
    ...     'road_type': ['road_art', 'road_loc'],
    ...     'geometry': [LineString([(0, 0), (100, 0)]), LineString([(200, 0), (300, 0)])]
    ... })
    >>> result = create_off_grid_areas(blocks, roads)
    >>> len(result)
    2
    """
    # Test each block
    off_grids = []
    for idx, block_row in blocks.iterrows():
        block = block_row.geometry
        block_id = idx

        off_grid = create_off_grid_inner_layer(
            block,
            roads,
            part_art_d=part_art_d,
            part_sec_d=part_sec_d,
            part_loc_d=part_loc_d,
        )

        if off_grid:
            reduction = ((block.area - off_grid.area) / block.area) * 100
            off_grids.append(
                {
                    "geometry": off_grid,
                    "block_id": block_id,
                    "original_area": block.area,
                    "off_grid_area": off_grid.area,
                    "reduction_pct": reduction,
                }
            )
        else:
            off_grids.append(
                {
                    "geometry": block,  # Keep original
                    "block_id": block_id,
                    "original_area": block.area,
                    "off_grid_area": 0,
                    "reduction_pct": 0,
                }
            )
    return off_grids


def extract_off_grid_inner_layer(
        output_path: Path,
        roads_layer_name: str,
        off_grid_layer_name: str,
        output_layer_name: str,
        part_art_d: float = 40.0,
        part_sec_d: float = 30.0,
        part_loc_d: float = 20.0,
) -> str:
    """
    Extract off-grid areas from a GeoPackage and save the results to a new layer.

    This function reads off-grid blocks and roads from a GeoPackage file, processes them
    to create off-grid areas by buffering edges based on adjacent road types, and saves
    the results back to the GeoPackage in a new layer.

    :param output_path: Path to the GeoPackage file containing input layers and where output will be saved
    :type output_path: Path
    :param roads_layer_name: Name of the layer in the GeoPackage containing road data
    :type roads_layer_name: str
    :param off_grid_layer_name: Name of the layer in the GeoPackage containing off-grid blocks
    :type off_grid_layer_name: str
    :param output_layer_name: Name of the output layer to create in the GeoPackage
    :type output_layer_name: str
    :param part_art_d: Buffer distance (in meters) for edges adjacent to arterial roads
    :type part_art_d: float
    :param part_sec_d: Buffer distance (in meters) for edges adjacent to secondary roads
    :type part_sec_d: float
    :param part_loc_d: Buffer distance (in meters) for edges adjacent to local roads
    :type part_loc_d: float
    :return: Name of the output layer that was created
    :rtype: str

    :Example:

    >>> from pathlib import Path
    >>> output_layer = extract_off_grid_inner_layer(
    ...     output_path=Path("outputs/results.gpkg"),
    ...     roads_layer_name="roads",
    ...     off_grid_layer_name="off_grid_blocks",
    ...     output_layer_name="off_grid_areas",
    ...     part_art_d=40.0,
    ...     part_sec_d=30.0,
    ...     part_loc_d=20.0
    ... )
    >>> print(f"Created layer: {output_layer}")
    """
    off_grid_layer = gpd.read_file(
        output_path, layer=off_grid_layer_name
    )
    roads_layer = gpd.read_file(
        output_path, layer=roads_layer_name
    )
    off_grid_cluster_layer = create_off_grid_inner_layers(
        off_grid_layer, roads_layer, part_art_d, part_sec_d, part_loc_d
    )

    gdf_out = gpd.GeoDataFrame(off_grid_cluster_layer, crs=off_grid_layer.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")

    print(
        f"  Created {len(off_grid_cluster_layer)} "
        f"perpendicular lines from guide points"
    )
    print(f"  Saved to layer: {output_layer_name}")

    return output_layer_name
