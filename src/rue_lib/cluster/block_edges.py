# src/rue_lib/cluster/block_edges.py
"""Extract edges from block polygons and assign road types."""

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Point, Polygon


def extract_block_edges(
    blocks_gdf: gpd.GeoDataFrame,
    roads_gdf: gpd.GeoDataFrame,
    tolerance: float = 5.0
) -> gpd.GeoDataFrame:
    """
    Extract edges (boundary segments) from blocks and assign road types.

    This function decomposes block polygons into individual edge segments
    and determines which road type each edge is adjacent to based on
    proximity to roads.

    Args:
        blocks_gdf: GeoDataFrame with block polygons (must have 'block_id' column)
        roads_gdf: GeoDataFrame with roads (must have 'road_type' column)
        tolerance: Distance tolerance for matching edges to roads (default: 5m)

    Returns:
        GeoDataFrame with edges as LineStrings, containing:
            - geometry: LineString for each edge
            - block_id: ID of the parent block
            - road_type: Type of adjacent road ('road_art', 'road_sec', 'road_loc', or None)
            - edge_index: Index of edge within the block (0-based)

    Example:
        >>> blocks = gpd.read_file('blocks.geojson')
        >>> roads = gpd.read_file('roads.geojson')
        >>> edges = extract_block_edges(blocks, roads, tolerance=5.0)
    """
    all_edges = []

    for idx, block_row in blocks_gdf.iterrows():
        block = block_row.geometry
        block_id = block_row.get('block_id', idx)

        if not isinstance(block, Polygon):
            continue

        # Get coordinates of the exterior ring
        coords = list(block.exterior.coords)

        # Create edge segments (skip the last duplicate point)
        for i in range(len(coords) - 1):
            edge_geom = LineString([coords[i], coords[i + 1]])

            # Find the nearest road to this edge
            road_type = _find_nearest_road_type(
                edge_geom, roads_gdf, tolerance
            )

            edge_data = {
                'geometry': edge_geom,
                'block_id': block_id,
                'road_type': road_type,
                'edge_index': i
            }

            # Copy over other block attributes
            for col in blocks_gdf.columns:
                if col not in ['geometry', 'block_id', 'edge_index', 'road_type']:
                    edge_data[col] = block_row.get(col)

            all_edges.append(edge_data)

    # Create GeoDataFrame
    if all_edges:
        edges_gdf = gpd.GeoDataFrame(all_edges, crs=blocks_gdf.crs)
        return edges_gdf
    else:
        return gpd.GeoDataFrame(
            columns=['geometry', 'block_id', 'road_type', 'edge_index'],
            crs=blocks_gdf.crs
        )


def _find_nearest_road_type(
    edge: LineString,
    roads_gdf: gpd.GeoDataFrame,
    tolerance: float
) -> str:
    """
    Find the road type of the nearest road to an edge.

    Args:
        edge: Edge LineString
        roads_gdf: GeoDataFrame with roads
        tolerance: Maximum distance to consider

    Returns:
        Road type string or None if no road within tolerance
    """
    if roads_gdf.empty or 'road_type' not in roads_gdf.columns:
        return None

    # Get edge midpoint for distance calculation
    edge_midpoint = edge.interpolate(0.5, normalized=True)

    min_distance = float('inf')
    nearest_road_type = None

    for _, road_row in roads_gdf.iterrows():
        road_geom = road_row.geometry

        # Calculate distance from edge to road
        distance = edge.distance(road_geom)

        # Also check distance from midpoint (sometimes more accurate for adjacency)
        midpoint_distance = edge_midpoint.distance(road_geom)

        # Use the minimum of the two distances
        effective_distance = min(distance, midpoint_distance)

        if effective_distance < min_distance:
            min_distance = effective_distance
            nearest_road_type = road_row.get('road_type')

    # Only return road type if within tolerance
    if min_distance <= tolerance:
        return nearest_road_type
    else:
        return None


def extract_block_edges_simple(
    blocks_gdf: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Extract edges from blocks without road type assignment (simplified version).

    This is useful when you already have road type information in the blocks
    themselves and just need to decompose them into edges.

    Args:
        blocks_gdf: GeoDataFrame with block polygons

    Returns:
        GeoDataFrame with edges as LineStrings

    Example:
        >>> blocks = gpd.read_file('blocks.geojson')
        >>> edges = extract_block_edges_simple(blocks)
    """
    all_edges = []

    for idx, block_row in blocks_gdf.iterrows():
        block = block_row.geometry
        block_id = block_row.get('block_id', idx)

        if not isinstance(block, Polygon):
            continue

        # Get coordinates of the exterior ring
        coords = list(block.exterior.coords)

        # Create edge segments
        for i in range(len(coords) - 1):
            edge_geom = LineString([coords[i], coords[i + 1]])

            edge_data = {
                'geometry': edge_geom,
                'block_id': block_id,
                'edge_index': i
            }

            # Copy over block attributes
            for col in blocks_gdf.columns:
                if col not in ['geometry', 'edge_index']:
                    edge_data[col] = block_row.get(col)

            all_edges.append(edge_data)

    # Create GeoDataFrame
    if all_edges:
        edges_gdf = gpd.GeoDataFrame(all_edges, crs=blocks_gdf.crs)
        return edges_gdf
    else:
        return gpd.GeoDataFrame(
            columns=['geometry', 'block_id', 'edge_index'],
            crs=blocks_gdf.crs
        )


def assign_road_types_to_edges(
    edges_gdf: gpd.GeoDataFrame,
    roads_gdf: gpd.GeoDataFrame,
    tolerance: float = 5.0
) -> gpd.GeoDataFrame:
    """
    Assign road types to existing edges based on proximity to roads.

    Args:
        edges_gdf: GeoDataFrame with edge LineStrings
        roads_gdf: GeoDataFrame with roads (must have 'road_type' column)
        tolerance: Distance tolerance for matching edges to roads (default: 5m)

    Returns:
        GeoDataFrame with 'road_type' column added/updated

    Example:
        >>> edges = extract_block_edges_simple(blocks)
        >>> edges_with_roads = assign_road_types_to_edges(edges, roads, tolerance=5.0)
    """
    edges_gdf = edges_gdf.copy()

    road_types = []
    for _, edge_row in edges_gdf.iterrows():
        edge_geom = edge_row.geometry
        road_type = _find_nearest_road_type(edge_geom, roads_gdf, tolerance)
        road_types.append(road_type)

    edges_gdf['road_type'] = road_types

    return edges_gdf