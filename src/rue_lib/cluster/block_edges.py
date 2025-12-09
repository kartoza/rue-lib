# src/rue_lib/cluster/block_edges.py
"""Extract edges from block polygons and assign road types."""

import geopandas as gpd
from shapely.geometry import LineString, Polygon

from rue_lib.cluster.helpers import find_closest_road_type


def extract_block_edges(
        blocks_gdf: gpd.GeoDataFrame,
        roads_gdf: gpd.GeoDataFrame,
        tolerance
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
        >>> edges = extract_block_edges(blocks, roads, tolerance)
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
            road_type = find_closest_road_type(
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
                if col not in ['geometry', 'block_id', 'edge_index',
                               'road_type']:
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
        road_type = find_closest_road_type(edge_geom, roads_gdf, tolerance)
        road_types.append(road_type)

    edges_gdf['road_type'] = road_types

    return edges_gdf
