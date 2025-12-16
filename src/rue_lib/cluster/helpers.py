# src/rue_lib/cluster/helpers.py
"""Helper functions for cluster generation and geometry operations."""

from typing import Optional

import geopandas as gpd
import numpy as np
from osgeo import ogr
from shapely.geometry import LineString, MultiPolygon, Point, Polygon

from rue_lib.core.definitions import BlockTypes, ClusterTypes, PropertyKeys, RoadTypes


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

    vec0_norm = vec0 / np.linalg.norm(vec0) if np.linalg.norm(vec0) > 0 else vec0
    vec1_norm = vec1 / np.linalg.norm(vec1) if np.linalg.norm(vec1) > 0 else vec1

    return np.dot(vec0_norm, vec1_norm)


def convert_to_quadrilateral(polygon: Polygon, min_area: float = 100.0) -> Optional[Polygon]:
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


def find_closest_road_type(
    edge: LineString | Point,
    roads: gpd.GeoDataFrame,
    max_distance: float = 1,
    default_type: Optional[RoadTypes] = RoadTypes.Local,
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
    if roads.empty or PropertyKeys.RoadType not in roads.columns:
        return None

    # Use the center point of the edge
    if isinstance(edge, Point):
        edge_center = edge
    else:
        edge_center = edge.interpolate(0.5, normalized=True)

    min_distance = float("inf")
    closest_type = default_type

    for _, road in roads.iterrows():
        dist = edge_center.distance(road.geometry)
        if dist < min_distance and dist <= max_distance:
            min_distance = dist
            new_type = road[PropertyKeys.RoadType]
            if (
                new_type != RoadTypes.Local
                or new_type == RoadTypes.Local
                and closest_type not in [RoadTypes.Artery, RoadTypes.Secondary]
            ):
                closest_type = road[PropertyKeys.RoadType]

    # Return None if the closest road is too far
    if min_distance > max_distance:
        return None

    if closest_type == "road_art":
        closest_type = RoadTypes.Artery
    elif closest_type == "road_sec":
        closest_type = RoadTypes.Secondary
    elif closest_type == "road_loc":
        closest_type = RoadTypes.Local
    return closest_type


def get_edge_angle(coords: list, vertex_idx: int) -> float:
    """
    Calculate the internal angle at a vertex in a polygon.

    Computes the angle between the incoming and outgoing edges at the
    specified vertex using the dot product of normalized edge vectors.

    Args:
        coords: List of polygon coordinates
        vertex_idx: Index of vertex to calculate angle at

    Returns:
        Internal angle in degrees (0-180)
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


def convert_polygonz_to_polygon(gpkg_path: str, layer_name: str):
    """
    Convert PolygonZ geometries to Polygon (remove Z coordinate).

    Args:
        gpkg_path: Path to GeoPackage file
        layer_name: Name of layer to convert
    """
    # Read the layer
    gdf = gpd.read_file(gpkg_path, layer=layer_name)

    # Convert geometries from PolygonZ to Polygon (drop Z coordinate)
    def drop_z(geom):
        if geom.has_z:
            # Extract 2D coordinates only
            if geom.geom_type == "Polygon":
                return Polygon([(x, y) for x, y, *_ in geom.exterior.coords])
            elif geom.geom_type == "MultiPolygon":
                return MultiPolygon(
                    [Polygon([(x, y) for x, y, *_ in poly.exterior.coords]) for poly in geom.geoms]
                )
        return geom

    gdf["geometry"] = gdf["geometry"].apply(drop_z)

    # Save back to the same layer (overwrite)
    gdf.to_file(gpkg_path, layer=layer_name, driver="GPKG")


def assign_cluster_type(gpkg_path: str, block_layer_name: str, final_layer_name: str):
    """
    Assign cluster_type from block layer to final layer based on geometry match.

    For each feature in final layer, finds matching feature in block layer
    and assigns its type as cluster_type.

    Args:
        gpkg_path: Path to GeoPackage containing the layers
        block_layer_name: Name of the layer containing block types
        final_layer_name: Name of the final layer to update with cluster_type
    """
    print(f"Assigning cluster types from {block_layer_name} to {final_layer_name}...")

    # Open the GeoPackage for reading and writing
    ds = ogr.Open(gpkg_path, 1)
    if ds is None:
        raise ValueError(f"Could not open {gpkg_path} for writing")

    final_layer = ds.GetLayerByName(final_layer_name)
    if final_layer is None:
        raise ValueError(f"Layer {final_layer_name} not found")

    block_layer = ds.GetLayerByName(block_layer_name)
    if block_layer is None:
        raise ValueError(f"Layer {block_layer_name} not found")

    # Create cluster_type field if it doesn't exist
    layer_defn = final_layer.GetLayerDefn()
    if layer_defn.GetFieldIndex("cluster_type") == -1:
        final_layer.CreateField(ogr.FieldDefn("cluster_type", ogr.OFTString))
        print("  Created cluster_type field")

    updated_count = 0

    # Loop through final layer
    for final_feat in final_layer:
        final_geom = final_feat.GetGeometryRef()
        if final_geom is None:
            continue

        # Find matching feature in block layer by centroid
        final_geom_centroid = final_geom.Centroid()
        block_layer.SetSpatialFilter(final_geom_centroid)
        block_type = None
        for block_feat in block_layer:
            block_geom = block_feat.GetGeometryRef()

            # Check if final geometry centroid is within block geometry
            if block_geom and final_geom_centroid.Within(block_geom):
                block_type = block_feat.GetField("type")
                break

        block_layer.ResetReading()

        # Assign cluster_type if found
        if block_type is not None:
            # Collect off-grid part
            if block_type == BlockTypes.ON_GRID_ART:
                block_type = ClusterTypes.ON_GRID_ART
            elif block_type == BlockTypes.ON_GRID_SEC:
                block_type = ClusterTypes.ON_GRID_SEC
            elif block_type == BlockTypes.OFF_GRID:
                block_type = ClusterTypes.ON_GRID_LOC

            final_feat.SetField("cluster_type", block_type)
            final_layer.SetFeature(final_feat)
            updated_count += 1

    print(f"  Updated {updated_count} features with cluster_type")

    # Clean up
    ds = None

    print(f"Cluster type assignment completed for layer: {final_layer_name}")
