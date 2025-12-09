# src/rue_lib/cluster/cold/sub_div_on_grid_corners_convex.py
"""Subdivide on-grid parts at convex corners.

This module identifies convex corners along off-grid edges of on-grid parts
and creates corner subdivisions. Convex corners are vertices with angles > 200°
where both adjacent edges are longer than half the plot width.
"""

from pathlib import Path
from typing import Optional

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, MultiPolygon, Point, Polygon
from shapely.ops import unary_union

from rue_lib.cluster.block_edges import extract_block_edges
from rue_lib.cluster.helpers import compute_angle_dot
from rue_lib.core.definitions import RoadTypes


def calculate_vertex_angle(polygon: Polygon, vertex_idx: int) -> float:
    """
    Calculate the internal angle at a vertex in degrees.

    Uses the angle between normalized edge vectors to determine the
    internal angle at the vertex.

    Args:
        polygon: Input polygon
        vertex_idx: Index of vertex to calculate angle at

    Returns:
        Angle in degrees (0-360)
    """
    coords = list(polygon.exterior.coords)[:-1]
    n = len(coords)

    if n < 2:
        return 0.0

    # Get vertex position indices
    idx_prev = (vertex_idx - 1) % n
    idx_curr = vertex_idx % n
    idx_next = (vertex_idx + 1) % n

    # Get coordinates
    p_prev = np.array(coords[idx_prev])
    p_curr = np.array(coords[idx_curr])
    p_next = np.array(coords[idx_next])

    # Calculate edge vectors
    vec0 = p_curr - p_prev  # incoming edge
    vec1 = p_next - p_curr  # outgoing edge

    # Normalize vectors
    vec0_norm = vec0 / (np.linalg.norm(vec0) + 1e-10)
    vec1_norm = vec1 / (np.linalg.norm(vec1) + 1e-10)

    # Reverse vec0 to get outward direction
    vec0_rev = -vec0_norm

    # Calculate angle using atan2 for proper quadrant
    # This gives the angle from vec0_rev to vec1_norm
    cross = np.cross(np.append(vec1_norm, 0), np.append(vec0_rev, 0))[2]
    dot = np.dot(vec1_norm, vec0_rev)

    angle_rad = np.arctan2(cross, dot)
    angle_deg = np.degrees(angle_rad)

    # Convert to 0-360 range
    if angle_deg < 0:
        angle_deg += 360

    return angle_deg


def make_convex_corner(
        vertex_point: Point,
        edge_vec0: np.ndarray,
        edge_vec1: np.ndarray,
        part_depth: float
) -> Polygon:
    """
    Create a convex corner polygon at a vertex.

    Creates a triangular or quadrilateral corner piece by:
    1. Offsetting perpendicular to both edges
    2. Finding intersection of offset rays
    3. Creating polygon from vertex and offset points

    Args:
        vertex_point: The vertex point
        edge_vec0: Vector of first edge (incoming)
        edge_vec1: Vector of second edge (outgoing)
        part_depth: Offset distance for the corner

    Returns:
        Convex corner polygon
    """
    vert_xyz = np.array([vertex_point.x, vertex_point.y, 0])

    # Normalize edge vectors
    vec0 = edge_vec0[:2] / (np.linalg.norm(edge_vec0[:2]) + 1e-10)
    vec1 = edge_vec1[:2] / (np.linalg.norm(edge_vec1[:2]) + 1e-10)

    # Create perpendicular vectors (rotate 90° counter-clockwise)
    vec0_perp = np.array([-vec0[1], vec0[0], 0]) * (part_depth + 1)
    vec1_perp = np.array([-vec1[1], vec1[0], 0]) * (part_depth + 1)

    # Offset points
    xyz0 = vert_xyz + vec0_perp
    xyz1 = vert_xyz + vec1_perp

    # Create rays from offset points along edge directions
    # Find intersection of the two rays
    # Ray 0: starts at xyz0, direction vec0
    # Ray 1: starts at xyz1, direction vec1

    # Solve for intersection using parametric form:
    # xyz0 + t0 * vec0 = xyz1 + t1 * vec1
    # This becomes a 2D line intersection problem

    # Convert to 2D for intersection
    A = np.array([vec0[:2], -vec1[:2]]).T
    b = xyz1[:2] - xyz0[:2]

    try:
        params = np.linalg.solve(A, b)
        t0 = params[0]
        xyz_cor = xyz0 + t0 * np.append(vec0, 0)
    except np.linalg.LinAlgError:
        # Parallel lines, use midpoint
        xyz_cor = (xyz0 + xyz1) / 2

    # Create corner polygon
    corner_coords = [
        (vert_xyz[0], vert_xyz[1]),
        (xyz1[0], xyz1[1]),
        (xyz_cor[0], xyz_cor[1]),
        (xyz0[0], xyz0[1]),
        (vert_xyz[0], vert_xyz[1])
    ]

    return Polygon(corner_coords)


def fix_pinched_polygons(polygons: list[Polygon]) -> list[Polygon]:
    """
    Fix polygons that have self-touching vertices (pinched).

    When a polygon touches itself at a vertex, split it into multiple
    separate polygons at the pinch points.

    Args:
        polygons: List of potentially pinched polygons

    Returns:
        List of fixed polygons (split at pinch points)
    """
    all_new_polygons = []

    for pgon in polygons:
        if not isinstance(pgon, Polygon) or pgon.is_empty:
            continue

        coords = list(pgon.exterior.coords)[:-1]  # Remove duplicate last point

        # Find duplicate positions (pinch points)
        seen_coords = {}
        pinch_coords = []

        for i, coord in enumerate(coords):
            coord_tuple = (round(coord[0], 6), round(coord[1], 6))
            if coord_tuple in seen_coords:
                pinch_coords.append(i)
            else:
                seen_coords[coord_tuple] = i

        if not pinch_coords:
            # No pinching, keep original
            all_new_polygons.append(pgon)
            continue

        # Split at pinch points
        new_polygon_coords = [[]]
        used_pinches = set()

        for i, coord in enumerate(coords):
            coord_tuple = (round(coord[0], 6), round(coord[1], 6))

            if i in pinch_coords:
                if coord_tuple in used_pinches:
                    # Second occurrence, close current polygon
                    new_polygon_coords[-1].append(coord)
                    # Start new polygon
                    new_polygon_coords.append([])
                else:
                    # First occurrence, mark as used
                    used_pinches.add(coord_tuple)
                    new_polygon_coords[-1].append(coord)
            else:
                new_polygon_coords[-1].append(coord)

        # Create polygons from coordinate groups
        for poly_coords in new_polygon_coords:
            if len(poly_coords) >= 3:
                # Close the polygon
                poly_coords.append(poly_coords[0])
                try:
                    new_pgon = Polygon(poly_coords)
                    if new_pgon.is_valid and new_pgon.area > 1.0:
                        all_new_polygons.append(new_pgon)
                except Exception:
                    pass

    return all_new_polygons


def find_touching_edge(
        edges_gdf: gpd.GeoDataFrame,
        point: Point,
        tolerance: float = 0.1
) -> Optional[int]:
    """
    Find an edge that is very close to a given point.

    Args:
        edges_gdf: GeoDataFrame containing edges
        point: Point to check
        tolerance: Distance tolerance

    Returns:
        Index of touching edge or None
    """
    for idx, edge in edges_gdf.iterrows():
        if not isinstance(edge.geometry, LineString):
            continue

        distance = point.distance(edge.geometry)
        if distance < tolerance:
            return idx

    return None


def transfer_edge_attributes(
        from_edges_gdf: gpd.GeoDataFrame,
        to_edges_gdf: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Transfer road_type attributes from source edges to target edges that touch.

    Args:
        from_edges_gdf: Source edges with attributes
        to_edges_gdf: Target edges to update

    Returns:
        Updated target edges GeoDataFrame
    """
    for idx, to_edge in to_edges_gdf.iterrows():
        centroid = to_edge.geometry.centroid
        touching_idx = find_touching_edge(from_edges_gdf, centroid)

        if touching_idx is not None:
            from_edge = from_edges_gdf.loc[touching_idx]
            if 'road_type' in from_edge and from_edge['road_type'] is not None:
                to_edges_gdf.at[idx, 'road_type'] = from_edge['road_type']

    return to_edges_gdf


def copy_attributes(
        source_row: gpd.GeoSeries,
        target_gdf: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Copy block attributes from source to all target features.

    Args:
        source_row: Source feature with attributes
        target_gdf: Target features to update

    Returns:
        Updated target GeoDataFrame
    """
    for attr in ['block_id', 'block_type', 'site']:
        if attr in source_row.index and source_row[attr] is not None:
            target_gdf[attr] = source_row[attr]

    return target_gdf


def subdivide_part_at_convex_corners(
        part_row: gpd.GeoSeries,
        part_edges_gdf: gpd.GeoDataFrame,
        part_sec_d: float,
        part_loc_d: float,
        plot_sec_w: float,
        plot_loc_w: float
) -> tuple[list[Polygon], list[Polygon]]:
    """
    Subdivide a part by creating corner pieces at convex vertices.

    Identifies convex corners (angle > 200°) along off-grid edges and creates
    corner subdivisions, leaving strips for the remaining area.

    Args:
        part_row: Part feature with geometry and type attribute
        part_edges_gdf: GeoDataFrame with part edges
        part_sec_d: Secondary road part depth
        part_loc_d: Local road part depth
        plot_sec_w: Secondary plot width
        plot_loc_w: Local plot width

    Returns:
        Tuple of (strip_polygons, corner_polygons)
    """
    part = part_row.geometry
    part_type = part_row.get('type', 'sec')

    # Determine part depth and plot width based on type
    if part_type == 'loc':
        part_depth = part_loc_d
        plot_width = plot_loc_w
    else:
        part_depth = part_sec_d
        plot_width = plot_sec_w

    # Get edges by type
    if part_edges_gdf.empty:
        return [part], []

    # Find off-grid edges
    og_edges = part_edges_gdf[
        part_edges_gdf.get('road_type', '') == 'off_grid'
    ]

    if og_edges.empty:
        return [part], []

    # Process vertices on off-grid edges to find convex corners
    convex_corners = []
    coords = list(part.exterior.coords)[:-1]

    # For each off-grid edge, check vertices
    for _, og_edge in og_edges.iterrows():
        if not isinstance(og_edge.geometry, LineString):
            continue

        edge_coords = list(og_edge.geometry.coords)

        # Find this edge's vertices in the polygon
        for edge_coord in edge_coords:
            # Find matching vertex in polygon
            for i, poly_coord in enumerate(coords):
                if (abs(poly_coord[0] - edge_coord[0]) < 0.01 and
                    abs(poly_coord[1] - edge_coord[1]) < 0.01):

                    # Calculate angle at this vertex
                    angle = calculate_vertex_angle(part, i)

                    # Check if convex (angle > 200°)
                    ang_threshold = 200
                    if angle > ang_threshold:
                        # Check edge lengths at this vertex
                        idx_prev = (i - 1) % len(coords)
                        idx_next = (i + 1) % len(coords)

                        p_curr = np.array(coords[i])
                        p_prev = np.array(coords[idx_prev])
                        p_next = np.array(coords[idx_next])

                        len0 = np.linalg.norm(p_curr - p_prev)
                        len1 = np.linalg.norm(p_next - p_curr)

                        # Only create corner if both edges are long enough
                        if len0 > (plot_loc_w / 2) and len1 > (plot_loc_w / 2):
                            # Create convex corner
                            vec0 = p_curr - p_prev
                            vec1 = p_next - p_curr

                            vec0_3d = np.append(vec0, 0)
                            vec1_3d = np.append(vec1, 0)

                            corner = make_convex_corner(
                                Point(p_curr),
                                vec0_3d,
                                vec1_3d,
                                part_depth
                            )

                            if corner.is_valid and corner.area > 1.0:
                                convex_corners.append(corner)

    if not convex_corners:
        return [part], []

    # Union all corners
    corners_union = unary_union(convex_corners)

    # Intersect with part to get trimmed corners
    corners_trim = corners_union.intersection(part)

    # Subtract corners from part to get strips
    strips = part.difference(corners_union)

    # Convert to lists of polygons
    strip_polygons = []
    if isinstance(strips, Polygon):
        strip_polygons = [strips]
    elif isinstance(strips, MultiPolygon):
        strip_polygons = list(strips.geoms)

    corner_polygons = []
    if isinstance(corners_trim, Polygon):
        corner_polygons = [corners_trim]
    elif isinstance(corners_trim, MultiPolygon):
        corner_polygons = list(corners_trim.geoms)

    # Fix any pinched polygons
    strip_polygons = fix_pinched_polygons(strip_polygons)
    corner_polygons = fix_pinched_polygons(corner_polygons)

    # Filter by area
    strip_polygons = [p for p in strip_polygons if p.area > 1.0]
    corner_polygons = [p for p in corner_polygons if p.area > 1.0]

    return strip_polygons, corner_polygons


def process_part(
        part_row: gpd.GeoSeries,
        roads_gdf: gpd.GeoDataFrame,
        part_sec_d: float,
        part_loc_d: float,
        plot_sec_w: float,
        plot_loc_w: float
) -> list[dict]:
    """
    Process a single part to subdivide at convex corners.

    Args:
        part_row: Part GeoSeries with geometry and attributes
        roads_gdf: GeoDataFrame with roads
        part_sec_d: Secondary part depth
        part_loc_d: Local part depth
        plot_sec_w: Secondary plot width
        plot_loc_w: Local plot width

    Returns:
        List of dictionaries with subdivided part geometries and attributes
    """
    part = part_row.geometry
    part_type = part_row.get('type', 'sec')

    # Extract part edges (simplified - would need proper implementation)
    # For now, create a simple edge representation
    part_gdf = gpd.GeoDataFrame([part_row], geometry='geometry', crs=roads_gdf.crs)

    # Create edges from part boundary
    coords = list(part.exterior.coords)
    edges = []
    for i in range(len(coords) - 1):
        edge = LineString([coords[i], coords[i + 1]])
        edges.append({
            'geometry': edge,
            'road_type': 'off_grid'  # Default - would need proper classification
        })

    part_edges_gdf = gpd.GeoDataFrame(edges, crs=roads_gdf.crs)

    # Subdivide at convex corners
    strips, corners = subdivide_part_at_convex_corners(
        part_row,
        part_edges_gdf,
        part_sec_d,
        part_loc_d,
        plot_sec_w,
        plot_loc_w
    )

    # Create result list
    results = []

    # Add strips with original type
    for strip in strips:
        result = {
            'geometry': strip,
            'type': part_type,
            'class': 'part'
        }
        # Copy block attributes
        for attr in ['block_id', 'block_type', 'site']:
            if attr in part_row.index:
                result[attr] = part_row[attr]
        results.append(result)

    # Add corners with doubled type (e.g., 'sec_sec')
    corner_type = f"{part_type}_{part_type}"
    for corner in corners:
        result = {
            'geometry': corner,
            'type': corner_type,
            'class': 'part'
        }
        # Copy block attributes
        for attr in ['block_id', 'block_type', 'site']:
            if attr in part_row.index:
                result[attr] = part_row[attr]
        results.append(result)

    return results


def subdivide_on_grid_corners_convex(
        output_path: Path,
        parts_layer_name: str,
        roads_layer_name: str,
        part_sec_d: float,
        part_loc_d: float,
        plot_sec_w: float,
        plot_loc_w: float,
        output_layer_name: str
):
    """
    Subdivide on-grid parts at convex corners.

    Processes all on-grid parts (excluding off-grid) and creates corner
    subdivisions at convex vertices along off-grid edges. This separates
    corner plots from strip plots.

    Args:
        output_path: Path to GeoPackage file
        parts_layer_name: Name of parts layer
        roads_layer_name: Name of roads layer
        part_sec_d: Secondary part depth (typically 30m)
        part_loc_d: Local part depth (typically 20m)
        plot_sec_w: Secondary plot width
        plot_loc_w: Local plot width
        output_layer_name: Name for output layer

    Returns:
        Name of output layer
    """
    # Read input layers
    parts_layer = gpd.read_file(output_path, layer=parts_layer_name)
    roads_layer = gpd.read_file(output_path, layer=roads_layer_name)

    if parts_layer.empty:
        return output_layer_name

    # Filter to only on-grid parts (exclude off_grid)
    on_grid_parts = parts_layer[parts_layer['type'] != 'off_grid']

    if on_grid_parts.empty:
        return output_layer_name

    # Process each part
    all_subdivided = []

    for idx, part_row in on_grid_parts.iterrows():
        subdivided = process_part(
            part_row,
            roads_layer,
            part_sec_d,
            part_loc_d,
            plot_sec_w,
            plot_loc_w
        )
        all_subdivided.extend(subdivided)

    # Create output GeoDataFrame
    if all_subdivided:
        subdivided_gdf = gpd.GeoDataFrame(all_subdivided, crs=parts_layer.crs)
        subdivided_gdf.to_file(output_path, layer=output_layer_name, driver="GPKG")

    return output_layer_name