# src/rue_lib/cluster/runner_warm.py
"""Generate warm blocks with off-grid subdivision and partitioning."""

from __future__ import annotations

import math
from pathlib import Path

import geopandas as gpd
from shapely.geometry import LineString, Point

from rue_lib.cluster.cold.cluster_on_grid import (
    create_off_grid_cold_clusters,
    create_off_grid_zero_clusters,
    create_perpendicular_lines_from_front_points,
    extract_off_grid_adjacent_lines,
    extract_vertices_from_lines,
    merge_vertices_into_lines_by_angle,
    sample_points_along_front_lines,
)
from rue_lib.cluster.cold.clusters import (
    merge_and_classify_off_grid_clusters,
    merge_and_classify_on_grid_clusters,
    merge_final_cold_clusters,
)
from rue_lib.cluster.cold.expand_roads_buffer import (
    create_buffered_lines_from_boundary_lines,
)
from rue_lib.cluster.cold.subdiv_at_convex_corner import (
    create_clusters_from_convex_points,
    find_convex_points,
)
from rue_lib.cluster.cold.subdiv_block import (
    find_concave_points,
    subdivide_blocks_by_concave_points,
)
from rue_lib.cluster.config import ClusterConfig
from rue_lib.cluster.helpers import convert_polygonz_to_polygon
from rue_lib.core.definitions import BlockTypes
from rue_lib.core.geometry import merge_connected_lines
from rue_lib.streets.operations import extract_by_expression
from rue_lib.streets.runner_utils import subtract_layer


def generate_cold(
    cfg: ClusterConfig, output_gpkg: Path, input_blocks_layer_name: str, roads_layer_name: str
):
    """
    Generate cold blocks with off-grid subdivision and partitioning.

    This function processes blocks to create:
    - Inner off-grid areas by offsetting edges inward
    - Frame parts (perimeter around off-grid)
    - Corner and side parts from block decomposition
    - Subdivided plots within off-grid areas
    - On-grid parts for arterial and secondary roads

    Args:
        cfg: ClusterConfig with partition settings
        output_gpkg: Path to output GeoPackage
        input_blocks_layer_name: Name of input blocks layer
        roads_layer_name: Name of roads layer
    """
    _part_art_d = cfg.on_grid_partition_depth_arterial_roads
    _part_sec_d = cfg.on_grid_partition_depth_secondary_roads
    _part_loc_d = cfg.on_grid_partition_depth_local_roads

    # TODO:
    #  We use the same width for off-grid and on-grid plots for now.
    _part_og_w = cfg.off_grid_cluster_width

    output_path = str(output_gpkg)
    print("==============================================================")
    print("COLD BLOCK")
    print("==============================================================")
    print("Step 1: Extract cold grid blocks...")
    cold_grid_layer_name = "200_cold_grid"
    extract_by_expression(
        output_path,
        input_blocks_layer_name,
        f"type = '{BlockTypes.COLD_GRID}'",
        output_path,
        cold_grid_layer_name,
    )
    # Convert from PolygonZ to Polygon
    convert_polygonz_to_polygon(output_path, cold_grid_layer_name)

    erased_layer_name = cold_grid_layer_name

    print("\nStep 2: Extract boundary lines adjacent to roads...")
    boundary_points_layer_name = "202_cold_boundary_points"
    boundary_lines_layer_name = "202_cold_boundary_lines"
    boundary_points_layer_name, boundary_lines_layer_name = extract_road_adjacent_vertices(
        output_path,
        erased_layer_name,
        roads_layer_name,
        output_path,
        boundary_points_layer_name,
    )
    merge_connected_lines(output_path, boundary_lines_layer_name)

    print("\nStep 4: Subtract local road buffer from cold grid...")
    cold_grid_updated_layer_name = "203_cold_grid_updated"
    subtract_layer(
        output_gpkg,
        cold_grid_layer_name,
        boundary_lines_layer_name,
        output_gpkg,
        cold_grid_updated_layer_name,
        cfg.on_grid_partition_depth_local_roads,
    )

    print("\nStep 5: Extract and subtract secondary road buffer...")
    extract_by_expression(
        output_path,
        "002_input_roads_buffer",
        "type = 'road_secondary'",
        output_path,
        "002_input_roads_buffer_sec",
    )

    subtract_layer(
        output_gpkg,
        cold_grid_updated_layer_name,
        "002_input_roads_buffer_sec",
        output_gpkg,
        cold_grid_updated_layer_name,
        cfg.on_grid_partition_depth_secondary_roads,
    )

    print("\nStep 6: Extract and subtract arterial road buffer...")
    extract_by_expression(
        output_path,
        "002_input_roads_buffer",
        "type = 'road_arterial'",
        output_path,
        "002_input_roads_buffer_art",
    )

    subtract_layer(
        output_gpkg,
        cold_grid_updated_layer_name,
        "002_input_roads_buffer_art",
        output_gpkg,
        cold_grid_updated_layer_name,
        cfg.on_grid_partition_depth_arterial_roads,
    )

    print("\nStep 7: Merge boundary vertices into lines by angle...")
    boundary_lines_from_vertices = "204_boundary_lines_from_vertices"
    merge_vertices_into_lines_by_angle(
        output_path,
        boundary_points_layer_name,
        output_path,
        boundary_lines_from_vertices,
    )

    print("\nStep 8: Create buffered lines from boundary lines...")
    buffered_lines_layer_name = "205_buffered_lines"
    create_buffered_lines_from_boundary_lines(
        output_path,
        "002_input_roads_buffer",
        boundary_lines_from_vertices,
        output_path,
        buffered_lines_layer_name,
        cfg,
    )

    print("\nStep 9: Clip buffered lines to cold grid...")
    clipped_lines_layer_name = "207_clipped_buffered_lines"
    subtract_layer(
        output_gpkg,
        cold_grid_layer_name,
        cold_grid_updated_layer_name,
        output_gpkg,
        clipped_lines_layer_name,
        0,
    )

    print("\nStep 10: Find concave points from boundary...")
    concave_points_layer_name = "208_concave_points"
    find_concave_points(
        output_path,
        erased_layer_name,
        boundary_points_layer_name,
        output_path,
        concave_points_layer_name,
    )

    print("\nStep 11: Subdivide blocks at concave corners...")
    cutting_lines_layer_name = "209_subdivided"
    subdivide_blocks_by_concave_points(
        output_path,
        erased_layer_name,
        concave_points_layer_name,
        boundary_points_layer_name,
        output_path,
        cutting_lines_layer_name,
        cfg.road_local_width_m / 2.0,
        clipped_lines_layer_name,
    )

    print("\nStep 12: Extract off-grid boundaries adjacent to on-grid blocks...")
    off_grid_block = "209_subdivided_blocks_off_grid"
    on_grid_block = "209_subdivided_blocks_on_grid"

    adjacent_off_grid_lines = "210_off_grid_adjacent_lines"
    extract_off_grid_adjacent_lines(
        output_path,
        off_grid_block,
        on_grid_block,
        output_path,
        adjacent_off_grid_lines,
    )

    print("\nStep 13: Extract vertices from off-grid adjacent lines...")
    extract_vertices_layer_name = "211_off_grid_adjacent_vertices"
    extract_vertices_from_lines(
        output_path,
        adjacent_off_grid_lines,
        extract_vertices_layer_name,
    )

    print("\nStep 14: Merge off-grid vertices into lines by angle...")
    lines_from_vertices_layer = "212_lines_from_vertices"
    merge_vertices_into_lines_by_angle(
        output_path,
        extract_vertices_layer_name,
        output_path,
        lines_from_vertices_layer,
    )

    print("\nStep 15: Find convex points from boundary...")
    convex_points_layer_name = "213_convex_points"
    find_convex_points(
        output_path,
        extract_vertices_layer_name,
        output_path,
        convex_points_layer_name,
    )

    print("\nStep 16: Create clusters from convex points...")
    convex_clusters_layer = "214_convex_clusters"
    max_partition_depth = max(
        cfg.on_grid_partition_depth_arterial_roads,
        cfg.on_grid_partition_depth_secondary_roads,
    )
    create_clusters_from_convex_points(
        output_path,
        convex_points_layer_name,
        on_grid_block,
        output_path,
        convex_clusters_layer,
        max_partition_depth,
    )

    print("\nStep 17: Sample points along front lines...")
    cold_clusters_points_layer = "215_off_grid_cold_clusters_points"
    sample_points_along_front_lines(
        output_path,
        lines_from_vertices_layer,
        None,
        output_path,
        cold_clusters_points_layer,
        width_m=float(cfg.off_grid_cluster_width),
    )

    print("\nStep 18: Create perpendicular lines from front line points...")
    perpendicular_lines_layer = "216_off_grid_perpendicular_lines"

    create_perpendicular_lines_from_front_points(
        output_path,
        cold_clusters_points_layer,
        lines_from_vertices_layer,
        off_grid_block,
        output_path,
        perpendicular_lines_layer,
    )

    print("\nStep 19: Add depth points to perpendicular lines...")
    depth_points_layer = "217_off_grid_depth_points"
    sample_points_along_front_lines(
        output_path,
        perpendicular_lines_layer,
        cold_clusters_points_layer,
        output_path,
        depth_points_layer,
        width_m=float(cfg.off_grid_cluster_depth),
        max_depth=1,
    )

    print("\nStep 20: Split off-grid blocks into preliminary clusters (off_grid0)...")
    off_grid0_layer = "218_off_grid0_clusters"
    create_off_grid_zero_clusters(
        output_path,
        off_grid_block,
        perpendicular_lines_layer,
        off_grid0_layer,
    )

    print("\nStep 21: Create cluster polygons from depth points...")
    clusters_layer = "219_off_grid_cold_clusters"
    create_off_grid_cold_clusters(
        output_path,
        off_grid0_layer,
        depth_points_layer,
        output_path,
        clusters_layer,
        perpendicular_lines_layer,
        buffer_distance=cfg.off_grid_cluster_width * 0.75,
        target_area_m2=cfg.off_grid_cluster_width * cfg.off_grid_cluster_depth,
    )

    print("\nStep 22: Merge and classify off-grid cold clusters...")
    off_grid_final_layer = "220_final_cold_off_grid_clusters"
    merge_and_classify_off_grid_clusters(
        output_path,
        off_grid_block,
        clusters_layer,
        output_path,
        off_grid_final_layer,
    )

    print("\nStep 23: Merge and classify on-grid cold clusters...")
    on_grid_final_layer = "221_final_cold_on_grid_clusters"
    merge_and_classify_on_grid_clusters(
        output_path,
        on_grid_block,
        convex_clusters_layer,
        concave_points_layer_name,
        "002_input_roads_buffer",
        output_path,
        on_grid_final_layer,
    )

    print("\nStep 24: Merge final on-grid and off-grid cold clusters...")
    final_clusters_layer = "222_final_cold_clusters"
    merge_final_cold_clusters(
        output_path,
        on_grid_final_layer,
        off_grid_final_layer,
        "209_subdivided_blocks",
        output_path,
        final_clusters_layer,
    )

    return final_clusters_layer


def get_priority_road_type(road_types: list[str]) -> str:
    """Get the highest priority road type from a list.

    Priority order: local > secondary > arterial
    """
    if "local" in road_types:
        return "local"
    elif "secondary" in road_types:
        return "secondary"
    elif "arterial" in road_types:
        return "arterial"
    else:
        return road_types[0] if road_types else "unknown"


def extract_road_adjacent_vertices(
    input_gpkg: str,
    erased_grid_layer_name: str,
    roads_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> tuple[str, str]:
    """
    Extract vertices from erased cold grid boundary that intersect with the roads buffer.

    Writes three layers:
    - A points layer containing the road-adjacent vertices
    - A lines layer containing the boundary lines that touch the roads buffer

    Returns:
        (points_layer_name, lines_layer_name)
    """
    gdf_grid = gpd.read_file(input_gpkg, layer=erased_grid_layer_name)
    gdf_roads = gpd.read_file(input_gpkg, layer=roads_layer_name)

    if gdf_grid.empty:
        raise ValueError(f"Layer {erased_grid_layer_name} is empty")
    if gdf_roads.empty:
        raise ValueError(f"Layer {roads_layer_name} is empty")

    print(f"  Processing {len(gdf_grid)} blocks...")

    vertices_to_write = []
    lines_layer_name = (
        output_layer_name.replace("points", "lines")
        if "points" in output_layer_name
        else f"{output_layer_name}_lines"
    )
    block_boundaries = {}

    for block_idx, grid_row in gdf_grid.iterrows():
        block_id = grid_row.get("id", block_idx)
        grid_geom = grid_row.geometry

        if grid_geom is None or grid_geom.is_empty:
            continue
        boundary = grid_geom.boundary
        if boundary is None or boundary.is_empty:
            continue
        if boundary.geom_type == "LineString":
            lines = [boundary]
        elif boundary.geom_type == "MultiLineString":
            lines = list(boundary.geoms)
        else:
            continue
        for line in lines:
            coords = list(line.coords)

            # Store all boundary coordinates for this block
            block_boundaries[block_id] = grid_geom

            for vertex_id, (x, y) in enumerate(coords):
                point = Point(x, y)
                point_buffer = point.buffer(0.1)

                # Check which roads this vertex touches
                for _, road_row in gdf_roads.iterrows():
                    if point_buffer.intersects(road_row.geometry):
                        vertices_to_write.append(
                            {
                                "x": x,
                                "y": y,
                                "block_id": block_id,
                                "orig_id": block_id,
                                "vertex_id": vertex_id,
                                "road_type": road_row["type"],
                            }
                        )
                        break

    # Group vertices by block_id only (not by road_type)
    # For each block, determine the priority road_type
    vertices_grouped: dict[int, list[dict]] = {}
    for vertex in vertices_to_write:
        block_id = vertex["block_id"]
        vertices_grouped.setdefault(block_id, []).append(vertex)

    # Build line segments by connecting adjacent road-touching vertices
    # Create individual line segments one at a time between consecutive vertices
    lines_to_write = []
    points_to_write = []

    for block_id, verts in vertices_grouped.items():
        road_types = list({v["road_type"] for v in verts})
        priority_road_type = get_priority_road_type(road_types)
        verts_sorted = sorted(verts, key=lambda v: v["vertex_id"])

        for i in range(len(verts_sorted) - 1):
            v1 = verts_sorted[i]
            v2 = verts_sorted[i + 1]

            if v2["vertex_id"] == v1["vertex_id"] + 1:
                coords = [(v1["x"], v1["y"]), (v2["x"], v2["y"])]
                line = LineString(coords)
                line_centroid = line.centroid.buffer(1)
                centroid_intersects_road = any(
                    line_centroid.intersects(road_row.geometry)
                    for _, road_row in gdf_roads.iterrows()
                )
                if centroid_intersects_road:
                    lines_to_write.append(
                        {
                            "geometry": line,
                            "block_id": block_id,
                            "orig_id": block_id,
                            "road_type": priority_road_type,
                        }
                    )

        boundary_geom = block_boundaries.get(block_id)
        if boundary_geom is None:
            continue

        boundary = boundary_geom.boundary
        if boundary is None or boundary.is_empty:
            continue

        if boundary.geom_type == "LineString":
            boundary_coords = list(boundary.coords)
        elif boundary.geom_type == "MultiLineString":
            boundary_coords = list(boundary.geoms[0].coords)
        else:
            continue

        num_boundary_verts = len(boundary_coords)

        for v in verts_sorted:
            vertex_id = v["vertex_id"]
            if num_boundary_verts >= 3 and vertex_id < num_boundary_verts:
                min_distance = 5.0  # meters
                current_x, current_y = boundary_coords[vertex_id]
                prev_vertex_id = (vertex_id - 1) % num_boundary_verts
                steps = 0
                while steps < num_boundary_verts - 2:
                    prev_x, prev_y = boundary_coords[prev_vertex_id]
                    dist_to_prev = math.sqrt((current_x - prev_x) ** 2 + (current_y - prev_y) ** 2)
                    if dist_to_prev >= min_distance:
                        break
                    prev_vertex_id = (prev_vertex_id - 1) % num_boundary_verts
                    steps += 1

                next_vertex_id = (vertex_id + 1) % num_boundary_verts
                steps = 0
                while steps < num_boundary_verts - 2:
                    next_x, next_y = boundary_coords[next_vertex_id]
                    dist_to_next = math.sqrt((current_x - next_x) ** 2 + (current_y - next_y) ** 2)
                    if dist_to_next >= min_distance:
                        break
                    next_vertex_id = (next_vertex_id + 1) % num_boundary_verts
                    steps += 1

                prev_x, prev_y = boundary_coords[prev_vertex_id]
                next_x, next_y = boundary_coords[next_vertex_id]

                cross_line = LineString([(prev_x, prev_y), (next_x, next_y)])

                is_interior = False
                line_length = cross_line.length
                if line_length > 0:
                    check_fractions = [0.15, 0.5, 0.85]
                    for fraction in check_fractions:
                        point_along_line = cross_line.interpolate(fraction, normalized=True)
                        point_buffer = point_along_line.buffer(0.1)
                        if point_buffer.intersects(boundary_geom):
                            is_interior = True
                            break

                v1_x = current_x - prev_x
                v1_y = current_y - prev_y

                v2_x = next_x - current_x
                v2_y = next_y - current_y

                dot = v1_x * v2_x + v1_y * v2_y
                cross = v1_x * v2_y - v1_y * v2_x
                angle_rad = math.atan2(cross, dot)
                angle_deg = math.degrees(angle_rad)
                if is_interior:
                    angle_deg = abs(angle_deg)
                else:
                    angle_deg = -abs(angle_deg)
            else:
                angle_deg = 0.0

            points_to_write.append(
                {
                    "geometry": Point(v["x"], v["y"]),
                    "block_id": block_id,
                    "vertex_id": v["vertex_id"],
                    "vertex_idx": v["vertex_id"],
                    "orig_id": block_id,
                    "road_type": priority_road_type,
                    "line_id": block_id,
                    "angle_deg": angle_deg,
                    "next_x": next_x,
                    "next_y": next_y,
                    "prev_x": prev_x,
                    "prev_y": prev_y,
                }
            )

    if lines_to_write:
        gdf_lines = gpd.GeoDataFrame(lines_to_write, geometry="geometry", crs=gdf_grid.crs)
    else:
        gdf_lines = gpd.GeoDataFrame([], geometry=[], crs=gdf_grid.crs)
    gdf_lines.to_file(output_gpkg, layer=lines_layer_name, driver="GPKG")

    if points_to_write:
        gdf_points = gpd.GeoDataFrame(points_to_write, geometry="geometry", crs=gdf_grid.crs)
    else:
        gdf_points = gpd.GeoDataFrame([], geometry=[], crs=gdf_grid.crs)
    gdf_points.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")

    total_vertices = len(points_to_write)
    total_lines = len(lines_to_write)
    print(f"  Extracted {total_vertices} vertices from road-adjacent boundaries")
    print(f"  Created points layer: {output_layer_name}")
    print(f"  Created lines layer: {lines_layer_name} ({total_lines} features)")
    return output_layer_name, lines_layer_name
