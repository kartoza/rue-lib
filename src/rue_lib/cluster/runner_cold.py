# src/rue_lib/cluster/runner_warm.py
"""Generate warm blocks with off-grid subdivision and partitioning."""

from __future__ import annotations

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
    print("Step 1: Generate inner part of off grid blocks...")
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

    print("\nStep 2: Erase roads buffer from cold grid...")
    erased_layer_name = cold_grid_layer_name

    print("\nStep 3: Extract boundary lines adjacent to roads...")
    boundary_points_layer_name = "202_cold_boundary_points"
    boundary_points_layer_name, boundary_lines_layer_name = extract_road_adjacent_vertices(
        output_path,
        erased_layer_name,
        roads_layer_name,
        output_path,
        boundary_points_layer_name,
    )

    subtract_layer(
        output_gpkg,
        "200_cold_grid",
        "202_cold_boundary_lines",
        output_gpkg,
        "200_cold_grid_updated",
        cfg.on_grid_partition_depth_local_roads,
    )

    extract_by_expression(
        output_path,
        "002_input_roads_buffer",
        "type = 'road_secondary'",
        output_path,
        "002_input_roads_buffer_sec",
    )

    subtract_layer(
        output_gpkg,
        "200_cold_grid_updated",
        "002_input_roads_buffer_sec",
        output_gpkg,
        "200_cold_grid_updated",
        cfg.on_grid_partition_depth_secondary_roads,
    )

    extract_by_expression(
        output_path,
        "002_input_roads_buffer",
        "type = 'road_arterial'",
        output_path,
        "002_input_roads_buffer_art",
    )

    subtract_layer(
        output_gpkg,
        "200_cold_grid_updated",
        "002_input_roads_buffer_art",
        output_gpkg,
        "200_cold_grid_updated",
        cfg.on_grid_partition_depth_arterial_roads,
    )

    boundary_lines_from_vertices = "202_boundary_lines_from_vertices"
    merge_vertices_into_lines_by_angle(
        output_path,
        boundary_points_layer_name,
        output_path,
        boundary_lines_from_vertices,
    )

    print("\nStep 6: Create buffered lines from boundary lines...")
    buffered_lines_layer_name = "206_buffered_lines"
    create_buffered_lines_from_boundary_lines(
        output_path,
        "002_input_roads_buffer",
        boundary_lines_from_vertices,
        output_path,
        buffered_lines_layer_name,
        cfg,
    )

    print("\nStep 7: Clip buffered lines to cold grid...")
    clipped_lines_layer_name = "207_clipped_buffered_lines"
    subtract_layer(
        output_gpkg,
        "200_cold_grid",
        "200_cold_grid_updated",
        output_gpkg,
        clipped_lines_layer_name,
        0,
    )

    # clip_buffered_lines_to_cold_grid(
    #     output_path,
    #     buffered_lines_layer_name,
    #     erased_layer_name,
    #     output_path,
    #     clipped_lines_layer_name,
    # )

    print("\nStep 4: Find concave points from boundary...")
    concave_points_layer_name = "203_concave_points"
    find_concave_points(
        output_path,
        erased_layer_name,
        boundary_points_layer_name,
        output_path,
        concave_points_layer_name,
    )

    print("\nStep 5: Subdivide blocks at concave corners...")
    cutting_lines_layer_name = "204_subdivided"
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

    print("\nStep 6: Create cluster on on-grid parts...")
    off_grid_block = "204_subdivided_blocks_off_grid"
    on_grid_block = "204_subdivided_blocks_on_grid"

    print("  Extracting off-grid boundaries adjacent to on-grid blocks...")
    adjacent_off_grid_lines = "205_off_grid_adjacent_lines"
    extract_off_grid_adjacent_lines(
        output_path,
        off_grid_block,
        on_grid_block,
        output_path,
        adjacent_off_grid_lines,
    )

    extract_vertices_layer_name = "208_off_grid_adjacent_vertices"
    extract_vertices_from_lines(
        output_path,
        adjacent_off_grid_lines,
        extract_vertices_layer_name,
    )

    lines_from_vertices_layer = "209_lines_from_vertices"
    merge_vertices_into_lines_by_angle(
        output_path,
        extract_vertices_layer_name,
        output_path,
        lines_from_vertices_layer,
    )

    print("\nStep 6b: Find convex points from boundary...")
    convex_points_layer_name = "208b_convex_points"
    find_convex_points(
        output_path,
        extract_vertices_layer_name,
        output_path,
        convex_points_layer_name,
    )

    print("\nStep 6c: Create clusters from convex points...")
    convex_clusters_layer = "208c_convex_clusters"
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

    print("\nStep 7: Sample points along front lines...")
    cold_clusters_points_layer = "210_off_grid_cold_clusters_points"

    sample_points_along_front_lines(
        output_path,
        lines_from_vertices_layer,
        None,
        output_path,
        cold_clusters_points_layer,
        width_m=float(cfg.off_grid_cluster_width),
    )

    print("\nStep 8: Create perpendicular lines from front line points...")
    perpendicular_lines_layer = "211_off_grid_perpendicular_lines"

    create_perpendicular_lines_from_front_points(
        output_path,
        cold_clusters_points_layer,
        lines_from_vertices_layer,
        off_grid_block,
        output_path,
        perpendicular_lines_layer,
    )

    print("\nStep 10: Add depth points to perpendicular lines for off-grid clustering...")
    depth_points_layer = "212_off_grid_depth_points"
    sample_points_along_front_lines(
        output_path,
        perpendicular_lines_layer,
        cold_clusters_points_layer,
        output_path,
        depth_points_layer,
        width_m=float(cfg.off_grid_cluster_depth),
        max_depth=1,
    )

    print("\nStep 9: Split off-grid blocks into preliminary clusters (off_grid0)...")
    off_grid0_layer = "213_off_grid0_clusters"
    create_off_grid_zero_clusters(
        output_path,
        off_grid_block,
        perpendicular_lines_layer,
        off_grid0_layer,
    )

    print("\nStep 11: Create cluster polygons from depth points...")
    clusters_layer = "214_off_grid_cold_clusters"
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

    print("\nStep 12: Merge and classify off-grid cold clusters...")
    off_grid_final_layer = "215_final_cold_off_grid_clusters"
    merge_and_classify_off_grid_clusters(
        output_path,
        off_grid_block,
        clusters_layer,
        output_path,
        off_grid_final_layer,
    )

    print("\nStep 13: Merge and classify on-grid cold clusters...")
    on_grid_final_layer = "216_final_cold_on_grid_clusters"
    merge_and_classify_on_grid_clusters(
        output_path,
        on_grid_block,
        convex_clusters_layer,
        concave_points_layer_name,
        "002_input_roads_buffer",
        output_path,
        on_grid_final_layer,
    )

    print("\nStep 14: Merge final on-grid and off-grid cold clusters...")
    final_clusters_layer = "217_final_cold_clusters"
    merge_final_cold_clusters(
        output_path,
        on_grid_final_layer,
        off_grid_final_layer,
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

    Writes both:
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
    lines_to_write = []
    points_to_write = []

    for block_id, verts in vertices_grouped.items():
        road_types = list({v["road_type"] for v in verts})
        priority_road_type = get_priority_road_type(road_types)
        verts_sorted = sorted(verts, key=lambda v: v["vertex_id"])

        current = []
        prev_id = None

        for v in verts_sorted:
            if prev_id is None or v["vertex_id"] == prev_id + 1:
                current.append(v)
            else:
                if len(current) >= 2:
                    coords = [(pt["x"], pt["y"]) for pt in current]
                    line = LineString(coords)
                    lines_to_write.append(
                        {
                            "geometry": line,
                            "block_id": block_id,
                            "orig_id": block_id,
                            "road_type": priority_road_type,
                        }
                    )
                current = [v]
            prev_id = v["vertex_id"]

        # Flush remaining vertices
        if len(current) >= 2:
            coords = [(pt["x"], pt["y"]) for pt in current]
            line = LineString(coords)
            lines_to_write.append(
                {
                    "geometry": line,
                    "block_id": block_id,
                    "orig_id": block_id,
                    "road_type": priority_road_type,
                }
            )

        for v in verts:
            points_to_write.append(
                {
                    "geometry": Point(v["x"], v["y"]),
                    "block_id": block_id,
                    "vertex_id": v["vertex_id"],
                    "vertex_idx": v["vertex_id"],
                    "orig_id": block_id,
                    "road_type": priority_road_type,
                    "line_id": block_id,
                    "angle_deg": 0.0,
                }
            )

    # Write lines layer
    if lines_to_write:
        gdf_lines = gpd.GeoDataFrame(lines_to_write, geometry="geometry", crs=gdf_grid.crs)
        gdf_lines.to_file(output_gpkg, layer=lines_layer_name, driver="GPKG")

    # Write points layer
    if points_to_write:
        gdf_points = gpd.GeoDataFrame(points_to_write, geometry="geometry", crs=gdf_grid.crs)
        gdf_points.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")

    total_vertices = len(points_to_write)
    total_lines = len(lines_to_write)
    print(f"  Extracted {total_vertices} vertices from road-adjacent boundaries")
    print(f"  Created lines layer: {lines_layer_name} ({total_lines} features)")
    print(f"  Created points layer: {output_layer_name}")
    return output_layer_name, lines_layer_name
