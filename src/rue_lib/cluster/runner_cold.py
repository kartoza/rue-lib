# src/rue_lib/cluster/runner_warm.py
"""Generate warm blocks with off-grid subdivision and partitioning."""

from __future__ import annotations

from pathlib import Path

from rue_lib.cluster.cold.add_on_grid_strips_art_sec_loc import generate_on_grid_strips_art_sec_loc
from rue_lib.cluster.cold.sub_div_on_grid_corners_convex import subdivide_on_grid_corners_convex
from rue_lib.cluster.cold.subdiv_at_concave_corner import subdivide_blocks_at_concave_corners
from rue_lib.cluster.config import ClusterConfig
from rue_lib.cluster.helpers import convert_polygonz_to_polygon
from rue_lib.core.definitions import BlockTypes
from rue_lib.streets.operations import extract_by_expression


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
    part_art_d = cfg.on_grid_partition_depth_arterial_roads
    part_sec_d = cfg.on_grid_partition_depth_secondary_roads
    part_loc_d = cfg.on_grid_partition_depth_local_roads

    # TODO:
    #  We use the same width for off-grid and on-grid plots for now.
    part_og_w = cfg.off_grid_cluster_width

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

    cold_grid_subdive_at_concave_layer_name = "201_cold_grid_subdive_at_concave"
    subdivide_blocks_at_concave_corners(
        output_path=output_path,
        input_layer_name=cold_grid_layer_name,
        output_layer_name=cold_grid_subdive_at_concave_layer_name,
        roads_layer_name=roads_layer_name,
    )
    cold_grid_block_layer_name = "202_cold_block"
    extract_by_expression(
        output_path,
        cold_grid_subdive_at_concave_layer_name,
        "type = 'block'",
        output_path,
        cold_grid_block_layer_name,
    )
    cold_grid_block_corner_layer_name = "202_cold_corner"
    extract_by_expression(
        output_path,
        cold_grid_subdive_at_concave_layer_name,
        "type = 'block_corner'",
        output_path,
        cold_grid_block_corner_layer_name,
    )

    cold_grid_block_strip_layer_name = "203_cold_block_strip"
    generate_on_grid_strips_art_sec_loc(
        output_path=output_path,
        blocks_layer_name=cold_grid_block_layer_name,
        roads_layer_name=roads_layer_name,
        part_art_d=part_art_d,
        part_sec_d=part_sec_d,
        part_loc_d=part_loc_d,
        output_layer_name=cold_grid_block_strip_layer_name,
    )

    subdivide_on_grid_corners_convex_layer_name = "204_subdivide_on_grid_corners_convex_layer_name"
    subdivide_on_grid_corners_convex(
        output_path=output_path,
        parts_layer_name=cold_grid_block_strip_layer_name,
        roads_layer_name=roads_layer_name,
        part_sec_d=part_sec_d,
        part_loc_d=part_loc_d,
        plot_sec_w=part_og_w,
        plot_loc_w=part_og_w,
        output_layer_name=subdivide_on_grid_corners_convex_layer_name,
    )

    return subdivide_on_grid_corners_convex_layer_name
