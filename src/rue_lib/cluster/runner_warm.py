# src/rue_lib/cluster/runner.py
"""Main runner for cluster/partition generation."""

from __future__ import annotations

from pathlib import Path

from rue_lib.cluster.art_sec_no_offgrid import (
    generate_art_sec_parts_no_offgrid
)
from rue_lib.cluster.block_parts import extract_block_parts_from_off_grid
from rue_lib.cluster.config import ClusterConfig
from rue_lib.cluster.frame import extract_frame
from rue_lib.cluster.off_grid import extract_off_grid_inner_layer
from rue_lib.cluster.off_grid_subdivision import extract_off_grid_cluster
from rue_lib.core.helpers import merge_gpkg_layers
from rue_lib.streets.operations import (
    extract_by_expression
)


def generate_warm(
        cfg: ClusterConfig, output_gpkg: Path, input_blocks_layer_name: str,
        roads_layer_name: str
):
    """
    Generate warm blocks from parcels.
    """
    output_path = str(output_gpkg)
    print("==============================================================")
    print("WARM BLOCK")
    print("==============================================================")
    print("Step 1: Generate inner part of off grid blocks...")
    warm_grid_layer_name = "103_warm_grid"
    extract_by_expression(
        output_path, input_blocks_layer_name,
        "type = 'on_grid_art' OR type = 'on_grid_sec' OR type = 'off_grid'",
        output_path,
        warm_grid_layer_name
    )
    off_grids_inner_layer_name = "104_off_grids_inner_layer"
    extract_off_grid_inner_layer(
        output_path=output_gpkg,
        roads_layer_name=roads_layer_name,
        off_grid_layer_name=warm_grid_layer_name,
        output_layer_name=off_grids_inner_layer_name,
        part_art_d=cfg.part_art_d,
        part_sec_d=cfg.part_sec_d,
        part_loc_d=cfg.part_loc_d
    )
    print("Step 2: Generate frame parts of off grid blocks...")
    off_grid_frame_layer_name = "105_off_grid_frame"
    extract_frame(
        output_path=output_gpkg,
        off_grid_layer_name=warm_grid_layer_name,
        off_grids_inside_layer_name=off_grids_inner_layer_name,
        output_layer_name=off_grid_frame_layer_name
    )
    print("Step 3: Extract all parts of blocks...")
    off_grid_corners_layer_name = "106_corners"
    off_grid_sides_layer_name = "107_sides"
    off_grid_off_grid_layer_name = "108_off_grid"
    output_not_generated_block_layer_name = "109_not_generated_block"
    extract_block_parts_from_off_grid(
        output_path=output_gpkg,
        warm_grid_layer_name=warm_grid_layer_name,
        off_grids_inner_layer_name=off_grids_inner_layer_name,
        off_grid_frame_layer_name=off_grid_frame_layer_name,
        output_not_generated_block_layer_name=output_not_generated_block_layer_name,
        angle_threshold=155.0,
        corner_distance=50.0,
        output_corner_layer_name=off_grid_corners_layer_name,
        output_sides_layer_name=off_grid_sides_layer_name,
        output_off_grid_layer_name=off_grid_off_grid_layer_name
    )
    print(
        "Step 4: Get blocks that is not in off grid and merge them with off grid blocks...")
    merge_gpkg_layers(
        gpkg_path=output_gpkg,
        layer_names=[
            off_grid_corners_layer_name, off_grid_sides_layer_name,
            off_grid_off_grid_layer_name, output_not_generated_block_layer_name
        ],
        output_layer_name="110_generated_part_with_off_grid_checkpoint",
    )
    print("Step 5: Subdiv inner off grid into parts...")
    off_grid_inner_cluster_layer_name = "111_inner_grid_subdiv"
    extract_off_grid_cluster(
        output_path=output_gpkg,
        off_grids_layer_name=off_grid_off_grid_layer_name,
        part_og_w=cfg.part_og_w,
        part_og_d=cfg.part_og_d,
        output_layer_name=off_grid_inner_cluster_layer_name,
        off_grid_plot_threshold=cfg.off_grid_plot_threshold,
        min_plot_area=cfg.part_og_w * cfg.part_og_d * cfg.off_grid_plot_threshold,
    )
    print("Step 6: Subdiv side off grid into parts...")
    off_grid_side_cluster_layer_name = "112_side_grid_subdiv"
    extract_off_grid_cluster(
        output_path=output_gpkg,
        off_grids_layer_name=off_grid_sides_layer_name,
        part_og_w=cfg.part_og_w,
        part_og_d=cfg.part_og_d,
        output_layer_name=off_grid_side_cluster_layer_name,
        off_grid_plot_threshold=cfg.off_grid_plot_threshold,
        min_plot_area=cfg.part_loc_d * cfg.plot_loc_w * cfg.off_grid_plot_threshold,
    )
    print("Step 7: Merge subdiv grids")
    merge_gpkg_layers(
        gpkg_path=output_gpkg,
        layer_names=[
            off_grid_inner_cluster_layer_name,
            off_grid_side_cluster_layer_name,
            off_grid_corners_layer_name
        ],
        output_layer_name="113_subdiv_into_parts_checkpoint",
    )

    print("Step 8: generate art sec parts no offgrid ")
    generate_art_sec_parts_no_offgrid(
        output_path=output_path,
        blocks_layer_name=output_not_generated_block_layer_name,
        roads_layer_name=roads_layer_name,
        part_art_d=cfg.part_art_d,
        part_sec_d=cfg.part_sec_d,
        part_loc_d=cfg.part_loc_d,
        output_layer_name="114_generate_art_sec_parts_no_offgrid",
    )
