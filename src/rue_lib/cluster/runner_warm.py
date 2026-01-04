# src/rue_lib/cluster/runner_warm.py
"""Generate warm blocks with off-grid subdivision and partitioning."""

from __future__ import annotations

from pathlib import Path

from rue_lib.cluster.art_sec_no_offgrid import generate_art_sec_parts_no_offgrid
from rue_lib.cluster.block_parts import extract_block_parts_from_off_grid
from rue_lib.cluster.config import ClusterConfig
from rue_lib.cluster.frame import extract_frame
from rue_lib.cluster.off_grid import extract_off_grid_inner_layer
from rue_lib.cluster.off_grid_subdivision import extract_off_grid_cluster
from rue_lib.core.definitions import BlockTypes
from rue_lib.core.helpers import merge_gpkg_layers
from rue_lib.streets.operations import extract_by_expression


def generate_warm(
    cfg: ClusterConfig, output_gpkg: Path, input_blocks_layer_name: str, roads_layer_name: str
) -> str:
    """
    Generate warm blocks with off-grid subdivision and partitioning.

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

    part_og_w = cfg.off_grid_cluster_width
    part_og_d = cfg.off_grid_cluster_depth

    output_path = str(output_gpkg)
    print("==============================================================")
    print("WARM BLOCK")
    print("==============================================================")
    print("Step 1: Generate inner part of off grid blocks...")
    warm_grid_layer_name = "100_warm_grid"
    extract_by_expression(
        output_path,
        input_blocks_layer_name,
        (
            f"type = '{BlockTypes.ON_GRID_ART}' OR "
            f"type = '{BlockTypes.ON_GRID_SEC}' OR "
            f"type = '{BlockTypes.OFF_GRID}'"
        ),
        output_path,
        warm_grid_layer_name,
    )
    warm_grid_off_grid_layer_name = "100_warm_grid_off_grid"
    extract_by_expression(
        output_path,
        warm_grid_layer_name,
        "type = 'off_grid'",
        output_path,
        warm_grid_off_grid_layer_name,
    )
    warm_grid_on_grid_layer_name = "100_warm_grid_on_grid"
    extract_by_expression(
        output_path,
        warm_grid_layer_name,
        "type != 'off_grid'",
        output_path,
        warm_grid_on_grid_layer_name,
    )
    off_grids_inner_layer_name = "101_off_grids_inner_layer"
    off_grids_inner_layer_name = extract_off_grid_inner_layer(
        output_path=output_gpkg,
        roads_layer_name=roads_layer_name,
        off_grid_layer_name=warm_grid_off_grid_layer_name,
        output_layer_name=off_grids_inner_layer_name,
        part_art_d=part_art_d,
        part_sec_d=part_sec_d,
        part_loc_d=part_loc_d,
    )
    if not off_grids_inner_layer_name:
        print("   No off grid found, as site is too small")
    else:
        print("Step 2: Generate frame parts of off grid blocks...")
        off_grid_frame_layer_name = "102_off_grid_frame"
        extract_frame(
            output_path=output_gpkg,
            off_grid_layer_name=warm_grid_layer_name,
            off_grids_inside_layer_name=off_grids_inner_layer_name,
            output_layer_name=off_grid_frame_layer_name,
        )
        print("Step 3: Extract all parts of blocks...")
        off_grid_corners_layer_name = "103_corners"
        off_grid_sides_layer_name = "104_sides"
        off_grid_off_grid_layer_name = "105_off_grid"
        extract_block_parts_from_off_grid(
            output_path=output_gpkg,
            roads_layer_name=roads_layer_name,
            warm_grid_layer_name=warm_grid_layer_name,
            off_grids_inner_layer_name=off_grids_inner_layer_name,
            off_grid_frame_layer_name=off_grid_frame_layer_name,
            angle_threshold=155.0,
            corner_distance=50.0,
            output_corner_layer_name=off_grid_corners_layer_name,
            output_sides_layer_name=off_grid_sides_layer_name,
            output_off_grid_layer_name=off_grid_off_grid_layer_name,
            part_art_d=part_art_d,
            part_sec_d=part_sec_d,
            part_loc_d=part_loc_d,
        )
        print("Step 4: Get blocks that is not in off grid and merge them with off grid blocks...")
        merge_gpkg_layers(
            gpkg_path=output_gpkg,
            layer_names=[
                off_grid_corners_layer_name,
                off_grid_sides_layer_name,
                off_grid_off_grid_layer_name,
                warm_grid_on_grid_layer_name,
            ],
            output_layer_name="107_generated_part_with_off_grid_checkpoint",
        )
        print("Step 5: Subdiv inner off grid into parts...")
        off_grid_inner_cluster_layer_name = "108_inner_grid_subdiv"
        extract_off_grid_cluster(
            output_path=output_gpkg,
            off_grids_layer_name=off_grid_off_grid_layer_name,
            part_og_w=part_og_w,
            part_og_d=part_og_d,
            output_layer_name=off_grid_inner_cluster_layer_name,
            min_plot_area=part_og_w * part_og_d * cfg.off_grid_plot_threshold,
        )
        print("Step 6: Subdiv side off grid into parts...")
        off_grid_side_cluster_layer_name = "109_side_grid_subdiv"
        extract_off_grid_cluster(
            output_path=output_gpkg,
            off_grids_layer_name=off_grid_sides_layer_name,
            part_og_w=part_og_w,
            part_og_d=part_og_d,
            output_layer_name=off_grid_side_cluster_layer_name,
            min_plot_area=part_loc_d * part_og_w * cfg.off_grid_plot_threshold,
        )

    print("Step 8: generate art sec parts no offgrid ")
    generate_art_sec_parts_no_offgrid_layer_name = generate_art_sec_parts_no_offgrid(
        output_path=output_path,
        blocks_layer_name=warm_grid_on_grid_layer_name,
        roads_layer_name="002_input_roads",
        roads_buffered_layer_name=roads_layer_name,
        part_art_d=part_art_d,
        art_road_width_m=cfg.road_arterial_width_m / 2,
        part_sec_d=part_sec_d,
        sec_road_width_m=cfg.road_secondary_width_m / 2,
        part_loc_d=part_loc_d,
        loc_road_width_m=cfg.road_local_width_m / 2,
        output_layer_name="110_generate_art_sec_parts_no_offgrid",
    )

    if off_grids_inner_layer_name:
        print("Step 7: Merge subdiv grids")
        return merge_gpkg_layers(
            gpkg_path=output_gpkg,
            layer_names=[
                off_grid_inner_cluster_layer_name,
                off_grid_side_cluster_layer_name,
                off_grid_corners_layer_name,
                generate_art_sec_parts_no_offgrid_layer_name,
            ],
            output_layer_name="111_warm_block_final",
        )
    else:
        print("Step 7: Merge subdiv grids")
        return merge_gpkg_layers(
            gpkg_path=output_gpkg,
            layer_names=[
                generate_art_sec_parts_no_offgrid_layer_name,
            ],
            output_layer_name="111_warm_block_final",
        )
