# src/rue_lib/cluster/runner.py
"""Main runner for cluster/partition generation."""

from __future__ import annotations

import os
from pathlib import Path

from osgeo import ogr

from rue_lib.cluster.block_parts import extract_block_parts_from_off_grid
from rue_lib.cluster.config import ClusterConfig
from rue_lib.cluster.frame import extract_frame
from rue_lib.cluster.off_grid import extract_off_grid_inner_layer
from rue_lib.cluster.off_grid_subdivision import (
    extract_off_grid_inner_cluster, extract_off_grid_side_cluster
)
from rue_lib.core.geometry import (
    get_utm_zone_from_layer,
    reproject_layer
)
from rue_lib.streets.operations import extract_by_expression


def generate_clusters(cfg: ClusterConfig) -> Path:
    """
    Generate clusters/partitions from blocks.

    This function takes blocks (output from streets module) and partitions them
    into clusters based on:
    - Proximity to different road types (on-grid clusters)
    - Off-grid cluster dimensions
    - Plot size requirements
    - Public space and amenity allocations

    Args:
        cfg: ClusterConfig with all settings

    Returns:
        Path to output directory containing cluster GeoJSON and summary

    Example:
        >>> config = ClusterConfig(
        ...     site_path="site.geojson",
        ...     roads_path="roads.geojson",
        ...     blocks_path="outputs/streets/all_grids_merged.geojson",
        ...     part_og_w=30.0,
        ...     part_og_d=35.0
        ... )
        >>> output = generate_clusters(config)
    """
    print("=" * 60)
    print("CLUSTER/PARTITION GENERATION")
    print("=" * 60)

    out_dir = Path(cfg.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    output_gpkg = out_dir / "outputs.gpkg"
    output_path = str(output_gpkg)
    if os.path.exists(output_path):
        os.remove(output_path)

    # Step 1: Read input data
    print("Step 1: Determining UTM zone...")
    site_ds = ogr.Open(cfg.blocks_path)
    site_layer = site_ds.GetLayer()
    utm_epsg = get_utm_zone_from_layer(site_layer)
    print(f"  Using UTM EPSG: {utm_epsg}")

    print("Step 2: Reprojecting layers to UTM...")
    roads_layer_name = reproject_layer(
        cfg.roads_path, output_path, utm_epsg, layer_name="01_roads"
    )
    blocks_layer_name = reproject_layer(
        cfg.blocks_path, output_path, utm_epsg, layer_name="02_blocks"
    )

    print("=" * 60)
    print("WARM BLOCK")
    print("=" * 60)
    print("Step 3: Generate inner part of off grid blocks...")
    off_grid_layer_name = "03_off_grid"
    extract_by_expression(
        output_path, blocks_layer_name,
        "type = 'off_grid'",
        output_path,
        off_grid_layer_name
    )
    off_grids_inner_layer_name = extract_off_grid_inner_layer(
        output_path=output_gpkg,
        roads_layer_name=roads_layer_name,
        off_grid_layer_name=off_grid_layer_name,
        output_layer_name="04_off_grids_inner_layer",
        part_art_d=cfg.part_art_d,
        part_sec_d=cfg.part_sec_d,
        part_loc_d=cfg.part_loc_d
    )
    print("Step 4: Generate frame parts of off grid blocks...")
    off_grid_frame_layer_name = "05_off_grid_frame"
    extract_frame(
        output_path=output_gpkg,
        off_grid_layer_name=off_grid_layer_name,
        off_grids_inside_layer_name=off_grids_inner_layer_name,
        output_layer_name=off_grid_frame_layer_name
    )
    print("Step 5: Extract all parts of blocks...")
    off_grid_corners_layer_name = "06_corners"
    off_grid_sides_layer_name = "07_sides"
    extract_block_parts_from_off_grid(
        output_path=output_gpkg,
        off_grid_layer_name=off_grid_layer_name,
        off_grids_inner_layer_name=off_grids_inner_layer_name,
        off_grid_frame_layer_name=off_grid_frame_layer_name,
        angle_threshold=155.0,
        corner_distance=50.0,
        output_corner_layer_name=off_grid_corners_layer_name,
        output_sides_layer_name=off_grid_sides_layer_name
    )
    print("Step 6: Subdiv inner off grid into parts...")
    off_grid_inner_cluster_layer_name = "08_inner_grid_subdiv"
    extract_off_grid_inner_cluster(
        output_path=output_gpkg,
        off_grids_inner_layer_name=off_grids_inner_layer_name,
        part_og_w=cfg.part_og_w,
        part_og_d=cfg.part_og_d,
        output_layer_name=off_grid_inner_cluster_layer_name,
        off_grid_plot_threshold=cfg.off_grid_plot_threshold
    )
    print("Step 7: Subdiv side off grid into parts...")
    off_grid_inner_cluster_layer_name = "09_side_grid_subdiv"
    extract_off_grid_side_cluster(
        output_path=output_gpkg,
        off_grid_sides_layer_name=off_grid_sides_layer_name,
        part_og_w=cfg.part_og_w,
        part_og_d=cfg.part_og_d,
        output_layer_name=off_grid_inner_cluster_layer_name,
        off_grid_plot_threshold=cfg.off_grid_plot_threshold
    )

    print("" + "=" * 60)
    print("CLUSTER GENERATION COMPLETE")
    print("=" * 60)

    return out_dir
