# src/rue_lib/cluster/runner.py
"""Main runner for cluster/partition generation."""

from __future__ import annotations

import os
from pathlib import Path

from osgeo import ogr

from rue_lib.cluster.config import ClusterConfig
from rue_lib.core.geometry import (
    get_utm_zone_from_layer,
    reproject_layer
)
from rue_lib.streets.operations import (
    extract_by_geometry_type
)
from .runner_cold import generate_cold
from .runner_warm import generate_warm
from ..core import merge_gpkg_layers
from ..core.roads import extract_roads_buffer


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
    """
    print("==============================================================")
    print("CLUSTER/PARTITION GENERATION")
    print("==============================================================")

    out_dir = Path(cfg.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    output_gpkg = out_dir / "outputs.gpkg"
    output_path = str(output_gpkg)
    if os.path.exists(output_path):
        os.remove(output_path)

    # Step 1: Read input data
    print("Step 1: Determining UTM zone...")
    site_ds = ogr.Open(cfg.input_path)
    site_layer = site_ds.GetLayer()
    utm_epsg = get_utm_zone_from_layer(site_layer)
    print(f"  Using UTM EPSG: {utm_epsg}")

    print("Step 2: Reprojecting layers to UTM...")
    input_layer_name = reproject_layer(
        cfg.input_path, output_path, utm_epsg, layer_name="000_input"
    )
    roads_layer_name = reproject_layer(
        cfg.roads_path, output_path, utm_epsg, layer_name="001_roads"
    )
    input_blocks_layer_name = "002_input_blocks"
    extract_by_geometry_type(
        output_path, input_layer_name,
        ["POLYGON", "MULTIPOLYGON"],
        output_path,
        input_blocks_layer_name
    )

    input_roads_layer_name = "002_input_roads"
    extract_by_geometry_type(
        output_path, input_layer_name,
        ["LINESTRING", "MULTILINESTRING"],
        output_path,
        input_roads_layer_name
    )

    input_roads_buffer_layer_name = "002_input_roads_buffer"
    extract_roads_buffer(
        input_path=output_path,
        input_layer_name=input_roads_layer_name,
        output_path=output_path,
        output_layer_name=input_roads_buffer_layer_name,
        # TODO:
        #  We add half of local as currently it contains locals
        road_arterial_width_m=(
                cfg.road_arterial_width_m + cfg.road_local_width_m
        ),
        road_secondary_width_m=(
                cfg.road_secondary_width_m + cfg.road_local_width_m
        ),
        road_local_width_m=cfg.road_local_width_m
    )

    # # Warm block generation
    warm_final_layer_name = generate_warm(
        cfg, output_gpkg, input_blocks_layer_name,
        input_roads_buffer_layer_name
    )

    # Cold block generation
    cold_final_layer_name = generate_cold(
        cfg, output_gpkg, input_blocks_layer_name,
        input_roads_buffer_layer_name
    )
    print("Final step: Merge all")
    merge_gpkg_layers(
        gpkg_path=output_gpkg,
        layer_names=[
            warm_final_layer_name,
            cold_final_layer_name
        ],
        output_layer_name="300_final",
    )
    print("" + "=" * 60)
    print("CLUSTER GENERATION COMPLETE")
    print("=" * 60)

    return out_dir
