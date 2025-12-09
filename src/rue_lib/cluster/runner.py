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
from .runner_warm import generate_warm


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
    utm_epsg = 3857
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
    generate_warm(cfg, output_gpkg, input_blocks_layer_name, roads_layer_name)

    print("" + "=" * 60)
    print("CLUSTER GENERATION COMPLETE")
    print("=" * 60)

    return out_dir
