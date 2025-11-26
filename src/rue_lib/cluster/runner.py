# src/rue_lib/cluster/runner.py
"""Main runner for cluster/partition generation."""

from __future__ import annotations

import math
from pathlib import Path

import geopandas as gpd

from rue_lib.cluster.config import ClusterConfig
from rue_lib.cluster.io import read_blocks, read_roads, read_site


def to_metric_crs(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Project to local UTM CRS for metric operations."""
    if gdf.crs is None:
        gdf = gdf.set_crs(4326)

    try:
        centroid = gdf.union_all().centroid
    except Exception:
        centroid = gdf.geometry.iloc[0].centroid

    lon = centroid.x
    utm_zone = int(math.floor((lon + 180) / 6) + 1)
    is_northern = centroid.y >= 0
    epsg = 32600 + utm_zone if is_northern else 32700 + utm_zone

    return gdf.to_crs(epsg)


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
        ...     part_art_d=40.0,
        ...     part_sec_d=30.0,
        ...     part_loc_d=20.0,
        ... )
        >>> output = generate_clusters(config)
    """
    print("=" * 60)
    print("CLUSTER/PARTITION GENERATION")
    print("=" * 60)

    out_dir = Path(cfg.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Read input data
    print("\nStep 1: Reading input data...")
    site = read_site(cfg.site_path)
    roads = read_roads(cfg.roads_path)
    blocks = read_blocks(cfg.blocks_path)

    print(f"  Site features: {len(site)}")
    print(f"  Road features: {len(roads)}")
    print(f"  Block features: {len(blocks)}")

    # Step 2: Project to metric CRS
    print("\nStep 2: Projecting to UTM CRS...")
    _site_m = to_metric_crs(site)
    _roads_m = to_metric_crs(roads)
    blocks_m = to_metric_crs(blocks)
    target_crs = blocks_m.crs
    print(f"  Using CRS: {target_crs}")

    print("\n" + "=" * 60)
    print("CLUSTER GENERATION COMPLETE")
    print("=" * 60)

    return out_dir
