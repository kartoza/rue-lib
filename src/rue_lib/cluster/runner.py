# src/rue_lib/cluster/runner.py
"""Main runner for cluster/partition generation."""

from __future__ import annotations

import os
from pathlib import Path

import geopandas as gpd
from osgeo import ogr
from shapely.validation import make_valid

from rue_lib.cluster.config import ClusterConfig
from rue_lib.core.geometry import get_utm_zone_from_layer, reproject_layer
from rue_lib.streets.operations import (
    export_layer_to_geojson,
    extract_by_expression,
    extract_by_geometry_type,
)

from ..core import merge_gpkg_layers
from ..core.definitions import BlockTypes
from ..core.roads import extract_roads_buffer
from .helpers import assign_cluster_type
from .runner_cold import generate_cold
from .runner_warm import generate_warm


def validate_and_clean_layer(gpkg_path: Path, layer_name: str) -> None:
    """
    Validate and clean geometries in a layer to prevent issues with viewers.

    This function:
    - Removes features with None/empty geometries
    - Fixes invalid geometries using make_valid
    - Ensures all geometries are valid Polygons or MultiPolygons

    Args:
        gpkg_path: Path to the GeoPackage file
        layer_name: Name of the layer to validate and clean
    """
    print(f"Validating and cleaning layer: {layer_name}...")

    gdf = gpd.read_file(gpkg_path, layer=layer_name)
    initial_count = len(gdf)

    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()]

    invalid_mask = ~gdf.geometry.is_valid
    invalid_count = invalid_mask.sum()

    if invalid_count > 0:
        print(f"  Found {invalid_count} invalid geometries, fixing...")
        gdf.loc[invalid_mask, "geometry"] = gdf.loc[invalid_mask, "geometry"].apply(make_valid)

    valid_geom_types = ["Polygon", "MultiPolygon"]
    gdf = gdf[gdf.geometry.geom_type.isin(valid_geom_types)]

    final_count = len(gdf)
    removed_count = initial_count - final_count

    if removed_count > 0:
        print(f"  Removed {removed_count} invalid/empty features")

    print(f"  Final feature count: {final_count}")
    gdf.to_file(gpkg_path, layer=layer_name, driver="GPKG")
    print(f"  Layer {layer_name} cleaned and saved")


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
    reproject_layer(cfg.roads_path, output_path, utm_epsg, layer_name="001_roads")
    input_blocks_layer_name = "002_input_blocks"
    extract_by_geometry_type(
        output_path,
        input_layer_name,
        ["POLYGON", "MULTIPOLYGON"],
        output_path,
        input_blocks_layer_name,
    )
    extract_by_expression(
        output_path,
        input_blocks_layer_name,
        (
            f"type = '{BlockTypes.COLD_GRID}' or "
            f"type = '{BlockTypes.ON_GRID_ART}' or "
            f"type = '{BlockTypes.ON_GRID_SEC}' or "
            f"type = '{BlockTypes.OFF_GRID}'"
        ),
        output_path,
        input_blocks_layer_name,
    )

    input_roads_layer_name = "002_input_roads"
    extract_by_geometry_type(
        output_path,
        input_layer_name,
        ["LINESTRING", "MULTILINESTRING"],
        output_path,
        input_roads_layer_name,
    )

    input_roads_buffer_layer_name = "002_input_roads_buffer"
    extract_roads_buffer(
        input_path=output_path,
        input_layer_name=input_roads_layer_name,
        output_path=output_path,
        output_layer_name=input_roads_buffer_layer_name,
        # TODO:
        #  We add half of local as currently it contains locals
        road_arterial_width_m=cfg.road_arterial_width_m,
        road_secondary_width_m=cfg.road_secondary_width_m,
        road_local_width_m=cfg.road_local_width_m,
    )

    # Warm block generation
    warm_final_layer_name = generate_warm(
        cfg, output_gpkg, input_blocks_layer_name, input_roads_buffer_layer_name
    )

    # Cold block generation
    cold_final_layer_name = generate_cold(
        cfg, output_gpkg, input_blocks_layer_name, input_roads_buffer_layer_name
    )

    input_roads_local_buffer_layer_name = "002_input_roads_buffer"
    extract_by_expression(
        output_path,
        input_roads_buffer_layer_name,
        "type = 'road_local'",
        output_path,
        input_roads_local_buffer_layer_name,
    )
    print("Final step: Merge all")
    final_layer_name = "300_final"
    merge_gpkg_layers(
        gpkg_path=output_gpkg,
        layer_names=[
            warm_final_layer_name,
            cold_final_layer_name,
            input_roads_local_buffer_layer_name,
        ],
        output_layer_name=final_layer_name,
    )
    assign_cluster_type(
        gpkg_path=output_gpkg,
        block_layer_name=input_blocks_layer_name,
        final_layer_name=final_layer_name,
    )

    validate_and_clean_layer(output_gpkg, final_layer_name)

    print("Exporting merged grids to GeoJSON...")
    output_geojson = out_dir / "outputs.geojson"
    export_layer_to_geojson(
        str(output_gpkg),
        final_layer_name,
        str(output_geojson),
    )

    print("" + "=" * 60)
    print("CLUSTER GENERATION COMPLETE")
    print("=" * 60)

    return out_dir
