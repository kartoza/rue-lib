"""Main runner for grid-based parcel generation."""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd

from rue_lib.core.geometry import merge_connected_lines
from rue_lib.core.helpers import remove_layer_from_gpkg
from rue_lib.geo import to_metric_crs
from rue_lib.site.config import SiteConfig
from rue_lib.site.financial import FinancialSite
from rue_lib.site.io import read_roads, read_site, save_geojson
from rue_lib.site.roads import buffer_roads
from rue_lib.streets.operations import extract_by_expression
from rue_lib.streets.runner_utils import subtract_layer


def generate_parcels(cfg: SiteConfig) -> Path:
    """
    Generate parcels.

    Args:
        cfg: SiteConfig with all settings

    Returns:
        Path to output GeoJSON file containing parcel polygons

    Example:
        >>> config = SiteConfig(
        ...     site_path="site.geojson",
        ...     roads_path="roads.geojson",
        ...     road_arterial_width_m=20,
        ...     road_secondary_width_m=15
        ... )
        >>> output = generate_parcels(config)
    """
    out_dir = Path(cfg.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gpkg_path = out_dir / "outputs.gpkg"
    gpkg_path = str(gpkg_path)

    site = read_site(cfg.site_path)
    roads = read_roads(cfg.roads_path)

    roads_layer_name = "roads"

    if site.crs and site.crs.is_projected:
        site_m = site
    else:
        site_m = to_metric_crs(site)

    if roads.crs and roads.crs.is_projected:
        roads_m = roads
    else:
        roads_m = to_metric_crs(roads)

    site_m.to_file(gpkg_path, layer="site", driver="GPKG")
    roads_m.to_file(gpkg_path, layer="roads", driver="GPKG")

    roads_buffer_layer = buffer_roads(
        roads_m, cfg.road_arterial_width_m, cfg.road_secondary_width_m
    )
    roads_buf_m = buffer_roads(
        roads_m,
        cfg.road_arterial_width_m - cfg.road_local_width_m,
        cfg.road_secondary_width_m - cfg.road_local_width_m,
    )
    extract_by_expression(
        gpkg_path, roads_layer_name, "type = 'road_art'", gpkg_path, "02_arterial_roads"
    )
    extract_by_expression(
        gpkg_path, roads_layer_name, "type = 'road_sec'", gpkg_path, "03_secondary_roads"
    )

    # Merge connected lines in both arterial and secondary roads before subtraction
    merge_connected_lines(gpkg_path, "02_arterial_roads")
    merge_connected_lines(gpkg_path, "03_secondary_roads")

    if not roads_buf_m.empty:
        roads_buf_m.to_file(gpkg_path, layer="roads_buffered", driver="GPKG")
        subtract_layer(
            input_gpkg=str(gpkg_path),
            base_layer_name="site",
            erase_layer_name="02_arterial_roads",
            output_gpkg=str(gpkg_path),
            output_layer_name="parcels",
            buffer_distance=(cfg.road_arterial_width_m - cfg.road_local_width_m) / 2,
            buffer_cap_style="flat",
        )
        subtract_layer(
            input_gpkg=str(gpkg_path),
            base_layer_name="parcels",
            erase_layer_name="03_secondary_roads",
            output_gpkg=str(gpkg_path),
            output_layer_name="parcels",
            buffer_distance=(cfg.road_secondary_width_m - cfg.road_local_width_m) / 2,
            buffer_cap_style="flat",
        )
    else:
        site_m.to_file(gpkg_path, layer="parcels", driver="GPKG")

    parcels_m = gpd.read_file(gpkg_path, layer="parcels")

    parcels_exploded = parcels_m.explode(index_parts=False, ignore_index=True)

    parcels_exploded["area_m2"] = parcels_exploded.geometry.area
    parcels_exploded["id"] = range(1, len(parcels_exploded) + 1)

    parcels_exploded.to_file(gpkg_path, layer="sites", driver="GPKG")

    remove_layer_from_gpkg(gpkg_path, "parcels")

    parcels_final = (
        parcels_exploded.to_crs(4326) if not parcels_exploded.empty else parcels_exploded
    )

    out_geojson = out_dir / "outputs.geojson"
    save_geojson(parcels_final, out_geojson)

    roads_buffer_layer = (
        roads_buffer_layer.to_crs(4326) if not roads_buffer_layer.empty else roads_buffer_layer
    )

    out_geojson = out_dir / "roads.geojson"
    save_geojson(roads_buffer_layer, out_geojson)

    print("Generating financial data")
    FinancialSite(config=cfg)

    return out_geojson
