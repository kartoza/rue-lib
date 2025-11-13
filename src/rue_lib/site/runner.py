"""Main runner for grid-based parcel generation."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import geopandas as gpd
from shapely.validation import make_valid

from rue_lib.site.generator import create_parcel_grid
from rue_lib.site.io import read_roads, read_site, save_geojson
from rue_lib.site.roads import buffer_roads


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


@dataclass
class SiteConfig:
    """Configuration for grid-based parcel generation."""

    site_path: str
    roads_path: str
    output_dir: str = "outputs/parcels"
    rows: int = 4
    cols: int = 4
    pad_m: float = 50.0
    min_parcel_area_m2: float = 5.0
    subtract_roads: bool = False


def generate_parcels(cfg: SiteConfig) -> Path:
    """
    Generate parcels using a grid-based approach.

    This creates ownership parcels by overlaying a rectangular grid
    on the site and optionally subtracting road corridors.

    Args:
        cfg: SiteConfig with all settings

    Returns:
        Path to output GeoJSON file containing parcel polygons

    Example:
        >>> config = SiteConfig(
        ...     site_path="site.geojson",
        ...     roads_path="roads.geojson",
        ...     rows=5,
        ...     cols=5
        ... )
        >>> output = generate_parcels(config)
    """
    out_dir = Path(cfg.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Read and project to metric CRS
    site = read_site(cfg.site_path)
    roads = read_roads(cfg.roads_path)
    site_m = to_metric_crs(site)
    roads_m = to_metric_crs(roads)

    # Buffer roads for optional subtraction
    roads_buf_m = buffer_roads(roads_m)

    # Create parcels as grid intersection with site
    parcels_m = create_parcel_grid(site_m, rows=cfg.rows, cols=cfg.cols, pad=cfg.pad_m)

    # Optionally subtract roads after grid creation
    if cfg.subtract_roads and not parcels_m.empty and not roads_buf_m.empty:
        roads_geom = make_valid(roads_buf_m.union_all()).buffer(0)
        parcels_m["geometry"] = parcels_m.geometry.apply(
            lambda g: make_valid(g).buffer(0).difference(roads_geom)
        )
        parcels_m["geometry"] = parcels_m.buffer(0)
        parcels_m["area_m2"] = parcels_m.geometry.area

    # Filter micro-slivers
    if not parcels_m.empty and cfg.min_parcel_area_m2 is not None:
        parcels_m = parcels_m[parcels_m["area_m2"] >= float(cfg.min_parcel_area_m2)].copy()

    # Back to WGS84
    parcels = parcels_m.to_crs(4326) if not parcels_m.empty else parcels_m

    out_path = out_dir / "parcels.geojson"
    save_geojson(parcels, out_path)

    return out_path
