"""Main runner for grid-based parcel generation."""

from __future__ import annotations

from pathlib import Path

from shapely.validation import make_valid

from rue_lib.geo import to_metric_crs
from rue_lib.site.config import SiteConfig
from rue_lib.site.io import read_roads, read_site, save_geojson
from rue_lib.site.roads import buffer_roads


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
        ...     road_arterial_width_m=20,
        ...     road_secondary_width_m=15
        ... )
        >>> output = generate_parcels(config)
    """
    out_dir = Path(cfg.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    site = read_site(cfg.site_path)
    roads = read_roads(cfg.roads_path)
    site_m = to_metric_crs(site)
    roads_m = to_metric_crs(roads)

    roads_buf_m = buffer_roads(roads_m, cfg.road_arterial_width_m, cfg.road_secondary_width_m)

    if not roads_buf_m.empty:
        roads_geom = make_valid(roads_buf_m.union_all()).buffer(0)
        site_m["geometry"] = site_m.geometry.apply(
            lambda g: make_valid(g).buffer(0).difference(roads_geom)
        )
        site_m["geometry"] = site_m.buffer(0)
        site_m["area_m2"] = site_m.geometry.area

    parcels = site_m.to_crs(4326) if not site_m.empty else site_m

    out_path = out_dir / "parcels.geojson"
    save_geojson(parcels, out_path)

    return out_path
