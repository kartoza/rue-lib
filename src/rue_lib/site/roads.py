"""Road processing for Step 1 - Site."""

import geopandas as gpd
import pandas as pd
from shapely.validation import make_valid


def buffer_roads(
    roads_m: gpd.GeoDataFrame, road_arterial_width_m: float = 16, road_secondary_width_m: float = 10
) -> gpd.GeoDataFrame:
    """
    Create road corridors as polygons.

    Args:
        roads_m: Roads in metric CRS
        road_arterial_width_m: Width of arterial roads (meters)
        road_secondary_width_m: Width of secondary roads (meters)

    Returns:
        Buffered road polygons
    """
    if "type" not in roads_m.columns:
        roads_m["type"] = None

    art = roads_m[roads_m["type"] == "road_art"].copy()
    sec = roads_m[roads_m["type"].isin(["road_sec", "road_sec_new"])].copy()

    parts = []

    if not art.empty:
        art["geometry"] = art.geometry.buffer(road_arterial_width_m / 2, cap_style=2, join_style=2)
        parts.append(art)

    if not sec.empty:
        sec["geometry"] = sec.geometry.buffer(road_secondary_width_m / 2, cap_style=2, join_style=2)
        parts.append(sec)

    if parts:
        return gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=roads_m.crs)

    return gpd.GeoDataFrame(columns=list(roads_m.columns), geometry=[], crs=roads_m.crs)


def subdivide_site_by_roads(
    site_m: gpd.GeoDataFrame, roads_buf_m: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """
    Boolean difference: site minus road buffers.

    Args:
        site_m: Site geometry in metric CRS
        roads_buf_m: Buffered road polygons

    Returns:
        Site split by roads
    """
    site_geom = make_valid(site_m.union_all()).buffer(0)

    if not roads_buf_m.empty:
        roads_geom = make_valid(roads_buf_m.union_all()).buffer(0)
        diff = site_geom.difference(roads_geom)
    else:
        diff = site_geom

    if diff.is_empty:
        return gpd.GeoDataFrame(geometry=[], crs=site_m.crs)

    out = gpd.GeoDataFrame(geometry=[diff], crs=site_m.crs).explode(
        index_parts=False, ignore_index=True
    )

    out["geometry"] = out.buffer(0)

    return out
