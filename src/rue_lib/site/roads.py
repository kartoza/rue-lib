"""Road processing for Step 1 - Site."""

import geopandas as gpd
import pandas as pd
from shapely.validation import make_valid


def buffer_roads(roads_m: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Create road corridors as polygons.

    Arterial roads: 8m half-buffer (16m total)
    Secondary roads: 5m half-buffer (10m total)

    Args:
        roads_m: Roads in metric CRS

    Returns:
        Buffered road polygons
    """
    if "type" not in roads_m.columns:
        roads_m["type"] = None

    art = roads_m[roads_m["type"] == "road_art"].copy()
    sec = roads_m[roads_m["type"].isin(["road_sec", "road_sec_new"])].copy()

    parts = []

    if not art.empty:
        art["geometry"] = art.geometry.buffer(8.0, cap_style=2, join_style=2)
        parts.append(art)

    if not sec.empty:
        sec["geometry"] = sec.geometry.buffer(5.0, cap_style=2, join_style=2)
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
