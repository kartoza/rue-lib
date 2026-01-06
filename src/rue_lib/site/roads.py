"""Road processing for Step 1 - Site."""

import geopandas as gpd
import pandas as pd
from shapely import wkt


def buffer_geometry(geom_wkt: str, distance: float, cap_style: int = 2, join_style: int = 2) -> str:
    """
    Buffer a geometry

    Args:
        geom_wkt: Geometry as WKT string
        distance: Buffer distance
        cap_style: 1=round, 2=flat, 3=square
        join_style: 1=round, 2=mitre, 3=bevel
    Returns:
        Buffered geometry as WKT string
    """
    try:
        geom = wkt.loads(geom_wkt)
        buffered = geom.buffer(distance, cap_style=cap_style, join_style=join_style)
        return buffered.wkt
    except Exception:
        return geom_wkt


def buffer_roads(
    roads_m: gpd.GeoDataFrame,
    road_arterial_width_m: float = 16,
    road_secondary_width_m: float = 10,
    road_type_key: str = "type",
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
    if road_type_key not in roads_m.columns:
        roads_m[road_type_key] = None

    art = roads_m[roads_m[road_type_key] == "road_art"].copy()
    sec = roads_m[roads_m[road_type_key].isin(["road_sec", "road_sec_new"])].copy()

    parts = []

    if not art.empty:
        art["geometry"] = art.geometry.apply(
            lambda g: wkt.loads(
                buffer_geometry(g.wkt, road_arterial_width_m / 2, cap_style=2, join_style=2)
            )
        )
        parts.append(art)

    if not sec.empty:
        sec["geometry"] = sec.geometry.apply(
            lambda g: wkt.loads(
                buffer_geometry(g.wkt, road_secondary_width_m / 2, cap_style=2, join_style=2)
            )
        )
        parts.append(sec)

    if parts:
        return gpd.GeoDataFrame(pd.concat(parts, ignore_index=True), crs=roads_m.crs)

    return gpd.GeoDataFrame(columns=list(roads_m.columns), geometry=[], crs=roads_m.crs)
