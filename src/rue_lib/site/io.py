"""I/O operations for Step 1 - Site."""

from pathlib import Path

import geopandas as gpd


def read_site(path: str | Path) -> gpd.GeoDataFrame:
    """Read site polygon(s) from file."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Site file not found: {p}")

    gdf = gpd.read_file(p)
    if gdf.empty:
        raise ValueError("Provided site file is empty.")

    return gdf


def read_roads(path: str | Path) -> gpd.GeoDataFrame:
    """Read roads from file and standardize column names."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Roads file not found: {p}")

    gdf = gpd.read_file(p)
    if gdf.empty:
        raise ValueError("Provided roads file is empty.")

    # Standardize column names
    if "type" not in gdf.columns and "road_type" in gdf.columns:
        gdf = gdf.rename(columns={"road_type": "type"})
    if "share" not in gdf.columns and "road_pcent" in gdf.columns:
        gdf = gdf.rename(columns={"road_pcent": "share"})

    return gdf


def save_geojson(gdf: gpd.GeoDataFrame, out_path: Path) -> None:
    """Save GeoDataFrame to GeoJSON."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(out_path, driver="GeoJSON")
