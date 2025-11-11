# src/rue_lib/streets/io.py
import json
from pathlib import Path

import geopandas as gpd


def read_parcels(path: str | Path) -> gpd.GeoDataFrame:
    """Read parcels from GeoJSON file."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Parcels file not found: {p}")

    gdf = gpd.read_file(p)
    if gdf.empty:
        raise ValueError("Parcels file is empty.")

    return gdf


def save_blocks(blocks, output_path):
    """Save blocks to a GeoJSON file.

    Accepts either a GeoDataFrame or an iterable of GeoJSON Feature dicts
    (as returned by `generate_blocks`). Handles empty inputs gracefully.
    """
    from pathlib import Path

    import geopandas as gpd

    # Already a GeoDataFrame
    if hasattr(blocks, "to_file"):
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        blocks.to_file(output_path, driver="GeoJSON")
        return

    # Empty iterable ⇒ write an empty FeatureCollection
    if not blocks:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump({"type": "FeatureCollection", "features": []}, f)
        return

    # Normal path
    try:
        gdf = gpd.GeoDataFrame.from_features(blocks)
        # Ensure geometry+CRS exist
        if "geometry" not in gdf.columns:
            raise ValueError("Missing geometry column after from_features")
        if gdf.crs is None:
            gdf.set_crs("EPSG:4326", inplace=True, allow_override=True)
    except Exception:
        # Fallback: build manually
        from shapely.geometry import shape as _shape

        rows, geoms = [], []
        for feat in blocks:
            props = dict(feat.get("properties", {}))
            geom = feat.get("geometry")
            rows.append(props)
            geoms.append(_shape(geom) if geom else None)
        gdf = gpd.GeoDataFrame(rows, geometry=geoms, crs="EPSG:4326")

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_file(output_path, driver="GeoJSON")
