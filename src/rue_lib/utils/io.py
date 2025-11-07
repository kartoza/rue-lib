"""Input/output utilities."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from shapely.geometry import mapping
from shapely.geometry.base import BaseGeometry


def load_geojson(path: Path) -> dict:
    """Load a GeoJSON file."""
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def save_geojson(path: Path, features: list[dict]) -> None:
    """Save features to a GeoJSON file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump({"type": "FeatureCollection", "features": features}, f)


def save_json(path: Path, data: dict[str, Any]) -> None:
    """Save data to a JSON file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def to_feature(geom: BaseGeometry, props: dict[str, Any] | None = None) -> dict:
    """Convert a Shapely geometry to a GeoJSON feature."""
    return {"type": "Feature", "properties": props or {}, "geometry": mapping(geom)}
