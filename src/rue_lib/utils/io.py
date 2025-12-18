"""Input/output utilities."""

from __future__ import annotations

import json
import shutil
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
    geometry = mapping(geom)
    # Ensure coordinates are lists, not tuples (GeoJSON spec requirement)
    coordinates = geometry.get("coordinates")
    if isinstance(coordinates, tuple):
        geometry["coordinates"] = list(coordinates)
    return {"type": "Feature", "properties": props or {}, "geometry": geometry}


def prepare_geopackage(output_path: str | Path, template_path: str | Path | None = None) -> Path:
    """
    Prepare a geopackage for writing by copying a template if it doesn't exist.

    This function checks if a geopackage file already exists at the output path.
    If it doesn't exist, it copies a template geopackage to the output location.
    The template contains pre-configured layers with QGIS styling for beautiful
    cartographic output.

    Args:
        output_path: Path where the geopackage should be created
        template_path: Optional path to custom template. If None, uses built-in template

    Returns:
        Path object pointing to the prepared geopackage

    Raises:
        FileNotFoundError: If template_path is specified but doesn't exist

    Example:
        >>> output_gpkg = prepare_geopackage("outputs/results.gpkg")
        >>> print(f"Geopackage ready at: {output_gpkg}")
    """
    output_path = Path(output_path)

    # If geopackage already exists, return its path
    if output_path.exists():
        return output_path

    # Ensure parent directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Determine template path
    if template_path is None:
        # Use built-in template
        current_file = Path(__file__)
        template_path = current_file.parent.parent / "templates" / "template.gpkg"
    else:
        template_path = Path(template_path)

    # Check template exists
    if not template_path.exists():
        raise FileNotFoundError(f"Template geopackage not found: {template_path}")

    # Copy template to output location
    shutil.copy2(template_path, output_path)

    return output_path
