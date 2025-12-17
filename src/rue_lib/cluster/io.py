# src/rue_lib/cluster/io.py
"""Input/output operations for cluster generation."""

from pathlib import Path

import geopandas as gpd

from ..utils.io import prepare_geopackage


def read_site(path: str) -> gpd.GeoDataFrame:
    """
    Read site boundary from file.

    Args:
        path: Path to GeoJSON or GeoPackage file containing site boundary

    Returns:
        GeoDataFrame containing site boundary geometry
    """
    return gpd.read_file(path)


def read_roads(path: str) -> gpd.GeoDataFrame:
    """
    Read roads from file.

    Args:
        path: Path to GeoJSON or GeoPackage file containing road network

    Returns:
        GeoDataFrame containing road geometries with road_type attributes
    """
    return gpd.read_file(path)


def read_blocks(path: str, layer: str = None) -> gpd.GeoDataFrame:
    """
    Read blocks from file.

    Args:
        path: Path to file (GeoJSON or GeoPackage)
        layer: Layer name (for GeoPackage files)

    Returns:
        GeoDataFrame containing blocks
    """
    if layer:
        return gpd.read_file(path, layer=layer)
    return gpd.read_file(path)


def save_geojson(gdf: gpd.GeoDataFrame, path: Path) -> None:
    """
    Save GeoDataFrame to GeoJSON file.

    Args:
        gdf: GeoDataFrame to save
        path: Output path for GeoJSON file
    """
    gdf.to_file(path, driver="GeoJSON")


def explode_polygon_geodataframe(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Explode MultiPolygon geometries into individual Polygon features.

    Args:
        gdf: Input GeoDataFrame that may contain MultiPolygon geometries

    Returns:
        GeoDataFrame with all geometries as single Polygons
    """
    # Only process if there are any MultiPolygon geometries
    if gdf.empty or not any(gdf.geom_type == "MultiPolygon"):
        return gdf

    # Use explode() to convert MultiPolygons to individual Polygons
    exploded = gdf.explode(index_parts=False)

    # Reset index to ensure unique indices
    exploded = exploded.reset_index(drop=True)

    return exploded


def save_geopackage(gdf: gpd.GeoDataFrame, path: Path, layer: str) -> None:
    """
    Save GeoDataFrame to GeoPackage file.

    This function prepares the geopackage using a template if it doesn't exist,
    ensuring consistent styling and structure for QGIS visualization.
    All MultiPolygon geometries are automatically exploded into individual Polygons.

    Args:
        gdf: GeoDataFrame to save
        path: Output path for GeoPackage file
        layer: Name of layer to create/update in GeoPackage
    """
    # Prepare the geopackage (copy from template if needed)
    prepare_geopackage(path)

    # Explode MultiPolygons into individual Polygons if needed
    gdf_exploded = explode_polygon_geodataframe(gdf)

    # Save the GeoDataFrame to the prepared geopackage
    gdf_exploded.to_file(path, layer=layer, driver="GPKG")
