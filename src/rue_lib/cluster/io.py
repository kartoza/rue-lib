# src/rue_lib/cluster/io.py
"""Input/output operations for cluster generation."""

from pathlib import Path

import geopandas as gpd


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


def save_geopackage(gdf: gpd.GeoDataFrame, path: Path, layer: str) -> None:
    """
    Save GeoDataFrame to GeoPackage file.

    Args:
        gdf: GeoDataFrame to save
        path: Output path for GeoPackage file
        layer: Name of layer to create/update in GeoPackage
    """
    gdf.to_file(path, layer=layer, driver="GPKG")
