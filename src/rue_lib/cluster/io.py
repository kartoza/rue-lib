# src/rue_lib/cluster/io.py
"""Input/output operations for cluster generation."""

from pathlib import Path

import geopandas as gpd

from ..utils.io import prepare_geopackage


def safe_geodataframe(data=None, geometry=None, crs=None, **kwargs) -> gpd.GeoDataFrame:
    """
    Safely create a GeoDataFrame, only assigning CRS if geometry is present.

    Args:
        data: DataFrame data
        geometry: Geometry column
        crs: Coordinate reference system
        **kwargs: Additional arguments for GeoDataFrame

    Returns:
        GeoDataFrame with CRS only if geometry is present
    """
    if data is None:
        data = []

    # If data is a list, check if it's empty or contains geometry
    if isinstance(data, list):
        if len(data) == 0 or not any("geometry" in item for item in data if isinstance(item, dict)):
            # No geometry data, create empty GeoDataFrame without CRS
            return gpd.GeoDataFrame(data, **kwargs)

    # Check if geometry parameter is empty
    if geometry is not None and hasattr(geometry, "__len__") and len(geometry) == 0:
        # Empty geometry, create GeoDataFrame without CRS
        return gpd.GeoDataFrame(data, geometry=[], **kwargs)

    # Has geometry data, can safely assign CRS
    return gpd.GeoDataFrame(data, geometry=geometry, crs=crs, **kwargs)


def explode_polygon_geodataframe(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Explode MultiPolygon geometries into individual Polygons.

    Args:
        gdf: Input GeoDataFrame that may contain MultiPolygon geometries

    Returns:
        GeoDataFrame with all MultiPolygons exploded to individual Polygons
    """
    if gdf.empty:
        return gdf

    # Check if there's a geometry column
    if "geometry" not in gdf.columns or gdf.geometry.empty:
        return gdf

    # Explode MultiPolygons into individual Polygons
    return gdf.explode(index_parts=False, ignore_index=True)


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

    # Check if we have a valid GeoDataFrame with geometry before saving
    if gdf_exploded.empty or "geometry" not in gdf_exploded.columns:
        print(f"Warning: Cannot save empty or non-geometric data to layer '{layer}'. Skipping.")
        return

    # Ensure it's actually a GeoDataFrame
    if not isinstance(gdf_exploded, gpd.GeoDataFrame):
        print(f"Warning: Data for layer '{layer}' is not a GeoDataFrame. Converting...")
        gdf_exploded = gpd.GeoDataFrame(gdf_exploded)

    # Save the GeoDataFrame to the prepared geopackage
    gdf_exploded.to_file(path, layer=layer, driver="GPKG")
