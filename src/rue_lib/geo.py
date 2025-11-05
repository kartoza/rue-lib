"""Geospatial functionality for rue-lib."""

import geopandas as gpd
from osgeo import gdal, ogr


def get_gdal_version() -> str:
    """
    Get the GDAL version.

    Returns:
        GDAL version string.
    """
    return gdal.VersionInfo("VERSION_NUM")


def create_sample_geodataframe() -> gpd.GeoDataFrame:
    """
    Create a sample GeoDataFrame.

    Returns:
        A simple GeoDataFrame with sample data.
    """
    from shapely.geometry import Point

    # Create sample data
    data = {
        "city": ["London", "Paris", "Berlin"],
        "population": [8982000, 2161000, 3645000],
        "geometry": [
            Point(-0.1276, 51.5074),  # London
            Point(2.3522, 48.8566),  # Paris
            Point(13.4050, 52.5200),  # Berlin
        ],
    }

    return gpd.GeoDataFrame(data, crs="EPSG:4326")


def get_driver_count() -> int:
    """
    Get the number of OGR drivers available.

    Returns:
        Number of available OGR drivers.
    """
    return ogr.GetDriverCount()
