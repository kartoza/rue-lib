"""Geospatial functionality for rue-lib.

Note: This module requires GDAL and GeoPandas to be installed.
Install with: pip install rue-lib[geo]
Or use the Nix development environment: nix develop
"""

import math

from shapely import GEOSException

try:
    import geopandas as gpd
    from osgeo import gdal, ogr

    HAS_GEO = True
except ImportError:
    HAS_GEO = False
    gdal = None
    ogr = None
    gpd = None


def _check_geo_available():
    """Check if geospatial dependencies are available."""
    if not HAS_GEO:
        raise ImportError(
            "Geospatial dependencies not installed. Install with: pip install rue-lib[geo]"
        )


def get_gdal_version() -> str:
    """
    Get the GDAL version.

    Returns:
        GDAL version string.

    Raises:
        ImportError: If GDAL is not installed.
    """
    _check_geo_available()
    return gdal.VersionInfo("VERSION_NUM")


def create_sample_geodataframe() -> "gpd.GeoDataFrame":
    """
    Create a sample GeoDataFrame.

    Returns:
        A simple GeoDataFrame with sample data.

    Raises:
        ImportError: If GeoPandas is not installed.
    """
    _check_geo_available()
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

    Raises:
        ImportError: If GDAL/OGR is not installed.
    """
    _check_geo_available()
    return ogr.GetDriverCount()


def to_metric_crs(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Project to local UTM CRS for metric operations."""
    if gdf.crs is None:
        gdf = gdf.set_crs(4326)

    try:
        centroid = gdf.union_all().centroid
    except (GEOSException, ValueError, TypeError):
        centroid = gdf.geometry.iloc[0].centroid

    lon = centroid.x
    utm_zone = int(math.floor((lon + 180) / 6) + 1)
    is_northern = centroid.y >= 0
    epsg = 32600 + utm_zone if is_northern else 32700 + utm_zone

    return gdf.to_crs(epsg)
