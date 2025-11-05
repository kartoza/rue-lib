"""Tests for geospatial functionality."""

import geopandas as gpd

from rue_lib.geo import create_sample_geodataframe, get_driver_count, get_gdal_version


def test_get_gdal_version():
    """Test GDAL version retrieval."""
    version = get_gdal_version()
    assert version is not None
    assert len(version) > 0


def test_create_sample_geodataframe():
    """Test sample GeoDataFrame creation."""
    gdf = create_sample_geodataframe()
    assert isinstance(gdf, gpd.GeoDataFrame)
    assert len(gdf) == 3
    assert "city" in gdf.columns
    assert "population" in gdf.columns
    assert "geometry" in gdf.columns
    assert gdf.crs.to_string() == "EPSG:4326"


def test_get_driver_count():
    """Test OGR driver count."""
    count = get_driver_count()
    assert count > 0
