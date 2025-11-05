"""Tests for geospatial functionality."""

import pytest

from rue_lib.geo import (
    HAS_GEO,
    create_sample_geodataframe,
    get_driver_count,
    get_gdal_version,
)

# Skip all tests in this module if geo dependencies are not available
pytestmark = pytest.mark.skipif(not HAS_GEO, reason="Geo dependencies not installed")


def test_get_gdal_version():
    """Test GDAL version retrieval."""
    version = get_gdal_version()
    assert version is not None
    assert len(version) > 0


def test_create_sample_geodataframe():
    """Test sample GeoDataFrame creation."""
    import geopandas as gpd

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
