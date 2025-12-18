"""Tests for geospatial functionality."""

from unittest.mock import patch

import pytest

from rue_lib.geo import (
    HAS_GEO,
    _check_geo_available,
    create_sample_geodataframe,
    get_driver_count,
    get_gdal_version,
    to_metric_crs,
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

    # Test specific data content
    cities = gdf["city"].tolist()
    assert "London" in cities
    assert "Paris" in cities
    assert "Berlin" in cities

    # Test populations are integers
    assert all(isinstance(pop, int) for pop in gdf["population"])

    # Test geometries are Points
    from shapely.geometry import Point

    assert all(isinstance(geom, Point) for geom in gdf.geometry)


def test_get_driver_count():
    """Test OGR driver count."""
    count = get_driver_count()
    assert count > 0
    assert isinstance(count, int)


class TestToMetricCRS:
    """Test the to_metric_crs function."""

    def test_to_metric_crs_london(self):
        """Test UTM conversion for London (northern hemisphere)."""
        import geopandas as gpd
        from shapely.geometry import Point

        # London coordinates
        gdf = gpd.GeoDataFrame(
            {
                "id": [1],
                "geometry": [Point(-0.1276, 51.5074)],  # London
            },
            crs="EPSG:4326",
        )

        result = to_metric_crs(gdf)

        # Should be UTM Zone 30N (EPSG:32630)
        assert result.crs.to_epsg() == 32630
        assert len(result) == 1

    def test_to_metric_crs_sydney(self):
        """Test UTM conversion for Sydney (southern hemisphere)."""
        import geopandas as gpd
        from shapely.geometry import Point

        # Sydney coordinates
        gdf = gpd.GeoDataFrame(
            {
                "id": [1],
                "geometry": [Point(151.2093, -33.8688)],  # Sydney
            },
            crs="EPSG:4326",
        )

        result = to_metric_crs(gdf)

        # Should be UTM Zone 56S (EPSG:32756)
        assert result.crs.to_epsg() == 32756
        assert len(result) == 1

    def test_to_metric_crs_no_crs(self):
        """Test UTM conversion when input has no CRS."""
        import geopandas as gpd
        from shapely.geometry import Point

        # Create GDF without CRS
        gdf = gpd.GeoDataFrame(
            {"id": [1], "geometry": [Point(-0.1276, 51.5074)]}
        )  # No CRS specified

        result = to_metric_crs(gdf)

        # Should assume EPSG:4326 and convert to UTM
        assert result.crs.to_epsg() == 32630

    def test_to_metric_crs_multiple_points(self):
        """Test UTM conversion with multiple points."""
        import geopandas as gpd
        from shapely.geometry import Point

        # Multiple points in same UTM zone
        gdf = gpd.GeoDataFrame(
            {
                "id": [1, 2, 3],
                "geometry": [
                    Point(2.3522, 48.8566),  # Paris
                    Point(2.1734, 41.4036),  # Barcelona
                    Point(4.9041, 52.3676),  # Amsterdam
                ],
            },
            crs="EPSG:4326",
        )

        result = to_metric_crs(gdf)

        # All should be in UTM Zone 31N
        assert result.crs.to_epsg() == 32631
        assert len(result) == 3

    def test_to_metric_crs_edge_case_centroid_calculation(self):
        """Test UTM conversion when union_all fails."""
        import geopandas as gpd
        from shapely import GEOSException
        from shapely.geometry import Point

        gdf = gpd.GeoDataFrame({"id": [1], "geometry": [Point(-0.1276, 51.5074)]}, crs="EPSG:4326")

        # Mock union_all to raise exception
        with patch.object(gdf, "union_all", side_effect=GEOSException("Mock error")):
            result = to_metric_crs(gdf)

            # Should fall back to first geometry centroid
            assert result.crs.to_epsg() == 32630
            assert len(result) == 1

    def test_to_metric_crs_utm_zone_calculation(self):
        """Test UTM zone calculation for various longitudes."""
        import geopandas as gpd
        from shapely.geometry import Point

        test_cases = [
            (-177, 0, 32601),  # Zone 1N
            (-3, 0, 32630),  # Zone 30N
            (0, 0, 32631),  # Zone 31N
            (3, 0, 32631),  # Zone 31N
            (177, 0, 32660),  # Zone 60N
            (0, -45, 32731),  # Zone 31S
        ]

        for lon, lat, expected_epsg in test_cases:
            gdf = gpd.GeoDataFrame({"id": [1], "geometry": [Point(lon, lat)]}, crs="EPSG:4326")

            result = to_metric_crs(gdf)
            assert result.crs.to_epsg() == expected_epsg


class TestGeoAvailabilityChecks:
    """Test geo availability checking functions."""

    def test_check_geo_available_when_available(self):
        """Test _check_geo_available when geo is available."""
        # Should not raise when HAS_GEO is True
        if HAS_GEO:
            _check_geo_available()  # Should not raise

    @pytest.mark.skipif(HAS_GEO, reason="Only test when geo is not available")
    def test_check_geo_available_when_not_available(self):
        """Test _check_geo_available when geo is not available."""
        with pytest.raises(ImportError, match="Geospatial dependencies not installed"):
            _check_geo_available()

    def test_has_geo_flag(self):
        """Test that HAS_GEO flag is boolean."""
        assert isinstance(HAS_GEO, bool)

    @patch("rue_lib.geo.HAS_GEO", False)
    def test_functions_raise_when_geo_unavailable(self):
        """Test that functions raise ImportError when geo is unavailable."""
        with pytest.raises(ImportError):
            get_gdal_version()

        with pytest.raises(ImportError):
            get_driver_count()

        with pytest.raises(ImportError):
            create_sample_geodataframe()
