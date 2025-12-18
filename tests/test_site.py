"""Tests for site module."""

import tempfile
from pathlib import Path
from unittest.mock import patch

import geopandas as gpd
import pytest
from shapely.geometry import Polygon

from rue_lib.site.config import SiteConfig
from rue_lib.site.runner import generate_parcels


class TestSiteConfig:
    """Test SiteConfig dataclass."""

    def test_site_config_required_fields(self):
        """Test SiteConfig creation with required fields."""
        config = SiteConfig(site_path="site.geojson", roads_path="roads.geojson")

        assert config.site_path == "site.geojson"
        assert config.roads_path == "roads.geojson"
        # Check default values
        assert config.output_dir == "outputs"
        assert config.geopackage_path == "outputs/output.gpkg"

    def test_site_config_custom_values(self):
        """Test SiteConfig with custom values."""
        config = SiteConfig(
            site_path="/path/to/site.geojson",
            roads_path="/path/to/roads.geojson",
            output_dir="/custom/output",
            geopackage_path="/custom/data.gpkg",
            road_arterial_width_m=25.0,
            road_secondary_width_m=18.0,
            road_local_width_m=12.0,
        )

        assert config.site_path == "/path/to/site.geojson"
        assert config.roads_path == "/path/to/roads.geojson"
        assert config.output_dir == "/custom/output"
        assert config.geopackage_path == "/custom/data.gpkg"
        assert config.road_arterial_width_m == 25.0
        assert config.road_secondary_width_m == 18.0
        assert config.road_local_width_m == 12.0

    def test_site_config_default_road_widths(self):
        """Test SiteConfig default road width values."""
        config = SiteConfig(site_path="site.geojson", roads_path="roads.geojson")

        assert config.road_arterial_width_m == 20.0
        assert config.road_secondary_width_m == 15.0
        assert config.road_local_width_m == 10.0


class TestGenerateParcels:
    """Test the generate_parcels function."""

    @pytest.fixture
    def sample_site_gdf(self):
        """Create sample site GeoDataFrame."""
        polygon = Polygon([(0, 0), (100, 0), (100, 100), (0, 100), (0, 0)])
        return gpd.GeoDataFrame({"geometry": [polygon]}, crs="EPSG:4326")

    @pytest.fixture
    def sample_roads_gdf(self):
        """Create sample roads GeoDataFrame."""
        from shapely.geometry import LineString

        road = LineString([(10, 50), (90, 50)])
        return gpd.GeoDataFrame({"geometry": [road], "type": ["arterial"]}, crs="EPSG:4326")

    @pytest.fixture
    def sample_config(self):
        """Create sample configuration."""
        return SiteConfig(
            site_path="test_site.geojson",
            roads_path="test_roads.geojson",
            output_dir="test_output",
            geopackage_path="test_output/test.gpkg",
        )

    @patch("rue_lib.site.runner.read_site")
    @patch("rue_lib.site.runner.read_roads")
    @patch("rue_lib.site.runner.prepare_geopackage")
    @patch("rue_lib.site.runner.buffer_roads")
    @patch("rue_lib.site.runner.erase_layer")
    @patch("rue_lib.site.runner.remove_layer_from_gpkg")
    @patch("rue_lib.site.runner.save_geojson")
    @patch("rue_lib.site.runner.FinancialSite")
    @patch("geopandas.read_file")
    def test_generate_parcels_with_roads(
        self,
        mock_read_file,
        mock_financial,
        mock_save_geojson,
        mock_remove_layer,
        mock_erase_layer,
        mock_buffer_roads,
        mock_prepare_geopackage,
        mock_read_roads,
        mock_read_site,
        sample_site_gdf,
        sample_roads_gdf,
        sample_config,
    ):
        """Test generate_parcels with roads buffering."""
        # Setup mocks
        mock_read_site.return_value = sample_site_gdf
        mock_read_roads.return_value = sample_roads_gdf

        # Mock buffered roads (non-empty)
        buffered_roads = sample_roads_gdf.copy()
        buffered_roads.geometry = buffered_roads.geometry.buffer(5)
        mock_buffer_roads.return_value = buffered_roads

        # Mock parcels after erase operation
        parcels_gdf = sample_site_gdf.copy()
        parcels_gdf["area_m2"] = [10000.0]
        parcels_gdf["id"] = [1]
        mock_read_file.return_value = parcels_gdf

        with tempfile.TemporaryDirectory() as temp_dir:
            config = SiteConfig(
                site_path="test_site.geojson",
                roads_path="test_roads.geojson",
                output_dir=temp_dir,
                geopackage_path=f"{temp_dir}/test.gpkg",
            )

            result = generate_parcels(config)

            # Verify calls
            mock_read_site.assert_called_once_with("test_site.geojson")
            mock_read_roads.assert_called_once_with("test_roads.geojson")
            mock_prepare_geopackage.assert_called_once()
            mock_buffer_roads.assert_called()
            mock_erase_layer.assert_called_once()
            mock_save_geojson.assert_called()
            mock_financial.assert_called_once()

            # Check return value
            assert isinstance(result, Path)
            assert str(result).endswith("roads.geojson")

    @patch("rue_lib.site.runner.read_site")
    @patch("rue_lib.site.runner.read_roads")
    @patch("rue_lib.site.runner.prepare_geopackage")
    @patch("rue_lib.site.runner.buffer_roads")
    @patch("rue_lib.site.runner.remove_layer_from_gpkg")
    @patch("rue_lib.site.runner.save_geojson")
    @patch("rue_lib.site.runner.FinancialSite")
    @patch("geopandas.read_file")
    def test_generate_parcels_empty_roads(
        self,
        mock_read_file,
        mock_financial,
        mock_save_geojson,
        mock_remove_layer,
        mock_buffer_roads,
        mock_prepare_geopackage,
        mock_read_roads,
        mock_read_site,
        sample_site_gdf,
        sample_config,
    ):
        """Test generate_parcels with empty roads."""
        # Setup mocks
        mock_read_site.return_value = sample_site_gdf

        # Empty roads
        empty_roads = gpd.GeoDataFrame({"geometry": []}, crs="EPSG:4326")
        mock_read_roads.return_value = empty_roads

        # Mock empty buffered roads
        mock_buffer_roads.return_value = empty_roads

        # Mock parcels (same as site since no roads to erase)
        parcels_gdf = sample_site_gdf.copy()
        parcels_gdf["area_m2"] = [10000.0]
        parcels_gdf["id"] = [1]
        mock_read_file.return_value = parcels_gdf

        with tempfile.TemporaryDirectory() as temp_dir:
            config = SiteConfig(
                site_path="test_site.geojson",
                roads_path="test_roads.geojson",
                output_dir=temp_dir,
                geopackage_path=f"{temp_dir}/test.gpkg",
            )

            result = generate_parcels(config)

            # Verify no erase operation when roads are empty
            mock_read_site.assert_called_once()
            mock_read_roads.assert_called_once()
            mock_buffer_roads.assert_called()

            # Should write site directly to parcels layer when roads are empty
            assert isinstance(result, Path)

    @patch("rue_lib.site.runner.read_site")
    @patch("rue_lib.site.runner.read_roads")
    @patch("rue_lib.site.runner.to_metric_crs")
    def test_generate_parcels_coordinate_transformation(
        self, mock_to_metric_crs, mock_read_roads, mock_read_site, sample_site_gdf, sample_roads_gdf
    ):
        """Test coordinate transformation logic."""
        # Mock geographic CRS data
        geo_site = sample_site_gdf.copy()
        geo_site.crs = "EPSG:4326"  # Geographic CRS
        geo_roads = sample_roads_gdf.copy()
        geo_roads.crs = "EPSG:4326"

        mock_read_site.return_value = geo_site
        mock_read_roads.return_value = geo_roads

        # Mock metric CRS versions
        metric_site = geo_site.to_crs("EPSG:3857")
        metric_roads = geo_roads.to_crs("EPSG:3857")
        mock_to_metric_crs.side_effect = [metric_site, metric_roads]

        with (
            patch("rue_lib.site.runner.prepare_geopackage"),
            patch("rue_lib.site.runner.buffer_roads") as mock_buffer,
            patch("rue_lib.site.runner.save_geojson"),
            patch("rue_lib.site.runner.FinancialSite"),
            patch("geopandas.read_file"),
        ):
            mock_buffer.return_value = gpd.GeoDataFrame({"geometry": []}, crs="EPSG:3857")

            config = SiteConfig(site_path="test_site.geojson", roads_path="test_roads.geojson")

            generate_parcels(config)

            # Should transform both datasets to metric CRS
            assert mock_to_metric_crs.call_count == 2

    def test_generate_parcels_creates_output_directory(self, sample_config):
        """Test that generate_parcels creates output directory."""
        with (
            patch("rue_lib.site.runner.read_site"),
            patch("rue_lib.site.runner.read_roads"),
            patch("rue_lib.site.runner.prepare_geopackage"),
            patch("rue_lib.site.runner.buffer_roads"),
            patch("rue_lib.site.runner.save_geojson"),
            patch("rue_lib.site.runner.FinancialSite"),
            patch("geopandas.read_file"),
            patch("pathlib.Path.mkdir") as mock_mkdir,
        ):
            # Mock empty roads to avoid erase operation
            empty_gdf = gpd.GeoDataFrame({"geometry": []}, crs="EPSG:4326")
            with (
                patch("rue_lib.site.runner.read_site", return_value=empty_gdf),
                patch("rue_lib.site.runner.read_roads", return_value=empty_gdf),
                patch("rue_lib.site.runner.buffer_roads", return_value=empty_gdf),
            ):
                generate_parcels(sample_config)

                mock_mkdir.assert_called_with(parents=True, exist_ok=True)
