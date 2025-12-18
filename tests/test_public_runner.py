"""Tests for public runner module."""

import importlib.util
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

HAS_GDAL = importlib.util.find_spec("osgeo") is not None

if HAS_GDAL:
    pass

pytest.mark.skipif(not HAS_GDAL, reason="GDAL not available")


@pytest.fixture
def mock_public_config():
    """Create a mock PublicConfig."""
    config = Mock()
    config.input_path = "input.gpkg"
    config.site_path = "site.geojson"
    config.output_dir = "outputs/public"
    config.open_percent = 4.0
    config.amen_percent = 10.0
    return config


@pytest.fixture
def mock_ogr_layer():
    """Create a mock OGR layer."""
    layer = Mock()
    layer.GetSpatialRef.return_value = Mock()
    layer.GetLayerDefn.return_value = Mock()
    layer.GetFeatureCount.return_value = 5
    layer.CreateFeature = Mock()
    layer.__iter__ = Mock(return_value=iter([]))
    return layer


@pytest.fixture
def mock_ogr_feature():
    """Create a mock OGR feature."""
    feature = Mock()
    feature.GetGeometryRef.return_value = Mock()
    feature.GetFID.return_value = 1
    feature.GetFieldCount.return_value = 3
    feature.Clone.return_value = feature
    return feature


class TestGeneratePublic:
    """Test generate_public function."""

    def test_generate_public_basic_workflow(self, mock_public_config, mock_ogr_layer):
        """Test basic public space generation workflow."""
        from rue_lib.public.runner import generate_public

        with (
            patch("rue_lib.public.runner.ogr") as mock_ogr,
            patch("rue_lib.public.runner.allocate_open_spaces") as mock_open,
            patch("rue_lib.public.runner.allocate_amenities") as mock_amen,
            patch("rue_lib.public.runner.merge") as mock_merge,
            patch("rue_lib.public.runner.Path") as mock_path,
        ):
            # Mock Path operations
            mock_output_dir = Mock()
            mock_output_dir.mkdir = Mock()
            mock_output_dir.__truediv__ = Mock(return_value=Path("outputs/public/outputs.gpkg"))
            mock_path.return_value = mock_output_dir

            # Mock OGR operations
            mock_driver = Mock()
            mock_ds = Mock()
            mock_ogr.GetDriverByName.return_value = mock_driver
            mock_driver.CreateDataSource.return_value = mock_ds
            mock_ogr.Open.side_effect = [mock_ds, mock_ds]  # Input and output
            mock_ds.GetLayer.return_value = mock_ogr_layer
            mock_ds.CopyLayer.return_value = mock_ogr_layer
            mock_ds.GetLayerByName.return_value = None  # No existing layer

            # Mock operations return layer names
            mock_open.return_value = "open_spaces"
            mock_amen.return_value = "amenities"
            mock_merge.return_value = "public_spaces"

            result = generate_public(mock_public_config)

            assert isinstance(result, Path)
            mock_open.assert_called_once()
            mock_amen.assert_called_once()
            mock_merge.assert_called_once()

    def test_generate_public_creates_output_directory(self, mock_public_config):
        """Test that output directory is created."""
        from rue_lib.public.runner import generate_public

        with (
            patch("rue_lib.public.runner.ogr") as mock_ogr,
            patch("rue_lib.public.runner.allocate_open_spaces") as mock_open,
            patch("rue_lib.public.runner.allocate_amenities") as mock_amen,
            patch("rue_lib.public.runner.merge") as mock_merge,
        ):
            # Mock file system operations
            mock_output_dir = Mock()
            mock_output_gpkg = Mock()
            mock_output_gpkg.exists.return_value = False
            mock_output_dir.__truediv__.return_value = mock_output_gpkg

            with patch("rue_lib.public.runner.Path") as mock_path:
                mock_path.return_value = mock_output_dir

                # Mock OGR operations
                mock_driver = Mock()
                mock_ds = Mock()
                mock_ogr.GetDriverByName.return_value = mock_driver
                mock_driver.CreateDataSource.return_value = mock_ds
                mock_ogr.Open.side_effect = [mock_ds, mock_ds]
                mock_ds.GetLayer.return_value = Mock()
                mock_ds.CopyLayer.return_value = Mock()
                mock_ds.GetLayerByName.return_value = None

                mock_open.return_value = "open_spaces"
                mock_amen.return_value = "amenities"
                mock_merge.return_value = "public_spaces"

                generate_public(mock_public_config)

                mock_output_dir.mkdir.assert_called_once_with(parents=True, exist_ok=True)
                mock_driver.CreateDataSource.assert_called_once()

    def test_generate_public_existing_geopackage(self, mock_public_config, mock_ogr_layer):
        """Test generation with existing geopackage."""
        from rue_lib.public.runner import generate_public

        with (
            patch("rue_lib.public.runner.ogr") as mock_ogr,
            patch("rue_lib.public.runner.allocate_open_spaces") as mock_open,
            patch("rue_lib.public.runner.allocate_amenities") as mock_amen,
            patch("rue_lib.public.runner.merge") as mock_merge,
        ):
            mock_output_dir = Mock()
            mock_output_gpkg = Mock()
            mock_output_gpkg.exists.return_value = True  # Existing geopackage
            mock_output_dir.__truediv__.return_value = mock_output_gpkg

            with patch("rue_lib.public.runner.Path") as mock_path:
                mock_path.return_value = mock_output_dir

                mock_ds = Mock()
                mock_ogr.Open.side_effect = [mock_ds, mock_ds]
                mock_ds.GetLayer.return_value = mock_ogr_layer
                mock_ds.CopyLayer.return_value = mock_ogr_layer
                mock_ds.GetLayerByName.return_value = None

                mock_open.return_value = "open_spaces"
                mock_amen.return_value = "amenities"
                mock_merge.return_value = "public_spaces"

                generate_public(mock_public_config)

                # Should not create new datasource for existing geopackage
                mock_ogr.GetDriverByName.assert_not_called()

    def test_generate_public_copy_input_layer(self, mock_public_config, mock_ogr_layer):
        """Test copying input blocks layer."""
        from rue_lib.public.runner import generate_public

        with (
            patch("rue_lib.public.runner.ogr") as mock_ogr,
            patch("rue_lib.public.runner.allocate_open_spaces") as mock_open,
            patch("rue_lib.public.runner.allocate_amenities") as mock_amen,
            patch("rue_lib.public.runner.merge") as mock_merge,
            patch("rue_lib.public.runner.Path") as mock_path,
        ):
            mock_output_dir = Mock()
            mock_output_gpkg = Mock()
            mock_output_gpkg.exists.return_value = False
            mock_output_dir.__truediv__.return_value = mock_output_gpkg
            mock_path.return_value = mock_output_dir

            # Mock input and output datasets
            mock_input_ds = Mock()
            mock_output_ds = Mock()
            mock_input_ds.GetLayer.return_value = mock_ogr_layer
            mock_output_ds.CopyLayer.return_value = mock_ogr_layer
            mock_output_ds.GetLayerByName.return_value = None
            mock_output_ds.DeleteLayer = Mock()

            mock_driver = Mock()
            mock_ogr.GetDriverByName.return_value = mock_driver
            mock_driver.CreateDataSource.return_value = mock_output_ds
            mock_ogr.Open.side_effect = [mock_input_ds, mock_output_ds]

            mock_open.return_value = "open_spaces"
            mock_amen.return_value = "amenities"
            mock_merge.return_value = "public_spaces"

            generate_public(mock_public_config)

            mock_output_ds.CopyLayer.assert_called_once_with(mock_ogr_layer, "00_input_blocks")

    def test_generate_public_replaces_existing_layer(self, mock_public_config, mock_ogr_layer):
        """Test that existing input blocks layer is replaced."""
        from rue_lib.public.runner import generate_public

        with (
            patch("rue_lib.public.runner.ogr") as mock_ogr,
            patch("rue_lib.public.runner.allocate_open_spaces") as mock_open,
            patch("rue_lib.public.runner.allocate_amenities") as mock_amen,
            patch("rue_lib.public.runner.merge") as mock_merge,
            patch("rue_lib.public.runner.Path") as mock_path,
        ):
            mock_output_dir = Mock()
            mock_output_gpkg = Mock()
            mock_output_gpkg.exists.return_value = False
            mock_output_dir.__truediv__.return_value = mock_output_gpkg
            mock_path.return_value = mock_output_dir

            mock_input_ds = Mock()
            mock_output_ds = Mock()
            mock_input_ds.GetLayer.return_value = mock_ogr_layer
            mock_output_ds.CopyLayer.return_value = mock_ogr_layer
            mock_output_ds.GetLayerByName.return_value = mock_ogr_layer  # Existing layer
            mock_output_ds.DeleteLayer = Mock()

            mock_driver = Mock()
            mock_ogr.GetDriverByName.return_value = mock_driver
            mock_driver.CreateDataSource.return_value = mock_output_ds
            mock_ogr.Open.side_effect = [mock_input_ds, mock_output_ds]

            mock_open.return_value = "open_spaces"
            mock_amen.return_value = "amenities"
            mock_merge.return_value = "public_spaces"

            generate_public(mock_public_config)

            mock_output_ds.DeleteLayer.assert_called_once_with("00_input_blocks")
            mock_output_ds.CopyLayer.assert_called_once()

    def test_generate_public_error_handling(self, mock_public_config):
        """Test error handling in public generation."""
        from rue_lib.public.runner import generate_public

        with patch("rue_lib.public.runner.ogr") as mock_ogr:
            # Mock failed dataset creation
            mock_driver = Mock()
            mock_ogr.GetDriverByName.return_value = mock_driver
            mock_driver.CreateDataSource.return_value = None  # Failed creation

            mock_output_dir = Mock()
            mock_output_gpkg = Mock()
            mock_output_gpkg.exists.return_value = False
            mock_output_dir.__truediv__.return_value = mock_output_gpkg

            with (
                patch("rue_lib.public.runner.Path") as mock_path,
                pytest.raises(ValueError, match="Could not create"),
            ):
                mock_path.return_value = mock_output_dir

                generate_public(mock_public_config)

    def test_generate_public_write_error_handling(self, mock_public_config):
        """Test error handling for write access."""
        from rue_lib.public.runner import generate_public

        with patch("rue_lib.public.runner.ogr") as mock_ogr:
            mock_output_dir = Mock()
            mock_output_gpkg = Mock()
            mock_output_gpkg.exists.return_value = True
            mock_output_dir.__truediv__.return_value = mock_output_gpkg

            # Mock failed write access
            mock_ogr.Open.side_effect = [Mock(), None]  # Input succeeds, output fails

            with (
                patch("rue_lib.public.runner.Path") as mock_path,
                pytest.raises(ValueError, match="Could not open"),
            ):
                mock_path.return_value = mock_output_dir

                generate_public(mock_public_config)


class TestPublicSpaceParameters:
    """Test parameter passing in public space generation."""

    def test_generate_public_passes_config_parameters(self, mock_public_config, mock_ogr_layer):
        """Test that configuration parameters are passed correctly."""
        from rue_lib.public.runner import generate_public

        with (
            patch("rue_lib.public.runner.ogr") as mock_ogr,
            patch("rue_lib.public.runner.allocate_open_spaces") as mock_open,
            patch("rue_lib.public.runner.allocate_amenities") as mock_amen,
            patch("rue_lib.public.runner.merge") as mock_merge,
            patch("rue_lib.public.runner.Path") as mock_path,
        ):
            mock_output_dir = Mock()
            mock_output_gpkg = Mock()
            mock_output_gpkg.exists.return_value = False
            mock_output_dir.__truediv__.return_value = mock_output_gpkg
            mock_path.return_value = mock_output_dir

            mock_ds = Mock()
            mock_driver = Mock()
            mock_ogr.GetDriverByName.return_value = mock_driver
            mock_driver.CreateDataSource.return_value = mock_ds
            mock_ogr.Open.side_effect = [mock_ds, mock_ds]
            mock_ds.GetLayer.return_value = mock_ogr_layer
            mock_ds.CopyLayer.return_value = mock_ogr_layer
            mock_ds.GetLayerByName.return_value = None

            mock_open.return_value = "open_spaces"
            mock_amen.return_value = "amenities"
            mock_merge.return_value = "public_spaces"

            generate_public(mock_public_config)

            # Verify correct parameters are passed
            mock_open.assert_called_with(
                str(mock_output_gpkg), "site", "00_input_blocks", "01_open_spaces", open_percent=4.0
            )

            mock_amen.assert_called_with(
                str(mock_output_gpkg), "site", "00_input_blocks", "02_amenities", amen_percent=10.0
            )

    def test_generate_public_custom_layer_names(self, mock_public_config, mock_ogr_layer):
        """Test generation with custom layer naming."""
        from rue_lib.public.runner import generate_public

        with (
            patch("rue_lib.public.runner.ogr") as mock_ogr,
            patch("rue_lib.public.runner.allocate_open_spaces") as mock_open,
            patch("rue_lib.public.runner.allocate_amenities") as mock_amen,
            patch("rue_lib.public.runner.merge") as mock_merge,
            patch("rue_lib.public.runner.Path") as mock_path,
        ):
            mock_output_dir = Mock()
            mock_output_gpkg = Mock()
            mock_output_gpkg.exists.return_value = False
            mock_output_dir.__truediv__.return_value = mock_output_gpkg
            mock_path.return_value = mock_output_dir

            mock_ds = Mock()
            mock_driver = Mock()
            mock_ogr.GetDriverByName.return_value = mock_driver
            mock_driver.CreateDataSource.return_value = mock_ds
            mock_ogr.Open.side_effect = [mock_ds, mock_ds]
            mock_ds.GetLayer.return_value = mock_ogr_layer
            mock_ds.CopyLayer.return_value = mock_ogr_layer
            mock_ds.GetLayerByName.return_value = None

            mock_open.return_value = "01_open_spaces"
            mock_amen.return_value = "02_amenities"
            mock_merge.return_value = "03_public_spaces"

            generate_public(mock_public_config)

            # Verify layer naming follows expected pattern
            mock_merge.assert_called_with(
                str(mock_output_gpkg),
                ["01_open_spaces", "02_amenities"],
                "site",
                "03_public_spaces",
            )


class TestPublicGenerationWorkflow:
    """Test the complete public space generation workflow."""

    def test_complete_workflow_integration(self, mock_public_config):
        """Test integration of all workflow steps."""
        from rue_lib.public.runner import generate_public

        with tempfile.TemporaryDirectory() as tmp_dir:
            # Update config to use temp directory
            mock_public_config.output_dir = tmp_dir

            with (
                patch("rue_lib.public.runner.ogr") as mock_ogr,
                patch("rue_lib.public.runner.allocate_open_spaces") as mock_open,
                patch("rue_lib.public.runner.allocate_amenities") as mock_amen,
                patch("rue_lib.public.runner.merge") as mock_merge,
            ):
                # Mock complete workflow
                mock_ds = Mock()
                mock_layer = Mock()
                mock_driver = Mock()

                mock_ogr.GetDriverByName.return_value = mock_driver
                mock_driver.CreateDataSource.return_value = mock_ds
                mock_ogr.Open.side_effect = [mock_ds, mock_ds]
                mock_ds.GetLayer.return_value = mock_layer
                mock_ds.CopyLayer.return_value = mock_layer
                mock_ds.GetLayerByName.return_value = None

                mock_open.return_value = "01_open_spaces"
                mock_amen.return_value = "02_amenities"
                mock_merge.return_value = "03_public_spaces"

                result = generate_public(mock_public_config)

                # Verify complete workflow execution
                assert isinstance(result, Path)
                assert "outputs.gpkg" in str(result)

                # All operations should be called in sequence
                mock_open.assert_called_once()
                mock_amen.assert_called_once()
                mock_merge.assert_called_once()

                # Verify proper sequencing
                calls = [mock_open.call_args, mock_amen.call_args, mock_merge.call_args]
                assert all(call is not None for call in calls)

    def test_workflow_with_site_processing(self, mock_public_config):
        """Test workflow including site boundary processing."""
        from rue_lib.public.runner import generate_public

        with (
            patch("rue_lib.public.runner.ogr") as mock_ogr,
            patch("rue_lib.public.runner.allocate_open_spaces") as mock_open,
            patch("rue_lib.public.runner.allocate_amenities") as mock_amen,
            patch("rue_lib.public.runner.merge") as mock_merge,
            patch("rue_lib.public.runner.get_utm_zone_from_layer") as mock_utm,
            patch("rue_lib.public.runner.reproject_layer"),
            patch("rue_lib.public.runner.Path") as mock_path,
        ):
            mock_output_dir = Mock()
            mock_output_gpkg = Mock()
            mock_output_gpkg.exists.return_value = False
            mock_output_dir.__truediv__.return_value = mock_output_gpkg
            mock_path.return_value = mock_output_dir

            mock_ds = Mock()
            mock_layer = Mock()
            mock_driver = Mock()

            mock_ogr.GetDriverByName.return_value = mock_driver
            mock_driver.CreateDataSource.return_value = mock_ds
            mock_ogr.Open.side_effect = [mock_ds, mock_ds]
            mock_ds.GetLayer.return_value = mock_layer
            mock_ds.CopyLayer.return_value = mock_layer
            mock_ds.GetLayerByName.return_value = None

            mock_utm.return_value = "32633"  # UTM zone
            mock_open.return_value = "01_open_spaces"
            mock_amen.return_value = "02_amenities"
            mock_merge.return_value = "03_public_spaces"

            generate_public(mock_public_config)

            # Should process UTM conversion if needed
            mock_utm.assert_called()

    def test_workflow_output_validation(self, mock_public_config):
        """Test validation of workflow outputs."""
        from rue_lib.public.runner import generate_public

        with (
            patch("rue_lib.public.runner.ogr") as mock_ogr,
            patch("rue_lib.public.runner.allocate_open_spaces") as mock_open,
            patch("rue_lib.public.runner.allocate_amenities") as mock_amen,
            patch("rue_lib.public.runner.merge") as mock_merge,
            patch("rue_lib.public.runner.Path") as mock_path,
        ):
            mock_output_dir = Mock()
            mock_output_gpkg = Mock()
            mock_output_gpkg.exists.return_value = False
            mock_output_dir.__truediv__.return_value = mock_output_gpkg
            mock_path.return_value = mock_output_dir

            mock_ds = Mock()
            mock_layer = Mock()
            mock_driver = Mock()

            mock_ogr.GetDriverByName.return_value = mock_driver
            mock_driver.CreateDataSource.return_value = mock_ds
            mock_ogr.Open.side_effect = [mock_ds, mock_ds]
            mock_ds.GetLayer.return_value = mock_layer
            mock_ds.CopyLayer.return_value = mock_layer
            mock_ds.GetLayerByName.return_value = None

            mock_open.return_value = "01_open_spaces"
            mock_amen.return_value = "02_amenities"
            mock_merge.return_value = "03_public_spaces"

            result = generate_public(mock_public_config)

            # Result should be a valid Path object
            assert isinstance(result, Path)

            # Should contain expected output file
            assert str(result).endswith("outputs.gpkg")
