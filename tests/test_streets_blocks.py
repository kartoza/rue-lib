"""Tests for streets blocks module."""

import tempfile
from unittest.mock import Mock, patch

import pytest

try:
    from osgeo import ogr

    HAS_GDAL = True
except ImportError:
    HAS_GDAL = False

pytest.mark.skipif(not HAS_GDAL, reason="GDAL not available")


@pytest.fixture
def mock_street_config():
    """Create a mock StreetConfig."""
    config = Mock()
    config.geopackage_path = "test.gpkg"
    config.road_arterial_width_m = 20.0
    config.road_secondary_width_m = 15.0
    config.road_locals_width_m = 10.0
    config.sidewalk_width_m = 3.0
    config.on_grid_partition_depth_arterial_roads = 60.0
    config.on_grid_partition_depth_secondary_roads = 60.0
    config.perpendicular_line_length = 1000.0
    config.optimize_grid_rotation = True
    config.grid_rotation_angle_step = 5.0
    config.use_ternary_search = False
    config.clip_to_boundary = True
    config.tolerance_area_ratio = 0.70
    config.tolerance_boundary_distance = 10.0
    return config


@pytest.fixture
def mock_ogr_geometry():
    """Create a mock OGR geometry."""
    geom = Mock()
    geom.Clone.return_value = geom
    geom.GetGeometryType.return_value = ogr.wkbPolygon25D
    geom.GetGeometryCount.return_value = 2
    geom.GetGeometryRef.return_value = geom
    geom.ExportToWkb.return_value = b"\x01\x03\x00\x00\x00"
    geom.ExportToWkt.return_value = "POLYGON((0 0,1 0,1 1,0 1,0 0))"
    geom.GetEnvelope.return_value = (0, 1, 0, 1)  # (minX, maxX, minY, maxY)
    geom.Area.return_value = 1.0
    geom.Intersection.return_value = geom
    geom.Difference.return_value = geom
    geom.Union.return_value = geom
    return geom


class TestExplodeMultipolygon:
    """Test explode_multipolygon function in blocks module."""

    def test_explode_single_polygon(self, mock_ogr_geometry):
        """Test exploding a single polygon."""
        from rue_lib.streets.blocks import explode_multipolygon

        mock_ogr_geometry.GetGeometryType.return_value = ogr.wkbPolygon25D

        result = explode_multipolygon(mock_ogr_geometry)

        assert len(result) == 1
        assert result[0] == mock_ogr_geometry

    def test_explode_multipolygon(self, mock_ogr_geometry):
        """Test exploding a multipolygon."""
        from rue_lib.streets.blocks import explode_multipolygon

        mock_ogr_geometry.GetGeometryType.return_value = ogr.wkbMultiPolygon25D
        mock_ogr_geometry.GetGeometryCount.return_value = 2

        sub_geom = Mock()
        sub_geom.GetGeometryType.return_value = ogr.wkbPolygon25D
        sub_geom.Clone.return_value = sub_geom
        mock_ogr_geometry.GetGeometryRef.return_value = sub_geom

        result = explode_multipolygon(mock_ogr_geometry)

        assert len(result) == 2

    def test_explode_2d_geometry_support(self, mock_ogr_geometry):
        """Test that 2D geometries are supported."""
        from rue_lib.streets.blocks import explode_multipolygon

        mock_ogr_geometry.GetGeometryType.return_value = ogr.wkbPolygon

        result = explode_multipolygon(mock_ogr_geometry)

        assert len(result) == 1


class TestMergeLines:
    """Test merge_lines function."""

    def test_merge_lines_basic(self, mock_street_config):
        """Test basic line merging functionality."""
        from rue_lib.streets.blocks import merge_lines

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.blocks.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()
                mock_feature = Mock()
                mock_geom = Mock()

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer
                mock_feature.GetGeometryRef.return_value = mock_geom
                mock_feature.GetFieldCount.return_value = 0
                mock_layer.__iter__ = Mock(return_value=iter([mock_feature]))
                mock_layer.GetLayerDefn.return_value = Mock()

                with patch("rue_lib.streets.blocks.create_layer") as mock_create:
                    mock_output_layer = Mock()
                    mock_create.return_value = mock_output_layer

                    merge_lines(
                        gpkg_path,
                        "perpendicular_lines",
                        "arterial_edges",
                        "street_blocks",
                        "secondary_setback",
                        "arterial_setback",
                        "merged_lines",
                    )

                    mock_create.assert_called()
                    mock_output_layer.CreateFeature.assert_called()

    def test_merge_lines_with_line_boundary_extraction(self):
        """Test line merging with polygon boundary extraction."""
        from rue_lib.streets.blocks import merge_lines

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.blocks.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()
                mock_feature = Mock()
                mock_polygon = Mock()
                mock_boundary = Mock()

                # Mock polygon with boundary
                mock_polygon.GetGeometryType.return_value = ogr.wkbPolygon25D
                mock_polygon.GetBoundary.return_value = mock_boundary
                mock_feature.GetGeometryRef.return_value = mock_polygon
                mock_feature.GetFieldCount.return_value = 0

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer
                mock_layer.__iter__ = Mock(return_value=iter([mock_feature]))
                mock_layer.GetLayerDefn.return_value = Mock()

                with (
                    patch("rue_lib.streets.blocks.create_layer") as mock_create,
                    patch("rue_lib.streets.blocks.break_linestring_by_angle") as mock_break,
                ):
                    mock_output_layer = Mock()
                    mock_create.return_value = mock_output_layer
                    mock_break.return_value = [mock_boundary]

                    merge_lines(
                        gpkg_path,
                        "perpendicular_lines",
                        "arterial_edges",
                        "street_blocks",
                        "secondary_setback",
                        "arterial_setback",
                        "merged_lines",
                    )

                    mock_polygon.GetBoundary.assert_called()
                    mock_break.assert_called()

    def test_merge_lines_handles_empty_layers(self):
        """Test line merging handles empty layers gracefully."""
        from rue_lib.streets.blocks import merge_lines

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.blocks.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer
                mock_layer.__iter__ = Mock(return_value=iter([]))  # Empty layer
                mock_layer.GetLayerDefn.return_value = Mock()

                with patch("rue_lib.streets.blocks.create_layer") as mock_create:
                    mock_output_layer = Mock()
                    mock_create.return_value = mock_output_layer

                    merge_lines(
                        gpkg_path,
                        "perpendicular_lines",
                        "arterial_edges",
                        "street_blocks",
                        "secondary_setback",
                        "arterial_setback",
                        "merged_lines",
                    )

                    mock_create.assert_called()
                    # Should not crash with empty layers


class TestCreateStreetBlocks:
    """Test create_street_blocks function."""

    def test_create_street_blocks_basic(self, mock_street_config):
        """Test basic street blocks creation."""
        from rue_lib.streets.blocks import create_street_blocks

        with (
            patch("rue_lib.streets.blocks.ogr") as mock_ogr,
            patch("rue_lib.streets.blocks.merge_lines") as mock_merge,
            patch("rue_lib.streets.blocks.create_grid_from_on_grid") as mock_grid,
        ):
            mock_ds = Mock()
            mock_layer = Mock()

            mock_ogr.Open.return_value = mock_ds
            mock_ds.GetLayerByName.return_value = mock_layer

            create_street_blocks(
                mock_street_config.geopackage_path,
                "perpendicular_lines",
                "arterial_edges",
                "on_grid_blocks",
                "secondary_setback",
                "arterial_setback",
                "street_blocks",
            )

            mock_merge.assert_called_once()
            mock_grid.assert_called_once()

    def test_create_street_blocks_with_cleanup(self, mock_street_config):
        """Test street blocks creation with cleanup step."""
        from rue_lib.streets.blocks import create_street_blocks

        with (
            patch("rue_lib.streets.blocks.ogr") as mock_ogr,
            patch("rue_lib.streets.blocks.merge_lines"),
            patch("rue_lib.streets.blocks.create_grid_from_on_grid"),
            patch("rue_lib.streets.blocks.cleanup_grid_blocks") as mock_cleanup,
        ):
            mock_ds = Mock()
            mock_layer = Mock()

            mock_ogr.Open.return_value = mock_ds
            mock_ds.GetLayerByName.return_value = mock_layer

            create_street_blocks(
                mock_street_config.geopackage_path,
                "perpendicular_lines",
                "arterial_edges",
                "on_grid_blocks",
                "secondary_setback",
                "arterial_setback",
                "street_blocks",
            )

            mock_cleanup.assert_called_once()


class TestGenerateStreetBlocks:
    """Test generate_street_blocks main function."""

    def test_generate_street_blocks_full_workflow(self, mock_street_config):
        """Test full street blocks generation workflow."""
        from rue_lib.streets.blocks import generate_street_blocks

        with (
            patch("rue_lib.streets.blocks.extract_arterial_edge_lines") as mock_extract,
            patch("rue_lib.streets.blocks.create_street_blocks") as mock_create,
        ):
            generate_street_blocks(
                mock_street_config.geopackage_path,
                "perpendicular_lines",
                "on_grid_blocks",
                "secondary_setback",
                "arterial_setback",
            )

            mock_extract.assert_called_once()
            mock_create.assert_called_once()

    def test_generate_street_blocks_custom_layer_names(self, mock_street_config):
        """Test street blocks generation with custom layer names."""
        from rue_lib.streets.blocks import generate_street_blocks

        with (
            patch("rue_lib.streets.blocks.extract_arterial_edge_lines") as mock_extract,
            patch("rue_lib.streets.blocks.create_street_blocks") as mock_create,
        ):
            generate_street_blocks(
                mock_street_config.geopackage_path,
                "custom_perpendicular",
                "custom_on_grid",
                "custom_secondary",
                "custom_arterial",
                arterial_edges_layer_name="custom_arterial_edges",
                merged_lines_layer_name="custom_merged",
                street_blocks_layer_name="custom_blocks",
            )

            mock_extract.assert_called_once_with(
                mock_street_config.geopackage_path,
                "arterial_roads",  # Default arterial roads layer
                "custom_arterial",
                "custom_secondary",
                "custom_arterial_edges",
            )

            mock_create.assert_called_once_with(
                mock_street_config.geopackage_path,
                "custom_perpendicular",
                "custom_arterial_edges",
                "custom_on_grid",
                "custom_secondary",
                "custom_arterial",
                "custom_blocks",
            )


class TestFilterClassifiedBlocks:
    """Test filter_classified_blocks function."""

    def test_filter_classified_blocks(self):
        """Test filtering classified blocks."""
        from rue_lib.streets.blocks import filter_classified_blocks

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.blocks.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer

                # Mock features
                mock_feature = Mock()
                mock_feature.GetGeometryRef.return_value = Mock()
                mock_layer.__iter__ = Mock(return_value=iter([mock_feature]))

                filter_classified_blocks(gpkg_path, "input_blocks", "filtered_blocks")

                mock_ds.GetLayerByName.assert_called()


class TestGenerateOnGridBlocks:
    """Test generate_on_grid_blocks function."""

    def test_generate_on_grid_blocks(self):
        """Test generating on-grid blocks."""
        from rue_lib.streets.blocks import generate_on_grid_blocks

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.blocks.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer

                # Mock features
                mock_feature = Mock()
                mock_feature.GetGeometryRef.return_value = Mock()
                mock_layer.__iter__ = Mock(return_value=iter([mock_feature]))

                generate_on_grid_blocks(gpkg_path, "perpendicular_lines", "on_grid_blocks")

                mock_ds.GetLayerByName.assert_called()


class TestPolygonizeAndClassifyBlocks:
    """Test polygonize_and_classify_blocks function."""

    def test_polygonize_and_classify_blocks(self):
        """Test polygonizing and classifying blocks."""
        from rue_lib.streets.blocks import polygonize_and_classify_blocks

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.blocks.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer

                # Mock features
                mock_feature = Mock()
                mock_feature.GetGeometryRef.return_value = Mock()
                mock_layer.__iter__ = Mock(return_value=iter([mock_feature]))

                polygonize_and_classify_blocks(
                    gpkg_path, "merged_lines", "classified_blocks", "parcel_bounds"
                )

                mock_ds.GetLayerByName.assert_called()
