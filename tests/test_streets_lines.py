"""Tests for streets lines module."""

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
def mock_ogr_geometry():
    """Create a mock OGR geometry."""
    geom = Mock()
    geom.Clone.return_value = geom
    geom.GetGeometryType.return_value = ogr.wkbLineString25D
    geom.GetPointCount.return_value = 5
    geom.GetPoint.side_effect = [(0, 0, 0), (1, 1, 0), (2, 0, 0), (3, 1, 0), (4, 0, 0)]
    geom.ExportToWkt.return_value = "LINESTRING(0 0,1 1,2 0,3 1,4 0)"
    geom.Length.return_value = 4.0
    geom.Intersection.return_value = geom
    geom.Difference.return_value = geom
    geom.Union.return_value = geom
    geom.Buffer.return_value = geom
    return geom


@pytest.fixture
def mock_ogr_polygon():
    """Create a mock OGR polygon geometry."""
    geom = Mock()
    geom.Clone.return_value = geom
    geom.GetGeometryType.return_value = ogr.wkbPolygon25D
    geom.GetEnvelope.return_value = (0, 10, 0, 10)
    geom.Area.return_value = 100.0
    geom.Intersection.return_value = geom
    geom.Contains.return_value = True
    geom.Intersects.return_value = True
    return geom


class TestExtractArterialEdgeLines:
    """Test extract_arterial_edge_lines function."""

    def test_extract_arterial_edge_lines_basic(self):
        """Test basic arterial edge line extraction."""
        from rue_lib.streets.lines import extract_arterial_edge_lines

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.lines.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_arterial_layer = Mock()
                mock_arterial_setback_layer = Mock()
                mock_secondary_setback_layer = Mock()
                mock_output_layer = Mock()

                # Mock arterial roads features
                mock_arterial_feature = Mock()
                mock_arterial_geom = Mock()
                mock_arterial_feature.GetGeometryRef.return_value = mock_arterial_geom
                mock_arterial_layer.__iter__ = Mock(return_value=iter([mock_arterial_feature]))

                # Mock setback features
                mock_setback_feature = Mock()
                mock_setback_geom = Mock()
                mock_setback_feature.GetGeometryRef.return_value = mock_setback_geom
                mock_arterial_setback_layer.__iter__ = Mock(
                    return_value=iter([mock_setback_feature])
                )
                mock_secondary_setback_layer.__iter__ = Mock(
                    return_value=iter([mock_setback_feature])
                )

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.side_effect = [
                    mock_arterial_layer,
                    mock_arterial_setback_layer,
                    mock_secondary_setback_layer,
                ]
                mock_ds.CreateLayer.return_value = mock_output_layer

                # Mock union operations
                mock_arterial_geom.Union.return_value = mock_arterial_geom
                mock_setback_geom.Union.return_value = mock_setback_geom
                mock_setback_geom.Difference.return_value = mock_setback_geom

                extract_arterial_edge_lines(
                    gpkg_path,
                    "arterial_roads",
                    "arterial_setback",
                    "secondary_setback",
                    "arterial_edges",
                )

                mock_ds.CreateLayer.assert_called_once()
                mock_arterial_geom.Union.assert_called()
                mock_setback_geom.Union.assert_called()

    def test_extract_arterial_edge_lines_custom_parameters(self):
        """Test arterial edge extraction with custom parameters."""
        from rue_lib.streets.lines import extract_arterial_edge_lines

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.lines.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()
                mock_feature = Mock()
                mock_geom = Mock()

                mock_feature.GetGeometryRef.return_value = mock_geom
                mock_layer.__iter__ = Mock(return_value=iter([mock_feature]))
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer
                mock_ds.CreateLayer.return_value = mock_layer

                mock_geom.Union.return_value = mock_geom
                mock_geom.Difference.return_value = mock_geom

                extract_arterial_edge_lines(
                    gpkg_path,
                    "arterial_roads",
                    "arterial_setback",
                    "secondary_setback",
                    "arterial_edges",
                    clip_buffer=0.5,
                    sample_distance=10.0,
                )

                mock_ds.CreateLayer.assert_called()

    def test_extract_arterial_edge_lines_empty_layers(self):
        """Test arterial edge extraction with empty input layers."""
        from rue_lib.streets.lines import extract_arterial_edge_lines

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.lines.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()

                # Empty layers
                mock_layer.__iter__ = Mock(return_value=iter([]))
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer
                mock_ds.CreateLayer.return_value = mock_layer

                extract_arterial_edge_lines(
                    gpkg_path,
                    "arterial_roads",
                    "arterial_setback",
                    "secondary_setback",
                    "arterial_edges",
                )

                mock_ds.CreateLayer.assert_called_once()
                # Should handle empty layers without crashing

    def test_extract_arterial_edge_lines_geometry_operations(self):
        """Test geometry operations in arterial edge extraction."""
        from rue_lib.streets.lines import extract_arterial_edge_lines

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.lines.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()
                mock_output_layer = Mock()

                # Mock complex geometry operations
                mock_feature = Mock()
                mock_geom = Mock()
                mock_union_geom = Mock()
                mock_diff_geom = Mock()

                mock_feature.GetGeometryRef.return_value = mock_geom
                mock_layer.__iter__ = Mock(return_value=iter([mock_feature]))
                mock_geom.Union.return_value = mock_union_geom
                mock_union_geom.Difference.return_value = mock_diff_geom

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer
                mock_ds.CreateLayer.return_value = mock_output_layer

                extract_arterial_edge_lines(
                    gpkg_path,
                    "arterial_roads",
                    "arterial_setback",
                    "secondary_setback",
                    "arterial_edges",
                )

                # Verify geometry operations were called
                mock_geom.Union.assert_called()
                mock_union_geom.Difference.assert_called()


class TestCreateDivisionPoints:
    """Test create_division_points function."""

    def test_create_division_points(self):
        """Test creating division points."""
        from rue_lib.streets.lines import create_division_points

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.lines.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer

                # Mock features
                mock_feature = Mock()
                mock_feature.GetGeometryRef.return_value = Mock()
                mock_layer.__iter__ = Mock(return_value=iter([mock_feature]))

                create_division_points(
                    gpkg_path, "arterial_roads", "division_points", distance=100.0
                )

                mock_ds.GetLayerByName.assert_called()


class TestCreatePerpendicularLines:
    """Test create_perpendicular_lines function."""

    def test_create_perpendicular_lines(self):
        """Test creating perpendicular lines."""
        from rue_lib.streets.lines import create_perpendicular_lines

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.lines.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer

                # Mock features
                mock_feature = Mock()
                mock_feature.GetGeometryRef.return_value = Mock()
                mock_layer.__iter__ = Mock(return_value=iter([mock_feature]))

                create_perpendicular_lines(
                    gpkg_path, "division_points", "perpendicular_lines", line_length=200.0
                )

                mock_ds.GetLayerByName.assert_called()
