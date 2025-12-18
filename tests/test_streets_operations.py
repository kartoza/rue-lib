"""Tests for streets operations module."""

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
    geom.GetGeometryType.return_value = ogr.wkbPolygon25D
    geom.GetGeometryCount.return_value = 2
    geom.GetGeometryRef.return_value = geom
    geom.ExportToWkb.return_value = b"\x01\x03\x00\x00\x00"
    geom.ExportToWkt.return_value = "POLYGON((0 0,1 0,1 1,0 1,0 0))"
    return geom


@pytest.fixture
def mock_shapely_polygon():
    """Create a mock shapely polygon."""
    from unittest.mock import Mock

    polygon = Mock()
    polygon.area = 1.0
    polygon.bounds = (0, 0, 1, 1)
    polygon.is_valid = True
    return polygon


class TestExplodeMultipolygon:
    """Test explode_multipolygon function."""

    def test_explode_single_polygon(self, mock_ogr_geometry):
        """Test exploding a single polygon."""
        from rue_lib.streets.operations import explode_multipolygon

        mock_ogr_geometry.GetGeometryType.return_value = ogr.wkbPolygon25D

        result = explode_multipolygon(mock_ogr_geometry)

        assert len(result) == 1
        assert result[0] == mock_ogr_geometry

    def test_explode_multipolygon(self, mock_ogr_geometry):
        """Test exploding a multipolygon."""
        from rue_lib.streets.operations import explode_multipolygon

        mock_ogr_geometry.GetGeometryType.return_value = ogr.wkbMultiPolygon25D
        mock_ogr_geometry.GetGeometryCount.return_value = 2

        sub_geom = Mock()
        sub_geom.GetGeometryType.return_value = ogr.wkbPolygon25D
        sub_geom.Clone.return_value = sub_geom
        mock_ogr_geometry.GetGeometryRef.return_value = sub_geom

        result = explode_multipolygon(mock_ogr_geometry)

        assert len(result) == 2

    def test_explode_2d_polygon(self, mock_ogr_geometry):
        """Test exploding 2D polygon (should still work)."""
        from rue_lib.streets.operations import explode_multipolygon

        mock_ogr_geometry.GetGeometryType.return_value = ogr.wkbPolygon

        result = explode_multipolygon(mock_ogr_geometry)

        assert len(result) == 1


class TestOgrToShapely:
    """Test _ogr_to_shapely conversion function."""

    def test_ogr_to_shapely_none_input(self):
        """Test conversion with None input."""
        from rue_lib.streets.operations import _ogr_to_shapely

        result = _ogr_to_shapely(None)
        assert result is None

    def test_ogr_to_shapely_wkb_success(self, mock_ogr_geometry):
        """Test successful WKB conversion."""
        from rue_lib.streets.operations import _ogr_to_shapely

        mock_ogr_geometry.ExportToWkb.return_value = b"\x01\x03\x00\x00\x00"

        with patch("rue_lib.streets.operations.wkb") as mock_wkb:
            mock_polygon = Mock()
            mock_wkb.loads.return_value = mock_polygon

            result = _ogr_to_shapely(mock_ogr_geometry)

            assert result == mock_polygon
            mock_wkb.loads.assert_called_once()

    def test_ogr_to_shapely_wkb_memoryview(self, mock_ogr_geometry):
        """Test WKB conversion with memoryview."""
        from rue_lib.streets.operations import _ogr_to_shapely

        # Mock memoryview WKB data
        mock_ogr_geometry.ExportToWkb.return_value = memoryview(b"\x01\x03\x00\x00\x00")

        with patch("rue_lib.streets.operations.wkb") as mock_wkb:
            mock_polygon = Mock()
            mock_wkb.loads.return_value = mock_polygon

            result = _ogr_to_shapely(mock_ogr_geometry)

            assert result == mock_polygon

    def test_ogr_to_shapely_fallback_to_wkt(self, mock_ogr_geometry):
        """Test fallback to WKT when WKB fails."""
        from rue_lib.streets.operations import _ogr_to_shapely

        mock_ogr_geometry.ExportToWkb.side_effect = TypeError("Invalid WKB")
        mock_ogr_geometry.ExportToWkt.return_value = "POLYGON((0 0,1 0,1 1,0 1,0 0))"

        with patch("rue_lib.streets.operations.wkt") as mock_wkt:
            mock_polygon = Mock()
            mock_wkt.loads.return_value = mock_polygon

            result = _ogr_to_shapely(mock_ogr_geometry)

            assert result == mock_polygon
            mock_wkt.loads.assert_called_once_with("POLYGON((0 0,1 0,1 1,0 1,0 0))")

    def test_ogr_to_shapely_both_fail(self, mock_ogr_geometry):
        """Test when both WKB and WKT conversion fail."""
        from rue_lib.streets.operations import _ogr_to_shapely

        mock_ogr_geometry.ExportToWkb.side_effect = TypeError("Invalid WKB")
        mock_ogr_geometry.ExportToWkt.side_effect = RuntimeError("Invalid WKT")

        result = _ogr_to_shapely(mock_ogr_geometry)

        assert result is None


class TestBreakMultipartFeatures:
    """Test break_multipart_features function."""

    def test_break_multipart_features(self):
        """Test breaking multipart features."""
        from rue_lib.streets.operations import break_multipart_features

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_layer = Mock()

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_layer

                # Mock multipart geometry
                mock_feature = Mock()
                mock_geom = Mock()
                mock_geom.GetGeometryType.return_value = ogr.wkbMultiPolygon25D
                mock_feature.GetGeometryRef.return_value = mock_geom
                mock_layer.__iter__ = Mock(return_value=iter([mock_feature]))

                with patch("rue_lib.streets.operations.explode_multipolygon") as mock_explode:
                    mock_explode.return_value = [mock_geom]

                    break_multipart_features(gpkg_path, "input", "output", ogr.wkbPolygon25D)

                    mock_explode.assert_called_once_with(mock_geom)


class TestClipLayer:
    """Test clip_layer function."""

    def test_clip_layer_basic(self):
        """Test basic layer clipping."""
        from rue_lib.streets.operations import clip_layer

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_source_layer = Mock()
                mock_clip_layer = Mock()

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.side_effect = [mock_source_layer, mock_clip_layer]

                # Mock features
                mock_feature = Mock()
                mock_feature.GetGeometryRef.return_value = Mock()
                mock_source_layer.__iter__ = Mock(return_value=iter([mock_feature]))
                mock_clip_layer.__iter__ = Mock(return_value=iter([mock_feature]))

                clip_layer(gpkg_path, "source", "clip", "output", ogr.wkbPolygon25D)

                mock_ds.GetLayerByName.assert_called()


class TestEraseLayer:
    """Test erase_layer function."""

    def test_erase_layer(self):
        """Test layer erase operation."""
        from rue_lib.streets.operations import erase_layer

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.streets.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_input_layer = Mock()
                mock_erase_layer = Mock()

                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.side_effect = [mock_input_layer, mock_erase_layer]

                erase_layer(gpkg_path, "input", "erase", "output", ogr.wkbPolygon25D)

                mock_ds.GetLayerByName.assert_called()
