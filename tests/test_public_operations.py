"""Tests for public operations module."""

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
    geom.GetEnvelope.return_value = (0, 10, 0, 10)
    geom.Area.return_value = 100.0
    geom.Centroid.return_value = geom
    geom.GetX.return_value = 5.0
    geom.GetY.return_value = 5.0
    geom.Distance.return_value = 10.0
    geom.Intersection.return_value = geom
    geom.Contains.return_value = True
    return geom


@pytest.fixture
def mock_ogr_feature():
    """Create a mock OGR feature."""
    feature = Mock()
    feature.GetGeometryRef.return_value = Mock()
    feature.GetFID.return_value = 1
    feature.GetFieldCount.return_value = 3
    feature.GetFieldDefnRef.return_value = Mock()
    feature.GetField.return_value = "test_value"
    feature.SetField = Mock()
    feature.Clone.return_value = feature
    return feature


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


class TestAllocateOpenSpaces:
    """Test allocate_open_spaces function."""

    def test_allocate_open_spaces_basic(self, mock_ogr_layer, mock_ogr_feature, mock_ogr_geometry):
        """Test basic open space allocation."""
        from rue_lib.public.operations import allocate_open_spaces

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.public.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.side_effect = [mock_ogr_layer, mock_ogr_layer]
                mock_ds.CreateLayer.return_value = mock_ogr_layer
                mock_ds.DeleteLayer = Mock()

                # Mock features for iteration
                mock_ogr_feature.GetGeometryRef.return_value = mock_ogr_geometry
                mock_ogr_layer.__iter__ = Mock(return_value=iter([mock_ogr_feature]))
                mock_ogr_layer.GetLayerDefn.return_value = Mock()

                result = allocate_open_spaces(
                    gpkg_path, "parcels", "blocks", "open_spaces", open_percent=4.0
                )

                assert result == "open_spaces"
                mock_ds.CreateLayer.assert_called_once()
                mock_ogr_layer.CreateFeature.assert_called()

    def test_allocate_open_spaces_custom_percentage(
        self, mock_ogr_layer, mock_ogr_feature, mock_ogr_geometry
    ):
        """Test open space allocation with custom percentage."""
        from rue_lib.public.operations import allocate_open_spaces

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.public.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.side_effect = [mock_ogr_layer, mock_ogr_layer]
                mock_ds.CreateLayer.return_value = mock_ogr_layer
                mock_ds.DeleteLayer = Mock()

                mock_ogr_feature.GetGeometryRef.return_value = mock_ogr_geometry
                mock_ogr_layer.__iter__ = Mock(return_value=iter([mock_ogr_feature]))
                mock_ogr_layer.GetLayerDefn.return_value = Mock()

                result = allocate_open_spaces(
                    gpkg_path, "parcels", "blocks", "open_spaces", open_percent=8.0
                )

                assert result == "open_spaces"

    def test_allocate_open_spaces_empty_layers(self, mock_ogr_layer):
        """Test open space allocation with empty input layers."""
        from rue_lib.public.operations import allocate_open_spaces

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.public.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.side_effect = [mock_ogr_layer, mock_ogr_layer]
                mock_ds.CreateLayer.return_value = mock_ogr_layer
                mock_ds.DeleteLayer = Mock()

                # Empty layers
                mock_ogr_layer.__iter__ = Mock(return_value=iter([]))
                mock_ogr_layer.GetLayerDefn.return_value = Mock()

                result = allocate_open_spaces(gpkg_path, "parcels", "blocks", "open_spaces")

                assert result == "open_spaces"
                mock_ds.CreateLayer.assert_called_once()

    def test_allocate_open_spaces_centroid_calculation(
        self, mock_ogr_layer, mock_ogr_feature, mock_ogr_geometry
    ):
        """Test centroid calculation in open space allocation."""
        from rue_lib.public.operations import allocate_open_spaces

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.public.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.side_effect = [mock_ogr_layer, mock_ogr_layer]
                mock_ds.CreateLayer.return_value = mock_ogr_layer
                mock_ds.DeleteLayer = Mock()

                # Mock multiple features for centroid calculation
                mock_feature1 = Mock()
                mock_geom1 = Mock()
                mock_geom1.Centroid.return_value = mock_geom1
                mock_geom1.GetX.return_value = 2.0
                mock_geom1.GetY.return_value = 2.0
                mock_feature1.GetGeometryRef.return_value = mock_geom1

                mock_feature2 = Mock()
                mock_geom2 = Mock()
                mock_geom2.Centroid.return_value = mock_geom2
                mock_geom2.GetX.return_value = 8.0
                mock_geom2.GetY.return_value = 8.0
                mock_feature2.GetGeometryRef.return_value = mock_geom2

                mock_ogr_layer.__iter__ = Mock(return_value=iter([mock_feature1, mock_feature2]))
                mock_ogr_layer.GetLayerDefn.return_value = Mock()

                result = allocate_open_spaces(gpkg_path, "parcels", "blocks", "open_spaces")

                # Should calculate centroids for distance sorting
                mock_geom1.Centroid.assert_called()
                mock_geom2.Centroid.assert_called()
                assert result == "open_spaces"


class TestAllocateAmenities:
    """Test allocate_amenities function."""

    def test_allocate_amenities_basic(self, mock_ogr_layer, mock_ogr_feature, mock_ogr_geometry):
        """Test basic amenity allocation."""
        from rue_lib.public.operations import allocate_amenities

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.public.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.side_effect = [mock_ogr_layer, mock_ogr_layer]
                mock_ds.CreateLayer.return_value = mock_ogr_layer
                mock_ds.DeleteLayer = Mock()

                mock_ogr_feature.GetGeometryRef.return_value = mock_ogr_geometry
                mock_ogr_layer.__iter__ = Mock(return_value=iter([mock_ogr_feature]))
                mock_ogr_layer.GetLayerDefn.return_value = Mock()

                result = allocate_amenities(
                    gpkg_path, "parcels", "blocks", "amenities", amen_percent=10.0
                )

                assert result == "amenities"
                mock_ds.CreateLayer.assert_called_once()
                mock_ogr_layer.CreateFeature.assert_called()

    def test_allocate_amenities_high_percentage(
        self, mock_ogr_layer, mock_ogr_feature, mock_ogr_geometry
    ):
        """Test amenity allocation with high percentage requirement."""
        from rue_lib.public.operations import allocate_amenities

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.public.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.side_effect = [mock_ogr_layer, mock_ogr_layer]
                mock_ds.CreateLayer.return_value = mock_ogr_layer
                mock_ds.DeleteLayer = Mock()

                mock_ogr_feature.GetGeometryRef.return_value = mock_ogr_geometry
                mock_ogr_layer.__iter__ = Mock(return_value=iter([mock_ogr_feature]))
                mock_ogr_layer.GetLayerDefn.return_value = Mock()

                result = allocate_amenities(
                    gpkg_path, "parcels", "blocks", "amenities", amen_percent=25.0
                )

                assert result == "amenities"

    def test_allocate_amenities_area_prioritization(self, mock_ogr_layer, mock_ogr_geometry):
        """Test amenity allocation prioritizes larger areas."""
        from rue_lib.public.operations import allocate_amenities

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.public.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.side_effect = [mock_ogr_layer, mock_ogr_layer]
                mock_ds.CreateLayer.return_value = mock_ogr_layer
                mock_ds.DeleteLayer = Mock()

                # Mock features with different areas
                mock_feature_small = Mock()
                mock_geom_small = Mock()
                mock_geom_small.Area.return_value = 50.0
                mock_feature_small.GetGeometryRef.return_value = mock_geom_small

                mock_feature_large = Mock()
                mock_geom_large = Mock()
                mock_geom_large.Area.return_value = 500.0
                mock_feature_large.GetGeometryRef.return_value = mock_geom_large

                mock_ogr_layer.__iter__ = Mock(
                    return_value=iter([mock_feature_small, mock_feature_large])
                )
                mock_ogr_layer.GetLayerDefn.return_value = Mock()

                result = allocate_amenities(gpkg_path, "parcels", "blocks", "amenities")

                # Should prioritize larger areas
                assert result == "amenities"
                mock_ogr_layer.CreateFeature.assert_called()


class TestMergePublicSpaces:
    """Test merge function for public spaces."""

    def test_merge_public_spaces_basic(self, mock_ogr_layer, mock_ogr_feature, mock_ogr_geometry):
        """Test basic merging of public spaces."""
        from rue_lib.public.operations import merge

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.public.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_ogr_layer
                mock_ds.CreateLayer.return_value = mock_ogr_layer
                mock_ds.DeleteLayer = Mock()

                mock_ogr_feature.GetGeometryRef.return_value = mock_ogr_geometry
                mock_ogr_layer.__iter__ = Mock(return_value=iter([mock_ogr_feature]))
                mock_ogr_layer.GetLayerDefn.return_value = Mock()

                result = merge(gpkg_path, ["open_spaces", "amenities"], "parcels", "public_spaces")

                assert result == "public_spaces"
                mock_ds.CreateLayer.assert_called()
                mock_ogr_layer.CreateFeature.assert_called()

    def test_merge_public_spaces_multiple_layers(
        self, mock_ogr_layer, mock_ogr_feature, mock_ogr_geometry
    ):
        """Test merging multiple public space layers."""
        from rue_lib.public.operations import merge

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.public.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_ogr_layer
                mock_ds.CreateLayer.return_value = mock_ogr_layer
                mock_ds.DeleteLayer = Mock()

                mock_ogr_feature.GetGeometryRef.return_value = mock_ogr_geometry
                mock_ogr_layer.__iter__ = Mock(return_value=iter([mock_ogr_feature]))
                mock_ogr_layer.GetLayerDefn.return_value = Mock()

                result = merge(
                    gpkg_path,
                    ["open_spaces", "amenities", "playgrounds", "community_centers"],
                    "parcels",
                    "all_public_spaces",
                )

                assert result == "all_public_spaces"
                # Should iterate through all input layers
                assert mock_ds.GetLayerByName.call_count >= 4

    def test_merge_public_spaces_empty_layers(self, mock_ogr_layer):
        """Test merging with empty input layers."""
        from rue_lib.public.operations import merge

        with tempfile.NamedTemporaryFile(suffix=".gpkg") as tmp_file:
            gpkg_path = tmp_file.name

            with patch("rue_lib.public.operations.ogr") as mock_ogr:
                mock_ds = Mock()
                mock_ogr.Open.return_value = mock_ds
                mock_ds.GetLayerByName.return_value = mock_ogr_layer
                mock_ds.CreateLayer.return_value = mock_ogr_layer
                mock_ds.DeleteLayer = Mock()

                # Empty layers
                mock_ogr_layer.__iter__ = Mock(return_value=iter([]))
                mock_ogr_layer.GetLayerDefn.return_value = Mock()

                result = merge(gpkg_path, ["open_spaces", "amenities"], "parcels", "public_spaces")

                assert result == "public_spaces"
                mock_ds.CreateLayer.assert_called_once()
