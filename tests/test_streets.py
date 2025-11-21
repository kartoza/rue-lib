"""Tests for streets operations."""

import tempfile
from pathlib import Path

import pytest
from osgeo import ogr

from rue_lib.streets.operations import create_local_streets_zone


class TestCreateLocalStreetsZone:
    """Tests for create_local_streets_zone function."""

    @pytest.fixture
    def temp_gpkg(self):
        """Create a temporary GeoPackage with test grid blocks."""
        with tempfile.TemporaryDirectory() as tmpdir:
            gpkg_path = str(Path(tmpdir) / "test.gpkg")

            # Create test GeoPackage with sample grid blocks
            driver = ogr.GetDriverByName("GPKG")
            ds = driver.CreateDataSource(gpkg_path)

            # Create spatial reference (UTM zone 33N)
            srs = ogr.osr.SpatialReference()
            srs.ImportFromEPSG(32633)

            # Create test layer with grid blocks
            layer = ds.CreateLayer("test_grid", srs, ogr.wkbPolygon)

            # Create a simple square grid block (100m x 100m)
            feature = ogr.Feature(layer.GetLayerDefn())
            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(0, 0)
            ring.AddPoint(100, 0)
            ring.AddPoint(100, 100)
            ring.AddPoint(0, 100)
            ring.AddPoint(0, 0)
            poly = ogr.Geometry(ogr.wkbPolygon)
            poly.AddGeometry(ring)
            feature.SetGeometry(poly)
            layer.CreateFeature(feature)

            # Create another grid block
            feature2 = ogr.Feature(layer.GetLayerDefn())
            ring2 = ogr.Geometry(ogr.wkbLinearRing)
            ring2.AddPoint(150, 0)
            ring2.AddPoint(250, 0)
            ring2.AddPoint(250, 100)
            ring2.AddPoint(150, 100)
            ring2.AddPoint(150, 0)
            poly2 = ogr.Geometry(ogr.wkbPolygon)
            poly2.AddGeometry(ring2)
            feature2.SetGeometry(poly2)
            layer.CreateFeature(feature2)

            ds = None

            yield gpkg_path

    def test_creates_both_layers(self, temp_gpkg):
        """Test that both inner and outer layers are created."""
        outer_layer, inner_layer = create_local_streets_zone(
            temp_gpkg,
            "test_grid",
            temp_gpkg,
            "local_streets",
            sidewalk_width_m=3.0,
            road_width_m=10.0,
        )

        assert outer_layer == "local_streets_outer"
        assert inner_layer == "local_streets_inner"

        # Verify layers exist in the GeoPackage
        ds = ogr.Open(temp_gpkg)
        assert ds.GetLayerByName(outer_layer) is not None
        assert ds.GetLayerByName(inner_layer) is not None
        ds = None

    def test_inner_area_smaller_than_outer(self, temp_gpkg):
        """Test that inner area is smaller than outer area."""
        create_local_streets_zone(
            temp_gpkg,
            "test_grid",
            temp_gpkg,
            "local_streets",
            sidewalk_width_m=3.0,
            road_width_m=10.0,
        )

        ds = ogr.Open(temp_gpkg)
        inner_layer = ds.GetLayerByName("local_streets_inner")
        outer_layer = ds.GetLayerByName("local_streets_outer")

        inner_feature = inner_layer.GetNextFeature()
        outer_feature = outer_layer.GetNextFeature()

        inner_area = inner_feature.GetGeometryRef().GetArea()
        outer_area = outer_feature.GetGeometryRef().GetArea()

        assert inner_area < outer_area
        ds = None

    def test_attributes_are_set(self, temp_gpkg):
        """Test that all attributes are correctly set."""
        create_local_streets_zone(
            temp_gpkg,
            "test_grid",
            temp_gpkg,
            "local_streets",
            sidewalk_width_m=3.0,
            road_width_m=10.0,
        )

        ds = ogr.Open(temp_gpkg)
        inner_layer = ds.GetLayerByName("local_streets_inner")
        inner_feature = inner_layer.GetNextFeature()

        # Check all fields exist and have values
        assert inner_feature.GetField("area_m2") > 0
        assert inner_feature.GetField("sidewalk_area") > 0
        assert inner_feature.GetField("buffer_dist") > 0
        assert inner_feature.GetField("sidewalk_w") == 3.0
        assert inner_feature.GetField("road_w") == 10.0
        assert inner_feature.GetField("zone_type") == "buildable"

        outer_layer = ds.GetLayerByName("local_streets_outer")
        outer_feature = outer_layer.GetNextFeature()

        assert outer_feature.GetField("area_m2") > 0
        assert outer_feature.GetField("sidewalk_area") > 0
        assert outer_feature.GetField("buffer_dist") > 0
        assert outer_feature.GetField("sidewalk_w") == 3.0
        assert outer_feature.GetField("road_w") == 10.0
        assert outer_feature.GetField("zone_type") == "street_sidewalk"

        ds = None

    def test_sidewalk_area_matches_difference(self, temp_gpkg):
        """Test that sidewalk_area equals outer_area - inner_area."""
        create_local_streets_zone(
            temp_gpkg,
            "test_grid",
            temp_gpkg,
            "local_streets",
            sidewalk_width_m=3.0,
            road_width_m=10.0,
        )

        ds = ogr.Open(temp_gpkg)
        inner_layer = ds.GetLayerByName("local_streets_inner")
        outer_layer = ds.GetLayerByName("local_streets_outer")

        inner_feature = inner_layer.GetNextFeature()
        outer_feature = outer_layer.GetNextFeature()

        inner_area = inner_feature.GetField("area_m2")
        outer_area = outer_feature.GetField("area_m2")
        sidewalk_area = inner_feature.GetField("sidewalk_area")

        # Should match within floating point precision
        assert abs((outer_area - inner_area) - sidewalk_area) < 0.01

        ds = None

    def test_geometries_are_dissolved(self, temp_gpkg):
        """Test that multiple input features are dissolved into one."""
        create_local_streets_zone(
            temp_gpkg,
            "test_grid",
            temp_gpkg,
            "local_streets",
            sidewalk_width_m=3.0,
            road_width_m=10.0,
        )

        ds = ogr.Open(temp_gpkg)
        inner_layer = ds.GetLayerByName("local_streets_inner")
        outer_layer = ds.GetLayerByName("local_streets_outer")

        # Should have exactly 1 feature (dissolved from 2 input features)
        assert inner_layer.GetFeatureCount() == 1
        assert outer_layer.GetFeatureCount() == 1

        ds = None

    def test_with_different_parameters(self, temp_gpkg):
        """Test with different sidewalk and road width parameters."""
        create_local_streets_zone(
            temp_gpkg,
            "test_grid",
            temp_gpkg,
            "local_streets_v2",
            sidewalk_width_m=5.0,
            road_width_m=15.0,
        )

        ds = ogr.Open(temp_gpkg)
        inner_layer = ds.GetLayerByName("local_streets_v2_inner")
        inner_feature = inner_layer.GetNextFeature()

        assert inner_feature.GetField("sidewalk_w") == 5.0
        assert inner_feature.GetField("road_w") == 15.0

        ds = None

    def test_layer_replacement(self, temp_gpkg):
        """Test that existing layers are replaced when running twice."""
        # Run once
        create_local_streets_zone(
            temp_gpkg,
            "test_grid",
            temp_gpkg,
            "local_streets",
            sidewalk_width_m=3.0,
            road_width_m=10.0,
        )

        # Run again with different parameters
        create_local_streets_zone(
            temp_gpkg,
            "test_grid",
            temp_gpkg,
            "local_streets",
            sidewalk_width_m=5.0,
            road_width_m=12.0,
        )

        ds = ogr.Open(temp_gpkg)
        inner_layer = ds.GetLayerByName("local_streets_inner")
        inner_feature = inner_layer.GetNextFeature()

        # Should have the new parameters
        assert inner_feature.GetField("sidewalk_w") == 5.0
        assert inner_feature.GetField("road_w") == 12.0

        # Should still have only 1 feature
        assert inner_layer.GetFeatureCount() == 1

        ds = None

    def test_geometries_are_valid(self, temp_gpkg):
        """Test that output geometries are valid."""
        create_local_streets_zone(
            temp_gpkg,
            "test_grid",
            temp_gpkg,
            "local_streets",
            sidewalk_width_m=3.0,
            road_width_m=10.0,
        )

        ds = ogr.Open(temp_gpkg)
        inner_layer = ds.GetLayerByName("local_streets_inner")
        outer_layer = ds.GetLayerByName("local_streets_outer")

        # Get features and clone their geometries before closing dataset
        inner_feature = inner_layer.GetNextFeature()
        outer_feature = outer_layer.GetNextFeature()

        assert inner_feature is not None
        assert outer_feature is not None

        inner_geom = inner_feature.GetGeometryRef().Clone()
        outer_geom = outer_feature.GetGeometryRef().Clone()

        ds = None

        # Validate after closing dataset
        assert inner_geom.IsValid()
        assert outer_geom.IsValid()
        assert not inner_geom.IsEmpty()
        assert not outer_geom.IsEmpty()
