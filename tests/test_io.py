"""Tests for I/O utilities."""

import json
import tempfile
from pathlib import Path

import geopandas as gpd
import pytest
from shapely.geometry import Point

from rue_lib.utils.io import (
    load_geojson,
    prepare_geopackage,
    save_geojson,
    save_json,
    to_feature,
)


class TestPrepareGeopackage:
    """Tests for prepare_geopackage function."""

    def test_creates_geopackage_from_template_when_not_exists(self):
        """Test that geopackage is created from template when it doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            output_path = tmpdir_path / "test_output.gpkg"

            # Ensure output doesn't exist initially
            assert not output_path.exists()

            # Call prepare_geopackage
            result_path = prepare_geopackage(output_path)

            # Verify geopackage was created
            assert result_path.exists()
            assert result_path == output_path
            assert result_path.suffix == ".gpkg"

    def test_returns_existing_geopackage_when_exists(self):
        """Test that existing geopackage is returned without modification."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            output_path = tmpdir_path / "existing.gpkg"

            # Create an existing geopackage with custom content
            test_gdf = gpd.GeoDataFrame(
                {"id": [1], "name": ["test"]}, geometry=[Point(0, 0)], crs="EPSG:4326"
            )
            test_gdf.to_file(output_path, driver="GPKG")

            original_size = output_path.stat().st_size

            # Call prepare_geopackage
            result_path = prepare_geopackage(output_path)

            # Verify original file is returned unchanged
            assert result_path == output_path
            assert result_path.exists()
            assert output_path.stat().st_size == original_size

    def test_creates_parent_directories(self):
        """Test that parent directories are created when they don't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            output_path = tmpdir_path / "nested" / "deep" / "test.gpkg"

            # Ensure parent directories don't exist
            assert not output_path.parent.exists()

            # Call prepare_geopackage
            result_path = prepare_geopackage(output_path)

            # Verify directories were created
            assert output_path.parent.exists()
            assert result_path.exists()

    def test_custom_template_path(self):
        """Test using a custom template path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create a custom template
            template_path = tmpdir_path / "custom_template.gpkg"
            custom_gdf = gpd.GeoDataFrame(
                {"custom_field": [1, 2], "value": ["A", "B"]},
                geometry=[Point(1, 1), Point(2, 2)],
                crs="EPSG:4326",
            )
            custom_gdf.to_file(template_path, driver="GPKG")

            # Use custom template
            output_path = tmpdir_path / "output.gpkg"
            result_path = prepare_geopackage(output_path, template_path)

            # Verify custom template was used
            assert result_path.exists()
            result_gdf = gpd.read_file(result_path)
            assert "custom_field" in result_gdf.columns
            assert len(result_gdf) == 2

    def test_template_not_found_raises_error(self):
        """Test that FileNotFoundError is raised when template doesn't exist."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)
            output_path = tmpdir_path / "test.gpkg"
            nonexistent_template = tmpdir_path / "missing_template.gpkg"

            with pytest.raises(FileNotFoundError, match="Template geopackage not found"):
                prepare_geopackage(output_path, nonexistent_template)

    def test_string_paths_accepted(self):
        """Test that string paths are accepted and converted to Path objects."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path_str = str(Path(tmpdir) / "test.gpkg")

            result_path = prepare_geopackage(output_path_str)

            assert isinstance(result_path, Path)
            assert result_path.exists()

    def test_built_in_template_exists(self):
        """Test that the built-in template exists and is valid."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.gpkg"

            # This should not raise an error if built-in template exists
            result_path = prepare_geopackage(output_path)
            assert result_path.exists()

            # Verify it's a valid geopackage by reading with geopandas
            # Should not raise an error if it's a valid geopackage
            try:
                gpd.read_file(result_path)
            except Exception:
                # If reading fails, check what layers are available
                import fiona

                layers = fiona.listlayers(str(result_path))
                assert len(layers) > 0  # Should have at least some layers


class TestGeojsonFunctions:
    """Tests for GeoJSON utility functions."""

    def test_load_geojson(self):
        """Test loading GeoJSON files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            geojson_path = Path(tmpdir) / "test.geojson"
            test_data = {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "properties": {"name": "test"},
                        "geometry": {"type": "Point", "coordinates": [0, 0]},
                    }
                ],
            }

            with geojson_path.open("w") as f:
                json.dump(test_data, f)

            loaded_data = load_geojson(geojson_path)
            assert loaded_data == test_data

    def test_save_geojson(self):
        """Test saving features to GeoJSON."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.geojson"
            features = [
                {
                    "type": "Feature",
                    "properties": {"id": 1},
                    "geometry": {"type": "Point", "coordinates": [1, 1]},
                }
            ]

            save_geojson(output_path, features)

            assert output_path.exists()
            with output_path.open("r") as f:
                saved_data = json.load(f)

            assert saved_data["type"] == "FeatureCollection"
            assert saved_data["features"] == features

    def test_save_json(self):
        """Test saving JSON data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "test.json"
            test_data = {"key": "value", "number": 42}

            save_json(output_path, test_data)

            assert output_path.exists()
            with output_path.open("r") as f:
                loaded_data = json.load(f)

            assert loaded_data == test_data

    def test_to_feature(self):
        """Test converting Shapely geometry to GeoJSON feature."""
        point = Point(1, 2)
        props = {"name": "test_point", "value": 42}

        feature = to_feature(point, props)

        assert feature["type"] == "Feature"
        assert feature["properties"] == props
        assert feature["geometry"]["type"] == "Point"
        assert feature["geometry"]["coordinates"] == [1.0, 2.0]

    def test_to_feature_no_properties(self):
        """Test converting Shapely geometry to GeoJSON feature without properties."""
        point = Point(3, 4)

        feature = to_feature(point)

        assert feature["type"] == "Feature"
        assert feature["properties"] == {}
        assert feature["geometry"]["type"] == "Point"
        assert feature["geometry"]["coordinates"] == [3.0, 4.0]
