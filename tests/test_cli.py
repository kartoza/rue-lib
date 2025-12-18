"""Tests for CLI modules."""

import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

from rue_lib.cli.step1 import main


class TestStep1CLI:
    """Test the step1 CLI functionality."""

    def test_main_missing_site_file(self, capsys):
        """Test CLI exits when site file is missing."""
        with patch("sys.argv", ["step1", "--site", "nonexistent.geojson"]):
            with pytest.raises(SystemExit) as excinfo:
                main()

            assert excinfo.value.code == 1
            captured = capsys.readouterr()
            assert "Site file not found" in captured.err

    def test_main_missing_roads_file(self, capsys):
        """Test CLI exits when roads file is missing."""
        with tempfile.NamedTemporaryFile(suffix=".geojson", delete=False) as site_file:
            site_file.write(b'{"type": "FeatureCollection", "features": []}')
            site_path = site_file.name

        try:
            with patch(
                "sys.argv", ["step1", "--site", site_path, "--roads", "nonexistent.geojson"]
            ):
                with pytest.raises(SystemExit) as excinfo:
                    main()

                assert excinfo.value.code == 1
                captured = capsys.readouterr()
                assert "Roads file not found" in captured.err
        finally:
            Path(site_path).unlink()

    @patch("rue_lib.cli.step1.generate_parcels")
    def test_main_success(self, mock_generate, capsys):
        """Test successful CLI execution."""
        with tempfile.NamedTemporaryFile(suffix=".geojson", delete=False) as site_file:
            site_file.write(b'{"type": "FeatureCollection", "features": []}')
            site_path = site_file.name

        with tempfile.NamedTemporaryFile(suffix=".geojson", delete=False) as roads_file:
            roads_file.write(b'{"type": "FeatureCollection", "features": []}')
            roads_path = roads_file.name

        try:
            mock_generate.return_value = "test_output.geojson"

            with patch(
                "sys.argv",
                [
                    "step1",
                    "--site",
                    site_path,
                    "--roads",
                    roads_path,
                    "--rows",
                    "3",
                    "--cols",
                    "3",
                    "--pad",
                    "15.0",
                    "--trim-by-roads",
                    "--min_owner_area",
                    "10.0",
                    "--out",
                    "test_output",
                ],
            ):
                main()

            captured = capsys.readouterr()
            assert "Final landowners GeoJSON: test_output.geojson" in captured.out

            # Verify generate_parcels was called with correct config
            mock_generate.assert_called_once()
            config = mock_generate.call_args[0][0]
            assert config.site_path == site_path
            assert config.roads_path == roads_path
            assert config.rows == 3
            assert config.cols == 3
            assert config.pad_m == 15.0
            assert config.output_dir == "test_output"

        finally:
            Path(site_path).unlink()
            Path(roads_path).unlink()

    def test_main_default_arguments(self):
        """Test CLI with default arguments."""
        with (
            patch("rue_lib.cli.step1.Path.cwd") as mock_cwd,
            patch("rue_lib.cli.step1.Path.exists") as mock_exists,
            patch("rue_lib.cli.step1.generate_parcels") as mock_generate,
            patch("sys.argv", ["step1"]),
        ):
            mock_cwd.return_value = Path("/test/dir")
            mock_exists.return_value = True
            mock_generate.return_value = "output.geojson"

            main()

            # Verify default paths were used
            config = mock_generate.call_args[0][0]
            assert config.site_path == "/test/dir/site.geojson"
            assert config.roads_path == "/test/dir/roads.geojson"
            assert config.rows == 4
            assert config.cols == 4
            assert config.pad_m == 10.0
            assert config.output_dir == "outputs/site_mob_phase1"
