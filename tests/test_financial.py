"""Tests for financial module."""

import json
import tempfile
from pathlib import Path

from rue_lib.financial import FinancialModel


class TestFinancialModel:
    """Test FinancialModel class."""

    def test_data_in_json_empty(self):
        """Test data_in_json returns instance dictionary."""
        model = FinancialModel()
        result = model.data_in_json()
        assert isinstance(result, dict)
        # Should return empty dict for basic instance
        assert result == {}

    def test_data_in_json_with_attributes(self):
        """Test data_in_json returns all instance attributes."""
        model = FinancialModel()
        model.cost = 100.0
        model.revenue = 150.0
        model.description = "test model"

        result = model.data_in_json()
        assert result == {"cost": 100.0, "revenue": 150.0, "description": "test model"}

    def test_save_creates_directory(self):
        """Test save creates output directory if it doesn't exist."""
        model = FinancialModel()
        model.test_value = 42

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir) / "subdir" / "output"
            model.save(str(output_path))

            # Check directory was created
            assert output_path.exists()
            assert output_path.is_dir()

            # Check file was created
            json_file = output_path / "financial.json"
            assert json_file.exists()

    def test_save_writes_json_file(self):
        """Test save writes correct JSON content."""
        model = FinancialModel()
        model.total_cost = 1000.0
        model.units = 50
        model.cost_per_unit = 20.0

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir)
            model.save(str(output_path))

            json_file = output_path / "financial.json"
            assert json_file.exists()

            # Read and verify JSON content
            with json_file.open() as f:
                data = json.load(f)

            assert data == {"total_cost": 1000.0, "units": 50, "cost_per_unit": 20.0}

    def test_save_json_formatting(self):
        """Test save writes properly formatted JSON."""
        model = FinancialModel()
        model.nested = {"key": "value", "number": 123}
        model.list_data = [1, 2, 3]

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir)
            model.save(str(output_path))

            json_file = output_path / "financial.json"
            content = json_file.read_text()

            # Should be formatted with indent=2
            assert "{\n  " in content  # Check indentation
            assert '"nested":' in content
            assert '"list_data":' in content

    def test_save_overwrites_existing_file(self):
        """Test save overwrites existing financial.json file."""
        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir)
            json_file = output_path / "financial.json"

            # Create initial file
            json_file.write_text('{"old": "data"}')

            # Save new data
            model = FinancialModel()
            model.new_data = "updated"
            model.save(str(output_path))

            # Verify file was overwritten
            with json_file.open() as f:
                data = json.load(f)

            assert data == {"new_data": "updated"}
            assert "old" not in data

    def test_save_with_path_object(self):
        """Test save works with Path objects."""
        model = FinancialModel()
        model.path_test = True

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir)
            model.save(output_path)  # Pass Path object directly

            json_file = output_path / "financial.json"
            assert json_file.exists()

            with json_file.open() as f:
                data = json.load(f)

            assert data == {"path_test": True}

    def test_save_complex_data_types(self):
        """Test save handles various data types correctly."""
        model = FinancialModel()
        model.integer = 42
        model.float_val = 3.14159
        model.string = "test string"
        model.boolean = True
        model.none_val = None
        model.list_val = [1, "two", 3.0]
        model.dict_val = {"nested": {"deep": "value"}}

        with tempfile.TemporaryDirectory() as temp_dir:
            output_path = Path(temp_dir)
            model.save(str(output_path))

            json_file = output_path / "financial.json"
            with json_file.open() as f:
                data = json.load(f)

            assert data["integer"] == 42
            assert data["float_val"] == 3.14159
            assert data["string"] == "test string"
            assert data["boolean"] is True
            assert data["none_val"] is None
            assert data["list_val"] == [1, "two", 3.0]
            assert data["dict_val"] == {"nested": {"deep": "value"}}
