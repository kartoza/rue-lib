"""Tests for public config module."""

from pathlib import Path

from rue_lib.public.config import PublicConfig


class TestPublicConfig:
    """Test PublicConfig dataclass."""

    def test_public_config_required_fields(self):
        """Test PublicConfig creation with required fields."""
        config = PublicConfig(input_path="input.geojson")

        assert config.input_path == "input.geojson"

    def test_public_config_default_values(self):
        """Test PublicConfig default values."""
        config = PublicConfig(input_path="input.geojson")

        assert config.site_path == "outputs/step1_parcels/parcels.geojson"
        assert config.output_dir == "outputs/public"
        assert config.open_percent == 4.0
        assert config.amen_percent == 10.0

    def test_public_config_custom_values(self):
        """Test PublicConfig with custom values."""
        config = PublicConfig(
            input_path="/custom/input.geojson",
            site_path="/custom/site.geojson",
            output_dir="/custom/output",
            open_percent=6.0,
            amen_percent=15.0,
        )

        assert config.input_path == "/custom/input.geojson"
        assert config.site_path == "/custom/site.geojson"
        assert config.output_dir == "/custom/output"
        assert config.open_percent == 6.0
        assert config.amen_percent == 15.0

    def test_public_config_percentage_validation(self):
        """Test that percentage values are reasonable."""
        config = PublicConfig(input_path="input.geojson")

        # Percentages should be positive
        assert config.open_percent > 0
        assert config.amen_percent > 0

        # Percentages should be reasonable (not over 100%)
        assert config.open_percent <= 100
        assert config.amen_percent <= 100

    def test_public_config_percentage_relationship(self):
        """Test logical relationship between percentages."""
        config = PublicConfig(input_path="input.geojson")

        # Amenity percentage should typically be higher than open space
        # (amenities include more diverse facilities)
        assert config.amen_percent >= config.open_percent

    def test_public_config_file_paths_are_strings(self):
        """Test that file paths are strings."""
        config = PublicConfig(input_path="input.geojson")

        assert isinstance(config.input_path, str)
        assert isinstance(config.site_path, str)
        assert isinstance(config.output_dir, str)

    def test_public_config_numeric_types(self):
        """Test that numeric values have correct types."""
        config = PublicConfig(input_path="input.geojson")

        assert isinstance(config.open_percent, float)
        assert isinstance(config.amen_percent, float)

    def test_public_config_with_path_objects(self):
        """Test PublicConfig with Path objects."""
        input_path = Path("test_input.geojson")
        site_path = Path("test_site.geojson")
        output_dir = Path("test_output")

        config = PublicConfig(
            input_path=str(input_path), site_path=str(site_path), output_dir=str(output_dir)
        )

        # Should convert to strings
        assert config.input_path == "test_input.geojson"
        assert config.site_path == "test_site.geojson"
        assert config.output_dir == "test_output"

    def test_public_config_realistic_percentages(self):
        """Test that default percentages are realistic for urban planning."""
        config = PublicConfig(input_path="input.geojson")

        # Open space requirements are typically 2-10% for urban areas
        assert 2.0 <= config.open_percent <= 10.0

        # Amenity requirements are typically 5-20% for urban areas
        assert 5.0 <= config.amen_percent <= 20.0

    def test_public_config_output_paths(self):
        """Test output path construction."""
        config = PublicConfig(input_path="input.geojson", output_dir="custom_output")

        assert config.output_dir == "custom_output"

        # Verify default site path structure
        assert "parcels" in config.site_path
        assert config.site_path.endswith(".geojson")

    def test_public_config_validation_edge_cases(self):
        """Test edge cases for configuration validation."""
        # Minimum valid percentages
        config_min = PublicConfig(input_path="input.geojson", open_percent=0.1, amen_percent=0.1)
        assert config_min.open_percent == 0.1
        assert config_min.amen_percent == 0.1

        # High percentages (still valid)
        config_high = PublicConfig(input_path="input.geojson", open_percent=50.0, amen_percent=60.0)
        assert config_high.open_percent == 50.0
        assert config_high.amen_percent == 60.0

    def test_public_config_equality(self):
        """Test PublicConfig equality comparison."""
        config1 = PublicConfig(input_path="input.geojson", open_percent=5.0, amen_percent=12.0)

        config2 = PublicConfig(input_path="input.geojson", open_percent=5.0, amen_percent=12.0)

        config3 = PublicConfig(input_path="different.geojson", open_percent=5.0, amen_percent=12.0)

        assert config1 == config2
        assert config1 != config3

    def test_public_config_repr(self):
        """Test PublicConfig string representation."""
        config = PublicConfig(input_path="test.geojson")

        repr_str = repr(config)

        assert "PublicConfig" in repr_str
        assert "test.geojson" in repr_str
        assert "4.0" in repr_str  # open_percent
        assert "10.0" in repr_str  # amen_percent
