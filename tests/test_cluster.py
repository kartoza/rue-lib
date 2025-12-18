"""Tests for cluster modules."""

from pathlib import Path

from rue_lib.cluster.config import ClusterConfig


class TestClusterConfig:
    """Test ClusterConfig dataclass."""

    def test_cluster_config_required_fields(self):
        """Test ClusterConfig creation with required fields."""
        config = ClusterConfig(roads_path="roads.geojson", input_path="input.geojson")

        assert config.roads_path == "roads.geojson"
        assert config.input_path == "input.geojson"

    def test_cluster_config_default_values(self):
        """Test ClusterConfig default values."""
        config = ClusterConfig(roads_path="roads.geojson", input_path="input.geojson")

        # Output settings
        assert config.output_dir == "outputs"
        assert config.geopackage_path == "outputs/output.gpkg"

        # Thresholds
        assert config.off_grid_plot_threshold == 0.5

        # Road widths
        assert config.road_arterial_width_m == 20.0
        assert config.road_secondary_width_m == 15.0
        assert config.road_local_width_m == 10.0

        # On-grid partition depths
        assert config.on_grid_partition_depth_arterial_roads == 40.0
        assert config.on_grid_partition_depth_secondary_roads == 30.0
        assert config.on_grid_partition_depth_local_roads == 20.0

        # Off-grid cluster dimensions
        assert config.off_grid_cluster_depth == 45.0
        assert config.off_grid_cluster_width == 30.0

        # Infrastructure
        assert config.sidewalk_width_m == 3.0

    def test_cluster_config_custom_values(self):
        """Test ClusterConfig with custom values."""
        config = ClusterConfig(
            roads_path="/custom/roads.geojson",
            input_path="/custom/input.geojson",
            output_dir="/custom/output",
            geopackage_path="/custom/data.gpkg",
            off_grid_plot_threshold=0.6,
            road_arterial_width_m=25.0,
            road_secondary_width_m=18.0,
            road_local_width_m=12.0,
            on_grid_partition_depth_arterial_roads=50.0,
            on_grid_partition_depth_secondary_roads=35.0,
            on_grid_partition_depth_local_roads=25.0,
            off_grid_cluster_depth=50.0,
            off_grid_cluster_width=35.0,
            sidewalk_width_m=4.0,
        )

        assert config.roads_path == "/custom/roads.geojson"
        assert config.input_path == "/custom/input.geojson"
        assert config.output_dir == "/custom/output"
        assert config.geopackage_path == "/custom/data.gpkg"
        assert config.off_grid_plot_threshold == 0.6
        assert config.road_arterial_width_m == 25.0
        assert config.road_secondary_width_m == 18.0
        assert config.road_local_width_m == 12.0
        assert config.on_grid_partition_depth_arterial_roads == 50.0
        assert config.on_grid_partition_depth_secondary_roads == 35.0
        assert config.on_grid_partition_depth_local_roads == 25.0
        assert config.off_grid_cluster_depth == 50.0
        assert config.off_grid_cluster_width == 35.0
        assert config.sidewalk_width_m == 4.0

    def test_cluster_config_road_width_hierarchy(self):
        """Test logical hierarchy of road widths."""
        config = ClusterConfig(roads_path="roads.geojson", input_path="input.geojson")

        # Arterial should be widest
        assert config.road_arterial_width_m > config.road_secondary_width_m
        assert config.road_arterial_width_m > config.road_local_width_m

        # Secondary should be wider than local
        assert config.road_secondary_width_m > config.road_local_width_m

    def test_cluster_config_partition_depth_hierarchy(self):
        """Test logical hierarchy of partition depths."""
        config = ClusterConfig(roads_path="roads.geojson", input_path="input.geojson")

        # Arterial partitions should be deepest
        assert (
            config.on_grid_partition_depth_arterial_roads
            > config.on_grid_partition_depth_secondary_roads
        )
        assert (
            config.on_grid_partition_depth_arterial_roads
            > config.on_grid_partition_depth_local_roads
        )

        # Secondary should be deeper than local
        assert (
            config.on_grid_partition_depth_secondary_roads
            > config.on_grid_partition_depth_local_roads
        )

    def test_cluster_config_positive_values(self):
        """Test that all numeric values are positive."""
        config = ClusterConfig(roads_path="roads.geojson", input_path="input.geojson")

        # Road widths should be positive
        assert config.road_arterial_width_m > 0
        assert config.road_secondary_width_m > 0
        assert config.road_local_width_m > 0
        assert config.sidewalk_width_m > 0

        # Partition depths should be positive
        assert config.on_grid_partition_depth_arterial_roads > 0
        assert config.on_grid_partition_depth_secondary_roads > 0
        assert config.on_grid_partition_depth_local_roads > 0

        # Cluster dimensions should be positive
        assert config.off_grid_cluster_depth > 0
        assert config.off_grid_cluster_width > 0

    def test_cluster_config_threshold_range(self):
        """Test off-grid plot threshold is in valid range."""
        config = ClusterConfig(roads_path="roads.geojson", input_path="input.geojson")

        # Threshold should be between 0 and 1 (percentage)
        assert 0 <= config.off_grid_plot_threshold <= 1

    def test_cluster_config_file_paths_are_strings(self):
        """Test that file paths are strings."""
        config = ClusterConfig(roads_path="roads.geojson", input_path="input.geojson")

        assert isinstance(config.roads_path, str)
        assert isinstance(config.input_path, str)
        assert isinstance(config.output_dir, str)
        assert isinstance(config.geopackage_path, str)

    def test_cluster_config_numeric_types(self):
        """Test that numeric values have correct types."""
        config = ClusterConfig(roads_path="roads.geojson", input_path="input.geojson")

        # Should all be floats
        assert isinstance(config.off_grid_plot_threshold, float)
        assert isinstance(config.road_arterial_width_m, float)
        assert isinstance(config.road_secondary_width_m, float)
        assert isinstance(config.road_local_width_m, float)
        assert isinstance(config.on_grid_partition_depth_arterial_roads, float)
        assert isinstance(config.on_grid_partition_depth_secondary_roads, float)
        assert isinstance(config.on_grid_partition_depth_local_roads, float)
        assert isinstance(config.off_grid_cluster_depth, float)
        assert isinstance(config.off_grid_cluster_width, float)
        assert isinstance(config.sidewalk_width_m, float)

    def test_cluster_config_with_path_objects(self):
        """Test ClusterConfig with Path objects."""
        roads_path = Path("test_roads.geojson")
        input_path = Path("test_input.geojson")

        config = ClusterConfig(roads_path=str(roads_path), input_path=str(input_path))

        # Should convert to strings
        assert config.roads_path == "test_roads.geojson"
        assert config.input_path == "test_input.geojson"

    def test_cluster_config_realistic_dimensions(self):
        """Test that default dimensions are realistic for urban planning."""
        config = ClusterConfig(roads_path="roads.geojson", input_path="input.geojson")

        # Road widths should be reasonable for urban context (5-30m)
        assert 5 <= config.road_local_width_m <= 30
        assert 5 <= config.road_secondary_width_m <= 30
        assert 5 <= config.road_arterial_width_m <= 30

        # Sidewalks should be reasonable (1-5m)
        assert 1 <= config.sidewalk_width_m <= 5

        # Partition depths should be reasonable for city blocks (10-100m)
        assert 10 <= config.on_grid_partition_depth_local_roads <= 100
        assert 10 <= config.on_grid_partition_depth_secondary_roads <= 100
        assert 10 <= config.on_grid_partition_depth_arterial_roads <= 100

        # Cluster dimensions should be reasonable (20-60m)
        assert 20 <= config.off_grid_cluster_depth <= 60
        assert 20 <= config.off_grid_cluster_width <= 60
