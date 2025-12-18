"""Tests for streets config module."""

from rue_lib.streets.config import StreetConfig


class TestStreetConfig:
    """Test StreetConfig dataclass."""

    def test_street_config_required_fields(self):
        """Test StreetConfig creation with required fields."""
        config = StreetConfig(parcel_path="parcels.geojson", roads_path="roads.geojson")

        assert config.parcel_path == "parcels.geojson"
        assert config.roads_path == "roads.geojson"

    def test_street_config_default_values(self):
        """Test StreetConfig default values."""
        config = StreetConfig(parcel_path="parcels.geojson", roads_path="roads.geojson")

        # Basic defaults
        assert config.output_dir == "outputs"
        assert config.geopackage_path == "outputs/output.gpkg"
        assert config.road_arterial_width_m == 20.0
        assert config.road_secondary_width_m == 15.0

        # Grid configuration defaults
        assert config.on_grid_partition_depth_arterial_roads == 60.0
        assert config.on_grid_partition_depth_secondary_roads == 60.0
        assert config.off_grid_partitions_preferred_depth == 140.0
        assert config.off_grid_partitions_preferred_width == 140.0

        # Line generation defaults
        assert config.perpendicular_line_length == 1000.0
        assert config.optimize_grid_rotation is True
        assert config.grid_rotation_angle_step == 5.0
        assert config.use_ternary_search is False
        assert config.clip_to_boundary is True

        # Tolerance defaults
        assert config.tolerance_area_ratio == 0.70
        assert config.tolerance_boundary_distance == 10.0

        # Infrastructure defaults
        assert config.sidewalk_width_m == 3.0
        assert config.road_locals_width_m == 10.0

    def test_street_config_cluster_geometry_defaults(self):
        """Test cluster geometry configuration defaults."""
        config = StreetConfig(parcel_path="parcels.geojson", roads_path="roads.geojson")

        # On-grid cluster depths
        assert config.part_art_d == 40.0
        assert config.part_sec_d == 30.0
        assert config.part_loc_d == 20.0

        # Off-grid cluster dimensions
        assert config.off_grid_cluster_depth == 45.0
        assert config.off_grid_cluster_width == 30.0

    def test_street_config_cluster_counts_defaults(self):
        """Test cluster count configuration defaults."""
        config = StreetConfig(parcel_path="parcels.geojson", roads_path="roads.geojson")

        # Off-grid cluster counts
        assert config.off_grid_arterial_clusters_depth == 0
        assert config.off_grid_secondary_clusters_depth == 0
        assert config.off_grid_local_clusters_depth == 2
        assert config.off_grid_local_clusters_width == 3

    def test_street_config_custom_values(self):
        """Test StreetConfig with custom values."""
        config = StreetConfig(
            parcel_path="/custom/parcels.geojson",
            roads_path="/custom/roads.geojson",
            output_dir="/custom/output",
            geopackage_path="/custom/data.gpkg",
            road_arterial_width_m=25.0,
            road_secondary_width_m=18.0,
            on_grid_partition_depth_arterial_roads=70.0,
            on_grid_partition_depth_secondary_roads=65.0,
            off_grid_partitions_preferred_depth=150.0,
            off_grid_partitions_preferred_width=160.0,
            perpendicular_line_length=1200.0,
            optimize_grid_rotation=False,
            grid_rotation_angle_step=10.0,
            use_ternary_search=True,
            clip_to_boundary=False,
            tolerance_area_ratio=0.80,
            tolerance_boundary_distance=15.0,
            sidewalk_width_m=4.0,
            road_locals_width_m=12.0,
            part_art_d=50.0,
            part_sec_d=35.0,
            part_loc_d=25.0,
            off_grid_cluster_depth=50.0,
            off_grid_cluster_width=35.0,
            off_grid_arterial_clusters_depth=1,
            off_grid_secondary_clusters_depth=1,
            off_grid_local_clusters_depth=3,
            off_grid_local_clusters_width=4,
        )

        # Verify all custom values
        assert config.parcel_path == "/custom/parcels.geojson"
        assert config.roads_path == "/custom/roads.geojson"
        assert config.output_dir == "/custom/output"
        assert config.geopackage_path == "/custom/data.gpkg"
        assert config.road_arterial_width_m == 25.0
        assert config.road_secondary_width_m == 18.0
        assert config.on_grid_partition_depth_arterial_roads == 70.0
        assert config.on_grid_partition_depth_secondary_roads == 65.0
        assert config.off_grid_partitions_preferred_depth == 150.0
        assert config.off_grid_partitions_preferred_width == 160.0
        assert config.perpendicular_line_length == 1200.0
        assert config.optimize_grid_rotation is False
        assert config.grid_rotation_angle_step == 10.0
        assert config.use_ternary_search is True
        assert config.clip_to_boundary is False
        assert config.tolerance_area_ratio == 0.80
        assert config.tolerance_boundary_distance == 15.0
        assert config.sidewalk_width_m == 4.0
        assert config.road_locals_width_m == 12.0
        assert config.part_art_d == 50.0
        assert config.part_sec_d == 35.0
        assert config.part_loc_d == 25.0
        assert config.off_grid_cluster_depth == 50.0
        assert config.off_grid_cluster_width == 35.0
        assert config.off_grid_arterial_clusters_depth == 1
        assert config.off_grid_secondary_clusters_depth == 1
        assert config.off_grid_local_clusters_depth == 3
        assert config.off_grid_local_clusters_width == 4

    def test_street_config_road_width_relationships(self):
        """Test logical relationships between road width settings."""
        config = StreetConfig(parcel_path="parcels.geojson", roads_path="roads.geojson")

        # Arterial roads should be wider than secondary roads
        assert config.road_arterial_width_m > config.road_secondary_width_m

        # Secondary roads should be wider than local roads
        assert config.road_secondary_width_m > config.road_locals_width_m

    def test_street_config_cluster_depth_relationships(self):
        """Test logical relationships between cluster depths."""
        config = StreetConfig(parcel_path="parcels.geojson", roads_path="roads.geojson")

        # Arterial clusters should be deepest
        assert config.part_art_d > config.part_sec_d
        assert config.part_art_d > config.part_loc_d

        # Secondary clusters should be deeper than local
        assert config.part_sec_d > config.part_loc_d

    def test_street_config_positive_values(self):
        """Test that numeric configuration values are positive."""
        config = StreetConfig(parcel_path="parcels.geojson", roads_path="roads.geojson")

        # All dimensions should be positive
        assert config.road_arterial_width_m > 0
        assert config.road_secondary_width_m > 0
        assert config.road_locals_width_m > 0
        assert config.sidewalk_width_m > 0

        # All cluster dimensions should be positive
        assert config.part_art_d > 0
        assert config.part_sec_d > 0
        assert config.part_loc_d > 0
        assert config.off_grid_cluster_depth > 0
        assert config.off_grid_cluster_width > 0

        # All partition dimensions should be positive
        assert config.on_grid_partition_depth_arterial_roads > 0
        assert config.on_grid_partition_depth_secondary_roads > 0
        assert config.off_grid_partitions_preferred_depth > 0
        assert config.off_grid_partitions_preferred_width > 0

        # Other positive values
        assert config.perpendicular_line_length > 0
        assert config.grid_rotation_angle_step > 0
        assert config.tolerance_area_ratio > 0
        assert config.tolerance_boundary_distance > 0

    def test_street_config_cluster_counts_non_negative(self):
        """Test that cluster count values are non-negative integers."""
        config = StreetConfig(parcel_path="parcels.geojson", roads_path="roads.geojson")

        assert config.off_grid_arterial_clusters_depth >= 0
        assert config.off_grid_secondary_clusters_depth >= 0
        assert config.off_grid_local_clusters_depth >= 0
        assert config.off_grid_local_clusters_width >= 0

        # Should be integers
        assert isinstance(config.off_grid_arterial_clusters_depth, int)
        assert isinstance(config.off_grid_secondary_clusters_depth, int)
        assert isinstance(config.off_grid_local_clusters_depth, int)
        assert isinstance(config.off_grid_local_clusters_width, int)
