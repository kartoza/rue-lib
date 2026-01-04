# src/rue_lib/streets/config.py
from dataclasses import dataclass


@dataclass
class StreetConfig:
    """Configuration for street generation."""

    parcel_path: str  # Output generated from step 1
    roads_path: str
    output_dir: str = "outputs"
    geopackage_path: str = "outputs/output.gpkg"
    road_arterial_width_m: float = 20.0
    road_secondary_width_m: float = 15.0
    on_grid_partition_depth_arterial_roads: float = 60.0
    on_grid_partition_depth_secondary_roads: float = 60.0
    off_grid_partitions_preferred_depth: float = 140.0
    off_grid_partitions_preferred_width: float = 140.0
    perpendicular_line_length: float = 1000.0  # Length of perpendicular lines
    optimize_grid_rotation: bool = True
    grid_rotation_angle_step: float = 5.0
    use_ternary_search: bool = False
    clip_to_boundary: bool = True
    tolerance_area_ratio: float = 0.70
    tolerance_boundary_distance: float = 10.0
    sidewalk_width_m: float = 3.0
    road_locals_width_m: float = 10.0

    # Cluster / partition geometry
    # Depth of on-grid clusters along different road types (meters)
    part_art_d: float = 40.0  # Depth of on-grid clusters along arteries
    part_sec_d: float = 30.0  # Depth of on-grid clusters along secondaries
    part_loc_d: float = 20.0  # Depth of on-grid clusters along local roads

    off_grid_cluster_depth: float = 45.0  # Depth of off-grid clusters (meters)
    off_grid_cluster_width: float = 30.0  # Width of off-grid clusters (meters)

    # Number of clusters along depth for arterial road off-grid partitions
    off_grid_arterial_clusters_depth: int = 0
    # Number of clusters along depth for secondary road off-grid partitions
    off_grid_secondary_clusters_depth: int = 0
    # Number of clusters along depth for local road off-grid partitions
    off_grid_local_clusters_depth: int = 2
    # Number of clusters along width for local road off-grid partitions
    off_grid_local_clusters_width: int = 3

    # Streets / dead end boundary
    dead_end_buffer_distance: float = 15.0
