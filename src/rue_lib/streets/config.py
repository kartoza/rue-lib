# src/rue_lib/streets/config.py
from dataclasses import dataclass


@dataclass
class StreetConfig:
    """Configuration for street generation."""

    parcel_path: str  # Output generated from step 1
    roads_path: str
    road_arterial_width_m: float = 20.0
    road_secondary_width_m: float = 15.0
    on_grid_partition_depth_arterial_roads: float = 40.0
    on_grid_partition_depth_secondary_roads: float = 30.0
    off_grid_partitions_preferred_depth: float = 140.0
    off_grid_partitions_preferred_width: float = 140.0
    arterial_setback_depth: float = 60.0  # Depth of arterial road setback zone
    secondary_setback_depth: float = 60.0  # Depth of secondary road setback zone
    perpendicular_line_length: float = 1000.0  # Length of perpendicular lines
    output_dir: str = "outputs/streets"
    road_local_width_m: float = 12.0
    optimize_grid_rotation: bool = True
    grid_rotation_angle_step: float = 5.0
    use_ternary_search: bool = False
    clip_to_boundary: bool = True
    tolerance_area_ratio: float = (
        0.70  # Allow cells with area ratio >= this (0.0-1.0). 0.0=perfect only, 0.95=allow 95%+
    )
    tolerance_boundary_distance: float = (
        10.0  # Allow cells within this distance of boundary (meters). 0.0=no touching
    )
    sidewalk_width_m: float = 3.0
    road_locals_width_m: float = 10.0
