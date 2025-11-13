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
