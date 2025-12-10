# src/rue_lib/cluster/config.py
from dataclasses import dataclass


@dataclass
class ClusterConfig:
    """
    Configuration for cluster/partition generation.

    This dataclass holds all configuration parameters for generating clusters/partitions
    from street blocks, including partition dimensions, plot sizes, road widths, and
    output settings.

    Attributes:
        roads_path: Path to roads network (GeoJSON or GeoPackage)
        input_path: Path to streets output or previous step output
        output_dir: Directory for output files
        road_arterial_width_m: Width of arterial roads (meters)
        road_secondary_width_m: Width of secondary roads (meters)
        road_local_width_m: Width of local roads (meters)
        off_grid_plot_threshold: Off-grid plot threshold (0.5 = 50% of target)
    """
    # Input paths
    roads_path: str  # Path to roads network
    input_path: str  # Path to streets output, or previous step output
    output_dir: str  # Path to output directory

    # Fixed data
    # Off-grid plot threshold (0.5 = 50% of target, below it is park)
    off_grid_plot_threshold: float = 0.5

    # Neighborhood / public roads
    road_arterial_width_m: float = 20.0
    road_secondary_width_m: float = 15.0
    road_local_width_m: float = 10.0

    # Neighbourhood / on-grid partitions
    on_grid_partition_depth_arterial_roads: float = 40.0
    on_grid_partition_depth_secondary_roads: float = 30.0
    on_grid_partition_depth_local_roads: float = 20.0

    # Neighbourhood / off-grid partitions
    off_grid_cluster_depth: float = 45.0
    off_grid_cluster_width: float = 30.0

    # Neighborhood / public spaces
    sidewalk_width_m: float = 3.0
