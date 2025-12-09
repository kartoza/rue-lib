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
        part_art_d: Depth of on-grid partitions along arterial roads (meters)
        part_sec_d: Depth of on-grid partitions along secondary roads (meters)
        part_loc_d: Depth of on-grid partitions along local roads (meters)
        part_og_d: Depth of off-grid partitions (meters)
        part_og_w: Width of off-grid strips (meters)
        plot_art_w: Typical plot width on arterial roads (meters)
        plot_sec_w: Typical plot width on secondary roads (meters)
        plot_loc_w: Typical plot width on local roads (meters)
        blk_art_num_og_d: Number of off-grid layers in depth behind arterial roads
        blk_sec_num_og_d: Number of off-grid layers in depth behind secondary roads
        blk_loc_num_og_d: Number of off-grid layers in depth behind local roads
        blk_art_num_og_w: Number of off-grid layers in width behind arterial roads
        blk_sec_num_og_w: Number of off-grid layers in width behind secondary roads
        blk_loc_num_og_w: Number of off-grid layers in width behind local roads
        open_percent: Target percent of site area for open space
        amen_percent: Target percent of site area for social amenities
        output_dir: Directory for output files
        road_arterial_width_m: Width of arterial roads (meters)
        road_secondary_width_m: Width of secondary roads (meters)
        road_local_width_m: Width of local roads (meters)
        off_grid_plot_threshold: Off-grid plot threshold (0.5 = 50% of target)
    """
    # Input paths
    roads_path: str  # Path to roads network
    input_path: str  # Path to streets output, or previous step output

    # Cluster / partition geometry
    # Depth of on-grid clusters along different road types (meters)
    part_art_d: float = 40.0  # Depth of on-grid clusters along arteries
    part_sec_d: float = 30.0  # Depth of on-grid clusters along secondaries
    part_loc_d: float = 20.0  # Depth of on-grid clusters along local roads

    # Off-grid cluster dimensions (meters)
    part_og_d: float = 140.0  # Depth of off-grid clusters
    part_og_w: float = 140.0  # Width of off-grid strips

    # Plot sizes (meters)
    # Typical plot widths on different road types
    plot_art_w: float = 20.0  # Typical plot width on an arterial
    plot_sec_w: float = 15.0  # Typical plot width on a secondary
    plot_loc_w: float = 12.0  # Typical plot width on a local road

    # Block → off-grid configuration
    # Number of off-grid layers in depth behind each road type
    blk_art_num_og_d: int = 2  # Number of off-grid layers (depth) behind arterials
    blk_sec_num_og_d: int = 2  # Number of off-grid layers (depth) behind secondaries
    blk_loc_num_og_d: int = 1  # Number of off-grid layers (depth) behind local roads

    # Number of off-grid layers in width behind each road type
    blk_art_num_og_w: int = 2  # Number of off-grid layers (width) behind arterials
    blk_sec_num_og_w: int = 2  # Number of off-grid layers (width) behind secondaries
    blk_loc_num_og_w: int = 1  # Number of off-grid layers (width) behind local roads

    # Public space & amenities (percentages)
    open_percent: float = 15.0  # Target percent of site area for open space
    amen_percent: float = 10.0  # Target percent of site area for social amenities

    # Output configuration
    output_dir: str = "outputs/clusters"

    # Road widths
    road_arterial_width_m: float = 20.0
    road_secondary_width_m: float = 15.0
    road_local_width_m: float = 12.0

    # Fixed data
    off_grid_plot_threshold: float = 0.5  # Off-grid plot threshold (0.5 = 50% of target, below it is park)
