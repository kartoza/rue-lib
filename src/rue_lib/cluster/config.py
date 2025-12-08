# src/rue_lib/cluster/config.py
from dataclasses import dataclass


@dataclass
class ClusterConfig:
    """Configuration for cluster/partition generation."""

    # Input paths
    site_path: str  # Path to site boundary
    roads_path: str  # Path to roads network
    blocks_path: str  # Path to blocks (output from streets module)

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
