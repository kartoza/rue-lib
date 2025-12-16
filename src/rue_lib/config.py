class MainConfig:
    """
    Configuration for street network generation and urban block clustering.

    This class defines parameters for:
    - Road widths for different road types (arterial, secondary, local)
    - On-grid partition depths along different road types
    - Off-grid cluster dimensions (depth and width)
    - Off-grid cluster counts for different road hierarchies
    - Public space dimensions (sidewalks)

    Attributes:
        road_arterial_width_m: Width of arterial roads in meters
        road_secondary_width_m: Width of secondary roads in meters
        road_local_width_m: Width of local roads in meters
        on_grid_partition_depth_arterial_roads: Depth of parcels along arterial roads in meters
        on_grid_partition_depth_secondary_roads: Depth of parcels along secondary roads in meters
        on_grid_partition_depth_local_roads: Depth of parcels along local roads in meters
        off_grid_cluster_depth: Depth of off-grid clusters in meters
        off_grid_cluster_width: Width of off-grid clusters in meters
        off_grid_arterial_clusters_depth: Number of cluster rows along arterial roads
        off_grid_arterial_clusters_width: Number of cluster columns along arterial roads
        off_grid_secondary_clusters_depth: Number of cluster rows along secondary roads
        off_grid_secondary_clusters_width: Number of cluster columns along secondary roads
        off_grid_local_clusters_depth: Number of cluster rows along local roads
        off_grid_local_clusters_width: Number of cluster columns along local roads
        sidewalk_width_m: Width of sidewalks in meters

    Notes:
        - All dimensional values are in meters
        - Road widths represent the full width of the road
        - Cluster counts are integers representing the number of clusters
        - These values can be adjusted based on local planning standards and regulations
    """

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

    # Neighbourhood / urban block structure
    off_grid_arterial_clusters_depth: int = 0
    off_grid_arterial_clusters_width: int = 3
    off_grid_secondary_clusters_depth: int = 0
    off_grid_secondary_clusters_width: int = 3
    off_grid_local_clusters_depth: int = 2
    off_grid_local_clusters_width: int = 3

    # Neighborhood / public spaces
    open_percent: float = 6.0
    amen_percent: float = 5.0
    sidewalk_width_m: float = 3.0
