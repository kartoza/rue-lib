# examples/step2_generate_streets.py

from rue_lib.config import MainConfig
from rue_lib.streets.runner import StreetConfig, generate_streets


def main():
    """
    Step 2: Generate street blocks from parcels and roads.
    """

    config = StreetConfig(
        parcel_path="outputs/step1_parcels/outputs.geojson",
        roads_path="data/roads.geojson",
        output_dir="outputs/step2_streets",
        # Neighborhood / public roads
        road_arterial_width_m=MainConfig.road_arterial_width_m,
        road_secondary_width_m=MainConfig.road_secondary_width_m,
        road_locals_width_m=MainConfig.road_local_width_m,
        part_art_d=MainConfig.on_grid_partition_depth_arterial_roads,
        part_sec_d=MainConfig.on_grid_partition_depth_secondary_roads,
        part_loc_d=MainConfig.on_grid_partition_depth_local_roads,
        # Neighbourhood / on-grid partitions
        on_grid_partition_depth_arterial_roads=MainConfig.on_grid_partition_depth_arterial_roads,
        on_grid_partition_depth_secondary_roads=MainConfig.on_grid_partition_depth_secondary_roads,
        # Neighbourhood / off-grid partitions
        off_grid_cluster_depth=MainConfig.off_grid_cluster_depth,
        off_grid_cluster_width=MainConfig.off_grid_cluster_width,
        # Neighbourhood / urban block structure
        off_grid_arterial_clusters_depth=MainConfig.off_grid_arterial_clusters_depth,
        off_grid_secondary_clusters_depth=MainConfig.off_grid_secondary_clusters_depth,
        off_grid_local_clusters_depth=MainConfig.off_grid_local_clusters_depth,
        off_grid_local_clusters_width=MainConfig.off_grid_local_clusters_width,
        # Neighborhood / public spaces
        sidewalk_width_m=MainConfig.sidewalk_width_m,
    )

    print("=" * 60)
    print("STEP 2: Generating Street Blocks")
    print("=" * 60)
    print(f"Parcels: {config.parcel_path}")
    print(f"Roads: {config.roads_path}")
    print(f"Output: {config.output_dir}")
    print()

    output_path = generate_streets(config)

    print()
    print("✔ Step 2 completed successfully!")
    print(f"  Output: {output_path}")


if __name__ == "__main__":
    main()
