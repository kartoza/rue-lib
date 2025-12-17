# examples/step3_generate_clusters.py
import argparse

from rue_lib.cluster.runner import ClusterConfig, generate_clusters
from rue_lib.config import MainConfig


def main():
    """
    Step 3: Generate clusters/partitions from street blocks and roads.
    """
    parser = argparse.ArgumentParser(description="Generate clusters/partitions from street blocks")
    parser.add_argument(
        "--input",
        default="outputs/step2_streets/all_grids_merged.geojson",
        help="Path to street blocks geojson file",
    )
    parser.add_argument(
        "--output-dir", default="outputs/step3_clusters", help="Output directory for clusters"
    )
    parser.add_argument(
        "--geopackage", default="outputs/output.gpkg", help="Path to output geopackage file"
    )
    args = parser.parse_args()

    config = ClusterConfig(
        roads_path="data/roads.geojson",
        input_path=args.input,
        output_dir=args.output_dir,
        geopackage_path=args.geopackage,
        # Neighborhood / public roads
        road_arterial_width_m=MainConfig.road_arterial_width_m,
        road_secondary_width_m=MainConfig.road_secondary_width_m,
        road_local_width_m=MainConfig.road_local_width_m,
        # Neighbourhood / on-grid partitions
        on_grid_partition_depth_arterial_roads=MainConfig.on_grid_partition_depth_arterial_roads,
        on_grid_partition_depth_secondary_roads=MainConfig.on_grid_partition_depth_secondary_roads,
        # Neighbourhood / off-grid partitions
        off_grid_cluster_depth=MainConfig.off_grid_cluster_depth,
        off_grid_cluster_width=MainConfig.off_grid_cluster_width,
        # Neighborhood / public spaces
        sidewalk_width_m=MainConfig.sidewalk_width_m,
    )

    print("=" * 60)
    print("STEP 3: Generating Clusters/Partitions")
    print("=" * 60)
    print(f"Input: {config.input_path}")
    print(f"Output: {config.output_dir}")
    print()

    output_path = generate_clusters(config)

    print()
    print("✔ Step 3 completed successfully!")
    print(f"  Output: {output_path}")


if __name__ == "__main__":
    main()
