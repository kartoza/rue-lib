# examples/step3_generate_clusters.py

from rue_lib.cluster.runner import ClusterConfig, generate_clusters


def main():
    """
    Step 2: Generate street blocks from parcels and roads.
    """

    config = ClusterConfig(
        roads_path="data/roads.geojson",
        input_path="outputs/step2_streets/all_grids_merged.geojson",
        output_dir="outputs/step3_clusters",
        part_art_d=40,
        part_sec_d=30,
        part_loc_d=20,

        part_og_d=45.0,
        part_og_w=30.0,

        # Road widths
        road_arterial_width_m=40.0,
        road_secondary_width_m=30.0,
        road_local_width_m=10.0
    )

    print("=" * 60)
    print("STEP 3: Generating Clusters/Partitions")
    print("=" * 60)
    print(f"Input: {config.input_path}")
    print(f"Output: {config.output_dir}")
    print()

    output_path = generate_clusters(config)

    print()
    print("✔ Step 2 completed successfully!")
    print(f"  Output: {output_path}")


if __name__ == "__main__":
    main()
