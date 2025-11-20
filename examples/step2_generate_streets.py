# examples/step2_generate_streets.py

from rue_lib.streets.runner import StreetConfig, generate_streets


def main():
    """
    Step 2: Generate street blocks from parcels and roads.
    """

    config = StreetConfig(
        parcel_path="outputs/step1_parcels/subsites.geojson",
        roads_path="data/roads.geojson",
        on_grid_partition_depth_arterial_roads=40.0,
        on_grid_partition_depth_secondary_roads=30.0,
        off_grid_partitions_preferred_depth=100,
        off_grid_partitions_preferred_width=100,
        output_dir="outputs/step2_streets",
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
