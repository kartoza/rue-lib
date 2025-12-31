# examples/step2_generate_streets_with_local_roads.py
"""
Example: Run the streets pipeline and replace local streets with a provided GeoJSON.

The local streets centerlines are taken from data/local_streets.geojson and buffered
according to StreetConfig. Outputs are written to outputs/step2_streets_local.
"""

from rue_lib.config import MainConfig
from rue_lib.streets import StreetConfig, generate_streets_with_local_roads


def main():
    cfg = StreetConfig(
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
        dead_end_buffer_distance=MainConfig.dead_end_buffer_distance,
    )

    local_roads_geojson = "data/local_streets.geojson"

    print("=" * 60)
    print("STEP 2 (with provided local streets): Generating Streets")
    print("=" * 60)
    print(f"Parcels: {cfg.parcel_path}")
    print(f"Roads: {cfg.roads_path}")
    print(f"Local streets (provided): {local_roads_geojson}")
    print(f"Output: {cfg.output_dir}")
    print()

    output = generate_streets_with_local_roads(cfg, local_roads_geojson)

    print()
    print("✔ Streets generation completed with provided local streets")
    print(f"  GeoPackage: {output}")


if __name__ == "__main__":
    main()
