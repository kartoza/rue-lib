# src/rue_lib/streets/runner.py
from pathlib import Path

from osgeo import gdal, ogr

from rue_lib.core.geometry import buffer_layer, get_utm_zone_from_layer, reproject_layer

from .blocks import (
    filter_offgrid_blocks,
    merge_lines,
    polygonize_and_classify_blocks,
)
from .config import StreetConfig
from .lines import (
    create_division_points,
    create_perpendicular_lines,
    extract_arterial_edge_lines,
)
from .operations import (
    calculate_required_rings,
    cleanup_intermediate_layers,
    clip_layer,
    extract_by_expression,
    multiring_buffer,
)

# Enable GDAL exceptions
gdal.UseExceptions()


def generate_streets(cfg: StreetConfig) -> Path:
    """
    Generate street blocks from roads and parcels

    Returns:
        Path to output blocks file
    """
    # Create output directory
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_gpkg = output_dir / "outputs.gpkg"
    output_path = str(output_gpkg)

    print("Step 1: Determining UTM zone...")
    # Get UTM zone from site layer
    site_ds = ogr.Open(cfg.parcel_path)
    site_layer = site_ds.GetLayer()
    utm_epsg = get_utm_zone_from_layer(site_layer)
    site_ds = None
    print(f"Using UTM EPSG: {utm_epsg}")

    print("Step 2: Reprojecting layers to UTM...")
    site_layer_name = reproject_layer(cfg.parcel_path, output_path, utm_epsg)
    roads_layer_name = reproject_layer(cfg.roads_path, output_path, utm_epsg)

    print("Step 3: Extracting arterial roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_art'", output_path, "arterial_roads"
    )

    print("Step 4: Extracting secondary roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_sec'", output_path, "secondary_roads"
    )

    print("Step 5: Creating arterial road setback zone...")
    # Use the larger setback depth for the arterial setback
    buffer_layer(
        output_path,
        "arterial_roads",
        cfg.on_grid_partition_depth_arterial_roads,
        output_path,
        "arterial_buffered_large",
        dissolve=True,
    )

    print("Step 6: Clipping arterial setback by site...")
    clip_layer(
        output_path,
        "arterial_buffered_large",
        output_path,
        site_layer_name,
        output_path,
        "arterial_road_setback",
    )

    print("Step 7: Creating street block rings...")

    num_rings = calculate_required_rings(
        output_path, site_layer_name, cfg.off_grid_partitions_preferred_depth
    )
    multiring_buffer(
        output_path,
        "arterial_buffered_large",
        num_rings,
        cfg.off_grid_partitions_preferred_depth,
        output_path,
        "arterial_offset_buffer_rings",
    )

    print("Step 8: Clipping block rings by site...")
    clip_layer(
        output_path,
        "arterial_offset_buffer_rings",
        output_path,
        site_layer_name,
        output_path,
        "street_blocks",
    )

    print("Step 9: Creating secondary road setback zone...")
    # Use the larger setback depth for secondary roads
    buffer_layer(
        output_path,
        "secondary_roads",
        cfg.secondary_setback_depth,
        output_path,
        "secondary_roads_buffered_large",
        dissolve=True,
    )

    print("Step 10: Clipping secondary setback by site...")
    clip_layer(
        output_path,
        "secondary_roads_buffered_large",
        output_path,
        site_layer_name,
        output_path,
        "secondary_road_setback",
    )

    print("Step 11: Extracting arterial edge lines...")
    extract_arterial_edge_lines(
        output_path,
        "arterial_roads",
        "arterial_road_setback",
        "secondary_road_setback",
        "arterial_edge_lines",
        clip_buffer=0.1,
    )

    print("Step 12: Creating division points along arterial edges...")
    create_division_points(
        output_path,
        "arterial_edge_lines",
        "division_points",
        cfg.off_grid_partitions_preferred_width,
    )

    print("Step 13: Creating perpendicular lines from division points...")
    create_perpendicular_lines(
        output_path,
        "division_points",
        "arterial_edge_lines",
        site_layer_name,
        "perpendicular_lines",
        cfg.perpendicular_line_length,
    )

    print("Step 14: Merging all lines...")
    merge_lines(
        output_path,
        "perpendicular_lines",
        "arterial_edge_lines",
        "street_blocks",
        "secondary_road_setback",
        "arterial_road_setback",
        "merged_lines",
    )

    print("Step 15: Classifying blocks by setback adjacency...")
    polygonize_and_classify_blocks(
        output_path,
        "merged_lines",
        "arterial_road_setback",
        "secondary_road_setback",
        "classified_blocks",
    )

    print("Step 16: Filtering off-grid blocks by shape and size...")
    filter_offgrid_blocks(
        output_path,
        "classified_blocks",
        "filtered_blocks",
        cfg.off_grid_partitions_preferred_width,
        cfg.off_grid_partitions_preferred_depth,
        cfg.on_grid_partition_depth_arterial_roads,
        cfg.on_grid_partition_depth_secondary_roads,
        arterial_preferred_width=cfg.off_grid_partitions_preferred_width,
        secondary_preferred_width=cfg.off_grid_partitions_preferred_width,
        area_threshold=0.6,
        squareness_threshold=0.7,
    )

    # Clean up intermediate layers
    print("Cleaning up intermediate layers...")
    final_layers = [
        "arterial_roads",
        "arterial_road_setback",
        "street_blocks",
        "secondary_road_setback",
        "arterial_edge_lines",
        "division_points",
        "perpendicular_lines",
        "merged_lines",
        "classified_blocks",
        "filtered_blocks",
    ]
    cleanup_intermediate_layers(output_path, final_layers)

    print(f"\nProcessing complete! Output saved to: {output_gpkg}")
    print("\nFinal layers:")
    print(
        f"  - arterial_road_setback: {cfg.arterial_setback_depth}m "
        f"buffer zone around arterial roads"
    )
    print(
        f"  - street_blocks: Concentric block rings ({cfg.off_grid_partitions_preferred_depth}m "
        f"spacing)"
    )
    print(
        f"  - secondary_road_setback: {cfg.secondary_setback_depth}m buffer "
        f"zone around secondary roads"
    )

    return output_gpkg
