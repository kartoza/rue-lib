# src/rue_lib/streets/runner.py
from pathlib import Path

from osgeo import gdal, ogr

from rue_lib.core.geometry import buffer_layer, get_utm_zone_from_layer, reproject_layer
from rue_lib.streets.blocks import generate_on_grid_blocks
from rue_lib.streets.grids import (
    create_cold_boundaries,
    extract_grid_lines_in_buffer,
    generate_local_streets,
    grids_from_site,
)
from rue_lib.streets.runner_utils import (
    apply_inner_buffer_to_cells,
    create_guide_points_from_site_boundary,
    create_perpendicular_lines_from_guide_points,
    create_perpendicular_lines_inside_buffer_from_points,
)

from .cell import fix_grid_cells_with_perpendicular_lines, remove_dead_end_cells
from .config import StreetConfig
from .operations import (
    create_on_grid_cells_from_perpendiculars,
    erase_layer,
    export_layer_to_geojson,
    extract_by_expression,
    extract_site_boundary_lines,
    merge_grid_layers_with_type,
    merge_layers_without_overlaps,
)

gdal.UseExceptions()


def generate_streets(cfg: StreetConfig) -> Path:
    """
    Generate street blocks from roads and parcels (version 3)

    Returns:
        Path to output blocks file
    """
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_gpkg = output_dir / "outputs.gpkg"
    output_path = str(output_gpkg)

    print("Step 1: Determining UTM zone...")
    site_ds = ogr.Open(cfg.parcel_path)
    site_layer = site_ds.GetLayer()
    utm_epsg = get_utm_zone_from_layer(site_layer)
    site_ds = None
    print(f"Using UTM EPSG: {utm_epsg}")

    print("Step 2: Reprojecting layers to UTM...")
    site_layer_name = reproject_layer(cfg.parcel_path, output_path, utm_epsg)
    roads_layer_name = reproject_layer(cfg.roads_path, output_path, utm_epsg)

    print("Step 4: Extracting arterial roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_art'", output_path, "04_arterial_roads"
    )

    print("Step 5: Extracting secondary roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_sec'", output_path, "05_secondary_roads"
    )

    print("Step 6: Creating arterial road setback zone...")
    buffer_layer(
        output_path,
        "04_arterial_roads",
        cfg.on_grid_partition_depth_arterial_roads,
        output_path,
        "06_arterial_setback",
        dissolve=True,
    )

    print("Step 7: Creating secondary road setback zone...")
    buffer_layer(
        output_path,
        "05_secondary_roads",
        cfg.on_grid_partition_depth_secondary_roads,
        output_path,
        "07_secondary_setback",
        dissolve=True,
    )

    print("Step 8: Removing arterial setback from site...")
    erase_layer(
        output_path,
        site_layer_name,
        output_path,
        "06_arterial_setback",
        output_path,
        "08_site_minus_arterial_setback",
    )

    print("Step 9: Removing secondary setback from site...")
    erase_layer(
        output_path,
        "08_site_minus_arterial_setback",
        output_path,
        "07_secondary_setback",
        output_path,
        "09_site_minus_all_setbacks",
    )

    print("Step 10: Generating on-grid blocks...")
    generate_on_grid_blocks(output_gpkg, site_layer_name, cfg)

    print("Step 13: Extracting site boundary lines")
    extract_site_boundary_lines(
        output_gpkg,
        "09_site_minus_all_setbacks",
        "10a_arterial_setback_clipped",
        "10b_secondary_setback_clipped",
        output_gpkg,
        "13_site_boundary_lines",
    )

    print("Step 14: Generating off-grid blocks...")
    grids_from_site(
        output_gpkg,
        "09_site_minus_all_setbacks",
        "13_site_boundary_lines",
        cfg.off_grid_partitions_preferred_width,
        cfg.off_grid_partitions_preferred_depth,
        grid_layer_name="14_off_grid_cells",
        point_layer_name="14_off_grid_points",
    )

    print("Step 14a: Creating inner buffer zone from site boundary...")
    _buffered_layer = apply_inner_buffer_to_cells(
        output_gpkg,
        "09_site_minus_all_setbacks",
        cfg.part_loc_d + cfg.sidewalk_width_m + cfg.road_locals_width_m / 2.0,
    )

    print("Step 14b: Extracting grid lines inside buffer zone...")
    lines_layer = extract_grid_lines_in_buffer(
        output_gpkg,
        "14_off_grid_cells",
        "14a_site_boundary_inner_buffer",
    )
    fixed_cells_layer = None
    if lines_layer is not None:
        perp_inside_layer = create_perpendicular_lines_inside_buffer_from_points(
            output_gpkg,
            lines_layer,
            "14a_site_boundary_inner_buffer",
            line_length=cfg.part_loc_d * 2,
        )

        if perp_inside_layer is not None:
            print("Step 14d: Fixing grid cells with perpendicular lines...")
            fixed_cells_layer = fix_grid_cells_with_perpendicular_lines(
                output_gpkg,
                "14_off_grid_cells",
                perp_inside_layer,
                "14a_site_boundary_inner_buffer",
                target_area=cfg.off_grid_partitions_preferred_width
                * cfg.off_grid_partitions_preferred_depth,
            )

    cleaned_cells_layer = ""
    if fixed_cells_layer is not None:
        print("Step 14e: Removing dead end cells")
        cleaned_cells_layer = remove_dead_end_cells(
            output_gpkg,
            fixed_cells_layer,
            "13_site_boundary_lines",
            "09_site_minus_all_setbacks",
        )

    create_cold_boundaries(
        output_gpkg,
        site_layer_name,
        cleaned_cells_layer,
        arterial_setback_layer="10a_arterial_setback_clipped",
        secondary_setback_layer="10b_secondary_setback_clipped",
        output_layer_name="14_cold_boundaries",
    )

    print("Step 15: Generating on-grid cells")

    print("Step 15a: Creating guide points from site boundary...")
    _guide_points_layer = create_guide_points_from_site_boundary(
        output_gpkg,
        "13_site_boundary_lines",
        perp_inside_layer,
        output_layer_name="13_site_boundary_points",
    )

    print("Step 15b: Creating perpendicular lines from guide points...")
    _guide_perp_layer = create_perpendicular_lines_from_guide_points(
        output_gpkg,
        "13_site_boundary_lines",
        "09_site_minus_all_setbacks",
        "13_site_boundary_points",
        line_length=max(
            cfg.on_grid_partition_depth_secondary_roads, cfg.on_grid_partition_depth_arterial_roads
        )
        * 1.05,
        output_layer_name="13_site_boundary_perp_from_points",
    )

    print("Merge arterial and secondary setbacks with overlaps resolved...")
    merge_layers_without_overlaps(
        output_path,
        ["10a_arterial_setback_clipped", "10b_secondary_setback_clipped"],
        output_path,
        "10_setback_clipped_merged",
    )

    print("Step 16: Creating on-grid cells from merged setbacks and perpendiculars...")
    _on_grid_cells_layer = create_on_grid_cells_from_perpendiculars(
        output_gpkg,
        "10_setback_clipped_merged",
        "13_site_boundary_perp_from_points",
        output_gpkg,
        "16_on_grid_cells",
    )

    on_grid_outer_layer, on_grid_inner_layer = generate_local_streets(
        output_gpkg, cfg, "16_on_grid_cells"
    )
    off_grid_outer_layer, off_grid_inner_layer = generate_local_streets(
        output_gpkg, cfg, cleaned_cells_layer
    )

    print("Step 17: Merging all grid layers with grid_type information...")
    merge_grid_layers_with_type(
        str(output_gpkg),
        str(output_gpkg),
        "17_all_grids_merged",
        [
            (cleaned_cells_layer, "off_grid"),
            (off_grid_inner_layer, "off_grid"),
            ("16_on_grid_cells", "on_grid"),
            (on_grid_inner_layer, "on_grid"),
            ("14_cold_boundaries", "cold_boundaries"),
        ],
    )

    print("Step 18: Exporting merged grids to GeoJSON...")
    output_geojson = output_dir / "all_grids_merged.geojson"
    export_layer_to_geojson(
        str(output_gpkg),
        "17_all_grids_merged",
        str(output_geojson),
    )

    print(f"\nProcessing complete! Output saved to: {output_gpkg}")
    print("\nFinal layers:")
    print(f"  - {site_layer}: Site polygon with roads subtracted")
    print("  - 17_all_grids_merged: Merged grid cells with grid_type classification")
    print("\nGeoJSON export:")
    print(f"  - {output_geojson}: Merged grids with grid_type classification")

    return output_gpkg
