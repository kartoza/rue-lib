# src/rue_lib/streets/runner.py
from pathlib import Path

import geopandas as gpd
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
    create_dead_end_boundary_lines,
    create_guide_points_from_site_boundary,
    create_perpendicular_lines_from_guide_points,
    create_perpendicular_lines_inside_buffer_from_points,
    polygons_to_lines_layer,
    subtract_layer,
)
from rue_lib.utils.io import prepare_geopackage

from ..site.roads import buffer_roads
from .cell import (
    fix_grid_cells_with_perpendicular_lines,
    merge_small_cells_with_neighbors,
    remove_dead_end_cells,
)
from .config import StreetConfig
from .financial import FinancialStreet
from .operations import (
    classify_on_grid_cells_by_setback,
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

    output_gpkg = Path(cfg.geopackage_path)

    # Prepare geopackage with template (creates from template if needed)
    prepare_geopackage(output_gpkg)

    output_path = str(output_gpkg)

    print("Step 1: Determining UTM zone...")
    site_ds = ogr.Open(cfg.parcel_path)
    site_layer = site_ds.GetLayer()
    utm_epsg = get_utm_zone_from_layer(site_layer)
    site_ds = None
    print(f"Using UTM EPSG: {utm_epsg}")

    print("Step 2: Reprojecting layers to UTM...")
    site_layer_name = reproject_layer(
        cfg.parcel_path, output_path, utm_epsg, layer_name="00_site", is_append_epsg=False
    )
    roads_layer_name = reproject_layer(
        cfg.roads_path, output_path, utm_epsg, layer_name="00_roads", is_append_epsg=False
    )
    input_roads_buffer_layer_name = "00_roads_buffer"
    roads_m = gpd.read_file(output_path, layer=roads_layer_name)
    roads_buf_m = buffer_roads(
        roads_m, cfg.road_arterial_width_m, cfg.road_secondary_width_m, road_type_key="road_type"
    )
    roads_buf_m.to_file(output_gpkg, layer=input_roads_buffer_layer_name, driver="GPKG")

    print("Step 4: Extracting arterial roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_art'", output_path, "04_arterial_roads"
    )

    print("Step 5: Extracting secondary roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_sec'", output_path, "05_secondary_roads"
    )

    print("Step 5a: Extracting local roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_loc'", output_path, "05_local_roads"
    )

    preferred_depth_on_grid_arterial = (
        cfg.part_art_d
        + (cfg.off_grid_arterial_clusters_depth * cfg.off_grid_cluster_depth)
        + cfg.road_locals_width_m
        + cfg.part_loc_d
    )

    preferred_depth_on_grid_secondary = (
        cfg.part_sec_d
        + (cfg.off_grid_secondary_clusters_depth * cfg.off_grid_cluster_depth)
        + cfg.road_locals_width_m
        + cfg.part_loc_d
    )

    preferred_depth_off_cluster_grid = (
        cfg.part_loc_d * 2
        + (cfg.off_grid_cluster_depth * cfg.off_grid_local_clusters_depth)
        + cfg.road_locals_width_m
    )

    preferred_width_off_cluster_grid = (
        cfg.part_loc_d * 2
        + (cfg.off_grid_cluster_width * cfg.off_grid_local_clusters_width)
        + cfg.road_locals_width_m
    )

    print("Step 6: Creating arterial road setback zone...")
    buffer_layer(
        output_path,
        "04_arterial_roads",
        preferred_depth_on_grid_arterial,
        output_path,
        "06_arterial_setback",
        dissolve=True,
    )

    print("Step 7: Creating secondary road setback zone...")
    buffer_layer(
        output_path,
        "05_secondary_roads",
        preferred_depth_on_grid_secondary,
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

    print("Step 13a: Deriving dead-end boundary lines (site edges minus boundary lines)...")
    dead_end_lines_layer = create_dead_end_boundary_lines(
        output_gpkg,
        "09_site_minus_all_setbacks",
        "13_site_boundary_lines",
        output_layer_name="13a_dead_end_lines",
        buffer_distance=min(
            cfg.on_grid_partition_depth_arterial_roads, cfg.on_grid_partition_depth_secondary_roads
        )
        * 0.4,
    )

    print("Step 14: Generating off-grid blocks...")
    grids_from_site(
        output_gpkg,
        "09_site_minus_all_setbacks",
        "13_site_boundary_lines",
        preferred_width_off_cluster_grid,
        preferred_depth_off_cluster_grid,
        grid_layer_name="14_off_grid_cells",
        point_layer_name="14_off_grid_points",
        dead_end_lines_layer=dead_end_lines_layer,
    )

    print("Step 14a: Creating inner buffer zone from site boundary...")
    _buffered_layer = apply_inner_buffer_to_cells(
        output_gpkg,
        "09_site_minus_all_setbacks",
        cfg.part_loc_d + cfg.road_locals_width_m / 2.0,
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
                target_area=preferred_depth_off_cluster_grid * preferred_width_off_cluster_grid,
            )

    print("Step 14e: Removing dead end cells")
    cleaned_cells_layer = remove_dead_end_cells(
        output_gpkg, fixed_cells_layer, "13a_dead_end_lines_buffered"
    )

    print("Step 14f: Merging small cells with neighbors")
    cleaned_cells_layer = merge_small_cells_with_neighbors(
        output_gpkg,
        cleaned_cells_layer,
        target_area=preferred_depth_off_cluster_grid * preferred_width_off_cluster_grid,
        area_threshold_ratio=0.5,
    )

    print("Step 15: Generating on-grid cells")

    print("Step 15a: Creating guide points from site boundary...")
    _guide_points_layer = create_guide_points_from_site_boundary(
        output_gpkg,
        "13_site_boundary_lines",
        perp_inside_layer,
        output_layer_name="13_site_boundary_points",
        lines_without_points_layer="13c_lines_without_points",
        min_line_length_threshold=preferred_width_off_cluster_grid,
    )

    print("Step 15b: Creating perpendicular lines from guide points...")
    _guide_perp_layer = create_perpendicular_lines_from_guide_points(
        output_gpkg,
        "13_site_boundary_lines",
        "09_site_minus_all_setbacks",
        "13_site_boundary_points",
        line_length=max(preferred_depth_on_grid_secondary, preferred_depth_on_grid_arterial) * 1.05,
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

    print("Step 16a: Remove dead-end on-grid cells")
    cleaned_on_grid_cells_layer = remove_dead_end_cells(
        output_gpkg, _on_grid_cells_layer, "13a_dead_end_lines_buffered"
    )

    print("Step 16b: Classifying on-grid cells by road type (arterial vs secondary)...")
    buffer_layer(
        output_path,
        "04_arterial_roads",
        preferred_depth_on_grid_arterial / 2,
        output_path,
        "06_arterial_buffer_for_classification",
        dissolve=True,
    )
    arterial_cells_layer, secondary_cells_layer = classify_on_grid_cells_by_setback(
        output_gpkg,
        cleaned_on_grid_cells_layer,
        "06_arterial_buffer_for_classification",
        output_gpkg,
        "16a_on_grid_arterial_cells",
        "16b_on_grid_secondary_cells",
    )

    print("Step 17: Creating cold boundaries...")
    create_cold_boundaries(
        output_gpkg,
        site_layer_name,
        cleaned_cells_layer,
        cleaned_on_grid_cells_layer,
        output_layer_name="14_cold_boundaries",
    )

    on_grid_arterial_outer_layer, on_grid_arterial_inner_layer = generate_local_streets(
        output_gpkg, cfg, arterial_cells_layer
    )
    on_grid_secondary_outer_layer, on_grid_secondary_inner_layer = generate_local_streets(
        output_gpkg, cfg, secondary_cells_layer
    )
    off_grid_outer_layer, off_grid_inner_layer = generate_local_streets(
        output_gpkg, cfg, cleaned_cells_layer
    )

    print("Exporting local roads from off-grid cells as linework...")
    local_roads_layer_name = "18_local_roads"
    polygons_to_lines_layer(
        output_gpkg,
        [
            cleaned_cells_layer,
            cleaned_on_grid_cells_layer,
        ],
        output_gpkg,
        local_roads_layer_name,
    )

    # Remove cold boundaries from local roads
    subtract_layer(
        output_gpkg,
        "14_cold_boundaries",
        cleaned_on_grid_cells_layer,
        output_gpkg,
        "14a_local_roads_minus_cold_boundaries",
        cfg.road_locals_width_m / 2.0,
    )

    subtract_layer(
        output_gpkg,
        "14a_local_roads_minus_cold_boundaries",
        cleaned_cells_layer,
        output_gpkg,
        "14a_local_roads_minus_cold_boundaries",
        cfg.road_locals_width_m / 2.0,
    )

    subtract_layer(
        output_gpkg,
        "14a_local_roads_minus_cold_boundaries",
        "05_secondary_roads",
        output_gpkg,
        "14a_local_roads_minus_cold_boundaries",
        cfg.road_secondary_width_m / 2.0,
    )

    subtract_layer(
        output_gpkg,
        "14a_local_roads_minus_cold_boundaries",
        "04_arterial_roads",
        output_gpkg,
        "14a_local_roads_minus_cold_boundaries",
        cfg.road_arterial_width_m / 2.0,
    )

    subtract_layer(
        output_gpkg,
        local_roads_layer_name,
        "05_secondary_roads",
        output_gpkg,
        local_roads_layer_name,
        cfg.road_secondary_width_m / 2.0,
    )

    subtract_layer(
        output_gpkg,
        local_roads_layer_name,
        "04_arterial_roads",
        output_gpkg,
        local_roads_layer_name,
        cfg.road_arterial_width_m / 2.0,
    )

    print("Step 17: Merging all grid layers with grid_type information...")
    all_grid_layer_name = "17_all_grids_merged"
    merge_grid_layers_with_type(
        str(output_gpkg),
        str(output_gpkg),
        all_grid_layer_name,
        [
            (off_grid_inner_layer, "off_grid"),
            (on_grid_arterial_inner_layer, "on_grid_art"),
            (on_grid_secondary_inner_layer, "on_grid_sec"),
            ("14a_local_roads_minus_cold_boundaries", "cold_boundaries"),
            ("04_arterial_roads", "road_arterial"),
            ("05_secondary_roads", "road_secondary"),
            (local_roads_layer_name, "road_local"),
        ],
    )

    # --------------------------------------
    # Extract local roads from merged grids
    # --------------------------------------
    local_road_grid_layer_name = "17_local_road_grid"
    extract_by_expression(
        output_path,
        all_grid_layer_name,
        (
            "type = 'on_grid_art_local_streets' OR "
            "type = 'on_grid_sec_local_streets' OR "
            "type = 'on_grid_art_local_streets' OR "
            "type = 'off_grid_local_streets' "
        ),
        output_path,
        local_road_grid_layer_name,
    )
    buildable_zone_layer_name = "17_buildable_zone"
    extract_by_expression(
        output_path,
        all_grid_layer_name,
        "zone_type = 'buildable'",
        output_path,
        buildable_zone_layer_name,
    )

    print("Cutting local roads with buildable zone to get roads outside buildable area...")
    erase_layer(
        output_path,
        local_road_grid_layer_name,
        output_path,
        buildable_zone_layer_name,
        output_path,
        "17_local_roads_buffer",
    )

    print("Step 18: Exporting merged grids to GeoJSON...")
    output_geojson = output_dir / "outputs.geojson"
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

    print("Step 19: Generating financial data")
    FinancialStreet(config=cfg)

    return output_gpkg
