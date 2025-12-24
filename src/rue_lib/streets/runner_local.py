# src/rue_lib/streets/runner_local.py
from pathlib import Path

import geopandas as gpd
from osgeo import ogr

from rue_lib.core.geometry import get_utm_zone_from_layer, reproject_layer
from rue_lib.site.roads import buffer_roads
from rue_lib.streets.financial import FinancialStreet
from rue_lib.streets.grids import generate_local_streets
from rue_lib.streets.runner_utils import polygons_to_lines_graph_based, subtract_layer

from .config import StreetConfig
from .operations import (
    export_layer_to_geojson,
    extract_by_expression,
    merge_grid_layers_with_type,
)
from .runner_local_operations import (
    build_base_polygons,
    classify_blocks_by_road_type,
    create_dead_end_lines,
)


def generate_streets_with_local_roads(cfg: StreetConfig, local_roads_geojson: str) -> Path:
    """Run the standard streets workflow, then ingest provided local roads GeoJSON and buffer it.

    Args:
        cfg: Street configuration (same as generate_streets)
        local_roads_geojson: Path to a GeoJSON containing local street centerlines

    Returns:
        StreetOutputs including buffered local streets GeoJSON path.
    """
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_gpkg = output_dir / "outputs.gpkg"
    output_path = str(output_gpkg)

    if output_gpkg.exists():
        output_gpkg.unlink()

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

    print("Step 6: Subtract roads buffer from site...")
    site_minus_roads_layer = "06_site_minus_roads"
    subtract_layer(
        output_path,
        site_layer_name,
        input_roads_buffer_layer_name,
        output_path,
        site_minus_roads_layer,
    )

    print("Step 7: Creating dead-end lines...")
    dead_end_lines_layer = "07_dead_end_lines"
    create_dead_end_lines(
        output_path,
        site_minus_roads_layer,
        input_roads_buffer_layer_name,
        dead_end_lines_layer,
        buffer_distance=0.1,
    )

    print("\nLoading provided local streets GeoJSON for buffering...")
    gdf_local = gpd.read_file(local_roads_geojson)
    if gdf_local.empty:
        raise RuntimeError(f"No geometries found in local_roads_geojson: {local_roads_geojson}")

    if utm_epsg and gdf_local.crs and utm_epsg != gdf_local.crs:
        print(f"  Reprojecting local roads to match streets CRS: {utm_epsg}")
        gdf_local = gdf_local.to_crs(utm_epsg)

    local_roads_layer = "06_updated_local_roads"
    gdf_local.to_file(output_gpkg, layer=local_roads_layer, driver="GPKG")

    # Step 7: build base polygons by merging site-minus-roads boundary with provided local roads
    print("Step 7: Creating base polygons from site minus roads and local roads...")
    base_polygons_layer = "07_base_polygons"
    build_base_polygons(
        output_gpkg,
        site_minus_roads_layer=site_layer_name,
        local_roads_layer=local_roads_layer,
        output_layer_name=base_polygons_layer,
        line_extension=cfg.road_locals_width_m * 0.75,
    )

    # Step 8: Classify blocks by road type
    print("Step 8: Classifying blocks by adjacent road type...")
    classified_blocks_layer = "08_classified_blocks"
    classify_blocks_by_road_type(
        output_gpkg,
        base_polygons_layer,
        input_roads_buffer_layer_name,
        dead_end_lines_layer,
        site_layer_name,
        classified_blocks_layer,
        road_type_column="road_type",
    )

    print("Step 9: Generating local streets zones")
    on_grid_layer_name = "09_on_grid_cells"
    extract_by_expression(
        output_path,
        classified_blocks_layer,
        "block_type IN ('on_grid_art', 'on_grid_sec', 'on_grid_art_sec')",
        output_path,
        on_grid_layer_name,
    )
    layers_for_local_roads = []
    # Check if layer exists
    on_grid_inner_layer = None
    if ogr.Open(str(output_gpkg)).GetLayerByName(on_grid_layer_name):
        _on_grid_outer_layer, on_grid_inner_layer = generate_local_streets(
            output_gpkg, cfg, on_grid_layer_name
        )
        layers_for_local_roads.append(on_grid_layer_name)
    else:
        _on_grid_outer_layer, on_grid_inner_layer = None, None

    off_grid_layer_name = "10_off_grid_cells"
    extract_by_expression(
        output_path,
        classified_blocks_layer,
        "block_type IN ('off_grid')",
        output_path,
        off_grid_layer_name,
    )
    off_grid_inner_layer = None
    if ogr.Open(str(output_gpkg)).GetLayerByName(off_grid_layer_name):
        _off_grid_outer_layer, off_grid_inner_layer = generate_local_streets(
            output_gpkg, cfg, off_grid_layer_name
        )
        layers_for_local_roads.append(off_grid_layer_name)
    else:
        _off_grid_outer_layer, off_grid_inner_layer = None, None

    local_roads_layer_name = "11_local_roads"
    if layers_for_local_roads:
        polygons_to_lines_graph_based(
            output_gpkg,
            layers_for_local_roads,
            output_gpkg,
            local_roads_layer_name,
            merge_distance=1,
        )
        subtract_layer(
            output_path,
            local_roads_layer_name,
            input_roads_buffer_layer_name,
            output_path,
            local_roads_layer_name,
        )
    else:
        print("  Warning: No valid layers for local roads generation, skipping...")

    cold_boundary_layer = "12_cold_boundary"
    extract_by_expression(
        output_path,
        classified_blocks_layer,
        "block_type IN ('cold_boundary')",
        output_path,
        cold_boundary_layer,
    )

    cold_boundary_subtracted_layer = "13_cold_boundary_subtracted"
    if ogr.Open(str(output_gpkg)).GetLayerByName(local_roads_layer_name):
        subtract_layer(
            output_gpkg,
            cold_boundary_layer,
            local_roads_layer_name,
            output_gpkg,
            cold_boundary_subtracted_layer,
            cfg.road_locals_width_m / 2.0,
        )
    else:
        cold_boundary_subtracted_layer = cold_boundary_layer

    subtract_layer(
        output_gpkg,
        cold_boundary_subtracted_layer,
        "05_secondary_roads",
        output_gpkg,
        cold_boundary_subtracted_layer,
        cfg.road_secondary_width_m / 2.0,
    )

    subtract_layer(
        output_gpkg,
        cold_boundary_layer,
        "04_arterial_roads",
        output_gpkg,
        cold_boundary_layer,
        cfg.road_arterial_width_m / 2.0,
    )

    layers_to_merge = []

    if off_grid_inner_layer is not None:
        layers_to_merge.append((off_grid_inner_layer, "off_grid"))
    if on_grid_inner_layer is not None:
        on_grid_art_layer = "14_on_grid_art"
        extract_by_expression(
            output_path,
            on_grid_inner_layer,
            "block_type IN ('on_grid_art', 'on_grid_art_sec')",
            output_path,
            on_grid_art_layer,
        )

        on_grid_sec_layer = "14_on_grid_sec"
        extract_by_expression(
            output_path,
            on_grid_inner_layer,
            "block_type IN ('on_grid_sec')",
            output_path,
            on_grid_sec_layer,
        )
        layers_to_merge.append((on_grid_art_layer, "on_grid_art"))
        layers_to_merge.append((on_grid_sec_layer, "on_grid_sec"))

    layers_to_merge.append((cold_boundary_subtracted_layer, "cold_boundaries"))

    layers_to_merge.extend(
        [
            ("04_arterial_roads", "road_arterial"),
            ("05_secondary_roads", "road_secondary"),
            (local_roads_layer_name, "road_local"),
        ]
    )
    all_grid_layer_name = "14_all_grids_merged"
    merge_grid_layers_with_type(
        str(output_gpkg),
        str(output_gpkg),
        all_grid_layer_name,
        layers_to_merge,
    )

    # Layers for financial calculations
    buildable_zone_layer_name = "15_buildable_zone"
    extract_by_expression(
        output_path,
        all_grid_layer_name,
        "zone_type = 'buildable'",
        output_path,
        buildable_zone_layer_name,
    )
    buffered_local_roads_layer_name = "16_local_roads_buffer"
    subtract_layer(
        output_path,
        site_minus_roads_layer,
        cold_boundary_subtracted_layer,
        output_path,
        buffered_local_roads_layer_name,
    )
    if on_grid_inner_layer:
        subtract_layer(
            output_path,
            buffered_local_roads_layer_name,
            on_grid_inner_layer,
            output_path,
            buffered_local_roads_layer_name,
        )
    if off_grid_inner_layer:
        subtract_layer(
            output_path,
            buffered_local_roads_layer_name,
            off_grid_inner_layer,
            output_path,
            buffered_local_roads_layer_name,
        )

    print("Step 18: Exporting merged grids to GeoJSON...")
    output_geojson = output_dir / "outputs.geojson"
    export_layer_to_geojson(
        str(output_gpkg),
        all_grid_layer_name,
        str(output_geojson),
    )

    print("Exporting local streets to GeoJSON for editor use...")
    local_streets_geojson = output_dir / "local_streets.geojson"
    export_layer_to_geojson(
        str(output_gpkg),
        local_roads_layer_name,
        str(local_streets_geojson),
    )

    # Overwrite the local streets GeoJSON with the provided input (projected data)
    local_streets_geojson = Path(output_gpkg).parent / "local_streets.geojson"
    gdf_local.to_file(local_streets_geojson, driver="GeoJSON")

    print(f"\nProcessing complete! Output saved to: {output_gpkg}")
    print("\nFinal layers:")
    print(f"  - {site_layer}: Site polygon with roads subtracted")
    print(f"  - {all_grid_layer_name}: Merged grid cells with grid_type classification")
    print("\nGeoJSON export:")
    print(f"  - {output_geojson}: Merged grids with grid_type classification")
    print(f"  - {local_streets_geojson}: Local streets centerlines")

    print("Step 17: Generating financial data")
    FinancialStreet(
        config=cfg,
        local_roads_layer=local_roads_layer_name,
        local_roads_buffer_layer=buffered_local_roads_layer_name,
    )

    return output_gpkg
