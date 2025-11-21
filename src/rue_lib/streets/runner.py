# src/rue_lib/streets/runner.py
from pathlib import Path

from osgeo import gdal, ogr

from rue_lib.core.geometry import buffer_layer, get_utm_zone_from_layer, reproject_layer

from .blocks_orthogonal import (
    clip_site_by_roads,
    create_grid_for_polygons,
)
from .config import StreetConfig
from .operations import (
    break_multipart_features,
    cleanup_grid_blocks,
    clip_layer,
    create_grid_from_on_grid,
    create_local_streets_zone,
    erase_layer,
    extract_by_expression,
)

gdal.UseExceptions()


def generate_on_grid_blocks(output_path: Path, cfg: StreetConfig) -> Path:
    print("Clip arterial setback to site boundary...")
    clip_layer(
        output_path,
        "arterial_setback",
        output_path,
        "site_clipped_by_roads",
        output_path,
        "arterial_setback_clipped",
    )

    print("Clip secondary setback to site boundary...")
    clip_layer(
        output_path,
        "secondary_setback",
        output_path,
        "site_clipped_by_roads",
        output_path,
        "secondary_setback_clipped",
    )

    print("Intersect arterial and secondary setbacks...")
    clip_layer(
        output_path,
        "arterial_setback_clipped",
        output_path,
        "secondary_setback_clipped",
        output_path,
        "intersected_setbacks",
    )

    print("Arterial setback without intersection...")
    erase_layer(
        output_path,
        "arterial_setback_clipped",
        output_path,
        "intersected_setbacks",
        output_path,
        "arterial_setback_final",
    )

    print("Secondary setback without intersection...")
    erase_layer(
        output_path,
        "secondary_setback_clipped",
        output_path,
        "intersected_setbacks",
        output_path,
        "secondary_setback_final",
    )

    print("Breaking arterial setback multipart features...")
    break_multipart_features(
        output_path,
        "arterial_setback_final",
        output_path,
        "arterial_setback_final",
    )

    print("Breaking secondary setback multipart features...")
    break_multipart_features(
        output_path,
        "secondary_setback_final",
        output_path,
        "secondary_setback_final",
    )

    print("Create grid from on-grid arterial setback...")
    create_grid_from_on_grid(
        output_path,
        "arterial_setback_final",
        "arterial_roads",
        cfg.off_grid_partitions_preferred_width,
        output_path,
        "arterial_setback_grid",
        road_buffer_distance=cfg.road_arterial_width_m,
    )

    print("Clean up grid blocks...")
    cleanup_grid_blocks(
        output_path,
        "arterial_setback_grid",
        output_path,
        "arterial_setback_grid_cleaned",
        cfg.arterial_setback_depth * cfg.off_grid_partitions_preferred_width * 0.5,
    )

    print("Create grid from on-grid secondary setback...")
    create_grid_from_on_grid(
        output_path,
        "secondary_setback_final",
        "secondary_roads",
        cfg.off_grid_partitions_preferred_width,
        output_path,
        "secondary_setback_grid",
        intersected_setbacks_layer_name="intersected_setbacks",
        road_buffer_distance=cfg.road_secondary_width_m,
    )

    print("Clean up grid blocks...")
    cleanup_grid_blocks(
        output_path,
        "secondary_setback_grid",
        output_path,
        "secondary_setback_grid_cleaned",
        cfg.secondary_setback_depth * cfg.off_grid_partitions_preferred_width * 0.5,
    )

    return output_path


def generate_local_streets(
    output_path: Path, cfg: StreetConfig, input_layer_name: str
) -> tuple[str, str]:
    """Generate local streets zone from grid blocks.

    Creates a zone for local streets with sidewalks by:
    1. Creating an inner buffer: -(sidewalk_width + road_width/2)
    2. Creating an outer rounded buffer: +sidewalk_width

    Both inner and outer zones are saved as separate layers.

    Args:
        output_path (Path): Path to the GeoPackage containing grid blocks.
        cfg (StreetConfig): Configuration with sidewalk and road width parameters.
        input_layer_name (str): Name of the input grid blocks layer.

    Returns:
        tuple[str, str]: Names of (outer_layer, inner_layer) created.
    """
    output_layer_name = f"{input_layer_name}_local_streets"

    print(f"Creating local streets zone from {input_layer_name}...")

    outer_layer, inner_layer = create_local_streets_zone(
        str(output_path),
        input_layer_name,
        str(output_path),
        output_layer_name,
        cfg.sidewalk_width_m,
        cfg.road_locals_width_m,
    )

    print(f"  Created layers: {outer_layer} (outer), {inner_layer} (inner)")

    return (outer_layer, inner_layer)


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

    print("Step 3: Creating site polygon clipped by roads...")
    clipped_site_layer = clip_site_by_roads(
        output_path,
        site_layer_name,
        roads_layer_name,
        "site_clipped_by_roads",
        cfg.road_arterial_width_m,
        cfg.road_secondary_width_m,
        cfg.road_local_width_m,
    )

    print("Step 4: Extracting arterial roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_art'", output_path, "arterial_roads"
    )

    print("Step 5: Extracting secondary roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_sec'", output_path, "secondary_roads"
    )

    print("Step 6: Creating arterial road setback zone...")
    buffer_layer(
        output_path,
        "arterial_roads",
        cfg.arterial_setback_depth,
        output_path,
        "arterial_setback",
        dissolve=True,
    )

    print("Step 7: Creating secondary road setback zone...")
    buffer_layer(
        output_path,
        "secondary_roads",
        cfg.secondary_setback_depth,
        output_path,
        "secondary_setback",
        dissolve=True,
    )

    print("Step 8: Removing arterial setback from site...")
    erase_layer(
        output_path,
        clipped_site_layer,
        output_path,
        "arterial_setback",
        output_path,
        "site_minus_arterial_setback",
    )

    print("Step 9: Removing secondary setback from site...")
    erase_layer(
        output_path,
        "site_minus_arterial_setback",
        output_path,
        "secondary_setback",
        output_path,
        "site_minus_all_setbacks",
    )

    print("Step 10: Creating grid for each site polygon...")
    if cfg.optimize_grid_rotation:
        search_method = (
            "ternary search"
            if cfg.use_ternary_search
            else f"linear search ({cfg.grid_rotation_angle_step}° step)"
        )
        print(f"  Optimizing grid rotation using {search_method}...")
    grid_layer = create_grid_for_polygons(
        output_path,
        "site_minus_all_setbacks",
        "site_grid",
        cfg.off_grid_partitions_preferred_width,
        cfg.off_grid_partitions_preferred_depth,
        optimize_rotation=cfg.optimize_grid_rotation,
        rotation_angle_step=cfg.grid_rotation_angle_step,
        use_ternary_search=cfg.use_ternary_search,
        clip_to_boundary=cfg.clip_to_boundary,
    )

    print("Step 11: Filter grid for each site polygon...")
    if cfg.optimize_grid_rotation:
        search_method = (
            "ternary search"
            if cfg.use_ternary_search
            else f"linear search ({cfg.grid_rotation_angle_step}° step)"
        )
        print(f"  Optimizing grid rotation using {search_method}...")
    grid_layer_filtered = create_grid_for_polygons(
        output_path,
        "site_minus_all_setbacks",
        "site_grid_filtered",
        cfg.off_grid_partitions_preferred_width,
        cfg.off_grid_partitions_preferred_depth,
        optimize_rotation=cfg.optimize_grid_rotation,
        rotation_angle_step=cfg.grid_rotation_angle_step,
        use_ternary_search=cfg.use_ternary_search,
        clip_to_boundary=False,
        tolerance_area_ratio=cfg.tolerance_area_ratio,
        tolerance_boundary_distance=cfg.tolerance_boundary_distance,
    )

    print("Step 12: Creating residual polygons (site minus grid)...")
    erase_layer(
        output_path,
        "site_minus_all_setbacks",
        output_path,
        grid_layer_filtered,
        output_path,
        "site_residual",
    )

    print("Step 13: Finding additional grid cells in residual areas...")
    print("  Searching for optimal rotations in leftover spaces...")
    create_grid_for_polygons(
        output_path,
        "site_residual",
        "site_residual_grid",
        cfg.off_grid_partitions_preferred_width,
        cfg.off_grid_partitions_preferred_depth,
        optimize_rotation=True,  # Always optimize for residual areas
        rotation_angle_step=cfg.grid_rotation_angle_step,
        use_ternary_search=cfg.use_ternary_search,
        clip_to_boundary=False,  # Only perfect interior cells for residual
    )

    generate_on_grid_blocks(output_gpkg, cfg)

    generate_local_streets(
        output_gpkg,
        cfg,
        "arterial_setback_grid_cleaned",
    )

    generate_local_streets(
        output_gpkg,
        cfg,
        "secondary_setback_grid_cleaned",
    )

    print(f"\nProcessing complete! Output saved to: {output_gpkg}")
    print("\nFinal layers:")
    print(f"  - {clipped_site_layer}: Site polygon with roads subtracted")
    print(
        f"  - {grid_layer}: Grid cells ({cfg.off_grid_partitions_preferred_width}m x "
        f"{cfg.off_grid_partitions_preferred_depth}m)"
    )

    return output_gpkg
