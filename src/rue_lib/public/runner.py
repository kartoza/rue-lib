# src/rue_lib/public/runner.py
from pathlib import Path

from osgeo import gdal, ogr

from .config import PublicConfig
from .operations import (
    allocate_amenities,
    allocate_open_spaces, merge,
)
from ..streets.operations import export_layer_to_geojson

gdal.UseExceptions()


def generate_public(cfg: PublicConfig) -> Path:
    """
    Generate public spaces from input data.

    Returns:
        Path to output file
    """
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_gpkg = output_dir / "outputs.gpkg"
    output_path = str(output_gpkg)

    # Create output GeoPackage if it doesn't exist
    if not output_gpkg.exists():
        driver = ogr.GetDriverByName("GPKG")
        ds = driver.CreateDataSource(output_path)
        if ds is None:
            raise ValueError(f"Could not create {output_path}")
        ds = None

    print("Step 1: Copying input blocks layer...")
    clusters_layer = "111_warm_block_final"

    # Copy input blocks to output geopackage
    input_ds = ogr.Open(cfg.input_path, 0)
    input_layer = input_ds.GetLayer()
    if input_layer is None:
        raise ValueError(
            f"Layer {clusters_layer} not found in {cfg.input_path}"
        )

    output_ds = ogr.Open(output_path, 1)
    if output_ds is None:
        raise ValueError(f"Could not open {output_path} for writing")

    # Remove layer if it already exists
    if output_ds.GetLayerByName("00_input_blocks"):
        output_ds.DeleteLayer("00_input_blocks")

    # Copy the layer
    output_ds.CopyLayer(input_layer, "00_input_blocks")
    print(f"  Copied {input_layer.GetFeatureCount()} blocks to output")

    # ---------------------------
    # Site
    # ---------------------------
    site_ds = ogr.Open(cfg.site_path, 0)
    site_layer = site_ds.GetLayer()
    # Remove layer if it already exists
    if output_ds.GetLayerByName("00_site_path"):
        output_ds.DeleteLayer("00_site_path")

    # Copy the layer
    output_ds.CopyLayer(site_layer, "00_site_path")
    print(
        f"  Copied site layer with {site_layer.GetFeatureCount()} features to output"
    )

    # Copy the layer
    final_layer_name = "04_final"
    if output_ds.GetLayerByName(final_layer_name):
        output_ds.DeleteLayer(final_layer_name)
    output_ds.CopyLayer(input_layer, final_layer_name)

    input_ds = None
    output_ds = None

    print("\nStep 2: Allocating open spaces...")

    open_spaces_layer_name = allocate_open_spaces(
        input_path=cfg.input_path,
        parcels_path=cfg.site_path,
        output_gpkg=output_path,
        output_layer_name="02_open_spaces",
        open_percent=cfg.open_percent,
    )

    print("\nStep 3: Allocating amenities...")

    amenities_layer_name = allocate_amenities(
        input_path=cfg.input_path,
        open_spaces_layer_name=open_spaces_layer_name,
        parcels_path=cfg.site_path,
        output_gpkg=output_path,
        output_layer_name="03_amenities",
        amen_percent=cfg.amen_percent,
    )

    merge(
        input_gpkg=output_path,
        final_layer_name=final_layer_name,
        open_spaces_layer_name=open_spaces_layer_name,
        amenities_layer_name=amenities_layer_name,
    )
    print(f"\nProcessing complete! Output saved to: {output_gpkg}")

    print("Step 18: Exporting merged grids to GeoJSON...")
    output_geojson = output_dir / "outputs.geojson"
    export_layer_to_geojson(
        str(output_gpkg),
        final_layer_name,
        str(output_geojson),
    )

    return output_gpkg
