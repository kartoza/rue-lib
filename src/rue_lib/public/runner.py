# src/rue_lib/public/runner.py
from pathlib import Path

from osgeo import gdal, ogr

from ..core.geometry import get_utm_zone_from_layer, reproject_layer
from ..streets.operations import export_layer_to_geojson
from .config import PublicConfig
from .financial import FinancialPublic
from .open_space import allocate_open_spaces
from .operations import (
    allocate_amenities,
    allocate_cluster_index,
    merge,
)

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

    # Copy input blocks to output geopackage
    input_ds = ogr.Open(cfg.input_path, 0)
    input_layer = input_ds.GetLayer()
    output_ds = ogr.Open(output_path, 1)
    if output_ds is None:
        raise ValueError(f"Could not open {output_path} for writing")

    # Remove layer if it already exists
    input_blocks_layer_name = "00_input_blocks"
    if output_ds.GetLayerByName(input_blocks_layer_name):
        output_ds.DeleteLayer(input_blocks_layer_name)

    utm_epsg = get_utm_zone_from_layer(input_layer)
    reproject_layer(
        cfg.input_path,
        output_path,
        utm_epsg,
        layer_name=input_blocks_layer_name,
        is_append_epsg=False,
    )

    # ---------------------------
    # Site
    # ---------------------------
    input_site_layer_name = "00_input_site"
    # Remove layer if it already exists
    if output_ds.GetLayerByName(input_site_layer_name):
        output_ds.DeleteLayer(input_site_layer_name)

    reproject_layer(
        cfg.site_path, output_path, utm_epsg, layer_name=input_site_layer_name, is_append_epsg=False
    )

    # Copy the layer
    final_layer_name = "04_final"
    if output_ds.GetLayerByName(final_layer_name):
        output_ds.DeleteLayer(final_layer_name)
    reproject_layer(
        cfg.input_path, output_path, utm_epsg, layer_name=final_layer_name, is_append_epsg=False
    )

    input_ds = None
    output_ds = None

    print("\nStep 1: Allocating cluster index for all clusters...")
    allocate_cluster_index(output_path, input_blocks_layer_name)
    allocate_cluster_index(output_path, final_layer_name)

    print("\nStep 2: Allocating open spaces...")

    open_spaces_layer_name = allocate_open_spaces(
        output_gpkg=output_path,
        parcel_layer_name=input_site_layer_name,
        block_layer_name=input_blocks_layer_name,
        output_layer_name="02_open_spaces",
        open_percent=cfg.open_percent,
    )

    print("\nStep 3: Allocating amenities...")

    amenities_layer_name = allocate_amenities(
        output_gpkg=output_path,
        parcel_layer_name=input_site_layer_name,
        block_layer_name=input_blocks_layer_name,
        open_spaces_layer_name=open_spaces_layer_name,
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

    print("Step 19: Generating financial data")
    FinancialPublic(config=cfg)

    return output_gpkg
