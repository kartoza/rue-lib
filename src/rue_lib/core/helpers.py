from pathlib import Path

import geopandas as gpd
from osgeo import ogr
from shapely import wkt


def feature_geom_to_shapely(feat):
    """
    Convert a feature geometry (QGIS or OGR) to Shapely.
    """
    if hasattr(feat, "geometry"):
        geom = feat.geometry()
        wkt_str = geom.asWkt() if hasattr(geom, "asWkt") else geom.ExportToWkt()
    else:
        geom = feat.GetGeometryRef()
        wkt_str = geom.ExportToWkt()
    return wkt.loads(wkt_str)


def remove_layer_from_gpkg(gpkg_path: Path, layer_name: str) -> None:
    """Delete a layer from a GeoPackage using OGR, if it exists."""
    ds = ogr.Open(str(gpkg_path), update=1)
    if ds is None:
        print(f"Warning: could not open {gpkg_path} to remove layer '{layer_name}'")
        return

    try:
        layer_count = ds.GetLayerCount()
        layer_index_to_delete = None

        for idx in range(layer_count):
            lyr = ds.GetLayerByIndex(idx)
            if lyr is not None and lyr.GetName() == layer_name:
                layer_index_to_delete = idx
                break

        if layer_index_to_delete is not None:
            res = ds.DeleteLayer(layer_index_to_delete)
            if res != 0:
                print(f"Warning: failed to delete layer '{layer_name}' from {gpkg_path}")
    finally:
        ds = None


def merge_gpkg_layers(
    gpkg_path: Path,
    layer_names: list[str],
    output_layer_name: str,
    add_source_column: bool = False,
    remove_source_layers: bool = False,
) -> str:
    """
    Merge multiple layers from a GeoPackage into a single layer.

    Args:
        gpkg_path: Path to the GeoPackage file
        layer_names: list of layer names to merge
        output_layer_name: Name of the output merged layer
        add_source_column: If True, add a 'source_layer' column indicating origin
        remove_source_layers: If True, delete the source layers after merging

    Returns:
        Name of the output layer

    Example:
        >>> merge_gpkg_layers(
        ...     Path("output.gpkg"),
        ...     ["06_corners", "07_sides", "08_off_grid"],
        ...     "09_all_parts"
        ... )
    """
    gdfs = []

    # Read all layers
    for layer_name in layer_names:
        try:
            gdf = gpd.read_file(gpkg_path, layer=layer_name)

            if gdf.empty:
                print(f"  Skipping empty layer: {layer_name}")
                continue

            # Add source column if requested
            if add_source_column:
                gdf["source_layer"] = layer_name

            gdfs.append(gdf)
            print(f"  Read layer '{layer_name}': {len(gdf)} features")

        except Exception as e:
            print(f"  Warning: Failed to read layer '{layer_name}': {e}")

    if not gdfs:
        print("  No valid layers found to merge")
        return output_layer_name

    # Merge all GeoDataFrames
    merged_gdf = gpd.GeoDataFrame(gpd.pd.concat(gdfs, ignore_index=True), crs=gdfs[0].crs)

    print(
        f"  Merged {len(gdfs)} layers into '{output_layer_name}': {len(merged_gdf)} total features"
    )

    # Write merged layer to GeoPackage
    merged_gdf.to_file(gpkg_path, layer=output_layer_name, driver="GPKG")

    # Optionally remove source layers
    if remove_source_layers:
        for layer_name in layer_names:
            remove_layer_from_gpkg(gpkg_path, layer_name)
            print(f"  Removed source layer: {layer_name}")

    return output_layer_name
