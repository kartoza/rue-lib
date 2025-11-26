from pathlib import Path

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
