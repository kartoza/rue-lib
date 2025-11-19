"""Basic geometry operations and cleaning utilities."""

from __future__ import annotations

import math
import os

from osgeo import ogr, osr
from shapely.geometry import Polygon

from rue_lib.core.exceptions import GeometryError


def clean_edges(poly: Polygon) -> Polygon:
    """
    Fix self-intersections and dangles using buffer(0).

    Args:
        poly: Input polygon

    Returns:
        Cleaned polygon

    Raises:
        GeometryError: If the polygon is empty or invalid
    """
    if poly.is_empty:
        raise GeometryError("Cannot clean edges of an empty polygon")
    try:
        return poly.buffer(0)
    except Exception as e:
        raise GeometryError(f"Failed to clean edges: {e}") from e


def clean_angles(poly: Polygon, eps: float) -> Polygon:
    """
    Remove tiny spikes and near-collinear vertices.

    Args:
        poly: Input polygon
        eps: Simplification tolerance in meters

    Returns:
        Simplified polygon
    """
    try:
        simplified = poly.simplify(eps, preserve_topology=True)
        return simplified.buffer(0)
    except Exception as e:
        raise GeometryError(f"Failed to clean angles: {e}") from e


def bounce(poly: Polygon, dist: float) -> Polygon:
    """
    Offset inwards and back out to tidy jaggies (miter/square caps).

    Args:
        poly: Input polygon
        dist: Offset distance in meters

    Returns:
        Bounced polygon
    """
    if dist <= 0:
        return poly

    try:
        return poly.buffer(-dist, join_style=2, cap_style=2).buffer(dist, join_style=2, cap_style=2)
    except Exception as e:
        raise GeometryError(f"Failed to bounce polygon: {e}") from e


def ring_clean_short_edges(
    coords: list[tuple[float, float]], min_len: float
) -> list[tuple[float, float]]:
    """
    Remove vertices that create edges shorter than min_len.

    Args:
        coords: Ring coordinates (should be closed)
        min_len: Minimum edge length in meters

    Returns:
        Cleaned coordinate list
    """
    if len(coords) < 4:
        return coords

    # Ensure closure
    if coords[0] != coords[-1]:
        coords = coords + [coords[0]]

    out: list[tuple[float, float]] = [coords[0]]

    for i in range(1, len(coords) - 1):
        x0, y0 = out[-1]
        x1, y1 = coords[i]
        seg_len = math.hypot(x1 - x0, y1 - y0)

        if seg_len >= min_len:
            out.append((x1, y1))

    # Close the ring
    if out[0] != out[-1]:
        out.append(out[0])

    if len(out) < 4:
        return coords

    return out


def ring_clean_collinear(
    coords: list[tuple[float, float]], dot_thresh: float = 0.9999
) -> list[tuple[float, float]]:
    """
    Remove nearly collinear vertices using dot product threshold.

    Args:
        coords: Ring coordinates (should be closed)
        dot_thresh: Collinearity threshold (higher = more aggressive)

    Returns:
        Cleaned coordinate list
    """
    if len(coords) < 4:
        return coords

    if coords[0] != coords[-1]:
        coords = coords + [coords[0]]

    keep: list[tuple[float, float]] = []
    n = len(coords) - 1

    for i in range(n):
        p_prev = coords[i - 1 if i > 0 else n - 1]
        p = coords[i]
        p_next = coords[(i + 1) % n]

        v0 = (p[0] - p_prev[0], p[1] - p_prev[1])
        v1 = (p_next[0] - p[0], p_next[1] - p[1])

        n0 = math.hypot(v0[0], v0[1])
        n1 = math.hypot(v1[0], v1[1])

        if n0 == 0 or n1 == 0:
            continue

        v0 = (v0[0] / n0, v0[1] / n0)
        v1 = (v1[0] / n1, v1[1] / n1)

        dot = abs(v0[0] * v1[0] + v0[1] * v1[1])

        if dot < dot_thresh:
            keep.append(p)

    if not keep:
        return coords

    if keep[0] != keep[-1]:
        keep.append(keep[0])

    if len(keep) < 4:
        return coords

    return keep


def get_utm_zone_from_layer(layer):
    """Return the WGS84 UTM EPSG code for a layer based on its extent center.

    Computes the longitude/latitude of the layer's bounding-box center, maps
    longitude to a UTM zone in [1..60], and returns EPSG 326## for the northern
    hemisphere or 327## for the southern hemisphere.

    Args:
        layer (ogr.Layer): OGR layer; its extent is used to estimate the zone.

    Returns:
        int: EPSG code for the suggested UTM CRS (e.g., 32633 or 32736).

    Notes:
        - Uses a simple zone formula; special-case rules (Norway/Svalbard) are
          not handled.
        - Assumes coordinates are in degrees (lon/lat). If the layer is already
          projected, compute from a geographic layer instead.
    """
    extent = layer.GetExtent()
    lon_center = (extent[0] + extent[1]) / 2
    lat_center = (extent[2] + extent[3]) / 2

    zone = int((lon_center + 180) / 6) + 1

    if lat_center >= 0:
        epsg_code = 32600 + zone
    else:
        epsg_code = 32700 + zone

    return epsg_code


def reproject_layer(input_path, output_path, target_epsg):
    """Reproject a vector layer to a target CRS using GDAL/OGR.

    Args:
        input_path (str): Path to the source dataset (e.g., .geojson, .gpkg).
            The function reads the first (default) layer from this dataset.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if it
            does not exist.
        target_epsg (int): EPSG code for the target CRS (e.g., 32651 for WGS84/
            UTM zone 51N).

    Returns:
        str: The name of the created layer inside the output GeoPackage.

    Raises:
        Exception: Propagated GDAL/OGR errors if opening/creating datasets,
            creating layers, or transforming geometries fails. (GDAL exceptions
            are enabled via `gdal.UseExceptions()`.)
    """
    source_ds = ogr.Open(input_path)
    source_layer = source_ds.GetLayer()

    source_srs = source_layer.GetSpatialRef()
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(target_epsg)

    transform = osr.CoordinateTransformation(source_srs, target_srs)

    layer_name = os.path.splitext(os.path.basename(input_path))[0] + f"_{target_epsg}"

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    for i in range(output_ds.GetLayerCount()):
        if output_ds.GetLayerByIndex(i).GetName() == layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(layer_name, target_srs, source_layer.GetGeomType())

    source_layer_defn = source_layer.GetLayerDefn()
    for i in range(source_layer_defn.GetFieldCount()):
        field_defn = source_layer_defn.GetFieldDefn(i)
        output_layer.CreateField(field_defn)

    for feature in source_layer:
        geom = feature.GetGeometryRef()
        geom.Transform(transform)

        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)

        for i in range(source_layer_defn.GetFieldCount()):
            out_feature.SetField(
                source_layer_defn.GetFieldDefn(i).GetNameRef(), feature.GetField(i)
            )

        output_layer.CreateFeature(out_feature)
        out_feature = None

    source_ds = None
    output_ds = None

    return layer_name


def buffer_layer(input_path, layer_name, distance, output_path, output_layer_name, dissolve=True):
    """Buffer a vector layer and optionally dissolve the result.

    Args:
        input_path (str): Path to the input dataset (e.g., .gpkg, .geojson).
        layer_name (str): Name of the layer within `input_path` to buffer.
        distance (float): Buffer distance in layer units (e.g., meters if projected).
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output polygon layer to create.
        dissolve (bool): If True, dissolve all buffered features into one geometry.

    Returns:
        None
    """
    from shapely import wkb as shapely_wkb

    source_ds = ogr.Open(input_path)
    source_layer = source_ds.GetLayerByName(layer_name)

    srs = source_layer.GetSpatialRef()

    geoms = []
    for feature in source_layer:
        geom = feature.GetGeometryRef().Clone()

        # Convert to Shapely for sharp-edged buffers
        wkb_data = geom.ExportToWkb()
        if isinstance(wkb_data, bytearray):
            wkb_data = bytes(wkb_data)
        shapely_geom = shapely_wkb.loads(wkb_data)

        # Buffer with sharp edges: join_style=2 (mitre), cap_style=2 (square)
        buffered_shapely = shapely_geom.buffer(distance, join_style=2, cap_style=2)

        # Convert back to OGR
        buffered = ogr.CreateGeometryFromWkb(buffered_shapely.wkb)
        geoms.append(buffered)

    source_ds = None

    if dissolve and geoms:
        union_geom = geoms[0]
        for geom in geoms[1:]:
            union_geom = union_geom.Union(geom)
        geoms = [union_geom]

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    for i in range(output_ds.GetLayerCount()):
        if output_ds.GetLayerByIndex(i).GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)

    for geom in geoms:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)
        output_layer.CreateFeature(out_feature)
        out_feature = None

    output_ds = None
