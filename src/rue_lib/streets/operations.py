# src/rue_lib/streets/operations.py
import math
import os

from osgeo import ogr
from shapely import geometry, wkb, wkt
from shapely.errors import GEOSException, WKTReadingError
from shapely.geometry import CAP_STYLE, JOIN_STYLE
from shapely.ops import split, unary_union

from .geometry_utils import (
    extract_edges_from_geom,
    merge_connected_edges,
)


def _ogr_to_shapely(ogr_geom):
    if ogr_geom is None:
        return None
    try:
        wkb_data = ogr_geom.ExportToWkb()
        if isinstance(wkb_data, memoryview):
            wkb_data = bytes(wkb_data)
        if isinstance(wkb_data, (bytes, bytearray)):
            return wkb.loads(wkb_data)
    except (TypeError, GEOSException, AttributeError, RuntimeError):
        # TypeError: Invalid WKB data type
        # GEOSException: GEOS topology errors
        # AttributeError: Invalid OGR geometry object
        # RuntimeError: GDAL/OGR errors
        pass
    try:
        return wkt.loads(ogr_geom.ExportToWkt())
    except (WKTReadingError, GEOSException, AttributeError, RuntimeError):
        # WKTReadingError: Invalid WKT string
        # GEOSException: GEOS topology errors
        # AttributeError: Invalid OGR geometry object
        # RuntimeError: GDAL/OGR errors
        return None


def extract_by_expression(input_path, layer_name, expression, output_path, output_layer_name):
    """Extract features by attribute expression and write to a GeoPackage layer.

    Args:
        input_path (str): Path to the input dataset (e.g., .gpkg, .geojson).
        layer_name (str): Name of the layer within `input_path` to filter.
        expression (str): OGR attribute filter expression
            (e.g., "\"road_type\" = 'road_art'").
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output layer to create.

    Returns:
        None

    Raises:
        Exception: Propagated GDAL/OGR errors (dataset access, layer creation,
            or feature writes).
    """
    source_ds = ogr.Open(input_path)
    source_layer = source_ds.GetLayerByName(layer_name)

    source_layer_defn = source_layer.GetLayerDefn()
    srs = source_layer.GetSpatialRef()
    geom_type = source_layer.GetGeomType()

    field_defs = []
    for i in range(source_layer_defn.GetFieldCount()):
        field_defs.append(source_layer_defn.GetFieldDefn(i))

    source_layer.SetAttributeFilter(expression)
    features_data = []
    for feature in source_layer:
        geom = feature.GetGeometryRef().Clone()
        field_values = {}
        for i in range(source_layer_defn.GetFieldCount()):
            field_name = source_layer_defn.GetFieldDefn(i).GetNameRef()
            field_values[field_name] = feature.GetField(i)
        features_data.append((geom, field_values))

    source_ds = None

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    for i in range(output_ds.GetLayerCount()):
        if output_ds.GetLayerByIndex(i).GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(output_layer_name, srs, geom_type)

    for field_def in field_defs:
        output_layer.CreateField(field_def)

    for geom, field_values in features_data:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)

        for field_name, value in field_values.items():
            out_feature.SetField(field_name, value)

        output_layer.CreateFeature(out_feature)
        out_feature = None

    output_ds = None


def extract_by_geometry_type(input_path, layer_name, geom_types, output_path, output_layer_name):
    """Extract features by attribute expression and geometry type, write to a GeoPackage layer.

    Args:
        input_path (str): Path to the input dataset (e.g., .gpkg, .geojson).
        layer_name (str): Name of the layer within `input_path` to filter.
        geom_types (list): List of geometry type names to include
            (e.g., ['POLYGON', 'MULTIPOLYGON']).
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output layer to create.

    Returns:
        None

    Raises:
        Exception: Propagated GDAL/OGR errors (dataset access, layer creation,
            or feature writes).
    """
    source_ds = ogr.Open(input_path)
    source_layer = source_ds.GetLayerByName(layer_name)

    source_layer_defn = source_layer.GetLayerDefn()
    srs = source_layer.GetSpatialRef()
    geom_type = source_layer.GetGeomType()

    field_defs = []
    for i in range(source_layer_defn.GetFieldCount()):
        field_defs.append(source_layer_defn.GetFieldDefn(i))

    features_data = []
    for feature in source_layer:
        geom = feature.GetGeometryRef().Clone()

        # Filter by geometry type
        geom_type_name = geom.GetGeometryName()
        if geom_type_name not in geom_types:
            continue

        field_values = {}
        for i in range(source_layer_defn.GetFieldCount()):
            field_name = source_layer_defn.GetFieldDefn(i).GetNameRef()
            field_values[field_name] = feature.GetField(i)
        features_data.append((geom, field_values))

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    for i in range(output_ds.GetLayerCount()):
        if output_ds.GetLayerByIndex(i).GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(output_layer_name, srs, geom_type)

    for field_def in field_defs:
        output_layer.CreateField(field_def)

    for geom, field_values in features_data:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)

        for field_name, value in field_values.items():
            out_feature.SetField(field_name, value)

        output_layer.CreateFeature(out_feature)


def clip_layer(
    input_path, input_layer_name, clip_path, clip_layer_name, output_path, output_layer_name
):
    """Clip features of one layer by another and write the result to a GeoPackage.

    Args:
        input_path (str): Path to the dataset containing the input layer.
        input_layer_name (str): Name of the layer to be clipped.
        clip_path (str): Path to the dataset containing the clip layer
            (may be the same as `input_path`).
        clip_layer_name (str): Name of the layer whose geometries define the clip area.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output layer to create.

    Returns:
        None

    Raises:
        Exception: Propagated GDAL/OGR errors (dataset access, layer creation,
            or geometry operations).
    """
    # Open input
    input_ds = ogr.Open(input_path)
    input_layer = input_ds.GetLayerByName(input_layer_name)

    srs = input_layer.GetSpatialRef()
    geom_type = input_layer.GetGeomType()

    input_layer_defn = input_layer.GetLayerDefn()
    field_defs = []
    for i in range(input_layer_defn.GetFieldCount()):
        field_defs.append(input_layer_defn.GetFieldDefn(i))

    input_features = []
    for feature in input_layer:
        geom = feature.GetGeometryRef().Clone()
        field_values = {}
        for i in range(input_layer_defn.GetFieldCount()):
            field_name = input_layer_defn.GetFieldDefn(i).GetNameRef()
            field_values[field_name] = feature.GetField(i)
        input_features.append((geom, field_values))

    input_ds = None

    if clip_path == input_path:
        # Same file, need to be careful
        clip_ds = ogr.Open(clip_path)
    else:
        clip_ds = ogr.Open(clip_path)

    clip_layer = clip_ds.GetLayerByName(clip_layer_name)

    clip_geoms = []
    for feature in clip_layer:
        clip_geoms.append(feature.GetGeometryRef().Clone())

    clip_ds = None

    clip_geom = clip_geoms[0]
    for geom in clip_geoms[1:]:
        clip_geom = clip_geom.Union(geom)

    # Validate and fix the clip geometry
    if not clip_geom.IsValid():
        print("Warning: Clip geometry is invalid, attempting to fix...")
        clip_geom = clip_geom.Buffer(0)

    clipped_features = []
    for geom, field_values in input_features:
        # Validate input geometry before intersection
        if not geom.IsValid():
            geom = geom.Buffer(0)

        clipped_geom = geom.Intersection(clip_geom)

        # Validate resulting geometry
        if not clipped_geom.IsEmpty():
            if not clipped_geom.IsValid():
                clipped_geom = clipped_geom.Buffer(0)

            if not clipped_geom.IsEmpty():
                clipped_features.append((clipped_geom, field_values))

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    # Remove layer if exists
    for i in range(output_ds.GetLayerCount()):
        if output_ds.GetLayerByIndex(i).GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(output_layer_name, srs, geom_type)

    id_field = ogr.FieldDefn("id", ogr.OFTInteger)
    output_layer.CreateField(id_field)

    for field_def in field_defs:
        try:
            output_layer.CreateField(field_def)
        except RuntimeError:
            # Field may already exist (if same name as 'id'), skip
            pass

    feature_id = 1
    for geom, field_values in clipped_features:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)
        out_feature.SetField("id", feature_id)

        for field_name, value in field_values.items():
            out_feature.SetField(field_name, value)

        output_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None


def erase_layer(
    input_path, input_layer_name, erase_path, erase_layer_name, output_path, output_layer_name
):
    """Erase (subtract) features of one layer from another using difference operation.

    Args:
        input_path (str): Path to the dataset containing the input layer.
        input_layer_name (str): Name of the layer to be erased from.
        erase_path (str): Path to the dataset containing the erase layer
            (may be the same as `input_path`).
        erase_layer_name (str): Name of the layer whose geometries will be subtracted.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output layer to create.

    Returns:
        None

    Raises:
        Exception: Propagated GDAL/OGR errors (dataset access, layer creation,
            or geometry operations).
    """
    # Open input
    input_ds = ogr.Open(input_path)
    input_layer = input_ds.GetLayerByName(input_layer_name)

    srs = input_layer.GetSpatialRef()
    geom_type = input_layer.GetGeomType()

    input_layer_defn = input_layer.GetLayerDefn()
    field_defs = []
    for i in range(input_layer_defn.GetFieldCount()):
        field_defs.append(input_layer_defn.GetFieldDefn(i))

    input_features = []
    for feature in input_layer:
        geom = feature.GetGeometryRef().Clone()
        field_values = {}
        for i in range(input_layer_defn.GetFieldCount()):
            field_name = input_layer_defn.GetFieldDefn(i).GetNameRef()
            field_values[field_name] = feature.GetField(i)
        input_features.append((geom, field_values))

    input_ds = None

    if erase_path == input_path:
        erase_ds = ogr.Open(erase_path)
    else:
        erase_ds = ogr.Open(erase_path)

    erase_layer = erase_ds.GetLayerByName(erase_layer_name)

    erase_geoms = []
    for feature in erase_layer:
        erase_geoms.append(feature.GetGeometryRef().Clone())

    erase_ds = None

    # Union all erase geometries
    erase_geom = erase_geoms[0]
    for geom in erase_geoms[1:]:
        erase_geom = erase_geom.Union(geom)

    # Validate and fix the erase geometry
    if not erase_geom.IsValid():
        print("Warning: Erase geometry is invalid, attempting to fix...")
        erase_geom = erase_geom.Buffer(0)

    # Add a small buffer to ensure complete removal
    # This helps handle floating-point precision issues
    buffer_distance = 0.0001  # Small buffer in layer units
    buffered_erase_geom = erase_geom.Buffer(buffer_distance)

    # Perform difference operation (subtract erase geometry from input)
    erased_features = []
    for geom, field_values in input_features:
        # Skip features that don't intersect the erase geometry
        if not geom.Intersects(buffered_erase_geom):
            erased_features.append((geom, field_values))
            continue

        # Validate input geometry before difference
        if not geom.IsValid():
            geom = geom.Buffer(0)

        # Use buffered erase geometry for more complete removal
        erased_geom = geom.Difference(buffered_erase_geom)

        # Validate resulting geometry
        if not erased_geom.IsEmpty():
            if not erased_geom.IsValid():
                erased_geom = erased_geom.Buffer(0)

            # Additional check: remove tiny fragments
            # For linestrings, check length; for polygons, check area
            if not erased_geom.IsEmpty():
                if geom_type in (ogr.wkbLineString, ogr.wkbMultiLineString):
                    # Remove line segments shorter than threshold
                    if erased_geom.Length() > buffer_distance * 2:
                        erased_features.append((erased_geom, field_values))
                else:
                    # For polygons or other geometries, check area
                    if (
                        hasattr(erased_geom, "Area")
                        and erased_geom.Area() > buffer_distance * buffer_distance
                    ):
                        erased_features.append((erased_geom, field_values))
                    elif not hasattr(erased_geom, "Area"):
                        # For point geometries or others without area
                        erased_features.append((erased_geom, field_values))

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    # Remove layer if exists
    for i in range(output_ds.GetLayerCount()):
        if output_ds.GetLayerByIndex(i).GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(output_layer_name, srs, geom_type)

    # Check if 'id' field already exists in the input layer
    has_id_field = any(field_def.GetNameRef() == "id" for field_def in field_defs)

    if not has_id_field:
        id_field = ogr.FieldDefn("id", ogr.OFTInteger)
        output_layer.CreateField(id_field)

    for field_def in field_defs:
        output_layer.CreateField(field_def)

    feature_id = 1
    for geom, field_values in erased_features:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)

        if not has_id_field:
            out_feature.SetField("id", feature_id)
        elif "id" in field_values:
            # Preserve existing id
            pass

        for field_name, value in field_values.items():
            out_feature.SetField(field_name, value)

        output_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None


def merge_layers_without_overlaps(input_path, layer_names, output_path, output_layer_name):
    """Merge polygon layers and dissolve overlaps into a single layer.

    Args:
        input_path (str): Path to the dataset containing the layers to merge.
        layer_names (list[str]): Names of the layers to merge. Layers must exist in `input_path`.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output layer to create.

    Returns:
        str: Name of the created output layer.
    """
    input_path = str(input_path)
    output_path = str(output_path)

    ds = ogr.Open(input_path)
    if ds is None:
        raise RuntimeError(f"Cannot open dataset: {input_path}")

    srs = None
    shapely_geoms = []

    for layer_name in layer_names:
        layer = ds.GetLayerByName(layer_name)
        if layer is None:
            ds = None
            raise RuntimeError(f"Cannot find layer: {layer_name}")

        if srs is None:
            srs = layer.GetSpatialRef()

        for feature in layer:
            geom = feature.GetGeometryRef()
            if geom is None:
                continue

            geom = geom.Clone()
            if not geom.IsValid():
                geom = geom.Buffer(0)

            shapely_geom = wkt.loads(geom.ExportToWkt())
            if not shapely_geom.is_valid:
                shapely_geom = shapely_geom.buffer(0)

            if not shapely_geom.is_empty:
                shapely_geoms.append(shapely_geom)

    ds = None

    if not shapely_geoms:
        raise RuntimeError(f"No geometries found to merge in layers: {layer_names}")

    merged_geom = unary_union(shapely_geoms)
    if merged_geom.is_empty:
        raise RuntimeError("Merged geometry is empty after union operation.")

    if not merged_geom.is_valid:
        merged_geom = merged_geom.buffer(0)

    def _flatten_to_polygons(geom):
        if geom.is_empty:
            return []
        if geom.geom_type == "Polygon":
            return [geom]
        if geom.geom_type == "MultiPolygon":
            return list(geom.geoms)
        if geom.geom_type == "GeometryCollection":
            polygons = []
            for sub_geom in geom.geoms:
                polygons.extend(_flatten_to_polygons(sub_geom))
            return polygons
        return []

    merged_polygons = _flatten_to_polygons(merged_geom)

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    if output_ds.GetLayerByName(output_layer_name):
        output_ds.DeleteLayer(output_layer_name)

    output_layer = output_ds.CreateLayer(output_layer_name, srs, geom_type=ogr.wkbPolygon)

    id_field = ogr.FieldDefn("id", ogr.OFTInteger)
    output_layer.CreateField(id_field)

    feature_id = 1
    for poly in merged_polygons:
        ogr_geom = ogr.CreateGeometryFromWkb(poly.wkb)
        if not ogr_geom.IsValid():
            ogr_geom = ogr_geom.Buffer(0)

        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(ogr_geom)
        out_feature.SetField("id", feature_id)
        output_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None

    return output_layer_name


def create_on_grid_cells_from_perpendiculars(
    input_path,
    setback_layer_name,
    perp_lines_layer_name,
    output_path,
    output_layer_name,
):
    """Create on-grid cells by splitting merged setback polygons with perpendicular lines.

    Lines are only used when their midpoint lies inside the merged setback polygons.
    Each setback polygon is split by its relevant lines, and resulting cells are written
    with basic attribution about the lines that touched them.
    """
    input_path = str(input_path)
    output_path = str(output_path)

    ds = ogr.Open(input_path)
    if ds is None:
        raise RuntimeError(f"Cannot open dataset: {input_path}")

    setback_layer = ds.GetLayerByName(setback_layer_name)
    perp_layer = ds.GetLayerByName(perp_lines_layer_name)

    if setback_layer is None:
        ds = None
        raise RuntimeError(f"Cannot find setback layer: {setback_layer_name}")
    if perp_layer is None:
        ds = None
        raise RuntimeError(f"Cannot find perpendicular lines layer: {perp_lines_layer_name}")

    srs = setback_layer.GetSpatialRef()

    setback_polygons = []
    for feature in setback_layer:
        geom = feature.GetGeometryRef()
        if geom is None:
            continue
        geom = geom.Clone()
        if not geom.IsValid():
            geom = geom.Buffer(0)

        shapely_geom = _ogr_to_shapely(geom)
        if shapely_geom is None:
            continue
        if not shapely_geom.is_valid:
            shapely_geom = shapely_geom.buffer(0)
        if shapely_geom.is_empty:
            continue

        setback_id = feature.GetField("id") if feature.GetFieldIndex("id") >= 0 else None
        setback_polygons.append((setback_id, shapely_geom))

    if not setback_polygons:
        ds = None
        raise RuntimeError("No setback polygons found to create on-grid cells.")

    setback_union = unary_union([geom for _sid, geom in setback_polygons])

    # Collect perpendicular lines whose midpoint falls inside the merged setback area
    candidate_lines = []
    line_counter = 1
    for feature in perp_layer:
        geom = feature.GetGeometryRef()
        if geom is None:
            continue
        geom = geom.Clone()
        if not geom.IsValid():
            geom = geom.Buffer(0)

        shapely_line = _ogr_to_shapely(geom)
        if shapely_line is None:
            continue
        if shapely_line.is_empty:
            continue
        if not shapely_line.is_valid:
            shapely_line = shapely_line.buffer(0)
        if shapely_line.is_empty:
            continue

        midpoint = shapely_line.interpolate(0.5, normalized=True)
        if not setback_union.contains(midpoint):
            continue

        line_id = feature.GetField("id") if feature.GetFieldIndex("id") >= 0 else None
        if line_id is None:
            line_id = line_counter
            line_counter += 1

        candidate_lines.append({"id": line_id, "geom": shapely_line, "midpoint": midpoint})

    ds = None

    if not candidate_lines:
        raise RuntimeError(
            "No perpendicular lines found with midpoint inside the merged setback layer."
        )

    cells = []

    for setback_id, poly in setback_polygons:
        relevant_lines = []
        for line in candidate_lines:
            if poly.contains(line["midpoint"]) and poly.intersects(line["geom"]):
                relevant_lines.append(line)

        if not relevant_lines:
            cells.append(
                {
                    "geom": poly,
                    "setback_id": setback_id,
                    "line_ids": [],
                    "area_m2": poly.area,
                }
            )
            continue

        split_lines = geometry.MultiLineString([line["geom"] for line in relevant_lines])

        try:
            split_result = split(poly, split_lines)
            split_parts = (
                list(split_result.geoms) if hasattr(split_result, "geoms") else [split_result]
            )
        except Exception as exc:
            print(f"Warning: split failed for polygon {setback_id}: {exc}")
            split_parts = [poly]

        for part in split_parts:
            if part.is_empty:
                continue
            if not part.is_valid:
                part = part.buffer(0)
            if part.is_empty:
                continue

            touching_ids = [line["id"] for line in relevant_lines if part.intersects(line["geom"])]

            cells.append(
                {
                    "geom": part,
                    "setback_id": setback_id,
                    "line_ids": touching_ids,
                    "area_m2": part.area,
                }
            )

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    if output_ds is None:
        raise RuntimeError(f"Cannot open output dataset: {output_path}")

    if output_ds.GetLayerByName(output_layer_name):
        output_ds.DeleteLayer(output_layer_name)

    output_layer = output_ds.CreateLayer(output_layer_name, srs, geom_type=ogr.wkbPolygon)

    id_field = ogr.FieldDefn("id", ogr.OFTInteger)
    output_layer.CreateField(id_field)

    setback_id_field = ogr.FieldDefn("setback_id", ogr.OFTInteger)
    output_layer.CreateField(setback_id_field)

    line_count_field = ogr.FieldDefn("line_count", ogr.OFTInteger)
    output_layer.CreateField(line_count_field)

    line_ids_field = ogr.FieldDefn("line_ids", ogr.OFTString)
    line_ids_field.SetWidth(255)
    output_layer.CreateField(line_ids_field)

    area_field = ogr.FieldDefn("area_m2", ogr.OFTReal)
    output_layer.CreateField(area_field)

    feature_id = 1
    for cell in cells:
        geom = cell["geom"]
        ogr_geom = ogr.CreateGeometryFromWkb(geom.wkb)
        if not ogr_geom.IsValid():
            ogr_geom = ogr_geom.Buffer(0)

        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(ogr_geom)
        out_feature.SetField("id", feature_id)

        if cell["setback_id"] is not None:
            out_feature.SetField("setback_id", int(cell["setback_id"]))

        out_feature.SetField("line_count", len(cell["line_ids"]))
        out_feature.SetField("line_ids", ",".join(str(lid) for lid in cell["line_ids"]))
        out_feature.SetField("area_m2", float(cell["area_m2"]))

        output_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None

    return output_layer_name


def classify_on_grid_cells_by_setback(
    input_path,
    cells_layer_name,
    arterial_setback_layer_name,
    output_path,
    arterial_output_layer_name,
    secondary_output_layer_name,
):
    """Classify on-grid cells into arterial and secondary based on setback intersection.

    For each cell:
    - If it intersects arterial setback → classify as arterial (priority)
    - Else → classify as secondary

    Args:
        input_path: Path to GeoPackage with cells and setback layers
        cells_layer_name: Name of the on-grid cells layer to classify
        arterial_setback_layer_name: Name of arterial setback layer
        output_path: Path to output GeoPackage
        arterial_output_layer_name: Name for arterial cells output layer
        secondary_output_layer_name: Name for secondary cells output layer

    Returns:
        Tuple of (arterial_layer_name, secondary_layer_name)
    """
    ds = ogr.Open(input_path, 1 if input_path == output_path else 0)
    if ds is None:
        raise RuntimeError(f"Cannot open dataset: {input_path}")

    cells_layer = ds.GetLayerByName(cells_layer_name)
    arterial_layer = ds.GetLayerByName(arterial_setback_layer_name)

    if cells_layer is None:
        ds = None
        raise RuntimeError(f"Cannot find cells layer: {cells_layer_name}")
    if arterial_layer is None:
        ds = None
        raise RuntimeError(f"Cannot find arterial setback layer: {arterial_setback_layer_name}")

    srs = cells_layer.GetSpatialRef()

    arterial_geoms = []
    for feat in arterial_layer:
        geom = feat.GetGeometryRef()
        if geom is not None:
            shapely_geom = _ogr_to_shapely(geom.Clone())
            if shapely_geom and not shapely_geom.is_empty:
                arterial_geoms.append(shapely_geom)

    arterial_cells = []
    secondary_cells = []

    cells_layer.ResetReading()
    for cell_feat in cells_layer:
        cell_geom_ogr = cell_feat.GetGeometryRef()
        if cell_geom_ogr is None:
            continue

        cell_geom = _ogr_to_shapely(cell_geom_ogr.Clone())
        if cell_geom is None or cell_geom.is_empty:
            continue

        intersects_arterial = any(cell_geom.intersects(art_geom) for art_geom in arterial_geoms)

        cell_data = {
            "geometry": cell_geom_ogr.Clone(),
            "fields": {},
        }

        for i in range(cell_feat.GetFieldCount()):
            field_name = cell_feat.GetFieldDefnRef(i).GetName()
            if cell_feat.IsFieldSet(i):
                cell_data["fields"][field_name] = cell_feat.GetField(i)

        if intersects_arterial:
            arterial_cells.append(cell_data)
        else:
            secondary_cells.append(cell_data)

    if input_path == output_path:
        output_ds = ds
    else:
        output_ds = ogr.Open(output_path, 1)
        if output_ds is None:
            ds = None
            raise RuntimeError(f"Cannot open output dataset: {output_path}")

    cells_layer_defn = cells_layer.GetLayerDefn()

    if output_ds.GetLayerByName(arterial_output_layer_name):
        output_ds.DeleteLayer(arterial_output_layer_name)
    arterial_output_layer = output_ds.CreateLayer(
        arterial_output_layer_name, srs, geom_type=ogr.wkbPolygon
    )

    for i in range(cells_layer_defn.GetFieldCount()):
        field_defn = cells_layer_defn.GetFieldDefn(i)
        arterial_output_layer.CreateField(field_defn)

    for cell_data in arterial_cells:
        feat = ogr.Feature(arterial_output_layer.GetLayerDefn())
        feat.SetGeometry(cell_data["geometry"])
        for field_name, field_value in cell_data["fields"].items():
            try:
                feat.SetField(field_name, field_value)
            except RuntimeError:
                pass
        arterial_output_layer.CreateFeature(feat)
        feat = None

    if output_ds.GetLayerByName(secondary_output_layer_name):
        output_ds.DeleteLayer(secondary_output_layer_name)
    secondary_output_layer = output_ds.CreateLayer(
        secondary_output_layer_name, srs, geom_type=ogr.wkbPolygon
    )

    for i in range(cells_layer_defn.GetFieldCount()):
        field_defn = cells_layer_defn.GetFieldDefn(i)
        secondary_output_layer.CreateField(field_defn)

    for cell_data in secondary_cells:
        feat = ogr.Feature(secondary_output_layer.GetLayerDefn())
        feat.SetGeometry(cell_data["geometry"])
        for field_name, field_value in cell_data["fields"].items():
            try:
                feat.SetField(field_name, field_value)
            except RuntimeError:
                pass
        secondary_output_layer.CreateFeature(feat)
        feat = None

    print(f"  Classified {len(arterial_cells)} arterial cells")
    print(f"  Classified {len(secondary_cells)} secondary cells")

    if output_ds is not None:
        output_ds.FlushCache()
    ds = None
    output_ds = None

    return arterial_output_layer_name, secondary_output_layer_name


def create_local_streets_zone(
    input_path,
    input_layer_name,
    output_path,
    output_layer_name,
    sidewalk_width_m,
    road_width_m,
):
    """Create local streets zone with sidewalks from grid blocks.

    This creates a zone for local streets by:
    1. Creating an inner (negative) buffer using sidewalk_width + half of road_width
    2. Creating an outer (positive) rounded buffer of sidewalk_width from the inner result

    The result represents the area where local streets and sidewalks will be placed.
    Both inner and outer buffer zones are saved as separate layers.

    Args:
        input_path (str): Path to the input dataset containing grid blocks.
        input_layer_name (str): Name of the layer with grid blocks.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Base name for the output layers.
        sidewalk_width_m (float): Width of sidewalk in meters.
        road_width_m (float): Width of local road in meters.

    Returns:
        tuple[str, str]: Names of (outer_layer, inner_layer) created.

    Raises:
        Exception: Propagated GDAL/OGR errors.
    """
    inner_buffer_distance = -(road_width_m / 2.0)
    outer_buffer_distance = sidewalk_width_m

    input_ds = ogr.Open(input_path)
    if input_ds is None:
        raise RuntimeError(f"Could not open input dataset: {input_path}")

    input_layer = input_ds.GetLayerByName(input_layer_name)
    if input_layer is None:
        raise RuntimeError(f"Could not find layer: {input_layer_name}")

    srs = input_layer.GetSpatialRef()

    inner_geoms = []
    outer_geoms = []
    block_types = []

    for feature in input_layer:
        geom = feature.GetGeometryRef()
        if geom is None:
            continue

        shapely_geom = _ogr_to_shapely(geom.Clone())
        if shapely_geom is None or shapely_geom.is_empty:
            continue

        if not shapely_geom.is_valid:
            shapely_geom = shapely_geom.buffer(0)

        inner_buffered = shapely_geom.buffer(
            inner_buffer_distance,
            join_style=JOIN_STYLE.mitre,
            cap_style=CAP_STYLE.square,
        )

        if inner_buffered.is_empty:
            continue

        if not inner_buffered.is_valid:
            inner_buffered = inner_buffered.buffer(0)
        if inner_buffered.is_empty:
            continue

        inner_buffered = inner_buffered.simplify(0.01, preserve_topology=True)
        outer_buffered = inner_buffered.buffer(
            outer_buffer_distance,
            join_style=JOIN_STYLE.round,
            cap_style=CAP_STYLE.round,
        )

        if outer_buffered.is_empty:
            continue

        if not outer_buffered.is_valid:
            outer_buffered = outer_buffered.buffer(0)
        if outer_buffered.is_empty:
            continue

        outer_buffered = outer_buffered.simplify(0.01, preserve_topology=True)

        inner_geoms.append(ogr.CreateGeometryFromWkb(inner_buffered.wkb))
        outer_geoms.append(ogr.CreateGeometryFromWkb(outer_buffered.wkb))

        # Store block_type if it exists in the feature
        block_type = ""
        if feature.GetFieldIndex("block_type") >= 0:
            block_type = feature.GetField("block_type") or ""
        block_types.append(block_type)

    input_ds = None

    if not inner_geoms:
        print(f"No valid geometries after inner buffer from {input_layer_name}")
        return (None, None)

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    inner_layer_name = f"{output_layer_name}_inner"
    outer_layer_name = f"{output_layer_name}_outer"

    for layer_name in [inner_layer_name, outer_layer_name]:
        for i in range(output_ds.GetLayerCount()):
            if output_ds.GetLayerByIndex(i).GetName() == layer_name:
                output_ds.DeleteLayer(i)
                break

    inner_layer = output_ds.CreateLayer(inner_layer_name, srs, ogr.wkbUnknown)
    inner_layer.CreateField(ogr.FieldDefn("area_m2", ogr.OFTReal))
    inner_layer.CreateField(ogr.FieldDefn("sidewalk_area", ogr.OFTReal))
    inner_layer.CreateField(ogr.FieldDefn("buffer_dist", ogr.OFTReal))
    inner_layer.CreateField(ogr.FieldDefn("sidewalk_w", ogr.OFTReal))
    inner_layer.CreateField(ogr.FieldDefn("road_w", ogr.OFTReal))
    inner_layer.CreateField(ogr.FieldDefn("zone_type", ogr.OFTString))
    inner_layer.CreateField(ogr.FieldDefn("block_type", ogr.OFTString))

    outer_layer = output_ds.CreateLayer(outer_layer_name, srs, ogr.wkbUnknown)
    outer_layer.CreateField(ogr.FieldDefn("area_m2", ogr.OFTReal))
    outer_layer.CreateField(ogr.FieldDefn("sidewalk_area", ogr.OFTReal))
    outer_layer.CreateField(ogr.FieldDefn("buffer_dist", ogr.OFTReal))
    outer_layer.CreateField(ogr.FieldDefn("sidewalk_w", ogr.OFTReal))
    outer_layer.CreateField(ogr.FieldDefn("road_w", ogr.OFTReal))
    outer_layer.CreateField(ogr.FieldDefn("zone_type", ogr.OFTString))

    for _idx, (inner_geom, outer_geom, block_type) in enumerate(
        zip(inner_geoms, outer_geoms, block_types)
    ):
        inner_area = inner_geom.GetArea()
        outer_area = outer_geom.GetArea()
        sidewalk_area = outer_area - inner_area

        inner_feature = ogr.Feature(inner_layer.GetLayerDefn())
        inner_feature.SetGeometry(inner_geom)
        inner_feature.SetField("area_m2", inner_area)
        inner_feature.SetField("sidewalk_area", sidewalk_area)
        inner_feature.SetField("buffer_dist", abs(inner_buffer_distance))
        inner_feature.SetField("sidewalk_w", sidewalk_width_m)
        inner_feature.SetField("road_w", road_width_m)
        inner_feature.SetField("zone_type", "buildable")
        inner_feature.SetField("block_type", block_type)
        inner_layer.CreateFeature(inner_feature)
        inner_feature = None

        outer_feature = ogr.Feature(outer_layer.GetLayerDefn())
        outer_feature.SetGeometry(outer_geom)
        outer_feature.SetField("area_m2", outer_area)
        outer_feature.SetField("sidewalk_area", sidewalk_area)
        outer_feature.SetField("buffer_dist", abs(outer_buffer_distance))
        outer_feature.SetField("sidewalk_w", sidewalk_width_m)
        outer_feature.SetField("road_w", road_width_m)
        outer_feature.SetField("zone_type", "street_sidewalk")
        outer_layer.CreateFeature(outer_feature)
        outer_feature = None

    if output_ds is not None:
        output_ds.FlushCache()
    output_ds = None

    return (outer_layer_name, inner_layer_name)


def export_geometry_vertices(
    input_path,
    input_layer_name,
    output_path,
    output_layer_name,
):
    """Export all vertices (points) from geometries in a layer for debugging.

    Args:
        input_path: Path to input GeoPackage
        input_layer_name: Name of layer to extract vertices from
        output_path: Path to output GeoPackage
        output_layer_name: Name for output points layer

    Returns:
        str: Name of created points layer
    """
    input_ds = ogr.Open(input_path)
    if input_ds is None:
        raise RuntimeError(f"Cannot open input dataset: {input_path}")

    input_layer = input_ds.GetLayerByName(input_layer_name)
    if input_layer is None:
        input_ds = None
        raise RuntimeError(f"Cannot find layer: {input_layer_name}")

    srs = input_layer.GetSpatialRef()

    # Collect all vertices
    vertices = []
    for feature in input_layer:
        geom = feature.GetGeometryRef()
        if geom is None:
            continue

        # Extract vertices based on geometry type
        if geom.GetGeometryName() in ["POLYGON", "MULTIPOLYGON"]:
            # Get all rings from polygon(s)
            for i in range(geom.GetGeometryCount()):
                ring = geom.GetGeometryRef(i)
                if ring.GetGeometryName() == "POLYGON":
                    # MultiPolygon case
                    for j in range(ring.GetGeometryCount()):
                        linear_ring = ring.GetGeometryRef(j)
                        for k in range(linear_ring.GetPointCount()):
                            x, y = linear_ring.GetPoint_2D(k)
                            vertices.append((x, y))
                else:
                    # Simple Polygon case - ring is LinearRing
                    for k in range(ring.GetPointCount()):
                        x, y = ring.GetPoint_2D(k)
                        vertices.append((x, y))

    input_ds = None

    if not vertices:
        print(f"No vertices found in layer {input_layer_name}")
        return None

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)
    for i in range(output_ds.GetLayerCount()):
        if output_ds.GetLayerByIndex(i).GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break
    points_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbPoint)
    points_layer.CreateField(ogr.FieldDefn("vertex_id", ogr.OFTInteger))
    for idx, (x, y) in enumerate(vertices):
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(x, y)

        feature = ogr.Feature(points_layer.GetLayerDefn())
        feature.SetGeometry(point)
        feature.SetField("vertex_id", idx)
        points_layer.CreateFeature(feature)
        feature = None

    if output_ds is not None:
        output_ds.FlushCache()
    output_ds = None

    print(f"Exported {len(vertices)} vertices to {output_layer_name}")
    return output_layer_name


def extract_site_boundary_lines(
    input_path,
    site_layer_name,
    arterial_setback_layer_name,
    secondary_setback_layer_name,
    output_path,
    output_layer_name,
):
    """Extract and merge boundary lines from site that touch arterial or secondary setbacks.

    This function extracts the exterior boundary lines of the site polygons,
    identifies which lines touch arterial or secondary setbacks, and merges
    connected edges of the same setback type into continuous linestrings.

    Args:
        input_path (str): Path to the input dataset.
        site_layer_name (str): Name of the site layer (e.g., "arterial_setback_final").
        arterial_setback_layer_name (str): Name of arterial setback layer.
        secondary_setback_layer_name (str): Name of secondary setback layer.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output layer for merged boundary lines.

    Returns:
        str: Name of the created output layer.

    Raises:
        Exception: Propagated GDAL/OGR errors.

    Output Attributes:
        - id: Feature ID
        - length_m: Length of the merged line in meters
        - setback_type: Type of setback ('arterial' or 'secondary')
    """
    from shapely import geometry as shapely_geom
    from shapely import wkb as shapely_wkb
    from shapely.ops import unary_union

    # Open input dataset
    input_ds = ogr.Open(input_path)
    if input_ds is None:
        raise RuntimeError(f"Could not open input dataset: {input_path}")

    site_layer = input_ds.GetLayerByName(site_layer_name)
    if site_layer is None:
        raise RuntimeError(f"Could not find site layer: {site_layer_name}")

    arterial_layer = input_ds.GetLayerByName(arterial_setback_layer_name)
    secondary_layer = input_ds.GetLayerByName(secondary_setback_layer_name)

    srs = site_layer.GetSpatialRef()

    # Load setback geometries as Shapely objects
    arterial_shapely = []
    if arterial_layer is not None:
        for feature in arterial_layer:
            geom = feature.GetGeometryRef()
            if geom is not None:
                wkb_data = geom.ExportToWkb()
                if isinstance(wkb_data, bytearray):
                    wkb_data = bytes(wkb_data)
                arterial_shapely.append(shapely_wkb.loads(wkb_data))

    secondary_shapely = []
    if secondary_layer is not None:
        for feature in secondary_layer:
            geom = feature.GetGeometryRef()
            if geom is not None:
                wkb_data = geom.ExportToWkb()
                if isinstance(wkb_data, bytearray):
                    wkb_data = bytes(wkb_data)
                secondary_shapely.append(shapely_wkb.loads(wkb_data))

    # Dissolve setback geometries for faster distance calculations
    arterial_union = unary_union(arterial_shapely) if arterial_shapely else None
    secondary_union = unary_union(secondary_shapely) if secondary_shapely else None

    # Extract boundary lines from site polygons
    boundary_lines = []
    tolerance = 0.1  # Buffer tolerance in meters

    for site_feature in site_layer:
        site_geom = site_feature.GetGeometryRef()
        if site_geom is None:
            continue

        site_id = site_feature.GetFID()

        # Convert to Shapely
        wkb_data = site_geom.ExportToWkb()
        if isinstance(wkb_data, bytearray):
            wkb_data = bytes(wkb_data)
        site_shapely = shapely_wkb.loads(wkb_data)

        # Get exterior boundary
        if isinstance(site_shapely, shapely_geom.Polygon):
            boundaries = [site_shapely.exterior]
        elif isinstance(site_shapely, shapely_geom.MultiPolygon):
            boundaries = [poly.exterior for poly in site_shapely.geoms]
        else:
            continue

        for ring_idx, boundary in enumerate(boundaries):
            # Split boundary into segments
            coords = list(boundary.coords)
            for seg_idx in range(len(coords) - 1):
                # Create segment as Shapely LineString
                segment = shapely_geom.LineString([coords[seg_idx], coords[seg_idx + 1]])

                # Get the center point of the segment
                center_point = segment.centroid

                # Calculate distances to setbacks using the center point
                min_arterial_dist = (
                    center_point.distance(arterial_union)
                    if arterial_union is not None
                    else float("inf")
                )
                min_secondary_dist = (
                    center_point.distance(secondary_union)
                    if secondary_union is not None
                    else float("inf")
                )

                # Determine if segment is within tolerance of setbacks
                touches_arterial = min_arterial_dist <= tolerance
                touches_secondary = min_secondary_dist <= tolerance

                # Only include segments that touch at least one setback
                if touches_arterial or touches_secondary:
                    # Determine setback type - prioritize by distance
                    if touches_arterial and not touches_secondary:
                        setback_type = "arterial"
                    elif touches_secondary and not touches_arterial:
                        setback_type = "secondary"
                    elif touches_arterial and touches_secondary:
                        # Both within tolerance - check which is closer
                        if min_arterial_dist < min_secondary_dist:
                            setback_type = "arterial"
                        else:
                            setback_type = "secondary"
                    else:
                        # Fallback (shouldn't reach here)
                        setback_type = "arterial" if touches_arterial else "secondary"

                    boundary_lines.append(
                        {
                            "segment": segment,
                            "site_id": site_id,
                            "ring_idx": ring_idx,
                            "segment_idx": seg_idx,
                            "length": segment.length,
                            "setback_type": setback_type,
                            "touches_arterial": touches_arterial,
                            "touches_secondary": touches_secondary,
                            "dist_arterial": min_arterial_dist,
                            "dist_secondary": min_secondary_dist,
                        }
                    )

    input_ds = None

    # Merge connected edges by setback_type
    from shapely.ops import linemerge

    # Group segments by setback_type
    segments_by_type = {"arterial": [], "secondary": []}

    for line_data in boundary_lines:
        setback_type = line_data["setback_type"]
        segment = line_data["segment"]
        segments_by_type[setback_type].append(segment)

    # Merge connected segments for each type
    merged_lines = []
    for setback_type, segments in segments_by_type.items():
        if segments:
            # Use linemerge to connect adjacent segments
            merged = linemerge(segments)

            # Handle both LineString and MultiLineString results
            if merged.geom_type == "LineString":
                merged_lines.append({"geometry": merged, "setback_type": setback_type})
            elif merged.geom_type == "MultiLineString":
                for line in merged.geoms:
                    merged_lines.append({"geometry": line, "setback_type": setback_type})
            else:
                # If merge didn't work, keep original segments
                for seg in segments:
                    merged_lines.append({"geometry": seg, "setback_type": setback_type})

    # Write to output
    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    # Remove existing layer if present
    for i in range(output_ds.GetLayerCount()):
        if output_ds.GetLayerByIndex(i).GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    # Create output layer
    output_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbLineString)

    # Create fields for merged lines
    output_layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("length_m", ogr.OFTReal))
    output_layer.CreateField(ogr.FieldDefn("setback_type", ogr.OFTString))

    # Write merged features
    feature_id = 1
    for line_data in merged_lines:
        geometry = line_data["geometry"]
        setback_type = line_data["setback_type"]

        # Convert Shapely geometry to OGR
        line_ogr = ogr.CreateGeometryFromWkb(geometry.wkb)

        feature = ogr.Feature(output_layer.GetLayerDefn())
        feature.SetGeometry(line_ogr)
        feature.SetField("id", feature_id)
        feature.SetField("length_m", geometry.length)
        feature.SetField("setback_type", setback_type)

        output_layer.CreateFeature(feature)
        feature = None
        feature_id += 1

    output_ds = None

    return output_layer_name


def break_multipart_features(input_path, input_layer_name, output_path, output_layer_name):
    """Break multipart geometries into single-part features.

    Args:
        input_path (str): Path to the dataset containing the input layer.
        input_layer_name (str): Name of the layer with multipart geometries.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output layer to create.

    Returns:
        None

    Raises:
        Exception: Propagated GDAL/OGR errors (dataset access, layer creation,
            or feature writes).
    """
    # Open input
    input_ds = ogr.Open(input_path)
    input_layer = input_ds.GetLayerByName(input_layer_name)

    srs = input_layer.GetSpatialRef()
    geom_type = input_layer.GetGeomType()

    # Get base geometry type (without Multi prefix)
    if geom_type == ogr.wkbMultiPoint:
        base_geom_type = ogr.wkbPoint
    elif geom_type == ogr.wkbMultiLineString:
        base_geom_type = ogr.wkbLineString
    elif geom_type == ogr.wkbMultiPolygon:
        base_geom_type = ogr.wkbPolygon
    else:
        # If not multipart, use original type
        base_geom_type = geom_type

    input_layer_defn = input_layer.GetLayerDefn()
    field_defs = []
    for i in range(input_layer_defn.GetFieldCount()):
        field_defs.append(input_layer_defn.GetFieldDefn(i))

    # Check if 'id' field exists
    has_id_field = any(field_def.GetNameRef() == "id" for field_def in field_defs)

    # Extract features and break multipart geometries
    broken_features = []
    for feature in input_layer:
        geom = feature.GetGeometryRef().Clone()
        field_values = {}
        for i in range(input_layer_defn.GetFieldCount()):
            field_name = input_layer_defn.GetFieldDefn(i).GetNameRef()
            field_values[field_name] = feature.GetField(i)

        # Check if geometry is multipart
        geom_name = geom.GetGeometryName()
        if geom_name.startswith("MULTI"):
            # Break into single parts
            for i in range(geom.GetGeometryCount()):
                single_geom = geom.GetGeometryRef(i).Clone()
                broken_features.append((single_geom, field_values.copy()))
        else:
            # Keep as-is if already single-part
            broken_features.append((geom, field_values))

    input_ds = None

    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    # Remove layer if exists
    for i in range(output_ds.GetLayerCount()):
        if output_ds.GetLayerByIndex(i).GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(output_layer_name, srs, base_geom_type)

    # Create fields
    for field_def in field_defs:
        output_layer.CreateField(field_def)

    # Write features with incrementing id
    feature_id = 1
    for geom, field_values in broken_features:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)

        for field_name, value in field_values.items():
            # Skip setting the id field here, we'll set it separately
            if field_name != "id":
                out_feature.SetField(field_name, value)

        # Set or reset the id field if it exists
        if has_id_field:
            out_feature.SetField("id", feature_id)

        output_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None


def create_arterial_setback_grid(
    input_path,
    setback_layer_name,
    arterial_road_layer_name,
    grid_width,
    output_path,
    output_layer_name,
    arterial_buffer_distance=15,
):
    """
    Generate all polygon edges (as LineStrings) for arterial setback polygons,
    optionally limited to edges that intersect a buffer around the arterial roads.

    For each edge, the distance to the (unbuffered) arterial road union and
    the edge length are stored as attributes.

    Args:
        input_path (str): Path to the dataset containing the setback & arterial layers.
        setback_layer_name (str): Name of the arterial setback layer (polygons).
        arterial_road_layer_name (str): Name of the arterial road layer (lines).
        grid_width (float): Kept for compatibility; currently unused.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output layer to create.
        arterial_buffer_distance (float | None): Optional buffer distance (same
            units as layer CRS). If provided and > 0, only edges that intersect
            the buffered arterial roads are written.

    Returns:
        None

    Raises:
        Exception: Propagated GDAL/OGR errors.
    """
    input_ds = ogr.Open(input_path)
    if input_ds is None:
        raise RuntimeError(f"Could not open input dataset: {input_path}")

    setback_layer = input_ds.GetLayerByName(setback_layer_name)
    if setback_layer is None:
        raise RuntimeError(f"Could not find setback layer: {setback_layer_name}")

    arterial_layer = input_ds.GetLayerByName(arterial_road_layer_name)
    if arterial_layer is None:
        raise RuntimeError(f"Could not find arterial layer: {arterial_road_layer_name}")

    srs = setback_layer.GetSpatialRef()

    arterial_geoms = []
    for feature in arterial_layer:
        geom = feature.GetGeometryRef()
        if geom is None:
            continue
        arterial_geoms.append(geom.Clone())

    if not arterial_geoms:
        raise RuntimeError("No arterial geometries found in arterial layer.")

    arterial_union = arterial_geoms[0]
    for geom in arterial_geoms[1:]:
        arterial_union = arterial_union.Union(geom)

    # Optional buffer geometry for selecting relevant edges
    arterial_buffer_geom = None
    if arterial_buffer_distance is not None and arterial_buffer_distance > 0:
        arterial_buffer_geom = arterial_union.Buffer(arterial_buffer_distance)

    selected_edges = []

    setback_layer.ResetReading()
    for setback_feature in setback_layer:
        setback_geom = setback_feature.GetGeometryRef()
        if setback_geom is None:
            continue

        setback_geom = setback_geom.Clone()

        if not setback_geom.IsValid():
            setback_geom = setback_geom.Buffer(0)

        setback_id = (
            setback_feature.GetField("id") if setback_feature.GetFieldIndex("id") >= 0 else None
        )

        polygon_edges = extract_edges_from_geom(setback_geom)

        for edge_geom in polygon_edges:
            if arterial_buffer_geom is not None and not edge_geom.Within(arterial_buffer_geom):
                continue

            distance_to_road = edge_geom.Distance(arterial_union)
            edge_length = edge_geom.Length()

            selected_edges.append((edge_geom, setback_id, distance_to_road, edge_length))

    merged_edges = merge_connected_edges(selected_edges)

    input_ds = None

    division_points = []
    for edge_geom, setback_id, _distance_to_road, _edge_length in merged_edges:
        n = edge_geom.GetPointCount()
        if n < 2:
            continue
        segment_distances = [0]
        for i in range(n - 1):
            x1, y1 = edge_geom.GetX(i), edge_geom.GetY(i)
            x2, y2 = edge_geom.GetX(i + 1), edge_geom.GetY(i + 1)
            segment_length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            segment_distances.append(segment_distances[-1] + segment_length)

        total_length = segment_distances[-1]

        current_distance = 0
        while current_distance <= total_length:
            # Find which segment this distance falls on
            segment_idx = 0
            for i in range(len(segment_distances) - 1):
                if segment_distances[i] <= current_distance <= segment_distances[i + 1]:
                    segment_idx = i
                    break

            # Calculate position within the segment
            segment_start_dist = segment_distances[segment_idx]
            segment_end_dist = segment_distances[segment_idx + 1]
            segment_length = segment_end_dist - segment_start_dist

            if segment_length > 0:
                t = (current_distance - segment_start_dist) / segment_length
            else:
                t = 0

            # Interpolate position
            x1, y1 = edge_geom.GetX(segment_idx), edge_geom.GetY(segment_idx)
            x2, y2 = edge_geom.GetX(segment_idx + 1), edge_geom.GetY(segment_idx + 1)

            point_x = x1 + t * (x2 - x1)
            point_y = y1 + t * (y2 - y1)

            # Create point geometry
            point_geom = ogr.Geometry(ogr.wkbPoint)
            point_geom.AddPoint(point_x, point_y)

            division_points.append((point_geom, setback_id, current_distance, distance_to_road))

            current_distance += grid_width

        # Always add the end point if not already included
        if current_distance - grid_width < total_length:
            point_geom = ogr.Geometry(ogr.wkbPoint)
            point_geom.AddPoint(edge_geom.GetX(n - 1), edge_geom.GetY(n - 1))
            division_points.append((point_geom, setback_id, total_length, distance_to_road))

    driver = ogr.GetDriverByName("GPKG")
    if driver is None:
        raise RuntimeError("Could not get GPKG driver.")

    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
        if output_ds is None:
            raise RuntimeError(f"Could not open existing output: {output_path}")
    else:
        output_ds = driver.CreateDataSource(output_path)
        if output_ds is None:
            raise RuntimeError(f"Could not create output: {output_path}")

    for i in range(output_ds.GetLayerCount()):
        layer = output_ds.GetLayerByIndex(i)
        if layer.GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbPoint)
    if output_layer is None:
        raise RuntimeError("Could not create output layer.")

    id_field = ogr.FieldDefn("id", ogr.OFTInteger)
    output_layer.CreateField(id_field)

    setback_id_field = ogr.FieldDefn("setback_id", ogr.OFTInteger)
    output_layer.CreateField(setback_id_field)

    distance_field = ogr.FieldDefn("distance", ogr.OFTReal)
    output_layer.CreateField(distance_field)

    dist_to_road_field = ogr.FieldDefn("dist_to_road", ogr.OFTReal)
    output_layer.CreateField(dist_to_road_field)

    feature_id = 1
    for (
        point_geom,
        setback_id,
        distance,
        distance_to_road,
    ) in division_points:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(point_geom)
        out_feature.SetField("id", feature_id)
        if setback_id is not None:
            out_feature.SetField("setback_id", int(setback_id))
        out_feature.SetField("distance", float(distance))
        out_feature.SetField("dist_to_road", float(distance_to_road))

        output_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None


def create_grid_from_division_points(
    input_path,
    setback_layer_name,
    points_layer_name,
    output_path,
    output_layer_name,
):
    """Create grid polygons by splitting setback polygons with perpendicular lines.

    Args:
        input_path (str): Path to the dataset containing layers.
        setback_layer_name (str): Name of the arterial setback layer (polygons).
        points_layer_name (str): Name of the division points layer.
        output_path (str): Path to the output GeoPackage (.gpkg).
        output_layer_name (str): Name of the output layer to create.

    Returns:
        None
    """
    input_ds = ogr.Open(input_path)
    if input_ds is None:
        raise RuntimeError(f"Could not open input dataset: {input_path}")

    setback_layer = input_ds.GetLayerByName(setback_layer_name)
    if setback_layer is None:
        raise RuntimeError(f"Could not find setback layer: {setback_layer_name}")

    points_layer = input_ds.GetLayerByName(points_layer_name)
    if points_layer is None:
        raise RuntimeError(f"Could not find points layer: {points_layer_name}")

    srs = setback_layer.GetSpatialRef()

    # Group points by setback_id
    points_by_setback = {}
    for point_feature in points_layer:
        point_geom = point_feature.GetGeometryRef()
        if point_geom is None:
            continue

        setback_id = point_feature.GetField("setback_id")
        if setback_id not in points_by_setback:
            points_by_setback[setback_id] = []

        point_x = point_geom.GetX()
        point_y = point_geom.GetY()
        distance = point_feature.GetField("distance")
        points_by_setback[setback_id].append((point_x, point_y, distance))

    # Sort points by distance for each setback
    for setback_id in points_by_setback:
        points_by_setback[setback_id].sort(key=lambda p: p[2])

    grid_polygons = []

    # Process each setback polygon
    setback_layer.ResetReading()
    for setback_feature in setback_layer:
        setback_geom = setback_feature.GetGeometryRef()
        if setback_geom is None:
            continue

        setback_geom = setback_geom.Clone()

        if not setback_geom.IsValid():
            setback_geom = setback_geom.Buffer(0)

        setback_id = setback_feature.GetField("id")

        if setback_id not in points_by_setback:
            continue

        points = points_by_setback[setback_id]

        # Get polygon extent for perpendicular line length
        env = setback_geom.GetEnvelope()
        max_extent = max(env[1] - env[0], env[3] - env[2]) * 2

        # Create perpendicular cutting lines from each point
        cutting_lines = []
        for i, (point_x, point_y, _distance) in enumerate(points):
            if i > 0:
                prev_x, prev_y, _ = points[i - 1]
                edge_dx = point_x - prev_x
                edge_dy = point_y - prev_y
            elif i < len(points) - 1:
                next_x, next_y, _ = points[i + 1]
                edge_dx = next_x - point_x
                edge_dy = next_y - point_y
            else:
                # Single point, skip
                continue

            edge_length = math.sqrt(edge_dx * edge_dx + edge_dy * edge_dy)
            if edge_length == 0:
                continue

            # Normalize
            edge_dx /= edge_length
            edge_dy /= edge_length

            # Perpendicular direction (rotate 90 degrees)
            perp_dx = -edge_dy
            perp_dy = edge_dx

            # Create perpendicular line
            perp_line = ogr.Geometry(ogr.wkbLineString)
            perp_line.AddPoint(point_x - max_extent * perp_dx, point_y - max_extent * perp_dy)
            perp_line.AddPoint(point_x + max_extent * perp_dx, point_y + max_extent * perp_dy)

            cutting_lines.append(perp_line)

        # Split the polygon using the cutting lines
        # We'll use a buffer approach to split the polygon
        current_polys = [setback_geom]

        for cutting_line in cutting_lines:
            # Create a thin buffer around the cutting line to act as a cutting blade
            cutting_buffer = cutting_line.Buffer(0.01)  # Small buffer to ensure clean cut

            new_polys = []
            for poly in current_polys:
                # Split the polygon by subtracting the cutting buffer
                try:
                    difference = poly.Difference(cutting_buffer)

                    if difference.IsEmpty():
                        continue

                    # Handle both single and multipart results
                    if difference.GetGeometryName() == "POLYGON":
                        new_polys.append(difference)
                    elif difference.GetGeometryName() == "MULTIPOLYGON":
                        for j in range(difference.GetGeometryCount()):
                            part = difference.GetGeometryRef(j).Clone()
                            if not part.IsEmpty() and part.Area() > 0.1:  # Filter tiny fragments
                                new_polys.append(part)
                    else:
                        # Keep original if split didn't work
                        new_polys.append(poly)
                except Exception as e:
                    # If splitting fails, keep the original polygon
                    print(f"Warning: Failed to split polygon: {e}")
                    new_polys.append(poly)

            current_polys = new_polys

        # Add all resulting polygons to the output
        for poly in current_polys:
            if not poly.IsEmpty() and poly.Area() > 0.1:
                grid_polygons.append((poly, setback_id))

    input_ds = None

    # Write output
    driver = ogr.GetDriverByName("GPKG")
    if driver is None:
        raise RuntimeError("Could not get GPKG driver.")

    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
        if output_ds is None:
            raise RuntimeError(f"Could not open existing output: {output_path}")
    else:
        output_ds = driver.CreateDataSource(output_path)
        if output_ds is None:
            raise RuntimeError(f"Could not create output: {output_path}")

    for i in range(output_ds.GetLayerCount()):
        layer = output_ds.GetLayerByIndex(i)
        if layer.GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)
    if output_layer is None:
        raise RuntimeError("Could not create output layer.")

    id_field = ogr.FieldDefn("id", ogr.OFTInteger)
    output_layer.CreateField(id_field)

    setback_id_field = ogr.FieldDefn("setback_id", ogr.OFTInteger)
    output_layer.CreateField(setback_id_field)

    feature_id = 1
    for polygon_geom, setback_id in grid_polygons:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(polygon_geom)
        out_feature.SetField("id", feature_id)
        if setback_id is not None:
            out_feature.SetField("setback_id", int(setback_id))

        output_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None


def create_grid_from_on_grid(
    input_path,
    setback_layer_name,
    arterial_road_layer_name,
    grid_width,
    output_path,
    output_layer_name,
    intersected_setbacks_layer_name=None,
    road_buffer_distance=15,
):
    """Create grid polygons from arterial setback by division points and splitting.

    This function combines the creation of division points and the splitting of polygons into
    a single operation. It generates division points along edges closest to arterial roads,
    then splits the setback polygons using perpendicular lines from these points.

    If intersected_setbacks is provided, division points start from the corner of the merged
    edge that is closest to the intersected setback, ensuring grids align properly with
    intersection areas.

    Args:
        input_path (str): Path to the dataset containing the setback & arterial layers.
        setback_layer_name (str): Name of the arterial setback layer (polygons).
        arterial_road_layer_name (str): Name of the arterial road layer (lines).
        grid_width (float): Width for division points spacing.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output layer to create.
        intersected_setbacks_layer_name (str | None): Optional name of intersected setbacks layer.
            If provided, division points start from the corner closest to intersected setbacks.
        road_buffer_distance (float | None): Optional buffer distance (same
            units as layer CRS). If provided and > 0, only edges that intersect
            the buffered roads are considered.

    Returns:
        None

    Raises:
        Exception: Propagated GDAL/OGR errors.
    """
    input_ds = ogr.Open(input_path)
    if input_ds is None:
        raise RuntimeError(f"Could not open input dataset: {input_path}")

    setback_layer = input_ds.GetLayerByName(setback_layer_name)
    if setback_layer is None:
        raise RuntimeError(f"Could not find setback layer: {setback_layer_name}")

    arterial_layer = input_ds.GetLayerByName(arterial_road_layer_name)
    if arterial_layer is None:
        raise RuntimeError(f"Could not find arterial layer: {arterial_road_layer_name}")

    srs = setback_layer.GetSpatialRef()

    # Build arterial road union
    arterial_geoms = []
    for feature in arterial_layer:
        geom = feature.GetGeometryRef()
        if geom is None:
            continue
        arterial_geoms.append(geom.Clone())

    if not arterial_geoms:
        raise RuntimeError("No arterial geometries found in arterial layer.")

    arterial_union = arterial_geoms[0]
    for geom in arterial_geoms[1:]:
        arterial_union = arterial_union.Union(geom)

    # Optional buffer geometry for selecting relevant edges
    arterial_buffer_geom = None
    if road_buffer_distance is not None and road_buffer_distance > 0:
        arterial_buffer_geom = arterial_union.Buffer(road_buffer_distance)

    # Load intersected setbacks geometries if provided
    intersected_geoms = []
    if intersected_setbacks_layer_name is not None:
        intersected_layer = input_ds.GetLayerByName(intersected_setbacks_layer_name)
        if intersected_layer is not None:
            for feature in intersected_layer:
                geom = feature.GetGeometryRef()
                if geom is None:
                    continue
                intersected_geoms.append(geom.Clone())

    # Extract and merge edges closest to arterial roads
    selected_edges = []
    setback_layer.ResetReading()
    for setback_feature in setback_layer:
        setback_geom = setback_feature.GetGeometryRef()
        if setback_geom is None:
            continue

        setback_geom = setback_geom.Clone()

        if not setback_geom.IsValid():
            setback_geom = setback_geom.Buffer(0)

        setback_id = (
            setback_feature.GetField("id") if setback_feature.GetFieldIndex("id") >= 0 else None
        )

        polygon_edges = extract_edges_from_geom(setback_geom)

        for edge_geom in polygon_edges:
            if arterial_buffer_geom is not None and not edge_geom.Within(arterial_buffer_geom):
                continue

            distance_to_road = edge_geom.Distance(arterial_union)
            edge_length = edge_geom.Length()

            selected_edges.append((edge_geom, setback_id, distance_to_road, edge_length))

    merged_edges = merge_connected_edges(selected_edges)

    # Write merged edges to debug layer for visualization
    driver = ogr.GetDriverByName("GPKG")
    if driver is None:
        raise RuntimeError("Could not get GPKG driver.")

    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
        if output_ds is None:
            raise RuntimeError(f"Could not open existing output: {output_path}")
    else:
        output_ds = driver.CreateDataSource(output_path)
        if output_ds is None:
            raise RuntimeError(f"Could not create output: {output_path}")

    # Delete existing debug layer if it exists
    debug_layer_name = f"{output_layer_name}_merged_edges_debug"
    for i in range(output_ds.GetLayerCount()):
        layer = output_ds.GetLayerByIndex(i)
        if layer.GetName() == debug_layer_name:
            output_ds.DeleteLayer(i)
            break

    # Create debug layer for merged edges
    debug_layer = output_ds.CreateLayer(debug_layer_name, srs, ogr.wkbLineString)
    if debug_layer is None:
        raise RuntimeError("Could not create debug layer.")

    # Add fields
    debug_layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    debug_layer.CreateField(ogr.FieldDefn("setback_id", ogr.OFTInteger))
    debug_layer.CreateField(ogr.FieldDefn("dist_to_road", ogr.OFTReal))
    debug_layer.CreateField(ogr.FieldDefn("length", ogr.OFTReal))

    # Write merged edges
    feature_id = 1
    for edge_geom, setback_id, distance_to_road, edge_length in merged_edges:
        out_feature = ogr.Feature(debug_layer.GetLayerDefn())
        out_feature.SetGeometry(edge_geom)
        out_feature.SetField("id", feature_id)
        if setback_id is not None:
            out_feature.SetField("setback_id", int(setback_id))
        out_feature.SetField("dist_to_road", distance_to_road)
        out_feature.SetField("length", edge_length)
        debug_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None

    points_by_setback = {}
    for edge_geom, setback_id, _distance_to_road, _edge_length in merged_edges:
        if setback_id not in points_by_setback:
            points_by_setback[setback_id] = []

        n = edge_geom.GetPointCount()
        if n < 2:
            continue

        # Find the closest corner (start or end) to the intersected setback if available
        if intersected_geoms:
            # Get the closest intersected setback for this setback_id
            # We need to find which setback geometry this edge belongs to
            # For now, find the closest intersected geometry to either endpoint
            start_point = ogr.Geometry(ogr.wkbPoint)
            start_point.AddPoint(edge_geom.GetX(0), edge_geom.GetY(0))
            end_point = ogr.Geometry(ogr.wkbPoint)
            end_point.AddPoint(edge_geom.GetX(n - 1), edge_geom.GetY(n - 1))

            min_start_dist = float("inf")
            min_end_dist = float("inf")
            for intersected_geom in intersected_geoms:
                start_dist = start_point.Distance(intersected_geom)
                end_dist = end_point.Distance(intersected_geom)
                if start_dist < min_start_dist:
                    min_start_dist = start_dist
                if end_dist < min_end_dist:
                    min_end_dist = end_dist

            # If the end point is closer to intersected setback, reverse the direction
            # by starting from total_length and going backwards
            reverse_direction = min_end_dist < min_start_dist
        else:
            reverse_direction = False

        # Calculate cumulative distances along the line
        segment_distances = [0]
        for i in range(n - 1):
            x1, y1 = edge_geom.GetX(i), edge_geom.GetY(i)
            x2, y2 = edge_geom.GetX(i + 1), edge_geom.GetY(i + 1)
            segment_length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
            segment_distances.append(segment_distances[-1] + segment_length)

        total_length = segment_distances[-1]

        # Generate division points starting from the closest corner
        if reverse_direction:
            # Start from the end and go backwards
            current_distance = total_length
            while current_distance >= 0:
                # Find which segment this distance falls on
                segment_idx = 0
                for i in range(len(segment_distances) - 1):
                    if segment_distances[i] <= current_distance <= segment_distances[i + 1]:
                        segment_idx = i
                        break

                # Calculate position within the segment
                segment_start_dist = segment_distances[segment_idx]
                segment_end_dist = segment_distances[segment_idx + 1]
                segment_length = segment_end_dist - segment_start_dist

                if segment_length > 0:
                    t = (current_distance - segment_start_dist) / segment_length
                else:
                    t = 0

                # Interpolate position
                x1, y1 = edge_geom.GetX(segment_idx), edge_geom.GetY(segment_idx)
                x2, y2 = edge_geom.GetX(segment_idx + 1), edge_geom.GetY(segment_idx + 1)

                point_x = x1 + t * (x2 - x1)
                point_y = y1 + t * (y2 - y1)

                # Store with distance from the chosen start (reversed)
                dist_from_start = total_length - current_distance
                points_by_setback[setback_id].append((point_x, point_y, dist_from_start))

                current_distance -= grid_width
        else:
            # Start from the beginning (original behavior)
            current_distance = 0
            while current_distance <= total_length:
                # Find which segment this distance falls on
                segment_idx = 0
                for i in range(len(segment_distances) - 1):
                    if segment_distances[i] <= current_distance <= segment_distances[i + 1]:
                        segment_idx = i
                        break

                # Calculate position within the segment
                segment_start_dist = segment_distances[segment_idx]
                segment_end_dist = segment_distances[segment_idx + 1]
                segment_length = segment_end_dist - segment_start_dist

                if segment_length > 0:
                    t = (current_distance - segment_start_dist) / segment_length
                else:
                    t = 0

                # Interpolate position
                x1, y1 = edge_geom.GetX(segment_idx), edge_geom.GetY(segment_idx)
                x2, y2 = edge_geom.GetX(segment_idx + 1), edge_geom.GetY(segment_idx + 1)

                point_x = x1 + t * (x2 - x1)
                point_y = y1 + t * (y2 - y1)

                points_by_setback[setback_id].append((point_x, point_y, current_distance))

                current_distance += grid_width

    # Sort points by distance for each setback
    for setback_id in points_by_setback:
        points_by_setback[setback_id].sort(key=lambda p: p[2])

    # Write division points to debug layer for visualization
    driver = ogr.GetDriverByName("GPKG")
    if driver is None:
        raise RuntimeError("Could not get GPKG driver.")

    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
        if output_ds is None:
            raise RuntimeError(f"Could not open existing output: {output_path}")
    else:
        output_ds = driver.CreateDataSource(output_path)
        if output_ds is None:
            raise RuntimeError(f"Could not create output: {output_path}")

    # Delete existing debug layer if it exists
    points_debug_layer_name = f"{output_layer_name}_division_points_debug"
    for i in range(output_ds.GetLayerCount()):
        layer = output_ds.GetLayerByIndex(i)
        if layer.GetName() == points_debug_layer_name:
            output_ds.DeleteLayer(i)
            break

    # Create debug layer for division points
    points_debug_layer = output_ds.CreateLayer(points_debug_layer_name, srs, ogr.wkbPoint)
    if points_debug_layer is None:
        raise RuntimeError("Could not create points debug layer.")

    # Add fields
    points_debug_layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    points_debug_layer.CreateField(ogr.FieldDefn("setback_id", ogr.OFTInteger))
    points_debug_layer.CreateField(ogr.FieldDefn("distance", ogr.OFTReal))

    # Write division points
    feature_id = 1
    for setback_id, points in points_by_setback.items():
        for point_x, point_y, distance in points:
            point_geom = ogr.Geometry(ogr.wkbPoint)
            point_geom.AddPoint(point_x, point_y)

            out_feature = ogr.Feature(points_debug_layer.GetLayerDefn())
            out_feature.SetGeometry(point_geom)
            out_feature.SetField("id", feature_id)
            if setback_id is not None:
                out_feature.SetField("setback_id", int(setback_id))
            out_feature.SetField("distance", distance)
            points_debug_layer.CreateFeature(out_feature)
            out_feature = None
            feature_id += 1

    output_ds = None

    grid_polygons = []

    # Process each setback polygon and split using division points
    setback_layer.ResetReading()
    for setback_feature in setback_layer:
        setback_geom = setback_feature.GetGeometryRef()
        if setback_geom is None:
            continue

        setback_geom = setback_geom.Clone()

        if not setback_geom.IsValid():
            setback_geom = setback_geom.Buffer(0)

        setback_id = setback_feature.GetField("id")

        if setback_id not in points_by_setback:
            continue

        points = points_by_setback[setback_id]

        # Get polygon extent for perpendicular line length
        env = setback_geom.GetEnvelope()
        max_extent = max(env[1] - env[0], env[3] - env[2]) * 2

        # Create perpendicular cutting lines from each point
        cutting_lines = []
        for i, (point_x, point_y, _distance) in enumerate(points):
            if i > 0:
                prev_x, prev_y, _ = points[i - 1]
                edge_dx = point_x - prev_x
                edge_dy = point_y - prev_y
            elif i < len(points) - 1:
                next_x, next_y, _ = points[i + 1]
                edge_dx = next_x - point_x
                edge_dy = next_y - point_y
            else:
                # Single point, skip
                continue

            edge_length = math.sqrt(edge_dx * edge_dx + edge_dy * edge_dy)
            if edge_length == 0:
                continue

            # Normalize
            edge_dx /= edge_length
            edge_dy /= edge_length

            # Perpendicular direction (rotate 90 degrees)
            perp_dx = -edge_dy
            perp_dy = edge_dx

            # Create perpendicular line
            perp_line = ogr.Geometry(ogr.wkbLineString)
            perp_line.AddPoint(point_x - max_extent * perp_dx, point_y - max_extent * perp_dy)
            perp_line.AddPoint(point_x + max_extent * perp_dx, point_y + max_extent * perp_dy)

            cutting_lines.append(perp_line)

        # Split the polygon using the cutting lines
        current_polys = [setback_geom]

        for cutting_line in cutting_lines:
            # Create a thin buffer around the cutting line to act as a cutting blade
            cutting_buffer = cutting_line.Buffer(0.01)

            new_polys = []
            for poly in current_polys:
                # Split the polygon by subtracting the cutting buffer
                try:
                    difference = poly.Difference(cutting_buffer)

                    if difference.IsEmpty():
                        continue

                    # Handle both single and multipart results
                    if difference.GetGeometryName() == "POLYGON":
                        new_polys.append(difference)
                    elif difference.GetGeometryName() == "MULTIPOLYGON":
                        for j in range(difference.GetGeometryCount()):
                            part = difference.GetGeometryRef(j).Clone()
                            if not part.IsEmpty() and part.Area() > 0.1:
                                new_polys.append(part)
                    else:
                        # Keep original if split didn't work
                        new_polys.append(poly)
                except Exception as e:
                    # If splitting fails, keep the original polygon
                    print(f"Warning: Failed to split polygon: {e}")
                    new_polys.append(poly)

            current_polys = new_polys

        # Add all resulting polygons to the output
        for poly in current_polys:
            if not poly.IsEmpty() and poly.Area() > 0.1:
                grid_polygons.append((poly, setback_id))

    input_ds = None

    # Write output
    driver = ogr.GetDriverByName("GPKG")
    if driver is None:
        raise RuntimeError("Could not get GPKG driver.")

    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
        if output_ds is None:
            raise RuntimeError(f"Could not open existing output: {output_path}")
    else:
        output_ds = driver.CreateDataSource(output_path)
        if output_ds is None:
            raise RuntimeError(f"Could not create output: {output_path}")

    for i in range(output_ds.GetLayerCount()):
        layer = output_ds.GetLayerByIndex(i)
        if layer.GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)
    if output_layer is None:
        raise RuntimeError("Could not create output layer.")

    id_field = ogr.FieldDefn("id", ogr.OFTInteger)
    output_layer.CreateField(id_field)

    setback_id_field = ogr.FieldDefn("setback_id", ogr.OFTInteger)
    output_layer.CreateField(setback_id_field)

    feature_id = 1
    for polygon_geom, setback_id in grid_polygons:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(polygon_geom)
        out_feature.SetField("id", feature_id)
        if setback_id is not None:
            out_feature.SetField("setback_id", int(setback_id))

        output_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None


def cleanup_grid_blocks(
    input_path,
    grid_layer_name,
    output_path,
    output_layer_name,
    min_area,
):
    """Clean up grid blocks by merging small blocks and removing isolated ones.

    Args:
        input_path (str): Path to the dataset containing the grid layer.
        grid_layer_name (str): Name of the grid polygon layer.
        output_path (str): Path to the output GeoPackage (.gpkg).
        output_layer_name (str): Name of the output layer to create.
        min_area (float): Minimum area threshold for blocks.

    Returns:
        None
    """
    input_ds = ogr.Open(input_path)
    if input_ds is None:
        raise RuntimeError(f"Could not open input dataset: {input_path}")

    grid_layer = input_ds.GetLayerByName(grid_layer_name)
    if grid_layer is None:
        raise RuntimeError(f"Could not find grid layer: {grid_layer_name}")

    srs = grid_layer.GetSpatialRef()

    # Load all blocks
    blocks = []
    for feature in grid_layer:
        geom = feature.GetGeometryRef()
        if geom is None:
            continue

        geom = geom.Clone()
        if not geom.IsValid():
            geom = geom.Buffer(0)

        setback_id = (
            feature.GetField("setback_id") if feature.GetFieldIndex("setback_id") >= 0 else None
        )
        area = geom.Area()

        blocks.append(
            {
                "geom": geom,
                "setback_id": setback_id,
                "area": area,
                "merged": False,
                "removed": False,
            }
        )

    input_ds = None

    # Process small blocks
    max_distance_threshold = 10  # Maximum distance to consider as "close neighbor"

    for i, block in enumerate(blocks):
        if block["removed"] or block["merged"]:
            continue

        # Check if block is too small
        if block["area"] >= min_area:
            continue

        # Find nearest neighbor with shared boundary or very close neighbor
        best_neighbor_idx = None
        best_shared_length = 0
        best_distance = float("inf")

        for j, other in enumerate(blocks):
            if i == j or other["removed"] or other["merged"]:
                continue

            # Calculate distance between blocks
            distance = block["geom"].Distance(other["geom"])

            # Only consider blocks within threshold distance
            if distance > max_distance_threshold:
                continue

            # Use a buffer to detect shared boundaries
            # Buffer size should be slightly larger than expected gap
            buffer_size = 0.5
            buffered_block = block["geom"].Buffer(buffer_size)
            intersection = buffered_block.Intersection(other["geom"])

            shared_length = 0
            if not intersection.IsEmpty():
                geom_name = intersection.GetGeometryName()
                if geom_name in ("LINESTRING", "MULTILINESTRING"):
                    shared_length = intersection.Length()
                elif geom_name in ("POLYGON", "MULTIPOLYGON"):
                    # If intersection is a polygon, they overlap
                    # Use the perimeter of the intersection as shared boundary
                    boundary = intersection.Boundary()
                    if boundary is not None:
                        if boundary.GetGeometryName() == "LINESTRING":
                            shared_length = boundary.Length()
                        elif boundary.GetGeometryName() == "MULTILINESTRING":
                            for k in range(boundary.GetGeometryCount()):
                                shared_length += boundary.GetGeometryRef(k).Length()

            # Prioritize blocks with shared boundaries
            if shared_length > 0:
                if shared_length > best_shared_length:
                    best_shared_length = shared_length
                    best_neighbor_idx = j
                    best_distance = distance
            else:
                # If no shared boundary found, track nearest by distance
                if best_shared_length == 0 and distance < best_distance:
                    best_distance = distance
                    best_neighbor_idx = j

        if best_neighbor_idx is not None:
            merged_geom = block["geom"].Union(blocks[best_neighbor_idx]["geom"])

            if not merged_geom.IsValid():
                merged_geom = merged_geom.Buffer(0)

            # Use Shapely to dissolve internal boundaries and convert to single polygon
            try:
                # Convert OGR geometry to Shapely
                shapely_geom = wkt.loads(merged_geom.ExportToIsoWkt())

                # If it's a MultiPolygon or GeometryCollection, use unary_union to dissolve
                if isinstance(shapely_geom, (geometry.MultiPolygon, geometry.GeometryCollection)):
                    # Unary union will merge overlapping/touching polygons into single polygons
                    dissolved = unary_union(shapely_geom)

                    # If still multipolygon after union, buffer slightly to merge touching polygons
                    if isinstance(dissolved, geometry.MultiPolygon):
                        dissolved = dissolved.buffer(1, join_style=2)
                        dissolved = unary_union(dissolved)
                        dissolved = dissolved.buffer(-1, join_style=2)

                    shapely_geom = dissolved

                # Final buffer out and in to clean up any remaining internal edges
                shapely_geom = shapely_geom.buffer(0.001, join_style=2)  # join_style=2 is mitre
                shapely_geom = shapely_geom.buffer(-0.001, join_style=2)

                # Convert back to OGR geometry
                merged_geom = ogr.CreateGeometryFromWkb(shapely_geom.wkb)

                if not merged_geom.IsValid():
                    merged_geom = merged_geom.Buffer(0)

            except Exception as e:
                print(f"Warning: Shapely conversion failed, using OGR buffer method: {e}")
                # Fallback to original buffer method
                merged_geom = merged_geom.Buffer(0.001)
                merged_geom = merged_geom.Buffer(-0.001)

                if not merged_geom.IsValid():
                    merged_geom = merged_geom.Buffer(0)

                # If still multipolygon, take largest
                if merged_geom.GetGeometryName() == "MULTIPOLYGON":
                    largest_area = 0
                    largest_poly = None
                    for k in range(merged_geom.GetGeometryCount()):
                        poly = merged_geom.GetGeometryRef(k)
                        if poly.Area() > largest_area:
                            largest_area = poly.Area()
                            largest_poly = poly.Clone()
                    if largest_poly is not None:
                        merged_geom = largest_poly

            blocks[best_neighbor_idx]["geom"] = merged_geom
            blocks[best_neighbor_idx]["area"] = merged_geom.Area()

            block["merged"] = True

            if best_shared_length > 0:
                print(
                    f"Merged small block {i} (area: {block['area']:.2f}) "
                    f"with adjacent neighbor {best_neighbor_idx} "
                    f"(shared edge: {best_shared_length:.2f}m)"
                )
            else:
                print(
                    f"Merged small block {i} (area: {block['area']:.2f}) "
                    f"with close neighbor {best_neighbor_idx} (distance: {best_distance:.2f}m)"
                )
        else:
            block["removed"] = True
            print(
                f"Removed isolated block {i} (area: {block['area']:.2f}) - "
                f"no adjacent neighbors found"
            )

    # Collect remaining blocks
    cleaned_blocks = []
    for block in blocks:
        if not block["removed"] and not block["merged"]:
            cleaned_blocks.append((block["geom"], block["setback_id"]))

    # Write output
    driver = ogr.GetDriverByName("GPKG")
    if driver is None:
        raise RuntimeError("Could not get GPKG driver.")

    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
        if output_ds is None:
            raise RuntimeError(f"Could not open existing output: {output_path}")
    else:
        output_ds = driver.CreateDataSource(output_path)
        if output_ds is None:
            raise RuntimeError(f"Could not create output: {output_path}")

    for i in range(output_ds.GetLayerCount()):
        layer = output_ds.GetLayerByIndex(i)
        if layer.GetName() == output_layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)
    if output_layer is None:
        raise RuntimeError("Could not create output layer.")

    id_field = ogr.FieldDefn("id", ogr.OFTInteger)
    output_layer.CreateField(id_field)

    setback_id_field = ogr.FieldDefn("setback_id", ogr.OFTInteger)
    output_layer.CreateField(setback_id_field)

    area_field = ogr.FieldDefn("area", ogr.OFTReal)
    output_layer.CreateField(area_field)

    feature_id = 1
    for polygon_geom, setback_id in cleaned_blocks:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(polygon_geom)
        out_feature.SetField("id", feature_id)
        if setback_id is not None:
            out_feature.SetField("setback_id", int(setback_id))
        out_feature.SetField("area", float(polygon_geom.Area()))

        output_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None

    print(
        f"Cleaned grid: {len(cleaned_blocks)} blocks remaining from {len(blocks)} original blocks"
    )


def cleanup_intermediate_layers(gpkg_path, layers_to_keep):
    """Remove all layers from a GeoPackage except the ones specified.

    Args:
        gpkg_path (str): Path to the GeoPackage to modify.
        layers_to_keep (Iterable[str]): Names of layers that must be preserved.

    Returns:
        None

    Raises:
        Exception: Propagated GDAL/OGR errors if the GeoPackage cannot be opened
            or layers cannot be deleted.

    Notes:
        * Deletion is done by name matching; missing names are ignored.
        * Enumerates names first to avoid index-shift issues during deletion.
    """
    ds = ogr.Open(gpkg_path, 1)

    # Get all layer names
    layer_names = []
    for i in range(ds.GetLayerCount()):
        layer_names.append(ds.GetLayerByIndex(i).GetName())

    # Delete layers not in keep list
    for layer_name in layer_names:
        if layer_name not in layers_to_keep:
            for i in range(ds.GetLayerCount()):
                if ds.GetLayerByIndex(i).GetName() == layer_name:
                    ds.DeleteLayer(i)
                    break

    ds = None


def merge_grid_layers_with_type(
    input_path,
    output_path,
    output_layer_name,
    layer_configs,
):
    """Merge multiple grid layers into one, adding a grid_type field to each.

    Args:
        input_path (str): Path to input GeoPackage containing source layers.
        output_path (str): Path to output GeoPackage.
        output_layer_name (str): Name for the merged output layer.
        layer_configs (list): List of tuples (layer_name, grid_type_value)
            Example: [
                ("intersected_setback", "on_grid_intersected"),
                ("arterial_setback_grid_cleaned", "on_grid_art"),
                ("secondary_setback_grid_cleaned", "on_grid_sec"),
                ("site_minus_all_setbacks_grid_cells", "off_grid"),
            ]

    Returns:
        str: Name of the created output layer.

    Raises:
        RuntimeError: If input GeoPackage cannot be opened or layers not found.
    """
    if input_path == output_path:
        ds = ogr.Open(input_path, 1)
        if ds is None:
            raise RuntimeError(f"Cannot open GeoPackage: {input_path}")

        first_layer_name = layer_configs[0][0]
        first_layer = ds.GetLayerByName(first_layer_name)
        if first_layer is None:
            ds = None
            raise RuntimeError(f"Cannot find layer: {first_layer_name}")

        srs = first_layer.GetSpatialRef()

        # Collect all unique field definitions from all source layers
        all_field_defs = {}  # field_name -> field_defn

        for layer_name, _ in layer_configs:
            source_layer = ds.GetLayerByName(layer_name)
            if source_layer is None:
                continue

            layer_defn = source_layer.GetLayerDefn()
            for i in range(layer_defn.GetFieldCount()):
                field_defn = layer_defn.GetFieldDefn(i)
                field_name = field_defn.GetName()
                # Keep the first occurrence of each field name
                if field_name not in all_field_defs:
                    all_field_defs[field_name] = field_defn

        all_features_data = []

        for layer_name, grid_type_value in layer_configs:
            source_layer = ds.GetLayerByName(layer_name)
            if source_layer is None:
                print(f"Warning: Layer '{layer_name}' not found, skipping...")
                continue

            source_layer.ResetReading()
            layer_features = []
            for source_feature in source_layer:
                geom = source_feature.GetGeometryRef()
                geom_wkt = geom.ExportToWkt() if geom else None

                field_values = {}
                for i in range(source_feature.GetFieldCount()):
                    field_name = source_feature.GetFieldDefnRef(i).GetName()
                    if source_feature.IsFieldSet(i):
                        field_values[field_name] = source_feature.GetField(i)

                layer_features.append((geom_wkt, field_values, grid_type_value))

            all_features_data.extend(layer_features)
            print(f"  Collected {len(layer_features)} features from '{layer_name}'")

        if ds.GetLayerByName(output_layer_name):
            ds.DeleteLayer(output_layer_name)

        output_layer = ds.CreateLayer(output_layer_name, srs, geom_type=ogr.wkbPolygon)

        # Create all fields in output layer
        for _field_name, field_defn in all_field_defs.items():
            output_layer.CreateField(field_defn)

        grid_type_field = ogr.FieldDefn("type", ogr.OFTString)
        grid_type_field.SetWidth(50)
        output_layer.CreateField(grid_type_field)

        feature_count = 0
        for geom_wkt, field_values, grid_type_value in all_features_data:
            out_feature = ogr.Feature(output_layer.GetLayerDefn())

            if geom_wkt:
                geom = ogr.CreateGeometryFromWkt(geom_wkt)
                out_feature.SetGeometry(geom)

            for field_name, field_value in field_values.items():
                try:
                    out_feature.SetField(field_name, field_value)
                except RuntimeError:
                    pass

            out_feature.SetField("type", grid_type_value)

            output_layer.CreateFeature(out_feature)
            out_feature = None
            feature_count += 1

        if ds is not None:
            ds.FlushCache()
        ds = None

    else:
        input_ds = ogr.Open(input_path, 0)
        if input_ds is None:
            raise RuntimeError(f"Cannot open input GeoPackage: {input_path}")

        output_ds = ogr.Open(output_path, 1)
        if output_ds is None:
            input_ds = None
            raise RuntimeError(f"Cannot open output GeoPackage: {output_path}")

        first_layer_name = layer_configs[0][0]
        first_layer = input_ds.GetLayerByName(first_layer_name)
        if first_layer is None:
            input_ds = None
            output_ds = None
            raise RuntimeError(f"Cannot find layer: {first_layer_name}")

        srs = first_layer.GetSpatialRef()

        if output_ds.GetLayerByName(output_layer_name):
            output_ds.DeleteLayer(output_layer_name)

        output_layer = output_ds.CreateLayer(output_layer_name, srs, geom_type=ogr.wkbPolygon)

        all_field_defs = {}

        for layer_name, _ in layer_configs:
            source_layer = input_ds.GetLayerByName(layer_name)
            if source_layer is None:
                continue

            layer_defn = source_layer.GetLayerDefn()
            for i in range(layer_defn.GetFieldCount()):
                field_defn = layer_defn.GetFieldDefn(i)
                field_name = field_defn.GetName()
                if field_name not in all_field_defs:
                    all_field_defs[field_name] = field_defn

        for _field_name, field_defn in all_field_defs.items():
            output_layer.CreateField(field_defn)

        grid_type_field = ogr.FieldDefn("type", ogr.OFTString)
        grid_type_field.SetWidth(50)
        output_layer.CreateField(grid_type_field)

        feature_count = 0
        for layer_name, grid_type_value in layer_configs:
            source_layer = input_ds.GetLayerByName(layer_name)
            if source_layer is None:
                print(f"Warning: Layer '{layer_name}' not found, skipping...")
                continue

            source_layer.ResetReading()
            for source_feature in source_layer:
                out_feature = ogr.Feature(output_layer.GetLayerDefn())

                geom = source_feature.GetGeometryRef()
                if geom:
                    out_feature.SetGeometry(geom.Clone())

                # Copy all fields that exist in source feature
                for i in range(source_feature.GetFieldCount()):
                    field_name = source_feature.GetFieldDefnRef(i).GetName()
                    if source_feature.IsFieldSet(i):
                        try:
                            out_feature.SetField(field_name, source_feature.GetField(i))
                        except RuntimeError:
                            # Field doesn't exist in output schema, skip
                            pass

                out_feature.SetField("type", grid_type_value)

                output_layer.CreateFeature(out_feature)
                out_feature = None
                feature_count += 1

            print(f"  Merged {source_layer.GetFeatureCount()} features from '{layer_name}'")

        if output_ds is not None:
            output_ds.FlushCache()
        input_ds = None
        output_ds = None

    print(f"Merged {feature_count} total features into '{output_layer_name}'")

    return output_layer_name


def export_layer_to_geojson_gpd(
    input_path,
    layer_name,
    output_geojson_path,
):
    """Export a GeoPackage layer to GeoJSON format in EPSG:4326, filtering to only valid
        polygons and lines.

    Args:
        input_path (str): Path to input GeoPackage.
        layer_name (str): Name of layer to export.
        output_geojson_path (str): Path for output GeoJSON file.

    Returns:
        str: Path to the created GeoJSON file.

    Raises:
        RuntimeError: If input GeoPackage cannot be opened or layer not found.
    """
    import os

    import geopandas as gpd
    from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon

    try:
        gdf = gpd.read_file(input_path, layer=layer_name)
    except Exception as e:
        print(f"Cannot read layer '{layer_name}' from {input_path}: {e}")
        return None

    if gdf.empty:
        print(f"Warning: Layer '{layer_name}' is empty")
        if gdf.crs is not None and gdf.crs != "EPSG:4326":
            gdf = gdf.to_crs("EPSG:4326")
        gdf.to_file(output_geojson_path, driver="GeoJSON")
        return output_geojson_path

    valid_geom_types = (Polygon, MultiPolygon, LineString, MultiLineString)

    def is_valid_geometry(geom):
        if geom is None or geom.is_empty:
            return False
        return isinstance(geom, valid_geom_types)

    gdf_filtered = gdf[gdf.geometry.apply(is_valid_geometry)].copy()

    invalid_mask = ~gdf_filtered.geometry.is_valid
    if invalid_mask.any():
        print(f"  Fixing {invalid_mask.sum()} invalid geometries")
        gdf_filtered.loc[invalid_mask, "geometry"] = gdf_filtered.loc[
            invalid_mask, "geometry"
        ].buffer(0)

    skipped_count = len(gdf) - len(gdf_filtered)
    if skipped_count > 0:
        print(f"  Skipped {skipped_count} features with invalid geometry types")

    if gdf_filtered.empty:
        print(f"Warning: No valid polygon or line geometries found in layer '{layer_name}'")
        if gdf_filtered.crs is not None and gdf_filtered.crs != "EPSG:4326":
            gdf_filtered = gdf_filtered.to_crs("EPSG:4326")
        gdf_filtered.to_file(output_geojson_path, driver="GeoJSON")
        return output_geojson_path

    if gdf_filtered.crs is not None and gdf_filtered.crs != "EPSG:4326":
        print(f"  Transforming from {gdf_filtered.crs} to EPSG:4326")
        gdf_filtered = gdf_filtered.to_crs("EPSG:4326")

    if os.path.exists(output_geojson_path):
        os.remove(output_geojson_path)

    gdf_filtered.to_file(output_geojson_path, driver="GeoJSON")
    print(f"Exported {len(gdf_filtered)} features to '{output_geojson_path}'")

    return output_geojson_path


def export_layer_to_geojson(
    input_path,
    layer_name,
    output_geojson_path,
):
    """Export a GeoPackage layer to GeoJSON format.

    Args:
        input_path (str): Path to input GeoPackage.
        layer_name (str): Name of layer to export.
        output_geojson_path (str): Path for output GeoJSON file.

    Returns:
        str: Path to the created GeoJSON file.

    Raises:
        RuntimeError: If input GeoPackage cannot be opened or layer not found.
    """
    input_ds = ogr.Open(input_path, 0)
    if input_ds is None:
        raise RuntimeError(f"Cannot open input GeoPackage: {input_path}")

    layer = input_ds.GetLayerByName(layer_name)
    if layer is None:
        input_ds = None
        raise RuntimeError(f"Cannot find layer: {layer_name}")

    driver = ogr.GetDriverByName("GeoJSON")
    if driver is None:
        input_ds = None
        raise RuntimeError("GeoJSON driver not available")

    import os

    if os.path.exists(output_geojson_path):
        os.remove(output_geojson_path)

    output_ds = driver.CreateDataSource(output_geojson_path)
    if output_ds is None:
        input_ds = None
        raise RuntimeError(f"Cannot create GeoJSON file: {output_geojson_path}")

    output_layer = output_ds.CopyLayer(layer, layer_name, ["RFC7946=YES"])

    if output_layer is None:
        input_ds = None
        output_ds = None
        raise RuntimeError("Failed to copy layer to GeoJSON")

    feature_count = output_layer.GetFeatureCount()

    input_ds = None
    output_ds = None

    print(f"Exported {feature_count} features to '{output_geojson_path}'")

    return output_geojson_path
