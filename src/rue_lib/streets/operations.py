# src/rue_lib/streets/operations.py
import math
import os

from osgeo import ogr
from shapely import geometry, wkt
from shapely.ops import unary_union

from .geometry_utils import (
    extract_edges_from_geom,
    merge_connected_edges,
)


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


def calculate_required_rings(gpkg_path, layer_name, ring_spacing):
    """Calculate the number of rings needed to cover a layer's extent.

    Args:
        gpkg_path (str): Path to the GeoPackage
        layer_name (str): Name of the layer to analyze
        ring_spacing (float): Distance between rings

    Returns:
        int: Number of rings needed to cover the extent
    """
    ds = ogr.Open(gpkg_path)
    layer = ds.GetLayerByName(layer_name)

    extent = layer.GetExtent()
    ds = None

    width = extent[1] - extent[0]
    height = extent[3] - extent[2]
    diagonal = (width**2 + height**2) ** 0.5

    rings = int(diagonal / ring_spacing) + 1

    return rings


def multiring_buffer(input_path, layer_name, rings, distance, output_path, output_layer_name):
    """Create concentric ring buffers around the union of input features.

    Args:
        input_path (str): Path to the input dataset (e.g., .gpkg, .geojson).
        layer_name (str): Name of the layer within `input_path` to buffer.
        rings (int): Number of concentric rings to create.
        distance (float): Step width between rings (in layer units).
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output polygon layer to create.

    Returns:
        None

    Raises:
        Exception: Propagated GDAL/OGR errors when reading, buffering, or writing.

    Notes:
        * Source geometries are dissolved before ring generation.
        * The output spatial reference matches the input layer.
        * Ring polygons are non-overlapping annuli.
    """
    source_ds = ogr.Open(input_path)
    source_layer = source_ds.GetLayerByName(layer_name)

    srs = source_layer.GetSpatialRef()

    source_geoms = []
    for feature in source_layer:
        source_geoms.append(feature.GetGeometryRef().Clone())

    source_ds = None

    base_geom = source_geoms[0]
    for geom in source_geoms[1:]:
        base_geom = base_geom.Union(geom)

    # Validate and fix the base geometry before creating rings
    if not base_geom.IsValid():
        print("Warning: Base geometry is invalid, attempting to fix...")
        base_geom = base_geom.Buffer(0)

    ring_geoms = []
    for ring_num in range(1, rings + 1):
        outer_distance = ring_num * distance
        inner_distance = (ring_num - 1) * distance

        outer_buffer = base_geom.Buffer(outer_distance)
        inner_buffer = base_geom.Buffer(inner_distance)

        # Apply buffer(0) to fix any self-intersections or invalid geometries
        if not outer_buffer.IsValid():
            outer_buffer = outer_buffer.Buffer(0)
        if not inner_buffer.IsValid():
            inner_buffer = inner_buffer.Buffer(0)

        ring_geom = outer_buffer.Difference(inner_buffer)

        # Validate and fix the resulting ring geometry
        if not ring_geom.IsValid():
            ring_geom = ring_geom.Buffer(0)

        if not ring_geom.IsEmpty():
            ring_geoms.append((ring_geom, ring_num))

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

    field_defn = ogr.FieldDefn("ring", ogr.OFTInteger)
    output_layer.CreateField(field_defn)

    for ring_geom, ring_num in ring_geoms:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(ring_geom)
        out_feature.SetField("ring", ring_num)
        output_layer.CreateFeature(out_feature)
        out_feature = None

    output_ds = None


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
    arterial_buffer_distance=15,
):
    """Create grid polygons from arterial setback by creating division points and splitting polygons

    This function combines the creation of division points and the splitting of polygons into
    a single operation. It generates division points along edges closest to arterial roads,
    then splits the setback polygons using perpendicular lines from these points.

    Args:
        input_path (str): Path to the dataset containing the setback & arterial layers.
        setback_layer_name (str): Name of the arterial setback layer (polygons).
        arterial_road_layer_name (str): Name of the arterial road layer (lines).
        grid_width (float): Width for division points spacing.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Name of the output layer to create.
        arterial_buffer_distance (float | None): Optional buffer distance (same
            units as layer CRS). If provided and > 0, only edges that intersect
            the buffered arterial roads are considered.

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
    if arterial_buffer_distance is not None and arterial_buffer_distance > 0:
        arterial_buffer_geom = arterial_union.Buffer(arterial_buffer_distance)

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

            points_by_setback[setback_id].append((point_x, point_y, current_distance))

            current_distance += grid_width

        # Always add the end point if not already included
        if current_distance - grid_width < total_length:
            points_by_setback[setback_id].append(
                (edge_geom.GetX(n - 1), edge_geom.GetY(n - 1), total_length)
            )

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
