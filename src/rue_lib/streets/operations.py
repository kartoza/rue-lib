# src/rue_lib/streets/operations.py
import os

from osgeo import ogr


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

    # Perform difference operation (subtract erase geometry from input)
    erased_features = []
    for geom, field_values in input_features:
        # Validate input geometry before difference
        if not geom.IsValid():
            geom = geom.Buffer(0)

        erased_geom = geom.Difference(erase_geom)

        # Validate resulting geometry
        if not erased_geom.IsEmpty():
            if not erased_geom.IsValid():
                erased_geom = erased_geom.Buffer(0)

            if not erased_geom.IsEmpty():
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

        # Set id field if we created it, or update it if it already existed
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
