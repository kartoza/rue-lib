# src/rue_lib/streets/runner.py

import os
from dataclasses import dataclass
from pathlib import Path

from osgeo import gdal, ogr, osr

# Enable GDAL exceptions
gdal.UseExceptions()


@dataclass
class StreetConfig:
    """Configuration for street generation."""

    parcel_path: str  # Output generated from step 1
    roads_path: str
    on_grid_partition_depth_arterial_roads: float = 40.0
    on_grid_partition_depth_secondary_roads: float = 30.0
    off_grid_partitions_preferred_depth: float = 44.0
    off_grid_partitions_preferred_width: float = 30.0
    arterial_setback_depth: float = 60.0  # Depth of arterial road setback zone
    secondary_setback_depth: float = 60.0  # Depth of secondary road setback zone
    perpendicular_line_length: float = 500.0  # Length of perpendicular lines
    output_dir: str = "outputs/streets"
    road_local_width_m: float = 12.0  # ROAD_LOC_W_


def get_utm_zone_from_layer(layer):
    """Determine appropriate UTM zone from layer extent."""
    extent = layer.GetExtent()
    lon_center = (extent[0] + extent[1]) / 2
    lat_center = (extent[2] + extent[3]) / 2

    # Calculate UTM zone
    zone = int((lon_center + 180) / 6) + 1

    # Determine if northern or southern hemisphere
    if lat_center >= 0:
        epsg_code = 32600 + zone  # Northern hemisphere
    else:
        epsg_code = 32700 + zone  # Southern hemisphere

    return epsg_code


def reproject_layer(input_path, output_path, target_epsg):
    """Reproject a layer to target CRS."""
    source_ds = ogr.Open(input_path)
    source_layer = source_ds.GetLayer()

    source_srs = source_layer.GetSpatialRef()
    target_srs = osr.SpatialReference()
    target_srs.ImportFromEPSG(target_epsg)

    transform = osr.CoordinateTransformation(source_srs, target_srs)

    layer_name = os.path.splitext(os.path.basename(input_path))[0] + f"_{target_epsg}"

    # Create output
    driver = ogr.GetDriverByName("GPKG")
    if os.path.exists(output_path):
        output_ds = driver.Open(output_path, 1)
    else:
        output_ds = driver.CreateDataSource(output_path)

    # Remove layer if exists
    for i in range(output_ds.GetLayerCount()):
        if output_ds.GetLayerByIndex(i).GetName() == layer_name:
            output_ds.DeleteLayer(i)
            break

    output_layer = output_ds.CreateLayer(layer_name, target_srs, source_layer.GetGeomType())

    # Copy field definitions
    source_layer_defn = source_layer.GetLayerDefn()
    for i in range(source_layer_defn.GetFieldCount()):
        field_defn = source_layer_defn.GetFieldDefn(i)
        output_layer.CreateField(field_defn)

    # Transform and copy features
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

    # Explicitly close datasets
    source_ds = None
    output_ds = None

    return layer_name


def extract_by_expression(input_path, layer_name, expression, output_path, output_layer_name):
    """Extract features matching an expression."""
    source_ds = ogr.Open(input_path)
    source_layer = source_ds.GetLayerByName(layer_name)

    # Copy field definitions and geometry type before closing
    source_layer_defn = source_layer.GetLayerDefn()
    srs = source_layer.GetSpatialRef()
    geom_type = source_layer.GetGeomType()

    # Get fields info
    field_defs = []
    for i in range(source_layer_defn.GetFieldCount()):
        field_defs.append(source_layer_defn.GetFieldDefn(i))

    # Apply filter and get features
    source_layer.SetAttributeFilter(expression)
    features_data = []
    for feature in source_layer:
        geom = feature.GetGeometryRef().Clone()
        field_values = {}
        for i in range(source_layer_defn.GetFieldCount()):
            field_name = source_layer_defn.GetFieldDefn(i).GetNameRef()
            field_values[field_name] = feature.GetField(i)
        features_data.append((geom, field_values))

    # Close source dataset
    source_ds = None

    # Now open output dataset
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

    # Create fields
    for field_def in field_defs:
        output_layer.CreateField(field_def)

    # Write features
    for geom, field_values in features_data:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)

        for field_name, value in field_values.items():
            out_feature.SetField(field_name, value)

        output_layer.CreateFeature(out_feature)
        out_feature = None

    output_ds = None


def buffer_layer(input_path, layer_name, distance, output_path, output_layer_name, dissolve=True):
    """Buffer a layer."""
    source_ds = ogr.Open(input_path)
    source_layer = source_ds.GetLayerByName(layer_name)

    srs = source_layer.GetSpatialRef()

    # Collect all geometries
    geoms = []
    for feature in source_layer:
        geom = feature.GetGeometryRef().Clone()
        buffered = geom.Buffer(distance)
        geoms.append(buffered)

    # Close source
    source_ds = None

    # Process geometries
    if dissolve and geoms:
        union_geom = geoms[0]
        for geom in geoms[1:]:
            union_geom = union_geom.Union(geom)
        geoms = [union_geom]

    # Open output
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

    output_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)

    for geom in geoms:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)
        output_layer.CreateFeature(out_feature)
        out_feature = None

    output_ds = None


def clip_layer(
    input_path, input_layer_name, clip_path, clip_layer_name, output_path, output_layer_name
):
    """Clip a layer by another layer."""
    # Open input
    input_ds = ogr.Open(input_path)
    input_layer = input_ds.GetLayerByName(input_layer_name)

    srs = input_layer.GetSpatialRef()
    geom_type = input_layer.GetGeomType()

    # Copy field definitions
    input_layer_defn = input_layer.GetLayerDefn()
    field_defs = []
    for i in range(input_layer_defn.GetFieldCount()):
        field_defs.append(input_layer_defn.GetFieldDefn(i))

    # Get input features
    input_features = []
    for feature in input_layer:
        geom = feature.GetGeometryRef().Clone()
        field_values = {}
        for i in range(input_layer_defn.GetFieldCount()):
            field_name = input_layer_defn.GetFieldDefn(i).GetNameRef()
            field_values[field_name] = feature.GetField(i)
        input_features.append((geom, field_values))

    input_ds = None

    # Open clip layer
    if clip_path == input_path:
        # Same file, need to be careful
        clip_ds = ogr.Open(clip_path)
    else:
        clip_ds = ogr.Open(clip_path)

    clip_layer = clip_ds.GetLayerByName(clip_layer_name)

    # Get clip geometry
    clip_geoms = []
    for feature in clip_layer:
        clip_geoms.append(feature.GetGeometryRef().Clone())

    clip_ds = None

    clip_geom = clip_geoms[0]
    for geom in clip_geoms[1:]:
        clip_geom = clip_geom.Union(geom)

    # Clip features
    clipped_features = []
    for geom, field_values in input_features:
        clipped_geom = geom.Intersection(clip_geom)
        if not clipped_geom.IsEmpty():
            clipped_features.append((clipped_geom, field_values))

    # Open output
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

    # Create fields
    for field_def in field_defs:
        output_layer.CreateField(field_def)

    # Write clipped features
    for geom, field_values in clipped_features:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)

        for field_name, value in field_values.items():
            out_feature.SetField(field_name, value)

        output_layer.CreateFeature(out_feature)
        out_feature = None

    output_ds = None


def multiring_buffer(input_path, layer_name, rings, distance, output_path, output_layer_name):
    """Create concentric ring buffers."""
    source_ds = ogr.Open(input_path)
    source_layer = source_ds.GetLayerByName(layer_name)

    srs = source_layer.GetSpatialRef()

    # Collect all source geometries
    source_geoms = []
    for feature in source_layer:
        source_geoms.append(feature.GetGeometryRef().Clone())

    source_ds = None

    # Union source geometries
    base_geom = source_geoms[0]
    for geom in source_geoms[1:]:
        base_geom = base_geom.Union(geom)

    # Create rings
    ring_geoms = []
    for ring_num in range(1, rings + 1):
        outer_distance = ring_num * distance
        inner_distance = (ring_num - 1) * distance

        outer_buffer = base_geom.Buffer(outer_distance)
        inner_buffer = base_geom.Buffer(inner_distance)

        ring_geom = outer_buffer.Difference(inner_buffer)
        ring_geoms.append((ring_geom, ring_num))

    # Open output
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

    output_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)

    # Add ring number field
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
    """Remove all layers except the ones specified."""
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


def generate_streets(cfg: StreetConfig) -> Path:
    """
    Generate street blocks from roads and parcels

    Returns:
        Path to output blocks file
    """
    # Create output directory
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_gpkg = output_dir / "outputs.gpkg"
    output_path = str(output_gpkg)

    print("Step 1: Determining UTM zone...")
    # Get UTM zone from site layer
    site_ds = ogr.Open(cfg.parcel_path)
    site_layer = site_ds.GetLayer()
    utm_epsg = get_utm_zone_from_layer(site_layer)
    site_ds = None
    print(f"Using UTM EPSG: {utm_epsg}")

    print("Step 2: Reprojecting layers to UTM...")
    site_layer_name = reproject_layer(cfg.parcel_path, output_path, utm_epsg)
    roads_layer_name = reproject_layer(cfg.roads_path, output_path, utm_epsg)

    print("Step 3: Extracting arterial roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_art'", output_path, "arterial_roads"
    )

    print("Step 4: Extracting secondary roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_sec'", output_path, "secondary_roads"
    )

    print("Step 5: Creating arterial road setback zone...")
    # Use the larger setback depth for the arterial setback
    buffer_layer(
        output_path,
        "arterial_roads",
        cfg.on_grid_partition_depth_arterial_roads,
        output_path,
        "arterial_buffered_large",
        dissolve=True,
    )

    print("Step 6: Clipping arterial setback by site...")
    clip_layer(
        output_path,
        "arterial_buffered_large",
        output_path,
        site_layer_name,
        output_path,
        "arterial_road_setback",
    )

    print("Step 7: Creating street block rings...")

    multiring_buffer(
        output_path,
        "arterial_buffered_large",
        10,  # Number of rings
        cfg.off_grid_partitions_preferred_depth,
        output_path,
        "arterial_offset_buffer_rings",
    )

    print("Step 8: Clipping block rings by site...")
    clip_layer(
        output_path,
        "arterial_offset_buffer_rings",
        output_path,
        site_layer_name,
        output_path,
        "street_blocks",
    )

    print("Step 9: Creating secondary road setback zone...")
    # Use the larger setback depth for secondary roads
    buffer_layer(
        output_path,
        "secondary_roads",
        cfg.secondary_setback_depth,
        output_path,
        "secondary_roads_buffered_large",
        dissolve=True,
    )

    print("Step 10: Clipping secondary setback by site...")
    clip_layer(
        output_path,
        "secondary_roads_buffered_large",
        output_path,
        site_layer_name,
        output_path,
        "secondary_road_setback",
    )

    # Clean up intermediate layers
    print("Cleaning up intermediate layers...")
    final_layers = [
        "arterial_road_setback",
        "street_blocks",
        "secondary_road_setback",
    ]
    cleanup_intermediate_layers(output_path, final_layers)

    print(f"\nProcessing complete! Output saved to: {output_gpkg}")
    print("\nFinal layers:")
    print(
        f"  - arterial_road_setback: {cfg.arterial_setback_depth}m "
        f"buffer zone around arterial roads"
    )
    print(
        f"  - street_blocks: Concentric block rings ({cfg.off_grid_partitions_preferred_depth}m "
        f"spacing)"
    )
    print(
        f"  - secondary_road_setback: {cfg.secondary_setback_depth}m buffer "
        f"zone around secondary roads"
    )

    return output_gpkg
