# src/rue_lib/streets/runner.py

import os
from dataclasses import dataclass
from pathlib import Path

from osgeo import gdal, ogr

from rue_lib.core.geometry import buffer_layer, get_utm_zone_from_layer, reproject_layer

# Enable GDAL exceptions
gdal.UseExceptions()


@dataclass
class StreetConfig:
    """Configuration for street generation."""

    parcel_path: str  # Output generated from step 1
    roads_path: str
    on_grid_partition_depth_arterial_roads: float = 40.0
    on_grid_partition_depth_secondary_roads: float = 30.0
    off_grid_partitions_preferred_depth: float = 80.0
    off_grid_partitions_preferred_width: float = 90.0
    arterial_setback_depth: float = 60.0  # Depth of arterial road setback zone
    secondary_setback_depth: float = 60.0  # Depth of secondary road setback zone
    perpendicular_line_length: float = 1000.0  # Length of perpendicular lines
    output_dir: str = "outputs/streets"
    road_local_width_m: float = 12.0  # ROAD_LOC_W_


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

    clipped_features = []
    for geom, field_values in input_features:
        clipped_geom = geom.Intersection(clip_geom)
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

    for field_def in field_defs:
        output_layer.CreateField(field_def)

    for geom, field_values in clipped_features:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)

        for field_name, value in field_values.items():
            out_feature.SetField(field_name, value)

        output_layer.CreateFeature(out_feature)
        out_feature = None

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

    ring_geoms = []
    for ring_num in range(1, rings + 1):
        outer_distance = ring_num * distance
        inner_distance = (ring_num - 1) * distance

        outer_buffer = base_geom.Buffer(outer_distance)
        inner_buffer = base_geom.Buffer(inner_distance)

        ring_geom = outer_buffer.Difference(inner_buffer)
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


def extract_arterial_edge_lines(
    gpkg_path,
    arterial_roads_layer_name,
    arterial_setback_layer_name,
    secondary_setback_layer_name,
    output_layer_name,
    clip_buffer=0.1,
    sample_distance=5.0,
):
    ds = ogr.Open(gpkg_path, 1)

    arterial_layer = ds.GetLayerByName(arterial_roads_layer_name)
    arterial_geoms = []
    for feature in arterial_layer:
        arterial_geoms.append(feature.GetGeometryRef().Clone())

    arterial_union = arterial_geoms[0]
    for geom in arterial_geoms[1:]:
        arterial_union = arterial_union.Union(geom)

    arterial_setback_layer = ds.GetLayerByName(arterial_setback_layer_name)
    arterial_setback_geoms = []
    for feature in arterial_setback_layer:
        arterial_setback_geoms.append(feature.GetGeometryRef().Clone())

    arterial_setback_union = arterial_setback_geoms[0]
    for geom in arterial_setback_geoms[1:]:
        arterial_setback_union = arterial_setback_union.Union(geom)

    secondary_setback_layer = ds.GetLayerByName(secondary_setback_layer_name)
    secondary_setback_geoms = []
    for feature in secondary_setback_layer:
        secondary_setback_geoms.append(feature.GetGeometryRef().Clone())

    secondary_setback_union = secondary_setback_geoms[0]
    for geom in secondary_setback_geoms[1:]:
        secondary_setback_union = secondary_setback_union.Union(geom)

    srs = arterial_setback_layer.GetSpatialRef()

    secondary_setback_buffered = secondary_setback_union.Buffer(clip_buffer)
    clipped_geom = arterial_setback_union.Difference(secondary_setback_buffered)

    if clipped_geom.IsEmpty():
        print("Warning: Clipped geometry is empty!")
        ds = None
        return

    def process_polygon(poly):
        exterior_ring = poly.GetGeometryRef(0)

        points = []
        for i in range(exterior_ring.GetPointCount()):
            pt = ogr.Geometry(ogr.wkbPoint)
            pt.AddPoint(exterior_ring.GetX(i), exterior_ring.GetY(i))
            points.append(pt)

        point_distances = [(pt, pt.Distance(arterial_union)) for pt in points]
        min_dist = min(dist for _, dist in point_distances)
        max_dist = max(dist for _, dist in point_distances)
        threshold = min_dist + (max_dist - min_dist) * 0.3

        close_point_indices = []
        for i in range(len(points)):
            if point_distances[i][1] <= threshold:
                close_point_indices.append(i)

        if len(close_point_indices) == 0:
            return None

        sequences = []
        current_seq = [close_point_indices[0]]

        for i in range(1, len(close_point_indices)):
            if close_point_indices[i] == close_point_indices[i - 1] + 1:
                current_seq.append(close_point_indices[i])
            else:
                sequences.append(current_seq)
                current_seq = [close_point_indices[i]]
        sequences.append(current_seq)

        longest_seq = max(sequences, key=len)

        line = ogr.Geometry(ogr.wkbLineString)
        for idx in longest_seq:
            line.AddPoint(exterior_ring.GetX(idx), exterior_ring.GetY(idx))

        return line

    result_lines = []
    geom_type = ogr.GT_Flatten(clipped_geom.GetGeometryType())

    if geom_type == ogr.wkbPolygon:
        result = process_polygon(clipped_geom)
        if result:
            result_lines.append(result)
    elif geom_type == ogr.wkbMultiPolygon:
        for i in range(clipped_geom.GetGeometryCount()):
            result = process_polygon(clipped_geom.GetGeometryRef(i))
            if result:
                result_lines.append(result)

    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbLineString)

    for line in result_lines:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(line)
        output_layer.CreateFeature(out_feature)
        out_feature = None

    ds = None
    print(f"Created {len(result_lines)} edge line(s)")


def create_division_points(gpkg_path, line_layer_name, output_layer_name, preferred_width):
    ds = ogr.Open(gpkg_path, 1)

    line_layer = ds.GetLayerByName(line_layer_name)
    srs = line_layer.GetSpatialRef()

    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbPoint)

    line_id_field = ogr.FieldDefn("line_id", ogr.OFTInteger)
    output_layer.CreateField(line_id_field)

    point_id_field = ogr.FieldDefn("point_id", ogr.OFTInteger)
    output_layer.CreateField(point_id_field)

    width_field = ogr.FieldDefn("width", ogr.OFTReal)
    output_layer.CreateField(width_field)

    n_bands_field = ogr.FieldDefn("n_bands", ogr.OFTInteger)
    output_layer.CreateField(n_bands_field)

    line_id = 0
    total_points = 0

    for feature in line_layer:
        line_geom = feature.GetGeometryRef()
        L = line_geom.Length()

        if L < preferred_width:
            continue

        best_n = None
        best_diff = float("inf")

        for n in range(3, 101, 2):
            w = L / n
            diff = abs(w - preferred_width)
            if diff < best_diff:
                best_diff = diff
                best_n = n

        w = L / best_n
        n_full_bands = int(L / w)

        for k in range(n_full_bands - 1):
            distance = (k + 1) * w

            point = line_geom.Value(distance)

            out_feature = ogr.Feature(output_layer.GetLayerDefn())
            out_feature.SetGeometry(point)
            out_feature.SetField("line_id", line_id)
            out_feature.SetField("point_id", k)
            out_feature.SetField("width", w)
            out_feature.SetField("n_bands", n_full_bands)
            output_layer.CreateFeature(out_feature)
            out_feature = None
            total_points += 1

        line_id += 1

    ds = None
    print(
        f"Created {total_points} points from {line_id} lines with average width {preferred_width}m"
    )


def create_perpendicular_lines(
    gpkg_path,
    points_layer_name,
    lines_layer_name,
    site_layer_name,
    output_layer_name,
    perpendicular_length,
):
    ds = ogr.Open(gpkg_path, 1)

    points_layer = ds.GetLayerByName(points_layer_name)
    lines_layer = ds.GetLayerByName(lines_layer_name)
    site_layer = ds.GetLayerByName(site_layer_name)

    srs = points_layer.GetSpatialRef()

    site_geoms = []
    for feature in site_layer:
        site_geoms.append(feature.GetGeometryRef().Clone())

    site_union = site_geoms[0]
    for geom in site_geoms[1:]:
        site_union = site_union.Union(geom)

    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbLineString)

    line_id_field = ogr.FieldDefn("line_id", ogr.OFTInteger)
    output_layer.CreateField(line_id_field)

    point_id_field = ogr.FieldDefn("point_id", ogr.OFTInteger)
    output_layer.CreateField(point_id_field)

    perp_count = 0

    for line_feature in lines_layer:
        line_geom = line_feature.GetGeometryRef()
        line_fid = line_feature.GetFID()

        points_layer.ResetReading()
        for point_feature in points_layer:
            point_geom = point_feature.GetGeometryRef()

            if line_geom.Distance(point_geom) < 1.0:
                point_id = point_feature.GetField("point_id")

                x = point_geom.GetX()
                y = point_geom.GetY()

                sample_dist = 1.0
                line_length = line_geom.Length()

                closest_dist = float("inf")
                closest_point_dist = 0

                step = max(1.0, line_length / 100)
                current = 0
                while current <= line_length:
                    test_point = line_geom.Value(current)
                    dist = ((test_point.GetX() - x) ** 2 + (test_point.GetY() - y) ** 2) ** 0.5
                    if dist < closest_dist:
                        closest_dist = dist
                        closest_point_dist = current
                    current += step

                point_before = line_geom.Value(max(0, closest_point_dist - sample_dist))
                point_after = line_geom.Value(min(line_length, closest_point_dist + sample_dist))

                dx = point_after.GetX() - point_before.GetX()
                dy = point_after.GetY() - point_before.GetY()

                length = (dx**2 + dy**2) ** 0.5
                if length < 0.001:
                    continue

                dx_norm = dx / length
                dy_norm = dy / length

                perp_dx = -dy_norm
                perp_dy = dx_norm

                x1 = x + perp_dx * perpendicular_length
                y1 = y + perp_dy * perpendicular_length
                x2 = x - perp_dx * perpendicular_length
                y2 = y - perp_dy * perpendicular_length

                perp_line = ogr.Geometry(ogr.wkbLineString)
                perp_line.AddPoint(x1, y1)
                perp_line.AddPoint(x2, y2)

                clipped_line = perp_line.Intersection(site_union)

                if not clipped_line.IsEmpty():
                    geom_type = ogr.GT_Flatten(clipped_line.GetGeometryType())

                    if geom_type == ogr.wkbLineString:
                        out_feature = ogr.Feature(output_layer.GetLayerDefn())
                        out_feature.SetGeometry(clipped_line)
                        out_feature.SetField("line_id", line_fid)
                        out_feature.SetField("point_id", point_id)
                        output_layer.CreateFeature(out_feature)
                        out_feature = None
                        perp_count += 1
                    elif geom_type == ogr.wkbMultiLineString:
                        for i in range(clipped_line.GetGeometryCount()):
                            line_part = clipped_line.GetGeometryRef(i)
                            out_feature = ogr.Feature(output_layer.GetLayerDefn())
                            out_feature.SetGeometry(line_part)
                            out_feature.SetField("line_id", line_fid)
                            out_feature.SetField("point_id", point_id)
                            output_layer.CreateFeature(out_feature)
                            out_feature = None
                            perp_count += 1

    ds = None
    print(f"Created {perp_count} perpendicular lines clipped by site")


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

    num_rings = calculate_required_rings(
        output_path, site_layer_name, cfg.off_grid_partitions_preferred_depth
    )
    multiring_buffer(
        output_path,
        "arterial_buffered_large",
        num_rings,
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

    print("Step 11: Extracting arterial edge lines...")
    extract_arterial_edge_lines(
        output_path,
        "arterial_roads",
        "arterial_road_setback",
        "secondary_road_setback",
        "arterial_edge_lines",
        clip_buffer=0.1,
    )

    print("Step 12: Creating division points along arterial edges...")
    create_division_points(
        output_path,
        "arterial_edge_lines",
        "division_points",
        cfg.off_grid_partitions_preferred_width,
    )

    print("Step 13: Creating perpendicular lines from division points...")
    create_perpendicular_lines(
        output_path,
        "division_points",
        "arterial_edge_lines",
        site_layer_name,
        "perpendicular_lines",
        cfg.perpendicular_line_length,
    )

    # Clean up intermediate layers
    print("Cleaning up intermediate layers...")
    final_layers = [
        "arterial_roads",
        "arterial_road_setback",
        "street_blocks",
        "secondary_road_setback",
        "arterial_edge_lines",
        "division_points",
        "perpendicular_lines",
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
