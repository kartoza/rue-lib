# src/rue_lib/streets/lines.py
from osgeo import ogr


def extract_arterial_edge_lines(
    gpkg_path,
    arterial_roads_layer_name,
    arterial_setback_layer_name,
    secondary_setback_layer_name,
    output_layer_name,
    clip_buffer=0.1,
    sample_distance=5.0,
):
    """Extract edge lines from arterial setback zones, excluding secondary setback overlap.

    Args:
        gpkg_path (str): Path to the GeoPackage
        arterial_roads_layer_name (str): Name of arterial roads layer
        arterial_setback_layer_name (str): Name of arterial setback layer
        secondary_setback_layer_name (str): Name of secondary setback layer
        output_layer_name (str): Name of output edge lines layer
        clip_buffer (float): Buffer distance for clipping (default 0.1)
        sample_distance (float): Sampling distance for edge detection (default 5.0)

    Returns:
        None
    """
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

        # Use a more generous threshold to capture more of the edge
        # Increased from 0.3 to 0.5 to get more of the edge line
        threshold = min_dist + (max_dist - min_dist) * 0.5

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

        # Handle wrapping: if first and last sequences are continuous across the ring boundary
        # (i.e., last index is n-1 and first index is 0), merge them
        if len(sequences) > 1:
            first_seq = sequences[0]
            last_seq = sequences[-1]
            if first_seq[0] == 0 and last_seq[-1] == len(points) - 1:
                # Merge: last sequence + first sequence
                sequences[-1] = last_seq + first_seq
                sequences.pop(0)

        # Return all sequences that are significant (not just the longest)
        # Filter out very short sequences (less than 3 points)
        result_lines = []
        for seq in sequences:
            if len(seq) >= 3:
                line = ogr.Geometry(ogr.wkbLineString)
                for idx in seq:
                    line.AddPoint(exterior_ring.GetX(idx), exterior_ring.GetY(idx))
                result_lines.append(line)

        return result_lines

    result_lines = []
    geom_type = ogr.GT_Flatten(clipped_geom.GetGeometryType())

    if geom_type == ogr.wkbPolygon:
        lines = process_polygon(clipped_geom)
        if lines:
            result_lines.extend(lines)
    elif geom_type == ogr.wkbMultiPolygon:
        for i in range(clipped_geom.GetGeometryCount()):
            lines = process_polygon(clipped_geom.GetGeometryRef(i))
            if lines:
                result_lines.extend(lines)

    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbLineString25D)

    for line in result_lines:
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(line)
        output_layer.CreateFeature(out_feature)
        out_feature = None

    ds = None
    print(f"Created {len(result_lines)} edge line(s)")


def create_division_points(gpkg_path, line_layer_name, output_layer_name, preferred_width):
    """Create division points along lines at regular intervals.

    Args:
        gpkg_path (str): Path to the GeoPackage
        line_layer_name (str): Name of the input line layer
        output_layer_name (str): Name of the output points layer
        preferred_width (float): Preferred spacing between division points

    Returns:
        None
    """
    ds = ogr.Open(gpkg_path, 1)

    line_layer = ds.GetLayerByName(line_layer_name)
    srs = line_layer.GetSpatialRef()

    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbPoint25D)

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

        # Draw first point
        point = line_geom.Value(0)
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(point)
        out_feature.SetField("line_id", line_id)
        out_feature.SetField("point_id", 0)
        out_feature.SetField("width", w)
        out_feature.SetField("n_bands", n_full_bands)
        output_layer.CreateFeature(out_feature)
        out_feature = None
        total_points += 1

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
    """Create perpendicular lines from division points, clipped to site boundary.

    Args:
        gpkg_path (str): Path to the GeoPackage
        points_layer_name (str): Name of the division points layer
        lines_layer_name (str): Name of the reference lines layer
        site_layer_name (str): Name of the site boundary layer
        output_layer_name (str): Name of the output perpendicular lines layer
        perpendicular_length (float): Length of perpendicular lines to create

    Returns:
        None
    """
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

    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbLineString25D)

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
                perp_line.AddPoint(x2, y2)
                perp_line.AddPoint(x1, y1)

                try:
                    clipped_line = perp_line.Intersection(site_union)
                except Exception as e:
                    print(f"Warning: Failed to clip perpendicular line at point {point_id}: {e}")
                    continue

                if not clipped_line.IsEmpty():
                    geom_type = ogr.GT_Flatten(clipped_line.GetGeometryType())

                    if geom_type == ogr.wkbLineString:
                        if clipped_line.GetPointCount() >= 2:
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
                            if line_part.GetPointCount() >= 2:
                                out_feature = ogr.Feature(output_layer.GetLayerDefn())
                                out_feature.SetGeometry(line_part)
                                out_feature.SetField("line_id", line_fid)
                                out_feature.SetField("point_id", point_id)
                                output_layer.CreateFeature(out_feature)
                                out_feature = None
                                perp_count += 1

    ds = None
    print(f"Created {perp_count} perpendicular lines clipped by site")
