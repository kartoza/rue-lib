import math

from osgeo import ogr

from rue_lib.core.helpers import create_or_replace_layer


def line_end_intersects_buffer(line: ogr.Geometry, buffer_geom: ogr.Geometry) -> bool:
    if buffer_geom is None or line is None:
        return False

    n_points = line.GetPointCount()
    if n_points < 2:
        return False
    end_x, end_y, _ = line.GetPoint(n_points - 1)
    end_pt = ogr.Geometry(ogr.wkbPoint)
    end_pt.AddPoint(end_x, end_y)
    return end_pt.Intersects(buffer_geom)


def split_polygon_by_line(poly: ogr.Geometry, line: ogr.Geometry, buffer_width: float = 0.00000001):
    """
    Split a polygon using a buffered line.

    Returns a list of polygon parts. If the line does not split the polygon,
    returns [poly].
    """
    if not poly.Intersects(line):
        return [poly.Clone()]

    cut_strip = line.Buffer(buffer_width)
    diff = poly.Difference(cut_strip)

    result_polys = []

    gtype = diff.GetGeometryType()
    if gtype == ogr.wkbPolygon:
        result_polys.append(diff.Clone())
    elif gtype in (ogr.wkbMultiPolygon, ogr.wkbGeometryCollection):
        for i in range(diff.GetGeometryCount()):
            g = diff.GetGeometryRef(i)
            if g is not None and g.GetGeometryType() == ogr.wkbPolygon:
                result_polys.append(g.Clone())
    if not result_polys:
        result_polys.append(poly.Clone())

    return result_polys


def find_concave_points(
    input_gpkg: str,
    erased_grid_layer_name: str,
    boundary_points_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> str:
    """
    Find concave points from the boundary vertices.

    A concave point is a vertex where the interior angle is greater than 180 degrees.
    This is determined by calculating the cross product of consecutive edge vectors.

    Args:
        input_gpkg: Path to input GeoPackage
        erased_grid_layer_name: Name of erased cold grid layer
        boundary_points_layer_name: Name of boundary points layer
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output concave points layer

    Returns:
        Name of the output layer
    """
    ds = ogr.Open(input_gpkg, 0)
    if ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    grid_layer = ds.GetLayerByName(erased_grid_layer_name)
    if grid_layer is None:
        raise ValueError(f"Layer {erased_grid_layer_name} not found")

    points_layer = ds.GetLayerByName(boundary_points_layer_name)
    if points_layer is None:
        raise ValueError(f"Layer {boundary_points_layer_name} not found")
    srs = grid_layer.GetSpatialRef()

    print(f"  Processing {grid_layer.GetFeatureCount()} blocks...")
    concave_points = []
    block_id = 0
    for grid_feat in grid_layer:
        block_id += 1
        grid_geom = grid_feat.GetGeometryRef()
        if grid_geom is None:
            continue
        boundary = grid_geom.GetBoundary()
        if boundary is None:
            continue
        if boundary.GetGeometryType() == ogr.wkbLineString:
            lines = [boundary]
        elif boundary.GetGeometryType() == ogr.wkbMultiLineString:
            lines = [boundary.GetGeometryRef(i) for i in range(boundary.GetGeometryCount())]
        else:
            continue
        for line in lines:
            point_count = line.GetPointCount()
            if point_count < 3:
                continue
            vertices = []
            for i in range(point_count):
                x = line.GetX(i)
                y = line.GetY(i)
                vertices.append((x, y))

            for i in range(point_count - 1):
                p0 = vertices[(i - 1) % (point_count - 1)]
                p1 = vertices[i]
                p2 = vertices[(i + 1) % (point_count - 1)]
                v1_x = p1[0] - p0[0]
                v1_y = p1[1] - p0[1]
                v2_x = p2[0] - p1[0]
                v2_y = p2[1] - p1[1]
                cross = v1_x * v2_y - v1_y * v2_x
                if cross > 15000:
                    point_geom = ogr.Geometry(ogr.wkbPoint)
                    point_geom.AddPoint(p1[0], p1[1])
                    points_layer.SetSpatialFilter(point_geom.Buffer(0.01))
                    is_boundary_point = False
                    for pt_feat in points_layer:
                        if pt_feat.GetField("block_id") == block_id:
                            is_boundary_point = True
                            break
                    points_layer.ResetReading()

                    if is_boundary_point:
                        concave_points.append(
                            {
                                "x": p1[0],
                                "y": p1[1],
                                "block_id": block_id,
                                "vertex_id": i,
                                "cross_product": cross,
                            }
                        )
    grid_layer = None
    points_layer = None
    ds = None
    out_ds = ogr.Open(output_gpkg, 1)
    if out_ds is None:
        raise ValueError(f"Could not open {output_gpkg} for writing")
    for i in range(out_ds.GetLayerCount()):
        layer = out_ds.GetLayerByIndex(i)
        if layer.GetName() == output_layer_name:
            out_ds.DeleteLayer(i)
            break
    out_layer = out_ds.CreateLayer(output_layer_name, srs, ogr.wkbPoint)
    out_layer.CreateField(ogr.FieldDefn("block_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("vertex_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("cross_product", ogr.OFTReal))
    for point in concave_points:
        pt_geom = ogr.Geometry(ogr.wkbPoint)
        pt_geom.AddPoint(point["x"], point["y"])

        out_feat = ogr.Feature(out_layer.GetLayerDefn())
        out_feat.SetGeometry(pt_geom)
        out_feat.SetField("block_id", point["block_id"])
        out_feat.SetField("vertex_id", point["vertex_id"])
        out_feat.SetField("cross_product", point["cross_product"])

        out_layer.CreateFeature(out_feat)
        out_feat = None
    out_layer = None
    out_ds = None

    total_concave = len(concave_points)
    print(f"  Found {total_concave} concave points")
    print(f"  Created layer: {output_layer_name}")
    return output_layer_name


def split_buffered_lines_by_subdivided_blocks(
    buffered_line_features: list[dict],
    subdivided_blocks: list[dict],
) -> list[dict]:
    """Split buffered lines by intersecting with subdivided blocks.

    Args:
        buffered_line_features: List of buffered line features, each containing:
            - geometry: OGR geometry of the buffered line
            - block_id: Original block ID
            - road_type: Type of road (arterial/secondary/local)
        subdivided_blocks: List of subdivided block features, each containing:
            - geometry: OGR geometry of the subdivided block
            - orig_block_id: Original block ID before subdivision
            - (other fields not used in this function)

    Returns:
        List of split buffered line features, each containing:
            - geometry: OGR geometry of the split buffered line
            - block_id: Original block ID
            - road_type: Type of road

    Note:
        - Handles Polygon, MultiPolygon, and GeometryCollection result types
        - MultiPolygons are split into individual polygons
        - Geometry collections are traversed to extract polygon geometries
        - If no intersection occurs, the original buffered line is kept
    """
    sdv_by_block = {}
    for _block in subdivided_blocks:
        block_id_for_line = _block["orig_block_id"]
        sdv_by_block.setdefault(block_id_for_line, []).append((_block["geometry"], _block["id"]))

    split_buffered_lines = []
    for buff_line in buffered_line_features:
        buff_geom = buff_line["geometry"]
        buff_block_id = buff_line["block_id"]
        buff_road_type = buff_line["road_type"]
        sdv_blocks_for_block = sdv_by_block.get(buff_block_id, [])
        if not sdv_blocks_for_block:
            split_buffered_lines.append(buff_line)
            continue
        for subdivided_block_geom, subdivided_block_id in sdv_blocks_for_block:
            if not buff_geom.Intersects(subdivided_block_geom):
                continue

            try:
                intersection = buff_geom.Intersection(subdivided_block_geom)
                if intersection and not intersection.IsEmpty():
                    geom_type = intersection.GetGeometryType()
                    if geom_type == ogr.wkbPolygon or geom_type == ogr.wkbPolygon25D:
                        split_buffered_lines.append(
                            {
                                "geometry": intersection.Clone(),
                                "block_id": buff_block_id,
                                "road_type": buff_road_type,
                                "id": subdivided_block_id,
                            }
                        )
                    elif geom_type == ogr.wkbMultiPolygon or geom_type == ogr.wkbMultiPolygon25D:
                        for i in range(intersection.GetGeometryCount()):
                            poly = intersection.GetGeometryRef(i)
                            if poly and not poly.IsEmpty():
                                split_buffered_lines.append(
                                    {
                                        "geometry": poly.Clone(),
                                        "block_id": buff_block_id,
                                        "road_type": buff_road_type,
                                        "id": subdivided_block_id,
                                    }
                                )
                    elif geom_type == ogr.wkbGeometryCollection:
                        # Extract polygons from geometry collection
                        for i in range(intersection.GetGeometryCount()):
                            geom = intersection.GetGeometryRef(i)
                            gt = geom.GetGeometryType()
                            if (
                                gt == ogr.wkbPolygon
                                or gt == ogr.wkbPolygon25D
                                or gt == ogr.wkbMultiPolygon
                                or gt == ogr.wkbMultiPolygon25D
                            ):
                                if gt == ogr.wkbMultiPolygon or gt == ogr.wkbMultiPolygon25D:
                                    for j in range(geom.GetGeometryCount()):
                                        poly = geom.GetGeometryRef(j)
                                        if poly and not poly.IsEmpty():
                                            split_buffered_lines.append(
                                                {
                                                    "geometry": poly.Clone(),
                                                    "block_id": buff_block_id,
                                                    "road_type": buff_road_type,
                                                    "id": subdivided_block_id,
                                                }
                                            )
                                else:
                                    split_buffered_lines.append(
                                        {
                                            "geometry": geom.Clone(),
                                            "block_id": buff_block_id,
                                            "road_type": buff_road_type,
                                            "id": subdivided_block_id,
                                        }
                                    )
            except RuntimeError as e:
                print(
                    f"    Warning: Failed to intersect buffered line "
                    f"with block {buff_block_id}: {e}"
                )
                continue

    return split_buffered_lines


def subdivide_blocks_by_concave_points(
    input_gpkg: str,
    erased_grid_layer_name: str,
    concave_points_layer_name: str,
    boundary_points_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
    offset_distance: float = 10.0,
    buffered_lines_layer_name: str = "",
) -> str:
    """
    Subdivide blocks at concave corners by creating perpendicular cutting lines.

    For each concave point in boundary_points:
    1. Find prev and next points in boundary_points sequence
    2. Move offset_distance (10m) from concave point towards prev and next
    3. Create perpendicular lines from those positions
    4. Output the cutting lines
    5. Optionally split clipped buffered lines using the same cutting lines

    Args:
        input_gpkg: Path to input GeoPackage
        erased_grid_layer_name: Name of erased cold grid layer
        concave_points_layer_name: Name of concave points layer
        boundary_points_layer_name: Name of boundary points layer
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output cutting lines layer
        offset_distance: Distance to move from concave point (default 10m)
        buffered_lines_layer_name: Optional layer name of clipped buffered lines to split
    Returns:
        Name of the output layer
    """
    ds = ogr.Open(input_gpkg, 0)
    if ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    grid_layer = ds.GetLayerByName(erased_grid_layer_name)
    if grid_layer is None:
        raise ValueError(f"Layer {erased_grid_layer_name} not found")
    concave_layer = ds.GetLayerByName(concave_points_layer_name)
    if concave_layer is None:
        raise ValueError(f"Layer {concave_points_layer_name} not found")
    boundary_layer = ds.GetLayerByName(boundary_points_layer_name)
    if boundary_layer is None:
        raise ValueError(f"Layer {boundary_points_layer_name} not found")
    srs = grid_layer.GetSpatialRef()
    print(f"  Processing {grid_layer.GetFeatureCount()} blocks...")

    concave_points_set = set()
    for concave_feat in concave_layer:
        block_id = concave_feat.GetField("block_id")
        concave_geom = concave_feat.GetGeometryRef()
        concave_x = concave_geom.GetX()
        concave_y = concave_geom.GetY()
        concave_points_set.add((block_id, concave_x, concave_y))

    boundary_points_by_block = {}
    for boundary_feat in boundary_layer:
        block_id = boundary_feat.GetField("block_id")
        vertex_id = boundary_feat.GetField("vertex_id")
        boundary_geom = boundary_feat.GetGeometryRef()
        x = boundary_geom.GetX()
        y = boundary_geom.GetY()

        if block_id not in boundary_points_by_block:
            boundary_points_by_block[block_id] = []

        boundary_points_by_block[block_id].append(
            {
                "x": x,
                "y": y,
                "vertex_id": vertex_id,
            }
        )

    for block_id in boundary_points_by_block:
        boundary_points_by_block[block_id].sort(key=lambda p: p["vertex_id"])

    grid_geometries = {}
    block_id = 0
    for grid_feat in grid_layer:
        block_id += 1
        grid_geom = grid_feat.GetGeometryRef()
        if grid_geom:
            grid_geometries[block_id] = grid_geom.Clone()

    grid_layer = None
    concave_layer = None
    boundary_layer = None
    ds = None

    cutting_lines_to_write = []
    line_id = 0

    for block_id, boundary_points in boundary_points_by_block.items():
        if block_id not in grid_geometries:
            continue

        block_geom = grid_geometries[block_id]

        for i, point in enumerate(boundary_points):
            point_x = point["x"]
            point_y = point["y"]

            is_concave = False
            for concave_block_id, concave_x, concave_y in concave_points_set:
                if (
                    concave_block_id == block_id
                    and abs(point_x - concave_x) < 0.01
                    and abs(point_y - concave_y) < 0.01
                ):
                    is_concave = True
                    break

            if not is_concave:
                continue

            num_points = len(boundary_points)
            prev_idx = (i - 1) % num_points
            next_idx = (i + 1) % num_points

            prev_point = boundary_points[prev_idx]
            next_point = boundary_points[next_idx]

            prev_x = prev_point["x"]
            prev_y = prev_point["y"]
            next_x = next_point["x"]
            next_y = next_point["y"]

            concave_point_geom = ogr.Geometry(ogr.wkbPoint)
            concave_point_geom.AddPoint(point_x, point_y)

            dir1_x = prev_x - point_x
            dir1_y = prev_y - point_y
            dir1_length = math.sqrt(dir1_x**2 + dir1_y**2)

            if dir1_length < 0.001:
                continue

            dir1_x /= dir1_length
            dir1_y /= dir1_length
            offset1_x = point_x + dir1_x * offset_distance
            offset1_y = point_y + dir1_y * offset_distance

            dir2_x = next_x - point_x
            dir2_y = next_y - point_y
            dir2_length = math.sqrt(dir2_x**2 + dir2_y**2)

            if dir2_length < 0.001:
                continue

            dir2_x /= dir2_length
            dir2_y /= dir2_length
            offset2_x = point_x + dir2_x * offset_distance
            offset2_y = point_y + dir2_y * offset_distance

            perp1_x = -dir1_y
            perp1_y = dir1_x

            test_line1_pos = ogr.Geometry(ogr.wkbLineString)
            test_line1_pos.AddPoint(offset1_x, offset1_y)
            test_line1_pos.AddPoint(offset1_x + perp1_x * 5.0, offset1_y + perp1_y * 5.0)

            test_line1_neg = ogr.Geometry(ogr.wkbLineString)
            test_line1_neg.AddPoint(offset1_x, offset1_y)
            test_line1_neg.AddPoint(offset1_x - perp1_x * 5.0, offset1_y - perp1_y * 5.0)

            direction1_x, direction1_y = None, None
            if block_geom and line_end_intersects_buffer(test_line1_pos, block_geom):
                direction1_x, direction1_y = perp1_x, perp1_y
            elif block_geom and line_end_intersects_buffer(test_line1_neg, block_geom):
                direction1_x, direction1_y = -perp1_x, -perp1_y

            if direction1_x is not None:
                current_length = 1.0
                while current_length < 1000.0:
                    test_line = ogr.Geometry(ogr.wkbLineString)
                    test_line.AddPoint(offset1_x, offset1_y)
                    test_line.AddPoint(
                        offset1_x + direction1_x * current_length,
                        offset1_y + direction1_y * current_length,
                    )

                    if not line_end_intersects_buffer(test_line, block_geom):
                        break

                    current_length += 5.0

                line1 = ogr.Geometry(ogr.wkbLineString)
                line1.AddPoint(offset1_x, offset1_y)
                line1.AddPoint(
                    offset1_x + direction1_x * max(current_length, 1.0),
                    offset1_y + direction1_y * max(current_length, 1.0),
                )

                cutting_lines_to_write.append(
                    {
                        "line_id": line_id,
                        "geometry": line1,
                        "block_id": block_id,
                        "concave_x": point_x,
                        "concave_y": point_y,
                        "line_type": "line1",
                    }
                )
                line_id += 1

            perp2_x = -dir2_y
            perp2_y = dir2_x

            test_line2_pos = ogr.Geometry(ogr.wkbLineString)
            test_line2_pos.AddPoint(offset2_x, offset2_y)
            test_line2_pos.AddPoint(offset2_x + perp2_x * 5.0, offset2_y + perp2_y * 5.0)

            test_line2_neg = ogr.Geometry(ogr.wkbLineString)
            test_line2_neg.AddPoint(offset2_x, offset2_y)
            test_line2_neg.AddPoint(offset2_x - perp2_x * 5.0, offset2_y - perp2_y * 5.0)

            direction2_x, direction2_y = None, None
            if block_geom and line_end_intersects_buffer(test_line2_pos, block_geom):
                direction2_x, direction2_y = perp2_x, perp2_y
            elif block_geom and line_end_intersects_buffer(test_line2_neg, block_geom):
                direction2_x, direction2_y = -perp2_x, -perp2_y

            if direction2_x is not None:
                current_length = 1.0
                while current_length < 200.0:
                    test_line = ogr.Geometry(ogr.wkbLineString)
                    test_line.AddPoint(offset2_x, offset2_y)
                    test_line.AddPoint(
                        offset2_x + direction2_x * current_length,
                        offset2_y + direction2_y * current_length,
                    )

                    if not line_end_intersects_buffer(test_line, block_geom):
                        break

                    current_length += 5.0

                line2 = ogr.Geometry(ogr.wkbLineString)
                line2.AddPoint(offset2_x, offset2_y)
                line2.AddPoint(
                    offset2_x + direction2_x * max(current_length, 1.0),
                    offset2_y + direction2_y * max(current_length, 1.0),
                )

                cutting_lines_to_write.append(
                    {
                        "line_id": line_id,
                        "geometry": line2,
                        "block_id": block_id,
                        "concave_x": point_x,
                        "concave_y": point_y,
                        "line_type": "line2",
                    }
                )
                line_id += 1

    concave_points_by_block = {}
    for c_block_id, concave_x, concave_y in concave_points_set:
        pt = ogr.Geometry(ogr.wkbPoint)
        pt.AddPoint(concave_x, concave_y)
        concave_points_by_block.setdefault(c_block_id, []).append(
            {"geom": pt, "x": concave_x, "y": concave_y}
        )

    working_blocks = {}
    for b_id, geom in grid_geometries.items():
        working_blocks[(b_id, 0)] = geom.Clone()

    subdivided_blocks = []
    for line_data in cutting_lines_to_write:
        line_geom = line_data["geometry"]
        block_id_for_line = line_data["block_id"]
        concave_x = line_data["concave_x"]
        concave_y = line_data["concave_y"]

        conc_pt = ogr.Geometry(ogr.wkbPoint)
        conc_pt.AddPoint(concave_x, concave_y)

        keys_for_block = [k for k in working_blocks.keys() if k[0] == block_id_for_line]

        for key in keys_for_block:
            orig_block_id, local_id = key
            poly = working_blocks.pop(key)

            if not poly.Intersects(line_geom):
                working_blocks[key] = poly
                continue

            parts = split_polygon_by_line(poly, line_geom)

            for part in parts:
                new_local_id = (
                    max(
                        [lid for (bid, lid) in working_blocks.keys() if bid == orig_block_id]
                        or [local_id]
                    )
                    + 1
                )
                working_blocks[(orig_block_id, new_local_id)] = part.Clone()

            poly = None
            break

    block_id = 0
    for (orig_block_id, _local_id), poly in working_blocks.items():
        is_concave = 0
        for concave_point_in_block in concave_points_by_block.get(orig_block_id, []):
            pt = concave_point_in_block["geom"]
            concave_x = concave_point_in_block["x"]
            concave_y = concave_point_in_block["y"]
            if poly.Intersects(pt):
                is_concave = 1
                break
        subdivided_blocks.append(
            {
                "id": block_id,
                "orig_block_id": orig_block_id,
                "geometry": poly.Clone(),
                "is_concave": is_concave,
                "concave_x": is_concave and concave_x or None,
                "concave_y": is_concave and concave_y or None,
                "type": is_concave and "corner" or "block",
                "block_type": "cold",
            }
        )
        block_id += 1

    out_ds = ogr.Open(output_gpkg, 1)
    if out_ds is None:
        raise ValueError(f"Could not open {output_gpkg} for writing")

    lines_layer = create_or_replace_layer(out_ds, output_layer_name, srs, ogr.wkbLineString)

    lines_layer.CreateField(ogr.FieldDefn("line_id", ogr.OFTInteger))
    lines_layer.CreateField(ogr.FieldDefn("block_id", ogr.OFTInteger))
    lines_layer.CreateField(ogr.FieldDefn("concave_x", ogr.OFTReal))
    lines_layer.CreateField(ogr.FieldDefn("concave_y", ogr.OFTReal))
    lines_layer.CreateField(ogr.FieldDefn("line_type", ogr.OFTString))

    for line_data in cutting_lines_to_write:
        out_feat = ogr.Feature(lines_layer.GetLayerDefn())
        out_feat.SetGeometry(line_data["geometry"])
        out_feat.SetField("line_id", line_data["line_id"])
        out_feat.SetField("block_id", line_data["block_id"])
        out_feat.SetField("concave_x", line_data["concave_x"])
        out_feat.SetField("concave_y", line_data["concave_y"])
        out_feat.SetField("line_type", line_data["line_type"])
        lines_layer.CreateFeature(out_feat)
        out_feat = None

    blocks_layer_name = f"{output_layer_name}_blocks"
    blocks_layer = create_or_replace_layer(out_ds, blocks_layer_name, srs, ogr.wkbPolygon)

    blocks_layer.CreateField(ogr.FieldDefn("orig_id", ogr.OFTInteger))
    blocks_layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
    blocks_layer.CreateField(ogr.FieldDefn("is_concave", ogr.OFTInteger))
    blocks_layer.CreateField(ogr.FieldDefn("concave_x", ogr.OFTReal))
    blocks_layer.CreateField(ogr.FieldDefn("concave_y", ogr.OFTReal))
    blocks_layer.CreateField(ogr.FieldDefn("type", ogr.OFTString))
    blocks_layer.CreateField(ogr.FieldDefn("block_type", ogr.OFTString))

    for blk in subdivided_blocks:
        feat = ogr.Feature(blocks_layer.GetLayerDefn())
        feat.SetGeometry(blk["geometry"])
        feat.SetField("orig_id", blk["orig_block_id"])
        feat.SetField("is_concave", blk["is_concave"])
        feat.SetField("type", blk["type"])
        feat.SetField("block_type", blk["block_type"])
        feat.SetField("id", blk["id"])

        if blk["concave_x"] is not None:
            feat.SetField("concave_x", blk["concave_x"])
        if blk["concave_y"] is not None:
            feat.SetField("concave_y", blk["concave_y"])

        blocks_layer.CreateFeature(feat)
        feat = None

    buffered_lines_layer = None
    for i in range(out_ds.GetLayerCount()):
        layer = out_ds.GetLayerByIndex(i)
        if layer.GetName() == buffered_lines_layer_name:
            buffered_lines_layer = layer
            break

    if buffered_lines_layer:
        buffered_line_features = []
        for feat in buffered_lines_layer:
            geom = feat.GetGeometryRef()
            if geom:
                buffered_line_features.append(
                    {
                        "geometry": geom.Clone(),
                        "block_id": feat.GetField("block_id"),
                        "road_type": feat.GetField("road_type"),
                    }
                )

        split_buffered_lines = split_buffered_lines_by_subdivided_blocks(
            buffered_line_features, subdivided_blocks
        )

        split_buffered_layer_name = f"{blocks_layer_name}_on_grid"
        split_buffered_layer = create_or_replace_layer(
            out_ds, split_buffered_layer_name, srs, ogr.wkbPolygon
        )

        split_buffered_layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))
        split_buffered_layer.CreateField(ogr.FieldDefn("block_id", ogr.OFTInteger))
        split_buffered_layer.CreateField(ogr.FieldDefn("road_type", ogr.OFTString))
        split_buffered_layer.CreateField(ogr.FieldDefn("type", ogr.OFTString))

        for buff_data in split_buffered_lines:
            feat = ogr.Feature(split_buffered_layer.GetLayerDefn())
            feat.SetGeometry(buff_data["geometry"])
            feat.SetField("id", buff_data["id"])
            feat.SetField("block_id", buff_data["block_id"])
            feat.SetField("road_type", buff_data["road_type"])
            feat.SetField("type", "loc_loc")
            split_buffered_layer.CreateFeature(feat)
            feat = None

        split_buffered_layer = None
        print(f"  Split buffered lines layer: {split_buffered_layer_name}")

        off_grid_layer_name = f"{blocks_layer_name}_off_grid"

        for i in range(out_ds.GetLayerCount()):
            layer = out_ds.GetLayerByIndex(i)
            if layer.GetName() == off_grid_layer_name:
                out_ds.DeleteLayer(i)
                break

        remaining_features = []
        for blk in subdivided_blocks:
            block_geom = blk["geometry"]
            orig_block_id = blk["orig_block_id"]

            intersecting_buffer = None
            for buff_data in buffered_line_features:
                if buff_data["block_id"] == orig_block_id:
                    intersecting_buffer = buff_data["geometry"]
                    break

            if not intersecting_buffer:
                continue

            remaining_geom = block_geom.Clone()
            remaining_geom = remaining_geom.Buffer(0).Difference(
                intersecting_buffer.Buffer(0.00001)
            )

            if remaining_geom and not remaining_geom.IsEmpty():
                remaining_features.append(
                    {
                        "id": orig_block_id,
                        "geometry": remaining_geom,
                        "orig_id": orig_block_id,
                        "is_concave": blk["is_concave"],
                        "concave_x": blk["concave_x"],
                        "concave_y": blk["concave_y"],
                        "type": blk["type"],
                        "block_type": blk["block_type"],
                    }
                )

        buffered_lines_layer = None

        remaining_layer = create_or_replace_layer(
            out_ds,
            off_grid_layer_name,
            srs,
            ogr.wkbPolygon,
            [
                ("id", ogr.OFTInteger),
                ("orig_id", ogr.OFTInteger),
                ("is_concave", ogr.OFTInteger),
                ("concave_x", ogr.OFTReal),
                ("concave_y", ogr.OFTReal),
                ("type", ogr.OFTString),
                ("block_type", ogr.OFTString),
            ],
        )
        id = 0
        for rem_data in remaining_features:
            feat = ogr.Feature(remaining_layer.GetLayerDefn())
            feat.SetGeometry(rem_data["geometry"])
            feat.SetField("id", id)
            id += 1
            feat.SetField("orig_id", rem_data["orig_id"])
            feat.SetField("is_concave", rem_data["is_concave"])
            feat.SetField("type", "concave_corner" if rem_data["is_concave"] else "off_grid")
            feat.SetField("block_type", rem_data["block_type"])

            if rem_data["concave_x"] is not None:
                feat.SetField("concave_x", rem_data["concave_x"])
            if rem_data["concave_y"] is not None:
                feat.SetField("concave_y", rem_data["concave_y"])

            remaining_layer.CreateFeature(feat)
            feat = None

        remaining_layer = None
        print(f"  Created off-grid cold boundary layer: {off_grid_layer_name}")

    # Clean up
    lines_layer = None
    blocks_layer = None
    out_ds = None

    total_lines = len(cutting_lines_to_write)
    print(f"  Created {total_lines} cutting lines")
    print(f"  Created layer: {output_layer_name}")
    print(f"  Created subdivided blocks layer: {blocks_layer_name}")
    return output_layer_name
