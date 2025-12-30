import math
from typing import Union

import geopandas as gpd
from osgeo import ogr
from shapely.geometry import LineString, MultiPolygon, Point, Polygon

from rue_lib.core.definitions import ClusterTypes, ColorTypes


def line_end_intersects_buffer_shapely(
    line: LineString, buffer_geom: Union[Polygon, MultiPolygon]
) -> bool:
    """Check if the end point of a line intersects with a buffer geometry (shapely version)."""
    if buffer_geom is None or line is None or line.is_empty:
        return False

    coords = list(line.coords)
    if len(coords) < 2:
        return False

    end_pt = Point(coords[-1])
    return buffer_geom.intersects(end_pt)


def split_polygon_by_line_shapely(
    poly: Polygon, line: LineString, buffer_width: float = 0.00000001
) -> list[Polygon]:
    """
    Split a polygon using a buffered line (shapely version).

    Returns a list of polygon parts. If the line does not split the polygon,
    returns [poly].
    """
    if not poly.intersects(line):
        return [poly]

    cut_strip = line.buffer(buffer_width)
    diff = poly.difference(cut_strip)

    result_polys = []

    if diff.is_empty:
        return [poly]

    if isinstance(diff, Polygon):
        result_polys.append(diff)
    elif isinstance(diff, MultiPolygon):
        result_polys.extend(list(diff.geoms))

    if not result_polys:
        result_polys.append(poly)

    return result_polys


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
    min_distance: float = 1.0,
) -> str:
    """
    Find concave points from the boundary vertices.

    A concave point is a vertex where the interior angle is between 80 and 100 degrees.
    This is determined by calculating the angle between vectors formed by
    prev->current and current->next points.

    Args:
        input_gpkg: Path to input GeoPackage
        erased_grid_layer_name: Name of erased cold grid layer (not used, kept for compatibility)
        boundary_points_layer_name: Name of boundary points layer
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output concave points layer
        min_distance: Minimum distance to previous/next points (default: 10.0 meters)

    Returns:
        Name of the output layer
    """
    gdf_points = gpd.read_file(input_gpkg, layer=boundary_points_layer_name)
    if gdf_points.empty:
        raise ValueError(f"Layer {boundary_points_layer_name} is empty")
    points_by_block = {}
    for _, row in gdf_points.iterrows():
        block_id = row["block_id"]
        vertex_id = row["vertex_id"]
        geom = row.geometry

        if block_id not in points_by_block:
            points_by_block[block_id] = []

        points_by_block[block_id].append(
            {
                "x": geom.x,
                "y": geom.y,
                "vertex_id": vertex_id,
            }
        )
    for block_id in points_by_block:
        points_by_block[block_id].sort(key=lambda p: p["vertex_id"])
    print(f"  Processing {len(points_by_block)} blocks...")
    concave_points = []
    skipped_too_close = 0
    for block_id, points in points_by_block.items():
        num_points = len(points)

        if num_points < 3:
            continue

        for i, current_point in enumerate(points):
            prev_point = points[(i - 1) % num_points]
            next_point = points[(i + 1) % num_points]
            v1_x = current_point["x"] - prev_point["x"]
            v1_y = current_point["y"] - prev_point["y"]

            v2_x = next_point["x"] - current_point["x"]
            v2_y = next_point["y"] - current_point["y"]

            # Calculate distances to previous and next points
            dist_to_prev = math.sqrt(v1_x * v1_x + v1_y * v1_y)
            dist_to_next = math.sqrt(v2_x * v2_x + v2_y * v2_y)

            # Skip if either distance is too small
            if dist_to_prev < min_distance or dist_to_next < min_distance:
                skipped_too_close += 1
                continue

            dot = v1_x * v2_x + v1_y * v2_y
            cross = v1_x * v2_y - v1_y * v2_x
            angle_rad = math.atan2(cross, dot)
            angle_deg = math.degrees(angle_rad)

            if angle_deg < -80:
                concave_points.append(
                    {
                        "geometry": Point(current_point["x"], current_point["y"]),
                        "block_id": block_id,
                        "vertex_id": current_point["vertex_id"],
                        "angle": angle_deg,
                    }
                )

    if concave_points:
        gdf_concave = gpd.GeoDataFrame(concave_points, geometry="geometry", crs=gdf_points.crs)
    else:
        gdf_concave = gpd.GeoDataFrame([], geometry=[], crs=gdf_points.crs)

    gdf_concave.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")

    total_concave = len(concave_points)
    print(f"  Found {total_concave} concave points")
    if skipped_too_close > 0:
        print(f"  Skipped {skipped_too_close} points (distance < {min_distance})")
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
    idx = 0
    for buff_line in buffered_line_features:
        buff_geom = buff_line["geometry"]
        buff_block_id = buff_line["block_id"]
        sdv_blocks_for_block = sdv_by_block.get(buff_block_id, [])
        if not sdv_blocks_for_block:
            split_buffered_lines.append(buff_line)
            continue
        for subdivided_block_geom, _subdivided_block_id in sdv_blocks_for_block:
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
                                "id": idx,
                            }
                        )
                        idx += 1
                    elif geom_type == ogr.wkbMultiPolygon or geom_type == ogr.wkbMultiPolygon25D:
                        for i in range(intersection.GetGeometryCount()):
                            poly = intersection.GetGeometryRef(i)
                            if poly and not poly.IsEmpty():
                                split_buffered_lines.append(
                                    {
                                        "geometry": poly.Clone(),
                                        "block_id": buff_block_id,
                                        "id": idx,
                                    }
                                )
                    elif geom_type == ogr.wkbGeometryCollection:
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
                                                    "id": idx,
                                                }
                                            )
                                            idx += 1
                                else:
                                    split_buffered_lines.append(
                                        {
                                            "geometry": geom.Clone(),
                                            "block_id": buff_block_id,
                                            "id": idx,
                                        }
                                    )
                                    idx += 1
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
    gdf_grid = gpd.read_file(input_gpkg, layer=erased_grid_layer_name)
    gdf_concave = gpd.read_file(input_gpkg, layer=concave_points_layer_name)
    gdf_boundary = gpd.read_file(input_gpkg, layer=boundary_points_layer_name)

    crs = gdf_grid.crs
    print(f"  Processing {len(gdf_grid)} blocks...")

    concave_points_set = set()
    for _, row in gdf_concave.iterrows():
        block_id = row["block_id"]
        concave_geom = row.geometry
        concave_points_set.add((block_id, concave_geom.x, concave_geom.y))

    # Build boundary points by block
    boundary_points_by_block = {}
    for _, row in gdf_boundary.iterrows():
        block_id = row["block_id"]
        vertex_id = row["vertex_id"]
        geom = row.geometry

        if block_id not in boundary_points_by_block:
            boundary_points_by_block[block_id] = []

        boundary_points_by_block[block_id].append(
            {
                "x": geom.x,
                "y": geom.y,
                "vertex_id": vertex_id,
            }
        )

    for block_id in boundary_points_by_block:
        boundary_points_by_block[block_id].sort(key=lambda p: p["vertex_id"])

    # Build grid geometries dict
    grid_geometries = {}
    for _, row in gdf_grid.iterrows():
        block_id = row["id"]
        grid_geometries[block_id] = row.geometry

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

            test_line1_pos = LineString(
                [(offset1_x, offset1_y), (offset1_x + perp1_x * 5.0, offset1_y + perp1_y * 5.0)]
            )
            test_line1_neg = LineString(
                [(offset1_x, offset1_y), (offset1_x - perp1_x * 5.0, offset1_y - perp1_y * 5.0)]
            )

            direction1_x, direction1_y = None, None
            if block_geom and line_end_intersects_buffer_shapely(test_line1_pos, block_geom):
                direction1_x, direction1_y = perp1_x, perp1_y
            elif block_geom and line_end_intersects_buffer_shapely(test_line1_neg, block_geom):
                direction1_x, direction1_y = -perp1_x, -perp1_y

            if direction1_x is not None:
                current_length = 1.0
                while True:
                    test_line = LineString(
                        [
                            (offset1_x, offset1_y),
                            (
                                offset1_x + direction1_x * current_length,
                                offset1_y + direction1_y * current_length,
                            ),
                        ]
                    )

                    if not line_end_intersects_buffer_shapely(test_line, block_geom):
                        break

                    current_length += 5.0

                line1 = LineString(
                    [
                        (offset1_x, offset1_y),
                        (
                            offset1_x + direction1_x * max(current_length, 1.0),
                            offset1_y + direction1_y * max(current_length, 1.0),
                        ),
                    ]
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

            test_line2_pos = LineString(
                [(offset2_x, offset2_y), (offset2_x + perp2_x * 5.0, offset2_y + perp2_y * 5.0)]
            )
            test_line2_neg = LineString(
                [(offset2_x, offset2_y), (offset2_x - perp2_x * 5.0, offset2_y - perp2_y * 5.0)]
            )

            direction2_x, direction2_y = None, None
            if block_geom and line_end_intersects_buffer_shapely(test_line2_pos, block_geom):
                direction2_x, direction2_y = perp2_x, perp2_y
            elif block_geom and line_end_intersects_buffer_shapely(test_line2_neg, block_geom):
                direction2_x, direction2_y = -perp2_x, -perp2_y

            if direction2_x is not None:
                current_length = 1.0
                while True:
                    test_line = LineString(
                        [
                            (offset2_x, offset2_y),
                            (
                                offset2_x + direction2_x * current_length,
                                offset2_y + direction2_y * current_length,
                            ),
                        ]
                    )

                    if not line_end_intersects_buffer_shapely(test_line, block_geom):
                        break

                    current_length += 5.0

                line2 = LineString(
                    [
                        (offset2_x, offset2_y),
                        (
                            offset2_x + direction2_x * max(current_length, 1.0),
                            offset2_y + direction2_y * max(current_length, 1.0),
                        ),
                    ]
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

    # Build concave points by block
    concave_points_by_block = {}
    for c_block_id, concave_x, concave_y in concave_points_set:
        pt = Point(concave_x, concave_y)
        concave_points_by_block.setdefault(c_block_id, []).append(
            {"geom": pt, "x": concave_x, "y": concave_y}
        )

    # Initialize working blocks
    working_blocks = {}
    for b_id, geom in grid_geometries.items():
        working_blocks[(b_id, 0)] = geom

    # Split blocks by cutting lines
    for line_data in cutting_lines_to_write:
        line_geom = line_data["geometry"]
        block_id_for_line = line_data["block_id"]

        keys_for_block = [k for k in working_blocks.keys() if k[0] == block_id_for_line]

        for key in keys_for_block:
            orig_block_id, local_id = key
            poly = working_blocks.pop(key)

            if not poly.intersects(line_geom):
                working_blocks[key] = poly
                continue

            parts = split_polygon_by_line_shapely(poly, line_geom)

            for part in parts:
                new_local_id = (
                    max(
                        [lid for (bid, lid) in working_blocks.keys() if bid == orig_block_id]
                        or [local_id]
                    )
                    + 1
                )
                working_blocks[(orig_block_id, new_local_id)] = part

            break

    # Build subdivided blocks list
    subdivided_blocks = []
    block_id = 0
    for (orig_block_id, _local_id), poly in working_blocks.items():
        is_concave = 0
        concave_x_val = None
        concave_y_val = None
        for concave_point_in_block in concave_points_by_block.get(orig_block_id, []):
            pt = concave_point_in_block["geom"]
            concave_x = concave_point_in_block["x"]
            concave_y = concave_point_in_block["y"]
            if poly.intersects(pt):
                is_concave = 1
                concave_x_val = concave_x
                concave_y_val = concave_y
                break
        subdivided_blocks.append(
            {
                "id": block_id,
                "orig_block_id": orig_block_id,
                "geometry": poly,
                "is_concave": is_concave,
                "concave_x": concave_x_val,
                "concave_y": concave_y_val,
                "type": "corner" if is_concave else "block",
                "block_type": "cold",
            }
        )
        block_id += 1

    # Write cutting lines to geopackage
    gdf_lines = gpd.GeoDataFrame(cutting_lines_to_write, geometry="geometry", crs=crs)
    gdf_lines.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")

    # Write subdivided blocks
    blocks_layer_name = f"{output_layer_name}_blocks"
    gdf_blocks = gpd.GeoDataFrame(subdivided_blocks, geometry="geometry", crs=crs)
    gdf_blocks.to_file(output_gpkg, layer=blocks_layer_name, driver="GPKG")

    # Handle buffered lines if provided
    if buffered_lines_layer_name:
        try:
            gdf_buffered = gpd.read_file(output_gpkg, layer=buffered_lines_layer_name)

            buffered_line_features = []
            for _, row in gdf_buffered.iterrows():
                buffered_line_features.append({"geometry": row.geometry, "block_id": row["id"]})

            # Split buffered lines by subdivided blocks
            split_buffered_lines = []
            idx = 0
            for buff_line in buffered_line_features:
                buff_geom = buff_line["geometry"]
                buff_block_id = buff_line["block_id"]

                for blk in subdivided_blocks:
                    if blk["orig_block_id"] != buff_block_id:
                        continue

                    subdivided_block_geom = blk["geometry"]

                    if not buff_geom.intersects(subdivided_block_geom):
                        continue

                    try:
                        intersection = buff_geom.intersection(subdivided_block_geom)
                        if intersection and not intersection.is_empty:
                            if isinstance(intersection, Polygon):
                                split_buffered_lines.append(
                                    {
                                        "geometry": intersection,
                                        "block_id": buff_block_id,
                                        "id": idx,
                                    }
                                )
                                idx += 1
                            elif isinstance(intersection, MultiPolygon):
                                for poly in intersection.geoms:
                                    if poly and not poly.is_empty:
                                        split_buffered_lines.append(
                                            {
                                                "geometry": poly,
                                                "block_id": buff_block_id,
                                                "id": idx,
                                            }
                                        )
                                        idx += 1
                    except Exception as e:
                        print(
                            f"    Warning: Failed to intersect buffered line "
                            f"with block {buff_block_id}: {e}"
                        )
                        continue

            # Write split buffered lines
            split_buffered_layer_name = f"{blocks_layer_name}_on_grid"
            if split_buffered_lines:
                gdf_split = gpd.GeoDataFrame(split_buffered_lines, geometry="geometry", crs=crs)
                gdf_split["type"] = "loc_loc"
                gdf_split.to_file(output_gpkg, layer=split_buffered_layer_name, driver="GPKG")
                print(f"  Split buffered lines layer: {split_buffered_layer_name}")

            # Create off-grid layer
            off_grid_layer_name = f"{blocks_layer_name}_off_grid"
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

                remaining_geom = block_geom.buffer(0).difference(
                    intersecting_buffer.buffer(0.00001)
                )

                if remaining_geom and not remaining_geom.is_empty:
                    _type = (
                        ClusterTypes.CONCAVE_CORNER
                        if blk["is_concave"]
                        else ClusterTypes.OFF_GRID_COLD
                    )
                    remaining_features.append(
                        {
                            "geometry": remaining_geom,
                            "orig_id": orig_block_id,
                            "is_concave": blk["is_concave"],
                            "concave_x": blk["concave_x"],
                            "concave_y": blk["concave_y"],
                            "type": _type,
                            "block_type": blk["block_type"],
                            "color": ColorTypes[_type],
                        }
                    )

            if remaining_features:
                gdf_remaining = gpd.GeoDataFrame(remaining_features, geometry="geometry", crs=crs)
                gdf_remaining.reset_index(drop=True, inplace=True)
                gdf_remaining["id"] = gdf_remaining.index
                gdf_remaining.to_file(output_gpkg, layer=off_grid_layer_name, driver="GPKG")
                print(f"  Created off-grid cold boundary layer: {off_grid_layer_name}")

        except Exception as e:
            print(f"  Warning: Could not process buffered lines layer: {e}")

    total_lines = len(cutting_lines_to_write)
    print(f"  Created {total_lines} cutting lines")
    print(f"  Created layer: {output_layer_name}")
    print(f"  Created subdivided blocks layer: {blocks_layer_name}")
    return output_layer_name
