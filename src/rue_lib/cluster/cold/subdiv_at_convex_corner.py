import math

import geopandas as gpd
from shapely.geometry import LineString, Point
from shapely.ops import split as shapely_split
from shapely.ops import unary_union

from rue_lib.core.definitions import ClusterTypes, ColorTypes


def find_convex_points(
    input_gpkg: str,
    points_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> str:
    """
    Find convex points from the boundary vertices.
    """
    gdf_points = gpd.read_file(input_gpkg, layer=points_layer_name)

    # Map boundary points by block_id for fast lookup
    points_by_block: dict[int, list[tuple[float, float]]] = {}
    for _, brow in gdf_points.iterrows():
        blk = brow.get("orig_id")
        geom = brow.geometry
        if blk is None or geom is None or not isinstance(geom, Point):
            continue
        points_by_block.setdefault(int(blk), []).append(
            {
                "id": brow.get("vertex_idx"),
                "coords": (geom.x, geom.y),
                "angle": brow.get("angle_deg"),
                "line_id": brow.get("line_id"),
            }
        )

    records = []
    geoms = []

    for block_id, boundary_pts in points_by_block.items():
        idx = 0
        for point in boundary_pts:
            pt = point["coords"]

            if idx == 0:
                idx += 1
                continue

            if idx == len(boundary_pts) - 1:
                idx += 1
                continue

            p0 = boundary_pts[idx - 1]["coords"]
            p1 = pt
            p2 = boundary_pts[idx + 1]["coords"]

            p0_dist = math.dist(p0, p1)
            p2_dist = math.dist(p1, p2)

            ang = point["angle"]

            if ang > 60 and ang < 100:
                records.append(
                    {
                        "is_convex": 1,
                        "block_id": int(block_id),
                        "vertex_id": point["id"],
                        "angle": ang,
                        "p0_dist": p0_dist,
                        "p2_dist": p2_dist,
                        "p0_x": p0[0],
                        "p0_y": p0[1],
                        "p2_x": p2[0],
                        "p2_y": p2[1],
                    }
                )
                geoms.append(Point(pt))
            idx += 1

    if not geoms:
        print("  No convex points found")
        return output_layer_name

    gdf_out = gpd.GeoDataFrame(records, geometry=geoms, crs=gdf_points.crs)
    gdf_out.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
    print(f"  Found {len(gdf_out)} convex points")
    return output_layer_name


def create_clusters_from_convex_points(
    input_gpkg: str,
    convex_points_layer_name: str,
    blocks_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
    max_partition_depth: float,
) -> str:
    """
    Create cluster polygons from convex points by splitting blocks with perpendicular lines.

    Args:
        input_gpkg: Path to input GeoPackage
        convex_points_layer_name: Name of convex points layer
        blocks_layer_name: Name of blocks layer
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output clusters layer
        max_partition_depth: Maximum partition depth value

    Returns:
        Name of the output layer
    """
    print("\nCreating clusters from convex points...")

    gdf_convex = gpd.read_file(input_gpkg, layer=convex_points_layer_name)
    gdf_blocks = gpd.read_file(input_gpkg, layer=blocks_layer_name)

    print(f"  Loaded {len(gdf_convex)} convex points")
    print(f"  Loaded {len(gdf_blocks)} blocks")

    if gdf_convex.empty or gdf_blocks.empty:
        print("  Warning: one or more input layers are empty")
        return output_layer_name

    blocks_by_id = {}
    for _, row in gdf_blocks.iterrows():
        block_id = row.get("id")
        if block_id is not None:
            blocks_by_id[int(block_id)] = row.geometry

    line_length = max_partition_depth * 1.25

    cluster_polygons = []
    cluster_records = []
    cluster_id = 0
    split_lines_by_convex = {}

    convex_by_block = {}
    for _, row in gdf_convex.iterrows():
        block_id = row.get("block_id")
        if block_id is not None:
            if block_id not in convex_by_block:
                convex_by_block[block_id] = []
            convex_by_block[block_id].append(row)

    print(f"  Processing {len(convex_by_block)} blocks with convex points...")

    convex_idx = -1
    for block_id, convex_points in convex_by_block.items():
        if block_id not in blocks_by_id:
            print(f"    Warning: Block {block_id} not found")
            continue

        block_geom = blocks_by_id[block_id]
        for convex_row in convex_points:
            convex_idx += 1
            curr_point = convex_row.geometry
            cx, cy = curr_point.x, curr_point.y

            p0_x = convex_row.get("p0_x")
            p0_y = convex_row.get("p0_y")
            p2_x = convex_row.get("p2_x")
            p2_y = convex_row.get("p2_y")

            center_x0 = (cx + p0_x) / 2
            center_y0 = (cy + p0_y) / 2

            dx0 = p0_x - cx
            dy0 = p0_y - cy
            len0 = math.hypot(dx0, dy0)
            if len0 > 0:
                dx0, dy0 = dx0 / len0, dy0 / len0
                perp_x0, perp_y0 = -dy0, dx0

                test_dist = 1
                end1_x = center_x0 + perp_x0 * test_dist
                end1_y = center_y0 + perp_y0 * test_dist
                end2_x = center_x0 - perp_x0 * test_dist
                end2_y = center_y0 - perp_y0 * test_dist

                test_pt1 = Point(end1_x, end1_y)
                test_pt2 = Point(end2_x, end2_y)

                if block_geom.contains(test_pt1):
                    line_end_x = cx + perp_x0 * line_length
                    line_end_y = cy + perp_y0 * line_length
                    line0 = LineString([(cx, cy), (line_end_x, line_end_y)])
                elif block_geom.contains(test_pt2):
                    line_end_x = cx - perp_x0 * line_length
                    line_end_y = cy - perp_y0 * line_length
                    line0 = LineString([(cx, cy), (line_end_x, line_end_y)])
                else:
                    continue
                split_lines_by_convex[convex_idx] = {
                    "line0": line0,
                    "convex": convex_row,
                    "block_id": block_id,
                }
            center_x2 = (cx + p2_x) / 2
            center_y2 = (cy + p2_y) / 2
            dx2 = p2_x - cx
            dy2 = p2_y - cy
            len2 = math.hypot(dx2, dy2)
            if len2 > 0:
                dx2, dy2 = dx2 / len2, dy2 / len2
                perp_x2, perp_y2 = -dy2, dx2
                test_dist = 1
                end1_x = center_x2 + perp_x2 * test_dist
                end1_y = center_y2 + perp_y2 * test_dist
                end2_x = center_x2 - perp_x2 * test_dist
                end2_y = center_y2 - perp_y2 * test_dist
                test_pt1 = Point(end1_x, end1_y)
                test_pt2 = Point(end2_x, end2_y)
                if block_geom.contains(test_pt1):
                    line_end_x = cx + perp_x2 * line_length
                    line_end_y = cy + perp_y2 * line_length
                    line2 = LineString([(cx, cy), (line_end_x, line_end_y)])
                elif block_geom.contains(test_pt2):
                    line_end_x = cx - perp_x2 * line_length
                    line_end_y = cy - perp_y2 * line_length
                    line2 = LineString([(cx, cy), (line_end_x, line_end_y)])
                else:
                    continue

                if convex_idx not in split_lines_by_convex:
                    continue

                split_lines_by_convex[convex_idx]["line2"] = line2

    for split_lines in split_lines_by_convex:
        split_lines_data = split_lines_by_convex[split_lines]
        block_geom = blocks_by_id[split_lines_data["block_id"]]
        try:
            lines_union = unary_union([split_lines_data["line0"], split_lines_data["line2"]])
            split_result = shapely_split(block_geom, lines_union)

            if hasattr(split_result, "geoms"):
                parts = list(split_result.geoms)
            else:
                parts = [split_result]

            min_area = 0.20 * (max_partition_depth**2)

            corner_parts_count = 0

            parts.sort(key=lambda g: g.centroid.distance(split_lines_data["convex"].geometry))

            if len(parts) > 2:
                part = parts[0]
                if part.is_empty or part.area <= 0:
                    continue

                part_area = float(part.area)

                if part_area < min_area:
                    continue

                cluster_polygons.append(part)
                cluster_records.append(
                    {
                        "id": cluster_id,
                        "block_id": int(split_lines_data["block_id"]),
                        "type": ClusterTypes.ON_GRID_LOC_LOC,
                        "color": ColorTypes[ClusterTypes.ON_GRID_LOC_LOC],
                        "on_grid": 1,
                        "block_type": "cold",
                        "area": part_area,
                        "distance_to_corner": float(
                            part.centroid.distance(split_lines_data["convex"].geometry)
                        ),
                    }
                )
                cluster_id += 1
                corner_parts_count += 1
            print(f"    Block {block_id}: created {corner_parts_count} corner clusters")

        except Exception as e:
            print(f"    Warning: Failed to split block {block_id}: {e}")

    if not cluster_polygons:
        print("  No clusters created")
        return output_layer_name

    gdf_clusters = gpd.GeoDataFrame(cluster_records, geometry=cluster_polygons, crs=gdf_blocks.crs)

    gdf_clusters.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
    print(f"  Created {len(gdf_clusters)} clusters from {len(gdf_convex)} convex points")

    return output_layer_name
