import math

import geopandas as gpd
from shapely.geometry import LineString, Point
from shapely.ops import split as shapely_split
from shapely.ops import unary_union


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

    convex_by_block = {}
    for _, row in gdf_convex.iterrows():
        block_id = row.get("block_id")
        if block_id is not None:
            if block_id not in convex_by_block:
                convex_by_block[block_id] = []
            convex_by_block[block_id].append(row)

    print(f"  Processing {len(convex_by_block)} blocks with convex points...")

    for block_id, convex_points in convex_by_block.items():
        if block_id not in blocks_by_id:
            print(f"    Warning: Block {block_id} not found")
            continue

        block_geom = blocks_by_id[block_id]
        split_lines = []

        for convex_row in convex_points:
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

                line0_extended = line0
                split_lines.append(line0_extended)

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

                line2_extended = line2
                split_lines.append(line2_extended)

        if split_lines:
            try:
                lines_union = unary_union(split_lines)
                split_result = shapely_split(block_geom, lines_union)

                if hasattr(split_result, "geoms"):
                    parts = list(split_result.geoms)
                else:
                    parts = [split_result]

                min_area = 10.0

                corner_parts_count = 0
                for part in parts:
                    if part.is_empty or part.area <= 0:
                        continue

                    part_area = float(part.area)

                    if part_area < min_area:
                        continue

                    cluster_polygons.append(part)
                    cluster_records.append(
                        {
                            "id": cluster_id,
                            "block_id": int(block_id),
                            "type": "convex_cluster",
                            "area": part_area,
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


def subtract_convex_clusters_from_blocks(
    input_gpkg: str,
    blocks_layer_name: str,
    convex_clusters_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> str:
    """
    Subtract convex corner clusters from on-grid blocks to create loc_loc clusters.

    Args:
        input_gpkg: Path to input GeoPackage
        blocks_layer_name: Name of the on-grid blocks layer
        convex_clusters_layer_name: Name of the convex corner clusters layer
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output layer with loc_loc clusters

    Returns:
        Name of the output layer
    """
    print("\nSubtracting convex clusters from blocks to create loc_loc clusters...")

    # Load layers
    gdf_blocks = gpd.read_file(input_gpkg, layer=blocks_layer_name)
    gdf_convex = gpd.read_file(input_gpkg, layer=convex_clusters_layer_name)

    print(f"  Loaded {len(gdf_blocks)} on-grid blocks")
    print(f"  Loaded {len(gdf_convex)} convex corner clusters")

    if gdf_blocks.empty:
        print("  Warning: blocks layer is empty")
        return output_layer_name

    if gdf_convex.empty:
        print("  No convex clusters to subtract, using original blocks")
        # Just copy blocks with loc_loc type
        gdf_blocks["type"] = "loc_loc"
        gdf_blocks.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
        return output_layer_name

    # Group convex clusters by block_id
    convex_by_block = {}
    for _, row in gdf_convex.iterrows():
        block_id = row.get("block_id")
        if block_id is not None:
            if block_id not in convex_by_block:
                convex_by_block[block_id] = []
            convex_by_block[block_id].append(row.geometry)

    print(f"  Grouped convex clusters by {len(convex_by_block)} blocks")

    # Subtract convex clusters from each block
    result_polygons = []
    result_records = []
    cluster_id = 0

    for _, block_row in gdf_blocks.iterrows():
        block_id = block_row.get("id")
        block_geom = block_row.geometry

        if block_geom is None or block_geom.is_empty:
            continue

        # Check if this block has convex clusters to subtract
        if block_id in convex_by_block:
            convex_geoms = convex_by_block[block_id]
            convex_union = unary_union(convex_geoms)

            try:
                # Subtract the convex clusters from the block
                remaining = block_geom.difference(convex_union)

                if remaining.is_empty or remaining.area <= 0:
                    print(f"    Block {block_id}: completely removed by convex clusters")
                    continue

                # Handle both single and multi-part results
                if remaining.geom_type == "Polygon":
                    parts = [remaining]
                elif remaining.geom_type == "MultiPolygon":
                    parts = list(remaining.geoms)
                else:
                    print(f"    Block {block_id}: unexpected geometry type {remaining.geom_type}")
                    continue

                # Create a loc_loc cluster for each part
                for part in parts:
                    if part.is_empty or part.area <= 10.0:  # Skip very small parts
                        continue

                    result_polygons.append(part)
                    result_records.append(
                        {
                            "id": cluster_id,
                            "block_id": int(block_id),
                            "type": "loc_loc",
                            "area": float(part.area),
                        }
                    )
                    cluster_id += 1

                print(f"    Block {block_id}: created {len(parts)} loc_loc cluster(s)")

            except Exception as e:
                print(f"    Warning: Failed to subtract from block {block_id}: {e}")
                # Keep the original block if subtraction fails
                result_polygons.append(block_geom)
                result_records.append(
                    {
                        "id": cluster_id,
                        "block_id": int(block_id),
                        "type": "loc_loc",
                        "area": float(block_geom.area),
                    }
                )
                cluster_id += 1
        else:
            # No convex clusters for this block, keep it as is
            result_polygons.append(block_geom)
            result_records.append(
                {
                    "id": cluster_id,
                    "block_id": int(block_id),
                    "type": "loc_loc",
                    "area": float(block_geom.area),
                }
            )
            cluster_id += 1

    if not result_polygons:
        print("  No loc_loc clusters created")
        return output_layer_name

    gdf_result = gpd.GeoDataFrame(result_records, geometry=result_polygons, crs=gdf_blocks.crs)
    gdf_result.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
    print(f"  Created {len(gdf_result)} loc_loc clusters")

    return output_layer_name
