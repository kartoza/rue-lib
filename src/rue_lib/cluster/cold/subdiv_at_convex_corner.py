import math

import geopandas as gpd
from shapely.geometry import Point


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
        boundary_pts = sorted(boundary_pts, key=lambda x: x["id"])
        for point in boundary_pts:
            pt = point["coords"]

            if idx == 0:
                idx += 1
                continue

            if idx == len(boundary_pts) - 1:
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
