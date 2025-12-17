import geopandas as gpd
import pandas as pd
from shapely.ops import unary_union

from rue_lib.core.definitions import ClusterTypes, ColorTypes


def merge_and_classify_on_grid_clusters(
    input_gpkg: str,
    blocks_layer_name: str,
    convex_clusters_layer_name: str,
    concave_points_layer_name: str,
    roads_buffer_layer_name: str,
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
    gdf_concave_points = gpd.read_file(input_gpkg, layer=concave_points_layer_name)
    gdf_roads_buffer = gpd.read_file(input_gpkg, layer=roads_buffer_layer_name)

    print(f"  Loaded {len(gdf_blocks)} on-grid blocks")
    print(f"  Loaded {len(gdf_convex)} convex corner clusters")

    if gdf_blocks.empty or gdf_convex.empty:
        gdf_result = gpd.GeoDataFrame([], geometry=[], crs=gdf_blocks.crs)
        gdf_result.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
        return output_layer_name

    convex_by_block = {}
    convex_records_by_block = {}
    for _, row in gdf_convex.iterrows():
        block_id = row.get("block_id")
        if block_id is not None:
            if block_id not in convex_by_block:
                convex_by_block[block_id] = []
            convex_by_block[block_id].append(row.geometry)
            convex_records_by_block.setdefault(block_id, []).append(row.to_dict())

    concave_points = []
    for _, row in gdf_concave_points.iterrows():
        concave_points.append(row.geometry)

    print(f"  Grouped convex clusters by {len(convex_by_block)} blocks")

    result_polygons = []
    result_records = []
    cluster_id = 0

    def _is_concave_point_in_geom(concave_points, geom):
        for point in concave_points:
            if geom.intersects(point.buffer(1)):
                return True
        return False

    for _, block_row in gdf_blocks.iterrows():
        block_id = block_row.get("id")
        block_geom = block_row.geometry

        if block_geom is None or block_geom.is_empty:
            continue

        if block_id in convex_by_block:
            convex_geoms = convex_by_block[block_id]
            convex_union = unary_union(convex_geoms)

            try:
                remaining = block_geom.difference(convex_union.buffer(0.000001))

                if remaining.is_empty or remaining.area <= 0:
                    print(f"    Block {block_id}: completely removed by convex clusters")
                    continue

                if remaining.geom_type == "Polygon":
                    parts = [remaining]
                elif remaining.geom_type == "MultiPolygon":
                    parts = list(remaining.geoms)
                else:
                    print(f"    Block {block_id}: unexpected geometry type {remaining.geom_type}")
                    continue

                result_polygons.extend(convex_geoms)
                result_records.extend(convex_records_by_block[block_id])

                for part in parts:
                    if part.is_empty or part.area <= 10.0:
                        continue
                    is_concave = _is_concave_point_in_geom(concave_points, part)
                    _type = is_concave and ClusterTypes.ON_GRID_LOC_LOC or ClusterTypes.ON_GRID_LOC

                    result_polygons.append(part)
                    result_records.append(
                        {
                            "id": cluster_id,
                            "block_id": int(block_id),
                            "type": _type,
                            "color": ColorTypes[_type],
                            "area": float(part.area),
                            "block_type": "cold",
                            "on_grid": 1,
                        }
                    )
                    cluster_id += 1

                print(f"    Block {block_id}: created {len(parts)} loc_loc cluster(s)")

            except Exception as e:
                print(f"    Warning: Failed to subtract from block {block_id}: {e}")
        else:
            result_polygons.append(block_geom)
            is_concave = _is_concave_point_in_geom(concave_points, block_geom)
            _type = is_concave and ClusterTypes.ON_GRID_LOC_LOC or ClusterTypes.ON_GRID_LOC
            result_records.append(
                {
                    "id": cluster_id,
                    "block_id": int(block_id),
                    "type": _type,
                    "color": ColorTypes[_type],
                    "area": float(block_geom.area),
                    "block_type": "cold",
                    "on_grid": 1,
                }
            )
            cluster_id += 1

    if not result_polygons:
        print("  No loc_loc clusters created")
        return output_layer_name

    def _get_road_type_priority(road_types, remove=False):
        """Get road type in priority order: arterial > secondary > local."""
        for road_type, short_name in [
            ("road_arterial", "art"),
            ("road_secondary", "sec"),
            ("road_local", "loc"),
        ]:
            if road_type in road_types:
                if remove:
                    road_types.remove(road_type)
                return short_name
        return "loc"

    cluster_id = 0
    for record in result_records:
        record["id"] = cluster_id
        road_types = []

        roads_touched = gdf_roads_buffer[
            gdf_roads_buffer.intersects(result_polygons[cluster_id].buffer(0.1))
        ]
        if not roads_touched.empty:
            road_types = roads_touched["type"].unique().tolist()
            road_types = [rt for rt in road_types if rt]

        is_corner = "_" in record.get("type")

        if is_corner:
            primary_type = _get_road_type_priority(road_types, remove=True)
            secondary_type = _get_road_type_priority(road_types, remove=False)
            record["type"] = f"{primary_type}_{secondary_type}"
        else:
            if "road_arterial" in road_types and "road_local" in road_types:
                record["type"] = "loc"
            else:
                record["type"] = _get_road_type_priority(road_types, remove=False)

        record["road_types"] = ",".join(road_types)

        cluster_id += 1

    gdf_result = gpd.GeoDataFrame(result_records, geometry=result_polygons, crs=gdf_blocks.crs)
    gdf_result.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
    print(f"  Created {len(gdf_result)} loc_loc clusters")

    return output_layer_name


def merge_and_classify_off_grid_clusters(
    input_gpkg: str,
    off_grid_blocks_layer_name: str,
    cold_clusters_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> str:
    """
    Merge and classify off-grid cold clusters.

    Args:
        input_gpkg: Path to input GeoPackage
        off_grid_blocks_layer_name: Name of the off-grid blocks layer
        cold_clusters_layer_name: Name of the cold clusters layer
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output layer with merged clusters

    Returns:
        Name of the output layer
    """
    gdf_blocks = gpd.read_file(input_gpkg, layer=off_grid_blocks_layer_name)
    gdf_cold = gpd.read_file(input_gpkg, layer=cold_clusters_layer_name)

    print(f"  Loaded {len(gdf_blocks)} off-grid blocks")
    print(f"  Loaded {len(gdf_cold)} cold clusters")

    if gdf_blocks.empty or gdf_cold.empty:
        return output_layer_name

    cold_by_orig_id = {}
    for _, row in gdf_cold.iterrows():
        orig_id = row.get("orig_id")
        if orig_id is not None:
            if orig_id not in cold_by_orig_id:
                cold_by_orig_id[orig_id] = []
            cold_by_orig_id[orig_id].append(row)

    print(f"  Grouped cold clusters by {len(cold_by_orig_id)} orig_ids")

    result_polygons = []
    result_records = []
    cluster_id = 0

    for _, block_row in gdf_blocks.iterrows():
        block_id = block_row.get("id")
        block_geom = block_row.geometry

        if block_geom is None or block_geom.is_empty:
            continue

        if block_id in cold_by_orig_id:
            for cold_cluster in cold_by_orig_id[block_id]:
                result_polygons.append(cold_cluster.geometry)
                cold_cluster_dict = cold_cluster.to_dict()
                cold_cluster_dict["id"] = cluster_id
                cold_cluster_dict["block_id"] = int(block_id)
                cold_cluster_dict["on_grid"] = 0
                cold_cluster_dict["block_type"] = "cold"
                result_records.append(cold_cluster_dict)
                cluster_id += 1
        else:
            result_polygons.append(block_geom)
            block_row_dict = block_row.to_dict()
            block_row_dict["id"] = cluster_id
            block_row_dict["block_id"] = int(block_id)
            block_row_dict["type"] = (
                ClusterTypes.CONCAVE_CORNER
                if block_row_dict["is_concave"]
                else ClusterTypes.OFF_GRID_COLD
            )
            block_row_dict["on_grid"] = 0
            block_row_dict["color"] = ColorTypes[block_row_dict["type"]]
            result_records.append(block_row_dict)
            cluster_id += 1

    if not result_polygons:
        print("  No loc_loc clusters created")
        return output_layer_name

    gdf_result = gpd.GeoDataFrame(result_records, geometry=result_polygons, crs=gdf_blocks.crs)
    gdf_result.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
    print(f"  Created {len(gdf_result)} loc_loc clusters from off-grid blocks")

    return output_layer_name


def merge_final_cold_clusters(
    input_gpkg: str,
    on_grid_clusters_layer_name: str,
    off_grid_clusters_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> str:
    """
    Merge final on-grid and off-grid cold clusters into a single layer.

    Args:
        input_gpkg: Path to input GeoPackage
        on_grid_clusters_layer_name: Name of the on-grid clusters layer
        off_grid_clusters_layer_name: Name of the off-grid clusters layer
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for merged output layer

    Returns:
        Name of the output layer
    """
    print("\nMerging final on-grid and off-grid cold clusters...")

    # Load both layers
    gdf_on_grid = gpd.read_file(input_gpkg, layer=on_grid_clusters_layer_name)
    gdf_off_grid = gpd.read_file(input_gpkg, layer=off_grid_clusters_layer_name)

    print(f"  Loaded {len(gdf_on_grid)} on-grid clusters")
    print(f"  Loaded {len(gdf_off_grid)} off-grid clusters")

    if gdf_on_grid.empty and gdf_off_grid.empty:
        print("  Warning: both input layers are empty")
        return output_layer_name

    if gdf_on_grid.empty:
        gdf_merged = gdf_off_grid.copy()
    elif gdf_off_grid.empty:
        gdf_merged = gdf_on_grid.copy()
    else:
        gdf_merged = gpd.GeoDataFrame(
            pd.concat([gdf_on_grid, gdf_off_grid], ignore_index=True),
            crs=gdf_on_grid.crs,
        )

    gdf_merged["id"] = range(len(gdf_merged))
    gdf_merged.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
    print(f"  Created {len(gdf_merged)} total cold clusters")

    return output_layer_name
