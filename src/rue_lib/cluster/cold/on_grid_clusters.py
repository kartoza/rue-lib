import geopandas as gpd
from shapely.ops import unary_union


def merge_and_classify_on_grid_clusters(
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

    if gdf_blocks.empty or gdf_convex.empty:
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

    print(f"  Grouped convex clusters by {len(convex_by_block)} blocks")

    result_polygons = []
    result_records = []
    cluster_id = 0

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

                    result_polygons.append(part)
                    result_records.append(
                        {
                            "id": cluster_id,
                            "block_id": int(block_id),
                            "type": "loc",
                            "area": float(part.area),
                            "block_type": "cold",
                            "on_grid": 1,
                        }
                    )
                    cluster_id += 1

                print(f"    Block {block_id}: created {len(parts)} loc_loc cluster(s)")

            except Exception as e:
                print(f"    Warning: Failed to subtract from block {block_id}: {e}")
                result_polygons.append(block_geom)
                result_records.append(
                    {
                        "id": cluster_id,
                        "block_id": int(block_id),
                        "type": "loc_loc",
                        "area": float(block_geom.area),
                        "block_type": "cold",
                        "on_grid": 1,
                    }
                )
                cluster_id += 1
        else:
            result_polygons.append(block_geom)
            result_records.append(
                {
                    "id": cluster_id,
                    "block_id": int(block_id),
                    "type": "loc_loc",
                    "area": float(block_geom.area),
                    "block_type": "cold",
                    "on_grid": 1,
                }
            )
            cluster_id += 1

    if not result_polygons:
        print("  No loc_loc clusters created")
        return output_layer_name

    # Update records id
    cluster_id = 0
    for record in result_records:
        record["id"] = cluster_id
        cluster_id += 1

    gdf_result = gpd.GeoDataFrame(result_records, geometry=result_polygons, crs=gdf_blocks.crs)
    gdf_result.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
    print(f"  Created {len(gdf_result)} loc_loc clusters")

    return output_layer_name
