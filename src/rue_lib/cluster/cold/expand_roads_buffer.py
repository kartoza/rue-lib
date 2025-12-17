import geopandas as gpd

from rue_lib.cluster.config import ClusterConfig


def create_buffered_lines_from_boundary_lines(
    input_gpkg: str,
    roads_layer_name: str,
    boundary_lines_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
    cfg: ClusterConfig,
) -> str:
    """
    Buffer road-adjacent boundary lines to create road strips.

    Uses the boundary lines layer (e.g. 202_cold_boundary_lines) instead of points
    to generate rectangular buffers per road type.
    """
    gdf_lines = gpd.read_file(input_gpkg, layer=boundary_lines_layer_name)
    gdf_roads = gpd.read_file(input_gpkg, layer=roads_layer_name)

    if gdf_lines.empty:
        raise ValueError(f"Layer {boundary_lines_layer_name} is empty")

    buffer_depths = {
        "road_local": cfg.on_grid_partition_depth_local_roads * 2,
        "road_sec": cfg.on_grid_partition_depth_secondary_roads * 2,
        "road_art": cfg.on_grid_partition_depth_arterial_roads * 2,
    }

    buffered_lines = []
    search_buffer = 1.0  # Small buffer to find intersecting roads

    for _, row in gdf_lines.iterrows():
        geom = row.geometry

        if geom is None or geom.is_empty or geom.length < 10.0:
            continue

        # Check road type by finding closest road in gdf_roads
        # Buffer the boundary line geometry slightly and find intersecting roads
        buffered_search = geom.centroid.buffer(search_buffer)

        # Find roads that intersect with the buffered boundary line
        road_type = None
        min_distance = float("inf")

        for _, road_row in gdf_roads.iterrows():
            road_geom = road_row.geometry
            if road_geom is None or road_geom.is_empty:
                continue

            # Check if the buffered boundary intersects with the road
            if buffered_search.intersects(road_geom):
                # Calculate distance to find the closest road
                distance = geom.distance(road_geom)
                if distance < min_distance:
                    min_distance = distance
                    road_type = road_row["road_type"]

        # If no road found by intersection, skip this boundary line
        if road_type is None:
            road_type = "road_local"  # Default to local if no road found

        depth = buffer_depths.get(road_type, 0)
        if depth <= 0:
            continue

        # Buffer with flat caps (cap_style=2) - only sides, no extension at endpoints
        # join_style=2 is mitered joins for sharp corners
        # buffered_geom = geom.buffer(depth / 2.0, join_style=2, cap_style=3)
        buffered_geom = geom

        if buffered_geom and not buffered_geom.is_empty:
            buffered_lines.append(
                {
                    "geometry": geom,
                    "road_type": road_type,
                    "lenght_m": geom.length,
                    "buffer_depth_m": depth,
                    "block_id": row.get("block_id", None),
                }
            )

    if buffered_lines:
        gdf_buffered = gpd.GeoDataFrame(buffered_lines, geometry="geometry", crs=gdf_lines.crs)
    else:
        gdf_buffered = gpd.GeoDataFrame([], geometry=[], crs=gdf_lines.crs)

    gdf_buffered.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")

    print(f"  Created buffered lines layer: {output_layer_name}")
    return output_layer_name


def clip_buffered_lines_to_cold_grid(
    input_gpkg: str,
    buffered_lines_layer_name: str,
    cold_grid_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> str:
    """
    Clip buffered lines to only keep geometries inside cold grid blocks.

    Args:
        input_gpkg: Path to input GeoPackage
        buffered_lines_layer_name: Name of buffered lines layer (e.g., 206_buffered_lines)
        cold_grid_layer_name: Name of cold grid erased layer (e.g., 201_cold_grid_erased)
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output clipped layer

    Returns:
        Name of the output layer
    """
    # Read input data using geopandas
    gdf_buffered = gpd.read_file(input_gpkg, layer=buffered_lines_layer_name)
    gdf_cold_grid = gpd.read_file(input_gpkg, layer=cold_grid_layer_name)

    if gdf_buffered.empty:
        raise ValueError(f"Layer {buffered_lines_layer_name} is empty")
    if gdf_cold_grid.empty:
        raise ValueError(f"Layer {cold_grid_layer_name} is empty")

    cold_grid_geoms = {}
    for idx, row in gdf_cold_grid.iterrows():
        block_id = row.get("block_id", idx + 1)
        cold_grid_geoms[block_id] = row.geometry

    clipped_geometries = []
    for _, row in gdf_buffered.iterrows():
        buffered_geom = row.geometry
        block_id = None
        road_type = row["road_type"]

        if buffered_geom is None or buffered_geom.is_empty:
            continue
        cold_grid_geom = cold_grid_geoms.get(block_id)
        if cold_grid_geom is None:
            continue
        try:
            clipped_geom = buffered_geom.intersection(cold_grid_geom)
            if clipped_geom and not clipped_geom.is_empty:
                clipped_geometries.append(
                    {
                        "geometry": clipped_geom,
                        "block_id": block_id,
                        "road_type": road_type,
                    }
                )
        except Exception as e:
            print(f"    Warning: Failed to clip geometry for block {block_id}: {e}")
            continue

    # Write output
    if clipped_geometries:
        gdf_clipped = gpd.GeoDataFrame(
            clipped_geometries, geometry="geometry", crs=gdf_buffered.crs
        )
    else:
        # Create empty layer with proper schema
        gdf_clipped = gpd.GeoDataFrame([], geometry=[], crs=gdf_buffered.crs)

    gdf_clipped.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")

    print(f"  Created clipped layer: {output_layer_name} ({len(clipped_geometries)} features)")
    return output_layer_name
