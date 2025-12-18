# src/rue_lib/public/open_space.py
import geopandas as gpd
import pandas as pd

from rue_lib.core.definitions import ClusterTypes, ColorTypes

MAX_OPEN_SPACE_OFF_GRID0_PER_GRID = 2


def _find_adjacent_blocks(
    central_off_grid0_block, blocks_in_site, cluster_type, min_shared_length: float = 1.0
):
    """
    Find blocks that are truly adjacent (sharing an edge) to the given central block.

    Diagonal touching (point intersection only) is not considered adjacent.

    Args:
        central_off_grid0_block: The central block row from GeoDataFrame
        blocks_in_site: GeoDataFrame of all blocks in the site
        cluster_type: Type of cluster to filter for
        min_shared_length: Minimum shared edge length to consider adjacent (default 1.0m)

    Returns:
        GeoDataFrame of adjacent blocks
    """

    central_geom = central_off_grid0_block.geometry
    filtered_blocks = blocks_in_site[blocks_in_site["type"] == cluster_type].copy()

    if filtered_blocks.empty:
        return filtered_blocks

    adjacent_indices = []

    for idx, block_row in filtered_blocks.iterrows():
        block_geom = block_row.geometry
        if not central_geom.intersects(block_geom):
            continue
        intersection = central_geom.intersection(block_geom)
        if intersection.geom_type == "LineString":
            shared_length = intersection.length
        elif intersection.geom_type == "MultiLineString":
            shared_length = sum(line.length for line in intersection.geoms)
        else:
            continue
        if shared_length >= min_shared_length:
            adjacent_indices.append(idx)

    return filtered_blocks.loc[adjacent_indices]


def _allocate_single_block(
    block,
    all_allocated_blocks,
    allocated_indices,
    current_area,
    group_id,
    site_id,
    group_rank,
    parcel_centroid,
    is_root=0,
):
    """
    Allocate a single block as open space with common metadata.

    Returns:
        Updated current_area
    """
    block_idx = block.name
    if block_idx in allocated_indices:
        return current_area

    block_data = block.copy()
    block_data["allocation"] = "open_space"
    block_data["group_id"] = group_id
    block_data["is_root"] = is_root
    block_data["site_id"] = site_id
    block_data["group_rank"] = group_rank
    block_data["dist_to_center"] = block["centroid"].distance(parcel_centroid)
    block_data["type"] = block_data["type"] + "_os"
    if "cluster_type" in block_data:
        block_data["cluster_type"] = block_data["cluster_type"] + "_os"
    block_data["color"] = ColorTypes[block_data["type"]]

    all_allocated_blocks.append(block_data)
    allocated_indices.add(block_idx)
    return current_area + block["block_area"]


def _allocate_adjacent_blocks_of_type(
    source_block,
    blocks_in_site,
    cluster_type,
    all_allocated_blocks,
    allocated_indices,
    current_area,
    required_open_area,
    group_id,
    site_id,
    group_rank,
    parcel_centroid,
):
    """
    Find and allocate adjacent blocks of a specific type.

    Returns:
        Tuple of (updated current_area, list of allocated blocks of this type)
    """
    adjacent_blocks = _find_adjacent_blocks(source_block, blocks_in_site, cluster_type)
    allocated_this_type = []

    for _, block in adjacent_blocks.iterrows():
        if current_area >= required_open_area:
            break

        block_idx = block.name
        if block_idx in allocated_indices:
            continue

        current_area = _allocate_single_block(
            block,
            all_allocated_blocks,
            allocated_indices,
            current_area,
            group_id,
            site_id,
            group_rank,
            parcel_centroid,
        )
        allocated_this_type.append(block)

    return current_area, allocated_this_type


def _allocate_with_cascade(
    source_block,
    blocks_in_site,
    all_allocated_blocks,
    allocated_indices,
    current_area,
    required_open_area,
    group_id,
    site_id,
    group_rank,
    parcel_centroid,
):
    """
    Allocate adjacent blocks with cascade: loc -> loc_loc pattern.

    Returns:
        Updated current_area
    """
    if current_area >= required_open_area:
        return current_area
    current_area, loc_blocks = _allocate_adjacent_blocks_of_type(
        source_block,
        blocks_in_site,
        ClusterTypes.ON_GRID_LOC,
        all_allocated_blocks,
        allocated_indices,
        current_area,
        required_open_area,
        group_id,
        site_id,
        group_rank,
        parcel_centroid,
    )
    if current_area >= required_open_area:
        return current_area

    # Allocate loc_loc blocks adjacent to allocated loc blocks
    for loc_block in loc_blocks:
        if current_area >= required_open_area:
            break

        current_area, _ = _allocate_adjacent_blocks_of_type(
            loc_block,
            blocks_in_site,
            ClusterTypes.ON_GRID_LOC_LOC,
            all_allocated_blocks,
            allocated_indices,
            current_area,
            required_open_area,
            group_id,
            site_id,
            group_rank,
            parcel_centroid,
        )

    return current_area


def allocate_open_spaces(
    output_gpkg: str,
    parcel_layer_name: str,
    block_layer_name: str,
    output_layer_name: str,
    open_percent: float = 4.0,
) -> str:
    """
    Allocate open spaces starting from central blocks.

    Implements Mobius logic:
    1. Find central block (closest to site centroid)
    2. Group blocks by off-grid roots with attached on-grid parts
    3. Allocate first 2 groups for open space
    4. Sort by distance to central block

    Args:
        output_gpkg: Path to output GeoPackage
        parcel_layer_name: Name for parcel layer
        block_layer_name: Name for blocks layer
        output_layer_name: Name for output open spaces layer
        open_percent: Percentage of site area to allocate as open space

    Returns:
        Name of the output layer
    """
    print(f"Allocating open spaces ({open_percent}% of site area)...")

    # Load input layers
    gdf_blocks = gpd.read_file(output_gpkg, layer=block_layer_name)
    gdf_parcels = gpd.read_file(output_gpkg, layer=parcel_layer_name)

    # Ensure we have an id column to track original feature IDs
    if "cluster_index" not in gdf_blocks.columns:
        gdf_blocks["cluster_index"] = range(len(gdf_blocks))

    all_allocated_blocks = []
    site_centroids = []
    site_id = 0

    # Process each parcel/site
    for _, parcel_row in gdf_parcels.iterrows():
        site_id += 1
        parcel_geom = parcel_row.geometry
        site_area = parcel_geom.area
        required_open_area = site_area * (open_percent / 100.0)

        print(f"  Site {site_id}: area={site_area:.2f}, required_open={required_open_area:.2f}")

        # Get all blocks that intersect this parcel
        blocks_in_site = gdf_blocks[gdf_blocks.intersects(parcel_geom)].copy()

        if blocks_in_site.empty:
            print(f"    No blocks found for site {site_id}")
            continue

        # Calculate centroids and distances
        parcel_centroid = parcel_geom.centroid

        # Store site centroid for debugging
        site_centroids.append(
            {
                "site_id": site_id,
                "geometry": parcel_centroid,
                "site_area": site_area,
                "required_open_area": required_open_area,
            }
        )

        blocks_in_site["centroid"] = blocks_in_site.geometry.centroid
        blocks_in_site["distance"] = blocks_in_site["centroid"].distance(parcel_centroid)
        blocks_in_site["block_area"] = blocks_in_site.geometry.area

        # Find off-grid0 (warm) blocks sorted by distance to parcel centroid
        off_grid0_blocks = blocks_in_site[
            blocks_in_site["type"] == ClusterTypes.OFF_GRID_WARM
        ].copy()

        # Initialize allocation tracking
        allocated_indices = set()
        block_ids = set()
        current_area = 0.0
        plot_index = None
        group_rank = 0

        if off_grid0_blocks.empty:
            print(f"    No off-grid0 blocks found for site {site_id}, will try concave corners")
        else:
            off_grid0_blocks = off_grid0_blocks.sort_values("distance")
            central_off_grid0_block = off_grid0_blocks.iloc[0]
            plot_index = central_off_grid0_block.get("plot_index")
            block_ids.add(central_off_grid0_block.get("block_id"))

            print(f"    Found {len(blocks_in_site)} blocks, central block identified")

            # Phase 1: Allocate central off-grid0 block and cascade to adjacent loc/loc_loc
            group_rank += 1
            current_area = _allocate_single_block(
                central_off_grid0_block,
                all_allocated_blocks,
                allocated_indices,
                current_area,
                plot_index,
                site_id,
                group_rank,
                parcel_centroid,
                is_root=1,
            )

            current_area = _allocate_with_cascade(
                central_off_grid0_block,
                blocks_in_site,
                all_allocated_blocks,
                allocated_indices,
                current_area,
                required_open_area,
                plot_index,
                site_id,
                group_rank,
                parcel_centroid,
            )

            # Phase 2: Allocate one adjacent off-grid0 if space remains and under limit
            group_rank += 1
            off_grid0_count = 1  # Already allocated central block
            if (
                current_area < required_open_area
                and off_grid0_count < MAX_OPEN_SPACE_OFF_GRID0_PER_GRID
            ):
                off_grid0_adjacents = _find_adjacent_blocks(
                    central_off_grid0_block, blocks_in_site, ClusterTypes.OFF_GRID_WARM
                )

                if not off_grid0_adjacents.empty:
                    off_grid0_adjacent = off_grid0_adjacents.iloc[0]
                    off_grid0_idx = off_grid0_adjacent.name

                    if off_grid0_idx not in allocated_indices:
                        off_grid0_count += 1
                        block_ids.add(off_grid0_adjacent.get("block_id"))

                        current_area = _allocate_single_block(
                            off_grid0_adjacent,
                            all_allocated_blocks,
                            allocated_indices,
                            current_area,
                            plot_index,
                            site_id,
                            group_rank,
                            parcel_centroid,
                        )

                        current_area = _allocate_with_cascade(
                            off_grid0_adjacent,
                            blocks_in_site,
                            all_allocated_blocks,
                            allocated_indices,
                            current_area,
                            required_open_area,
                            plot_index,
                            site_id,
                            group_rank,
                            parcel_centroid,
                        )

        # Phase 3: Allocate concave corner blocks with adjacent loc_loc if space remains
        group_rank += 1
        if current_area < required_open_area:
            print(f"    Phase {group_rank}: Finding concave corner blocks...")
            concave_corners = blocks_in_site[
                blocks_in_site["type"] == ClusterTypes.CONCAVE_CORNER
            ].copy()

            if not concave_corners.empty:
                concave_corners = concave_corners.sort_values("distance")

                for _, concave_block in concave_corners.iterrows():
                    if current_area >= required_open_area:
                        break

                    concave_idx = concave_block.name
                    if concave_idx in allocated_indices:
                        continue

                    if plot_index is None:
                        plot_index = concave_block.get("id")

                    current_area = _allocate_single_block(
                        concave_block,
                        all_allocated_blocks,
                        allocated_indices,
                        current_area,
                        plot_index,
                        site_id,
                        group_rank,
                        parcel_centroid,
                    )

                    loc_loc_adjacent = _find_adjacent_blocks(
                        concave_block, blocks_in_site, ClusterTypes.ON_GRID_LOC_LOC
                    )

                    if not loc_loc_adjacent.empty:
                        current_area = _allocate_single_block(
                            loc_loc_adjacent.iloc[0],
                            all_allocated_blocks,
                            allocated_indices,
                            current_area,
                            plot_index,
                            site_id,
                            group_rank,
                            parcel_centroid,
                        )

        # Phase 4: Allocate remaining off-grid0 blocks (furthest first) with one adjacent loc
        group_rank += 1
        if current_area < required_open_area and not off_grid0_blocks.empty:
            print(f"    Phase {group_rank}: Allocating remaining off-grid0 blocks...")

            while current_area < required_open_area:
                # Find remaining off-grid0 blocks not yet allocated and from different block_ids
                remaining_off_grid0 = off_grid0_blocks[
                    ~off_grid0_blocks.index.isin(allocated_indices)
                    & ~off_grid0_blocks["block_id"].isin(block_ids)
                ].copy()

                if remaining_off_grid0.empty:
                    break

                # Sort by distance (furthest first to distribute evenly)
                remaining_off_grid0 = remaining_off_grid0.sort_values("distance", ascending=False)
                remaining_block = remaining_off_grid0.iloc[0]
                block_ids.add(remaining_block.get("block_id"))

                # Allocate the off-grid0 block
                current_area = _allocate_single_block(
                    remaining_block,
                    all_allocated_blocks,
                    allocated_indices,
                    current_area,
                    plot_index,
                    site_id,
                    group_rank,
                    parcel_centroid,
                )

                # Try to allocate one adjacent loc block
                if current_area < required_open_area:
                    loc_adjacent = _find_adjacent_blocks(
                        remaining_block, blocks_in_site, ClusterTypes.ON_GRID_LOC
                    )

                    if not loc_adjacent.empty:
                        current_area = _allocate_single_block(
                            loc_adjacent.iloc[0],
                            all_allocated_blocks,
                            allocated_indices,
                            current_area,
                            plot_index,
                            site_id,
                            group_rank,
                            parcel_centroid,
                        )

        print(
            f"    Allocated {len(allocated_indices)} blocks as open space (area={current_area:.2f})"
        )

    # Create output GeoDataFrame
    if all_allocated_blocks:
        gdf_output = gpd.GeoDataFrame(all_allocated_blocks, crs=gdf_blocks.crs)
        cols_to_drop = ["centroid", "distance", "block_area"]
        gdf_output = gdf_output.drop(columns=[c for c in cols_to_drop if c in gdf_output.columns])
        gdf_output.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
        print(f"Open spaces saved to layer: {output_layer_name}")
    else:
        # Create empty layer with proper schema
        gdf_output = gdf_blocks.iloc[:0].copy()
        gdf_output["allocation"] = pd.Series(dtype="str")
        gdf_output["group_id"] = pd.Series(dtype="int")
        gdf_output["is_root"] = pd.Series(dtype="int")
        gdf_output["group_rank"] = pd.Series(dtype="int")
        gdf_output["dist_to_center"] = pd.Series(dtype="float")
        gdf_output["cluster_index"] = pd.Series(dtype="int")
        gdf_output.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
        print(f"No open spaces allocated, empty layer created: {output_layer_name}")

    # Save site centroids for debugging
    if site_centroids:
        gdf_centroids = gpd.GeoDataFrame(site_centroids, crs=gdf_blocks.crs)
        centroids_layer_name = f"{output_layer_name}_centroids"
        gdf_centroids.to_file(output_gpkg, layer=centroids_layer_name, driver="GPKG")
        print(f"Site centroids saved to layer: {centroids_layer_name}")

    return output_layer_name
