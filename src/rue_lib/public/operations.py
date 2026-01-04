# src/rue_lib/public/operations.py

import geopandas as gpd
import pandas as pd
from osgeo import ogr

from rue_lib.core.definitions import BlockTypes, ClusterTypes, ColorTypes
from rue_lib.public.open_space import find_adjacent_blocks


def allocate_cluster_index(
    output_gpkg: str,
    layer_name: str,
    *,
    id_field: str = "cluster_index",
    start_at: int = 0,
) -> str:
    """Ensure a layer has a unique cluster index column.

    This is used to carry a stable identifier through downstream public-space
    allocation steps (open spaces, amenities) so features can be merged back
    into a final layer.

    Args:
        output_gpkg: Path to the GeoPackage.
        layer_name: Layer name to update in-place.
        id_field: Column name to use for the unique identifier.
        start_at: Starting integer value.

    Returns:
        The updated layer name.
    """
    # Prefer reading with pyogrio (if available) so we can use the layer FID as
    # a stable fallback identifier when no explicit unique `id` column exists.
    gdf = gpd.read_file(output_gpkg, layer=layer_name)
    if gdf.empty:
        gdf[id_field] = pd.Series(dtype="int64")
        gdf.to_file(output_gpkg, layer=layer_name, driver="GPKG")
        return layer_name

    if "id" in gdf.columns and gdf["id"].notna().all() and gdf["id"].is_unique:
        gdf[id_field] = gdf["id"].astype(int)
    else:
        needs_reassign = (
            id_field not in gdf.columns or gdf[id_field].isna().any() or not gdf[id_field].is_unique
        )
        if needs_reassign:
            if getattr(gdf.index, "name", None) == "fid":
                fid_min = int(gdf.index.min())
                gdf[id_field] = (gdf.index.astype("int64") - fid_min + start_at).astype("int64")
            else:
                gdf[id_field] = pd.Series(range(start_at, start_at + len(gdf)), dtype="int64")
        else:
            gdf[id_field] = gdf[id_field].astype(int)

    gdf.to_file(output_gpkg, layer=layer_name, driver="GPKG")
    return layer_name


def _allocate_single_amenity(
    row,
    *,
    all_rows: list,
    allocated: set,
    current_area: float,
    group_id: int,
    site_id: int,
    group_rank: int,
    parcel_centroid,
    is_root: int = 0,
) -> float:
    """Allocate a single cluster as amenity and return updated area."""
    idx = row.name
    if idx in allocated:
        return current_area

    out = row.copy()
    out["allocation"] = "amenity"
    out["group_id"] = group_id
    out["is_root"] = is_root
    out["site_id"] = site_id
    out["group_rank"] = group_rank
    out["dist_to_center"] = row["centroid"].distance(parcel_centroid)

    out["type"] = f"{out['type']}_am"
    if "cluster_type" in out:
        out["cluster_type"] = f"{out['cluster_type']}_am"
    try:
        out["color"] = ColorTypes[out["type"]]
    except KeyError:
        # Keep the existing color if we don't have a mapping for this type.
        pass

    all_rows.append(out)
    allocated.add(idx)
    return current_area + row["block_area"]


def _order_block_ids_by_distance(
    block_ids: list[int], blocks_in_site_without_os: gpd.GeoDataFrame
) -> list[int]:
    """Order block_ids by the distance of their closest off_grid0 (fallback: closest any)."""
    scored: list[tuple[float, int]] = []
    for bid in block_ids:
        subset = blocks_in_site_without_os[blocks_in_site_without_os["block_id"] == bid]
        off0 = subset[subset["type"] == ClusterTypes.OFF_GRID_WARM]
        if not off0.empty:
            score = float(off0["distance"].min())
        else:
            score = float(subset["distance"].min())
        scored.append((score, bid))
    scored.sort(key=lambda t: t[0])
    return [bid for _, bid in scored]


def _allocate_in_block(
    block_id: int,
    cluster_type_start: ClusterTypes,
    *,
    blocks_in_site_without_os: gpd.GeoDataFrame,
    allocated: set,
    all_allocated_blocks: list,
    current_area: float,
    group_rank: int,
    required_amen_area: float,
    site_id: int,
    parcel_centroid,
) -> tuple[float, int]:
    """Allocate amenities within a specific block.

    Returns:
        Tuple of (updated_current_area, updated_group_rank)
    """
    in_block = blocks_in_site_without_os[blocks_in_site_without_os["block_id"] == block_id].copy()
    if in_block.empty:
        return current_area, group_rank

    budget_area = 0.3 * sum(in_block["block_area"])
    allocated_area_in_block = sum(in_block.loc[in_block.index.isin(allocated), "block_area"])

    allocated_in_block = set(in_block.index).intersection(allocated)

    if allocated_area_in_block > budget_area:
        return current_area, group_rank

    off0_candidates = in_block[
        (in_block["type"] == cluster_type_start) & ~in_block.index.isin(allocated)
    ].copy()
    if off0_candidates.empty:
        return current_area, group_rank

    group_rank += 1
    group_id = int(block_id) if block_id != -1 else group_rank

    off0_candidates = off0_candidates.sort_values("distance")
    seed_off_grid0 = off0_candidates.iloc[0]
    seed_loc = None

    if current_area < required_amen_area and allocated_area_in_block < budget_area:
        current_area = _allocate_single_amenity(
            seed_off_grid0,
            all_rows=all_allocated_blocks,
            allocated=allocated,
            current_area=current_area,
            group_id=group_id,
            site_id=site_id,
            group_rank=group_rank,
            parcel_centroid=parcel_centroid,
            is_root=1,
        )
        allocated_area_in_block += seed_off_grid0["block_area"]

    if current_area < required_amen_area and allocated_area_in_block < budget_area:
        loc_adj = find_adjacent_blocks(seed_off_grid0, in_block, ClusterTypes.ON_GRID_LOC)
        loc_adj = loc_adj[~loc_adj.index.isin(allocated)]
        if not loc_adj.empty:
            loc_adj = loc_adj.copy()
            loc_adj["d_seed"] = loc_adj["centroid"].distance(seed_off_grid0["centroid"])
            seed_loc = loc_adj.sort_values("d_seed").iloc[0]
            current_area = _allocate_single_amenity(
                seed_loc,
                all_rows=all_allocated_blocks,
                allocated=allocated,
                current_area=current_area,
                group_id=group_id,
                site_id=site_id,
                group_rank=group_rank,
                parcel_centroid=parcel_centroid,
            )
            allocated_area_in_block += seed_loc["block_area"]

    if current_area < required_amen_area and allocated_area_in_block < budget_area:
        off0_adj = find_adjacent_blocks(seed_off_grid0, in_block, cluster_type_start)
        off0_adj = off0_adj[~off0_adj.index.isin(allocated)]
        if not off0_adj.empty:
            current_area = _allocate_single_amenity(
                off0_adj.iloc[0],
                all_rows=all_allocated_blocks,
                allocated=allocated,
                current_area=current_area,
                group_id=group_id,
                site_id=site_id,
                group_rank=group_rank,
                parcel_centroid=parcel_centroid,
            )
            allocated_in_block.add(off0_adj.iloc[0].name)
            allocated_area_in_block += off0_adj.iloc[0]["block_area"]

    while current_area < required_amen_area and allocated_area_in_block < budget_area:
        next_off_grid0 = None
        for allocated_idx in list(allocated_in_block):
            adj = find_adjacent_blocks(in_block.loc[allocated_idx], in_block, cluster_type_start)
            adj = adj[~adj.index.isin(allocated)]
            if not adj.empty:
                next_off_grid0 = adj.sort_values("distance").iloc[0]
                break

        if next_off_grid0 is None:
            break

        if allocated_area_in_block >= budget_area:
            break

        current_area = _allocate_single_amenity(
            next_off_grid0,
            all_rows=all_allocated_blocks,
            allocated=allocated,
            current_area=current_area,
            group_id=group_id,
            site_id=site_id,
            group_rank=group_rank,
            parcel_centroid=parcel_centroid,
        )
        allocated_area_in_block += next_off_grid0["block_area"]

        if current_area >= required_amen_area or allocated_area_in_block >= budget_area:
            break

        locs = find_adjacent_blocks(next_off_grid0, in_block, ClusterTypes.ON_GRID_LOC)
        locs = locs[~locs.index.isin(allocated)]
        if locs.empty:
            continue

        if seed_loc is not None:
            locs = locs.copy()
            locs["d_loc"] = locs["centroid"].distance(seed_loc["centroid"])
            chosen_loc = locs.sort_values("d_loc", ascending=False).iloc[0]
        else:
            chosen_loc = locs.sort_values("distance").iloc[0]

        if allocated_area_in_block >= budget_area:
            break

        current_area = _allocate_single_amenity(
            chosen_loc,
            all_rows=all_allocated_blocks,
            allocated=allocated,
            current_area=current_area,
            group_id=group_id,
            site_id=site_id,
            group_rank=group_rank,
            parcel_centroid=parcel_centroid,
        )
        allocated_area_in_block += chosen_loc["block_area"]

    return current_area, group_rank


def allocate_amenities(
    output_gpkg: str,
    parcel_layer_name: str,
    block_layer_name: str,
    open_spaces_layer_name: str,
    output_layer_name: str,
    amen_percent: float = 10.0,
) -> str:
    """Allocate amenities from the remaining (non-open-space) clusters.

    Strategy (mirrors the intent of `allocate_open_spaces`, but for amenities):
    - For each parcel/site, allocate amenities by `block_id` groups.
    - Seed a block with: center `off_grid0` + exactly one adjacent `loc` (never `loc_loc`).
    - If still under the required area, add one adjacent `off_grid0` (no extra loc required).
    - If still under the required area, keep expanding inside that same `block_id` by:
      `off_grid0` adjacent to an already allocated cluster, then choose a `loc` adjacent to it
      that is furthest from the initial `loc` (to spread allocation).
    - Stop expanding that `block_id` after ~30% of its clusters are allocated, then move to
      another `block_id` (prefer blocks with no open space).
    """
    print(f"Allocating amenities ({amen_percent}% of site area)...")

    gdf_blocks = gpd.read_file(output_gpkg, layer=block_layer_name)
    gdf_parcels = gpd.read_file(output_gpkg, layer=parcel_layer_name)
    gdf_open_spaces = gpd.read_file(output_gpkg, layer=open_spaces_layer_name)

    if "cluster_index" not in gdf_blocks.columns:
        gdf_blocks["cluster_index"] = range(len(gdf_blocks))

    open_space_ids: set[int] = set()
    if "cluster_index" in gdf_open_spaces.columns:
        open_space_ids = set(gdf_open_spaces["cluster_index"].dropna().astype(int))

    all_allocated_blocks: list = []
    site_id = 0

    for _, parcel_row in gdf_parcels.iterrows():
        site_id += 1
        parcel_geom = parcel_row.geometry
        site_area = parcel_geom.area
        required_amen_area = site_area * (amen_percent / 100.0)

        print(f"  Site {site_id}: area={site_area:.2f}, required_amen={required_amen_area:.2f}")

        blocks_in_site = gdf_blocks[gdf_blocks.intersects(parcel_geom)].copy()
        warm_blocks = blocks_in_site[blocks_in_site["cluster_type"] != BlockTypes.COLD_GRID]
        cold_blocks = blocks_in_site[blocks_in_site["cluster_type"] == BlockTypes.COLD_GRID]
        blocks_in_site_without_os = blocks_in_site.copy()
        if open_space_ids:
            blocks_in_site_without_os = blocks_in_site[
                ~blocks_in_site["cluster_index"].isin(open_space_ids)
            ].copy()

        if blocks_in_site_without_os.empty:
            print(f"    No remaining blocks found for site {site_id}")
            continue

        print(f"    Remaining blocks found for site {site_id}: {len(blocks_in_site_without_os)}")
        parcel_centroid = parcel_geom.centroid
        blocks_in_site_without_os["centroid"] = blocks_in_site_without_os.geometry.centroid
        blocks_in_site_without_os["distance"] = blocks_in_site_without_os["centroid"].distance(
            parcel_centroid
        )
        blocks_in_site_without_os["block_area"] = blocks_in_site_without_os.geometry.area

        open_space_block_ids: set[int] = set()
        if "block_id" in gdf_open_spaces.columns and not gdf_open_spaces.empty:
            open_space_block_ids = set(
                gdf_open_spaces[gdf_open_spaces.intersects(parcel_geom)]["block_id"]
                .dropna()
                .astype(int)
            )

        allocated: set = set()
        group_rank = 0
        current_area = 0.0

        if "block_id" not in blocks_in_site_without_os.columns:
            raise ValueError(f"`{block_layer_name}` must contain a `block_id` column")

        processed_block_ids: set[int] = set()

        # Phase 1: allocate in the center block (even if it has open space).
        center_off0 = blocks_in_site_without_os[
            blocks_in_site_without_os["type"] == ClusterTypes.OFF_GRID_WARM
        ].copy()
        if not center_off0.empty:
            center_off0 = center_off0.sort_values("distance")
            center_block_id = int(center_off0.iloc[0]["block_id"])
            current_area, group_rank = _allocate_in_block(
                center_block_id,
                ClusterTypes.OFF_GRID_WARM,
                blocks_in_site_without_os=blocks_in_site_without_os,
                allocated=allocated,
                all_allocated_blocks=all_allocated_blocks,
                current_area=current_area,
                group_rank=group_rank,
                required_amen_area=required_amen_area,
                site_id=site_id,
                parcel_centroid=parcel_centroid,
            )
            processed_block_ids.add(center_block_id)

        # Phase 2: allocate for blocks WITHOUT open space.
        no_open_space_block_ids = sorted(
            set(warm_blocks["block_id"].dropna().astype(int)).difference(open_space_block_ids)
        )
        for bid in _order_block_ids_by_distance(no_open_space_block_ids, blocks_in_site_without_os):
            if current_area >= required_amen_area:
                break
            if bid in processed_block_ids:
                continue
            current_area, group_rank = _allocate_in_block(
                bid,
                ClusterTypes.OFF_GRID_WARM,
                blocks_in_site_without_os=blocks_in_site_without_os,
                allocated=allocated,
                all_allocated_blocks=all_allocated_blocks,
                current_area=current_area,
                group_rank=group_rank,
                required_amen_area=required_amen_area,
                site_id=site_id,
                parcel_centroid=parcel_centroid,
            )
            processed_block_ids.add(bid)

        no_open_space_block_ids = sorted(
            set(cold_blocks["block_id"].dropna().astype(int)).difference(open_space_block_ids)
        )
        for bid in _order_block_ids_by_distance(no_open_space_block_ids, blocks_in_site_without_os):
            if current_area >= required_amen_area:
                break
            if bid in processed_block_ids:
                continue
            current_area, group_rank = _allocate_in_block(
                bid,
                ClusterTypes.OFF_GRID_COLD,
                blocks_in_site_without_os=blocks_in_site_without_os,
                allocated=allocated,
                all_allocated_blocks=all_allocated_blocks,
                current_area=current_area,
                group_rank=group_rank,
                required_amen_area=required_amen_area,
                site_id=site_id,
                parcel_centroid=parcel_centroid,
            )
            processed_block_ids.add(bid)

        # Phase 3: if still short, allocate per-block for blocks WITH open space.
        with_open_space_block_ids = sorted(
            set(blocks_in_site["block_id"].dropna().astype(int)).intersection(open_space_block_ids)
        )
        for bid in _order_block_ids_by_distance(
            with_open_space_block_ids, blocks_in_site_without_os
        ):
            if current_area >= required_amen_area:
                break
            if bid in processed_block_ids:
                continue
            current_area, group_rank = _allocate_in_block(
                bid,
                ClusterTypes.OFF_GRID_WARM,
                blocks_in_site_without_os=blocks_in_site_without_os,
                allocated=allocated,
                all_allocated_blocks=all_allocated_blocks,
                current_area=current_area,
                group_rank=group_rank,
                required_amen_area=required_amen_area,
                site_id=site_id,
                parcel_centroid=parcel_centroid,
            )
            processed_block_ids.add(bid)

        print(f"    Allocated {len(allocated)} clusters as amenities (area={current_area:.2f})")

    if all_allocated_blocks:
        gdf_output = gpd.GeoDataFrame(all_allocated_blocks, crs=gdf_blocks.crs)
        cols_to_drop = ["centroid", "distance", "block_area"]
        gdf_output = gdf_output.drop(columns=[c for c in cols_to_drop if c in gdf_output.columns])
        gdf_output.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
        print(f"Amenities saved to layer: {output_layer_name}")
    else:
        gdf_output = gdf_blocks.iloc[:0].copy()
        gdf_output["allocation"] = pd.Series(dtype="str")
        gdf_output["group_id"] = pd.Series(dtype="int")
        gdf_output["is_root"] = pd.Series(dtype="int")
        gdf_output["site_id"] = pd.Series(dtype="int")
        gdf_output["group_rank"] = pd.Series(dtype="int")
        gdf_output["dist_to_center"] = pd.Series(dtype="float")
        if "cluster_index" not in gdf_output.columns:
            gdf_output["cluster_index"] = pd.Series(dtype="int")
        gdf_output.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
        print(f"No amenities allocated, empty layer created: {output_layer_name}")

    return output_layer_name


def merge(
    input_gpkg: str,
    final_layer_name: str,
    open_spaces_layer_name: str,
    amenities_layer_name: str,
):
    """
    Merge open spaces and amenities into the final layer.

    Replaces features in final_layer with matching block_id and plot_index
    from open_spaces and amenities layers.

    Args:
        input_gpkg: Path to GeoPackage containing all layers
        final_layer_name: Name of the final layer to update
        open_spaces_layer_name: Name of the open spaces layer
        amenities_layer_name: Name of the amenities layer
    """
    print(f"Merging layers into {final_layer_name}...")

    # Open the GeoPackage for reading and writing
    ds = ogr.Open(input_gpkg, 1)
    if ds is None:
        raise ValueError(f"Could not open {input_gpkg} for writing")

    final_layer = ds.GetLayerByName(final_layer_name)
    if final_layer is None:
        raise ValueError(f"Layer {final_layer_name} not found")

    open_spaces_layer = ds.GetLayerByName(open_spaces_layer_name)
    if open_spaces_layer is None:
        raise ValueError(f"Layer {open_spaces_layer_name} not found")

    amenities_layer = ds.GetLayerByName(amenities_layer_name)
    if amenities_layer is None:
        raise ValueError(f"Layer {amenities_layer_name} not found")

    # Process open spaces
    print(f"  Processing {open_spaces_layer_name}...")
    replaced_count = 0

    for open_space_feat in open_spaces_layer:
        cluster_index = open_space_feat.GetField("cluster_index")
        if cluster_index is None:
            continue
        try:
            cluster_index = int(cluster_index)
        except (TypeError, ValueError):
            continue

        # Find matching feature in final layer
        final_layer.ResetReading()
        for final_feat in final_layer:
            idx_field = final_layer.GetLayerDefn().GetFieldIndex("cluster_index")
            final_cluster_index = (
                final_feat.GetField("cluster_index") if idx_field != -1 else final_feat.GetFID()
            )
            try:
                final_cluster_index = int(final_cluster_index)
            except (TypeError, ValueError):
                continue

            if final_cluster_index == cluster_index:
                # Replace with open space feature
                fid = final_feat.GetFID()

                # Create updated feature
                updated_feat = ogr.Feature(final_layer.GetLayerDefn())
                updated_feat.SetFID(fid)
                updated_feat.SetGeometry(open_space_feat.GetGeometryRef())

                # Copy all fields from open space feature
                for i in range(open_space_feat.GetFieldCount()):
                    field_name = open_space_feat.GetFieldDefnRef(i).GetName()
                    if final_layer.GetLayerDefn().GetFieldIndex(field_name) >= 0:
                        value = open_space_feat.GetField(i)
                        updated_feat.SetField(field_name, value)

                final_layer.SetFeature(updated_feat)
                replaced_count += 1
                break

    print(f"    Replaced {replaced_count} features from open spaces")

    # Process amenities
    print(f"  Processing {amenities_layer_name}...")
    replaced_count = 0

    for amenity_feat in amenities_layer:
        cluster_index = amenity_feat.GetField("cluster_index")
        if cluster_index is None:
            continue
        try:
            cluster_index = int(cluster_index)
        except (TypeError, ValueError):
            continue

        # Find matching feature in final layer
        final_layer.ResetReading()
        for final_feat in final_layer:
            idx_field = final_layer.GetLayerDefn().GetFieldIndex("cluster_index")
            final_cluster_index = (
                final_feat.GetField("cluster_index") if idx_field != -1 else final_feat.GetFID()
            )
            try:
                final_cluster_index = int(final_cluster_index)
            except (TypeError, ValueError):
                continue

            if final_cluster_index == cluster_index:
                # Replace with amenity feature
                fid = final_feat.GetFID()

                # Create updated feature
                updated_feat = ogr.Feature(final_layer.GetLayerDefn())
                updated_feat.SetFID(fid)
                updated_feat.SetGeometry(amenity_feat.GetGeometryRef())

                # Copy all fields from amenity feature
                for i in range(amenity_feat.GetFieldCount()):
                    field_name = amenity_feat.GetFieldDefnRef(i).GetName()
                    if final_layer.GetLayerDefn().GetFieldIndex(field_name) >= 0:
                        value = amenity_feat.GetField(i)
                        updated_feat.SetField(field_name, value)

                final_layer.SetFeature(updated_feat)
                replaced_count += 1
                break

    print(f"    Replaced {replaced_count} features from amenities")

    # Clean up
    ds = None

    print(f"Merge completed for layer: {final_layer_name}")
