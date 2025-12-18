# src/rue_lib/public/operations.py
import geopandas as gpd
import pandas as pd
from osgeo import ogr

from rue_lib.core.definitions import ClusterTypes, ColorTypes


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
            if gdf.index.name == "fid":
                fid_min = int(gdf.index.min())
                gdf[id_field] = (gdf.index.astype("int64") - fid_min + start_at).astype("int64")
            else:
                gdf[id_field] = pd.Series(range(start_at, start_at + len(gdf)), dtype="int64")
        else:
            gdf[id_field] = gdf[id_field].astype(int)

    gdf.to_file(output_gpkg, layer=layer_name, driver="GPKG")
    return layer_name


def allocate_amenities(
    output_gpkg: str,
    parcel_layer_name: str,
    block_layer_name: str,
    open_spaces_layer_name: str,
    output_layer_name: str,
    amen_percent: float = 10.0,
) -> str:
    """
    Allocate amenities from remaining blocks after open space allocation.

    Implements Mobius logic:
    1. Get remaining blocks (not allocated as open space)
    2. Filter through parts (off_grid blocks and their attached on-grid parts)
    3. Sort by distance to site center
    4. Allocate until amenity percentage requirement is met

    Args:
        output_gpkg: Path to output GeoPackage
        parcel_layer_name: Name for parcel layer
        block_layer_name: Name for blocks layer
        open_spaces_layer_name: Name of the open spaces layer
        output_layer_name: Name for output amenities layer
        amen_percent: Percentage of site area to allocate as amenities

    Returns:
        Name of the output layer
    """
    print(f"Allocating amenities ({amen_percent}% of site area)...")

    # Load input layers
    gdf_blocks = gpd.read_file(output_gpkg, layer=block_layer_name)
    gdf_parcels = gpd.read_file(output_gpkg, layer=parcel_layer_name)
    gdf_open_spaces = gpd.read_file(output_gpkg, layer=open_spaces_layer_name)

    # Ensure we have an id column to track original feature IDs
    if "cluster_index" not in gdf_blocks.columns:
        gdf_blocks["cluster_index"] = range(len(gdf_blocks))

    # Get all open space block IDs to exclude them
    open_space_ids = set()
    if "cluster_index" in gdf_open_spaces.columns:
        open_space_ids = set(gdf_open_spaces["cluster_index"].dropna().astype(int))

    all_allocated_blocks = []
    site_id = 0

    # Process each parcel/site
    for _, parcel_row in gdf_parcels.iterrows():
        site_id += 1
        parcel_geom = parcel_row.geometry
        site_area = parcel_geom.area
        required_amen_area = site_area * (amen_percent / 100.0)

        print(f"  Site {site_id}: area={site_area:.2f}, required_amen={required_amen_area:.2f}")

        # Get all blocks that intersect this parcel (excluding open spaces)
        blocks_in_site = gdf_blocks[
            gdf_blocks.intersects(parcel_geom) & ~gdf_blocks["cluster_index"].isin(open_space_ids)
        ].copy()

        if blocks_in_site.empty:
            print(f"    No remaining blocks found for site {site_id}")
            continue

        print(f"    Remaining blocks found for site {site_id}: {len(blocks_in_site)}")

        parcel_centroid = parcel_geom.centroid
        blocks_in_site["centroid"] = blocks_in_site.geometry.centroid
        blocks_in_site["distance"] = blocks_in_site["centroid"].distance(parcel_centroid)
        blocks_in_site["block_area"] = blocks_in_site.geometry.area

        off_grid_types = [
            ClusterTypes.OFF_GRID_WARM,
            ClusterTypes.OFF_GRID_COLD,
            ClusterTypes.CONCAVE_CORNER,
        ]
        off_grid_blocks = blocks_in_site[blocks_in_site["type"].isin(off_grid_types)].copy()

        on_grid_block_types = ["loc", "loc_loc", "sec", "art", "off_grid2"]
        on_grid_blocks = blocks_in_site[blocks_in_site["type"].isin(on_grid_block_types)].copy()

        off_grid_groups = []
        for idx, og_block in off_grid_blocks.iterrows():
            group = {
                "root_idx": idx,
                "root": og_block,
                "attached_indices": [],
                "total_area": og_block["block_area"],
                "distance": og_block["distance"],
            }
            off_grid_groups.append(group)

        for idx, on_block in on_grid_blocks.iterrows():
            min_dist = float("inf")
            closest_group = None

            for group in off_grid_groups:
                root_block = group["root"]
                if root_block["block_id"] != on_block["block_id"]:
                    continue
                if not on_block.geometry.intersects(root_block.geometry):
                    continue

                dist = on_block["centroid"].distance(root_block["centroid"])
                if dist < min_dist:
                    min_dist = dist
                    closest_group = group

            if closest_group is not None:
                closest_group["attached_indices"].append(idx)
                closest_group["total_area"] += on_block["block_area"]

        off_grid_groups.sort(key=lambda x: x["distance"])

        print(f"    Found {len(blocks_in_site)} remaining blocks, {len(off_grid_groups)} groups")

        allocated_indices = []
        current_area = 0.0

        for group in off_grid_groups:
            if current_area >= required_amen_area:
                break

            root_idx = group["root_idx"]
            root_block = group["root"]
            allocated_indices.append(root_idx)
            current_area += root_block["block_area"]
            for attached_idx in group["attached_indices"]:
                if current_area >= required_amen_area:
                    break

                attached_block = blocks_in_site.loc[attached_idx]
                allocated_indices.append(attached_idx)
                current_area += attached_block["block_area"]

        print(
            f"    Allocated {len(allocated_indices)} blocks as amenities (area={current_area:.2f})"
        )

        for idx in allocated_indices:
            block_data = blocks_in_site.loc[idx].copy()
            block_data["allocation"] = "amenity"
            block_data["type"] = block_data["type"] + "_am"
            if "cluster_type" in block_data:
                block_data["cluster_type"] = block_data["cluster_type"] + "_am"
            block_data["color"] = ColorTypes[block_data["type"]]
            all_allocated_blocks.append(block_data)

    if all_allocated_blocks:
        gdf_output = gpd.GeoDataFrame(all_allocated_blocks, crs=gdf_blocks.crs)
        cols_to_drop = ["centroid", "distance", "block_area"]
        gdf_output = gdf_output.drop(columns=[c for c in cols_to_drop if c in gdf_output.columns])
        gdf_output.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")
        print(f"Amenities saved to layer: {output_layer_name}")
    else:
        gdf_output = gdf_blocks.iloc[:0].copy()
        gdf_output["allocation"] = pd.Series(dtype="str")
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
