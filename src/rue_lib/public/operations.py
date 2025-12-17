# src/rue_lib/public/operations.py
from osgeo import ogr, osr

from rue_lib.core.definitions import ClusterTypes, ColorTypes


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

    # Open input layers
    output_ds = ogr.Open(output_gpkg, 1)
    blocks_layer = output_ds.GetLayerByName(block_layer_name)
    parcels_layer = output_ds.GetLayerByName(parcel_layer_name)

    # Get spatial reference systems
    blocks_srs = blocks_layer.GetSpatialRef()

    # Create output layer
    if output_ds.GetLayerByName(output_layer_name):
        output_ds.DeleteLayer(output_layer_name)

    output_layer = output_ds.CreateLayer(output_layer_name, blocks_srs, ogr.wkbPolygon)

    # Create fields (copy from blocks layer plus type)
    layer_defn = blocks_layer.GetLayerDefn()
    for i in range(layer_defn.GetFieldCount()):
        field_defn = layer_defn.GetFieldDefn(i)
        output_layer.CreateField(field_defn)

    output_layer.CreateField(ogr.FieldDefn("allocation", ogr.OFTString))
    output_layer.CreateField(ogr.FieldDefn("group_id", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("is_root", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("group_rank", ogr.OFTInteger))

    dist_field = ogr.FieldDefn("dist_to_center", ogr.OFTReal)
    output_layer.CreateField(dist_field)

    cluster_index_field = ogr.FieldDefn("cluster_index", ogr.OFTInteger)
    output_layer.CreateField(cluster_index_field)

    site_id = 0

    # Process each parcel/site
    for parcel_feat in parcels_layer:
        site_id += 1
        parcel_geom = parcel_feat.GetGeometryRef().Clone()
        site_area = parcel_geom.Area()
        required_open_area = site_area * (open_percent / 100.0)

        print(f"  Site {site_id}: area={site_area:.2f}, required_open={required_open_area:.2f}")

        # Get all blocks for this site
        blocks_layer.SetSpatialFilter(parcel_geom)
        all_blocks = []
        parcel_centroid = parcel_geom.Centroid()

        for block_feat in blocks_layer:
            block_geom = block_feat.GetGeometryRef()
            if block_geom.Intersects(parcel_geom):
                block_type = block_feat.GetField("type")
                part_type = block_feat.GetField("part_type")
                is_corner = part_type == "corner" if part_type else False

                block_centroid = block_geom.Centroid()
                distance = parcel_centroid.Distance(block_centroid)
                all_blocks.append(
                    {
                        "block_id": block_feat.GetField("block_id"),
                        "feature": block_feat.Clone(),
                        "geometry": block_geom.Clone(),
                        "centroid": block_centroid,
                        "is_corner": is_corner,
                        "type": block_type,
                        "distance": distance,
                        "area": block_geom.Area(),
                    }
                )

        blocks_layer.ResetReading()

        if not all_blocks:
            print(f"    No blocks found for site {site_id}")
            continue

        # Find central block (closest off-grid block to site centroid)
        off_grid_blocks = [
            b
            for b in all_blocks
            if b["type"]
            in [ClusterTypes.OFF_GRID_WARM, ClusterTypes.OFF_GRID_COLD, ClusterTypes.CONCAVE_CORNER]
        ]

        if not off_grid_blocks:
            print(f"    No off-grid blocks found for site {site_id}")
            continue

        off_grid_blocks.sort(key=lambda x: x["distance"])
        central_block = off_grid_blocks[0]
        central_centroid = central_block["centroid"]
        print([group["type"] for group in off_grid_blocks])

        print(f"    Found {len(all_blocks)} blocks, central block identified")

        # Group blocks: off-grid roots with attached on-grid parts
        # Similar to reorderParts1 in Mobius
        off_grid_groups = []
        for og_block in off_grid_blocks:
            group = {
                "root": og_block,
                "attached": [],
                "total_area": og_block["area"],
                "distance_to_center": central_centroid.Distance(og_block["centroid"]),
            }
            off_grid_groups.append(group)

        # Attach on-grid blocks to nearest off-grid root
        on_grid_block_types = ["loc", "loc_loc", "sec", "art", "off_grid2"]
        on_grid_blocks = [b for b in all_blocks if b["type"] in on_grid_block_types]
        on_grid_blocks.sort(key=lambda x: x["distance"])
        for on_block in on_grid_blocks:
            min_dist = float("inf")
            closest_group = None
            for group in off_grid_groups:
                if group["root"]["block_id"] != on_block["block_id"]:
                    continue
                if not on_block["geometry"].Intersects(group["root"]["geometry"]):
                    continue
                dist = on_block["centroid"].Distance(group["root"]["centroid"])
                if dist < min_dist:
                    min_dist = dist
                    closest_group = group

            on_block["distance_to_center"] = central_centroid.Distance(on_block["centroid"])

            if closest_group is not None:
                closest_group["attached"].append(on_block)
                closest_group["total_area"] += on_block["area"]

        # Sort groups by distance to central block
        off_grid_groups.sort(key=lambda x: x["distance_to_center"])

        print(f"    Created {len(off_grid_groups)} off-grid groups")

        allocated_blocks: list[tuple[dict, int, int, int, float]] = []
        allocated_ids: set[tuple[int | None, int | None]] = set()
        current_area = 0.0

        def _block_key(b: dict) -> tuple[int | None, int | None]:
            """
            Use (block_id, FID) as uniqueness key so we don't double-allocate.
            """
            feat = b["feature"]
            return feat.GetFID()

        # Re-sort on-grid blocks:
        # 1. Keep first 2
        # 2. Reverse the remaining
        # 3. Merge
        if len(off_grid_groups) > 2:
            first_two = off_grid_groups[:2]
            remaining = off_grid_groups[2:]
            remaining.reverse()
            off_grid_groups = first_two + remaining

        print([group["root"]["type"] for group in off_grid_groups])
        for i, group in enumerate(off_grid_groups):
            if current_area >= required_open_area:
                break

            group_id = i + 1
            group_rank = i + 1
            dist_to_center = group["distance_to_center"]

            root_block = group["root"]
            key = _block_key(root_block)
            if key not in allocated_ids and current_area + root_block["area"] < required_open_area:
                allocated_blocks.append((root_block, group_id, 1, group_rank, dist_to_center))
                allocated_ids.add(key)
                current_area += root_block["area"]

                # This is for looping attached
                for attached in group["attached"]:
                    current_area += attached["area"]
                    if current_area >= required_open_area:
                        break
                    key = _block_key(attached)
                    if key in allocated_ids:
                        continue
                    allocated_blocks.append(
                        (attached, group_id, 0, group_rank, attached["distance_to_center"])
                    )
                    allocated_ids.add(key)

        # TODO: If still need more, add remaining off-grid blocks (sorted by distance)
        # Elsewhere in the site – open space
        # Consider only:
        #   concave_corner, off_grid0, off_grid1.
        # Prefer:
        # A combination of:
        # Concave corners that can “pull in” whole blocks,
        # Off-grid strips extruded from road edges,
        # Clusters of off_grid1.
        #   Spread them spatially using neighbor-avoidance.

        print(
            f"    Allocated {len(allocated_blocks)} blocks as open space (area={current_area:.2f})"
        )

        for block, group_id, is_root, group_rank, dist_to_center in allocated_blocks:
            out_feat = ogr.Feature(output_layer.GetLayerDefn())
            out_feat.SetGeometry(block["geometry"])

            src_feat = block["feature"]
            for i in range(src_feat.GetFieldCount()):
                value = src_feat.GetField(i)
                out_feat.SetField(i, value)

            new_type = out_feat.GetField("type") + "_os"
            new_cluster_type = out_feat.GetField("cluster_type") + "_os"
            out_feat.SetField("allocation", "open_space")
            out_feat.SetField("group_id", group_id)
            out_feat.SetField("group_rank", group_rank)
            out_feat.SetField("is_root", is_root)
            out_feat.SetField("dist_to_center", dist_to_center)
            out_feat.SetField("type", new_type)
            out_feat.SetField("cluster_type", new_cluster_type)
            out_feat.SetField("color", ColorTypes[new_type])
            out_feat.SetField("cluster_index", src_feat.GetFID())
            output_layer.CreateFeature(out_feat)
            out_feat = None

    # Clean up
    output_ds = None

    print(f"Open spaces saved to layer: {output_layer_name}")
    return output_layer_name


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

    # Open input layers
    output_ds = ogr.Open(output_gpkg, 1)
    blocks_layer = output_ds.GetLayerByName(block_layer_name)
    parcels_layer = output_ds.GetLayerByName(parcel_layer_name)

    open_spaces_layer = output_ds.GetLayerByName(open_spaces_layer_name)
    if open_spaces_layer is None:
        raise ValueError(f"Layer {open_spaces_layer_name} not found")

    # Get spatial reference systems
    blocks_srs = blocks_layer.GetSpatialRef()
    parcels_srs = parcels_layer.GetSpatialRef()

    # Create coordinate transformation
    transform = None
    if parcels_srs and blocks_srs:
        if not parcels_srs.IsSame(blocks_srs):
            transform = osr.CoordinateTransformation(parcels_srs, blocks_srs)

    # Get all open space block IDs to exclude them
    open_space_ids = set()
    for feat in open_spaces_layer:
        cluster_index = feat.GetField("cluster_index")
        if cluster_index is not None:
            open_space_ids.add(cluster_index)
    open_spaces_layer.ResetReading()

    # Create output layer
    if output_ds.GetLayerByName(output_layer_name):
        output_ds.DeleteLayer(output_layer_name)

    output_layer = output_ds.CreateLayer(output_layer_name, blocks_srs, ogr.wkbPolygon)

    # Create fields (copy from blocks layer plus type)
    layer_defn = blocks_layer.GetLayerDefn()
    for i in range(layer_defn.GetFieldCount()):
        field_defn = layer_defn.GetFieldDefn(i)
        output_layer.CreateField(field_defn)

    output_layer.CreateField(ogr.FieldDefn("allocation", ogr.OFTString))
    output_layer.CreateField(ogr.FieldDefn("cluster_index", ogr.OFTInteger))

    site_id = 0

    # Process each parcel/site
    for parcel_feat in parcels_layer:
        site_id += 1
        parcel_geom = parcel_feat.GetGeometryRef().Clone()

        if transform:
            parcel_geom.Transform(transform)

        site_area = parcel_geom.Area()
        required_amen_area = site_area * (amen_percent / 100.0)

        print(f"  Site {site_id}: area={site_area:.2f}, required_amen={required_amen_area:.2f}")

        # Get all blocks for this site (excluding open spaces)
        blocks_layer.SetSpatialFilter(parcel_geom)
        parcel_centroid = parcel_geom.Centroid()
        remaining_blocks = []

        for block_feat in blocks_layer:
            cluster_index = block_feat.GetFID()
            # Skip if already allocated as open space
            if cluster_index in open_space_ids:
                continue

            block_id = block_feat.GetField("block_id")
            block_geom = block_feat.GetGeometryRef()
            if block_geom.Intersects(parcel_geom):
                block_type = block_feat.GetField("type")
                block_centroid = block_geom.Centroid()
                distance = parcel_centroid.Distance(block_centroid)

                remaining_blocks.append(
                    {
                        "feature": block_feat.Clone(),
                        "geometry": block_geom.Clone(),
                        "centroid": block_centroid,
                        "type": block_type,
                        "distance": distance,
                        "area": block_geom.Area(),
                        "block_id": block_id,
                    }
                )

        blocks_layer.ResetReading()

        if not remaining_blocks:
            print(f"    No remaining blocks found for site {site_id}")
            continue
        print(f"    Remaining blocks found for site {site_id} : {len(remaining_blocks)}")

        # Filter through parts: get off-grid blocks with their attached on-grid parts
        off_grid_blocks = [
            b
            for b in remaining_blocks
            if b["type"]
            in [ClusterTypes.OFF_GRID_WARM, ClusterTypes.OFF_GRID_COLD, ClusterTypes.CONCAVE_CORNER]
        ]
        on_grid_block_types = ["loc", "loc_loc", "sec", "art", "off_grid2"]
        on_grid_blocks = [b for b in remaining_blocks if b["type"] in on_grid_block_types]

        # Group blocks similar to open space allocation
        off_grid_groups = []
        for og_block in off_grid_blocks:
            group = {
                "root": og_block,
                "attached": [],
                "total_area": og_block["area"],
                "distance": og_block["distance"],
            }
            off_grid_groups.append(group)

        # Attach on-grid blocks to nearest off-grid root
        for on_block in on_grid_blocks:
            min_dist = float("inf")
            closest_group = None
            for group in off_grid_groups:
                if group["root"]["block_id"] != on_block["block_id"]:
                    continue
                if not on_block["geometry"].Intersects(group["root"]["geometry"]):
                    continue
                dist = on_block["centroid"].Distance(group["root"]["centroid"])
                if dist < min_dist:
                    min_dist = dist
                    closest_group = group

            if closest_group:
                closest_group["attached"].append(on_block)
                closest_group["total_area"] += on_block["area"]

        # Sort groups by distance to site center
        off_grid_groups.sort(key=lambda x: x["distance"])

        print(f"    Found {len(remaining_blocks)} remaining blocks, {len(off_grid_groups)} groups")

        # Allocate groups until amenity requirement is met
        allocated_blocks = []
        current_area = 0.0

        for group in off_grid_groups:
            if current_area >= required_amen_area:
                break
            allocated_blocks.append(group["root"])
            current_area += group["root"]["area"]
            print(group["root"]["type"], len(group["attached"]))
            for attached in group["attached"]:
                current_area += attached["area"]
                if current_area >= required_amen_area:
                    break
                allocated_blocks.append(attached)

        print(
            f"    Allocated {len(allocated_blocks)} blocks as amenities (area={current_area:.2f})"
        )

        # Write allocated blocks to output
        for block in allocated_blocks:
            out_feat = ogr.Feature(output_layer.GetLayerDefn())
            out_feat.SetGeometry(block["geometry"])

            # Copy all fields from original
            src_feat = block["feature"]
            for i in range(src_feat.GetFieldCount()):
                value = src_feat.GetField(i)
                out_feat.SetField(i, value)

            new_type = out_feat.GetField("type") + "_am"
            new_cluster_type = out_feat.GetField("cluster_type") + "_am"
            out_feat.SetField("allocation", "amenity")
            out_feat.SetField("type", new_type)
            out_feat.SetField("cluster_type", new_cluster_type)
            out_feat.SetField("color", ColorTypes[new_type])
            out_feat.SetField("cluster_index", src_feat.GetFID())
            output_layer.CreateFeature(out_feat)
            out_feat = None

    # Clean up
    output_ds = None

    print(f"Amenities saved to layer: {output_layer_name}")
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

        # Find matching feature in final layer
        final_layer.ResetReading()
        for final_feat in final_layer:
            final_cluster_index = final_feat.GetFID()

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

        # Find matching feature in final layer
        final_layer.ResetReading()
        for final_feat in final_layer:
            final_cluster_index = final_feat.GetFID()

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
