# src/rue_lib/public/operations.py
from osgeo import ogr, osr

from rue_lib.core.definitions import ClusterTypes, ColorTypes


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
