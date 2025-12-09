# src/rue_lib/public/operations.py
from osgeo import ogr, osr


# TODO : We need to do this on the cluster module instead
def classify_cluster_blocks(
    input_gpkg: str,
    blocks_layer_name: str,
    parcels_path: str,
    output_gpkg: str,
    output_layer_name: str = "00_classified_blocks",
) -> str:
    """
    Classify cluster blocks as on-grid or off-grid based on neighbor presence.

    A block is considered "on_grid" if it sits on the edge of a cluster
    (missing a neighbor on any side). Blocks that have neighbors on all sides
    are tagged "off_grid". The result is written to a new layer with the
    added attributes ``parent_part`` and ``type`` containing the classification.
    """
    print("Classifying cluster blocks (on_grid vs off_grid)...")

    blocks_ds = ogr.Open(input_gpkg, 0)
    if blocks_ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    blocks_layer = blocks_ds.GetLayerByName(blocks_layer_name)
    if blocks_layer is None:
        raise ValueError(f"Layer {blocks_layer_name} not found in {input_gpkg}")

    parcels_ds = ogr.Open(parcels_path, 0)
    if parcels_ds is None:
        raise ValueError(f"Could not open {parcels_path}")

    parcels_layer = parcels_ds.GetLayer(0)

    blocks_srs = blocks_layer.GetSpatialRef()
    parcels_srs = parcels_layer.GetSpatialRef()

    transform = None
    if parcels_srs and blocks_srs and not parcels_srs.IsSame(blocks_srs):
        transform = osr.CoordinateTransformation(parcels_srs, blocks_srs)
        print(f"  Reprojecting parcels from {parcels_srs.GetName()} to {blocks_srs.GetName()}")

    driver = ogr.GetDriverByName("GPKG")
    output_ds = driver.Open(output_gpkg, 1)
    if output_ds is None:
        output_ds = driver.CreateDataSource(output_gpkg)
    if output_ds is None:
        raise ValueError(f"Could not create or open {output_gpkg} for writing")

    if output_ds.GetLayerByName(output_layer_name):
        output_ds.DeleteLayer(output_layer_name)

    geom_type = blocks_layer.GetGeomType()
    output_layer = output_ds.CreateLayer(output_layer_name, blocks_srs, geom_type)

    layer_defn = blocks_layer.GetLayerDefn()
    field_names = []
    for i in range(layer_defn.GetFieldCount()):
        field_defn = layer_defn.GetFieldDefn(i)
        field_names.append(field_defn.GetNameRef())
        output_layer.CreateField(field_defn)

    if output_layer.FindFieldIndex("parent_part", 1) == -1:
        parent_field = ogr.FieldDefn("parent_part", ogr.OFTString)
        parent_field.SetWidth(32)
        output_layer.CreateField(parent_field)

    output_layer.CreateField(ogr.FieldDefn("type", ogr.OFTString))

    def _adaptive_tolerance(envelope: tuple[float, float, float, float]) -> float:
        width = envelope[1] - envelope[0]
        height = envelope[3] - envelope[2]
        scale = max(width, height, 1.0)
        return max(scale * 1e-6, 1e-9)

    def _extract_boundary_edge_centers(geom: ogr.Geometry) -> list[tuple[float, float]]:
        """Extract center points of all boundary edges from polygonal geometry."""
        edge_centers: list[tuple[float, float]] = []

        def _process_ring(ring: ogr.Geometry):
            point_count = ring.GetPointCount()
            if point_count < 2:
                return
            for i in range(point_count - 1):
                x1, y1, _z1 = ring.GetPoint(i)
                x2, y2, _z2 = ring.GetPoint(i + 1)
                center_x = (x1 + x2) / 2.0
                center_y = (y1 + y2) / 2.0
                edge_centers.append((center_x, center_y))

        def _process_polygon(poly: ogr.Geometry):
            if poly.GetGeometryCount() > 0:
                ring = poly.GetGeometryRef(0)
                if ring is not None:
                    _process_ring(ring)

        gtype = geom.GetGeometryType()
        if gtype in (ogr.wkbPolygon, ogr.wkbPolygon25D):
            _process_polygon(geom)
        elif gtype in (ogr.wkbMultiPolygon, ogr.wkbMultiPolygon25D):
            for i in range(geom.GetGeometryCount()):
                sub = geom.GetGeometryRef(i)
                if sub is not None:
                    _process_polygon(sub)

        return edge_centers

    def _has_neighbors_all_sides(
        target_geom: ogr.Geometry,
        other_geoms: list[ogr.Geometry],
    ) -> bool:
        """
        Check if each boundary edge center point touches at least one neighboring geometry.

        Uses center points of boundary edges instead of vertices for more robust detection.
        """
        env = target_geom.GetEnvelope()
        tol = _adaptive_tolerance(env)
        probe_buffer = max(tol * 10, 1e-5)  # Slightly larger buffer for edge centers

        edge_centers = _extract_boundary_edge_centers(target_geom)
        if not edge_centers:
            return False

        for x, y in edge_centers:
            pt = ogr.Geometry(ogr.wkbPoint)
            pt.AddPoint(x, y)
            probe = pt.Buffer(probe_buffer)

            touches_neighbor = any(
                neighbor is not None and neighbor.Intersects(probe) for neighbor in other_geoms
            )
            if not touches_neighbor:
                return False

        return True

    site_id = 0
    total_blocks = 0

    for parcel_feat in parcels_layer:
        site_id += 1
        parcel_geom = parcel_feat.GetGeometryRef()
        if parcel_geom is None:
            continue

        parcel_geom = parcel_geom.Clone()
        if transform:
            parcel_geom.Transform(transform)

        blocks_layer.SetSpatialFilter(parcel_geom)

        block_entries = []
        for block_feat in blocks_layer:
            block_geom = block_feat.GetGeometryRef()
            if block_geom is None or block_geom.IsEmpty():
                continue

            geom_clone = block_geom.Clone()
            attrs = {name: block_feat.GetField(name) for name in field_names}

            block_entries.append(
                {
                    "geom": geom_clone,
                    "attrs": attrs,
                }
            )

        blocks_layer.SetSpatialFilter(None)
        blocks_layer.ResetReading()

        if not block_entries:
            print(f"  Site {site_id}: no blocks found")
            continue

        geoms = [entry["geom"] for entry in block_entries]
        print(f"  Site {site_id}: classifying {len(block_entries)} blocks")

        for idx, entry in enumerate(block_entries):
            other_geoms = geoms[:idx] + geoms[idx + 1 :]
            neighbors_all_sides = _has_neighbors_all_sides(entry["geom"], other_geoms)
            block_type = "off_grid" if neighbors_all_sides else "on_grid"

            out_feat = ogr.Feature(output_layer.GetLayerDefn())
            out_feat.SetGeometry(entry["geom"])

            for name, value in entry["attrs"].items():
                out_feat.SetField(name, value)

            out_feat.SetField("parent_part", block_type)
            out_feat.SetField("type", block_type)
            output_layer.CreateFeature(out_feat)
            out_feat = None
            total_blocks += 1

    blocks_ds = None
    parcels_ds = None
    output_ds = None

    print(f"Classification complete: wrote {total_blocks} blocks to layer {output_layer_name}")
    return output_layer_name


def allocate_open_spaces(
    input_gpkg: str,
    blocks_layer_name: str,
    parcels_path: str,
    output_gpkg: str,
    output_layer_name: str = "02_open_spaces",
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
        input_gpkg: Path to input GeoPackage containing all blocks
        blocks_layer_name: Name of the all blocks layer
        parcels_path: Path to parcels/sites GeoJSON or GeoPackage
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output open spaces layer
        open_percent: Percentage of site area to allocate as open space

    Returns:
        Name of the output layer
    """
    print(f"Allocating open spaces ({open_percent}% of site area)...")

    # Open input layers
    blocks_ds = ogr.Open(input_gpkg, 0)
    if blocks_ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    blocks_layer = blocks_ds.GetLayerByName(blocks_layer_name)
    if blocks_layer is None:
        raise ValueError(f"Layer {blocks_layer_name} not found")

    parcels_ds = ogr.Open(parcels_path, 0)
    if parcels_ds is None:
        raise ValueError(f"Could not open {parcels_path}")

    parcels_layer = parcels_ds.GetLayer(0)

    # Open output dataset
    output_ds = ogr.Open(output_gpkg, 1)
    if output_ds is None:
        raise ValueError(f"Could not open {output_gpkg} for writing")

    # Get spatial reference systems
    blocks_srs = blocks_layer.GetSpatialRef()
    parcels_srs = parcels_layer.GetSpatialRef()

    # Create coordinate transformation
    transform = None
    if parcels_srs and blocks_srs:
        if not parcels_srs.IsSame(blocks_srs):
            transform = osr.CoordinateTransformation(parcels_srs, blocks_srs)

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

    site_id = 0

    # Process each parcel/site
    for parcel_feat in parcels_layer:
        site_id += 1
        parcel_geom = parcel_feat.GetGeometryRef().Clone()

        if transform:
            parcel_geom.Transform(transform)

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
        off_grid_blocks = [b for b in all_blocks if b["type"] == "off_grid"]
        if not off_grid_blocks:
            print(f"    No off-grid blocks found for site {site_id}")
            continue

        off_grid_blocks.sort(key=lambda x: x["distance"])
        central_block = off_grid_blocks[0]
        central_centroid = central_block["centroid"]

        print(f"    Found {len(all_blocks)} blocks, central block identified")

        # Group blocks: off-grid roots with attached on-grid parts
        # Similar to reorderParts1 in Mobius
        off_grid_groups = []
        for og_block in off_grid_blocks:
            if og_block["block_id"] != central_block["block_id"]:
                continue
            group = {
                "root": og_block,
                "attached": [],
                "total_area": og_block["area"],
                "distance_to_center": central_centroid.Distance(og_block["centroid"]),
            }
            off_grid_groups.append(group)

        # Attach on-grid blocks to nearest off-grid root
        on_grid_blocks = [b for b in all_blocks if b["type"] == "on_grid"]
        on_grid_blocks.sort(key=lambda x: x["distance"])
        for on_block in on_grid_blocks:
            min_dist = float("inf")
            closest_group = None
            for group in off_grid_groups:
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
            fid = feat.GetFID()
            return (b.get("block_id"), fid)

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

        if current_area < required_open_area:
            for _i, group in enumerate(off_grid_groups):
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

            out_feat.SetField("allocation", "open_space")
            out_feat.SetField("group_id", group_id)
            out_feat.SetField("group_rank", group_rank)
            out_feat.SetField("is_root", is_root)
            out_feat.SetField("dist_to_center", dist_to_center)

            output_layer.CreateFeature(out_feat)
            out_feat = None

    # Clean up
    blocks_ds = None
    parcels_ds = None
    output_ds = None

    print(f"Open spaces saved to layer: {output_layer_name}")
    return output_layer_name


def allocate_amenities(
    input_gpkg: str,
    blocks_layer_name: str,
    open_spaces_layer_name: str,
    parcels_path: str,
    output_gpkg: str,
    output_layer_name: str = "03_amenities",
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
        input_gpkg: Path to input GeoPackage containing all blocks
        blocks_layer_name: Name of the all blocks layer
        open_spaces_layer_name: Name of the open spaces layer
        parcels_path: Path to parcels/sites GeoJSON or GeoPackage
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output amenities layer
        amen_percent: Percentage of site area to allocate as amenities

    Returns:
        Name of the output layer
    """
    print(f"Allocating amenities ({amen_percent}% of site area)...")

    # Open input layers
    blocks_ds = ogr.Open(input_gpkg, 0)
    if blocks_ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    blocks_layer = blocks_ds.GetLayerByName(blocks_layer_name)
    if blocks_layer is None:
        raise ValueError(f"Layer {blocks_layer_name} not found")

    open_spaces_layer = blocks_ds.GetLayerByName(open_spaces_layer_name)
    if open_spaces_layer is None:
        raise ValueError(f"Layer {open_spaces_layer_name} not found")

    parcels_ds = ogr.Open(parcels_path, 0)
    if parcels_ds is None:
        raise ValueError(f"Could not open {parcels_path}")

    parcels_layer = parcels_ds.GetLayer(0)

    # Open output dataset
    output_ds = ogr.Open(output_gpkg, 1)
    if output_ds is None:
        raise ValueError(f"Could not open {output_gpkg} for writing")

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
        block_id = feat.GetField("block_id")
        if block_id is not None:
            open_space_ids.add(block_id)
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
            block_id = block_feat.GetField("block_id")
            # Skip if already allocated as open space
            if block_id is not None and block_id in open_space_ids:
                continue

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

        # Filter through parts: get off-grid blocks with their attached on-grid parts
        off_grid_blocks = [b for b in remaining_blocks if b["type"] == "off_grid"]
        on_grid_blocks = [b for b in remaining_blocks if b["type"] == "on_grid"]

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
                if current_area >= required_amen_area:
                    break
                allocated_blocks.append(attached)
                current_area += attached["area"]

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

            out_feat.SetField("allocation", "amenity")

            output_layer.CreateFeature(out_feat)
            out_feat = None

    # Clean up
    blocks_ds = None
    parcels_ds = None
    output_ds = None

    print(f"Amenities saved to layer: {output_layer_name}")
    return output_layer_name
