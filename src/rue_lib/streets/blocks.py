# src/rue_lib/streets/blocks.py
import math

from osgeo import ogr

from .geometry_utils import break_linestring_by_angle


def merge_lines(
    gpkg_path,
    perpendicular_lines_layer,
    arterial_edge_lines_layer,
    street_blocks_layer,
    secondary_road_setback_layer,
    arterial_road_setback_layer,
    output_layer_name,
):
    """
    Merge multiple line sources into one layer, tagging setback provenance.

    Writes a unified line layer from:
      * `perpendicular_lines_layer`
      * `arterial_edge_lines_layer`
      * polygon boundaries of `street_blocks_layer`
      * polygon boundaries of `secondary_road_setback_layer`
      * polygon boundaries of `arterial_road_setback_layer`
      * (optional) site boundary

    Each output feature has:
      - source (str): producer of the line (e.g., 'perpendicular', 'arterial_edge',
        'street_blocks', 'secondary_setback', 'arterial_setback', 'site_boundary')
      - setback (str|NULL): 'secondary' or 'arterial' for setback boundaries; NULL otherwise
      - poly_id (int|NULL): FID of the source polygon when applicable
      - poly_area (float|NULL): Area of the source polygon when applicable

    Args:
        gpkg_path (str): GeoPackage path to read/write.
        perpendicular_lines_layer (str): Name of perpendicular lines layer.
        arterial_edge_lines_layer (str): Name of arterial edge lines layer.
        street_blocks_layer (str): Name of street blocks polygon layer.
        secondary_road_setback_layer (str): Name of secondary setback polygon layer.
        arterial_road_setback_layer (str): Name of arterial setback polygon layer.
        output_layer_name (str): Name of merged output line layer.

    Returns:
        None
    """
    ds = ogr.Open(gpkg_path, 1)
    line_id = 1

    perp_layer = ds.GetLayerByName(perpendicular_lines_layer)
    blocks_layer = ds.GetLayerByName(street_blocks_layer)
    secondary_layer = ds.GetLayerByName(secondary_road_setback_layer)
    arterial_setback_layer = ds.GetLayerByName(arterial_road_setback_layer)

    srs = perp_layer.GetSpatialRef()

    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbLineString)

    source_field = ogr.FieldDefn("source", ogr.OFTString)
    output_layer.CreateField(source_field)

    setback_field = ogr.FieldDefn("setback", ogr.OFTString)
    output_layer.CreateField(setback_field)

    poly_id_field = ogr.FieldDefn("line_id", ogr.OFTInteger)
    output_layer.CreateField(poly_id_field)

    poly_area_field = ogr.FieldDefn("poly_area", ogr.OFTReal)
    output_layer.CreateField(poly_area_field)

    line_count = 0

    for feature in perp_layer:
        geom = feature.GetGeometryRef()
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)
        out_feature.SetField("source", "perpendicular")
        out_feature.SetField("setback", None)
        out_feature.SetField("line_id", line_id)
        output_layer.CreateFeature(out_feature)
        out_feature = None
        line_count += 1
        line_id += 1

    # Add arterial_road_setback boundaries
    arterial_lines = []
    for feature in arterial_setback_layer:
        poly = feature.GetGeometryRef()
        boundary = poly.Boundary()
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(boundary)
        out_feature.SetField("source", "arterial_setback")
        out_feature.SetField("setback", "arterial")
        out_feature.SetField("line_id", line_id)
        output_layer.CreateFeature(out_feature)
        arterial_lines.append(boundary)
        out_feature = None
        line_count += 1
        line_id += 1

    for feature in blocks_layer:
        poly = feature.GetGeometryRef()
        boundary = poly.Boundary()
        for i in range(boundary.GetGeometryCount()):
            line = boundary.GetGeometryRef(i)
            segments = break_linestring_by_angle(line)
            for segment in segments:
                print(segment)
                out_feature = ogr.Feature(output_layer.GetLayerDefn())
                out_feature.SetGeometry(segment)
                out_feature.SetField("source", "street_blocks")
                out_feature.SetField("setback", None)
                out_feature.SetField("line_id", line_id)
                output_layer.CreateFeature(out_feature)
                out_feature = None
                line_count += 1
                line_id += 1

    for feature in secondary_layer:
        poly = feature.GetGeometryRef()
        boundary = poly.Boundary()

        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(boundary)
        out_feature.SetField("source", "secondary_setback")
        out_feature.SetField("setback", "secondary")
        out_feature.SetField("line_id", feature.GetFID())
        output_layer.CreateFeature(out_feature)
        out_feature = None
        line_count += 1
        line_id += 1

    ds = None
    print(f"Merged {line_count} lines into '{output_layer_name}'")


def polygonize_and_classify_blocks(
    gpkg_path,
    lines_layer_name,
    arterial_setback_layer_name,
    secondary_setback_layer_name,
    output_layer_name,
):
    """
    Create polygon blocks from merged lines and classify them by setback type.

    A block is classified based on spatial containment:
    - 'arterial_setback' if block centroid is inside arterial setback zone
    - 'secondary_setback' if block centroid is inside secondary setback zone (but not arterial)
    - 'off_grid' if block centroid is outside both setback zones

    Args:
        gpkg_path (str): Path to the GeoPackage
        lines_layer_name (str): Name of the merged lines layer
        arterial_setback_layer_name (str): Name of arterial setback polygon layer
        secondary_setback_layer_name (str): Name of secondary setback polygon layer
        output_layer_name (str): Name of the output classified blocks layer

    Returns:
        None
    """
    ds = ogr.Open(gpkg_path, 1)

    lines_layer = ds.GetLayerByName(lines_layer_name)
    arterial_layer = ds.GetLayerByName(arterial_setback_layer_name)
    secondary_layer = ds.GetLayerByName(secondary_setback_layer_name)
    srs = lines_layer.GetSpatialRef()

    # Create union of arterial setback polygons
    arterial_geoms = []
    for feature in arterial_layer:
        arterial_geoms.append(feature.GetGeometryRef().Clone())

    arterial_union = arterial_geoms[0]
    for geom in arterial_geoms[1:]:
        arterial_union = arterial_union.Union(geom)

    # Create union of secondary setback polygons
    secondary_geoms = []
    for feature in secondary_layer:
        secondary_geoms.append(feature.GetGeometryRef().Clone())

    secondary_union = secondary_geoms[0]
    for geom in secondary_geoms[1:]:
        secondary_union = secondary_union.Union(geom)

    # Collect all lines
    all_line_geoms = []
    for feature in lines_layer:
        geom = feature.GetGeometryRef().Clone()
        all_line_geoms.append(geom)

    # Create a union of all lines
    union_geom = all_line_geoms[0]
    for geom in all_line_geoms[1:]:
        union_geom = union_geom.Union(geom)

    # Polygonize
    polygons = []
    poly_geom = union_geom.Polygonize()

    if poly_geom is not None:
        geom_type = ogr.GT_Flatten(poly_geom.GetGeometryType())

        if geom_type == ogr.wkbPolygon:
            polygons.append(poly_geom.Clone())
        elif geom_type == ogr.wkbMultiPolygon or geom_type == ogr.wkbGeometryCollection:
            for i in range(poly_geom.GetGeometryCount()):
                sub_geom = poly_geom.GetGeometryRef(i)
                if ogr.GT_Flatten(sub_geom.GetGeometryType()) == ogr.wkbPolygon:
                    polygons.append(sub_geom.Clone())

    print(f"Created {len(polygons)} polygons from lines")

    # Remove existing output layer if it exists
    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    # Create output layer
    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)

    # Add classification fields
    block_id_field = ogr.FieldDefn("block_id", ogr.OFTInteger)
    output_layer.CreateField(block_id_field)

    classification_field = ogr.FieldDefn("block_type", ogr.OFTString)
    output_layer.CreateField(classification_field)

    area_field = ogr.FieldDefn("area", ogr.OFTReal)
    output_layer.CreateField(area_field)

    # Classify each polygon
    arterial_blocks = 0
    secondary_blocks = 0
    off_grid_blocks = 0

    for block_id, poly in enumerate(polygons):
        # Get centroid for classification
        centroid = poly.Centroid()

        # Classify based on spatial containment
        if arterial_union.Contains(centroid):
            block_type = "arterial_setback"
            arterial_blocks += 1
        elif secondary_union.Contains(centroid):
            block_type = "secondary_setback"
            secondary_blocks += 1
        else:
            block_type = "off_grid"
            off_grid_blocks += 1

        # Create output feature
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(poly)
        out_feature.SetField("block_id", block_id)
        out_feature.SetField("block_type", block_type)
        out_feature.SetField("area", poly.GetArea())

        output_layer.CreateFeature(out_feature)
        out_feature = None

    ds = None

    print("Block classification complete:")
    print(f"  - Arterial setback blocks: {arterial_blocks}")
    print(f"  - Secondary setback blocks: {secondary_blocks}")
    print(f"  - Off-grid blocks: {off_grid_blocks}")


def filter_offgrid_blocks(
    gpkg_path,
    classified_blocks_layer_name,
    output_layer_name,
    offgrid_preferred_width,
    offgrid_preferred_depth,
    arterial_setback_depth,
    secondary_setback_depth,
    arterial_preferred_width=None,  # Optional, defaults to offgrid_preferred_width
    secondary_preferred_width=None,  # Optional, defaults to offgrid_preferred_width
    area_threshold=0.5,  # 50% of preferred area
    squareness_threshold=0.6,  # How square-like (0-1, where 1 is perfect square)
):
    """
    Filter all blocks to remove those that don't meet shape and size requirements.

    Uses different size thresholds for different block types:
    - Arterial setback: arterial_preferred_width × arterial_setback_depth
    - Secondary setback: secondary_preferred_width × secondary_setback_depth
    - Off-grid: offgrid_preferred_width × offgrid_preferred_depth

    Args:
        gpkg_path (str): Path to the GeoPackage
        classified_blocks_layer_name (str): Name of the classified blocks layer
        output_layer_name (str): Name of the output filtered blocks layer
        offgrid_preferred_width (float): Preferred block width for off-grid blocks in meters
        offgrid_preferred_depth (float): Preferred block depth for off-grid blocks in meters
        arterial_setback_depth (float): Setback depth for arterial roads in meters
        secondary_setback_depth (float): Setback depth for secondary roads in meters
        arterial_preferred_width (float): Preferred width for arterial blocks
        (defaults to offgrid_preferred_width)
        secondary_preferred_width (float): Preferred width for secondary blocks
        (defaults to offgrid_preferred_width)
        area_threshold (float): Minimum area as fraction of preferred area (default 0.5 = 50%)
        squareness_threshold (float): Minimum squareness ratio (0-1, default 0.6)

    Returns:
        None
    """
    ds = ogr.Open(gpkg_path, 1)

    classified_layer = ds.GetLayerByName(classified_blocks_layer_name)
    srs = classified_layer.GetSpatialRef()

    # Set default widths if not provided
    if arterial_preferred_width is None:
        arterial_preferred_width = offgrid_preferred_width
    if secondary_preferred_width is None:
        secondary_preferred_width = offgrid_preferred_width

    # Calculate minimum required areas for each block type
    arterial_preferred_area = arterial_preferred_width * arterial_setback_depth
    arterial_min_area = arterial_preferred_area * area_threshold

    secondary_preferred_area = secondary_preferred_width * secondary_setback_depth
    secondary_min_area = secondary_preferred_area * area_threshold

    offgrid_preferred_area = offgrid_preferred_width * offgrid_preferred_depth
    offgrid_min_area = offgrid_preferred_area * area_threshold

    # Remove existing output layer if it exists
    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    # Create output layer with same schema
    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)

    # Copy field definitions from input layer
    classified_layer_defn = classified_layer.GetLayerDefn()
    for i in range(classified_layer_defn.GetFieldCount()):
        field_defn = classified_layer_defn.GetFieldDefn(i)
        output_layer.CreateField(field_defn)

    # Add filtering metadata fields
    squareness_field = ogr.FieldDefn("squareness", ogr.OFTReal)
    output_layer.CreateField(squareness_field)

    compactness_field = ogr.FieldDefn("compactness", ogr.OFTReal)
    output_layer.CreateField(compactness_field)

    vertices_field = ogr.FieldDefn("vertices", ogr.OFTInteger)
    output_layer.CreateField(vertices_field)

    filtered_field = ogr.FieldDefn("kept_reason", ogr.OFTString)
    output_layer.CreateField(filtered_field)

    def count_vertices(poly):
        """Count vertices in the exterior ring of a polygon."""
        exterior_ring = poly.GetGeometryRef(0)
        # Subtract 1 because the ring is closed (first point = last point)
        return exterior_ring.GetPointCount() - 1

    def calculate_compactness(poly):
        """
        Calculate compactness using isoperimetric quotient.

        Compactness = 4π × area / perimeter²

        For a circle: 1.0 (most compact)
        For a square: ~0.785
        For a triangle: ~0.60
        For elongated shapes: < 0.5
        """
        area = poly.GetArea()
        perimeter = poly.Boundary().Length()

        if perimeter == 0:
            return 0.0

        compactness = (4 * math.pi * area) / (perimeter**2)
        return compactness

    def calculate_squareness(poly):
        """
        Calculate how square-like a polygon is using multiple metrics.

        Returns value between 0 and 1, where 1 is a perfect square.
        """
        # 1. Check vertex count (squares should have ~4 vertices)
        vertices = count_vertices(poly)

        # Penalize heavily if less than 4 vertices (triangles, etc.)
        if vertices < 4:
            vertex_score = 0.0
        elif vertices == 4:
            vertex_score = 1.0
        else:
            # For more than 4 vertices, decay slowly (irregular squares are fine)
            vertex_score = max(0.3, 1.0 - (vertices - 4) * 0.05)

        # 2. Calculate compactness
        compactness = calculate_compactness(poly)
        # Normalize: squares have compactness ~0.785
        # Scale so that 0.785 -> 1.0, and values below that scale proportionally
        compactness_score = min(1.0, compactness / 0.785)

        # 3. Check aspect ratio of bounding box
        envelope = poly.GetEnvelope()
        bbox_width = envelope[1] - envelope[0]
        bbox_height = envelope[3] - envelope[2]

        if bbox_width == 0 or bbox_height == 0:
            aspect_score = 0.0
        else:
            aspect_ratio = min(bbox_width, bbox_height) / max(bbox_width, bbox_height)
            aspect_score = aspect_ratio

        # 4. Check how much of bounding box is filled
        bbox_area = bbox_width * bbox_height
        poly_area = poly.GetArea()

        if bbox_area == 0:
            extent_score = 0.0
        else:
            extent_ratio = poly_area / bbox_area
            extent_score = extent_ratio

        # Combined squareness metric with weights:
        # - Vertex count: 30% (important to filter triangles)
        # - Compactness: 30% (important for overall shape)
        # - Aspect ratio: 25% (should be roughly square)
        # - Extent ratio: 15% (should fill bbox reasonably)
        squareness = (
            vertex_score * 0.30
            + compactness_score * 0.30
            + aspect_score * 0.25
            + extent_score * 0.15
        )

        return squareness

    # Statistics by type
    stats = {
        "arterial_setback": {"kept": 0, "removed_size": 0, "removed_shape": 0},
        "secondary_setback": {"kept": 0, "removed_size": 0, "removed_shape": 0},
        "off_grid": {"kept": 0, "removed_size": 0, "removed_shape": 0},
    }

    for feature in classified_layer:
        poly = feature.GetGeometryRef()
        block_type = feature.GetField("block_type")
        area = poly.GetArea()

        # Get the appropriate minimum area for this block type
        if block_type == "arterial_setback":
            min_area = 0
        elif block_type == "secondary_setback":
            min_area = 0
        else:  # off_grid
            min_area = offgrid_min_area

        # Calculate shape metrics
        squareness = calculate_squareness(poly)
        compactness = calculate_compactness(poly)
        vertices = count_vertices(poly)

        keep_block = True
        reason = ""

        # Check area requirement
        if area < min_area:
            keep_block = False
            reason = f"area_too_small ({area:.0f} < {min_area:.0f})"
            stats[block_type]["removed_size"] += 1

        # Check squareness requirement (only if area is sufficient)
        elif squareness < squareness_threshold:
            keep_block = False
            reason = f"not_square (sq={squareness:.2f}, cmp={compactness:.2f}, v={vertices})"
            stats[block_type]["removed_shape"] += 1
        else:
            reason = "meets_requirements"
            stats[block_type]["kept"] += 1

        if keep_block:
            out_feature = ogr.Feature(output_layer.GetLayerDefn())
            out_feature.SetGeometry(poly)

            # Copy all fields
            for i in range(classified_layer_defn.GetFieldCount()):
                field_name = classified_layer_defn.GetFieldDefn(i).GetNameRef()
                out_feature.SetField(field_name, feature.GetField(i))

            out_feature.SetField("squareness", squareness)
            out_feature.SetField("compactness", compactness)
            out_feature.SetField("vertices", vertices)
            out_feature.SetField("kept_reason", reason)
            output_layer.CreateFeature(out_feature)
            out_feature = None

    ds = None

    # Print statistics
    print("\nBlock filtering complete:")
    print(f"\nArterial Setback Blocks (min area: {arterial_min_area:.0f} m²):")
    print(f"  - Kept: {stats['arterial_setback']['kept']}")
    print(f"  - Removed (too small): {stats['arterial_setback']['removed_size']}")
    print(f"  - Removed (not square): {stats['arterial_setback']['removed_shape']}")

    print(f"\nSecondary Setback Blocks (min area: {secondary_min_area:.0f} m²):")
    print(f"  - Kept: {stats['secondary_setback']['kept']}")
    print(f"  - Removed (too small): {stats['secondary_setback']['removed_size']}")
    print(f"  - Removed (not square): {stats['secondary_setback']['removed_shape']}")

    print(f"\nOff-Grid Blocks (min area: {offgrid_min_area:.0f} m²):")
    print(f"  - Kept: {stats['off_grid']['kept']}")
    print(f"  - Removed (too small): {stats['off_grid']['removed_size']}")
    print(f"  - Removed (not square): {stats['off_grid']['removed_shape']}")

    total_kept = sum(s["kept"] for s in stats.values())
    total_removed = sum(s["removed_size"] + s["removed_shape"] for s in stats.values())

    print("\nTotal Summary:")
    print(f"  - Total kept: {total_kept}")
    print(f"  - Total removed: {total_removed}")
    print(f"  - Squareness threshold: {squareness_threshold:.2f}")


def filter_classified_blocks(
    gpkg_path: str,
    input_layer_name: str,
    output_layer_name: str,
    allowed_types=("arterial_setback", "secondary_setback"),
):
    """
    Copy only features whose block_type is in allowed_types into output_layer_name.
    """
    ds = ogr.Open(gpkg_path, 1)
    in_lyr = ds.GetLayerByName(input_layer_name)
    if in_lyr is None:
        raise RuntimeError(f"Layer not found: {input_layer_name}")

    srs = in_lyr.GetSpatialRef()
    # drop output if exists
    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    out_lyr = ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)

    # Copy all fields from input
    in_defn = in_lyr.GetLayerDefn()
    for i in range(in_defn.GetFieldCount()):
        out_lyr.CreateField(in_defn.GetFieldDefn(i).Clone())

    kept = 0
    total = 0
    in_lyr.ResetReading()
    for feat in in_lyr:
        total += 1
        if feat.GetField("block_type") in allowed_types:
            out_f = ogr.Feature(out_lyr.GetLayerDefn())
            out_f.SetGeometry(feat.GetGeometryRef().Clone())
            # copy fields
            for i in range(in_defn.GetFieldCount()):
                name = in_defn.GetFieldDefn(i).GetNameRef()
                out_f.SetField(name, feat.GetField(name))
            out_lyr.CreateFeature(out_f)
            out_f = None
            kept += 1

    ds = None
    print(f"Filtered {kept}/{total} classified blocks -> '{output_layer_name}'")
