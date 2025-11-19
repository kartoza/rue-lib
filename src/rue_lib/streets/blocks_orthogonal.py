# src/rue_lib/streets/blocks_orthogonal.py
from __future__ import annotations

import math
from collections.abc import Iterable, Sequence
from dataclasses import dataclass

from osgeo import ogr
from shapely import ops as shapely_ops
from shapely import wkb as shapely_wkb
from shapely.geometry import GeometryCollection, Polygon


@dataclass
class BlockRecord:
    """Container for orthogonal block metrics."""

    block_id: int
    block_type: str
    geom: Polygon
    area: float
    quality: float = 0.0
    neighbors: int = 0
    vertex_count: int = 0
    compactness: float = 0.0
    aspect_ratio: float = 0.0
    extent_ratio: float = 0.0
    size_ratio: float = 0.0
    keep: bool = True


def _ogr_delete_layer(dataset: ogr.DataSource, layer_name: str) -> None:
    """Remove layer if it already exists."""

    for idx in range(dataset.GetLayerCount()):
        if dataset.GetLayerByIndex(idx).GetName() == layer_name:
            dataset.DeleteLayer(idx)
            break


def _iter_polygons(geom) -> Iterable[Polygon]:
    """Yield shapely polygons from any geometry container."""

    if geom is None:
        return

    if isinstance(geom, Polygon):
        if not geom.is_empty:
            yield geom
        return

    if isinstance(geom, GeometryCollection):
        for part in geom.geoms:
            yield from _iter_polygons(part)
        return

    if hasattr(geom, "geoms"):
        for part in geom.geoms:
            yield from _iter_polygons(part)


def _evaluate_grid_at_angle(
    polygon: Polygon,
    angle: float,
    cell_width: float,
    cell_depth: float,
    centroid_x: float,
    centroid_y: float,
) -> tuple[int, float]:
    """
    Helper function to evaluate grid coverage at a specific angle.

    Only counts cells that are:
    - Completely inside the polygon (not touching/intersecting boundary)
    - Perfect rectangles (not clipped)

    Returns:
        Tuple of (cell_count, total_area)
    """
    from shapely import affinity
    from shapely.geometry import box

    # Rotate the polygon by -angle
    rotated_polygon = affinity.rotate(polygon, -angle, origin=(centroid_x, centroid_y))

    # Get bounding box of rotated polygon
    minx, miny, maxx, maxy = rotated_polygon.bounds

    # Calculate grid dimensions
    num_cols = max(1, int(math.ceil((maxx - minx) / cell_width)))
    num_rows = max(1, int(math.ceil((maxy - miny) / cell_depth)))

    # Target area for a perfect cell
    perfect_cell_area = cell_width * cell_depth
    area_tolerance = 0.01  # 1% tolerance for floating point errors

    # Count cells that fit inside
    cell_count = 0
    total_area = 0.0

    for row in range(num_rows):
        for col in range(num_cols):
            cell_minx = minx + col * cell_width
            cell_miny = miny + row * cell_depth
            cell_maxx = cell_minx + cell_width
            cell_maxy = cell_miny + cell_depth

            cell = box(cell_minx, cell_miny, cell_maxx, cell_maxy)

            # Check if cell is completely within the polygon
            if not rotated_polygon.contains(cell):
                continue

            # Verify the cell is a perfect rectangle (not clipped)
            # by checking if its area matches the expected area
            if abs(cell.area - perfect_cell_area) > perfect_cell_area * area_tolerance:
                continue

            # Check that cell doesn't touch the polygon boundary
            if cell.touches(rotated_polygon.boundary):
                continue

            # This is a valid interior cell
            cell_count += 1
            total_area += cell.area

    return cell_count, total_area


def _get_polygon_edges(
    polygon: Polygon, min_length_ratio: float = 0.2
) -> list[tuple[float, float, float, float, float]]:
    """
    Extract longer edges from polygon for grid snapping.

    Args:
        polygon: Shapely polygon
        min_length_ratio: Minimum edge length as ratio of longest edge (default: 0.2)

    Returns:
        List of tuples: (x1, y1, x2, y2, length) for edges longer than threshold
        Sorted by length descending
    """
    import math

    if not hasattr(polygon, "exterior"):
        return []

    coords = list(polygon.exterior.coords)[:-1]  # Remove duplicate last point
    if len(coords) < 2:
        return []

    edges = []
    for i in range(len(coords)):
        p1 = coords[i]
        p2 = coords[(i + 1) % len(coords)]

        x1, y1 = p1
        x2, y2 = p2

        # Calculate edge length
        length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        edges.append((x1, y1, x2, y2, length))

    # Find longest edge
    max_length = max(edge[4] for edge in edges) if edges else 0
    if max_length == 0:
        return []

    # Filter edges longer than threshold
    threshold = max_length * min_length_ratio
    long_edges = [edge for edge in edges if edge[4] >= threshold]

    # Sort by length descending
    long_edges.sort(key=lambda e: e[4], reverse=True)

    return long_edges


def _find_best_grid_snap(
    rotated_polygon: Polygon, cell_width: float, cell_depth: float
) -> tuple[float, float]:
    """
    Find the best grid translation offset to align grid with polygon edges.

    After rotation, this shifts the grid to maximize alignment with polygon edges.

    Args:
        rotated_polygon: Polygon already rotated to desired angle
        cell_width: Width of grid cells
        cell_depth: Depth of grid cells

    Returns:
        Tuple of (offset_x, offset_y) to apply to grid origin
    """
    if not hasattr(rotated_polygon, "exterior"):
        return (0.0, 0.0)

    coords = list(rotated_polygon.exterior.coords)
    if len(coords) < 2:
        return (0.0, 0.0)

    # Get bounding box
    minx, miny, maxx, maxy = rotated_polygon.bounds

    # Find all horizontal and vertical edges in the rotated polygon
    horizontal_edges = []
    vertical_edges = []

    for i in range(len(coords) - 1):
        x1, y1 = coords[i]
        x2, y2 = coords[i + 1]

        dx = abs(x2 - x1)
        dy = abs(y2 - y1)

        # Consider edges that are mostly horizontal (within 5 degrees)
        if dx > dy * 10:  # tan(~5.7°) ≈ 0.1
            horizontal_edges.append((y1 + y2) / 2)  # Y coordinate of edge

        # Consider edges that are mostly vertical
        if dy > dx * 10:
            vertical_edges.append((x1 + x2) / 2)  # X coordinate of edge

    # Find best X offset to align with vertical edges
    offset_x = 0.0
    if vertical_edges:
        # Try to align grid with the leftmost significant vertical edge
        target_x = min(vertical_edges)
        # Calculate offset so that target_x falls on a grid line
        offset_x = (target_x - minx) % cell_width

    # Find best Y offset to align with horizontal edges
    offset_y = 0.0
    if horizontal_edges:
        # Try to align grid with the bottommost significant horizontal edge
        target_y = min(horizontal_edges)
        # Calculate offset so that target_y falls on a grid line
        offset_y = (target_y - miny) % cell_depth

    return (offset_x, offset_y)


def find_optimal_grid_rotation(
    polygon: Polygon,
    cell_width: float,
    cell_depth: float,
    angle_step: float = 5.0,
    use_ternary_search: bool = True,
) -> tuple[float, float, float]:
    """
    Find the rotation angle and snapping offset that maximizes grid cells fitting inside a polygon.

    Args:
        polygon: Shapely polygon to fit grid into
        cell_width: Width of each grid cell
        cell_depth: Depth (height) of each grid cell
        angle_step: Step size for testing angles in degrees (default: 5.0)
                   Only used if use_ternary_search is False
        use_ternary_search: Use ternary search optimization (default: True)

    Returns:
        Tuple of (rotation_angle, snap_offset_x, snap_offset_y) in degrees and units
    """
    from shapely import affinity

    # Get polygon centroid for rotation
    centroid = polygon.centroid
    cx, cy = centroid.x, centroid.y

    # Step 1: Find optimal rotation angle (without snapping)
    best_angle = 0.0

    if use_ternary_search:
        # Use ternary search to find optimal angle efficiently
        left = 0.0
        right = 180.0
        precision = 0.5

        while right - left > precision:
            mid1 = left + (right - left) / 3
            mid2 = right - (right - left) / 3

            count1, area1 = _evaluate_grid_at_angle(polygon, mid1, cell_width, cell_depth, cx, cy)
            count2, area2 = _evaluate_grid_at_angle(polygon, mid2, cell_width, cell_depth, cx, cy)

            score1 = count1 + area1 / 1e6
            score2 = count2 + area2 / 1e6

            if score1 < score2:
                left = mid1
            else:
                right = mid2

        best_angle = (left + right) / 2
        best_angle = round(best_angle, 1)

    else:
        # Linear search
        best_count = 0
        best_coverage = 0.0
        angles_to_test = [i * angle_step for i in range(int(180 / angle_step) + 1)]

        for angle in angles_to_test:
            cell_count, total_area = _evaluate_grid_at_angle(
                polygon, angle, cell_width, cell_depth, cx, cy
            )

            if cell_count > best_count or (cell_count == best_count and total_area > best_coverage):
                best_count = cell_count
                best_coverage = total_area
                best_angle = angle

    # Step 2: Find best snapping offset for the optimal angle
    rotated_poly = affinity.rotate(polygon, -best_angle, origin=(cx, cy))
    snap_offset_x, snap_offset_y = _find_best_grid_snap(rotated_poly, cell_width, cell_depth)

    return (best_angle, snap_offset_x, snap_offset_y)


def create_grid_for_polygons(
    gpkg_path: str,
    polygon_layer_name: str,
    output_layer_name: str,
    cell_width: float,
    cell_depth: float,
    optimize_rotation: bool = False,
    rotation_angle_step: float = 5.0,
    use_ternary_search: bool = True,
    clip_to_boundary: bool = False,
    tolerance_area_ratio: float = 0.0,
    tolerance_boundary_distance: float = 0.0,
) -> str:
    """
    Create a grid of rectangles for each polygon in the input layer.

    Args:
        gpkg_path: Path to the GeoPackage
        polygon_layer_name: Name of the layer containing polygons
        output_layer_name: Name for the output grid layer
        cell_width: Width of each grid cell
        cell_depth: Depth (height) of each grid cell
        optimize_rotation: If True, find optimal rotation angle for each polygon
        rotation_angle_step: Step size for testing rotation angles (default: 5.0 degrees)
                            Only used if use_ternary_search is False
        use_ternary_search: Use ternary search for faster optimization (default: True)
        clip_to_boundary: If True, include cells that intersect boundary and clip them.
                         If False (default), only include perfect interior cells.
        tolerance_area_ratio: Allow cells with area ratio >= this value (0.0-1.0).
                             E.g., 0.95 allows cells with ≥95% of perfect area.
                             Default: 0.0 (perfect only)
        tolerance_boundary_distance: Allow cells within this distance of boundary (meters).
                                     E.g., 0.5 allows cells up to 0.5m inside boundary.
                                     Default: 0.0 (no touching)

    Returns:
        Name of the created output layer
    """
    from shapely import affinity
    from shapely.geometry import box

    ds = ogr.Open(gpkg_path, 1)
    polygon_layer = ds.GetLayerByName(polygon_layer_name)
    srs = polygon_layer.GetSpatialRef()

    # Delete existing output layer if it exists
    _ogr_delete_layer(ds, output_layer_name)

    # Create output layer
    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)
    output_layer.CreateField(ogr.FieldDefn("grid_id", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("site_id", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("row", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("col", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("rotation", ogr.OFTReal))

    feature_defn = output_layer.GetLayerDefn()
    grid_id = 1

    # Process each polygon
    polygon_layer.ResetReading()
    for site_feature in polygon_layer:
        site_id = site_feature.GetFID()
        ogr_geom = site_feature.GetGeometryRef()

        if ogr_geom is None:
            continue

        # Convert to shapely
        wkb_data = ogr_geom.ExportToWkb()
        if wkb_data is None or isinstance(wkb_data, int):
            continue
        if isinstance(wkb_data, bytearray):
            wkb_data = bytes(wkb_data)

        site_geom = shapely_wkb.loads(wkb_data)

        rotation_angle = 0.0
        snap_offset_x = 0.0
        snap_offset_y = 0.0
        if optimize_rotation:
            rotation_angle, snap_offset_x, snap_offset_y = find_optimal_grid_rotation(
                site_geom, cell_width, cell_depth, rotation_angle_step, use_ternary_search
            )
            print(
                f"Site {site_id}: Optimal rotation = {rotation_angle:.1f}°,"
                f"snap offset = ({snap_offset_x:.2f}, {snap_offset_y:.2f})"
            )

        centroid = site_geom.centroid
        cx, cy = centroid.x, centroid.y

        working_geom = site_geom
        if rotation_angle != 0.0:
            working_geom = affinity.rotate(site_geom, -rotation_angle, origin=(cx, cy))

        minx, miny, maxx, maxy = working_geom.bounds

        grid_minx = minx + snap_offset_x
        grid_miny = miny + snap_offset_y

        num_cols = max(1, int(math.ceil((maxx - grid_minx) / cell_width)))
        num_rows = max(1, int(math.ceil((maxy - grid_miny) / cell_depth)))

        perfect_cell_area = cell_width * cell_depth
        area_tolerance = 0.01  # 1% tolerance for floating point errors

        for row in range(num_rows):
            for col in range(num_cols):
                cell_minx = grid_minx + col * cell_width
                cell_miny = grid_miny + row * cell_depth
                cell_maxx = cell_minx + cell_width
                cell_maxy = cell_miny + cell_depth

                cell = box(cell_minx, cell_miny, cell_maxx, cell_maxy)

                if clip_to_boundary:
                    # Mode: Include cells that intersect boundary, but clip them
                    if not working_geom.intersects(cell):
                        continue

                    # Clip cell to polygon boundary
                    clipped_cell = working_geom.intersection(cell)

                    # Skip if intersection is empty or not a polygon
                    if clipped_cell.is_empty:
                        continue
                    if clipped_cell.geom_type not in ["Polygon", "MultiPolygon"]:
                        continue

                    # Use the clipped cell as output
                    output_cell = clipped_cell
                    if rotation_angle != 0.0:
                        output_cell = affinity.rotate(clipped_cell, rotation_angle, origin=(cx, cy))

                else:
                    # Mode: Include cells that meet tolerance thresholds

                    # Check 1: Cell must be mostly within the polygon
                    # Calculate intersection area to determine overlap
                    if not working_geom.intersects(cell):
                        continue

                    intersection = working_geom.intersection(cell)
                    if intersection.is_empty:
                        continue

                    # Check 2: Area tolerance - cell must retain enough area
                    intersection_area = intersection.area
                    area_ratio = intersection_area / perfect_cell_area

                    # If tolerance is 0.0, require perfect area (within 1% for float errors)
                    # If tolerance > 0.0, allow cells with area_ratio >= tolerance_area_ratio
                    if tolerance_area_ratio > 0.0:
                        if area_ratio < tolerance_area_ratio:
                            continue
                    else:
                        # Strict mode: require perfect area
                        if abs(cell.area - perfect_cell_area) > perfect_cell_area * area_tolerance:
                            continue
                        if not working_geom.contains(cell):
                            continue

                    # Check 3: Boundary distance tolerance
                    # If tolerance_boundary_distance is 0.0, don't allow touching
                    # If > 0.0, allow cells whose edges are within tolerance
                    if tolerance_boundary_distance > 0.0:
                        if not working_geom.contains(cell.centroid):
                            continue
                    else:
                        # Cell must not touch boundary
                        if cell.touches(working_geom.boundary):
                            continue

                    if tolerance_area_ratio > 0.0 and area_ratio < 0.999:
                        # Cell is partially outside, use the clipped version
                        output_cell = intersection
                        if rotation_angle != 0.0:
                            output_cell = affinity.rotate(
                                intersection, rotation_angle, origin=(cx, cy)
                            )
                    else:
                        # Cell is perfect or nearly perfect, use original
                        output_cell = cell
                        if rotation_angle != 0.0:
                            output_cell = affinity.rotate(cell, rotation_angle, origin=(cx, cy))

                ogr_cell = ogr.CreateGeometryFromWkb(output_cell.wkb)
                out_feature = ogr.Feature(feature_defn)
                out_feature.SetGeometry(ogr_cell)
                out_feature.SetField("grid_id", grid_id)
                out_feature.SetField("site_id", site_id)
                out_feature.SetField("row", row)
                out_feature.SetField("col", col)
                out_feature.SetField("rotation", rotation_angle)
                output_layer.CreateFeature(out_feature)
                out_feature = None
                grid_id += 1

    ds = None
    return output_layer_name


def clip_site_by_roads(
    gpkg_path: str,
    site_layer_name: str,
    roads_layer_name: str,
    output_layer_name: str,
    arterial_width: float,
    secondary_width: float,
    local_width: float,
) -> str:
    """Subtract buffered roads from the site polygon."""

    ds = ogr.Open(gpkg_path, 1)
    site_layer = ds.GetLayerByName(site_layer_name)
    roads_layer = ds.GetLayerByName(roads_layer_name)
    srs = site_layer.GetSpatialRef()

    site_geoms: list[Polygon] = []
    site_layer.ResetReading()
    for feature in site_layer:
        ogr_geom = feature.GetGeometryRef()
        if ogr_geom is None:
            continue
        wkb_data = ogr_geom.ExportToWkb()
        if wkb_data is None or isinstance(wkb_data, int):
            print(f"Warning: Failed to export geometry to WKB for site feature {feature.GetFID()}")
            continue
        # Convert bytearray to bytes if needed
        if isinstance(wkb_data, bytearray):
            wkb_data = bytes(wkb_data)
        geom = shapely_wkb.loads(wkb_data)
        if not geom.is_valid:
            geom = geom.buffer(0)
        site_geoms.append(geom)

    site_union = shapely_ops.unary_union(site_geoms) if site_geoms else None

    road_buffers = []
    roads_layer.ResetReading()
    for feature in roads_layer:
        ogr_geom = feature.GetGeometryRef()
        if ogr_geom is None:
            continue
        wkb_data = ogr_geom.ExportToWkb()
        if wkb_data is None or isinstance(wkb_data, int):
            print(f"Warning: Failed to export geometry to WKB for road feature {feature.GetFID()}")
            continue
        # Convert bytearray to bytes if needed
        if isinstance(wkb_data, bytearray):
            wkb_data = bytes(wkb_data)
        road_geom = shapely_wkb.loads(wkb_data)
        # Try to get road type from various possible field names
        road_type = ""
        for field_name in ["road_type", "type", "TYPE", "ROAD_TYPE"]:
            try:
                value = feature.GetField(field_name)
                if value:
                    road_type = str(value).lower()
                    break
            except KeyError:
                continue
        width = local_width
        if road_type == "road_art":
            width = arterial_width
        elif road_type == "road_sec":
            width = secondary_width

        buffer_dist = max(width, 0.0) / 2.0
        buffered = road_geom.buffer(buffer_dist, cap_style=2, join_style=2)
        road_buffers.append(buffered)

    roads_union = shapely_ops.unary_union(road_buffers) if road_buffers else None

    if site_union is None:
        raise RuntimeError("Site layer does not contain any polygon geometry.")

    clipped = site_union
    if roads_union and not roads_union.is_empty:
        clipped = site_union.difference(roads_union)
        if clipped.is_empty:
            # Fall back to original site if subtraction removed everything
            clipped = site_union

    driver = ogr.GetDriverByName("GPKG")
    output_ds = driver.Open(gpkg_path, 1)
    _ogr_delete_layer(output_ds, output_layer_name)

    output_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)
    id_field = ogr.FieldDefn("site_id", ogr.OFTInteger)
    output_layer.CreateField(id_field)

    feature_defn = output_layer.GetLayerDefn()
    feature_id = 1
    for poly in _iter_polygons(clipped):
        ogr_geom = ogr.CreateGeometryFromWkb(poly.wkb)
        out_feature = ogr.Feature(feature_defn)
        out_feature.SetGeometry(ogr_geom)
        out_feature.SetField("site_id", feature_id)
        output_layer.CreateFeature(out_feature)
        out_feature = None
        feature_id += 1

    output_ds = None
    ds = None

    return output_layer_name


def merge_lines_with_site_boundary(
    gpkg_path: str,
    output_layer_name: str,
    site_layer_name: str,
    line_layer_names: Sequence[str],
) -> str:
    """Combine perpendicular lines, road centerlines, and site boundary lines."""

    ds = ogr.Open(gpkg_path, 1)
    base_layer = ds.GetLayerByName(line_layer_names[0])
    srs = base_layer.GetSpatialRef()

    _ogr_delete_layer(ds, output_layer_name)
    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbLineString)
    out_defn = output_layer.GetLayerDefn()

    def _write_geometry(geometry):
        geom_type = ogr.GT_Flatten(geometry.GetGeometryType())
        if geom_type == ogr.wkbLineString:
            out_feature = ogr.Feature(out_defn)
            out_feature.SetGeometry(geometry.Clone())
            output_layer.CreateFeature(out_feature)
            out_feature = None
        elif geom_type == ogr.wkbMultiLineString:
            for idx in range(geometry.GetGeometryCount()):
                _write_geometry(geometry.GetGeometryRef(idx))

    for name in line_layer_names:
        layer = ds.GetLayerByName(name)
        if layer is None:
            continue
        for feature in layer:
            geom = feature.GetGeometryRef()
            _write_geometry(geom)

    site_layer = ds.GetLayerByName(site_layer_name)
    for feature in site_layer:
        boundary = feature.GetGeometryRef().Boundary()
        _write_geometry(boundary)

    ds = None
    return output_layer_name


def polygonize_lines_to_blocks(
    gpkg_path: str,
    lines_layer_name: str,
    output_layer_name: str,
    block_type: str,
    min_area: float = 1.0,
) -> str:
    """Polygonize merged lines and persist polygons."""

    ds = ogr.Open(gpkg_path, 1)
    lines_layer = ds.GetLayerByName(lines_layer_name)
    srs = lines_layer.GetSpatialRef()

    geoms = []
    for feature in lines_layer:
        geom = shapely_wkb.loads(feature.GetGeometryRef().ExportToWkb())
        geoms.append(geom)

    polygons = list(shapely_ops.polygonize(geoms))

    _ogr_delete_layer(ds, output_layer_name)
    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)
    output_layer.CreateField(ogr.FieldDefn("block_id", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("block_type", ogr.OFTString))
    output_layer.CreateField(ogr.FieldDefn("area", ogr.OFTReal))

    feature_defn = output_layer.GetLayerDefn()
    block_id = 1
    for polygon in polygons:
        cleaned = polygon.buffer(0)
        if cleaned.is_empty or cleaned.area < min_area:
            continue
        ogr_geom = ogr.CreateGeometryFromWkb(cleaned.wkb)
        out_feature = ogr.Feature(feature_defn)
        out_feature.SetGeometry(ogr_geom)
        out_feature.SetField("block_id", block_id)
        out_feature.SetField("block_type", block_type)
        out_feature.SetField("area", cleaned.area)
        output_layer.CreateFeature(out_feature)
        out_feature = None
        block_id += 1

    ds = None
    return output_layer_name


def _compute_block_metrics(
    block: BlockRecord, preferred_width: float, preferred_depth: float
) -> None:
    """Populate derived metrics and quality scores."""

    geom = block.geom
    if not geom.is_valid:
        geom = geom.buffer(0)
        block.geom = geom

    area = geom.area
    perimeter = geom.length
    block.area = area

    vertices = len(geom.exterior.coords) - 1 if geom.exterior else 0
    block.vertex_count = max(vertices, 0)

    if vertices <= 3:
        vertex_score = 0.4
    elif vertices == 4:
        vertex_score = 1.0
    else:
        vertex_score = max(0.0, 1.0 - (vertices - 4) * 0.1)

    compactness = (4.0 * math.pi * area / (perimeter * perimeter)) if perimeter > 0 else 0.0
    block.compactness = compactness
    compactness_score = min(compactness / 0.785, 1.0) if compactness > 0 else 0.0

    minx, miny, maxx, maxy = geom.bounds
    width = maxx - minx
    height = maxy - miny
    if max(width, height) == 0:
        aspect_ratio = 0.0
    else:
        aspect_ratio = min(width, height) / max(width, height)
    block.aspect_ratio = aspect_ratio
    aspect_score = aspect_ratio

    bbox_area = width * height
    extent_ratio = min(area / bbox_area, 1.0) if bbox_area > 0 else 0.0
    block.extent_ratio = max(extent_ratio, 0.0)

    target_area = preferred_width * preferred_depth
    if target_area <= 0:
        size_ratio = 1.0
    else:
        size_ratio = area / target_area
    block.size_ratio = size_ratio

    if target_area <= 0:
        size_score = 1.0
    elif 0.5 <= size_ratio <= 1.5:
        size_score = 1.0
    else:
        size_score = max(0.0, 1.0 - (abs(size_ratio - 1.0) / 2.0))

    quality = (
        vertex_score * 0.25
        + compactness_score * 0.25
        + aspect_score * 0.20
        + block.extent_ratio * 0.15
        + size_score * 0.15
    )
    block.quality = max(0.0, min(quality, 1.0))


def _bounds_overlap(a: Sequence[float], b: Sequence[float]) -> bool:
    """Return True if bounding boxes intersect."""

    return not (a[2] < b[0] or a[0] > b[2] or a[3] < b[1] or a[1] > b[3])


def resolve_block_conflicts(
    gpkg_path: str,
    arterial_layer_name: str,
    secondary_layer_name: str,
    output_layer_name: str,
    preferred_width: float,
    preferred_depth: float,
) -> str:
    """Resolve overlaps between arterial and secondary blocks and persist winners."""

    ds = ogr.Open(gpkg_path, 1)
    arterial_layer = ds.GetLayerByName(arterial_layer_name)
    secondary_layer = ds.GetLayerByName(secondary_layer_name)
    srs = arterial_layer.GetSpatialRef()

    blocks: list[BlockRecord] = []
    block_id = 1
    for layer, block_type in ((arterial_layer, "arterial"), (secondary_layer, "secondary")):
        if layer is None:
            continue
        for feature in layer:
            geom = shapely_wkb.loads(feature.GetGeometryRef().ExportToWkb())
            if geom.is_empty or geom.area <= 0:
                continue
            record = BlockRecord(
                block_id=block_id, block_type=block_type, geom=geom, area=geom.area
            )
            _compute_block_metrics(record, preferred_width, preferred_depth)
            blocks.append(record)
            block_id += 1

    # Count neighbors using bounding-box prefiltering
    for idx, block in enumerate(blocks):
        count = 0
        for jdx, other in enumerate(blocks):
            if idx == jdx:
                continue
            if not _bounds_overlap(block.geom.bounds, other.geom.bounds):
                continue
            if block.geom.touches(other.geom) or block.geom.intersects(other.geom):
                count += 1
        block.neighbors = count

    max_neighbors = max((b.neighbors for b in blocks), default=0)

    for idx, block in enumerate(blocks):
        if not block.keep:
            continue
        for jdx in range(idx + 1, len(blocks)):
            other = blocks[jdx]
            if not other.keep or block.block_type == other.block_type:
                continue
            if not _bounds_overlap(block.geom.bounds, other.geom.bounds):
                continue
            if not (block.geom.intersects(other.geom) or block.geom.touches(other.geom)):
                continue

            score_block = block.quality * 0.7 + (
                (block.neighbors / max_neighbors) * 0.3 if max_neighbors else 0.0
            )
            score_other = other.quality * 0.7 + (
                (other.neighbors / max_neighbors) * 0.3 if max_neighbors else 0.0
            )

            if score_block >= score_other:
                other.keep = False
            else:
                block.keep = False
                break

    driver = ogr.GetDriverByName("GPKG")
    output_ds = driver.Open(gpkg_path, 1)
    _ogr_delete_layer(output_ds, output_layer_name)
    output_layer = output_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)

    output_layer.CreateField(ogr.FieldDefn("block_id", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("block_type", ogr.OFTString))
    output_layer.CreateField(ogr.FieldDefn("area", ogr.OFTReal))
    output_layer.CreateField(ogr.FieldDefn("quality", ogr.OFTReal))
    output_layer.CreateField(ogr.FieldDefn("neighbors", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("vertices", ogr.OFTInteger))
    output_layer.CreateField(ogr.FieldDefn("compactness", ogr.OFTReal))
    output_layer.CreateField(ogr.FieldDefn("aspect_ratio", ogr.OFTReal))
    output_layer.CreateField(ogr.FieldDefn("size_ratio", ogr.OFTReal))

    feature_defn = output_layer.GetLayerDefn()
    for record in blocks:
        if not record.keep:
            continue
        ogr_geom = ogr.CreateGeometryFromWkb(record.geom.buffer(0).wkb)
        out_feature = ogr.Feature(feature_defn)
        out_feature.SetGeometry(ogr_geom)
        out_feature.SetField("block_id", record.block_id)
        out_feature.SetField("block_type", record.block_type)
        out_feature.SetField("area", record.area)
        out_feature.SetField("quality", record.quality)
        out_feature.SetField("neighbors", record.neighbors)
        out_feature.SetField("vertices", record.vertex_count)
        out_feature.SetField("compactness", record.compactness)
        out_feature.SetField("aspect_ratio", record.aspect_ratio)
        out_feature.SetField("size_ratio", record.size_ratio)
        output_layer.CreateFeature(out_feature)
        out_feature = None

    output_ds = None
    ds = None

    return output_layer_name
