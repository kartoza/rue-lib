import math
from pathlib import Path
from typing import Optional

import geopandas as gpd
from shapely import LineString, Point, unary_union
from shapely.errors import TopologicalError

from rue_lib.core.geometry import buffer_layer

from .operations import clip_layer


def merge_setback_layers(
    output_path: Path,
    arterial_layer: str = "10a_arterial_setback_clipped",
    secondary_layer: str = "10b_secondary_setback_clipped",
    output_layer_name: str = "10c_setback_all",
) -> str:
    """Merge arterial and secondary setback layers into a single, non-overlapping layer.

    The result is a planar union of both layers – no overlapping polygons.

    Args:
        output_path: Path to the GeoPackage.
        arterial_layer: Name of the arterial setback layer to merge.
        secondary_layer: Name of the secondary setback layer to merge.
        output_layer_name: Name of the merged output layer.

    Returns:
        The name of the merged output layer.
    """
    print("\nMerging setback layers (planar union, no overlaps)...")
    print(f"  Arterial setback layer:  {arterial_layer}")
    print(f"  Secondary setback layer: {secondary_layer}")
    print(f"  Output merged layer:     {output_layer_name}")

    # Read input layers
    gdf_art = gpd.read_file(output_path, layer=arterial_layer)
    gdf_sec = gpd.read_file(output_path, layer=secondary_layer)

    print(f"  Loaded {len(gdf_art)} arterial setback feature(s)")
    print(f"  Loaded {len(gdf_sec)} secondary setback feature(s)")

    if gdf_art.empty and gdf_sec.empty:
        print("  Warning: both input setback layers are empty; nothing to merge")
        gdf_empty = gpd.GeoDataFrame(geometry=[], crs=None)
        gdf_empty.to_file(output_path, layer=output_layer_name, driver="GPKG")
        return output_layer_name

    frames = []
    if not gdf_art.empty:
        frames.append(gdf_art)
    if not gdf_sec.empty:
        frames.append(gdf_sec)

    gdf_all = gpd.GeoDataFrame(
        gpd.pd.concat(frames, ignore_index=True),
        crs=gdf_art.crs or gdf_sec.crs,
    )

    combined_geom = unary_union(gdf_all.geometry)

    if combined_geom.is_empty:
        print("  Warning: union geometry is empty after merge")
        gdf_empty = gpd.GeoDataFrame(geometry=[], crs=gdf_all.crs)
        gdf_empty.to_file(output_path, layer=output_layer_name, driver="GPKG")
        return output_layer_name

    if combined_geom.geom_type == "Polygon":
        geoms = [combined_geom]
    elif combined_geom.geom_type == "MultiPolygon":
        geoms = list(combined_geom.geoms)
    else:
        geoms = [combined_geom]

    gdf_merged = gpd.GeoDataFrame(geometry=geoms, crs=gdf_all.crs)

    print(f"  Merged setback polygons (non-overlapping): {len(gdf_merged)}")

    gdf_merged.to_file(output_path, layer=output_layer_name, driver="GPKG")
    print(f"  Saved merged setbacks layer as: {output_layer_name}")

    return output_layer_name


def apply_inner_buffer_to_cells(
    output_path: Path,
    site_boundary_layer: str,
    buffer_distance: float,
) -> str:
    """Apply inner buffer to site boundary lines and clip to site boundary.

    Args:
        output_path: Path to GeoPackage
        site_boundary_layer: Name of layer containing site boundary (polygon)
        buffer_distance: Buffer distance in meters (negative for inner buffer)

    Returns:
        Name of the buffered layer
    """
    print("\nCreating inner buffer zone from site boundary lines...")
    print(f"  Buffer distance: {buffer_distance}m")

    buffer_layer(
        str(output_path),
        "13_site_boundary_lines",
        buffer_distance,
        str(output_path),
        "14_site_boundary_buffer_temp",
        dissolve=True,
    )

    print("  Clipping buffer to site boundary...")
    clip_layer(
        output_path,
        "14_site_boundary_buffer_temp",
        output_path,
        site_boundary_layer,
        output_path,
        "14a_site_boundary_inner_buffer",
    )

    print("  Created buffered layer: 14a_site_boundary_inner_buffer")
    return "14a_site_boundary_inner_buffer"


def create_perpendicular_lines_inside_buffer_from_points(
    output_path: Path,
    points_layer: str,
    buffer_layer: str,
    line_length: float,
    output_layer_name: Optional[str] = None,
    tangent_step_fraction: float = 0.001,
) -> Optional[str]:
    """Create perpendicular lines inside the buffer, starting from points on the buffer boundary.

    The perpendicular is computed with respect to the buffer boundary:
      - For each point, compute the local tangent of the buffer boundary.
      - From the tangent, derive a normal (perpendicular) direction.
      - Choose the normal direction that goes INSIDE the buffer polygon.
      - Build a line segment from the point inward, clipped to the buffer polygon.

    Args:
        output_path: Path to the GeoPackage.
        points_layer: Name of the layer containing points on the buffer boundary.
        buffer_layer: Name of the layer containing the buffer polygon.
        line_length: Desired length of the perpendicular line (map units, e.g. meters).
        output_layer_name: Optional name for the output layer. If None, a default is used.
        tangent_step_fraction: Step along the boundary (fraction of boundary length)
                               used to approximate the local tangent.

    Returns:
        Name of the created perpendicular lines layer, or None on failure.
    """
    import geopandas as gpd

    print(f"\nCreating perpendicular lines inside buffer from points in layer: {points_layer}")

    # Load data
    gdf_pts = gpd.read_file(output_path, layer=points_layer)
    gdf_buffer = gpd.read_file(output_path, layer=buffer_layer)

    print(f"  Loaded {len(gdf_pts)} points")
    print(f"  Loaded {len(gdf_buffer)} buffer feature(s)")

    if gdf_pts.empty:
        print("  Warning: points_layer is empty")
        return None

    if gdf_buffer.empty:
        print("  Warning: buffer_layer is empty")
        return None

    # Dissolve buffer to single polygon/multipolygon and get its boundary
    buffer_geom = gdf_buffer.geometry.unary_union
    if buffer_geom is None or buffer_geom.is_empty:
        print("  Warning: buffer geometry is empty")
        return None

    try:
        buffer_boundary = buffer_geom.boundary
    except TopologicalError:
        print("  Warning: buffer geometry is invalid, attempting buffer(0) fix")
        buffer_geom = buffer_geom.buffer(0)
        buffer_boundary = buffer_geom.boundary

    boundary_len = buffer_boundary.length
    if boundary_len == 0:
        print("  Warning: buffer boundary length is zero")
        return None

    step = boundary_len * tangent_step_fraction
    if step == 0:
        step = line_length * 0.01

    if output_layer_name is None:
        output_layer_name = f"{points_layer}_perp_inside_buffer"

    half_eps_in = line_length * 0.01

    lines = []
    records = []

    for _idx, row in gdf_pts.iterrows():
        pt = row.geometry
        if pt is None or pt.is_empty:
            continue

        if not isinstance(pt, Point):
            continue

        x0, y0 = pt.x, pt.y

        s = buffer_boundary.project(Point(x0, y0))

        s0 = max(0.0, s - step)
        s1 = min(boundary_len, s + step)

        p0 = buffer_boundary.interpolate(s0)
        p1 = buffer_boundary.interpolate(s1)

        tx = p1.x - p0.x
        ty = p1.y - p0.y

        t_norm = math.hypot(tx, ty)
        if t_norm == 0:
            continue

        tx /= t_norm
        ty /= t_norm

        n1x, n1y = -ty, tx
        n2x, n2y = ty, -tx

        test1 = Point(x0 + n1x * half_eps_in, y0 + n1y * half_eps_in)
        test2 = Point(x0 + n2x * half_eps_in, y0 + n2y * half_eps_in)

        inside1 = buffer_geom.contains(test1)
        inside2 = buffer_geom.contains(test2)

        if inside1 and not inside2:
            nx, ny = n1x, n1y
        elif inside2 and not inside1:
            nx, ny = n2x, n2y
        elif inside1 and inside2:
            nx, ny = n1x, n1y
        else:
            continue

        end_x = x0 + nx * line_length
        end_y = y0 + ny * line_length

        candidate_line = LineString([(x0, y0), (end_x, end_y)])

        clipped = candidate_line.intersection(buffer_geom)
        if clipped.is_empty:
            continue

        if clipped.geom_type == "LineString":
            line_geom = clipped
        elif clipped.geom_type == "MultiLineString":
            segments = list(clipped.geoms)
            line_geom = min(segments, key=lambda seg: seg.distance(Point(x0, y0)))
        else:
            continue

        if line_geom.length == 0:
            continue

        lines.append(line_geom)
        rec = row.drop(labels="geometry").to_dict()
        records.append(rec)

    if not lines:
        print("  Warning: no perpendicular lines created")
        return None

    import geopandas as gpd

    gdf_out = gpd.GeoDataFrame(records, geometry=lines, crs=gdf_pts.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")

    print(f"  Saved perpendicular lines inside buffer to layer: {output_layer_name}")
    return output_layer_name


def create_guide_points_from_site_boundary(
    output_path: Path,
    site_boundary_lines_layer: str,
    perp_lines_layer: str,
    output_layer_name: str = "13_site_boundary_points",
) -> str:
    """Create guide points from site boundary lines and perpendicular lines.

    This function creates two types of points:
    1. Corner points where 2 site boundary lines meet
    2. Points where perpendicular lines intersect/touch site boundary lines

    Args:
        output_path: Path to the GeoPackage
        site_boundary_lines_layer: Name of layer containing site boundary lines
        perp_lines_layer: Name of layer containing perpendicular lines
        output_layer_name: Name for the output points layer

    Returns:
        Name of the created points layer
    """
    print("\nCreating guide points from site boundary...")
    print(f"  Site boundary lines: {site_boundary_lines_layer}")
    print(f"  Perpendicular lines: {perp_lines_layer}")

    # Load layers
    gdf_boundary = gpd.read_file(output_path, layer=site_boundary_lines_layer)
    gdf_perp = gpd.read_file(output_path, layer=perp_lines_layer)

    print(f"  Loaded {len(gdf_boundary)} boundary line(s)")
    print(f"  Loaded {len(gdf_perp)} perpendicular line(s)")

    points = []
    point_types = []
    print("  Extracting corner points...")

    boundary_lines = []
    for _, row in gdf_boundary.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        if geom.geom_type == "LineString":
            boundary_lines.append(geom)
        elif geom.geom_type == "MultiLineString":
            boundary_lines.extend(list(geom.geoms))

    endpoint_counts: dict[tuple[float, float], int] = {}

    for line in boundary_lines:
        coords = list(line.coords)
        if len(coords) >= 2:
            # Start and end endpoints
            for coord in (coords[0], coords[-1]):
                key = (round(coord[0], 6), round(coord[1], 6))
                endpoint_counts[key] = endpoint_counts.get(key, 0) + 1

    corner_points_set = {pt for pt, cnt in endpoint_counts.items() if cnt > 1}

    for pt_coords in corner_points_set:
        points.append(Point(pt_coords[0], pt_coords[1]))
        point_types.append("corner")

    print(f"  Found {len(corner_points_set)} non-dead-end corner points")

    print("  Finding perpendicular line intersections...")

    boundary_union = unary_union(boundary_lines)

    extension_factor = 1000.0

    intersection_count = 0
    for _, row in gdf_perp.iterrows():
        perp_geom = row.geometry
        if perp_geom is None or perp_geom.is_empty:
            continue

        try:
            from shapely.geometry import LineString

            coords = list(perp_geom.coords)
            if len(coords) < 2:
                continue

            start = coords[0]
            end = coords[-1]

            dx = end[0] - start[0]
            dy = end[1] - start[1]
            length = (dx**2 + dy**2) ** 0.5

            if length == 0:
                continue

            dx /= length
            dy /= length

            extended_start = (start[0] - dx * extension_factor, start[1] - dy * extension_factor)
            extended_end = (end[0] + dx * extension_factor, end[1] + dy * extension_factor)

            extended_line = LineString([extended_start, extended_end])

            intersection = extended_line.intersection(boundary_union)

            if intersection.is_empty:
                continue

            original_mid_x = (start[0] + end[0]) / 2
            original_mid_y = (start[1] + end[1]) / 2
            original_mid = Point(original_mid_x, original_mid_y)

            intersection_points = []
            if intersection.geom_type == "Point":
                intersection_points.append(intersection)
            elif intersection.geom_type == "MultiPoint":
                intersection_points.extend(list(intersection.geoms))
            elif intersection.geom_type == "LineString":
                coords = list(intersection.coords)
                if coords:
                    intersection_points.append(Point(coords[0]))
                    if len(coords) > 1:
                        intersection_points.append(Point(coords[-1]))
            elif intersection.geom_type == "MultiLineString":
                for line in intersection.geoms:
                    coords = list(line.coords)
                    if coords:
                        intersection_points.append(Point(coords[0]))
                        if len(coords) > 1:
                            intersection_points.append(Point(coords[-1]))

            if intersection_points:
                closest_pt = min(intersection_points, key=lambda pt: pt.distance(original_mid))
                pt_coords = (round(closest_pt.x, 6), round(closest_pt.y, 6))

                if pt_coords not in corner_points_set:
                    points.append(Point(pt_coords[0], pt_coords[1]))
                    point_types.append("perp_intersection")
                    intersection_count += 1

        except Exception as e:
            print(f"  Warning: Failed to compute intersection: {e}")
            continue

    print(f"  Found {intersection_count} perpendicular line intersection points")

    if not points:
        print("  Warning: No guide points created!")
        return output_layer_name

    gdf_points = gpd.GeoDataFrame(
        {"point_type": point_types},
        geometry=points,
        crs=gdf_boundary.crs,
    )

    gdf_points.to_file(output_path, layer=output_layer_name, driver="GPKG")

    print(f"  Created {len(points)} total guide points")
    print(f"    - {len(corner_points_set)} corner points")
    print(f"    - {intersection_count} perpendicular intersection points")
    print(f"  Saved to layer: {output_layer_name}")

    return output_layer_name


def create_perpendicular_lines_from_guide_points(
    output_path: Path,
    site_boundary_lines_layer: str,
    site_polygon_layer: str,
    guide_points_layer: str,
    line_length: float,
    output_layer_name: str = "13_site_boundary_perp_from_points",
) -> str | None:
    """Create perpendicular lines from guide points, oriented away from the site polygon.

    For each guide point:
      - Find boundary line segments that touch or are very close to the point.
      - If the point_type == "corner", create one perpendicular line for each touching segment.
      - Otherwise, create a perpendicular line from the closest boundary segment.
      - The perpendicular direction is chosen to go OUTSIDE the site polygon.
      - Any portion of the perpendicular line lying inside the polygon is removed,
        so the final line is only outside the boundary.

    Args:
        output_path: Path to the GeoPackage.
        site_boundary_lines_layer: Name of layer containing site boundary lines.
        site_polygon_layer: Name of layer containing the site polygon.
        guide_points_layer: Name of the layer containing guide points (with point_type field).
        line_length: Length of perpendicular lines (map units, e.g. meters).
        output_layer_name: Name of the output layer for perpendicular lines.

    Returns:
        Name of the created perpendicular line layer, or None if no lines are created.
    """
    print("\nCreating perpendicular lines from guide points...")
    print(f"  Boundary lines: {site_boundary_lines_layer}")
    print(f"  Site polygon:   {site_polygon_layer}")
    print(f"  Guide points:   {guide_points_layer}")

    # Load layers
    gdf_boundary = gpd.read_file(output_path, layer=site_boundary_lines_layer)
    gdf_site = gpd.read_file(output_path, layer=site_polygon_layer)
    gdf_points = gpd.read_file(output_path, layer=guide_points_layer)

    print(f"  Loaded {len(gdf_boundary)} boundary line(s)")
    print(f"  Loaded {len(gdf_site)} site polygon feature(s)")
    print(f"  Loaded {len(gdf_points)} guide point(s)")

    if gdf_boundary.empty or gdf_site.empty or gdf_points.empty:
        print("  Warning: one or more input layers are empty; no perpendiculars created")
        return None

    # Collect boundary line segments
    boundary_lines = []
    for _, row in gdf_boundary.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "LineString":
            boundary_lines.append(geom)
        elif geom.geom_type == "MultiLineString":
            boundary_lines.extend(list(geom.geoms))

    if not boundary_lines:
        print("  Warning: no boundary line segments found")
        return None

    # Site polygon geometry (to decide inside vs outside)
    site_geom = unary_union(gdf_site.geometry)
    if site_geom is None or site_geom.is_empty:
        print("  Warning: site polygon geometry is empty")
        return None

    lines = []
    records = []

    touch_tol = 1e-6  # distance to treat a point as touching a line
    tangent_step_fraction = 0.01  # fraction of line length for tangent calc
    inside_test_step = line_length * 0.01  # small step along normal to test inside/outside

    def boundary_segments_touching_point(pt: Point):
        """Return list of boundary line segments that touch/are close to the given point."""
        touching = []
        for line in boundary_lines:
            if line.distance(pt) <= touch_tol:
                touching.append(line)
        # Fallback: if nothing is within tol, take the closest line
        if not touching:
            min_d = float("inf")
            closest = None
            for line in boundary_lines:
                d = line.distance(pt)
                if d < min_d:
                    min_d = d
                    closest = line
            if closest is not None:
                touching.append(closest)
        return touching

    def perpendicular_line_from_segment_at_point(seg: LineString, pt: Point) -> LineString | None:
        """Build a perpendicular line from segment seg at pt, oriented away from site polygon."""
        if seg is None or seg.is_empty:
            return None

        # Parameter along segment
        seg_len = seg.length
        if seg_len == 0:
            return None

        s = seg.project(pt)
        step = max(seg_len * tangent_step_fraction, line_length * 0.01)

        s0 = max(0.0, s - step)
        s1 = min(seg_len, s + step)

        p0 = seg.interpolate(s0)
        p1 = seg.interpolate(s1)

        tx = p1.x - p0.x
        ty = p1.y - p0.y
        t_norm = math.hypot(tx, ty)
        if t_norm == 0:
            return None
        tx /= t_norm
        ty /= t_norm

        # Two candidate normals
        n1x, n1y = -ty, tx
        n2x, n2y = ty, -tx

        x0, y0 = pt.x, pt.y

        test1 = Point(x0 + n1x * inside_test_step, y0 + n1y * inside_test_step)
        test2 = Point(x0 + n2x * inside_test_step, y0 + n2y * inside_test_step)

        inside1 = site_geom.contains(test1)
        inside2 = site_geom.contains(test2)

        # We want the direction AWAY from the polygon, so pick the one that is NOT inside
        if inside1 and not inside2:
            nx, ny = n2x, n2y
        elif inside2 and not inside1:
            nx, ny = n1x, n1y
        elif not inside1 and not inside2:
            # Both directions considered outside (boundary may be ambiguous); just pick one
            nx, ny = n1x, n1y
        else:
            # Both directions go inside; skip
            return None

        end_x = x0 + nx * line_length
        end_y = y0 + ny * line_length

        raw_line = LineString([(x0, y0), (end_x, end_y)])

        # Remove any part of the line that lies inside the polygon, keep only outside piece(s)
        diff = raw_line.difference(site_geom)
        if diff.is_empty:
            return None

        if diff.geom_type == "LineString":
            return diff
        if diff.geom_type == "MultiLineString":
            # Return the piece that starts closest to the original point
            pieces = list(diff.geoms)
            return min(pieces, key=lambda seg2: seg2.distance(pt))

        return None

    for idx, row in gdf_points.iterrows():
        pt = row.geometry
        if pt is None or pt.is_empty or not isinstance(pt, Point):
            continue

        pt_type = row.get("point_type", "")

        touching_segments = boundary_segments_touching_point(pt)

        if not touching_segments:
            continue

        # Corner: create one perpendicular for each touching segment
        if pt_type == "corner":
            used = 0
            for seg in touching_segments:
                line = perpendicular_line_from_segment_at_point(seg, pt)
                if line is None or line.length == 0:
                    continue
                lines.append(line)
                rec = row.drop(labels="geometry").to_dict()
                rec["source_point_index"] = idx
                rec["source_point_type"] = pt_type
                records.append(rec)
                used += 1
            if used == 0:
                # Fallback: at least try with the closest segment
                seg = touching_segments[0]
                line = perpendicular_line_from_segment_at_point(seg, pt)
                if line is not None and line.length > 0:
                    lines.append(line)
                    rec = row.drop(labels="geometry").to_dict()
                    rec["source_point_index"] = idx
                    rec["source_point_type"] = pt_type
                    records.append(rec)
        else:
            # Non-corner: use only the closest touching segment
            closest_seg = None
            min_d = float("inf")
            for seg in touching_segments:
                d = seg.distance(pt)
                if d < min_d:
                    min_d = d
                    closest_seg = seg
            if closest_seg is None:
                continue
            line = perpendicular_line_from_segment_at_point(closest_seg, pt)
            if line is None or line.length == 0:
                continue
            lines.append(line)
            rec = row.drop(labels="geometry").to_dict()
            rec["source_point_index"] = idx
            rec["source_point_type"] = pt_type
            records.append(rec)

    if not lines:
        print("  Warning: no perpendicular lines created from guide points")
        return None

    gdf_out = gpd.GeoDataFrame(records, geometry=lines, crs=gdf_boundary.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")

    print(f"  Created {len(lines)} perpendicular lines from guide points")
    print(f"  Saved to layer: {output_layer_name}")

    return output_layer_name
