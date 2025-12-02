# src/rue_lib/streets/runner.py
from pathlib import Path
from typing import Optional

from osgeo import gdal, ogr

from rue_lib.core.geometry import buffer_layer, get_utm_zone_from_layer, reproject_layer
from rue_lib.streets.grids import grids_from_site

from .blocks_orthogonal import (
    clip_site_by_roads,
)
from .cell import fix_grid_cells_with_perpendicular_lines
from .config import StreetConfig
from .operations import (
    break_multipart_features,
    cleanup_grid_blocks,
    clip_layer,
    create_grid_from_on_grid,
    create_local_streets_zone,
    create_on_grid_cells_from_perpendiculars,
    erase_layer,
    export_layer_to_geojson,
    extract_by_expression,
    extract_site_boundary_lines,
    merge_grid_layers_with_type,
    merge_layers_without_overlaps,
)

gdal.UseExceptions()


def generate_on_grid_blocks(output_path: Path, cfg: StreetConfig) -> Path:
    print("Clip arterial setback to site boundary...")
    clip_layer(
        output_path,
        "06_arterial_setback",
        output_path,
        "site_clipped_by_roads",
        output_path,
        "10a_arterial_setback_clipped",
    )

    print("Clip secondary setback to site boundary...")
    clip_layer(
        output_path,
        "07_secondary_setback",
        output_path,
        "site_clipped_by_roads",
        output_path,
        "10b_secondary_setback_clipped",
    )

    print("Intersect arterial and secondary setbacks...")
    clip_layer(
        output_path,
        "10a_arterial_setback_clipped",
        output_path,
        "10b_secondary_setback_clipped",
        output_path,
        "10_intersected_setbacks",
    )

    print("Arterial setback without intersection...")
    erase_layer(
        output_path,
        "10a_arterial_setback_clipped",
        output_path,
        "10_intersected_setbacks",
        output_path,
        "11_arterial_setback_final",
    )

    print("Secondary setback without intersection...")
    erase_layer(
        output_path,
        "10b_secondary_setback_clipped",
        output_path,
        "10_intersected_setbacks",
        output_path,
        "12_secondary_setback_final",
    )

    print("Breaking arterial setback multipart features...")
    break_multipart_features(
        output_path,
        "11_arterial_setback_final",
        output_path,
        "11_arterial_setback_final",
    )

    print("Breaking secondary setback multipart features...")
    break_multipart_features(
        output_path,
        "12_secondary_setback_final",
        output_path,
        "12_secondary_setback_final",
    )

    print("Create grid from on-grid arterial setback...")
    create_grid_from_on_grid(
        output_path,
        "11_arterial_setback_final",
        "04_arterial_roads",
        cfg.off_grid_partitions_preferred_width,
        output_path,
        "11a_arterial_setback_grid",
        road_buffer_distance=cfg.road_arterial_width_m,
    )

    print("Clean up grid blocks...")
    cleanup_grid_blocks(
        output_path,
        "11a_arterial_setback_grid",
        output_path,
        "11_arterial_setback_grid_cleaned",
        cfg.arterial_setback_depth * cfg.off_grid_partitions_preferred_width * 0.5,
    )

    print("Create grid from on-grid secondary setback...")
    create_grid_from_on_grid(
        output_path,
        "12_secondary_setback_final",
        "05_secondary_roads",
        cfg.off_grid_partitions_preferred_width,
        output_path,
        "12a_secondary_setback_grid",
        intersected_setbacks_layer_name="10_intersected_setbacks",
        road_buffer_distance=cfg.road_secondary_width_m,
    )

    print("Clean up grid blocks...")
    cleanup_grid_blocks(
        output_path,
        "12a_secondary_setback_grid",
        output_path,
        "12_secondary_setback_grid_cleaned",
        cfg.secondary_setback_depth * cfg.off_grid_partitions_preferred_width * 0.5,
    )

    return output_path


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
    import geopandas as gpd
    from shapely.ops import unary_union

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

    # Step 1: Buffer the site_boundary_lines
    buffer_layer(
        str(output_path),
        "13_site_boundary_lines",
        buffer_distance,  # Use positive value, will invert direction below
        str(output_path),
        "14_site_boundary_buffer_temp",
        dissolve=True,
    )

    # Step 2: Clip the buffer to the site boundary
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


def extract_grid_lines_in_buffer(
    output_path: Path,
    grid_cells_layer: str,
    buffer_layer: str,
    coord_precision: int = 6,
) -> str | None:
    """Extract internal grid lines (shared edges between cells) that intersect the buffer boundary,
    and create points at the intersection between the grid line and the buffer boundary.

    The point moves along the grid line, so the final geometry is guaranteed to lie on BOTH:
      - the internal grid line
      - the buffer boundary line

    Attributes include:
      - line_id
      - line_length
      - dx, dy (unit direction vector of the grid line)

    Args:
        output_path: Path to GeoPackage
        grid_cells_layer: Name of layer containing grid cells
        buffer_layer: Name of layer containing the buffer zone (polygon)
        coord_precision: Decimal places for rounding coordinates when matching shared edges.

    Returns:
        Name of the extracted points layer, or None if no lines are found.
    """
    import math

    import geopandas as gpd
    from shapely.errors import TopologicalError
    from shapely.geometry import LineString, Point
    from shapely.ops import nearest_points

    print("\nExtracting internal grid lines on buffer boundary...")
    print(f"  Grid layer: {grid_cells_layer}")
    print(f"  Buffer layer: {buffer_layer}")

    # Read the grid cells and buffer zone
    gdf_cells = gpd.read_file(output_path, layer=grid_cells_layer)
    gdf_buffer = gpd.read_file(output_path, layer=buffer_layer)

    print(f"  Loaded {len(gdf_cells)} grid cells")
    print(f"  Loaded {len(gdf_buffer)} buffer feature(s)")

    if gdf_cells.empty:
        print("  Warning: grid_cells_layer is empty")
        return None

    if gdf_buffer.empty:
        print("  Warning: buffer_layer is empty")
        return None

    # Dissolve buffer to a single geometry and take its boundary
    buffer_geom = gdf_buffer.geometry.unary_union
    buffer_boundary = buffer_geom.boundary  # lines along which we want the intersections

    # Helpers to normalize points/edges with rounding
    def norm_point(pt):
        return (
            round(pt[0], coord_precision),
            round(pt[1], coord_precision),
        )

    def norm_edge(p1, p2):
        a = norm_point(p1)
        b = norm_point(p2)
        return (a, b) if a <= b else (b, a)

    edge_dict: dict[tuple, dict] = {}
    total_edges = 0

    for _idx, row in gdf_cells.iterrows():
        geom = row.geometry

        if geom is None or geom.is_empty:
            continue

        if not geom.is_valid:
            try:
                geom = geom.buffer(0)
            except TopologicalError:
                continue

        boundary = geom.boundary

        if boundary.geom_type == "LineString":
            lines = [boundary]
        elif boundary.geom_type == "MultiLineString":
            lines = list(boundary.geoms)
        else:
            continue

        for line in lines:
            coords = list(line.coords)
            for i in range(len(coords) - 1):
                p1, p2 = coords[i], coords[i + 1]
                if p1 == p2:
                    continue

                edge_key = norm_edge(p1, p2)
                total_edges += 1

                if edge_key not in edge_dict:
                    edge_dict[edge_key] = {
                        "count": 0,
                        "geom": LineString([p1, p2]),
                    }
                edge_dict[edge_key]["count"] += 1

    print(f"  Processed ~{total_edges} raw edges")
    print(f"  Unique normalized edges: {len(edge_dict)}")

    # Collect only shared edges that intersect the BUFFER BOUNDARY (not just buffer polygon)
    internal_lines: list[LineString] = []
    shared_edges = 0

    for _edge_key, info in edge_dict.items():
        if info["count"] > 1:
            edge = info["geom"]
            if buffer_boundary.intersects(edge):
                shared_edges += 1
                internal_lines.append(edge)

    print(f"  Shared edges (count > 1) intersecting buffer boundary: {shared_edges}")

    if len(internal_lines) == 0:
        print("  Warning: No internal lines intersect the buffer boundary!")
        return None

    # Helper to pick a single Point from intersection geometry
    def pick_intersection_point(edge: LineString, midpoint: Point):
        inter = edge.intersection(buffer_boundary)
        if inter.is_empty:
            # Fallback: project midpoint onto boundary
            _, snapped = nearest_points(midpoint, buffer_boundary)
            return snapped

        if inter.geom_type == "Point":
            return inter

        if inter.geom_type == "MultiPoint":
            points = list(inter.geoms)
            return min(points, key=lambda p: p.distance(midpoint))

        if inter.geom_type in ("LineString", "MultiLineString"):
            # Overlap: choose nearest point on overlap to the midpoint
            _, snapped = nearest_points(midpoint, inter)
            return snapped

        if inter.geom_type == "GeometryCollection":
            pts = [g for g in inter.geoms if g.geom_type == "Point"]
            if pts:
                return min(pts, key=lambda p: p.distance(midpoint))
            # fallback
            _, snapped = nearest_points(midpoint, inter)
            return snapped

        # Last resort
        _, snapped = nearest_points(midpoint, buffer_boundary)
        return snapped

    # Build points (intersection) + attributes, including line direction dx, dy
    points = []
    attrs = []

    for line_idx, edge in enumerate(internal_lines):
        # Midpoint purely used to pick the closest intersection when multiple
        midpoint = edge.interpolate(0.5, normalized=True)
        inter_point = pick_intersection_point(edge, midpoint)

        if not isinstance(inter_point, Point) or inter_point.is_empty:
            continue

        # Direction of the grid line (from first to last coordinate)
        x1, y1, z1 = edge.coords[0]
        x2, y2, z2 = edge.coords[-1]
        vx, vy = x2 - x1, y2 - y1
        length = math.hypot(vx, vy)
        if length == 0:
            continue

        dx = vx / length
        dy = vy / length

        points.append(inter_point)
        attrs.append(
            {
                "line_id": line_idx,
                "line_length": float(edge.length),
                "dx": dx,
                "dy": dy,
            }
        )

    print(f"  Created {len(points)} points at grid/buffer intersections")

    if not points:
        print("  Warning: No intersection points could be created!")
        return None

    # Create GeoDataFrame with points as geometry
    gdf_points = gpd.GeoDataFrame(
        attrs,
        geometry=points,
        crs=gdf_cells.crs,
    )

    output_layer_name = f"{grid_cells_layer}_internal_points_in_buffer"
    gdf_points.to_file(output_path, layer=output_layer_name, driver="GPKG")

    print(f"  Saved intersection points to layer: {output_layer_name}")
    return output_layer_name


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
    import math

    import geopandas as gpd
    from shapely.errors import TopologicalError
    from shapely.geometry import LineString, Point

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

    # Step along the boundary to approximate the tangent
    step = boundary_len * tangent_step_fraction
    if step == 0:
        step = line_length * 0.01  # small fallback

    if output_layer_name is None:
        output_layer_name = f"{points_layer}_perp_inside_buffer"

    half_eps_in = line_length * 0.01  # small step to test inside/outside

    lines = []
    records = []

    for _idx, row in gdf_pts.iterrows():
        pt = row.geometry
        if pt is None or pt.is_empty:
            continue

        if not isinstance(pt, Point):
            continue

        # Use only 2D coords (avoid mixing 2D/3D)
        x0, y0 = pt.x, pt.y

        # Position along buffer boundary (2D)
        s = buffer_boundary.project(Point(x0, y0))

        # Approximate local tangent by sampling slightly before and after s
        s0 = max(0.0, s - step)
        s1 = min(boundary_len, s + step)

        p0 = buffer_boundary.interpolate(s0)
        p1 = buffer_boundary.interpolate(s1)

        tx = p1.x - p0.x
        ty = p1.y - p0.y

        # If tangent is degenerate, skip
        t_norm = math.hypot(tx, ty)
        if t_norm == 0:
            continue

        tx /= t_norm
        ty /= t_norm

        # Two possible normals (perpendicular to tangent)
        n1x, n1y = -ty, tx
        n2x, n2y = ty, -tx

        # Choose the normal that goes INSIDE the buffer
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
            # Choose the segment that starts closest to the original point
            segments = list(clipped.geoms)
            line_geom = min(segments, key=lambda seg: seg.distance(Point(x0, y0)))
        else:
            continue

        if line_geom.length == 0:
            continue

        lines.append(line_geom)
        # Copy attributes except geometry
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
    import geopandas as gpd
    from shapely.geometry import Point
    from shapely.ops import unary_union

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

    # Count how many boundary segments share each endpoint
    endpoint_counts: dict[tuple[float, float], int] = {}

    for line in boundary_lines:
        coords = list(line.coords)
        if len(coords) >= 2:
            # Start and end endpoints
            for coord in (coords[0], coords[-1]):
                key = (round(coord[0], 6), round(coord[1], 6))
                endpoint_counts[key] = endpoint_counts.get(key, 0) + 1

    # Keep only non–dead-end corners: endpoints used by at least two segments
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
                # Take endpoints of all line segments
                for line in intersection.geoms:
                    coords = list(line.coords)
                    if coords:
                        intersection_points.append(Point(coords[0]))
                        if len(coords) > 1:
                            intersection_points.append(Point(coords[-1]))

            # Find the closest intersection point to the original perpendicular line midpoint
            if intersection_points:
                closest_pt = min(intersection_points, key=lambda pt: pt.distance(original_mid))
                pt_coords = (round(closest_pt.x, 6), round(closest_pt.y, 6))

                if pt_coords not in corner_points_set:  # Don't duplicate corner points
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

    # Create GeoDataFrame
    gdf_points = gpd.GeoDataFrame(
        {"point_type": point_types},
        geometry=points,
        crs=gdf_boundary.crs,
    )

    # Save to GeoPackage
    gdf_points.to_file(output_path, layer=output_layer_name, driver="GPKG")

    print(f"  Created {len(points)} total guide points")
    print(f"    - {len(corner_points_set)} corner points")
    print(f"    - {intersection_count} perpendicular intersection points")
    print(f"  Saved to layer: {output_layer_name}")

    return output_layer_name


# -------------------------------------------------------------------------
# New function: Create perpendicular lines from guide points, oriented away from site polygon
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
    import math

    import geopandas as gpd
    from shapely.geometry import LineString, Point
    from shapely.ops import unary_union

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

    # Tolerances
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


def generate_local_streets(
    output_path: Path, cfg: StreetConfig, input_layer_name: str
) -> tuple[str, str]:
    """Generate local streets zone from grid blocks.

    Creates a zone for local streets with sidewalks by:
    1. Creating an inner buffer: -(sidewalk_width + road_width/2)
    2. Creating an outer rounded buffer: +sidewalk_width

    Both inner and outer zones are saved as separate layers.

    Args:
        output_path (Path): Path to the GeoPackage containing grid blocks.
        cfg (StreetConfig): Configuration with sidewalk and road width parameters.
        input_layer_name (str): Name of the input grid blocks layer.

    Returns:
        tuple[str, str]: Names of (outer_layer, inner_layer) created.
    """
    output_layer_name = f"{input_layer_name}_local_streets"

    print(f"Creating local streets zone from {input_layer_name}...")

    outer_layer, inner_layer = create_local_streets_zone(
        str(output_path),
        input_layer_name,
        str(output_path),
        output_layer_name,
        cfg.sidewalk_width_m,
        cfg.road_locals_width_m,
    )

    print(f"  Created layers: {outer_layer} (outer), {inner_layer} (inner)")

    return (outer_layer, inner_layer)


def generate_streets(cfg: StreetConfig) -> Path:
    """
    Generate street blocks from roads and parcels (version 3)

    Returns:
        Path to output blocks file
    """
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_gpkg = output_dir / "outputs.gpkg"
    output_path = str(output_gpkg)

    print("Step 1: Determining UTM zone...")
    site_ds = ogr.Open(cfg.parcel_path)
    site_layer = site_ds.GetLayer()
    utm_epsg = get_utm_zone_from_layer(site_layer)
    site_ds = None
    print(f"Using UTM EPSG: {utm_epsg}")

    print("Step 2: Reprojecting layers to UTM...")
    site_layer_name = reproject_layer(cfg.parcel_path, output_path, utm_epsg)
    roads_layer_name = reproject_layer(cfg.roads_path, output_path, utm_epsg)

    print("Step 3: Creating site polygon clipped by roads...")
    clipped_site_layer = clip_site_by_roads(
        output_path,
        site_layer_name,
        roads_layer_name,
        "site_clipped_by_roads",
        cfg.road_arterial_width_m,
        cfg.road_secondary_width_m,
        cfg.road_local_width_m,
    )

    print("Step 4: Extracting arterial roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_art'", output_path, "04_arterial_roads"
    )

    print("Step 5: Extracting secondary roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_sec'", output_path, "05_secondary_roads"
    )

    print("Step 6: Creating arterial road setback zone...")
    buffer_layer(
        output_path,
        "04_arterial_roads",
        cfg.arterial_setback_depth,
        output_path,
        "06_arterial_setback",
        dissolve=True,
    )

    print("Step 7: Creating secondary road setback zone...")
    buffer_layer(
        output_path,
        "05_secondary_roads",
        cfg.secondary_setback_depth,
        output_path,
        "07_secondary_setback",
        dissolve=True,
    )

    print("Step 8: Removing arterial setback from site...")
    erase_layer(
        output_path,
        clipped_site_layer,
        output_path,
        "06_arterial_setback",
        output_path,
        "08_site_minus_arterial_setback",
    )

    print("Step 9: Removing secondary setback from site...")
    erase_layer(
        output_path,
        "08_site_minus_arterial_setback",
        output_path,
        "07_secondary_setback",
        output_path,
        "09_site_minus_all_setbacks",
    )

    print("Step 10: Generating on-grid blocks...")
    generate_on_grid_blocks(output_gpkg, cfg)

    print("Step 13: Extracting site boundary lines")
    extract_site_boundary_lines(
        output_gpkg,
        "09_site_minus_all_setbacks",
        "10a_arterial_setback_clipped",
        "10b_secondary_setback_clipped",
        output_gpkg,
        "13_site_boundary_lines",
    )

    print("Step 14: Generating off-grid blocks...")
    grids_from_site(
        output_gpkg,
        "09_site_minus_all_setbacks",
        "13_site_boundary_lines",
        cfg.off_grid_partitions_preferred_width,
        cfg.off_grid_partitions_preferred_depth,
        grid_layer_name="14_off_grid_cells",
        point_layer_name="14_off_grid_points",
    )

    print("Step 14a: Creating inner buffer zone from site boundary...")
    _buffered_layer = apply_inner_buffer_to_cells(
        output_gpkg,
        "09_site_minus_all_setbacks",
        cfg.part_loc_d,
    )

    print("Step 14b: Extracting grid lines inside buffer zone...")
    lines_layer = extract_grid_lines_in_buffer(
        output_gpkg,
        "14_off_grid_cells",
        "14a_site_boundary_inner_buffer",
    )
    if lines_layer is not None:
        perp_inside_layer = create_perpendicular_lines_inside_buffer_from_points(
            output_gpkg,
            lines_layer,
            "14a_site_boundary_inner_buffer",
            line_length=cfg.part_loc_d * 2,
        )

        if perp_inside_layer is not None:
            print("Step 14d: Fixing grid cells with perpendicular lines...")
            _fixed_cells_layer = fix_grid_cells_with_perpendicular_lines(
                output_gpkg,
                "14_off_grid_cells",
                perp_inside_layer,
                "14a_site_boundary_inner_buffer",
                target_area=cfg.off_grid_partitions_preferred_width
                * cfg.off_grid_partitions_preferred_depth,
            )

    print("Step 15: Generating on-grid cells")

    print("Step 15a: Creating guide points from site boundary...")
    _guide_points_layer = create_guide_points_from_site_boundary(
        output_gpkg,
        "13_site_boundary_lines",
        perp_inside_layer,
        output_layer_name="13_site_boundary_points",
    )

    print("Step 15b: Creating perpendicular lines from guide points...")
    _guide_perp_layer = create_perpendicular_lines_from_guide_points(
        output_gpkg,
        "13_site_boundary_lines",
        "09_site_minus_all_setbacks",
        "13_site_boundary_points",
        line_length=max(cfg.secondary_setback_depth, cfg.arterial_setback_depth) * 1.05,
        output_layer_name="13_site_boundary_perp_from_points",
    )

    print("Merge arterial and secondary setbacks with overlaps resolved...")
    merge_layers_without_overlaps(
        output_path,
        ["10a_arterial_setback_clipped", "10b_secondary_setback_clipped"],
        output_path,
        "10_setback_clipped_merged",
    )

    print("Step 16: Creating on-grid cells from merged setbacks and perpendiculars...")
    _on_grid_cells_layer = create_on_grid_cells_from_perpendiculars(
        output_gpkg,
        "10_setback_clipped_merged",
        "13_site_boundary_perp_from_points",
        output_gpkg,
        "16_on_grid_cells",
    )

    print("Step 17: Merging all grid layers with grid_type information...")
    merge_grid_layers_with_type(
        str(output_gpkg),
        str(output_gpkg),
        "17_all_grids_merged",
        [
            ("10_intersected_setbacks", "on_grid_intersected"),
            ("11_arterial_setback_grid_cleaned", "on_grid_art"),
            ("12_secondary_setback_grid_cleaned", "on_grid_sec"),
        ],
    )

    print("Step 18: Exporting merged grids to GeoJSON...")
    output_geojson = output_dir / "all_grids_merged.geojson"
    export_layer_to_geojson(
        str(output_gpkg),
        "17_all_grids_merged",
        str(output_geojson),
    )

    print(f"\nProcessing complete! Output saved to: {output_gpkg}")
    print("\nFinal layers:")
    print(f"  - {clipped_site_layer}: Site polygon with roads subtracted")
    print("  - 17_all_grids_merged: Merged grid cells with grid_type classification")
    print("\nGeoJSON export:")
    # print(f"  - {output_geojson}: Merged grids with grid_type classification")

    return output_gpkg
