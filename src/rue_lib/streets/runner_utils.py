import math
from pathlib import Path
from typing import Optional

import geopandas as gpd
import pandas as pd
from shapely import LineString, Point, unary_union
from shapely.errors import TopologicalError
from shapely.geometry import MultiLineString

from rue_lib.core.geometry import buffer_layer

from .operations import (
    clip_layer,
)


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


def polygons_to_lines_layer(
    input_path: Path,
    input_layer_names: str | list[str],
    output_path: Path,
    output_layer_name: str = "polygon_boundaries",
    dedupe_distance: float = 0.01,
) -> str:
    """Extract polygon boundaries from one or more layers and write them as LineStrings.

    This preserves individual polygon edges; shared edges are kept from both polygons
    instead of dissolving the polygons into a single union.

    Lines that are within dedupe_distance of each other are deduplicated.

    Args:
        input_path: Path to input GeoPackage
        input_layer_names: Layer name(s) to extract boundaries from
        output_path: Path to output GeoPackage
        output_layer_name: Name for output layer
        dedupe_distance: Distance threshold for removing duplicate lines (default: 0.01)
    """
    input_path = str(input_path)
    output_path = str(output_path)

    layer_names = (
        [input_layer_names] if isinstance(input_layer_names, str) else list(input_layer_names)
    )
    if not layer_names:
        raise RuntimeError("No input layer names provided for polygons_to_lines_layer")

    # Read all layers and combine
    all_gdfs = []
    for layer_name in layer_names:
        try:
            gdf = gpd.read_file(input_path, layer=layer_name)
            all_gdfs.append(gdf)
        except Exception as e:
            print(f"  Warning: cannot read layer '{layer_name}': {e}")
            return

    if not all_gdfs:
        raise RuntimeError("No layers could be read")

    gdf_combined = gpd.GeoDataFrame(pd.concat(all_gdfs, ignore_index=True))
    crs = gdf_combined.crs

    # Extract all boundary lines from polygons
    line_records = []
    feature_idx = 1

    for row in gdf_combined.itertuples():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        if not geom.is_valid:
            geom = geom.buffer(0)

        # Handle different geometry types
        if geom.geom_type == "Polygon":
            polygons = [geom]
        elif geom.geom_type == "MultiPolygon":
            polygons = list(geom.geoms)
        else:
            continue

        # Extract rings from each polygon
        for poly in polygons:
            rings = [poly.exterior] + list(poly.interiors)
            ring_id = 1
            for ring in rings:
                if ring.is_empty:
                    continue
                line_records.append(
                    {
                        "src_id": feature_idx,
                        "ring_id": ring_id,
                        "geometry": LineString(ring.coords),
                    }
                )
                ring_id += 1
            feature_idx += 1

    if not line_records:
        raise RuntimeError(f"No geometries found to convert in layers: {', '.join(layer_names)}")

    # Create GeoDataFrame of all lines
    gdf_lines = gpd.GeoDataFrame(line_records, geometry="geometry", crs=crs)

    # Deduplicate by removing close points, then reconstructing lines
    if dedupe_distance > 0:
        print(
            f"  Deduplicating points from {len(gdf_lines)} lines "
            f"(distance threshold: {dedupe_distance})..."
        )

        # Extract all unique points from all lines (normalize to 2D)
        all_points = []
        point_to_lines = {}  # Map points to the lines they belong to

        for idx, row in gdf_lines.iterrows():
            line_coords = list(row.geometry.coords)
            for coord in line_coords:
                # Normalize to 2D (take only x, y)
                coord_2d = (coord[0], coord[1])
                point_key = (round(coord_2d[0], 10), round(coord_2d[1], 10))  # Round for grouping
                if point_key not in point_to_lines:
                    all_points.append(coord_2d)
                    point_to_lines[point_key] = []
                point_to_lines[point_key].append(idx)

        print(f"  Extracted {len(all_points)} unique points")

        # Create GeoDataFrame of points for spatial indexing
        gdf_points = gpd.GeoDataFrame(geometry=[Point(coord) for coord in all_points], crs=crs)

        spatial_index = gdf_points.sindex
        keep_point_indices = set(range(len(all_points)))
        point_mapping = {}

        for idx in range(len(all_points)):
            if idx not in keep_point_indices:
                continue

            point = gdf_points.iloc[idx].geometry

            buffer_geom = point.buffer(dedupe_distance)
            possible_matches_idx = list(spatial_index.intersection(buffer_geom.bounds))

            for candidate_idx in possible_matches_idx:
                if candidate_idx <= idx or candidate_idx not in keep_point_indices:
                    continue

                candidate_point = gdf_points.iloc[candidate_idx].geometry

                if point.distance(candidate_point) < dedupe_distance:
                    point_mapping[candidate_idx] = idx
                    keep_point_indices.discard(candidate_idx)

        kept_points = [all_points[i] for i in sorted(keep_point_indices)]
        old_to_new_idx = {
            old_idx: new_idx for new_idx, old_idx in enumerate(sorted(keep_point_indices))
        }

        for removed_idx, kept_idx in point_mapping.items():
            old_to_new_idx[removed_idx] = old_to_new_idx[kept_idx]

        print(f"  After deduplication: {len(kept_points)} points")
        new_line_records = []
        for _idx, row in gdf_lines.iterrows():
            line_coords = list(row.geometry.coords)
            new_coords = []
            for coord in line_coords:
                coord_2d = (coord[0], coord[1])
                for point_idx, pt in enumerate(all_points):
                    if abs(pt[0] - coord_2d[0]) < 1e-9 and abs(pt[1] - coord_2d[1]) < 1e-9:
                        new_idx = old_to_new_idx.get(point_idx, point_idx)
                        if new_idx < len(kept_points):
                            new_coord = kept_points[new_idx]
                            if not new_coords or new_coords[-1] != new_coord:
                                new_coords.append(new_coord)
                        break

            if len(new_coords) >= 2:
                new_line_records.append(
                    {
                        "src_id": row["src_id"],
                        "ring_id": row["ring_id"],
                        "geometry": LineString(new_coords),
                    }
                )

        gdf_lines = gpd.GeoDataFrame(new_line_records, geometry="geometry", crs=crs)
        print(f"  After reconstruction: {len(gdf_lines)} lines")

    # Write output
    gdf_lines.to_file(output_path, layer=output_layer_name, driver="GPKG")

    return output_layer_name


def create_dead_end_boundary_lines(
    output_path: Path,
    site_minus_setbacks_layer: str,
    site_boundary_lines_layer: str,
    output_layer_name: str = "09_dead_end_lines",
    diff_buffer: float = 0.05,
    buffer_distance: float = 10.0,
) -> str:
    """Create dead-end boundary lines"""
    gdf_site = gpd.read_file(output_path, layer=site_minus_setbacks_layer)
    gdf_boundary = gpd.read_file(output_path, layer=site_boundary_lines_layer)

    dead_end_lines_buffered_layer = f"{output_layer_name}_buffered"

    # Convert polygons to exterior lines
    line_geoms = []
    for geom in gdf_site.geometry:
        if geom is None or geom.is_empty:
            continue
        if geom.geom_type == "Polygon":
            line_geoms.append(LineString(geom.exterior))
        elif geom.geom_type == "MultiPolygon":
            for poly in geom.geoms:
                line_geoms.append(LineString(poly.exterior))

    if not line_geoms:
        raise RuntimeError("No polygon boundaries found to convert to lines")

    site_lines_union = unary_union(line_geoms)
    boundary_union = unary_union(gdf_boundary.geometry) if not gdf_boundary.empty else None

    if boundary_union and not boundary_union.is_empty:
        site_lines_union = site_lines_union.difference(boundary_union.buffer(diff_buffer))

    # Flatten to LineString parts
    line_parts = []
    if site_lines_union.geom_type == "LineString":
        line_parts = [site_lines_union]
    elif site_lines_union.geom_type == "MultiLineString":
        line_parts = list(site_lines_union.geoms)
    elif site_lines_union.geom_type == "GeometryCollection":
        line_parts = [
            g for g in site_lines_union.geoms if g.geom_type in ("LineString", "MultiLineString")
        ]
        expanded = []
        for g in line_parts:
            if isinstance(g, MultiLineString):
                expanded.extend(list(g.geoms))
            else:
                expanded.append(g)
        line_parts = expanded

    line_parts_extended = []
    for line_part in line_parts:
        line_parts_extended.append(line_part.buffer(buffer_distance))

    if not line_parts:
        raise RuntimeError("No dead-end boundary lines found after subtraction")

    gdf_out = gpd.GeoDataFrame(geometry=line_parts, crs=gdf_site.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")

    gdf_out = gpd.GeoDataFrame(geometry=line_parts_extended, crs=gdf_site.crs)
    gdf_out.to_file(output_path, layer=dead_end_lines_buffered_layer, driver="GPKG")

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
    lines_without_points_layer: str = "13_lines_without_points",
    min_line_length_threshold: float = 100.0,
) -> str:
    """Create guide points from site boundary lines and perpendicular lines.

    This function creates multiple types of points (excluding dead-end corners):
    1. Non-dead-end corner points where 2+ lines meet
    2. Points where perpendicular lines intersect/touch site boundary lines
    3. Points at regular intervals along long segments (segments >= min_line_length_threshold)

    Args:
        output_path: Path to the GeoPackage
        site_boundary_lines_layer: Name of layer containing site boundary lines
        perp_lines_layer: Name of layer containing perpendicular lines
        output_layer_name: Name for the output points layer
        lines_without_points_layer: Name for layer containing long segments between points
        min_line_length_threshold: Minimum length for segments and interval for adding points
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
            for coord in (coords[0], coords[-1]):
                key = (round(coord[0], 6), round(coord[1], 6))
                endpoint_counts[key] = endpoint_counts.get(key, 0) + 1

    corner_points_set = {pt for pt, cnt in endpoint_counts.items() if cnt > 1}

    for pt_coords in corner_points_set:
        points.append(Point(pt_coords[0], pt_coords[1]))
        point_types.append("corner")

    print("  Finding perpendicular line intersections...")

    boundary_union = unary_union(boundary_lines)

    extension_factor = 1000.0

    intersection_count = 0
    for _, row in gdf_perp.iterrows():
        perp_geom = row.geometry
        if perp_geom is None or perp_geom.is_empty:
            continue

        try:
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

    print(f"\n  Breaking boundary lines at {len(points)} points...")
    breaking_points_geom = points

    broken_segments = []
    segment_lengths = []

    for line in boundary_lines:
        points_on_line = []
        snap_distance = 1.0

        for pt in breaking_points_geom:
            if line.distance(pt) < snap_distance:
                distance_along = line.project(pt)
                points_on_line.append(distance_along)

        points_on_line = sorted(set(points_on_line))

        if not points_on_line:
            broken_segments.append(line)
            segment_lengths.append(line.length)
        else:
            coords = list(line.coords)
            line_positions = [0.0] + points_on_line + [line.length]

            for i in range(len(line_positions) - 1):
                start_dist = line_positions[i]
                end_dist = line_positions[i + 1]

                if end_dist - start_dist < 0.1:
                    continue

                start_pt = line.interpolate(start_dist)
                end_pt = line.interpolate(end_dist)

                segment = LineString([start_pt, end_pt])
                broken_segments.append(segment)
                segment_lengths.append(segment.length)

    print(f"  Created {len(broken_segments)} segments from boundary lines")

    long_segments = []
    long_lengths = []

    for segment, length in zip(broken_segments, segment_lengths):
        if length >= min_line_length_threshold * 2:
            long_segments.append(segment)
            long_lengths.append(length)

    if long_segments:
        sorted_pairs = sorted(zip(long_segments, long_lengths), key=lambda x: x[1], reverse=True)
        long_segments = [seg for seg, _ in sorted_pairs]
        long_lengths = [length for _, length in sorted_pairs]

        print(f"  Found {len(long_segments)} segments longer than {min_line_length_threshold}m")

        gdf_long_segments = gpd.GeoDataFrame(
            {"length_m": long_lengths},
            geometry=long_segments,
            crs=gdf_boundary.crs,
        )

        gdf_long_segments.to_file(output_path, layer=lines_without_points_layer, driver="GPKG")
        print(f"  Saved longest segments to layer: {lines_without_points_layer}")

        for i, length in enumerate(long_lengths[:10]):
            print(f"    Segment {i + 1}: length = {length:.2f}m")
        if len(long_lengths) > 10:
            print(f"    ... and {len(long_lengths) - 10} more segments")

        print(f"\n  Adding points along long segments at {min_line_length_threshold}m intervals...")
        additional_points_count = 0
        for segment in long_segments:
            if segment.geom_type == "LineString":
                segment_length = segment.length
                num_intervals = int(segment_length / min_line_length_threshold)

                for i in range(1, num_intervals + 1):
                    distance = i * min_line_length_threshold
                    if distance < segment_length:
                        point = segment.interpolate(distance)
                        key = (round(point.x, 6), round(point.y, 6))

                        if key not in corner_points_set and not any(
                            abs(pt.x - point.x) < 0.001 and abs(pt.y - point.y) < 0.001
                            for pt in points
                        ):
                            points.append(Point(point.x, point.y))
                            point_types.append("long_segment_interval")
                            additional_points_count += 1

        print(f"  Added {additional_points_count} additional points along long segments")
        gdf_points = gpd.GeoDataFrame(
            {"point_type": point_types},
            geometry=points,
            crs=gdf_boundary.crs,
        )
        gdf_points.to_file(output_path, layer=output_layer_name, driver="GPKG")
        print(f"  Updated guide points layer with {len(points)} total points")

    else:
        print(f"  No segments longer than {min_line_length_threshold}m found")

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


try:
    # Shapely 2.x
    from shapely import make_valid as _make_valid
except Exception:
    _make_valid = None


def _fix_geom(g):
    """Fix invalid geometries, returning None if empty/invalid."""
    if g is None or g.is_empty:
        return None
    try:
        if _make_valid is not None:
            g2 = _make_valid(g)
        else:
            # classic "buffer(0)" fix
            g2 = g.buffer(0)
        if g2 is None or g2.is_empty:
            return None
        return g2
    except Exception:
        # last resort: try buffer(0)
        try:
            g2 = g.buffer(0)
            if g2 is None or g2.is_empty:
                return None
            return g2
        except Exception:
            return None


def subtract_layer(
    input_gpkg: str,
    base_layer_name: str,
    erase_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
    buffer_distance: float = 0.0,
) -> str:
    """
    Subtract one layer from another using geometric difference (GeoPandas/Shapely).

    Subtracts the geometries from erase_layer (optionally buffered) from base_layer.
    Both layers must be in the same input GeoPackage.

    Args:
        input_gpkg: Path to input GeoPackage containing both layers
        base_layer_name: Name of the layer to subtract from
        erase_layer_name: Name of the layer to subtract
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for the output layer
        buffer_distance: Optional buffer distance to apply to erase layer before subtraction

    Returns:
        Name of the output layer
    """
    gdf_base = gpd.read_file(input_gpkg, layer=base_layer_name)
    gdf_erase = gpd.read_file(input_gpkg, layer=erase_layer_name)

    if gdf_base.empty:
        raise ValueError(f"Layer '{base_layer_name}' is empty")
    if gdf_erase.empty:
        print("  Erase layer is empty. Writing base layer unchanged...")

    if gdf_erase.crs is not None and gdf_base.crs is not None and gdf_erase.crs != gdf_base.crs:
        gdf_erase = gdf_erase.to_crs(gdf_base.crs)

    gdf_base = gdf_base.copy()
    gdf_erase = gdf_erase.copy()

    gdf_base["geometry"] = gdf_base.geometry.apply(_fix_geom)
    gdf_erase["geometry"] = gdf_erase.geometry.apply(_fix_geom)

    gdf_base = gdf_base[gdf_base.geometry.notnull()].reset_index(drop=True)
    gdf_erase = gdf_erase[gdf_erase.geometry.notnull()].reset_index(drop=True)

    print(f"  Processing {len(gdf_base)} features from base layer...")

    erase_union = None
    if not gdf_erase.empty:
        erase_union = unary_union(list(gdf_erase.geometry))

    if erase_union is None:
        erased_geom = gdf_base.geometry
    else:
        # Buffer the erase geometry if buffer_distance is specified
        if buffer_distance > 0:
            erase_union = erase_union.buffer(buffer_distance, cap_style=3, join_style=2)
        else:
            erase_union = erase_union.buffer(0.001, cap_style=3, join_style=2)
        erased_geom = gdf_base.geometry.difference(erase_union)

    # Create a new GeoDataFrame with the erased geometry
    out = gdf_base.copy()
    out["geometry"] = erased_geom

    out["geometry"] = out.geometry.apply(_fix_geom)
    out = out[out.geometry.notnull()].reset_index(drop=True)

    # Ensure it's a GeoDataFrame before writing
    out = gpd.GeoDataFrame(out, geometry="geometry", crs=gdf_base.crs)
    out.to_file(output_gpkg, layer=output_layer_name, driver="GPKG")

    print(f"  Created layer: {output_layer_name}")
    return output_layer_name
