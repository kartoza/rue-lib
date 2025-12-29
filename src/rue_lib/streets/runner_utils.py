import math
from pathlib import Path
from typing import Optional

import geopandas as gpd
import pandas as pd
from shapely import LineString, Point, unary_union
from shapely.errors import TopologicalError
from shapely.geometry import CAP_STYLE, JOIN_STYLE, MultiLineString

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


def polygons_to_lines_graph_based(
    input_path: Path,
    input_layer_names: str | list[str],
    output_path: Path,
    output_layer_name: str = "polygon_boundaries_graph",
    merge_distance: float = 0.1,
) -> str:
    """Extract polygon boundaries using a graph-based approach.

    This function:
    1. Combines all polygons from input layers (preserves individual boundaries)
    2. Extracts vertices and builds a connectivity graph
    3. Merges nearby vertices (based on merge_distance)
    4. Creates lines from the vertex connectivity

    Args:
        input_path: Path to input GeoPackage
        input_layer_names: Layer name(s) to extract boundaries from
        output_path: Path to output GeoPackage
        output_layer_name: Name for output layer
        merge_distance: Distance threshold for merging nearby vertices (default: 0.1)
    """
    input_path = str(input_path)
    output_path = str(output_path)

    layer_names = (
        [input_layer_names] if isinstance(input_layer_names, str) else list(input_layer_names)
    )
    if not layer_names:
        raise RuntimeError("No input layer names provided")

    print("\nExtracting polygon boundaries using graph-based approach...")
    print(f"  Input layers: {', '.join(layer_names)}")

    # Read all layers and combine
    all_gdfs = []
    for layer_name in layer_names:
        try:
            gdf = gpd.read_file(input_path, layer=layer_name)
            all_gdfs.append(gdf)
            print(f"  Loaded {len(gdf)} features from '{layer_name}'")
        except Exception as e:
            print(f"  Warning: cannot read layer '{layer_name}': {e}")

    if not all_gdfs:
        raise RuntimeError("No layers could be read")

    # Determine common CRS (use the first layer's CRS)
    target_crs = all_gdfs[0].crs

    # Transform all layers to the common CRS if needed
    for i, gdf in enumerate(all_gdfs):
        if gdf.crs is not None and gdf.crs != target_crs:
            print(f"  Transforming layer {layer_names[i]} from {gdf.crs} to {target_crs}")
            all_gdfs[i] = gdf.to_crs(target_crs)

    gdf_combined = gpd.GeoDataFrame(pd.concat(all_gdfs, ignore_index=True), crs=target_crs)
    crs = gdf_combined.crs

    # Step 1: Collect all polygons (without merging)
    print("  Collecting polygons from all features...")
    all_polygons = []
    for geom in gdf_combined.geometry:
        if geom is None or geom.is_empty:
            continue

        if geom.geom_type == "Polygon":
            all_polygons.append(geom)
        elif geom.geom_type == "MultiPolygon":
            all_polygons.extend(list(geom.geoms))

    print(f"  Collected {len(all_polygons)} polygon(s)")

    # Step 2: Extract vertices and build connectivity graph
    print("  Building vertex connectivity graph...")

    # vertex_coords will store unique vertex coordinates
    # vertex_connections will store which vertices connect to each vertex
    vertex_coords = []
    vertex_id_map = {}
    vertex_connections = {}

    def get_or_create_vertex(coord):
        """Get existing vertex ID or create new one"""
        key = (round(coord[0], 10), round(coord[1], 10))
        if key not in vertex_id_map:
            vertex_id = len(vertex_coords)
            vertex_coords.append((coord[0], coord[1]))
            vertex_id_map[key] = vertex_id
            vertex_connections[vertex_id] = set()
        return vertex_id_map[key]

    def add_ring_to_graph(ring_coords):
        """Add a polygon ring to the connectivity graph"""
        if len(ring_coords) < 2:
            return
        vertex_ids = [get_or_create_vertex(coord) for coord in ring_coords]

        for i in range(len(vertex_ids)):
            current_id = vertex_ids[i]
            next_id = vertex_ids[(i + 1) % len(vertex_ids)]
            vertex_connections[current_id].add(next_id)
            vertex_connections[next_id].add(current_id)

    for poly in all_polygons:
        if poly is None or poly.is_empty:
            continue
        exterior_coords = list(poly.exterior.coords)
        add_ring_to_graph(exterior_coords)
        for interior in poly.interiors:
            interior_coords = list(interior.coords)
            add_ring_to_graph(interior_coords)

    print(f"  Extracted {len(vertex_coords)} unique vertices")

    # Step 3: Merge nearby vertices
    if merge_distance > 0:
        print(f"  Merging vertices within {merge_distance} units...")

        gdf_vertices = gpd.GeoDataFrame(geometry=[Point(coord) for coord in vertex_coords], crs=crs)
        spatial_index = gdf_vertices.sindex

        vertex_merge_map = {}
        merged_count = 0

        for vid in range(len(vertex_coords)):
            if vid in vertex_merge_map:
                continue

            vertex_point = gdf_vertices.iloc[vid].geometry
            buffer_geom = vertex_point.buffer(merge_distance)
            nearby_indices = list(spatial_index.intersection(buffer_geom.bounds))

            merge_candidates = []
            for candidate_vid in nearby_indices:
                if candidate_vid <= vid:
                    continue
                if candidate_vid in vertex_merge_map:
                    continue

                candidate_point = gdf_vertices.iloc[candidate_vid].geometry
                if vertex_point.distance(candidate_point) < merge_distance:
                    merge_candidates.append(candidate_vid)

            for candidate_vid in merge_candidates:
                vertex_merge_map[candidate_vid] = vid
                merged_count += 1

                for connected_vid in vertex_connections[candidate_vid]:
                    if connected_vid != candidate_vid and connected_vid != vid:
                        actual_connected = vertex_merge_map.get(connected_vid, connected_vid)
                        vertex_connections[vid].add(actual_connected)

                vertex_connections[candidate_vid].clear()

        print(f"  Merged {merged_count} vertices")

        for vid in vertex_connections:
            new_connections = set()
            for connected_vid in vertex_connections[vid]:
                actual_vid = vertex_merge_map.get(connected_vid, connected_vid)
                if actual_vid != vid:  # Don't create self-loops
                    new_connections.add(actual_vid)
            vertex_connections[vid] = new_connections

    debug_vertices = []
    debug_vertex_info = []
    for vid in range(len(vertex_coords)):
        if vid in vertex_merge_map:
            continue

        connections = vertex_connections[vid]
        if not connections:
            continue

        coord = vertex_coords[vid]
        debug_vertices.append(Point(coord[0], coord[1]))
        debug_vertex_info.append(
            {
                "vertex_id": vid,
                "num_connections": len(connections),
                "connected_to": ",".join(map(str, sorted(connections))),
                "x": float(coord[0]),
                "y": float(coord[1]),
            }
        )

    if debug_vertices:
        debug_layer_name = f"{output_layer_name}_vertices"
        gdf_debug_verts = gpd.GeoDataFrame(debug_vertex_info, geometry=debug_vertices, crs=crs)
        gdf_debug_verts.to_file(output_path, layer=debug_layer_name, driver="GPKG")
        print(f"  Saved {len(debug_vertices)} vertices to: {debug_layer_name}")

    print("  Checking for redundant collinear connections...")

    def are_collinear(p1, p2, p3, tolerance=0.01):
        """Check if three points are collinear using cross product

        Args:
            p1, p2, p3: Points as (x, y) tuples
            tolerance: Maximum cross product value to consider collinear
        """
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3

        cross = (y2 - y1) * (x3 - x2) - (y3 - y2) * (x2 - x1)
        return abs(cross) < tolerance

    def point_between(p1, p2, p3, tolerance=0.1):
        """Check if p2 is between p1 and p3

        Args:
            p1, p2, p3: Points as (x, y) tuples
            tolerance: Bounding box tolerance in meters
        """
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3

        min_x, max_x = (min(x1, x3), max(x1, x3))
        min_y, max_y = (min(y1, y3), max(y1, y3))

        return (
            min_x - tolerance <= x2 <= max_x + tolerance
            and min_y - tolerance <= y2 <= max_y + tolerance
        )

    redundant_edges = set()

    for vid in range(len(vertex_coords)):
        if vid in vertex_merge_map:
            continue

        connections = vertex_connections[vid]
        coord_a = vertex_coords[vid]
        conn_list = list(connections)
        for i, vid_b in enumerate(conn_list):
            if vid_b in vertex_merge_map:
                continue

            coord_b = vertex_coords[vid_b]

            for vid_c in conn_list[i + 1 :]:
                if vid_c in vertex_merge_map:
                    continue

                coord_c = vertex_coords[vid_c]
                if are_collinear(coord_a, coord_b, coord_c):
                    if vid_c in vertex_connections.get(vid_b, set()):
                        if point_between(coord_a, coord_b, coord_c):
                            edge_key = tuple(sorted([vid, vid_c]))
                            redundant_edges.add(edge_key)
                        elif point_between(coord_a, coord_c, coord_b):
                            edge_key = tuple(sorted([vid, vid_b]))
                            redundant_edges.add(edge_key)

    if redundant_edges:
        print(f"  Found {len(redundant_edges)} redundant collinear connections to remove")

        for vid in vertex_connections:
            new_connections = set()
            for connected_vid in vertex_connections[vid]:
                edge_key = tuple(sorted([vid, connected_vid]))
                if edge_key not in redundant_edges:
                    new_connections.add(connected_vid)
            vertex_connections[vid] = new_connections

    # Step 5: Create lines from connectivity graph
    print("  Creating lines from vertex connectivity...")

    lines = []
    line_info = []
    created_edges = set()

    for vid in range(len(vertex_coords)):
        if vid in vertex_merge_map:
            continue

        coord = vertex_coords[vid]
        connections = vertex_connections[vid]

        for connected_vid in connections:
            edge_key = tuple(sorted([vid, connected_vid]))

            if edge_key in created_edges:
                continue

            created_edges.add(edge_key)

            connected_coord = vertex_coords[connected_vid]
            line = LineString([coord, connected_coord])
            lines.append(line)
            line_info.append(
                {
                    "line_id": len(lines),
                    "vertex_from": vid,
                    "vertex_to": connected_vid,
                    "length": float(line.length),
                }
            )

    print(f"  Created {len(lines)} lines from connectivity graph")

    if not lines:
        print("  Warning: No lines created!")
        return output_layer_name

    # Save output
    gdf_lines = gpd.GeoDataFrame(line_info, geometry=lines, crs=crs)
    gdf_lines.to_file(output_path, layer=output_layer_name, driver="GPKG")
    print(f"  Saved to layer: {output_layer_name}")

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
        line_parts_extended.append(
            line_part.buffer(
                buffer_distance,
                join_style=JOIN_STYLE.mitre,
                cap_style=CAP_STYLE.round,
            )
        )

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
    site_boundary_lines_layer: str,
    line_length: float,
    output_layer_name: Optional[str] = None,
    tangent_step_fraction: float = 0.001,
    min_endpoint_distance: float = 10.0,
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
        min_endpoint_distance: Minimum distance from line endpoint to next/prev buffer points.
                               Lines with endpoints closer than this are skipped (default: 80.0).

    Returns:
        Name of the created perpendicular lines layer, or None on failure.
        Debug lines are saved to a layer named "{output_layer_name}_debug".
    """
    import geopandas as gpd

    print(f"\nCreating perpendicular lines inside buffer from points in layer: {points_layer}")

    # Load data
    gdf_pts = gpd.read_file(output_path, layer=points_layer)
    gdf_buffer = gpd.read_file(output_path, layer=buffer_layer)
    gdf_boundary_lines = gpd.read_file(output_path, layer=site_boundary_lines_layer)

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
    filtered_count = 0

    # Debug lines from endpoints to next/prev buffer points
    debug_lines = []
    debug_records = []

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

        line_coords = list(line_geom.coords)
        line_end = Point(line_coords[-1])
        line_start = Point(line_coords[0])

        # Initialize vertices from boundary line (will be set if endpoint intersects boundary)
        prev_vertex = None
        next_vertex = None

        # Check if line start point intersects with any site boundary line
        # Use small buffer for numerical tolerance
        if gdf_boundary_lines.intersects(line_start.buffer(0.1)).any():
            # Line starts on a boundary line
            pass

        if gdf_boundary_lines.intersects(line_end.buffer(0.1)).any():
            # Line ends on a boundary line - get the intersecting line
            intersecting_lines = gdf_boundary_lines[
                gdf_boundary_lines.intersects(line_end.buffer(0.1))
            ]

            if not intersecting_lines.empty:
                # Get the first intersecting line
                boundary_line = intersecting_lines.iloc[0].geometry

                # Project the line_end point onto the boundary line to find its position
                distance_along_line = boundary_line.project(line_end)

                # Get the vertices (coordinates) of the boundary line
                boundary_coords = list(boundary_line.coords)

                # Find which segment the point falls on
                cumulative_dist = 0
                prev_vertex = None
                next_vertex = None

                for i in range(len(boundary_coords) - 1):
                    p1 = Point(boundary_coords[i])
                    p2 = Point(boundary_coords[i + 1])
                    segment_length = p1.distance(p2)

                    if cumulative_dist <= distance_along_line <= cumulative_dist + segment_length:
                        # The point falls on this segment
                        prev_vertex = Point(boundary_coords[i])
                        next_vertex = Point(boundary_coords[i + 1])
                        break

                    cumulative_dist += segment_length

                # If point is at the start or end, handle edge cases
                if prev_vertex is None and next_vertex is None:
                    if distance_along_line <= 0:
                        # Point is at the start
                        prev_vertex = None
                        next_vertex = (
                            Point(boundary_coords[1]) if len(boundary_coords) > 1 else None
                        )
                    else:
                        # Point is at the end
                        prev_vertex = (
                            Point(boundary_coords[-2]) if len(boundary_coords) > 1 else None
                        )
                        next_vertex = None

        # Calculate distances to next and previous points
        # If line endpoint intersects a boundary line, use vertices from that line
        # Otherwise, use buffer boundary points
        if prev_vertex is not None and next_vertex is not None:
            # Use vertices from the boundary line
            next_pt = next_vertex
            prev_pt = prev_vertex
            dist_to_next = line_end.distance(next_pt)
            dist_to_prev = line_end.distance(prev_pt)
        else:
            # Fall back to buffer boundary logic
            search_distance = line_length * 2

            s_next = (s + search_distance) % boundary_len
            next_pt = buffer_boundary.interpolate(s_next)

            s_prev = (s - search_distance) % boundary_len
            prev_pt = buffer_boundary.interpolate(s_prev)

            dist_to_next = line_end.distance(next_pt)
            dist_to_prev = line_end.distance(prev_pt)

        # Create debug lines from endpoint to next/prev points
        debug_line_to_next = LineString([line_end, next_pt])
        debug_line_to_prev = LineString([line_end, prev_pt])

        debug_lines.append(debug_line_to_next)
        debug_records.append(
            {
                "source_idx": _idx,
                "line_type": "to_next",
                "distance": float(dist_to_next),
                "filtered": dist_to_next < min_endpoint_distance
                or dist_to_prev < min_endpoint_distance,
            }
        )

        debug_lines.append(debug_line_to_prev)
        debug_records.append(
            {
                "source_idx": _idx,
                "line_type": "to_prev",
                "distance": float(dist_to_prev),
                "filtered": dist_to_next < min_endpoint_distance
                or dist_to_prev < min_endpoint_distance,
            }
        )

        # Skip lines where endpoint is too close to next or prev buffer points
        if dist_to_next < min_endpoint_distance or dist_to_prev < min_endpoint_distance:
            filtered_count += 1
            continue

        lines.append(line_geom)
        rec = row.drop(labels="geometry").to_dict()
        rec["dist_end_to_next"] = float(dist_to_next)
        rec["dist_end_to_prev"] = float(dist_to_prev)
        rec["line_length"] = float(line_geom.length)
        records.append(rec)

    if not lines:
        print("  Warning: no perpendicular lines created")
        if filtered_count > 0:
            print(
                f"  {filtered_count} lines were filtered out (endpoint too close to buffer points)"
            )
        return None

    gdf_out = gpd.GeoDataFrame(records, geometry=lines, crs=gdf_pts.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")

    print(f"  Created {len(lines)} perpendicular lines inside buffer")
    if filtered_count > 0:
        print(f"  Filtered out {filtered_count} lines (endpoint too close to buffer points)")
    print(f"  Saved to layer: {output_layer_name}")

    # Save debug lines showing endpoint distances to next/prev buffer points
    if debug_lines:
        debug_layer_name = f"{output_layer_name}_debug"
        gdf_debug = gpd.GeoDataFrame(debug_records, geometry=debug_lines, crs=gdf_pts.crs)
        gdf_debug.to_file(output_path, layer=debug_layer_name, driver="GPKG")
        print(f"  Saved {len(debug_lines)} debug lines to layer: {debug_layer_name}")

    return output_layer_name


def create_guide_points_from_site_boundary(
    output_path: Path,
    site_boundary_lines_layer: str,
    perp_lines_layer: Optional[str],
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
        perp_lines_layer: Name of layer containing perpendicular lines (optional)
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

    # Load perpendicular lines if layer name is provided
    gdf_perp = None
    if perp_lines_layer is not None:
        try:
            gdf_perp = gpd.read_file(output_path, layer=perp_lines_layer)
        except Exception as e:
            print(f"  Warning: Could not load perpendicular lines layer '{perp_lines_layer}': {e}")
            gdf_perp = None

    print(f"  Loaded {len(gdf_boundary)} boundary line(s)")
    if gdf_perp is not None:
        print(f"  Loaded {len(gdf_perp)} perpendicular line(s)")
    else:
        print("  No perpendicular lines layer available")

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

    intersection_count = 0
    if gdf_perp is not None and not gdf_perp.empty:
        print("  Finding perpendicular line intersections...")

        boundary_union = unary_union(boundary_lines)

        extension_factor = 1000.0

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

                extended_start = (
                    start[0] - dx * extension_factor,
                    start[1] - dy * extension_factor,
                )
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
    simplify: bool = False,
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
            erase_union = erase_union.buffer(
                buffer_distance, cap_style=CAP_STYLE.square, join_style=JOIN_STYLE.mitre
            )

            if simplify:
                simplify_tolerance = buffer_distance * 0.1
                erase_union = erase_union.simplify(simplify_tolerance, preserve_topology=True)

            if not erase_union.is_valid:
                erase_union = erase_union.buffer(0)
        else:
            # Small buffer to clean geometry
            erase_union = erase_union.buffer(
                0.001, cap_style=CAP_STYLE.square, join_style=JOIN_STYLE.mitre
            )
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
