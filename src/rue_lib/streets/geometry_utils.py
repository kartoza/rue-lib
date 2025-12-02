# src/rue_lib/streets/geometry_utils.py
import math

from osgeo import ogr


def is_duplicate(new_geom, existing_geoms):
    """Check if geometry already exists in the output layer."""
    for existing_geom in existing_geoms:
        if new_geom.Equals(existing_geom):
            return True
    return False


def calculate_angle_change(x1, y1, x2, y2, x3, y3):
    """
    Calculate the angle change at point (x2, y2) between three consecutive points.

    Returns angle in degrees (0-180).
    """
    # Vector from p1 to p2
    dx1 = x2 - x1
    dy1 = y2 - y1

    # Vector from p2 to p3
    dx2 = x3 - x2
    dy2 = y3 - y2

    # Calculate angles of each vector
    angle1 = math.atan2(dy1, dx1)
    angle2 = math.atan2(dy2, dx2)

    # Calculate angle difference
    angle_diff = angle2 - angle1

    # Normalize to [-pi, pi]
    while angle_diff > math.pi:
        angle_diff -= 2 * math.pi
    while angle_diff < -math.pi:
        angle_diff += 2 * math.pi

    # Convert to degrees and return absolute value
    return abs(math.degrees(angle_diff))


def break_linestring_by_angle(linestring, angle_threshold=60.0):
    """
    Break a linestring at points where the angle change exceeds the threshold.

    Args:
        linestring: OGR LineString geometry
        angle_threshold: Angle threshold in degrees (default 30.0)

    Returns:
        List of OGR LineString geometries
    """
    point_count = linestring.GetPointCount()

    # If linestring has less than 3 points, return as is
    if point_count < 3:
        return [linestring.Clone()]

    # Find break points
    break_indices = [0]  # Start with first point

    for i in range(1, point_count - 1):
        x1 = linestring.GetX(i - 1)
        y1 = linestring.GetY(i - 1)
        x2 = linestring.GetX(i)
        y2 = linestring.GetY(i)
        x3 = linestring.GetX(i + 1)
        y3 = linestring.GetY(i + 1)

        angle_change = calculate_angle_change(x1, y1, x2, y2, x3, y3)

        if angle_change > angle_threshold:
            break_indices.append(i)

    break_indices.append(point_count - 1)  # End with last point

    # Create linestrings from break indices
    segments = []
    for i in range(len(break_indices) - 1):
        start_idx = break_indices[i]
        end_idx = break_indices[i + 1]

        # Create new linestring
        new_line = ogr.Geometry(ogr.wkbLineString)
        for j in range(start_idx, end_idx + 1):
            new_line.AddPoint(linestring.GetX(j), linestring.GetY(j))

        segments.append(new_line)

    return segments


def edges_from_ring(ring):
    """Extract edges from a polygon ring as individual LineString geometries.

    Args:
        ring: OGR LinearRing geometry

    Returns:
        List of OGR LineString geometries representing edges
    """
    edges = []
    if ring is None:
        return edges
    n = ring.GetPointCount()
    if n < 2:
        return edges
    for i in range(n - 1):
        x1, y1 = ring.GetX(i), ring.GetY(i)
        x2, y2 = ring.GetX(i + 1), ring.GetY(i + 1)
        line = ogr.Geometry(ogr.wkbLineString)
        line.AddPoint(x1, y1)
        line.AddPoint(x2, y2)
        edges.append(line)
    return edges


def extract_edges_from_geom(geom):
    """Extract all edges from a polygon or multipolygon geometry.

    Handles POLYGON, MULTIPOLYGON, GEOMETRYCOLLECTION and falls back to boundary

    Args:
        geom: OGR geometry

    Returns:
        List of OGR LineString geometries representing edges
    """
    edges = []
    if geom is None:
        return edges

    gname = geom.GetGeometryName().upper()

    if gname == "POLYGON":
        # outer + inner rings
        for i in range(geom.GetGeometryCount()):
            ring = geom.GetGeometryRef(i)
            edges.extend(edges_from_ring(ring))

    elif gname == "MULTIPOLYGON":
        for p in range(geom.GetGeometryCount()):
            poly = geom.GetGeometryRef(p)
            if poly is None:
                continue
            for i in range(poly.GetGeometryCount()):
                ring = poly.GetGeometryRef(i)
                edges.extend(edges_from_ring(ring))

    elif gname == "GEOMETRYCOLLECTION":
        # Recursively extract edges from each geometry in the collection
        for i in range(geom.GetGeometryCount()):
            sub_geom = geom.GetGeometryRef(i)
            if sub_geom is not None:
                edges.extend(extract_edges_from_geom(sub_geom))

    else:
        # Try to get boundary for other geometry types
        try:
            boundary = geom.GetBoundary()
            if boundary is not None:
                bname = boundary.GetGeometryName().upper()
                if bname == "LINESTRING":
                    edges.extend(edges_from_ring(boundary))
                elif bname == "MULTILINESTRING":
                    for i in range(boundary.GetGeometryCount()):
                        ring = boundary.GetGeometryRef(i)
                        edges.extend(edges_from_ring(ring))
        except Exception as e:
            # GetBoundary() may fail for some geometry types, just skip
            print(e)
            pass

    return edges


def get_line_endpoints(geom):
    """Get the start and end points of a LineString geometry.

    Args:
        geom: OGR LineString geometry

    Returns:
        Tuple of (start_point, end_point) as (x, y) tuples, or (None, None) if invalid
    """
    n = geom.GetPointCount()
    if n < 2:
        return None, None
    start = (round(geom.GetX(0), 6), round(geom.GetY(0), 6))
    end = (round(geom.GetX(n - 1), 6), round(geom.GetY(n - 1), 6))
    return start, end


def merge_connected_edges(edges_with_info):
    """Merge connected edges into continuous linestrings.

    Args:
        edges_with_info: List of tuples (edge_geom, setback_id, distance, length)

    Returns:
        List of tuples (merged_geom, setback_id, avg_distance, total_length)
    """
    if not edges_with_info:
        return []

    # Group by setback_id
    by_setback = {}
    for edge_geom, setback_id, dist, length in edges_with_info:
        if setback_id not in by_setback:
            by_setback[setback_id] = []
        by_setback[setback_id].append((edge_geom, dist, length))

    merged_results = []

    for setback_id, edges in by_setback.items():
        # Build edge data with endpoints
        edge_data = []
        for edge_geom, dist, length in edges:
            start, end = get_line_endpoints(edge_geom)
            if start is not None and end is not None:
                edge_data.append(
                    {
                        "geom": edge_geom,
                        "start": start,
                        "end": end,
                        "dist": dist,
                        "length": length,
                        "used": False,
                    }
                )

        # Build chains by connecting edges
        chains = []
        for edge in edge_data:
            if edge["used"]:
                continue

            # Start a new chain
            chain = [edge]
            edge["used"] = True

            # Try to extend forward
            current_end = edge["end"]
            extended = True
            while extended:
                extended = False
                for other in edge_data:
                    if other["used"]:
                        continue
                    if other["start"] == current_end:
                        chain.append(other)
                        other["used"] = True
                        current_end = other["end"]
                        extended = True
                        break
                    elif other["end"] == current_end:
                        # Reverse connection
                        chain.append(other)
                        other["used"] = True
                        current_end = other["start"]
                        extended = True
                        break

            # Try to extend backward
            current_start = chain[0]["start"]
            extended = True
            while extended:
                extended = False
                for other in edge_data:
                    if other["used"]:
                        continue
                    if other["end"] == current_start:
                        chain.insert(0, other)
                        other["used"] = True
                        current_start = other["start"]
                        extended = True
                        break
                    elif other["start"] == current_start:
                        # Reverse connection
                        chain.insert(0, other)
                        other["used"] = True
                        current_start = other["end"]
                        extended = True
                        break

            chains.append(chain)

        # Create merged linestrings from chains
        for chain in chains:
            merged_line = ogr.Geometry(ogr.wkbLineString)

            # Add points from first edge
            first_geom = chain[0]["geom"]
            for j in range(first_geom.GetPointCount()):
                merged_line.AddPoint(first_geom.GetX(j), first_geom.GetY(j))

            # Add remaining edges
            for k in range(1, len(chain)):
                edge_geom = chain[k]["geom"]
                prev_end = (
                    merged_line.GetX(merged_line.GetPointCount() - 1),
                    merged_line.GetY(merged_line.GetPointCount() - 1),
                )
                edge_start = (edge_geom.GetX(0), edge_geom.GetY(0))
                edge_end = (
                    edge_geom.GetX(edge_geom.GetPointCount() - 1),
                    edge_geom.GetY(edge_geom.GetPointCount() - 1),
                )

                # Check if we need to reverse
                start_dist = (
                    (prev_end[0] - edge_start[0]) ** 2 + (prev_end[1] - edge_start[1]) ** 2
                ) ** 0.5
                end_dist = (
                    (prev_end[0] - edge_end[0]) ** 2 + (prev_end[1] - edge_end[1]) ** 2
                ) ** 0.5

                if end_dist < start_dist:
                    # Add in reverse
                    for j in range(edge_geom.GetPointCount() - 2, -1, -1):
                        merged_line.AddPoint(edge_geom.GetX(j), edge_geom.GetY(j))
                else:
                    # Add normally, skip first point to avoid duplicate
                    for j in range(1, edge_geom.GetPointCount()):
                        merged_line.AddPoint(edge_geom.GetX(j), edge_geom.GetY(j))

            # Calculate average distance and total length
            total_length = sum(e["length"] for e in chain)
            avg_distance = (
                sum(e["dist"] * e["length"] for e in chain) / total_length
                if total_length > 0
                else 0
            )

            merged_results.append((merged_line, setback_id, avg_distance, total_length))

    return merged_results
