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
