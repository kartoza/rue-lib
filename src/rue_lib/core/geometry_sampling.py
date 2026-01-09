from __future__ import annotations

from typing import Union

from shapely import wkb, wkt
from shapely.geometry import LineString, MultiLineString, Point
from shapely.ops import linemerge

GeomInput = Union[str, bytes, bytearray, memoryview]


def extend_line(line: LineString, extension_m: float) -> LineString:
    """Extend a LineString by a given distance at both ends.

    Args:
        line: The LineString to extend (supports 2D and 3D coordinates)
        extension_m: Distance in meters to extend at each end

    Returns:
        A new LineString extended by extension_m at both ends
    """
    if not isinstance(line, LineString) or len(line.coords) < 2:
        return line

    coords = list(line.coords)

    has_z = line.has_z

    # Extend start
    if has_z:
        x1, y1, z1 = coords[0][0], coords[0][1], coords[0][2]
        x2, y2, z2 = coords[1][0], coords[1][1], coords[1][2]
    else:
        x1, y1 = coords[0][0], coords[0][1]
        x2, y2 = coords[1][0], coords[1][1]

    dx, dy = x2 - x1, y2 - y1
    length = (dx**2 + dy**2) ** 0.5
    if length > 0:
        dx, dy = dx / length, dy / length
        if has_z:
            new_start = (x1 - dx * extension_m, y1 - dy * extension_m, z1)
        else:
            new_start = (x1 - dx * extension_m, y1 - dy * extension_m)
    else:
        if has_z:
            new_start = (x1, y1, z1)
        else:
            new_start = (x1, y1)

    # Extend end
    if has_z:
        x1, y1, z1 = coords[-2][0], coords[-2][1], coords[-2][2]
        x2, y2, z2 = coords[-1][0], coords[-1][1], coords[-1][2]
    else:
        x1, y1 = coords[-2][0], coords[-2][1]
        x2, y2 = coords[-1][0], coords[-1][1]

    dx, dy = x2 - x1, y2 - y1
    length = (dx**2 + dy**2) ** 0.5
    if length > 0:
        dx, dy = dx / length, dy / length
        if has_z:
            new_end = (x2 + dx * extension_m, y2 + dy * extension_m, z2)
        else:
            new_end = (x2 + dx * extension_m, y2 + dy * extension_m)
    else:
        if has_z:
            new_end = (x2, y2, z2)
        else:
            new_end = (x2, y2)

    return LineString([new_start] + coords + [new_end])


def _load_line_geom(line_wkt_or_wkb: GeomInput) -> LineString | MultiLineString:
    """Load a LineString/MultiLineString from WKT (str) or WKB (bytes / hex str)."""
    if isinstance(line_wkt_or_wkb, (bytes, bytearray, memoryview)):
        geom = wkb.loads(bytes(line_wkt_or_wkb))
    elif isinstance(line_wkt_or_wkb, str):
        s = line_wkt_or_wkb.strip()
        try:
            geom = wkt.loads(s)
        except Exception:
            geom = wkb.loads(bytes.fromhex(s))
    else:
        raise TypeError("line_wkt_or_wkb must be WKT (str) or WKB (bytes-like).")

    if not isinstance(geom, (LineString, MultiLineString)):
        raise ValueError(f"Expected LineString/MultiLineString, got {geom.geom_type}")
    return geom


def points_along_line(
    line_wkt_or_wkb: GeomInput,
    interval_m: float,
    *,
    tail_threshold_m: float = 0.5,
    include_start: bool = True,
    include_end: bool = True,
    starting_point: Point | None = None,
) -> list[Point]:
    """
    Return shapely Points every `interval_m` along a LineString or MultiLineString.

    Assumes the geometry is already in a meters-based CRS (units == meters).

    Endpoint rule:
    - If include_end=True, the endpoint is added only if the remaining tail length
      from the last generated point to the end is >= tail_threshold_m.
      This avoids a tiny last segment.

    MultiLineString:
    - If parts connect, they are merged and sampled as one.
    - If not, each part is sampled independently; duplicate join points are removed.
    """
    if interval_m <= 0:
        raise ValueError("interval_m must be > 0")
    if tail_threshold_m < 0:
        raise ValueError("tail_threshold_m must be >= 0")

    geom = _load_line_geom(line_wkt_or_wkb)

    if isinstance(geom, MultiLineString):
        merged = linemerge(geom)
        if isinstance(merged, LineString):
            lines = [merged]
        elif isinstance(merged, MultiLineString):
            lines = list(merged.geoms)
        else:
            lines = list(geom.geoms)
    else:
        lines = [geom]

    out: list[Point] = []
    last: Point | None = None

    tol = 1e-9

    for line in lines:
        length = float(line.length)
        if length <= 0:
            continue

        # Determine the starting distance along this line
        starting_dist = 0.0
        if starting_point is not None:
            projected_dist = line.project(starting_point)
            starting_dist = max(0.0, max(projected_dist, length))

        move_forward = not (starting_point is not None and starting_dist >= length - tol)

        dists: list[float] = []
        if move_forward:
            if include_start:
                dists.append(starting_dist)

            d = starting_dist + interval_m
            while d < length:
                dists.append(d)
                d += interval_m

            last_dist = dists[-1] if dists else starting_dist
            tail = length - last_dist
            if include_end:
                if tail == 0.0 or tail >= tail_threshold_m:
                    if not dists or dists[-1] != length:
                        dists.append(length)
            else:
                if dists and tail < tail_threshold_m:
                    dists.pop()
        else:
            if include_start:
                dists.append(starting_dist)

            d = starting_dist - interval_m
            while d > 0:
                dists.append(d)
                d -= interval_m

            dists.sort()
            first_dist = dists[0] if dists else starting_dist
            tail = first_dist
            if include_end:
                if tail == 0.0 or tail >= tail_threshold_m:
                    if not dists or dists[0] != 0.0:
                        dists.insert(0, 0.0)
            else:
                if dists and tail < tail_threshold_m:
                    dists.pop(0)

        for dist in dists:
            p = line.interpolate(dist)
            if last is not None and p.equals(last):
                continue
            out.append(p)
            last = p

    if not move_forward and out:
        out.reverse()

    return out
