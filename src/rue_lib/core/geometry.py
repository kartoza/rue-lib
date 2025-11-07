"""Basic geometry operations and cleaning utilities."""

from __future__ import annotations

import math

from shapely.geometry import Polygon

from rue_lib.core.exceptions import GeometryError


def clean_edges(poly: Polygon) -> Polygon:
    """
    Fix self-intersections and dangles using buffer(0).

    Args:
        poly: Input polygon

    Returns:
        Cleaned polygon
    """
    try:
        return poly.buffer(0)
    except Exception as e:
        raise GeometryError(f"Failed to clean edges: {e}") from e


def clean_angles(poly: Polygon, eps: float) -> Polygon:
    """
    Remove tiny spikes and near-collinear vertices.

    Args:
        poly: Input polygon
        eps: Simplification tolerance in meters

    Returns:
        Simplified polygon
    """
    try:
        simplified = poly.simplify(eps, preserve_topology=True)
        return simplified.buffer(0)
    except Exception as e:
        raise GeometryError(f"Failed to clean angles: {e}") from e


def bounce(poly: Polygon, dist: float) -> Polygon:
    """
    Offset inwards and back out to tidy jaggies (miter/square caps).

    Args:
        poly: Input polygon
        dist: Offset distance in meters

    Returns:
        Bounced polygon
    """
    if dist <= 0:
        return poly

    try:
        return poly.buffer(-dist, join_style=2, cap_style=2).buffer(dist, join_style=2, cap_style=2)
    except Exception as e:
        raise GeometryError(f"Failed to bounce polygon: {e}") from e


def ring_clean_short_edges(
    coords: list[tuple[float, float]], min_len: float
) -> list[tuple[float, float]]:
    """
    Remove vertices that create edges shorter than min_len.

    Args:
        coords: Ring coordinates (should be closed)
        min_len: Minimum edge length in meters

    Returns:
        Cleaned coordinate list
    """
    if len(coords) < 4:
        return coords

    # Ensure closure
    if coords[0] != coords[-1]:
        coords = coords + [coords[0]]

    out: list[tuple[float, float]] = [coords[0]]

    for i in range(1, len(coords) - 1):
        x0, y0 = out[-1]
        x1, y1 = coords[i]
        seg_len = math.hypot(x1 - x0, y1 - y0)

        if seg_len >= min_len:
            out.append((x1, y1))

    # Close the ring
    if out[0] != out[-1]:
        out.append(out[0])

    if len(out) < 4:
        return coords

    return out


def ring_clean_collinear(
    coords: list[tuple[float, float]], dot_thresh: float = 0.9999
) -> list[tuple[float, float]]:
    """
    Remove nearly collinear vertices using dot product threshold.

    Args:
        coords: Ring coordinates (should be closed)
        dot_thresh: Collinearity threshold (higher = more aggressive)

    Returns:
        Cleaned coordinate list
    """
    if len(coords) < 4:
        return coords

    if coords[0] != coords[-1]:
        coords = coords + [coords[0]]

    keep: list[tuple[float, float]] = []
    n = len(coords) - 1

    for i in range(n):
        p_prev = coords[i - 1 if i > 0 else n - 1]
        p = coords[i]
        p_next = coords[(i + 1) % n]

        v0 = (p[0] - p_prev[0], p[1] - p_prev[1])
        v1 = (p_next[0] - p[0], p_next[1] - p[1])

        n0 = math.hypot(v0[0], v0[1])
        n1 = math.hypot(v1[0], v1[1])

        if n0 == 0 or n1 == 0:
            continue

        v0 = (v0[0] / n0, v0[1] / n0)
        v1 = (v1[0] / n1, v1[1] / n1)

        dot = abs(v0[0] * v1[0] + v0[1] * v1[1])

        if dot < dot_thresh:
            keep.append(p)

    if not keep:
        return coords

    if keep[0] != keep[-1]:
        keep.append(keep[0])

    if len(keep) < 4:
        return coords

    return keep
