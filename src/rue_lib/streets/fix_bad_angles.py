"""Fix bad cells with insufficient right angles by adjusting geometry.

This module handles cells that are too small and have only 2 right angles,
attempting to fix them by splitting non-orthogonal edges and adjusting them
to create proper rectangular cells.
"""

from __future__ import annotations

import math

from shapely.geometry import LineString
from shapely.ops import split as shapely_split

from rue_lib.streets.cell import Cell
from rue_lib.streets.grids import is_good_cell


def calculate_angle_degrees(v1: tuple[float, float], v2: tuple[float, float]) -> float:
    """Calculate angle between two vectors in degrees.

    Args:
        v1: First vector as (dx, dy)
        v2: Second vector as (dx, dy)

    Returns:
        Angle in degrees (0-180)
    """
    len1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2)
    len2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2)

    if len1 < 1e-6 or len2 < 1e-6:
        return 0.0

    dot = v1[0] * v2[0] + v1[1] * v2[1]
    cos_angle = max(-1.0, min(1.0, dot / (len1 * len2)))
    angle_rad = math.acos(cos_angle)
    return math.degrees(angle_rad)


def is_right_angle(
    v1: tuple[float, float], v2: tuple[float, float], tolerance: float = 0.1
) -> bool:
    """Check if angle between two vectors is approximately 90 degrees.

    Args:
        v1: First vector as (dx, dy)
        v2: Second vector as (dx, dy)
        tolerance: Cosine tolerance (0.1 corresponds to ~84-96 degrees)

    Returns:
        True if angle is approximately 90 degrees
    """
    len1 = math.sqrt(v1[0] ** 2 + v1[1] ** 2)
    len2 = math.sqrt(v2[0] ** 2 + v2[1] ** 2)

    if len1 < 1e-6 or len2 < 1e-6:
        return False

    dot = v1[0] * v2[0] + v1[1] * v2[1]
    cos_angle = dot / (len1 * len2)
    return abs(cos_angle) <= tolerance


def find_non_orthogonal_edges(coords: list[tuple[float, float]]) -> list[int]:
    """Find indices of corners where the angle is not 90 degrees.

    Args:
        coords: List of polygon coordinates (without closing point)

    Returns:
        List of corner indices where angle is not orthogonal
    """
    non_orthogonal = []

    for i in range(len(coords)):
        p0 = coords[i - 1]
        p1 = coords[i]
        p2 = coords[(i + 1) % len(coords)]

        v1 = (p1[0] - p0[0], p1[1] - p0[1])
        v2 = (p2[0] - p1[0], p2[1] - p1[1])

        if not is_right_angle(v1, v2):
            non_orthogonal.append(i)

    return non_orthogonal


def rotate_edge_to_parallel(
    edge_start: tuple[float, float],
    edge_end: tuple[float, float],
    reference_dir: tuple[float, float],
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Rotate an edge around its midpoint to make it parallel to reference direction.

    Args:
        edge_start: Start point of edge to rotate
        edge_end: End point of edge to rotate
        reference_dir: Reference direction vector (edge should be parallel to this)

    Returns:
        New (start, end) points of the rotated edge
    """
    # Find midpoint
    mid_x = (edge_start[0] + edge_end[0]) / 2.0
    mid_y = (edge_start[1] + edge_end[1]) / 2.0

    # Calculate edge length
    dx = edge_end[0] - edge_start[0]
    dy = edge_end[1] - edge_start[1]
    edge_len = math.sqrt(dx * dx + dy * dy)

    if edge_len < 1e-6:
        return edge_start, edge_end

    # Normalize reference direction
    ref_len = math.sqrt(reference_dir[0] ** 2 + reference_dir[1] ** 2)
    if ref_len < 1e-6:
        return edge_start, edge_end

    ref_unit = (reference_dir[0] / ref_len, reference_dir[1] / ref_len)

    # Create new edge aligned with reference direction, same length, centered at midpoint
    half_len = edge_len / 2.0
    new_start = (float(mid_x - ref_unit[0] * half_len), float(mid_y - ref_unit[1] * half_len))
    new_end = (float(mid_x + ref_unit[0] * half_len), float(mid_y + ref_unit[1] * half_len))

    return new_start, new_end


def fix_area_too_small_bad_angles(
    cell: Cell,
    target_area: float,
    all_cells: list[Cell] | None = None,
    cluster_width: float = 20.0,
) -> tuple[list[Cell], int | None]:
    """Fix cells with area_too_small and only 2 right angles.

    Args:
        cell: Cell to fix (must have quality["reason"] == "area_too_small")
        target_area: Target area in square meters
        all_cells: Optional list of all cells for neighbor-aware fixing
        cluster_width: Width of the cluster/split line in meters

    Returns:
        Tuple of (list containing fixed cell or original cell, neighbor index that was updated)
    """
    quality = cell.quality

    # Only handle area_too_small with exactly 2 right angles
    if quality.get("reason") != "area_too_small":
        return [cell], None

    if quality.get("right_angles", 0) != 2:
        return [cell], None

    print(f"[fix-angles] Attempting to fix cell with 2 right angles (area: {cell.geom.area:.1f}m²)")

    coords = list(cell.geom.exterior.coords)[:-1]

    non_ortho_indices = find_non_orthogonal_edges(coords)

    if len(non_ortho_indices) < 2:
        print(f"[fix-angles]   Expected 2 non-orthogonal corners, found {len(non_ortho_indices)}")
        return [cell], None

    # Check if it's a quadrilateral
    if len(coords) != 4:
        print(f"[fix-angles]   Cell is not a quadrilateral ({len(coords)} vertices)")
        return [cell], None

    # Check if good corners are adjacent
    right_angle_indices = [i for i in range(4) if i not in non_ortho_indices]

    if len(right_angle_indices) != 2:
        print(f"[fix-angles]   Expected 2 right-angle corners, found {len(right_angle_indices)}")
        return [cell], None

    good_idx_0 = right_angle_indices[0]
    good_idx_1 = right_angle_indices[1]
    diff = abs(good_idx_1 - good_idx_0)
    are_good_adjacent = (diff == 1) or (diff == 3)

    if not are_good_adjacent:
        print("[fix-angles]   Good corners are not adjacent (cannot fix with this method)")
        return [cell], None

    print(f"[fix-angles]   Non-orthogonal corners at indices: {non_ortho_indices}")
    print(f"[fix-angles]   Good corners at indices {right_angle_indices} are adjacent")

    # Get the longest edge to use for splitting
    edge_lengths = []
    for i in range(4):
        p1 = coords[i]
        p2 = coords[(i + 1) % 4]
        dx = p2[0] - p1[0]
        dy = p2[1] - p1[1]
        length = math.sqrt(dx * dx + dy * dy)
        edge_lengths.append((i, length, p1, p2))

    # Find longest edge
    longest_edge = max(edge_lengths, key=lambda x: x[1])
    edge_idx, edge_len, p1, p2 = longest_edge

    print(f"[fix-angles]   Longest edge is edge {edge_idx} with length {edge_len:.1f}m")

    # Calculate circular angles at both corners of the longest edge
    prev_corner = coords[(edge_idx - 1) % 4]
    next_corner = coords[(edge_idx + 2) % 4]
    angle_to_prev_p1 = math.atan2(prev_corner[1] - p1[1], prev_corner[0] - p1[0])
    angle_to_p2_from_p1 = math.atan2(p2[1] - p1[1], p2[0] - p1[0])
    circular_angle_p1 = angle_to_p2_from_p1 - angle_to_prev_p1
    if circular_angle_p1 < 0:
        circular_angle_p1 += 2 * math.pi
    circular_angle_p1_deg = math.degrees(circular_angle_p1)

    angle_to_p1_from_p2 = math.atan2(p1[1] - p2[1], p1[0] - p2[0])
    angle_to_next_p2 = math.atan2(next_corner[1] - p2[1], next_corner[0] - p2[0])
    circular_angle_p2 = angle_to_next_p2 - angle_to_p1_from_p2
    if circular_angle_p2 < 0:
        circular_angle_p2 += 2 * math.pi
    circular_angle_p2_deg = math.degrees(circular_angle_p2)

    # Use the corner with the smaller circular angle
    if circular_angle_p1_deg < circular_angle_p2_deg:
        corner = p1
        circular_angle_deg = circular_angle_p1_deg
        other_corner = p2
    else:
        corner = p2
        circular_angle_deg = circular_angle_p2_deg
        other_corner = p1

    print(f"[fix-angles]   Using circular angle: {circular_angle_deg:.1f}°")

    # Calculate direction along the longest edge
    edge_dir_x = (other_corner[0] - corner[0]) / edge_len
    edge_dir_y = (other_corner[1] - corner[1]) / edge_len

    # Create perpendicular direction to the longest edge
    # Perpendicular to (dx, dy) is (-dy, dx) or (dy, -dx)
    perp_dir_x = -edge_dir_y
    perp_dir_y = edge_dir_x

    # We need to find a point on the longest edge such that when we draw a perpendicular
    # line from it, the segment inside the cell has length = cluster_width
    # Using the circular angle and trigonometry:
    # distance_along_edge = cluster_width / tan(circular_angle_rad)
    circular_angle_rad = math.radians(circular_angle_deg)
    distance_along_edge = cluster_width / math.tan(circular_angle_rad)

    print(f"[fix-angles]   Distance along edge: {distance_along_edge:.1f}m")

    # Find point on longest edge at this distance from corner
    point_on_edge = (
        float(corner[0] + edge_dir_x * distance_along_edge),
        float(corner[1] + edge_dir_y * distance_along_edge),
    )

    # Extend the perpendicular line in both directions to ensure it crosses the polygon
    extension = cluster_width * 2
    extended_start = (
        point_on_edge[0] - perp_dir_x * extension,
        point_on_edge[1] - perp_dir_y * extension,
    )
    extended_end = (
        point_on_edge[0] + perp_dir_x * extension,
        point_on_edge[1] + perp_dir_y * extension,
    )

    split_line = LineString([extended_start, extended_end])

    # Check if split line intersects with a bad neighbor cell
    if all_cells is not None:
        bad_neighbor_idx = None
        for i, other_cell in enumerate(all_cells):
            if other_cell.id == cell.id:
                continue

            # Check if this neighbor is bad (not good)
            if other_cell.quality.get("is_good", False):
                continue

            # Check if split line intersects with neighbor
            if other_cell.geom.intersects(split_line):
                # Check if they share an edge (not just touch at a point)
                intersection = cell.geom.intersection(other_cell.geom)
                if intersection.geom_type == "LineString" and intersection.length > 1.0:
                    bad_neighbor_idx = i
                    break

        if bad_neighbor_idx is None:
            print("[fix-angles]   No bad neighbor found along split line, skipping fix")
            return [cell], None

        print(f"[fix-angles]   Found bad neighbor cell {all_cells[bad_neighbor_idx].id}")

    try:
        result = shapely_split(cell.geom, split_line)
        parts = list(result.geoms) if hasattr(result, "geoms") else [result]

        if len(parts) > 1:
            print(f"[fix-angles]   Split cell into {len(parts)} pieces")

            # Identify biggest and smallest parts
            biggest_part = max(parts, key=lambda g: g.area)
            smallest_part = min(parts, key=lambda g: g.area)

            # Create new cell with biggest part
            new_quality = is_good_cell(biggest_part, target_area)
            new_cell = Cell(id=cell.id, geom=biggest_part, quality=new_quality)

            print(
                f"[fix-angles]   Biggest piece: area={biggest_part.area:.1f}m², "
                f"right_angles={new_quality.get('right_angles', 0)}"
            )

            # If we have all_cells, merge smallest part with bad neighbor
            if all_cells is not None and bad_neighbor_idx is not None:
                # Merge geometries and ensure valid polygon
                neighbor_cell = all_cells[bad_neighbor_idx]
                merged_geom = neighbor_cell.geom.union(smallest_part)

                # If union creates a GeometryCollection or MultiPolygon, extract the main polygon
                if merged_geom.geom_type == "GeometryCollection":
                    # Get the largest polygon from the collection
                    polygons = [g for g in merged_geom.geoms if g.geom_type == "Polygon"]
                    if polygons:
                        merged_geom = max(polygons, key=lambda p: p.area)
                elif merged_geom.geom_type == "MultiPolygon":
                    # Get the largest polygon from the multipolygon
                    merged_geom = max(merged_geom.geoms, key=lambda p: p.area)

                # Simplify to remove unnecessary vertices and ensure valid geometry
                merged_geom = merged_geom.buffer(0)

                merged_quality = is_good_cell(merged_geom, target_area)
                all_cells[bad_neighbor_idx] = Cell(
                    id=neighbor_cell.id, geom=merged_geom, quality=merged_quality
                )
                print(
                    f"[fix-angles]   Merged small piece ({smallest_part.area:.1f}m²) "
                    f"with neighbor {neighbor_cell.id}, new area={merged_geom.area:.1f}m²"
                )

            return [new_cell], bad_neighbor_idx
        else:
            print("[fix-angles]   Split did not create multiple pieces")
            return [cell], None

    except Exception as e:
        print(f"[fix-angles]   Split failed: {e}")
        return [cell], None


def inspect_and_fix_bad_angles(
    cells: list[Cell],
    target_area: float,
    cluster_width: float = 20.0,
) -> list[Cell]:
    """Inspect and fix cells with bad angles (area_too_small with 2 right angles).

    Args:
        cells: List of cells to inspect and fix
        target_area: Target area in square meters
        cluster_width: Width for splitting operations in meters

    Returns:
        List of cells after fixing attempts
    """
    print("\n[fix-angles] Inspecting cells for bad angles...")

    # Count cells with area_too_small and 2 right angles
    bad_angle_cells = []
    for i, cell in enumerate(cells):
        quality = cell.quality
        if quality.get("reason") == "area_too_small" and quality.get("right_angles", 0) <= 2:
            bad_angle_cells.append((i, cell))

    if not bad_angle_cells:
        print("[fix-angles] No cells found with area_too_small and 2 right angles")
        return cells

    print(f"[fix-angles] Found {len(bad_angle_cells)} cells with area_too_small and 2 right angles")

    working_cells = list(cells)
    fixed_count = 0

    for original_idx, cell in bad_angle_cells:
        print(f"\n[fix-angles] Processing cell {original_idx}...")

        current_idx = None
        for i, c in enumerate(working_cells):
            if c.id == cell.id:
                current_idx = i
                break

        if current_idx is None:
            print("[fix-angles]   Cell not found in working list (may have been merged)")
            continue

        cell = working_cells[current_idx]

        result, updated_neighbor_idx = fix_area_too_small_bad_angles(
            cell,
            target_area,
            all_cells=working_cells,
            cluster_width=cluster_width,
        )

        # If a neighbor was updated, refresh our reference to it in bad_angle_cells
        if updated_neighbor_idx is not None:
            # Update any cells in bad_angle_cells that point to the updated neighbor
            updated_neighbor = working_cells[updated_neighbor_idx]
            for i, (orig_idx, bad_cell) in enumerate(bad_angle_cells):
                if bad_cell.id == updated_neighbor.id:
                    # Update the reference in the list
                    bad_angle_cells[i] = (orig_idx, updated_neighbor)
                    print(
                        f"[fix-angles]   Updated neighbor cell {updated_neighbor.id} in "
                        f"processing queue"
                    )
                    break

        # Update the current cell in working_cells
        if len(result) == 1:
            if result[0].quality.get("right_angles", 0) > cell.quality.get("right_angles", 0):
                fixed_count += 1
                working_cells[current_idx] = result[0]
                print("[fix-angles]   ✓ Cell improved")
            elif result[0].id == cell.id and result[0].geom == cell.geom:
                print("[fix-angles]   Cell unchanged (no bad neighbor or unfixable)")
            else:
                # Cell was updated but didn't improve right angles
                working_cells[current_idx] = result[0]
        elif len(result) != 1:
            print(f"[fix-angles]   Cell was split into {len(result)} pieces")
            working_cells.pop(current_idx)
            for split_cell in reversed(result):
                working_cells.insert(current_idx, split_cell)
            fixed_count += 1

    print("\n[fix-angles] Summary:")
    print(f"[fix-angles]   Processed: {len(bad_angle_cells)} cells")
    print(f"[fix-angles]   Fixed: {fixed_count} cells")

    return working_cells
