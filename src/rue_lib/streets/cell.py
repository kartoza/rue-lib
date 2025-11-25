from __future__ import annotations

from collections.abc import Sequence

from shapely.geometry import LineString, Polygon
from shapely.ops import split


class Cell:
    """Represents a grid cell with its geometry and quality information.

    Attributes:
        id (int): Unique identifier for the cell.
        quality (dict): Quality assessment dictionary containing metrics like
            area_ratio, num_vertices, right_angles, is_good, and reason.
        geom (Polygon): Shapely Polygon representing the cell's geometry.
    """

    id: int = 0
    quality: dict = {}
    geom: Polygon = None

    def __init__(self, id: int, geom: Polygon, quality: dict):
        """Initialize a Cell instance.

        Args:
            id (int): Unique identifier for the cell.
            geom (Polygon): Shapely Polygon representing the cell's geometry.
            quality (dict): Quality assessment dictionary.
        """
        self.id = id
        self.geom = geom
        self.quality = quality


def get_neighbors(cell: Cell, all_cells: list[Cell]) -> list[Cell]:
    """Find all neighboring cells that touch or intersect the given cell.

    A neighbor is defined as any cell whose geometry either touches or intersects
    the given cell's geometry, excluding the cell itself.

    Args:
        cell (Cell): The cell to find neighbors for.
        all_cells (list[Cell]): list of all cells in the grid to search.

    Returns:
        list[Cell]: list of neighboring cells. Returns empty list if no neighbors found.
    """
    neighbors = []
    for _cell in all_cells:
        _cell_geom = _cell.geom
        if cell.geom.equals(_cell_geom):
            continue
        if cell.geom.touches(_cell_geom) or cell.geom.intersects(_cell_geom):
            neighbors.append(_cell)
    return neighbors


def find_neighbor_corner_split_lines(
    bad_cell: Cell,
    all_cells: list[Cell],
) -> list[LineString]:
    """Find split lines based on corner points from neighboring good cells.

    Uses two strategies to generate potential split lines:
    1. Lines connecting pairs of neighbor corners across the bad cell
    2. Perpendicular lines extending from neighbor edge corners (creates rectangular splits)

    Each generated line is extended and tested to ensure it crosses the bad cell's
    interior with sufficient length (>1.0 units).

    Args:
        bad_cell (Cell): The oversized or malformed cell to generate split lines for.
        all_cells (list[Cell]): All cells in the grid, used to find neighbors.

    Returns:
        list[LineString]: list of potential split lines. Returns empty list if
            fewer than 2 neighbor corners are found or if no valid lines can be generated.
    """
    neighbors = get_neighbors(bad_cell, all_cells)

    neighbor_corners = []
    neighbor_cells = []

    for neighbor in neighbors:
        neighbor_geom = neighbor.geom
        coords = list(neighbor_geom.exterior.coords)[:-1]
        for coord in coords:
            neighbor_corners.append((coord[0], coord[1]))
        neighbor_cells.append((neighbor, coords))

    if len(neighbor_corners) < 2:
        return []

    unique_corners = list(set(neighbor_corners))

    split_lines = []
    bad_cell_buffered = bad_cell.geom.buffer(0.01)

    # STRATEGY 1: Lines connecting pairs of neighbor corners
    for i, corner1 in enumerate(unique_corners):
        for corner2 in unique_corners[i + 1 :]:
            dx = corner2[0] - corner1[0]
            dy = corner2[1] - corner1[1]
            length = (dx**2 + dy**2) ** 0.5

            if length < 1.0:
                continue

            extension = 200.0
            dir_x = dx / length
            dir_y = dy / length

            p1_ext = (corner1[0] - extension * dir_x, corner1[1] - extension * dir_y)
            p2_ext = (corner2[0] + extension * dir_x, corner2[1] + extension * dir_y)

            line = LineString([p1_ext, p2_ext])

            if bad_cell_buffered.intersects(line):
                intersection = bad_cell.geom.intersection(line)
                if not intersection.is_empty and intersection.length > 1.0:
                    split_lines.append(line)

    # STRATEGY 2: Perpendicular lines from neighbor corners
    for _neighbor_cell, neighbor_coords in neighbor_cells:
        for i in range(len(neighbor_coords)):
            p1 = neighbor_coords[i]
            p2 = neighbor_coords[(i + 1) % len(neighbor_coords)]

            dx = p2[0] - p1[0]
            dy = p2[1] - p1[1]
            edge_length = (dx**2 + dy**2) ** 0.5

            if edge_length < 1.0:
                continue

            edge_dir_x = dx / edge_length
            edge_dir_y = dy / edge_length

            # Calculate perpendicular direction (both ways)
            perp_dirs = [(-edge_dir_y, edge_dir_x), (edge_dir_y, -edge_dir_x)]

            # Try perpendicular lines from both corners of this edge
            for corner in [p1, p2]:
                for perp_dir_x, perp_dir_y in perp_dirs:
                    extension = 200.0
                    line_p1 = (
                        corner[0] - extension * perp_dir_x,
                        corner[1] - extension * perp_dir_y,
                    )
                    line_p2 = (
                        corner[0] + extension * perp_dir_x,
                        corner[1] + extension * perp_dir_y,
                    )

                    line = LineString([line_p1, line_p2])

                    # Check if this line crosses the bad cell
                    if bad_cell_buffered.intersects(line):
                        intersection = bad_cell.geom.intersection(line)
                        if not intersection.is_empty and intersection.length > 1.0:
                            split_lines.append(line)

    return split_lines


def fix_area_too_large(
    cell: Cell,
    target_area: float,
    all_cells: list[Cell] = None,
) -> list[Cell]:
    """Attempt to split an oversized cell into good cells using neighbor-aware split lines.

    Uses neighbor corner points to generate split lines that create
    rectangular cells consistent with the surrounding grid. The function:
    1. Finds split lines from neighboring good cell corners (if neighbors provided)
    2. Evaluates each split line using a multi-factor scoring system
    3. Selects the split that produces the highest quality cells
    4. Recursively splits any remaining oversized pieces

    The scoring system awards points for:
    - Creating fully "good" cells (+100 each)
    - Right angles (+10 per angle)
    - Quadrilateral shapes (+20)
    - Low area deviation (-50 × deviation)
    - Even splits across all pieces (+50 bonus)

    Args:
        cell (Cell): The oversized cell to fix. Must have quality["reason"] == "area_too_large".
        target_area (float): Target area for cells in square meters (e.g., 10000.0 for 100m × 100m).
        all_cells (list[Cell], optional): All cells in the grid. Required for neighbor-aware
            splitting. If None, function returns the cell unchanged.

    Returns:
        list[Cell]: list of split cells if successful. Returns [cell] unchanged if:
            - Cell's quality reason is not "area_too_large"
            - No valid split lines found
            - No split improves the quality score
    """
    from rue_lib.streets.grids import is_good_cell

    if cell.quality.get("reason") != "area_too_large":
        return [cell]

    neighbor_corner_lines = []
    best_splits = None

    if all_cells is not None:
        neighbor_corner_lines = find_neighbor_corner_split_lines(cell, all_cells)
        if neighbor_corner_lines:
            print(
                f"[cell-fix]   Found {len(neighbor_corner_lines)} "
                f"potential split lines from neighbor corners"
            )

    if neighbor_corner_lines:
        best_score = -float("inf")

        for split_line in neighbor_corner_lines:
            try:
                result = split(cell.geom, split_line)

                if result.is_empty or len(result.geoms) < 2:
                    continue

                split_cells = list(result.geoms)

                min_area_threshold = target_area * 0.1
                split_cells = [
                    c
                    for c in split_cells
                    if c.geom_type == "Polygon" and c.area >= min_area_threshold
                ]

                if len(split_cells) < 2:
                    continue

                # Calculate quality score for this split
                score = 0
                total_area_deviation = 0
                total_right_angles = 0
                good_cell_count = 0

                for split_cell in split_cells:
                    split_quality = is_good_cell(split_cell, target_area)

                    if split_quality.get("is_good"):
                        good_cell_count += 1
                        score += 100

                    ra = split_quality.get("right_angles", 0)
                    total_right_angles += ra
                    score += ra * 10

                    area_ratio = split_quality.get("area_ratio", 1.0)
                    area_deviation = abs(1.0 - area_ratio)
                    total_area_deviation += area_deviation
                    score -= area_deviation * 50

                    if split_quality.get("num_vertices") == 4:
                        score += 20

                # Average area deviation (lower is better)
                avg_area_deviation = total_area_deviation / len(split_cells)
                if avg_area_deviation < 0.2:
                    score += 50

                if score > best_score:
                    best_score = score
                    best_splits = split_cells
            except Exception as e:
                print(e)
                continue

        if best_splits and len(best_splits) > 1:
            min_area = min(c.area for c in best_splits)
            print(
                f"[cell-fix]   Successfully split using neighbor corner points "
                f"(score: {best_score:.1f}, "
                f"min area: {min_area:.1f}m²)"
            )
            final_cells = []
            for split_cell in best_splits:
                split_quality = is_good_cell(split_cell, target_area)
                _split_cell = Cell(id=len(all_cells), quality=split_quality, geom=split_cell)

                if split_quality.get("reason") == "area_too_large":
                    sub_splits = fix_area_too_large(
                        _split_cell,
                        target_area,
                        all_cells=None,
                    )
                    final_cells.extend(sub_splits)
                else:
                    final_cells.append(_split_cell)

            return final_cells

    if best_splits:
        final_cells = []

        for split_cell in best_splits:
            if split_cell.geom_type != "Polygon":
                continue

            split_quality = is_good_cell(split_cell, target_area)
            updated_split_cell = Cell(id=len(all_cells), quality=split_quality, geom=split_cell)

            if split_quality.get("reason") == "area_too_large":
                sub_splits = fix_area_too_large(updated_split_cell, target_area)
                final_cells.extend(sub_splits)
            else:
                final_cells.append(updated_split_cell)

        return final_cells
    else:
        return [cell]


def inspect_and_fix_cells(
    cells: Sequence[Cell],
    target_area: float,
    fix_oversized: bool = True,
) -> list[Cell]:
    """Inspect cells and attempt to fix cells using sequential fixing.

    This function performs a comprehensive analysis and fixing process for a grid of cells:
    1. Counts the number of bad cells (cells not meeting quality criteria)
    2. If fix_oversized=True, sequentially fixes oversized cells one at a time
    3. Uses neighbor-aware splitting to create consistent rectangular cells
    4. Updates the working cell list after each fix so later fixes benefit from earlier fixes
    5. Reports before/after statistics

    Sequential fixing is key: fixing cells in order allows later cells to use earlier
    fixed cells as good neighbors, dramatically improving split quality.

    Args:
        cells (Sequence[Cell]): Sequence of Cell objects to inspect and fix.
        target_area (float): Target area for cells in square meters (e.g., 10000.0 for 100m × 100m).
            Used to evaluate cell quality and determine split positions.
        fix_oversized (bool, optional): Whether to attempt fixing oversized cells. Defaults to True.
            If False, returns cells unchanged after counting bad cells.

    Returns:
        list[Cell]: list of cells after fixing. May contain more cells than input if oversized
            cells were successfully split. Cells that couldn't be fixed are marked with
            reason "area_too_large_unfixable" to prevent infinite loops.
    """
    from rue_lib.streets.grids import is_good_cell

    total_bad = 0

    for cell in cells:
        if not cell.quality.get("is_good", True):
            total_bad += 1

    if fix_oversized:
        print("\n[cell-fix] Attempting to fix oversized cells...")

        working_cells = list(cells)
        fixed_count = 0

        while True:
            bad_idx = None
            for i, q in enumerate(working_cells):
                try:
                    if q.quality.get("reason") == "area_too_large":
                        bad_idx = i
                        break
                except AttributeError:
                    pass

            if bad_idx is None:
                break

            print(f"\n[cell-fix] Fixing cell {bad_idx}...")

            cell = working_cells[bad_idx]

            split_result = fix_area_too_large(
                cell,
                target_area,
                all_cells=working_cells,
            )

            if len(split_result) > 1:
                fixed_count += 1
                print(f"[cell-fix]   Cell {bad_idx}: Split into {len(split_result)} pieces")

                working_cells.pop(bad_idx)

                for split_cell in reversed(split_result):
                    working_cells.insert(bad_idx, split_cell)
            else:
                print(f"[cell-fix]   Cell {bad_idx}: Unable to fix - marking as processed")
                working_cells[bad_idx].quality["reason"] = "area_too_large_unfixable"

        fixed_cells = working_cells

        print(f"[cell-fix] Fixed {fixed_count} oversized cells")

        fixed_qualities = [is_good_cell(cell.geom, target_area) for cell in fixed_cells]

        total_bad_after = sum(1 for q in fixed_qualities if not q.get("is_good", True))
        print("\n[cell-fix] After fixing:")
        print(f"[cell-fix]   Total cells: {len(fixed_cells)}")
        print(f"[cell-fix]   Bad cells: {total_bad_after}")
        print(f"[cell-fix]   Improvement: {total_bad - total_bad_after} fewer bad cells\n")

        return fixed_cells
    else:
        # No fixing requested, return inputs unchanged
        return list(cells)
