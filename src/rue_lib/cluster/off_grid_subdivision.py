# src/rue_lib/cluster/off_grid_subdivision.py
"""Subdivide off-grid areas into smaller plot clusters."""

from typing import Optional

import numpy as np
from shapely.geometry import Polygon


def compute_angle_dot(polygon: Polygon, vertex_idx: int) -> float:
    """
    Compute the dot product of normalized edge vectors at a vertex.

    This helps identify sharp corners (low dot product = sharp angle).

    Args:
        polygon: Input polygon
        vertex_idx: Index of vertex to check

    Returns:
        Dot product of normalized edge vectors (-1 to 1)
    """
    coords = list(polygon.exterior.coords)[:-1]  # Exclude duplicate last point
    n = len(coords)

    # Get three consecutive points
    p_prev = np.array(coords[(vertex_idx - 1) % n])
    p_curr = np.array(coords[vertex_idx % n])
    p_next = np.array(coords[(vertex_idx + 1) % n])

    # Edge vectors
    vec0 = p_curr - p_prev  # Incoming edge
    vec1 = p_next - p_curr  # Outgoing edge

    # Normalize
    vec0_norm = vec0 / np.linalg.norm(vec0) if np.linalg.norm(vec0) > 0 else vec0
    vec1_norm = vec1 / np.linalg.norm(vec1) if np.linalg.norm(vec1) > 0 else vec1

    # Dot product
    return np.dot(vec0_norm, vec1_norm)


def convert_to_quadrilateral(polygon: Polygon, min_area: float = 100.0) -> Optional[Polygon]:
    """
    Convert an off-grid polygon to a quadrilateral by removing vertices with sharpest angles.

    This simplifies complex polygons to rectangles/quads suitable for grid subdivision.

    Args:
        polygon: Input off-grid polygon
        min_area: Minimum area threshold (m²)

    Returns:
        Quadrilateral polygon or None if area too small or too few vertices
    """
    if polygon.area < min_area:
        return None

    coords = list(polygon.exterior.coords)[:-1]  # Exclude duplicate last point
    n_vertices = len(coords)

    if n_vertices < 4:
        return None

    # Compute angle dot products for all vertices
    # Lower dot product = sharper angle = keep this vertex
    # Higher dot product = smoother angle = remove this vertex
    vertex_dots = []
    for i in range(n_vertices):
        dot = compute_angle_dot(polygon, i)
        vertex_dots.append((i, dot))

    # Sort by dot product (descending) - vertices with highest dots are removed first
    sorted_vertices = sorted(vertex_dots, key=lambda x: x[1], reverse=True)

    # Keep only the 4 vertices with the lowest dot products (sharpest corners)
    vertices_to_keep = sorted([v[0] for v in sorted_vertices[-4:]])

    # Create new quadrilateral with these 4 vertices
    quad_coords = [coords[i] for i in vertices_to_keep]
    quad_coords.append(quad_coords[0])  # Close the polygon

    quad = Polygon(quad_coords)

    # Check if the worst corner (highest dot) is too smooth
    # If dot > -0.7, it's acceptable as a quad
    worst_dot = sorted_vertices[-1][1]

    if worst_dot > -0.7:
        return quad

    # If worst corner is too sharp, we might need to "chop" it
    # For now, return the quad and handle chopping later if needed
    return quad


def create_grid_positions(
    quad: Polygon,
    part_og_w: float,
    part_og_d: float,
    swap_orientation: bool = False,
) -> list[list[tuple[float, float]]]:
    """
    Create a grid of positions within a quadrilateral for plot subdivision.

    This generates a 2D grid of points that will be used to create individual plots.

    Args:
        quad: Quadrilateral polygon to subdivide
        part_og_w: Off-grid plot width (meters)
        part_og_d: Off-grid plot depth (meters)
        swap_orientation: If True, swap width and depth

    Returns:
        2D list of (x, y) positions: positions[i][j] = (x, y)
    """
    # Get the 4 edges and their lengths
    coords = list(quad.exterior.coords)[:-1]
    edges = []
    lengths = []

    for i in range(4):
        p0 = np.array(coords[i])
        p1 = np.array(coords[(i + 1) % 4])
        edge_vec = p1 - p0
        edge_len = np.linalg.norm(edge_vec)
        edges.append((p0, edge_vec))
        lengths.append(edge_len)

    # Determine orientation - find the edge most aligned with X-axis
    # This helps us orient the grid properly
    x_dots = []
    for i in range(4):
        edge_vec_norm = edges[i][1] / lengths[i] if lengths[i] > 0 else edges[i][1]
        dot = np.dot([1, 0], edge_vec_norm[:2])
        x_dots.append(dot)

    # Find the edge most aligned with +X direction
    sorted_indices = sorted(range(4), key=lambda i: x_dots[i])

    # Bottom edge (most aligned with +X)
    a = sorted_indices[-1]
    b = (a + 1) % 4
    c = (a + 2) % 4
    d = (a + 3) % 4

    # Get corner positions
    p0 = edges[a][0]  # Bottom-left
    p1 = edges[b][0]  # Bottom-right
    _p2 = edges[c][0]  # Top-right
    p3 = edges[d][0]  # Top-left

    # Width and depth dimensions
    wd = [part_og_w, part_og_d]
    if swap_orientation:
        wd = [part_og_d, part_og_w]

    # Calculate number of divisions
    # Use minimum of opposite edge lengths
    len_horizontal = min(lengths[a], lengths[c])
    len_vertical = min(lengths[b], lengths[d])

    num_x = int(np.floor((3 + len_horizontal) / wd[0])) + 1
    num_y = int(np.floor((3 + len_vertical) / wd[1])) + 1

    # Ensure at least 2 divisions in each direction
    if num_x == 1:
        num_x = 2
    if num_y == 1:
        num_y = 2

    # Create edge vectors
    vec_right0 = edges[a][1] / (num_x - 1)  # Bottom edge vector
    vec_right1 = -edges[c][1] / (num_x - 1)  # Top edge vector (reversed)

    # Extension for boundary (100m beyond edges)
    vec_right_exp = (
        vec_right1 / np.linalg.norm(vec_right1) * 100
        if np.linalg.norm(vec_right1) > 0
        else vec_right1
    )

    # Generate grid positions
    positions = []

    for i in range(num_x):
        column = []

        # Interpolate between bottom and top edges
        x0 = p0 + vec_right0 * i
        x1 = p3 + vec_right1 * i

        # Extend first and last columns beyond boundaries
        if i == 0:
            x0 = x0 - vec_right_exp
            x1 = x1 - vec_right_exp
        elif i == num_x - 1:
            x0 = x0 + vec_right_exp
            x1 = x1 + vec_right_exp

        # Vertical vector for this column
        vec_up = (x1 - x0) / (num_y - 1)
        vec_up_exp = vec_up / np.linalg.norm(vec_up) * 100 if np.linalg.norm(vec_up) > 0 else vec_up

        for j in range(num_y):
            xyz = x0 + vec_up * j

            # Extend first and last rows beyond boundaries
            if j == 0:
                xyz = xyz - vec_up_exp
            elif j == num_y - 1:
                xyz = xyz + vec_up_exp

            column.append(tuple(xyz[:2]))  # Store as (x, y)

        positions.append(column)

    return positions


def create_plot_grid_from_positions(
    positions: list[list[tuple[float, float]]],
    off_grid_boundary: Polygon,
) -> list[Polygon]:
    """
    Create individual plot polygons from grid positions.

    Creates rectangular plots by connecting adjacent grid points,
    then clips them to the off-grid boundary.

    Args:
        positions: 2D grid of positions from create_grid_positions
        off_grid_boundary: Boundary polygon to clip plots to

    Returns:
        List of plot polygons
    """
    plots = []

    num_x = len(positions)

    for i in range(num_x - 1):
        col0 = positions[i]
        col1 = positions[i + 1]
        num_y = len(col0)

        for j in range(num_y - 1):
            # Create quad from 4 adjacent points
            p0 = col0[j]
            p1 = col1[j]
            p2 = col1[j + 1]
            p3 = col0[j + 1]

            # Create polygon
            plot_coords = [p0, p1, p2, p3, p0]
            plot = Polygon(plot_coords)

            # Clip to off-grid boundary
            try:
                clipped = plot.intersection(off_grid_boundary)

                if not clipped.is_empty and clipped.area > 0:
                    if clipped.geom_type == "Polygon":
                        plots.append(clipped)
                    elif clipped.geom_type == "MultiPolygon":
                        # Take all valid polygons
                        plots.extend([p for p in clipped.geoms if p.area > 0])
            except Exception as e:
                print(f"Warning: Failed to clip plot at ({i}, {j}): {e}")
                continue

    return plots


def subdivide_off_grid(
    off_grid: Polygon,
    part_og_w: float = 140.0,
    part_og_d: float = 140.0,
    swap_orientation: bool = False,
    min_plot_area: float = None,
) -> list[Polygon]:
    """
    Subdivide an off-grid polygon into a grid of smaller plots.

    This is the main function that orchestrates the subdivision process:
    1. Simplify the off-grid to a quadrilateral
    2. Create a grid of positions
    3. Generate individual plot polygons
    4. Filter by minimum area

    Args:
        off_grid: Off-grid polygon to subdivide
        part_og_w: Plot width (meters)
        part_og_d: Plot depth (meters)
        swap_orientation: If True, swap width and depth
        min_plot_area: Minimum plot area (m²). Defaults to 30% of target plot size

    Returns:
        List of plot polygons

    Example:
        >>> off_grid = Polygon([(0, 0), (200, 0), (200, 200), (0, 200)])
        >>> plots = subdivide_off_grid(off_grid, part_og_w=50, part_og_d=50)
        >>> len(plots)  # Should create ~16 plots (4x4 grid)
        16
    """
    # Default minimum area: 30% of target plot size
    if min_plot_area is None:
        min_plot_area = part_og_w * part_og_d * 0.3

    # Step 1: Convert to quadrilateral
    quad = convert_to_quadrilateral(off_grid, min_area=100.0)

    if quad is None:
        # Off-grid too small or complex, return as single plot
        return [off_grid]

    # Step 2: Create grid positions
    try:
        positions = create_grid_positions(
            quad,
            part_og_w,
            part_og_d,
            swap_orientation,
        )
    except Exception as e:
        print(f"Warning: Failed to create grid positions: {e}")
        return [off_grid]

    # Step 3: Generate plot polygons
    try:
        plots = create_plot_grid_from_positions(positions, off_grid)
    except Exception as e:
        print(f"Warning: Failed to create plots from positions: {e}")
        return [off_grid]

    # Step 4: Filter by minimum area
    valid_plots = []
    small_plots = []

    for plot in plots:
        if plot.area >= min_plot_area:
            valid_plots.append(plot)
        else:
            small_plots.append(plot)

    # If no valid plots, return original
    if not valid_plots:
        return [off_grid]

    # For very small plots, could merge them into "park" areas
    # For now, we'll include them in the output
    # but mark them differently (caller can filter by area)

    return valid_plots + small_plots


def classify_plot_by_area(
    plot: Polygon,
    part_og_w: float,
    part_og_d: float,
    threshold: float = 0.3,
) -> str:
    """
    Classify a plot as 'plot' or 'park' based on area.

    Args:
        plot: Plot polygon
        part_og_w: Target plot width
        part_og_d: Target plot depth
        threshold: Area threshold (0.3 = 30% of target)

    Returns:
        'plot' or 'park'
    """
    target_area = part_og_w * part_og_d

    if plot.area >= target_area * threshold:
        return "plot"
    else:
        return "park"
