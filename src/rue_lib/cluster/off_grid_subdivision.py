# src/rue_lib/cluster/off_grid_subdivision.py
"""Subdivide off-grid areas into smaller plot clusters."""

from pathlib import Path

import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon

from rue_lib.cluster.classification import classify_plot_by_area
from rue_lib.cluster.helpers import convert_to_quadrilateral


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

    x_dots = []
    for i in range(4):
        edge_vec_norm = edges[i][1] / lengths[i] if lengths[i] > 0 else \
            edges[i][1]
        dot = np.dot([1, 0], edge_vec_norm[:2])
        x_dots.append(dot)

    sorted_indices = sorted(range(4), key=lambda i: x_dots[i])

    a = sorted_indices[-1]
    b = (a + 1) % 4
    c = (a + 2) % 4
    d = (a + 3) % 4

    # Get corner positions
    p0 = edges[a][0]
    p1 = edges[b][0]
    _p2 = edges[c][0]
    p3 = edges[d][0]

    wd = [part_og_w, part_og_d]
    if swap_orientation:
        wd = [part_og_d, part_og_w]

    len_horizontal = min(lengths[a], lengths[c])
    len_vertical = min(lengths[b], lengths[d])

    num_x = int(np.floor((3 + len_horizontal) / wd[0])) + 1
    num_y = int(np.floor((3 + len_vertical) / wd[1])) + 1

    if num_x == 1:
        num_x = 2
    if num_y == 1:
        num_y = 2

    vec_right0 = edges[a][1] / (num_x - 1)
    vec_right1 = -edges[c][1] / (num_x - 1)

    vec_right_exp = (
        vec_right1 / np.linalg.norm(vec_right1) * 100
        if np.linalg.norm(vec_right1) > 0
        else vec_right1
    )

    positions = []

    for i in range(num_x):
        column = []

        x0 = p0 + vec_right0 * i
        x1 = p3 + vec_right1 * i

        if i == 0:
            x0 = x0 - vec_right_exp
            x1 = x1 - vec_right_exp
        elif i == num_x - 1:
            x0 = x0 + vec_right_exp
            x1 = x1 + vec_right_exp

        vec_up = (x1 - x0) / (num_y - 1)
        vec_up_exp = vec_up / np.linalg.norm(vec_up) * 100 if np.linalg.norm(
            vec_up) > 0 else vec_up

        for j in range(num_y):
            xyz = x0 + vec_up * j

            if j == 0:
                xyz = xyz - vec_up_exp
            elif j == num_y - 1:
                xyz = xyz + vec_up_exp

            column.append(tuple(xyz[:2]))

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
            p0 = col0[j]
            p1 = col1[j]
            p2 = col1[j + 1]
            p3 = col0[j + 1]

            plot_coords = [p0, p1, p2, p3, p0]
            plot = Polygon(plot_coords)

            try:
                clipped = plot.intersection(off_grid_boundary)

                if not clipped.is_empty and clipped.area > 0:
                    if clipped.geom_type == "Polygon":
                        plots.append(clipped)
                    elif clipped.geom_type == "MultiPolygon":
                        plots.extend([p for p in clipped.geoms if p.area > 0])
            except Exception as e:
                print(f"Warning: Failed to clip plot at ({i}, {j}): {e}")
                continue

    return plots


def subdivide_off_grid(
        block_id: int,
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
    if min_plot_area is None:
        min_plot_area = part_og_w * part_og_d * 0.3
    quad = convert_to_quadrilateral(off_grid, min_area=100.0)

    if quad is None:
        print(f"Warning: No quad")
        return [off_grid]

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

    try:
        plots = create_plot_grid_from_positions(positions, off_grid)
    except Exception as e:
        print(f"Warning: Failed to create plots from positions: {e}")
        return [off_grid]

    valid_plots = []
    small_plots = []

    for plot in plots:
        if plot.area >= min_plot_area:
            valid_plots.append(plot)
        else:
            small_plots.append(plot)

    if not valid_plots:
        return [off_grid]

    return (valid_plots + small_plots)


def extract_off_grid_cluster(
        output_path: Path,
        off_grids_layer_name: str,
        part_og_w: float,
        part_og_d: float,
        output_layer_name: str,
        min_plot_area: float
):
    """
    Extract and subdivide off-grid areas into plot clusters.

    Reads off-grid geometries from a GeoPackage layer, subdivides each into
    individual plots, and classifies them by area. Prints progress information
    for each block processed.

    Args:
        output_path: Path to the GeoPackage file containing off-grid data
        off_grids_layer_name: Name of the layer containing off-grid geometries
        part_og_w: Target plot width in meters
        part_og_d: Target plot depth in meters
        output_layer_name: Layer name for the output plots

    Returns:
        None (currently collects plots but doesn't return them)
    """
    off_grids_layer = gpd.read_file(
        output_path, layer=off_grids_layer_name
    )

    all_plots = []

    for _idx, off_grid_part in off_grids_layer.iterrows():
        off_grid_geom = off_grid_part.geometry
        block_id = off_grid_part.get("block_id")
        parent_index = off_grid_part.get("part_index")

        print(f"  ---------------------------")
        print(f"  Block {block_id}:")
        print(f"    Off-grid area: {off_grid_geom.area:.2f} m²")

        try:
            # Subdivide the off-grid area using oriented approach
            plots = subdivide_off_grid(
                block_id,
                off_grid_geom,
                part_og_w=part_og_w,
                part_og_d=part_og_d,
                min_plot_area=min_plot_area
            )

            print(f"    ✓ Created {len(plots)} plots")

            # Collect plots
            for i, plot in enumerate(plots):
                feature = {
                    "geometry": plot,
                    "block_id": block_id,
                    "plot_type": classify_plot_by_area(
                        plot, part_og_w, part_og_d
                    ),
                    "plot_index": i,
                    "area": plot.area,
                    "parent_part": "off_grid",
                    "parent_index": parent_index,
                }
                try:
                    feature["color"] = off_grid_part["color"]
                    feature["type"] = off_grid_part["type"]
                except KeyError:
                    pass
                all_plots.append(feature)
        except Exception as e:
            print(f"    ✗ Error subdividing off-grid: {e}")

    gdf_out = gpd.GeoDataFrame(all_plots, crs=off_grids_layer.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")
    return output_layer_name
