import math
from pathlib import Path

import geopandas as gpd
import numpy as np
from osgeo import ogr
from scipy.spatial import Voronoi
from shapely import unary_union
from shapely.affinity import rotate
from shapely.errors import TopologicalError
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import nearest_points
from shapely.prepared import prep

from rue_lib.core.helpers import feature_geom_to_shapely
from rue_lib.streets.config import StreetConfig
from rue_lib.streets.operations import create_local_streets_zone


def is_good_cell(poly: Polygon, target_area) -> dict:
    """
    Check if a polygon is a "good cell" and return detailed quality information.

    Always checks all properties regardless of failures to provide complete information.

    Returns
    -------
    dict
        Dictionary with keys:
        - is_good (bool): True if cell meets all criteria
        - reason (str): Primary reason if not good, otherwise "good"
        - right_angles (int): Count of orthogonal corners (0-4)
        - num_vertices (int): Number of vertices
        - area_ratio (float): Actual area / target area
    """
    result = {
        "is_good": True,
        "reason": "good",
        "right_angles": 0,
        "num_vertices": 0,
        "area_ratio": 0.0,
    }

    if poly.geom_type != "Polygon":
        result["is_good"] = False
        result["reason"] = "not_polygon"
        return result

    area = poly.area
    result["area_ratio"] = area / target_area if target_area > 0 else 0.0

    max_area = target_area * 1.10
    min_area = target_area * 1

    if not (min_area <= area <= max_area):
        result["is_good"] = False
        if area < min_area:
            result["reason"] = "area_too_small"
        else:
            result["reason"] = "area_too_large"

    coords = list(poly.exterior.coords)[:-1]
    result["num_vertices"] = len(coords)

    if len(coords) != 4:
        result["is_good"] = False
        if result["reason"] == "good":
            result["reason"] = f"not_quadrilateral_{len(coords)}_vertices"
        if len(coords) < 3:
            return result

    orthogonal_corners = 0
    has_degenerate_edge = False

    for i in range(len(coords)):
        x0, y0 = coords[i - 1][0], coords[i - 1][1]
        x1, y1 = coords[i][0], coords[i][1]
        x2, y2 = coords[(i + 1) % len(coords)][0], coords[(i + 1) % len(coords)][1]
        v1 = (x1 - x0, y1 - y0)
        v2 = (x2 - x1, y2 - y1)
        len1 = (v1[0] ** 2 + v1[1] ** 2) ** 0.5
        len2 = (v2[0] ** 2 + v2[1] ** 2) ** 0.5

        if len1 < 1e-6 or len2 < 1e-6:
            has_degenerate_edge = True
            continue

        dot = v1[0] * v2[0] + v1[1] * v2[1]
        cos_angle = dot / (len1 * len2)
        if abs(cos_angle) <= 0.1:
            orthogonal_corners += 1

    result["right_angles"] = orthogonal_corners

    if has_degenerate_edge:
        result["is_good"] = False
        if result["reason"] == "good":
            result["reason"] = "degenerate_edge"

    if len(coords) == 4 and orthogonal_corners < 4:
        result["is_good"] = False
        if result["reason"] == "good":
            result["reason"] = f"insufficient_right_angles_{orthogonal_corners}"

    return result


def _voronoi_finite_polygons_2d(vor, radius=None):
    """Reconstruct infinite Voronoi regions in a 2D diagram to finite regions.

    Args:
        vor: Voronoi diagram from scipy.spatial.Voronoi
        radius: Optional radius for bounding infinite regions. If None, computed
            from point cloud extent.

    Returns:
        tuple: (regions, vertices) where regions is list of vertex indices per region
            and vertices is array of vertex coordinates.

    Raises:
        ValueError: If input is not 2D.
    """
    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = np.ptp(vor.points, axis=0).max() * 2

    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            new_regions.append(vertices)
            continue

        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                continue

            t = vor.points[p2] - vor.points[p1]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)


def _build_mesh_and_cells(
    start_x: float,
    start_y: float,
    maxx: float,
    maxy: float,
    grid_width: float,
    grid_depth: float,
    prepared_poly,
    polygon_rot: Polygon,
):
    """Build mesh points and Voronoi cells within a rotated polygon.

    Args:
        start_x: Starting x coordinate for mesh generation.
        start_y: Starting y coordinate for mesh generation.
        maxx: Maximum x bound of rotated polygon.
        maxy: Maximum y bound of rotated polygon.
        grid_width: Width of each grid cell.
        grid_depth: Depth of each grid cell.
        prepared_poly: Prepared geometry for efficient spatial queries.
        polygon_rot: Rotated polygon for clipping cells.

    Returns:
        tuple: (cells_rot, quality_rot, mesh_points_rot, good_cells, good_area)
            - cells_rot: List of clipped Voronoi polygons
            - quality_rot: List of quality dictionaries for each cell
            - mesh_points_rot: List of mesh points used to generate Voronoi
            - good_cells: Count of cells meeting quality criteria
            - good_area: Total area of good cells
    """
    mesh_points_rot = []

    y = start_y
    while y <= maxy + grid_depth:
        x = start_x
        while x <= maxx + grid_width:
            p = Point(x, y)
            if prepared_poly.context.oriented_envelope.contains(p):
                mesh_points_rot.append(p)
            x += grid_width
        y += grid_depth

    if len(mesh_points_rot) < 4:
        return [], [], [], 0, 0.0

    pts = np.array([[p.x, p.y] for p in mesh_points_rot])

    x_coords = pts[:, 0]
    y_coords = pts[:, 1]
    if np.allclose(x_coords, x_coords[0]) or np.allclose(y_coords, y_coords[0]):
        return [], [], [], 0, 0.0
    if np.ptp(x_coords) < 1e-6 or np.ptp(y_coords) < 1e-6:
        return [], [], [], 0, 0.0

    vor = Voronoi(pts)
    regions, vertices = _voronoi_finite_polygons_2d(vor)

    cells_rot = []
    quality_rot = []
    good_cells = 0
    good_area = 0.0
    target_area = grid_width * grid_depth

    for region in regions:
        if not region:
            continue
        coords = vertices[region]
        cell_rot = Polygon(coords)

        if cell_rot.is_empty:
            continue
        if not cell_rot.is_valid:
            cell_rot = cell_rot.buffer(0)
            if cell_rot.is_empty:
                continue

        clipped = cell_rot.intersection(polygon_rot)
        if clipped.is_empty or clipped.area <= 0:
            continue

        if clipped.geom_type == "Polygon":
            cells_rot.append(clipped)
            quality_info = is_good_cell(clipped, target_area)
            quality_rot.append(quality_info)
            if quality_info["is_good"]:
                good_cells += 1
                good_area += clipped.area
        elif clipped.geom_type == "MultiPolygon":
            for g in clipped.geoms:
                if g.is_empty or g.area <= 0:
                    continue
                cells_rot.append(g)
                quality_info = is_good_cell(g, target_area)
                quality_rot.append(quality_info)
                if quality_info["is_good"]:
                    good_cells += 1
                    good_area += g.area

    return cells_rot, quality_rot, mesh_points_rot, good_cells, good_area


def grids_from_polygon(
    polygon,
    arterial_line,
    grid_width: float = 100.0,
    grid_depth: float = 100.0,
    dead_end_lines=None,
):
    """Create a mesh-based grid (Voronoi cells from equidistant mesh points).

    The grid is aligned to the arterial line and clipped to the polygon boundary.

    Algorithm:
        1. If an arterial line exists, sample local tangent angles along the line
           from start to end and use them as candidate rotations.
        2. For each candidate angle:
            - Rotate polygon around a fixed origin (midpoint of arterial).
            - Align rows so the arterial sits roughly halfway between rows.
            - Search horizontally (left/right) for the best grid shift.
            - Count "good" cells with area in [95%, 105%] of target area.
            - Accumulate the total area of those good cells.
        3. Choose the angle+shift combination with:
            - The most good cells
            - If tied, the largest total good-cell area

    Args:
        polygon: Input polygon geometry (OGR or Shapely).
        arterial_line: Arterial line geometry for alignment (OGR or Shapely), or None.
        grid_width: Width of each grid cell in coordinate units. Default 100.0.
        grid_depth: Depth of each grid cell in coordinate units. Default 100.0.
        dead_end_lines: Optional shapely LineString/MultiLineString; intersections with
            bad cells are penalized in scoring.

    Returns:
        tuple: (grid_cells, mesh_points, cell_quality)
            - grid_cells: List of Polygon objects in original coordinates
            - mesh_points: List of Point objects in original coordinates
    """
    from rue_lib.streets.cell import Cell

    polygon_shply = feature_geom_to_shapely(polygon)
    mesh_padding = 0

    def _norm_angle_deg(angle: float) -> float:
        return ((angle + 90.0) % 180.0) - 90.0

    candidate_angles_deg: set[float] = set()

    if arterial_line is not None:
        arterial_line_shply = feature_geom_to_shapely(arterial_line)
        length = arterial_line_shply.length

        origin_point = arterial_line_shply.interpolate(length / 2.0)

        if length > 0:
            sample_step = max(grid_width, length / 20.0)
            d = 0.0
            while d < length:
                d2 = min(d + sample_step, length)
                if d2 <= d:
                    break
                p1 = arterial_line_shply.interpolate(d)
                p2 = arterial_line_shply.interpolate(d2)
                dx = p2.x - p1.x
                dy = p2.y - p1.y
                if not (np.isclose(dx, 0.0) and np.isclose(dy, 0.0)):
                    angle = np.degrees(np.arctan2(dy, dx))
                    candidate_angles_deg.add(round(_norm_angle_deg(angle), 3))
                d += sample_step

        if not candidate_angles_deg:
            candidate_angles_deg = {0.0}
    else:
        origin_point = polygon_shply.centroid
        candidate_angles_deg = {0.0}

    best_overall_cells = []
    best_overall_points = []
    best_overall_good = -1
    best_overall_good_area = -1.0
    best_overall_small = math.inf
    best_overall_large = math.inf
    best_overall_dead = math.inf
    best_angle_deg = 0.0

    def _score_tuple(good_cells, dead_cnt, small_cnt, large_cnt, good_area_val):
        return (good_cells, -dead_cnt, -small_cnt, -large_cnt, good_area_val)

    for angle_deg in candidate_angles_deg:
        polygon_rot = rotate(
            polygon_shply,
            -angle_deg,
            origin=(origin_point.x, origin_point.y),
            use_radians=False,
        )
        dead_end_rot = None
        dead_end_prepped = None
        if dead_end_lines is not None:
            dead_end_rot = rotate(
                dead_end_lines,
                -angle_deg,
                origin=(origin_point.x, origin_point.y),
                use_radians=False,
            )
            if not dead_end_rot.is_empty:
                dead_end_prepped = prep(dead_end_rot)

        minx, miny, maxx, maxy = polygon_rot.bounds
        prepared_poly = prep(polygon_rot)

        if arterial_line is not None:
            mid_y = origin_point.y
            dy = mid_y - miny
            remainder_y = dy % grid_depth
            start_y = miny + remainder_y - (grid_depth / 2.0)
        else:
            start_y = miny

        padded_start_x = minx + mesh_padding if minx + mesh_padding < maxx else minx

        if arterial_line is None:
            cells_rot, quality_rot, mesh_points_rot, best_good, best_good_area = (
                _build_mesh_and_cells(
                    padded_start_x,
                    start_y,
                    maxx,
                    maxy,
                    grid_width,
                    grid_depth,
                    prepared_poly,
                    polygon_rot,
                )
            )
            small_cnt = sum(1 for q in quality_rot if q["reason"] == "area_too_small")
            large_cnt = sum(1 for q in quality_rot if q["reason"] == "area_too_large")
            dead_cnt = 0
            if dead_end_prepped:
                for cell_geom, quality in zip(cells_rot, quality_rot):
                    if not quality["is_good"] and dead_end_prepped.intersects(cell_geom):
                        dead_cnt += 1
            best_small = small_cnt
            best_large = large_cnt
            best_dead = dead_cnt
        else:
            best_cells_rot = []
            best_mesh_points_rot = []
            best_good = -1
            best_good_area = -1.0
            best_small = math.inf
            best_large = math.inf
            best_dead = math.inf
            best_score = None

            shift_step = min(10, grid_width)
            shift = -grid_width
            while shift <= grid_width + 1e-6:
                small_cells = []
                large_cells = []

                start_x = padded_start_x + shift
                cells_cand, quality_cand, mesh_points_cand, good_cells, good_area = (
                    _build_mesh_and_cells(
                        start_x,
                        start_y,
                        maxx,
                        maxy,
                        grid_width,
                        grid_depth,
                        prepared_poly,
                        polygon_rot,
                    )
                )

                for quality in quality_cand:
                    if quality["reason"] == "area_too_small":
                        small_cells.append(quality)
                    elif quality["reason"] == "area_too_large":
                        large_cells.append(quality)
                dead_cnt = 0
                if dead_end_prepped:
                    for cell_geom, quality in zip(cells_cand, quality_cand):
                        if not quality["is_good"] and dead_end_prepped.intersects(cell_geom):
                            dead_cnt += 1

                cand_score = _score_tuple(
                    good_cells, dead_cnt, len(small_cells), len(large_cells), good_area
                )
                if (best_score is None) or (cand_score > best_score):
                    best_good = good_cells
                    best_good_area = good_area
                    best_cells_rot = cells_cand
                    best_mesh_points_rot = mesh_points_cand
                    best_small = len(small_cells)
                    best_large = len(large_cells)
                    best_dead = dead_cnt
                    best_score = cand_score

                shift += shift_step

            cells_rot = best_cells_rot
            mesh_points_rot = best_mesh_points_rot
            small_cnt = best_small
            large_cnt = best_large
            dead_cnt = best_dead

        overall_score = _score_tuple(best_good, dead_cnt, small_cnt, large_cnt, best_good_area)
        if (best_overall_points == []) or (
            overall_score
            > _score_tuple(
                best_overall_good,
                best_overall_dead,
                best_overall_small,
                best_overall_large,
                best_overall_good_area,
            )
        ):
            best_overall_good = best_good
            best_overall_good_area = best_good_area
            best_overall_small = small_cnt
            best_overall_large = large_cnt
            best_overall_dead = dead_cnt
            best_overall_cells = cells_rot
            best_overall_points = mesh_points_rot
            best_angle_deg = angle_deg

    if not best_overall_points:
        return [], []

    if best_angle_deg != 0.0:
        grid_cells = [
            rotate(
                c,
                best_angle_deg,
                origin=(origin_point.x, origin_point.y),
                use_radians=False,
            )
            for c in best_overall_cells
        ]
        mesh_points = [
            rotate(
                p,
                best_angle_deg,
                origin=(origin_point.x, origin_point.y),
                use_radians=False,
            )
            for p in best_overall_points
        ]
    else:
        grid_cells = best_overall_cells
        mesh_points = best_overall_points

    final_cells = []
    for i, c in enumerate(grid_cells):
        clipped = c.intersection(polygon_shply)
        if not clipped.is_empty and clipped.area > 0:
            if clipped.geom_type == "Polygon":
                _cell_quality = is_good_cell(clipped, grid_width * grid_depth)
                _cell = Cell(i, clipped, _cell_quality)
                final_cells.append(_cell)

    return final_cells, mesh_points


def grids_from_site(
    output_path: Path,
    site_name: str,
    site_boundary_line_name: str,
    grid_width: float = 100.0,
    grid_depth: float = 100.0,
    grid_layer_name: str | None = None,
    point_layer_name: str | None = None,
    dead_end_lines_layer: str | None = None,
):
    """Generate mesh-based grids for site polygons and save to GeoPackage.

    For each polygon in the site layer, finds its arterial boundary line,
    generates a mesh-based Voronoi grid, and saves both grid cells and
    mesh points into separate layers in the GeoPackage.

    Args:
        output_path: Path to the GeoPackage file.
        site_name: Name of the layer containing site polygons.
        site_boundary_line_name: Name of the layer containing boundary lines.
        grid_width: Width of each grid cell in coordinate units. Default 100.0.
        grid_depth: Depth of each grid cell in coordinate units. Default 100.0.
        grid_layer_name: Optional name for the output grid cells layer.
            If None, defaults to "{site_name}_grid_cells".
        point_layer_name: Optional name for the output mesh points layer.
            If None, defaults to "{site_name}_grid_points".

    Raises:
        RuntimeError: If the GeoPackage cannot be opened or required layers
            are not found.
    """
    ds = ogr.Open(str(output_path), 1)
    if ds is None:
        raise RuntimeError(f"Cannot open {output_path}")

    site_layer = ds.GetLayerByName(site_name)
    site_boundary_lines = ds.GetLayerByName(site_boundary_line_name)
    if site_layer is None or site_boundary_lines is None:
        raise RuntimeError("Could not find site or boundary line layers in GeoPackage")

    srs = site_layer.GetSpatialRef()

    dead_end_lines_geom = None
    if dead_end_lines_layer:
        gdf_dead = gpd.read_file(output_path, layer=dead_end_lines_layer)
        if not gdf_dead.empty:
            dead_end_lines_geom = unary_union(gdf_dead.geometry)

    if grid_layer_name is None:
        grid_layer_name = f"{site_name}_grid_cells"
    if point_layer_name is None:
        point_layer_name = f"{site_name}_grid_points"

    # Recreate layers
    if ds.GetLayerByName(grid_layer_name):
        ds.DeleteLayer(grid_layer_name)
    if ds.GetLayerByName(point_layer_name):
        ds.DeleteLayer(point_layer_name)

    grid_layer = ds.CreateLayer(grid_layer_name, srs, geom_type=ogr.wkbMultiPolygon25D)
    grid_layer.CreateField(ogr.FieldDefn("grid_id", ogr.OFTInteger))

    area_field = ogr.FieldDefn("area", ogr.OFTReal)
    area_field.SetPrecision(2)
    grid_layer.CreateField(area_field)

    grid_layer.CreateField(ogr.FieldDefn("is_good", ogr.OFTInteger))
    grid_layer.CreateField(ogr.FieldDefn("quality", ogr.OFTString))
    grid_layer.CreateField(ogr.FieldDefn("right_angles", ogr.OFTInteger))
    grid_layer.CreateField(ogr.FieldDefn("num_vertices", ogr.OFTInteger))

    area_ratio_field = ogr.FieldDefn("area_ratio", ogr.OFTReal)
    area_ratio_field.SetPrecision(4)
    grid_layer.CreateField(area_ratio_field)

    point_layer = ds.CreateLayer(point_layer_name, srs, geom_type=ogr.wkbPoint25D)
    point_layer.CreateField(ogr.FieldDefn("pt_id", ogr.OFTInteger))

    grid_id = 1
    pt_id = 1

    for polygon in site_layer:
        site_boundary_lines.ResetReading()

        arterial_line = None

        for site_boundary_line in site_boundary_lines:
            if site_boundary_line.GetGeometryRef().Intersects(polygon.GetGeometryRef()):
                if site_boundary_line.GetFieldAsString(2) == "arterial":
                    arterial_line = site_boundary_line
                    break

        print(f"\nDEBUG Generating grid for site polygon ID {polygon.GetFID()}")
        grid_cells, mesh_points = grids_from_polygon(
            polygon,
            arterial_line,
            grid_width,
            grid_depth,
            dead_end_lines=dead_end_lines_geom,
        )

        for _i, _cell in enumerate(grid_cells):
            cell = _cell.geom
            quality_info = _cell.quality
            feat = ogr.Feature(grid_layer.GetLayerDefn())
            feat.SetField("grid_id", grid_id)
            feat.SetField("area", cell.area)

            feat.SetField("is_good", 1 if quality_info["is_good"] else 0)
            feat.SetField("quality", quality_info["reason"])
            feat.SetField("right_angles", quality_info["right_angles"])
            feat.SetField("num_vertices", quality_info["num_vertices"])
            feat.SetField("area_ratio", quality_info["area_ratio"])

            grid_id += 1

            geom = ogr.CreateGeometryFromWkt(cell.wkt)
            # Ensure polygon is converted to multipolygon for compatibility
            if geom.GetGeometryType() == ogr.wkbPolygon:
                multi_geom = ogr.Geometry(ogr.wkbMultiPolygon)
                multi_geom.AddGeometry(geom)
                feat.SetGeometry(multi_geom)
            else:
                feat.SetGeometry(geom)
            grid_layer.CreateFeature(feat)
            feat = None

        for p in mesh_points:
            feat = ogr.Feature(point_layer.GetLayerDefn())
            feat.SetField("pt_id", pt_id)
            pt_id += 1

            geom = ogr.CreateGeometryFromWkt(p.wkt)
            feat.SetGeometry(geom)
            point_layer.CreateFeature(feat)
            feat = None

    ds = None


def extract_grid_lines_in_buffer(
    output_path: Path,
    grid_cells_layer: str,
    buffer_layer: str,
    coord_precision: int = 6,
) -> str | None:
    """Extract internal grid lines (shared edges between cells) that intersect the buffer boundary,
    and create points at the intersection between the grid line and the buffer boundary.

    The point moves along the grid line, so the final geometry is guaranteed to lie on BOTH:
      - the internal grid line
      - the buffer boundary line

    Attributes include:
      - line_id
      - line_length
      - dx, dy (unit direction vector of the grid line)

    Args:
        output_path: Path to GeoPackage
        grid_cells_layer: Name of layer containing grid cells
        buffer_layer: Name of layer containing the buffer zone (polygon)
        coord_precision: Decimal places for rounding coordinates when matching shared edges.

    Returns:
        Name of the extracted points layer, or None if no lines are found.
    """

    print("\nExtracting internal grid lines on buffer boundary...")
    print(f"  Grid layer: {grid_cells_layer}")
    print(f"  Buffer layer: {buffer_layer}")

    # Read the grid cells and buffer zone
    gdf_cells = gpd.read_file(output_path, layer=grid_cells_layer)
    gdf_buffer = gpd.read_file(output_path, layer=buffer_layer)

    print(f"  Loaded {len(gdf_cells)} grid cells")
    print(f"  Loaded {len(gdf_buffer)} buffer feature(s)")

    if gdf_cells.empty:
        print("  Warning: grid_cells_layer is empty")
        return None

    if gdf_buffer.empty:
        print("  Warning: buffer_layer is empty")
        return None

    # Dissolve buffer to a single geometry and take its boundary
    buffer_geom = gdf_buffer.geometry.unary_union
    buffer_boundary = buffer_geom.boundary

    def norm_point(pt):
        return (
            round(pt[0], coord_precision),
            round(pt[1], coord_precision),
        )

    def norm_edge(p1, p2):
        a = norm_point(p1)
        b = norm_point(p2)
        return (a, b) if a <= b else (b, a)

    edge_dict: dict[tuple, dict] = {}
    total_edges = 0

    for _idx, row in gdf_cells.iterrows():
        geom = row.geometry

        if geom is None or geom.is_empty:
            continue

        if not geom.is_valid:
            try:
                geom = geom.buffer(0)
            except TopologicalError:
                continue

        boundary = geom.boundary

        if boundary.geom_type == "LineString":
            lines = [boundary]
        elif boundary.geom_type == "MultiLineString":
            lines = list(boundary.geoms)
        else:
            continue

        for line in lines:
            coords = list(line.coords)
            for i in range(len(coords) - 1):
                p1, p2 = coords[i], coords[i + 1]
                if p1 == p2:
                    continue

                edge_key = norm_edge(p1, p2)
                total_edges += 1

                if edge_key not in edge_dict:
                    edge_dict[edge_key] = {
                        "count": 0,
                        "geom": LineString([p1, p2]),
                    }
                edge_dict[edge_key]["count"] += 1

    print(f"  Processed ~{total_edges} raw edges")
    print(f"  Unique normalized edges: {len(edge_dict)}")

    internal_lines: list[LineString] = []
    shared_edges = 0

    for _edge_key, info in edge_dict.items():
        if info["count"] > 1:
            edge = info["geom"]
            if buffer_boundary.intersects(edge):
                shared_edges += 1
                internal_lines.append(edge)

    print(f"  Shared edges (count > 1) intersecting buffer boundary: {shared_edges}")

    if len(internal_lines) == 0:
        print("  Warning: No internal lines intersect the buffer boundary!")
        return None

    # Helper to pick a single Point from intersection geometry
    def pick_intersection_point(edge: LineString, midpoint: Point):
        inter = edge.intersection(buffer_boundary)
        if inter.is_empty:
            _, snapped = nearest_points(midpoint, buffer_boundary)
            return snapped

        if inter.geom_type == "Point":
            return inter

        if inter.geom_type == "MultiPoint":
            points = list(inter.geoms)
            return min(points, key=lambda p: p.distance(midpoint))

        if inter.geom_type in ("LineString", "MultiLineString"):
            _, snapped = nearest_points(midpoint, inter)
            return snapped

        if inter.geom_type == "GeometryCollection":
            pts = [g for g in inter.geoms if g.geom_type == "Point"]
            if pts:
                return min(pts, key=lambda p: p.distance(midpoint))
            _, snapped = nearest_points(midpoint, inter)
            return snapped

        _, snapped = nearest_points(midpoint, buffer_boundary)
        return snapped

    points = []
    attrs = []

    for line_idx, edge in enumerate(internal_lines):
        midpoint = edge.interpolate(0.5, normalized=True)
        inter_point = pick_intersection_point(edge, midpoint)

        if not isinstance(inter_point, Point) or inter_point.is_empty:
            continue

        coord_start = edge.coords[0]
        coord_end = edge.coords[-1]
        x1, y1 = coord_start[0], coord_start[1]
        x2, y2 = coord_end[0], coord_end[1]
        vx, vy = x2 - x1, y2 - y1
        length = math.hypot(vx, vy)
        if length == 0:
            continue

        dx = vx / length
        dy = vy / length

        points.append(inter_point)
        attrs.append(
            {
                "line_id": line_idx,
                "line_length": float(edge.length),
                "dx": dx,
                "dy": dy,
            }
        )

    print(f"  Created {len(points)} points at grid/buffer intersections")

    if not points:
        print("  Warning: No intersection points could be created!")
        return None

    # Create GeoDataFrame with points as geometry
    gdf_points = gpd.GeoDataFrame(
        attrs,
        geometry=points,
        crs=gdf_cells.crs,
    )

    output_layer_name = f"{grid_cells_layer}_internal_points_in_buffer"
    gdf_points.to_file(output_path, layer=output_layer_name, driver="GPKG")

    print(f"  Saved intersection points to layer: {output_layer_name}")
    return output_layer_name


def generate_local_streets(
    output_path: Path, cfg: StreetConfig, input_layer_name: str
) -> tuple[str, str]:
    """Generate local streets zone from grid blocks.

    Creates a zone for local streets with sidewalks by:
    1. Creating an inner buffer: -(sidewalk_width + road_width/2)
    2. Creating an outer rounded buffer: +sidewalk_width

    Both inner and outer zones are saved as separate layers.

    Args:
        output_path (Path): Path to the GeoPackage containing grid blocks.
        cfg (StreetConfig): Configuration with sidewalk and road width parameters.
        input_layer_name (str): Name of the input grid blocks layer.

    Returns:
        tuple[str, str]: Names of (outer_layer, inner_layer) created.
    """
    output_layer_name = f"{input_layer_name}_local_streets"

    print(f"Creating local streets zone from {input_layer_name}...")

    outer_layer, inner_layer = create_local_streets_zone(
        str(output_path),
        input_layer_name,
        str(output_path),
        output_layer_name,
        cfg.sidewalk_width_m,
        cfg.road_locals_width_m,
    )

    print(f"  Created layers: {outer_layer} (outer), {inner_layer} (inner)")

    return (outer_layer, inner_layer)


def create_cold_boundaries(
    output_path: Path,
    subsites_layer: str,
    off_grid_cells_layer: str,
    arterial_setback_layer: str,
    secondary_setback_layer: str,
    output_layer_name: str,
) -> str:
    """Create cold boundaries by erasing subsites from grid cells and setback areas.

    Cold boundaries represent areas that are not suitable for development, calculated by
    taking the subsite polygons and removing:
    1. Off-grid cells (development cells)
    2. Arterial setback areas
    3. Secondary setback areas

    Args:
        output_path: Path to GeoPackage.
        subsites_layer: Name of layer containing subsite polygons.
        off_grid_cells_layer: Name of layer containing off-grid cells
            (e.g., "14_off_grid_cells_fixed_by_perp_lines_no_dead_ends").
        arterial_setback_layer: Name of layer containing arterial setback polygons
            (e.g., "10a_arterial_setback_clipped").
        secondary_setback_layer: Name of layer containing secondary setback polygons
            (e.g., "10b_secondary_setback_clipped").
        output_layer_name: Name of layer to save cleaned cells to.

    Returns:
        Name of the created cold boundaries layer.
    """
    print("\nCreating cold boundaries...")

    gdf_subsites = gpd.read_file(output_path, layer=subsites_layer)
    gdf_off_grid = gpd.read_file(output_path, layer=off_grid_cells_layer)
    gdf_arterial_setback = gpd.read_file(output_path, layer=arterial_setback_layer)
    gdf_secondary_setback = gpd.read_file(output_path, layer=secondary_setback_layer)

    print(f"  Loaded {len(gdf_subsites)} subsite(s)")
    print(f"  Loaded {len(gdf_off_grid)} off-grid cells")
    print(f"  Loaded {len(gdf_arterial_setback)} arterial setback polygon(s)")
    print(f"  Loaded {len(gdf_secondary_setback)} secondary setback polygon(s)")

    if gdf_subsites.empty:
        print("  Warning: No subsites found")
        return None

    subsites_union = gdf_subsites.geometry.unary_union
    print("  Merged subsites into unified geometry")

    layers_to_erase = []
    gap_fill_buffer = 1.0

    if not gdf_off_grid.empty:
        layers_to_erase.append(gdf_off_grid.geometry.unary_union)
        print("  Added off-grid cells to erase list")

    if not gdf_arterial_setback.empty:
        layers_to_erase.append(gdf_arterial_setback.geometry.unary_union)
        print("  Added arterial setback to erase list")

    if not gdf_secondary_setback.empty:
        layers_to_erase.append(gdf_secondary_setback.geometry.unary_union)
        print("  Added secondary setback to erase list")

    if not layers_to_erase:
        print("  Warning: No layers to erase, returning subsites as cold boundaries")
        cold_result = subsites_union
    else:
        erase_union = unary_union(layers_to_erase)
        print(f"  Merged {len(layers_to_erase)} layer(s) to erase")

        erase_buffered = erase_union.buffer(gap_fill_buffer, join_style=2, cap_style=2)
        erase_final = erase_buffered.buffer(-gap_fill_buffer, join_style=2, cap_style=2)
        print(f"  Applied gap closure (buffer={gap_fill_buffer})")

        cold_result = subsites_union.difference(erase_final)
        print("  Erased development areas from subsites")

    # Handle the result geometry
    if cold_result.is_empty:
        print("  Warning: No cold boundaries remain after erase operation")
        return None

    if cold_result.geom_type == "Polygon":
        cold_polygons = [cold_result]
    elif cold_result.geom_type == "MultiPolygon":
        cold_polygons = list(cold_result.geoms)
    else:
        print(f"  Warning: Unexpected geometry type: {cold_result.geom_type}")
        return None

    print(f"  Created {len(cold_polygons)} initial cold boundary polygon(s)")

    cold_boundaries = []
    cold_attrs = []

    for poly_id, poly in enumerate(cold_polygons):
        if poly.is_empty or poly.area < 0.1:
            continue

        cold_boundaries.append(poly)
        cold_attrs.append(
            {
                "id": poly_id,
                "cold_boundary": True,
                "area": poly.area,
            }
        )

    if not cold_boundaries:
        print("  Warning: No valid cold boundaries after filtering")
        return None

    gdf_cold = gpd.GeoDataFrame(cold_attrs, geometry=cold_boundaries, crs=gdf_subsites.crs)

    if not output_layer_name:
        output_layer_name = "cold_boundaries"
    gdf_cold.to_file(output_path, layer=output_layer_name, driver="GPKG")

    print(f"  Saved {len(cold_boundaries)} cold boundary polygon(s) to layer: {output_layer_name}")
    return output_layer_name
