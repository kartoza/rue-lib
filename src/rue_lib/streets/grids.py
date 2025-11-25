from pathlib import Path

import numpy as np
from osgeo import ogr
from scipy.spatial import Voronoi
from shapely.affinity import rotate
from shapely.geometry import Point, Polygon
from shapely.prepared import prep

from rue_lib.core.helpers import feature_geom_to_shapely


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
    min_area = target_area * 0.95

    if not (min_area <= area < max_area):
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

    Returns:
        tuple: (grid_cells, mesh_points, cell_quality)
            - grid_cells: List of Polygon objects in original coordinates
            - mesh_points: List of Point objects in original coordinates
            - cell_quality: List of quality dictionaries for each cell
    """
    polygon_shply = feature_geom_to_shapely(polygon)

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
    best_overall_quality = []
    best_overall_points = []
    best_overall_good = -1
    best_overall_good_area = -1.0
    best_angle_deg = 0.0

    for angle_deg in candidate_angles_deg:
        polygon_rot = rotate(
            polygon_shply,
            -angle_deg,
            origin=(origin_point.x, origin_point.y),
            use_radians=False,
        )

        minx, miny, maxx, maxy = polygon_rot.bounds
        prepared_poly = prep(polygon_rot)

        if arterial_line is not None:
            mid_y = origin_point.y
            dy = mid_y - miny
            remainder_y = dy % grid_depth
            start_y = miny + remainder_y - (grid_depth / 2.0)
        else:
            start_y = miny

        if arterial_line is None:
            cells_rot, quality_rot, mesh_points_rot, best_good, best_good_area = (
                _build_mesh_and_cells(
                    minx,
                    start_y,
                    maxx,
                    maxy,
                    grid_width,
                    grid_depth,
                    prepared_poly,
                    polygon_rot,
                )
            )
        else:
            best_cells_rot = []
            best_quality_rot = []
            best_mesh_points_rot = []
            best_good = -1
            best_good_area = -1.0

            shift_step = min(10.0, grid_width)
            shift = -grid_width
            while shift <= grid_width + 1e-6:
                start_x = minx + shift
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

                if good_cells > best_good or (
                    good_cells == best_good and good_area > best_good_area
                ):
                    best_good = good_cells
                    best_good_area = good_area
                    best_cells_rot = cells_cand
                    best_quality_rot = quality_cand
                    best_mesh_points_rot = mesh_points_cand

                shift += shift_step

            cells_rot = best_cells_rot
            quality_rot = best_quality_rot
            mesh_points_rot = best_mesh_points_rot

        if best_good > best_overall_good or (
            best_good == best_overall_good and best_good_area > best_overall_good_area
        ):
            best_overall_good = best_good
            best_overall_good_area = best_good_area
            best_overall_cells = cells_rot
            best_overall_quality = quality_rot
            best_overall_points = mesh_points_rot
            best_angle_deg = angle_deg

    if not best_overall_points:
        return [], [], []

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

    cell_quality = best_overall_quality

    final_cells = []
    final_quality = []
    for i, c in enumerate(grid_cells):
        clipped = c.intersection(polygon_shply)
        if not clipped.is_empty and clipped.area > 0:
            if clipped.geom_type == "Polygon":
                final_cells.append(clipped)
                final_quality.append(cell_quality[i])
            elif clipped.geom_type == "MultiPolygon":
                for g in clipped.geoms:
                    if not g.is_empty and g.area > 0:
                        final_cells.append(g)
                        final_quality.append(cell_quality[i])

    return final_cells, mesh_points, final_quality


def grids_from_site(
    output_path: Path,
    site_name: str,
    site_boundary_line_name: str,
    grid_width: float = 100.0,
    grid_depth: float = 100.0,
    grid_layer_name: str | None = None,
    point_layer_name: str | None = None,
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

    if grid_layer_name is None:
        grid_layer_name = f"{site_name}_grid_cells"
    if point_layer_name is None:
        point_layer_name = f"{site_name}_grid_points"

    # Recreate layers
    if ds.GetLayerByName(grid_layer_name):
        ds.DeleteLayer(grid_layer_name)
    if ds.GetLayerByName(point_layer_name):
        ds.DeleteLayer(point_layer_name)

    grid_layer = ds.CreateLayer(grid_layer_name, srs, geom_type=ogr.wkbPolygon)
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

    point_layer = ds.CreateLayer(point_layer_name, srs, geom_type=ogr.wkbPoint)
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

        grid_cells, mesh_points, cell_quality = grids_from_polygon(
            polygon, arterial_line, grid_width, grid_depth
        )

        for i, cell in enumerate(grid_cells):
            feat = ogr.Feature(grid_layer.GetLayerDefn())
            feat.SetField("grid_id", grid_id)
            feat.SetField("area", cell.area)

            quality_info = cell_quality[i]
            feat.SetField("is_good", 1 if quality_info["is_good"] else 0)
            feat.SetField("quality", quality_info["reason"])
            feat.SetField("right_angles", quality_info["right_angles"])
            feat.SetField("num_vertices", quality_info["num_vertices"])
            feat.SetField("area_ratio", quality_info["area_ratio"])

            grid_id += 1

            geom = ogr.CreateGeometryFromWkt(cell.wkt)
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
