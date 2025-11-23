from pathlib import Path

import numpy as np
from osgeo import ogr
from scipy.spatial import Voronoi
from shapely.affinity import rotate
from shapely.geometry import Point, Polygon
from shapely.prepared import prep

from rue_lib.core.helpers import feature_geom_to_shapely


def _voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite Voronoi regions in a 2D diagram to finite regions.
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


def grids_from_polygon(
    polygon,
    arterial_line,
    grid_width: float = 100.0,
    grid_depth: float = 100.0,
):
    """Create a mesh-based grid (Voronoi cells from equidistant mesh points)
    aligned to the arterial line and clipped to ``polygon``.

    Algorithm:
    - If an arterial line exists, sample local tangent angles along the line
      from start to end and use them as candidate rotations.
    - For each candidate angle:
        * Rotate polygon around a fixed origin (midpoint of arterial).
        * Align rows so the arterial sits roughly halfway between rows.
        * Search horizontally (left/right) for the best grid shift.
        * Count "good" cells with area in [width*depth, 1.10 * width*depth).
          Also accumulate the total area of those good cells.
    - Choose the angle+shift combination with:
        1. the most good cells;
        2. if tied, the largest total good-cell area
           (equivalent to highest good-area / polygon-area ratio).

    Returns
    -------
    (grid_cells, mesh_points)
        grid_cells: list[Polygon] in original coordinates
        mesh_points: list[Point] in original coordinates
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

        def build_mesh_and_cells(
            start_x: float,
            *,
            start_y=start_y,
            maxx=maxx,
            maxy=maxy,
            prepared_poly=prepared_poly,
            polygon_rot=polygon_rot,
        ):
            mesh_points_rot = []

            y = start_y
            while y <= maxy + grid_depth:
                x = start_x
                while x <= maxx + grid_width:
                    p = Point(x, y)
                    if prepared_poly.contains(p):
                        mesh_points_rot.append(p)
                    x += grid_width
                y += grid_depth

            if len(mesh_points_rot) < 4:
                return [], [], 0, 0.0

            pts = np.array([[p.x, p.y] for p in mesh_points_rot])

            x_coords = pts[:, 0]
            y_coords = pts[:, 1]
            if np.allclose(x_coords, x_coords[0]) or np.allclose(y_coords, y_coords[0]):
                return [], [], 0, 0.0
            if np.ptp(x_coords) < 1e-6 or np.ptp(y_coords) < 1e-6:
                return [], [], 0, 0.0

            vor = Voronoi(pts)
            regions, vertices = _voronoi_finite_polygons_2d(vor)

            cells_rot = []
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
                    area = clipped.area
                    cells_rot.append(clipped)
                    if target_area <= area < target_area * 1.10:
                        good_cells += 1
                        good_area += area
                elif clipped.geom_type == "MultiPolygon":
                    for g in clipped.geoms:
                        if g.is_empty or g.area <= 0:
                            continue
                        area = g.area
                        cells_rot.append(g)
                        if target_area <= area < target_area * 1.10:
                            good_cells += 1
                            good_area += area

            return cells_rot, mesh_points_rot, good_cells, good_area

        if arterial_line is None:
            cells_rot, mesh_points_rot, best_good, best_good_area = build_mesh_and_cells(minx)
        else:
            best_cells_rot = []
            best_mesh_points_rot = []
            best_good = -1
            best_good_area = -1.0

            shift_step = min(10.0, grid_width)
            shift = -grid_width
            while shift <= grid_width + 1e-6:
                start_x = minx + shift
                cells_cand, mesh_points_cand, good_cells, good_area = build_mesh_and_cells(start_x)

                if good_cells > best_good or (
                    good_cells == best_good and good_area > best_good_area
                ):
                    best_good = good_cells
                    best_good_area = good_area
                    best_cells_rot = cells_cand
                    best_mesh_points_rot = mesh_points_cand

                shift += shift_step

            cells_rot = best_cells_rot
            mesh_points_rot = best_mesh_points_rot

        if best_good > best_overall_good or (
            best_good == best_overall_good and best_good_area > best_overall_good_area
        ):
            best_overall_good = best_good
            best_overall_good_area = best_good_area
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
    for c in grid_cells:
        clipped = c.intersection(polygon_shply)
        if not clipped.is_empty and clipped.area > 0:
            if clipped.geom_type == "Polygon":
                final_cells.append(clipped)
            elif clipped.geom_type == "MultiPolygon":
                for g in clipped.geoms:
                    if not g.is_empty and g.area > 0:
                        final_cells.append(g)

    return final_cells, mesh_points


def grids_from_site(
    output_path: Path,
    site_name: str,
    site_boundary_line_name: str,
    grid_width: float = 100.0,
    grid_depth: float = 100.0,
    grid_layer_name: str | None = None,
    point_layer_name: str | None = None,
):
    """
    For each polygon in `site_name`, find its arterial boundary line,
    generate mesh-based grid, and save both grid cells and mesh points
    into layers in the same GeoPackage (`output_path`).
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

        grid_cells, mesh_points = grids_from_polygon(polygon, arterial_line, grid_width, grid_depth)

        for cell in grid_cells:
            feat = ogr.Feature(grid_layer.GetLayerDefn())
            feat.SetField("grid_id", grid_id)
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
