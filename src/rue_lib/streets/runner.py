# src/rue_lib/streets/runner.py
from pathlib import Path
from typing import Optional

from osgeo import gdal, ogr

from rue_lib.core.geometry import buffer_layer, get_utm_zone_from_layer, reproject_layer
from rue_lib.streets.grids import grids_from_site, is_good_cell

from .blocks_orthogonal import (
    clip_site_by_roads,
)
from .cell import Cell, fix_grid_cells_with_perpendicular_lines
from .config import StreetConfig
from .fix_bad_angles import inspect_and_fix_bad_angles
from .operations import (
    break_multipart_features,
    cleanup_grid_blocks,
    clip_layer,
    create_grid_from_on_grid,
    create_local_streets_zone,
    erase_layer,
    export_layer_to_geojson,
    extract_by_expression,
    extract_site_boundary_lines,
    merge_grid_layers_with_type,
)

gdal.UseExceptions()


def generate_on_grid_blocks(output_path: Path, cfg: StreetConfig) -> Path:
    print("Clip arterial setback to site boundary...")
    clip_layer(
        output_path,
        "arterial_setback",
        output_path,
        "site_clipped_by_roads",
        output_path,
        "arterial_setback_clipped",
    )

    print("Clip secondary setback to site boundary...")
    clip_layer(
        output_path,
        "secondary_setback",
        output_path,
        "site_clipped_by_roads",
        output_path,
        "secondary_setback_clipped",
    )

    print("Intersect arterial and secondary setbacks...")
    clip_layer(
        output_path,
        "arterial_setback_clipped",
        output_path,
        "secondary_setback_clipped",
        output_path,
        "intersected_setbacks",
    )

    print("Arterial setback without intersection...")
    erase_layer(
        output_path,
        "arterial_setback_clipped",
        output_path,
        "intersected_setbacks",
        output_path,
        "arterial_setback_final",
    )

    print("Secondary setback without intersection...")
    erase_layer(
        output_path,
        "secondary_setback_clipped",
        output_path,
        "intersected_setbacks",
        output_path,
        "secondary_setback_final",
    )

    print("Breaking arterial setback multipart features...")
    break_multipart_features(
        output_path,
        "arterial_setback_final",
        output_path,
        "arterial_setback_final",
    )

    print("Breaking secondary setback multipart features...")
    break_multipart_features(
        output_path,
        "secondary_setback_final",
        output_path,
        "secondary_setback_final",
    )

    print("Create grid from on-grid arterial setback...")
    create_grid_from_on_grid(
        output_path,
        "arterial_setback_final",
        "arterial_roads",
        cfg.off_grid_partitions_preferred_width,
        output_path,
        "arterial_setback_grid",
        road_buffer_distance=cfg.road_arterial_width_m,
    )

    print("Clean up grid blocks...")
    cleanup_grid_blocks(
        output_path,
        "arterial_setback_grid",
        output_path,
        "arterial_setback_grid_cleaned",
        cfg.arterial_setback_depth * cfg.off_grid_partitions_preferred_width * 0.5,
    )

    print("Create grid from on-grid secondary setback...")
    create_grid_from_on_grid(
        output_path,
        "secondary_setback_final",
        "secondary_roads",
        cfg.off_grid_partitions_preferred_width,
        output_path,
        "secondary_setback_grid",
        intersected_setbacks_layer_name="intersected_setbacks",
        road_buffer_distance=cfg.road_secondary_width_m,
    )

    print("Clean up grid blocks...")
    cleanup_grid_blocks(
        output_path,
        "secondary_setback_grid",
        output_path,
        "secondary_setback_grid_cleaned",
        cfg.secondary_setback_depth * cfg.off_grid_partitions_preferred_width * 0.5,
    )

    return output_path


def fix_grid_cells(
    output_path: Path,
    layer_name: str,
    target_area: float,
    cluster_width: float = 20.0,
) -> str:
    """Fix bad grid cells by adjusting angles and splitting oversized cells.

    Args:
        output_path: Path to GeoPackage
        layer_name: Name of layer containing grid cells
        target_area: Target area for cells in square meters
        cluster_width: Width for splitting operations in meters

    Returns:
        Name of the fixed layer
    """
    import geopandas as gpd

    print(f"\nFixing bad cells in layer: {layer_name}")

    # Read grid cells
    gdf = gpd.read_file(output_path, layer=layer_name)
    print(f"  Loaded {len(gdf)} cells")

    # Convert to Cell objects
    cells = []
    for idx, row in gdf.iterrows():
        quality = is_good_cell(row.geometry, target_area)
        cell = Cell(id=idx, geom=row.geometry, quality=quality)
        cells.append(cell)

    # Count initial bad cells
    initial_bad = sum(1 for c in cells if not c.quality.get("is_good", False))
    print(f"  Initial bad cells: {initial_bad}")

    # Fix bad angles
    cells = inspect_and_fix_bad_angles(cells, target_area, cluster_width=cluster_width)

    # Count final bad cells
    final_bad = sum(1 for c in cells if not c.quality.get("is_good", False))
    improvement = initial_bad - final_bad

    print(f"\n  Final bad cells: {final_bad}")
    if improvement > 0:
        print(f"  ✓ Improved {improvement} cells")

    # Convert back to GeoDataFrame
    data = []
    for cell in cells:
        data.append(
            {
                "id": cell.id,
                "geometry": cell.geom,
                "area": cell.geom.area,
                "is_good": cell.quality.get("is_good", False),
                "quality": cell.quality.get("reason", ""),
                "right_angles": cell.quality.get("right_angles", 0),
                "num_vertices": cell.quality.get("num_vertices", 0),
                "area_ratio": cell.quality.get("area_ratio", 0.0),
            }
        )

    gdf_fixed = gpd.GeoDataFrame(data, crs=gdf.crs)

    # Save to new layer
    fixed_layer_name = f"fixed_{layer_name}"
    gdf_fixed.to_file(output_path, layer=fixed_layer_name, driver="GPKG")
    print(f"  Saved fixed cells to layer: {fixed_layer_name}")

    return fixed_layer_name


def apply_inner_buffer_to_cells(
    output_path: Path,
    site_boundary_layer: str,
    buffer_distance: float,
) -> str:
    """Apply inner buffer to site boundary lines and clip to site boundary.

    Args:
        output_path: Path to GeoPackage
        site_boundary_layer: Name of layer containing site boundary (polygon)
        buffer_distance: Buffer distance in meters (negative for inner buffer)

    Returns:
        Name of the buffered layer
    """
    print("\nCreating inner buffer zone from site boundary lines...")
    print(f"  Buffer distance: {buffer_distance}m")

    # Step 1: Buffer the site_boundary_lines
    buffer_layer(
        str(output_path),
        "site_boundary_lines",
        buffer_distance,  # Use positive value, will invert direction below
        str(output_path),
        "site_boundary_buffer_temp",
        dissolve=True,
    )

    # Step 2: Clip the buffer to the site boundary
    print("  Clipping buffer to site boundary...")
    clip_layer(
        output_path,
        "site_boundary_buffer_temp",
        output_path,
        site_boundary_layer,
        output_path,
        "site_boundary_inner_buffer",
    )

    print("  Created buffered layer: site_boundary_inner_buffer")
    return "site_boundary_inner_buffer"


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
    import math

    import geopandas as gpd
    from shapely.errors import TopologicalError
    from shapely.geometry import LineString, Point
    from shapely.ops import nearest_points

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
    buffer_boundary = buffer_geom.boundary  # lines along which we want the intersections

    # Helpers to normalize points/edges with rounding
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

    # Collect only shared edges that intersect the BUFFER BOUNDARY (not just buffer polygon)
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
            # Fallback: project midpoint onto boundary
            _, snapped = nearest_points(midpoint, buffer_boundary)
            return snapped

        if inter.geom_type == "Point":
            return inter

        if inter.geom_type == "MultiPoint":
            points = list(inter.geoms)
            return min(points, key=lambda p: p.distance(midpoint))

        if inter.geom_type in ("LineString", "MultiLineString"):
            # Overlap: choose nearest point on overlap to the midpoint
            _, snapped = nearest_points(midpoint, inter)
            return snapped

        if inter.geom_type == "GeometryCollection":
            pts = [g for g in inter.geoms if g.geom_type == "Point"]
            if pts:
                return min(pts, key=lambda p: p.distance(midpoint))
            # fallback
            _, snapped = nearest_points(midpoint, inter)
            return snapped

        # Last resort
        _, snapped = nearest_points(midpoint, buffer_boundary)
        return snapped

    # Build points (intersection) + attributes, including line direction dx, dy
    points = []
    attrs = []

    for line_idx, edge in enumerate(internal_lines):
        # Midpoint purely used to pick the closest intersection when multiple
        midpoint = edge.interpolate(0.5, normalized=True)
        inter_point = pick_intersection_point(edge, midpoint)

        if not isinstance(inter_point, Point) or inter_point.is_empty:
            continue

        # Direction of the grid line (from first to last coordinate)
        x1, y1, z1 = edge.coords[0]
        x2, y2, z2 = edge.coords[-1]
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


def create_perpendicular_lines_inside_buffer_from_points(
    output_path: Path,
    points_layer: str,
    buffer_layer: str,
    line_length: float,
    output_layer_name: Optional[str] = None,
    tangent_step_fraction: float = 0.001,
) -> Optional[str]:
    """Create perpendicular lines inside the buffer, starting from points on the buffer boundary.

    The perpendicular is computed with respect to the buffer boundary:
      - For each point, compute the local tangent of the buffer boundary.
      - From the tangent, derive a normal (perpendicular) direction.
      - Choose the normal direction that goes INSIDE the buffer polygon.
      - Build a line segment from the point inward, clipped to the buffer polygon.

    Args:
        output_path: Path to the GeoPackage.
        points_layer: Name of the layer containing points on the buffer boundary.
        buffer_layer: Name of the layer containing the buffer polygon.
        line_length: Desired length of the perpendicular line (map units, e.g. meters).
        output_layer_name: Optional name for the output layer. If None, a default is used.
        tangent_step_fraction: Step along the boundary (fraction of boundary length)
                               used to approximate the local tangent.

    Returns:
        Name of the created perpendicular lines layer, or None on failure.
    """
    import math

    import geopandas as gpd
    from shapely.errors import TopologicalError
    from shapely.geometry import LineString, Point

    print(f"\nCreating perpendicular lines inside buffer from points in layer: {points_layer}")

    # Load data
    gdf_pts = gpd.read_file(output_path, layer=points_layer)
    gdf_buffer = gpd.read_file(output_path, layer=buffer_layer)

    print(f"  Loaded {len(gdf_pts)} points")
    print(f"  Loaded {len(gdf_buffer)} buffer feature(s)")

    if gdf_pts.empty:
        print("  Warning: points_layer is empty")
        return None

    if gdf_buffer.empty:
        print("  Warning: buffer_layer is empty")
        return None

    # Dissolve buffer to single polygon/multipolygon and get its boundary
    buffer_geom = gdf_buffer.geometry.unary_union
    if buffer_geom is None or buffer_geom.is_empty:
        print("  Warning: buffer geometry is empty")
        return None

    try:
        buffer_boundary = buffer_geom.boundary
    except TopologicalError:
        print("  Warning: buffer geometry is invalid, attempting buffer(0) fix")
        buffer_geom = buffer_geom.buffer(0)
        buffer_boundary = buffer_geom.boundary

    boundary_len = buffer_boundary.length
    if boundary_len == 0:
        print("  Warning: buffer boundary length is zero")
        return None

    # Step along the boundary to approximate the tangent
    step = boundary_len * tangent_step_fraction
    if step == 0:
        step = line_length * 0.01  # small fallback

    if output_layer_name is None:
        output_layer_name = f"{points_layer}_perp_inside_buffer"

    half_eps_in = line_length * 0.01  # small step to test inside/outside

    lines = []
    records = []

    for _idx, row in gdf_pts.iterrows():
        pt = row.geometry
        if pt is None or pt.is_empty:
            continue

        if not isinstance(pt, Point):
            continue

        # Use only 2D coords (avoid mixing 2D/3D)
        x0, y0 = pt.x, pt.y

        # Position along buffer boundary (2D)
        s = buffer_boundary.project(Point(x0, y0))

        # Approximate local tangent by sampling slightly before and after s
        s0 = max(0.0, s - step)
        s1 = min(boundary_len, s + step)

        p0 = buffer_boundary.interpolate(s0)
        p1 = buffer_boundary.interpolate(s1)

        tx = p1.x - p0.x
        ty = p1.y - p0.y

        # If tangent is degenerate, skip
        t_norm = math.hypot(tx, ty)
        if t_norm == 0:
            continue

        tx /= t_norm
        ty /= t_norm

        # Two possible normals (perpendicular to tangent)
        n1x, n1y = -ty, tx
        n2x, n2y = ty, -tx

        # Choose the normal that goes INSIDE the buffer
        test1 = Point(x0 + n1x * half_eps_in, y0 + n1y * half_eps_in)
        test2 = Point(x0 + n2x * half_eps_in, y0 + n2y * half_eps_in)

        inside1 = buffer_geom.contains(test1)
        inside2 = buffer_geom.contains(test2)

        if inside1 and not inside2:
            nx, ny = n1x, n1y
        elif inside2 and not inside1:
            nx, ny = n2x, n2y
        elif inside1 and inside2:
            # Both directions considered inside (e.g. very thick buffer); pick one
            nx, ny = n1x, n1y
        else:
            # Neither direction goes inside (point may not be exactly on boundary)
            continue

        # Build candidate line from boundary point inward **in 2D**
        end_x = x0 + nx * line_length
        end_y = y0 + ny * line_length

        # <<< FIXED PART: use coordinate tuples, not Point geometries >>>
        candidate_line = LineString([(x0, y0), (end_x, end_y)])

        # Clip line to buffer polygon so we stay inside
        clipped = candidate_line.intersection(buffer_geom)
        if clipped.is_empty:
            continue

        if clipped.geom_type == "LineString":
            line_geom = clipped
        elif clipped.geom_type == "MultiLineString":
            # Choose the segment that starts closest to the original point
            segments = list(clipped.geoms)
            line_geom = min(segments, key=lambda seg: seg.distance(Point(x0, y0)))
        else:
            continue

        if line_geom.length == 0:
            continue

        lines.append(line_geom)
        # Copy attributes except geometry
        rec = row.drop(labels="geometry").to_dict()
        records.append(rec)

    if not lines:
        print("  Warning: no perpendicular lines created")
        return None

    import geopandas as gpd

    gdf_out = gpd.GeoDataFrame(records, geometry=lines, crs=gdf_pts.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")

    print(f"  Saved perpendicular lines inside buffer to layer: {output_layer_name}")
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


def generate_streets(cfg: StreetConfig) -> Path:
    """
    Generate street blocks from roads and parcels (version 3)

    Returns:
        Path to output blocks file
    """
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    output_gpkg = output_dir / "outputs.gpkg"
    output_path = str(output_gpkg)

    print("Step 1: Determining UTM zone...")
    site_ds = ogr.Open(cfg.parcel_path)
    site_layer = site_ds.GetLayer()
    utm_epsg = get_utm_zone_from_layer(site_layer)
    site_ds = None
    print(f"Using UTM EPSG: {utm_epsg}")

    print("Step 2: Reprojecting layers to UTM...")
    site_layer_name = reproject_layer(cfg.parcel_path, output_path, utm_epsg)
    roads_layer_name = reproject_layer(cfg.roads_path, output_path, utm_epsg)

    print("Step 3: Creating site polygon clipped by roads...")
    clipped_site_layer = clip_site_by_roads(
        output_path,
        site_layer_name,
        roads_layer_name,
        "site_clipped_by_roads",
        cfg.road_arterial_width_m,
        cfg.road_secondary_width_m,
        cfg.road_local_width_m,
    )

    print("Step 4: Extracting arterial roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_art'", output_path, "arterial_roads"
    )

    print("Step 5: Extracting secondary roads...")
    extract_by_expression(
        output_path, roads_layer_name, "road_type = 'road_sec'", output_path, "secondary_roads"
    )

    print("Step 6: Creating arterial road setback zone...")
    buffer_layer(
        output_path,
        "arterial_roads",
        cfg.arterial_setback_depth,
        output_path,
        "arterial_setback",
        dissolve=True,
    )

    print("Step 7: Creating secondary road setback zone...")
    buffer_layer(
        output_path,
        "secondary_roads",
        cfg.secondary_setback_depth,
        output_path,
        "secondary_setback",
        dissolve=True,
    )

    print("Step 8: Removing arterial setback from site...")
    erase_layer(
        output_path,
        clipped_site_layer,
        output_path,
        "arterial_setback",
        output_path,
        "site_minus_arterial_setback",
    )

    print("Step 9: Removing secondary setback from site...")
    erase_layer(
        output_path,
        "site_minus_arterial_setback",
        output_path,
        "secondary_setback",
        output_path,
        "site_minus_all_setbacks",
    )

    print("Step 10: Generating on-grid blocks...")
    generate_on_grid_blocks(output_gpkg, cfg)

    print("Step 11: Generating local streets zones for arterial setback...")
    generate_local_streets(
        output_gpkg,
        cfg,
        "arterial_setback_grid_cleaned",
    )

    print("Step 12: Generating local streets zones for secondary setback...")
    generate_local_streets(
        output_gpkg,
        cfg,
        "secondary_setback_grid_cleaned",
    )

    print("Step 13: Extracting site boundary lines")
    extract_site_boundary_lines(
        output_gpkg,
        "site_minus_all_setbacks",
        "arterial_setback_final",
        "secondary_setback_final",
        output_gpkg,
        "site_boundary_lines",
    )

    print("Step 14: Generating off-grid blocks...")
    grids_from_site(
        output_gpkg,
        "site_minus_all_setbacks",
        "site_boundary_lines",
        cfg.off_grid_partitions_preferred_width,
        cfg.off_grid_partitions_preferred_depth,
    )

    print("Step 14a: Creating inner buffer zone from site boundary...")
    _buffered_layer = apply_inner_buffer_to_cells(
        output_gpkg,
        "site_minus_all_setbacks",
        cfg.part_loc_d,
    )

    print("Step 14b: Extracting grid lines inside buffer zone...")
    lines_layer = extract_grid_lines_in_buffer(
        output_gpkg,
        "site_minus_all_setbacks_grid_cells",
        "site_boundary_inner_buffer",
    )
    if lines_layer is not None:
        perp_inside_layer = create_perpendicular_lines_inside_buffer_from_points(
            output_gpkg,
            lines_layer,
            "site_boundary_inner_buffer",
            line_length=cfg.part_loc_d * 2,
        )

        if perp_inside_layer is not None:
            print("Step 14c: Fixing grid cells with perpendicular lines...")
            _fixed_cells_layer = fix_grid_cells_with_perpendicular_lines(
                output_gpkg,
                "site_minus_all_setbacks_grid_cells",
                perp_inside_layer,
                "site_boundary_inner_buffer",
                target_area=cfg.off_grid_partitions_preferred_width
                * cfg.off_grid_partitions_preferred_depth,
            )

    print("Step 17: Merging all grid layers with grid_type information...")
    merge_grid_layers_with_type(
        str(output_gpkg),
        str(output_gpkg),
        "all_grids_merged",
        [
            ("intersected_setbacks", "on_grid_intersected"),
            ("arterial_setback_grid_cleaned", "on_grid_art"),
            ("secondary_setback_grid_cleaned", "on_grid_sec"),
        ],
    )

    print("Step 18: Exporting merged grids to GeoJSON...")
    output_geojson = output_dir / "all_grids_merged.geojson"
    export_layer_to_geojson(
        str(output_gpkg),
        "all_grids_merged",
        str(output_geojson),
    )

    print(f"\nProcessing complete! Output saved to: {output_gpkg}")
    print("\nFinal layers:")
    print(f"  - {clipped_site_layer}: Site polygon with roads subtracted")
    print("  - all_grids_merged: Merged grid cells with grid_type classification")
    print("\nGeoJSON export:")
    # print(f"  - {output_geojson}: Merged grids with grid_type classification")

    return output_gpkg
