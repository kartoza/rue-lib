import math

from osgeo import ogr
from shapely import ops as shapely_ops

from rue_lib.core.helpers import create_or_replace_layer


def extract_off_grid_adjacent_lines(
    input_gpkg: str,
    off_grid_layer_name: str,
    on_grid_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> str:
    """Extract off-grid boundary line segments that touch on-grid blocks.

    Args:
        input_gpkg: Path to the GeoPackage containing the off/on grid layers.
        off_grid_layer_name: Layer name for off-grid polygons.
        on_grid_layer_name: Layer name for on-grid polygons.
        output_gpkg: Path to the GeoPackage where the output is written.
        output_layer_name: Name for the output lines layer.
        split_output_layer_name: Optional name for the split lines layer. If provided,
            lines are further split at bends over the angle threshold.
        angle_threshold_deg: Minimum turn angle (degrees) to trigger a split.

    Returns:
        The name of the created lines layer.
    """
    ds = ogr.Open(input_gpkg, 0)
    if ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    off_layer = ds.GetLayerByName(off_grid_layer_name)
    if off_layer is None:
        raise ValueError(f"Layer {off_grid_layer_name} not found")

    on_layer = ds.GetLayerByName(on_grid_layer_name)
    if on_layer is None:
        raise ValueError(f"Layer {on_grid_layer_name} not found")

    srs = off_layer.GetSpatialRef()

    print(f"  Found {off_layer.GetFeatureCount()} off-grid features")
    print(f"  Found {on_layer.GetFeatureCount()} on-grid features")

    # Collect all on-grid geometries and their boundaries
    on_grid_geoms = []
    for feat in on_layer:
        geom = feat.GetGeometryRef()
        if geom and not geom.IsEmpty():
            on_grid_geoms.append(geom.Clone())

    if not on_grid_geoms:
        print("  Warning: No on-grid geometries found")
        ds = None
        return output_layer_name

    # Create union of on-grid boundaries for intersection checking
    on_union = None
    for geom in on_grid_geoms:
        geom_fixed = geom.Buffer(0)
        on_union = geom_fixed if on_union is None else on_union.Union(geom_fixed)

    if on_union is None:
        print("  Warning: Failed to create on-grid union")
        ds = None
        return output_layer_name

    lines_to_write = []
    next_id = 0

    def collect_lines(geom: ogr.Geometry):
        gtype = geom.GetGeometryType()
        if gtype in (ogr.wkbLineString, ogr.wkbLineString25D):
            yield geom
        elif gtype in (ogr.wkbMultiLineString, ogr.wkbGeometryCollection):
            for i in range(geom.GetGeometryCount()):
                sub_geom = geom.GetGeometryRef(i)
                if sub_geom:
                    yield from collect_lines(sub_geom)

    processed = 0
    found_adjacent = 0

    for feat in off_layer:
        off_geom = feat.GetGeometryRef()

        if feat.GetField("is_concave") == 1:
            # Skip concave off-grid blocks
            continue

        if off_geom is None or off_geom.IsEmpty():
            continue

        processed += 1
        boundary = off_geom.GetBoundary()
        if boundary is None or boundary.IsEmpty():
            continue
        on_union_buffered = on_union.Buffer(0.01)

        if not boundary.Intersects(on_union_buffered):
            continue

        found_adjacent += 1
        try:
            shared_boundary = boundary.Intersection(on_union_buffered)
        except RuntimeError as e:
            print(f"    Warning: Failed to intersect boundaries: {e}")
            continue

        if shared_boundary is None or shared_boundary.IsEmpty():
            continue

        orig_id_idx = feat.GetFieldIndex("orig_id")
        orig_id_val = feat.GetField(orig_id_idx) if orig_id_idx != -1 else None
        geom_type = shared_boundary.GetGeometryType()
        if geom_type not in (
            ogr.wkbLineString,
            ogr.wkbLineString25D,
            ogr.wkbMultiLineString,
            ogr.wkbGeometryCollection,
        ):
            continue
        for line_part in collect_lines(shared_boundary):
            if line_part is None:
                continue
            wkt = line_part.ExportToWkt()
            if wkt:
                wkt_2d = wkt.replace(" Z ", " ").replace("LINESTRING Z", "LINESTRING")
                line_2d = ogr.CreateGeometryFromWkt(wkt_2d)
                if line_2d and not line_2d.IsEmpty() and line_2d.Length() > 0.01:
                    lines_to_write.append(
                        {
                            "geometry": line_2d,
                            "orig_id": orig_id_val if orig_id_val is not None else next_id,
                        }
                    )
                    next_id += 1
    print(f"  Processed {processed} off-grid blocks")
    print(f"  Found {found_adjacent} blocks adjacent to on-grid")
    print(f"  Collected {len(lines_to_write)} line segments")

    off_layer = None
    on_layer = None
    ds = None

    out_ds = ogr.Open(output_gpkg, 1)
    if out_ds is None:
        raise ValueError(f"Could not open {output_gpkg} for writing")
    merged_lines = []
    grouped = {}
    for line in lines_to_write:
        grouped.setdefault(line["orig_id"], []).append(line["geometry"])

    merge_successful = False
    for orig_id, geoms in grouped.items():
        shapely_lines = []
        for g in geoms:
            if g is None or g.IsEmpty():
                continue
            try:
                wkt = g.ExportToWkt()
                if wkt:
                    from shapely import wkt as shapely_wkt

                    shapely_lines.append(shapely_wkt.loads(wkt))
            except (RuntimeError, TypeError):
                continue

        if not shapely_lines:
            for geom in geoms:
                if geom and not geom.IsEmpty():
                    merged_lines.append({"geometry": geom, "orig_id": orig_id})
            continue

        try:
            dissolved = shapely_ops.unary_union(shapely_lines)
            try:
                merged = shapely_ops.linemerge(dissolved)
            except ValueError:
                merged = dissolved

            if merged.is_empty:
                for geom in geoms:
                    if geom and not geom.IsEmpty():
                        merged_lines.append({"geometry": geom, "orig_id": orig_id})
                continue

            if merged.geom_type == "LineString":
                merged_parts = [merged]
            elif merged.geom_type == "MultiLineString":
                merged_parts = list(merged.geoms)
            else:
                for geom in geoms:
                    if geom and not geom.IsEmpty():
                        merged_lines.append({"geometry": geom, "orig_id": orig_id})
                continue

            for part in merged_parts:
                ogr_geom = ogr.CreateGeometryFromWkb(part.wkb)
                if ogr_geom:
                    merged_lines.append({"geometry": ogr_geom, "orig_id": orig_id})
                    merge_successful = True
        except (RuntimeError, TypeError):
            for geom in geoms:
                if geom and not geom.IsEmpty():
                    merged_lines.append({"geometry": geom, "orig_id": orig_id})

    if merge_successful:
        print(f"  Merged {len(lines_to_write)} segments into {len(merged_lines)} features")
    else:
        print(f"  Writing {len(merged_lines)} line segments (merging skipped)")

    out_layer = create_or_replace_layer(out_ds, output_layer_name, srs, ogr.wkbLineString)
    out_layer.CreateField(ogr.FieldDefn("orig_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("id", ogr.OFTInteger))

    line_id = 1
    for line in merged_lines:
        feat = ogr.Feature(out_layer.GetLayerDefn())
        feat.SetGeometry(line["geometry"])
        feat.SetField("orig_id", line["orig_id"])
        feat.SetField("id", line_id)
        out_layer.CreateFeature(feat)
        feat = None
        line_id += 1

    out_layer = None
    out_ds = None

    print(f"  Created off-grid adjacent lines layer: {output_layer_name}")
    return output_layer_name


def extract_vertices_from_lines(
    input_gpkg: str,
    input_layer_name: str,
    output_layer_name: str,
    tol: float = 1,
) -> str:
    ds = ogr.Open(input_gpkg, 0)
    if ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    in_layer = ds.GetLayerByName(input_layer_name)
    if in_layer is None:
        raise ValueError(f"Layer {input_layer_name} not found")

    srs = in_layer.GetSpatialRef()

    points: list[dict] = []
    line_counter = 0
    seen_by_orig: dict[int, set[tuple[int, int]]] = {}

    def qkey(x: float, y: float) -> tuple[int, int]:
        return (int(round(x / tol)), int(round(y / tol)))

    def angle_deg(
        prev_xy: tuple[float, float], cur_xy: tuple[float, float], next_xy: tuple[float, float]
    ) -> float:
        px, py = prev_xy
        cx, cy = cur_xy
        nx, ny = next_xy

        ax = px - cx
        ay = py - cy
        bx = nx - cx
        by = ny - cy

        la = math.hypot(ax, ay)
        lb = math.hypot(bx, by)
        if la <= 1e-12 or lb <= 1e-12:
            return 0.0

        dot = ax * bx + ay * by
        c = dot / (la * lb)
        if c > 1.0:
            c = 1.0
        elif c < -1.0:
            c = -1.0

        return math.degrees(math.acos(c))

    for feat in in_layer:
        geom = feat.GetGeometryRef()
        if geom is None or geom.IsEmpty():
            continue

        orig_id_idx = feat.GetFieldIndex("orig_id")
        orig_id_val = feat.GetField(orig_id_idx) if orig_id_idx != -1 else None
        if orig_id_val is None:
            orig_id_val = int(feat.GetFID())

        seen = seen_by_orig.setdefault(int(orig_id_val), set())

        def handle_line(
            line_geom: ogr.Geometry, *, current_seen: set[tuple[int, int]], current_orig_id: int
        ):
            nonlocal line_counter

            n = line_geom.GetPointCount()
            if n < 2:
                line_counter += 1
                return

            coords = [(line_geom.GetX(i), line_geom.GetY(i)) for i in range(n)]
            x0, y0 = coords[0]
            xn, yn = coords[-1]
            if abs(x0 - xn) <= tol and abs(y0 - yn) <= tol:
                coords = coords[:-1]
            deduped: list[tuple[float, float]] = []
            for x, y in coords:
                if not deduped:
                    deduped.append((x, y))
                    continue
                lx, ly = deduped[-1]
                if abs(x - lx) <= tol and abs(y - ly) <= tol:
                    continue
                deduped.append((x, y))

            m = len(deduped)
            vertex_out_idx = 1

            for i, (x, y) in enumerate(deduped):
                k = qkey(x, y)
                if k in current_seen:
                    continue
                current_seen.add(k)
                if i > 0 and i < m - 1:
                    ang = angle_deg(deduped[i - 1], deduped[i], deduped[i + 1])
                else:
                    ang = 0.0

                pt_geom = ogr.Geometry(ogr.wkbPoint)
                pt_geom.AddPoint(x, y)
                points.append(
                    {
                        "geometry": pt_geom,
                        "orig_id": current_orig_id,
                        "line_id": line_counter,
                        "vertex_idx": vertex_out_idx,
                        "x": x,
                        "y": y,
                        "angle_deg": ang,
                    }
                )
                vertex_out_idx += 1

            line_counter += 1

        gtype = geom.GetGeometryType()
        if gtype in (ogr.wkbLineString, ogr.wkbLineString25D):
            handle_line(geom, current_seen=seen, current_orig_id=int(orig_id_val))
        elif gtype in (ogr.wkbMultiLineString, ogr.wkbGeometryCollection):
            for i in range(geom.GetGeometryCount()):
                sub = geom.GetGeometryRef(i)
                if sub and sub.GetGeometryType() in (ogr.wkbLineString, ogr.wkbLineString25D):
                    handle_line(sub, current_seen=seen, current_orig_id=int(orig_id_val))

    ds = None

    out_ds = ogr.Open(input_gpkg, 1)
    if out_ds is None:
        raise ValueError(f"Could not open {input_gpkg} for writing")

    out_layer = create_or_replace_layer(out_ds, output_layer_name, srs, ogr.wkbPoint)
    out_layer.CreateField(ogr.FieldDefn("orig_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("line_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("vertex_idx", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("x", ogr.OFTReal))
    out_layer.CreateField(ogr.FieldDefn("y", ogr.OFTReal))
    out_layer.CreateField(ogr.FieldDefn("angle_deg", ogr.OFTReal))

    for pt in points:
        out_feat = ogr.Feature(out_layer.GetLayerDefn())
        out_feat.SetGeometry(pt["geometry"])
        out_feat.SetField("orig_id", pt["orig_id"])
        out_feat.SetField("line_id", pt["line_id"])
        out_feat.SetField("vertex_idx", pt["vertex_idx"])
        out_feat.SetField("x", pt["x"])
        out_feat.SetField("y", pt["y"])
        out_feat.SetField("angle_deg", float(pt["angle_deg"]))
        out_layer.CreateFeature(out_feat)
        out_feat = None

    out_layer = None
    out_ds = None

    print(f"  Created vertices layer: {output_layer_name} ({len(points)} points)")
    return output_layer_name


def merge_vertices_into_lines_by_angle(
    input_gpkg: str,
    vertices_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
    angle_threshold_deg: float = 100.0,
) -> str:
    """Build line segments between consecutive vertices and merge them by group.

    Reads a point layer produced by `extract_vertices_from_lines` (must include:
    `orig_id`, `line_id`, `vertex_idx`, `angle_deg`, and point geometry).

    For each (orig_id, line_id) sequence (sorted by vertex_idx), this function:
      - Creates a line segment between every consecutive pair of vertices.
      - Assigns a `group_id` to each segment. The group continues until a "break"
        happens; a break is when the vertex angle <= angle_threshold_deg.
      - Merges all segments with the same group_id into a single LineString.
      - Calculates the total length of each merged line.

    Args:
        input_gpkg: Path to the GeoPackage containing the vertices layer.
        vertices_layer_name: Name of the point layer created by `extract_vertices_from_lines`.
        output_gpkg: Path to the GeoPackage to write the output layer to.
        output_layer_name: Name of the output line layer to create/replace.
        angle_threshold_deg: Break threshold; break occurs when vertex angle <= this value.

    Returns:
        The name of the created output layer.

    Raises:
        ValueError: If the datasets/layers cannot be opened.
    """
    ds = ogr.Open(input_gpkg, 0)
    if ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    v_layer = ds.GetLayerByName(vertices_layer_name)
    if v_layer is None:
        raise ValueError(f"Layer {vertices_layer_name} not found")

    srs = v_layer.GetSpatialRef()

    # Field indices
    defn = v_layer.GetLayerDefn()
    idx_orig = defn.GetFieldIndex("orig_id")
    idx_line = defn.GetFieldIndex("line_id")
    idx_vidx = defn.GetFieldIndex("vertex_idx")
    idx_ang = defn.GetFieldIndex("angle_deg")

    if idx_orig == -1 or idx_line == -1 or idx_vidx == -1:
        raise ValueError("Vertices layer must contain fields: orig_id, line_id, vertex_idx")

    # Group vertices by (orig_id, line_id)
    groups: dict[tuple[int, int], list[dict]] = {}

    for feat in v_layer:
        geom = feat.GetGeometryRef()
        if geom is None or geom.IsEmpty():
            continue

        orig_id = feat.GetField(idx_orig)
        line_id = feat.GetField(idx_line)
        vertex_idx = feat.GetField(idx_vidx)

        if orig_id is None or line_id is None or vertex_idx is None:
            continue

        ang = 0.0
        if idx_ang != -1:
            val = feat.GetField(idx_ang)
            if val is not None:
                ang = float(val)

        groups.setdefault((int(orig_id), int(line_id)), []).append(
            {
                "vertex_idx": int(vertex_idx),
                "x": float(geom.GetX()),
                "y": float(geom.GetY()),
                "angle_deg": float(ang),
            }
        )

    ds = None

    # Create segments and assign group IDs
    segments: list[dict] = []
    global_group_id = 1

    for (orig_id, line_id), pts in groups.items():
        pts.sort(key=lambda p: p["vertex_idx"])
        if len(pts) < 2:
            continue

        current_group_id = global_group_id
        last_break = False

        for i in range(len(pts) - 1):
            a = pts[i]
            b = pts[i + 1]

            angle_sum = a["angle_deg"] + b["angle_deg"]
            target_angle = b["angle_deg"]
            is_break = 1 if target_angle <= angle_threshold_deg else 0

            line = ogr.Geometry(ogr.wkbLineString)
            line.AddPoint(a["x"], a["y"])
            line.AddPoint(b["x"], b["y"])

            segments.append(
                {
                    "geometry": line,
                    "orig_id": orig_id,
                    "line_id": line_id,
                    "group_id": current_group_id,
                    "a_vidx": a["vertex_idx"],
                    "b_vidx": b["vertex_idx"],
                    "angle_sum": float(angle_sum),
                    "is_break": is_break,
                }
            )

            if is_break:
                last_break = True
            if last_break:
                global_group_id += 1
                current_group_id = global_group_id
                last_break = False

        # Move to next group for next line
        global_group_id += 1

    # Merge segments by group_id
    grouped_segments: dict[int, list[dict]] = {}
    for seg in segments:
        grouped_segments.setdefault(seg["group_id"], []).append(seg)

    merged_lines: list[dict] = []
    for group_id, segs in grouped_segments.items():
        # Sort segments by their vertex indices to maintain order
        segs.sort(key=lambda s: (s["orig_id"], s["line_id"], s["a_vidx"]))

        # Collect all coordinates in order
        coords = []
        for seg in segs:
            geom = seg["geometry"]
            if not coords:
                # Add first point of first segment
                coords.append((geom.GetX(0), geom.GetY(0)))
            # Add second point of each segment
            coords.append((geom.GetX(1), geom.GetY(1)))

        # Create merged LineString
        merged_line = ogr.Geometry(ogr.wkbLineString)
        for x, y in coords:
            merged_line.AddPoint(x, y)

        # Calculate total length
        total_length = merged_line.Length()

        # Get representative values from first segment
        first_seg = segs[0]
        merged_lines.append(
            {
                "geometry": merged_line,
                "orig_id": first_seg["orig_id"],
                "line_id": first_seg["line_id"],
                "group_id": group_id,
                "length": total_length,
                "num_segments": len(segs),
            }
        )

    out_ds = ogr.Open(output_gpkg, 1)
    if out_ds is None:
        raise ValueError(f"Could not open {output_gpkg} for writing")

    out_layer = create_or_replace_layer(out_ds, output_layer_name, srs, ogr.wkbLineString)
    out_layer.CreateField(ogr.FieldDefn("orig_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("line_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("group_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("length", ogr.OFTReal))
    out_layer.CreateField(ogr.FieldDefn("num_segments", ogr.OFTInteger))

    for row in merged_lines:
        feat = ogr.Feature(out_layer.GetLayerDefn())
        feat.SetGeometry(row["geometry"])
        feat.SetField("orig_id", row["orig_id"])
        feat.SetField("line_id", row["line_id"])
        feat.SetField("group_id", row["group_id"])
        feat.SetField("length", row["length"])
        feat.SetField("num_segments", row["num_segments"])
        out_layer.CreateFeature(feat)
        feat = None

    out_layer = None
    out_ds = None

    print(
        f"  Created lines layer: {output_layer_name} ({len(segments)} segments merged "
        f"into {len(merged_lines)} lines)"
    )
    return output_layer_name
