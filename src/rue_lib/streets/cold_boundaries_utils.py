from __future__ import annotations

from pathlib import Path

import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Point, Polygon


def extract_line_vertices(
    input_path: Path,
    line_layer_name: str,
    output_path: Path,
    output_layer_name: str = "18_local_roads_vertices",
    dedup_rounding: int = 6,
) -> str | None:
    """Extract unique vertices from a line layer and write them as points."""
    gdf_lines = gpd.read_file(input_path, layer=line_layer_name)
    if gdf_lines.empty:
        print(f"  Warning: line layer '{line_layer_name}' is empty; no vertices extracted")
        return None

    seen: set[tuple[float, float]] = set()
    points: list[Point] = []

    for geom in gdf_lines.geometry:
        if geom is None or geom.is_empty:
            continue
        if isinstance(geom, LineString):
            lines = [geom]
        elif isinstance(geom, MultiLineString):
            lines = list(geom.geoms)
        else:
            continue

        for line in lines:
            for x, y, *_ in line.coords:
                key = (round(x, dedup_rounding), round(y, dedup_rounding))
                if key in seen:
                    continue
                seen.add(key)
                points.append(Point(x, y))

    if not points:
        print(f"  Warning: no vertices extracted from '{line_layer_name}'")
        return None

    gdf_vertices = gpd.GeoDataFrame(
        {"vertex_id": list(range(len(points)))}, geometry=points, crs=gdf_lines.crs
    )
    gdf_vertices.to_file(output_path, layer=output_layer_name, driver="GPKG")
    print(
        f"  Extracted {len(points)} vertices from '{line_layer_name}' \n"
        f"  into layer '{output_layer_name}'"
    )
    return output_layer_name


def extract_cold_boundary_lines_from_vertices(
    output_path: Path,
    cold_boundaries_layer: str,
    vertices_layer: str,
    local_roads_layer: str,
    output_layer_name: str = "14_cold_boundary_lines",
    snap_distance: float = 1.0,
) -> str | None:
    """Split cold boundary polygons into lines using nearby vertices as break points."""
    gdf_cold = gpd.read_file(output_path, layer=cold_boundaries_layer)
    gdf_vertices = gpd.read_file(output_path, layer=vertices_layer)
    gdf_local_roads = gpd.read_file(output_path, layer=local_roads_layer)

    if gdf_cold.empty:
        print(f"  Warning: cold boundaries layer '{cold_boundaries_layer}' is empty")
        return None
    if gdf_vertices.empty:
        print(f"  Warning: vertices layer '{vertices_layer}' is empty")
        return None

    vertices_per_cold = {}

    for idx, row in gdf_cold.iterrows():
        geom = row.geometry
        cold_id = row.get("id", idx)
        if geom is None or geom.is_empty:
            continue
        cold_vertices = []
        for _vertex_idx, vrow in gdf_vertices.iterrows():
            vgeom = vrow.geometry
            if vgeom is None or vgeom.is_empty:
                continue
            if vgeom.buffer(0.1).intersects(geom):
                cold_vertices.append(vrow)
        vertices_per_cold[cold_id] = cold_vertices

    out_lines: list[LineString] = []
    records: list[dict] = []
    for cold_id, vertices in vertices_per_cold.items():
        print(f"  Cold boundary ID {cold_id} has {len(vertices)} associated vertices")

        if len(vertices) < 2:
            print(f"    Skipping cold boundary {cold_id}: insufficient vertices ({len(vertices)})")
            continue

        # Get the cold boundary polygon
        cold_geom = gdf_cold[gdf_cold.get("id", gdf_cold.index) == cold_id].iloc[0].geometry

        if cold_geom is None or cold_geom.is_empty:
            print(f"    Skipping cold boundary {cold_id}: invalid geometry")
            continue

        # Get the boundary of the polygon (exterior ring)
        if isinstance(cold_geom, Polygon):
            boundary = cold_geom.exterior
        else:
            print(f"    Skipping cold boundary {cold_id}: not a polygon")
            continue

        # Sort vertices by their position along the boundary
        vertices_with_position = []
        for v in vertices:
            if v.geometry is None or v.geometry.is_empty:
                continue
            # Project the vertex onto the boundary to get its position
            position = boundary.project(v.geometry)
            vertices_with_position.append((position, v.geometry.coords[0]))

        if len(vertices_with_position) < 2:
            print(f"    Skipping cold boundary {cold_id}: insufficient valid vertices")
            continue

        vertices_with_position.sort(key=lambda x: x[0])

        coords = [coord for _, coord in vertices_with_position]

        is_closing = False
        for _, road_row in gdf_local_roads.iterrows():
            road_geom = road_row.geometry
            if road_geom is None or road_geom.is_empty:
                continue
            line = LineString([coords[0], coords[-1]])
            line_length = line.length
            line_buffered = line.buffer(1.0)
            if line_buffered.intersects(road_geom):
                # Check if the intersection is significant (at least 70% of line length)
                intersection = line_buffered.intersection(road_geom)
                if not intersection.is_empty:
                    # Get the length of the intersection
                    if hasattr(intersection, "length"):
                        intersection_length = intersection.length
                    else:
                        # For non-linear geometries (like polygons), use the line's intersection
                        line_intersection = line.intersection(road_geom.buffer(1.0))
                        intersection_length = (
                            line_intersection.length if hasattr(line_intersection, "length") else 0
                        )

                    # Check if intersection is at least 70% of original line length
                    if intersection_length >= 0.7 * line_length:
                        is_closing = True
                        break

        if is_closing:
            if coords[0] != coords[-1]:
                coords.append(coords[0])

        line_string = LineString(coords)
        out_lines.append(line_string)
        records.append(
            {
                "cold_id": cold_id,
                "segment_idx": 0,
                "segment_length": float(line_string.length),
                "num_vertices": len(coords),
            }
        )

    if not out_lines:
        print("  Warning: no cold boundary lines were created from vertices")
        return None

    gdf_out = gpd.GeoDataFrame(records, geometry=out_lines, crs=gdf_cold.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")
    print(f"  Saved {len(out_lines)} cold boundary line segments to layer: {output_layer_name}")
    return output_layer_name
