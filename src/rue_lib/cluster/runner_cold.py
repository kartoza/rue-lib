# src/rue_lib/cluster/runner_warm.py
"""Generate warm blocks with off-grid subdivision and partitioning."""

from __future__ import annotations

from pathlib import Path

from osgeo import ogr

from rue_lib.cluster.cold.expand_roads_buffer import (
    clip_buffered_lines_to_cold_grid,
    create_buffered_lines_from_boundary_lines,
)
from rue_lib.cluster.cold.subdiv_block import (
    find_concave_points,
    subdivide_blocks_by_concave_points,
)
from rue_lib.cluster.config import ClusterConfig
from rue_lib.cluster.helpers import convert_polygonz_to_polygon
from rue_lib.core.definitions import BlockTypes
from rue_lib.streets.operations import extract_by_expression


def generate_cold(
    cfg: ClusterConfig, output_gpkg: Path, input_blocks_layer_name: str, roads_layer_name: str
):
    """
    Generate cold blocks with off-grid subdivision and partitioning.

    This function processes blocks to create:
    - Inner off-grid areas by offsetting edges inward
    - Frame parts (perimeter around off-grid)
    - Corner and side parts from block decomposition
    - Subdivided plots within off-grid areas
    - On-grid parts for arterial and secondary roads

    Args:
        cfg: ClusterConfig with partition settings
        output_gpkg: Path to output GeoPackage
        input_blocks_layer_name: Name of input blocks layer
        roads_layer_name: Name of roads layer
    """
    _part_art_d = cfg.on_grid_partition_depth_arterial_roads
    _part_sec_d = cfg.on_grid_partition_depth_secondary_roads
    _part_loc_d = cfg.on_grid_partition_depth_local_roads

    # TODO:
    #  We use the same width for off-grid and on-grid plots for now.
    _part_og_w = cfg.off_grid_cluster_width

    output_path = str(output_gpkg)
    print("==============================================================")
    print("COLD BLOCK")
    print("==============================================================")
    print("Step 1: Generate inner part of off grid blocks...")
    cold_grid_layer_name = "200_cold_grid"
    extract_by_expression(
        output_path,
        input_blocks_layer_name,
        f"type = '{BlockTypes.COLD_GRID}'",
        output_path,
        cold_grid_layer_name,
    )
    # Convert from PolygonZ to Polygon
    convert_polygonz_to_polygon(output_path, cold_grid_layer_name)

    print("\nStep 2: Erase roads buffer from cold grid...")
    erased_layer_name = "201_cold_grid_erased"
    erase_roads_from_cold_grid(
        output_path,
        cold_grid_layer_name,
        roads_layer_name,
        output_path,
        erased_layer_name,
    )

    print("\nStep 3: Extract boundary lines adjacent to roads...")
    boundary_points_layer_name = "202_cold_boundary_points"
    boundary_points_layer_name, boundary_lines_layer_name = extract_road_adjacent_vertices(
        output_path,
        erased_layer_name,
        roads_layer_name,
        output_path,
        boundary_points_layer_name,
    )

    print("\nStep 6: Create buffered lines from boundary lines...")
    buffered_lines_layer_name = "206_buffered_lines"
    create_buffered_lines_from_boundary_lines(
        output_path,
        boundary_lines_layer_name,
        erased_layer_name,
        output_path,
        buffered_lines_layer_name,
        cfg,
    )

    print("\nStep 7: Clip buffered lines to cold grid...")
    clipped_lines_layer_name = "207_clipped_buffered_lines"
    clip_buffered_lines_to_cold_grid(
        output_path,
        buffered_lines_layer_name,
        erased_layer_name,
        output_path,
        clipped_lines_layer_name,
    )

    print("\nStep 4: Find concave points from boundary...")
    concave_points_layer_name = "203_concave_points"
    find_concave_points(
        output_path,
        erased_layer_name,
        boundary_points_layer_name,
        output_path,
        concave_points_layer_name,
    )

    print("\nStep 5: Subdivide blocks at concave corners...")
    cutting_lines_layer_name = "204_subdivided"
    subdivide_blocks_by_concave_points(
        output_path,
        erased_layer_name,
        concave_points_layer_name,
        boundary_points_layer_name,
        output_path,
        cutting_lines_layer_name,
        cfg.road_local_width_m,
        clipped_lines_layer_name,
    )

    return cold_grid_layer_name


def erase_roads_from_cold_grid(
    input_gpkg: str,
    cold_grid_layer_name: str,
    roads_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> str:
    """
    Erase roads buffer from cold grid blocks.

    Args:
        input_gpkg: Path to input GeoPackage
        cold_grid_layer_name: Name of cold grid layer
        roads_layer_name: Name of roads buffer layer
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output erased layer

    Returns:
        Name of the output layer
    """
    # Step 1: Read input data and collect geometries
    ds = ogr.Open(input_gpkg, 0)
    if ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    cold_layer = ds.GetLayerByName(cold_grid_layer_name)
    if cold_layer is None:
        raise ValueError(f"Layer {cold_grid_layer_name} not found")

    roads_layer = ds.GetLayerByName(roads_layer_name)
    if roads_layer is None:
        raise ValueError(f"Layer {roads_layer_name} not found")

    # Get SRS and layer definition
    srs = cold_layer.GetSpatialRef()
    cold_defn = cold_layer.GetLayerDefn()

    # Collect all road geometries into one union
    road_union = None
    for road_feat in roads_layer:
        road_geom = road_feat.GetGeometryRef()
        if road_geom:
            # Fix invalid geometries using Buffer(0)
            road_geom_fixed = road_geom.Buffer(0)
            if road_union is None:
                road_union = road_geom_fixed
            else:
                road_union = road_union.Union(road_geom_fixed)

    print(f"  Processing {cold_layer.GetFeatureCount()} cold grid blocks...")

    # Collect cold grid features with geometries
    cold_features = []
    for cold_feat in cold_layer:
        cold_geom = cold_feat.GetGeometryRef()
        if cold_geom is None:
            continue

        # Clone geometry and attributes
        cold_features.append(
            {
                "geometry": cold_geom.Clone(),
                "attributes": {i: cold_feat.GetField(i) for i in range(cold_defn.GetFieldCount())},
            }
        )

    # Close the input dataset
    cold_layer = None
    roads_layer = None
    ds = None

    # Step 2: Open for writing and create output
    out_ds = ogr.Open(output_gpkg, 1)
    if out_ds is None:
        raise ValueError(f"Could not open {output_gpkg} for writing")

    # Delete existing layer if it exists
    for i in range(out_ds.GetLayerCount()):
        layer = out_ds.GetLayerByIndex(i)
        if layer.GetName() == output_layer_name:
            out_ds.DeleteLayer(i)
            break

    # Create output layer
    out_layer = out_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)

    # Copy fields
    for i in range(cold_defn.GetFieldCount()):
        field_defn = cold_defn.GetFieldDefn(i)
        out_layer.CreateField(field_defn)

    # Process and write features
    for feat_data in cold_features:
        cold_geom = feat_data["geometry"]

        # Fix invalid geometry using Buffer(0)
        cold_geom_fixed = cold_geom.Buffer(0)

        # Erase roads from cold grid block
        if road_union:
            try:
                erased_geom = cold_geom_fixed.Difference(road_union)
            except RuntimeError as e:
                print(f"    Warning: Failed to erase roads from block, keeping original: {e}")
                erased_geom = cold_geom_fixed
        else:
            erased_geom = cold_geom_fixed

        if erased_geom and not erased_geom.IsEmpty():
            # Create output feature
            out_feat = ogr.Feature(out_layer.GetLayerDefn())
            out_feat.SetGeometry(erased_geom)

            # Copy attributes
            for i, value in feat_data["attributes"].items():
                out_feat.SetField(i, value)

            out_layer.CreateFeature(out_feat)
            out_feat = None

    # Clean up
    out_layer = None
    out_ds = None

    print(f"  Created layer: {output_layer_name}")
    return output_layer_name


def extract_road_adjacent_vertices(
    input_gpkg: str,
    erased_grid_layer_name: str,
    roads_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> tuple[str, str]:
    """
    Extract vertices from erased cold grid boundary that intersect with the roads buffer.

    Writes both:
    - A points layer containing the road-adjacent vertices
    - A lines layer containing the boundary lines that touch the roads buffer

    Returns:
        (points_layer_name, lines_layer_name)
    """
    ds = ogr.Open(input_gpkg, 0)
    if ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    grid_layer = ds.GetLayerByName(erased_grid_layer_name)
    if grid_layer is None:
        raise ValueError(f"Layer {erased_grid_layer_name} not found")

    roads_layer = ds.GetLayerByName(roads_layer_name)
    if roads_layer is None:
        raise ValueError(f"Layer {roads_layer_name} not found")

    srs = grid_layer.GetSpatialRef()

    roads_data = []
    for road_feat in roads_layer:
        road_geom = road_feat.GetGeometryRef()
        road_type = road_feat.GetField("type")
        if road_geom:
            roads_data.append({"geometry": road_geom.Clone(), "type": road_type})

    print(f"  Processing {grid_layer.GetFeatureCount()} blocks...")

    vertices_to_write = []
    lines_layer_name = (
        output_layer_name.replace("points", "lines")
        if "points" in output_layer_name
        else f"{output_layer_name}_lines"
    )
    block_id = 0
    for grid_feat in grid_layer:
        block_id += 1
        grid_geom = grid_feat.GetGeometryRef()
        if grid_geom is None:
            continue

        boundary = grid_geom.GetBoundary()
        if boundary is None:
            continue

        if boundary.GetGeometryType() == ogr.wkbLineString:
            lines = [boundary]
        elif boundary.GetGeometryType() == ogr.wkbMultiLineString:
            lines = [boundary.GetGeometryRef(i) for i in range(boundary.GetGeometryCount())]
        else:
            continue

        for line in lines:
            point_count = line.GetPointCount()
            for i in range(point_count):
                x = line.GetX(i)
                y = line.GetY(i)

                point = ogr.Geometry(ogr.wkbPoint)
                point.AddPoint(x, y)

                for road_data in roads_data:
                    if point.Buffer(0.1).Intersects(road_data["geometry"]):
                        vertices_to_write.append(
                            {
                                "x": x,
                                "y": y,
                                "block_id": block_id,
                                "vertex_id": i,
                                "road_type": road_data["type"],
                            }
                        )
                        break

    grid_layer = None
    roads_layer = None
    ds = None

    # Build line segments by connecting adjacent road-touching vertices (per block and road type)
    lines_to_write = []
    vertices_grouped: dict[tuple[int, str], list[dict]] = {}
    for vertex in vertices_to_write:
        key = (vertex["block_id"], vertex["road_type"])
        vertices_grouped.setdefault(key, []).append(vertex)

    for (block_id, road_type), verts in vertices_grouped.items():
        verts_sorted = sorted(verts, key=lambda v: v["vertex_id"])
        current = []
        prev_id = None

        def flush(seq, b_id=block_id, r_type=road_type):
            if len(seq) < 2:
                return
            line = ogr.Geometry(ogr.wkbLineString)
            for pt in seq:
                line.AddPoint(pt["x"], pt["y"])
            lines_to_write.append(
                {
                    "geometry": line,
                    "block_id": b_id,
                    "road_type": r_type,
                }
            )

        for v in verts_sorted:
            if prev_id is None or v["vertex_id"] == prev_id + 1:
                current.append(v)
            else:
                flush(current)
                current = [v]
            prev_id = v["vertex_id"]

        flush(current)

    out_ds = ogr.Open(output_gpkg, 1)
    if out_ds is None:
        raise ValueError(f"Could not open {output_gpkg} for writing")

    # Delete existing layers if present
    for i in range(out_ds.GetLayerCount()):
        layer = out_ds.GetLayerByIndex(i)
        if layer.GetName() == lines_layer_name:
            out_ds.DeleteLayer(i)
            break

    for i in range(out_ds.GetLayerCount()):
        layer = out_ds.GetLayerByIndex(i)
        if layer.GetName() == output_layer_name:
            out_ds.DeleteLayer(i)
            break

    # Write lines layer
    lines_layer = out_ds.CreateLayer(lines_layer_name, srs, ogr.wkbLineString)
    lines_layer.CreateField(ogr.FieldDefn("block_id", ogr.OFTInteger))
    lines_layer.CreateField(ogr.FieldDefn("road_type", ogr.OFTString))

    for line in lines_to_write:
        out_feat = ogr.Feature(lines_layer.GetLayerDefn())
        out_feat.SetGeometry(line["geometry"])
        out_feat.SetField("block_id", line["block_id"])
        out_feat.SetField("road_type", line["road_type"])
        lines_layer.CreateFeature(out_feat)
        out_feat = None

    lines_layer = None

    # Write points layer
    for i in range(out_ds.GetLayerCount()):
        layer = out_ds.GetLayerByIndex(i)
        if layer.GetName() == output_layer_name:
            out_ds.DeleteLayer(i)
            break

    out_layer = out_ds.CreateLayer(output_layer_name, srs, ogr.wkbPoint)

    out_layer.CreateField(ogr.FieldDefn("block_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("vertex_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("road_type", ogr.OFTString))

    for vertex in vertices_to_write:
        point = ogr.Geometry(ogr.wkbPoint)
        point.AddPoint(vertex["x"], vertex["y"])

        out_feat = ogr.Feature(out_layer.GetLayerDefn())
        out_feat.SetGeometry(point)
        out_feat.SetField("block_id", vertex["block_id"])
        out_feat.SetField("vertex_id", vertex["vertex_id"])
        out_feat.SetField("road_type", vertex["road_type"])

        out_layer.CreateFeature(out_feat)
        out_feat = None

    out_layer = None
    out_ds = None

    total_vertices = len(vertices_to_write)
    total_lines = len(lines_to_write)
    print(f"  Extracted {total_vertices} vertices from road-adjacent boundaries")
    print(f"  Created lines layer: {lines_layer_name} ({total_lines} features)")
    print(f"  Created points layer: {output_layer_name}")
    return output_layer_name, lines_layer_name
