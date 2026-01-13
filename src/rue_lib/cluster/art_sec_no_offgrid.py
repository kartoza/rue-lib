# src/rue_lib/cluster/art_sec_no_offgrid.py
"""Generate arterial and secondary parts without off-grid subdivision.

This module creates block parts (corner and side parts) for arterial and
secondary road blocks by offsetting roads inward and performing boolean
operations to create the different part types.
"""

from pathlib import Path
from typing import Optional

import geopandas as gpd
import numpy as np
import pandas as pd
import pyogrio.errors
from shapely import unary_union
from shapely.errors import GEOSException
from shapely.geometry import LineString, MultiLineString, Polygon

from rue_lib.cluster.block_edges import extract_block_edges
from rue_lib.cluster.off_grid import extend_line
from rue_lib.core.definitions import ColorTypes, RoadTypes
from rue_lib.streets.geometry_utils import break_linestring_by_angle_shapely


def generate_art_sec_parts_no_offgrid(
    output_path: Path,
    blocks_layer_name: str,
    roads_layer_name: str,
    roads_buffered_layer_name: str,
    part_art_d: float,
    art_road_width_m: float,
    part_sec_d: float,
    sec_road_width_m: float,
    part_loc_d: float,
    loc_road_width_m: float,
    output_layer_name,
    ortho_direction: Optional[np.ndarray] = None,
):
    """
    Generate arterial, secondary, and local block parts without off-grid subdivision.

    This function processes blocks and creates parts based on a three-stage splitting
    process using arterial, secondary, and local roads. It performs geometric operations
    to create corner parts, side parts, and other part types by offsetting roads inward
    and performing boolean operations. Results are saved to a GeoPackage file.

    Splitting stages:
    1. Arterial roads: Blocks are split by buffered arterial roads
    2. Secondary roads: Parts from stage 1 are further split by buffered secondary roads
    3. Local roads: Parts from stage 2 are split by buffered local roads, but only if
       both resulting parts are large enough (reverts to original part if too small)

    Roads are extracted from roads_layer based on road_type attribute:
    - "road_art" for arterial roads
    - "road_sec" for secondary roads
    - "road_loc" for local roads

    Args:
        output_path: Path to the GeoPackage file containing input layers and
            where output will be saved.
        blocks_layer_name: Name of the layer containing block polygons.
        roads_layer_name: Name of the layer containing road geometries with
            road_type attributes ("road_art", "road_sec", "road_loc").
        roads_buffered_layer_name: Name of the layer containing buffered road
            geometries used for edge extraction and part classification.
        part_art_d: Offset depth for arterial road parts in meters (typically 40m).
        art_road_width_m: Width of arterial roads in meters.
        part_sec_d: Offset depth for secondary road parts in meters (typically 30m).
        sec_road_width_m: Width of secondary roads in meters.
        part_loc_d: Offset depth for local road parts in meters (typically 20m).
        loc_road_width_m: Width of local roads in meters.
        output_layer_name: Name for the output layer containing generated parts.
        ortho_direction: Optional orthogonal direction vector for perpendicular
            checks. Defaults to [0, 1, 0] if not provided.

    Returns:
        str: The output layer name.
    """
    try:
        blocks_layer = gpd.read_file(output_path, layer=blocks_layer_name)
    except pyogrio.errors.DataLayerError:
        print("No off grid on on grid.")
        return

    roads_layer = gpd.read_file(output_path, layer=roads_layer_name)
    roads_buffered_layer = gpd.read_file(output_path, layer=roads_buffered_layer_name)

    # Extract arterial, secondary, and local roads from roads_layer
    roads_art_layer = roads_layer[roads_layer.get("road_type", "") == "road_art"].copy()
    roads_sec_layer = roads_layer[roads_layer.get("road_type", "") == "road_sec"].copy()
    roads_loc_layer = roads_layer[roads_layer.get("type", "") == "road_local"].copy()

    if ortho_direction is None:
        ortho_direction = np.array([0, 1, 0])

    if blocks_layer.empty:
        return gpd.GeoDataFrame(columns=["geometry", "class", "type", "block_id"])

    all_parts = []
    edges_for_blocks = []
    edges_for_parts = []

    sec_roads_cleaned = roads_sec_layer.copy()
    road_art_layer_buffered = roads_art_layer.buffer(
        art_road_width_m / 2, cap_style="flat", join_style="mitre"
    )
    for idx, sec_road in sec_roads_cleaned.iterrows():
        sec_geom = sec_road.geometry
        for art_buffer in road_art_layer_buffered:
            try:
                sec_geom = sec_geom.difference(art_buffer)
            except (GEOSException, AttributeError, ValueError):
                pass
        sec_roads_cleaned.at[idx, "geometry"] = sec_geom

    art_roads_buffered = roads_art_layer.buffer(
        part_art_d + art_road_width_m, cap_style="flat", join_style="mitre"
    )

    sec_roads_buffered = []
    for sec_geom in sec_roads_cleaned.geometry:
        buffered = sec_geom.buffer(
            part_sec_d + sec_road_width_m, cap_style="flat", join_style="mitre"
        )
        sec_roads_buffered.append(buffered)

    loc_roads_broken = []
    segment_id = 0
    for _idx, loc_road in roads_loc_layer.iterrows():
        loc_geom = loc_road.geometry

        if isinstance(loc_geom, LineString):
            coords = list(loc_geom.coords)
            for i in range(len(coords) - 1):
                segment = LineString([coords[i], coords[i + 1]])
                loc_roads_broken.append({"geometry": segment, "segment_id": segment_id})
                segment_id += 1
        else:
            loc_roads_broken.append({"geometry": loc_geom, "segment_id": segment_id})
            segment_id += 1

    loc_segments_buffered = []
    for segment_data in loc_roads_broken:
        segment_geom = segment_data["geometry"]
        seg_id = segment_data["segment_id"]
        buffered = segment_geom.buffer(loc_road_width_m + 0.1, cap_style="flat")
        loc_segments_buffered.append({"geometry": buffered, "segment_id": seg_id})

    block_lines_intersecting = []

    all_blocks = []
    for idx, block_row in blocks_layer.iterrows():
        block = block_row.geometry
        block_id = block_row.get("block_id", idx)
        if pd.isna(block_id):
            block_id = idx
        if isinstance(block, Polygon):
            all_blocks.append({"geometry": block, "block_id": block_id})
        else:
            union = unary_union(block.geoms)
            if isinstance(union, Polygon):
                all_blocks.append({"geometry": union, "block_id": block_id})
            else:
                for _idx, sub_block in enumerate(block.geoms):
                    if isinstance(sub_block, Polygon):
                        all_blocks.append(
                            {"geometry": sub_block, "block_id": block_id + len(blocks_layer) + _idx}
                        )

    for idx, block_row in enumerate(all_blocks):
        block = block_row["geometry"]
        block_id = block_row.get("block_id", idx)
        if pd.isna(block_id):
            block_id = idx

        if isinstance(block, Polygon):
            boundary = block.boundary
        else:
            union = unary_union(block.geoms)
            if isinstance(union, Polygon):
                boundary = union.boundary
            else:
                continue

        line_groups = []

        if isinstance(boundary, LineString):
            coords = list(boundary.coords)
            current_group = []

            for i in range(len(coords) - 1):
                segment = LineString([coords[i], coords[i + 1]])
                intersects_any = False
                for seg_data in loc_segments_buffered:
                    seg_buffer = seg_data["geometry"]
                    if segment.centroid.intersects(seg_buffer):
                        intersects_any = True
                        break

                if intersects_any:
                    if not current_group:
                        current_group.append(coords[i])
                    current_group.append(coords[i + 1])
                else:
                    if len(current_group) >= 2:
                        line_groups.append(current_group)
                    current_group = []

            if len(current_group) >= 2:
                line_groups.append(current_group)

        if line_groups:
            merge_threshold = 1.0
            merged = False

            while True:
                merged = False
                for i in range(len(line_groups)):
                    if i >= len(line_groups):
                        break
                    for j in range(i + 1, len(line_groups)):
                        if j >= len(line_groups):
                            break

                        group1 = line_groups[i]
                        group2 = line_groups[j]

                        dist1 = (
                            (group1[-1][0] - group2[0][0]) ** 2
                            + (group1[-1][1] - group2[0][1]) ** 2
                        ) ** 0.5
                        dist2 = (
                            (group1[-1][0] - group2[-1][0]) ** 2
                            + (group1[-1][1] - group2[-1][1]) ** 2
                        ) ** 0.5
                        dist3 = (
                            (group1[0][0] - group2[0][0]) ** 2 + (group1[0][1] - group2[0][1]) ** 2
                        ) ** 0.5
                        dist4 = (
                            (group1[0][0] - group2[-1][0]) ** 2
                            + (group1[0][1] - group2[-1][1]) ** 2
                        ) ** 0.5

                        if dist1 <= merge_threshold:
                            line_groups[i] = group1 + group2[1:]
                            line_groups.pop(j)
                            merged = True
                            break
                        elif dist2 <= merge_threshold:
                            line_groups[i] = group1 + list(reversed(group2))[1:]
                            line_groups.pop(j)
                            merged = True
                            break
                        elif dist3 <= merge_threshold:
                            line_groups[i] = list(reversed(group1)) + group2[1:]
                            line_groups.pop(j)
                            merged = True
                            break
                        elif dist4 <= merge_threshold:
                            line_groups[i] = group2 + group1[1:]
                            line_groups.pop(j)
                            merged = True
                            break

                    if merged:
                        break

                if not merged:
                    break

            if len(line_groups) == 1:
                line_geom = LineString(line_groups[0])
                block_lines_intersecting.append({"geometry": line_geom, "block_id": block_id})
            elif len(line_groups) > 1:
                lines = [LineString(group) for group in line_groups if len(group) >= 2]
                if len(lines) == 1:
                    line_geom = lines[0]
                else:
                    line_geom = MultiLineString(lines)
                block_lines_intersecting.append({"geometry": line_geom, "block_id": block_id})

    local_lines_by_block_id = {}
    if block_lines_intersecting:
        for block_line in block_lines_intersecting:
            block_id = block_line["block_id"]
            broken_lines = break_linestring_by_angle_shapely(
                block_line["geometry"], angle_threshold=60
            )
            local_lines_by_block_id[block_id] = broken_lines

    for idx, block_row in enumerate(all_blocks):
        block = block_row["geometry"]

        # Get block ID
        block_id = block_row.get("block_id", idx)
        if pd.isna(block_id):
            block_id = idx

        parts_to_add = []

        intersecting_art_buffers = []
        for art_buffer in art_roads_buffered:
            if block.intersects(art_buffer):
                intersecting_art_buffers.append(art_buffer)

        intersecting_sec_buffers = []
        for sec_buffer in sec_roads_buffered:
            if block.intersects(sec_buffer):
                intersecting_sec_buffers.append(sec_buffer)

        if intersecting_art_buffers or intersecting_sec_buffers:
            try:
                intermediate_parts = []
                remaining_block = block
                if intersecting_art_buffers:
                    union_art_buffer = intersecting_art_buffers[0]
                    for buffer in intersecting_art_buffers[1:]:
                        union_art_buffer = union_art_buffer.union(buffer)

                    art_buffer_part = remaining_block.intersection(union_art_buffer)
                    remaining_block = remaining_block.difference(union_art_buffer)
                    intermediate_parts.append({"geometry": art_buffer_part, "is_art_or_sec": "art"})

                intermediate_parts.append({"geometry": remaining_block, "is_art_or_sec": None})

                if intersecting_sec_buffers:
                    union_sec_buffer = intersecting_sec_buffers[0]
                    for buffer in intersecting_sec_buffers[1:]:
                        union_sec_buffer = union_sec_buffer.union(buffer)

                    for part_data in intermediate_parts:
                        part = part_data["geometry"]
                        is_art = part_data["is_art_or_sec"] == "art"

                        if not part.intersects(union_sec_buffer):
                            parts_to_add.append(part_data)
                            continue

                        sec_buffer_part = part.intersection(union_sec_buffer)
                        part_remaining = part.difference(union_sec_buffer)

                        if isinstance(sec_buffer_part, Polygon) and sec_buffer_part.area > 1.0:
                            # If original part was art, this is art_sec; otherwise it's sec
                            sec_type = "art" if is_art else "sec"
                            parts_to_add.append(
                                {"geometry": sec_buffer_part, "is_art_or_sec": sec_type}
                            )
                        if isinstance(part_remaining, Polygon) and part_remaining.area > 1.0:
                            # Remaining part keeps original type
                            remaining_type = part_data["is_art_or_sec"]
                            parts_to_add.append(
                                {"geometry": part_remaining, "is_art_or_sec": remaining_type}
                            )
                else:
                    parts_to_add = intermediate_parts

            except (GEOSException, AttributeError, ValueError):
                parts_to_add = [{"geometry": block, "is_art_or_sec": None}]
        else:
            parts_to_add = [{"geometry": block, "is_art_or_sec": None}]

        if block_id in local_lines_by_block_id:
            local_lines = local_lines_by_block_id[block_id]

            final_parts = []
            is_corner_block = len(parts_to_add) >= 4
            min_area_factor = 4 if is_corner_block else 3

            # Get art/sec parts for local line filtering
            art_sec_parts = [
                part_data["geometry"]
                for part_data in parts_to_add
                if part_data["is_art_or_sec"] is not None
            ]

            filtered_local_lines = []
            if art_sec_parts:
                for local_line in local_lines:
                    line_centroid = local_line.centroid
                    for part_geom in art_sec_parts:
                        if line_centroid.buffer(1).intersects(part_geom):
                            filtered_local_lines.append(local_line)
                            break

                # Sort by length
                local_lines_sorted = sorted(
                    filtered_local_lines, key=lambda line: line.length, reverse=True
                )

                local_lines = local_lines_sorted
            else:
                local_lines = []

            for part_data in parts_to_add:
                current_part = part_data["geometry"]
                split_parts = []

                for local_line in local_lines[:2]:
                    local_off = extend_line(local_line, loc_road_width_m).buffer(
                        part_loc_d, join_style="mitre", cap_style="flat"
                    )

                    buffer_part_for_calculation = block_row["geometry"].intersection(local_off)
                    if (
                        buffer_part_for_calculation.area * min_area_factor
                        > block_row["geometry"].area
                    ):
                        continue

                    try:
                        buffer_part = current_part.intersection(local_off)
                        remaining_part = current_part.difference(local_off)
                        split_parts.append({"geometry": buffer_part, "is_art_or_sec": None})
                        current_part = remaining_part
                    except (GEOSException, AttributeError, ValueError):
                        continue

                final_parts.extend(split_parts)
                if current_part.area > 1.0:
                    final_parts.append({"geometry": current_part, "is_art_or_sec": None})

            parts_to_add = final_parts

        index = 0
        for part_data in parts_to_add:
            part_geom = part_data["geometry"]
            part_gdf = gpd.GeoDataFrame(
                [{"geometry": part_geom, "block_id": block_id}], crs=blocks_layer.crs
            )
            index += 1
            part_edges = extract_block_edges(part_gdf, roads_buffered_layer, default_type=None)

            unique_types = set(part_edges["road_type"])
            local_count = (part_edges["road_type"] == RoadTypes.Local).sum()

            if RoadTypes.Artery in unique_types:
                part_type = "art_sec" if RoadTypes.Secondary in unique_types else "art"
            elif RoadTypes.Secondary in unique_types:
                part_type = "sec_loc" if RoadTypes.Local in unique_types else "sec"
            else:
                if local_count >= 2 and len(part_edges) >= 4:
                    local_edges = part_edges[part_edges["road_type"] == RoadTypes.Local].copy()

                    edges_to_merge = set()
                    min_edge_length = 5.0
                    local_indices = local_edges.index.tolist()
                    for i, idx in enumerate(local_indices):
                        edge = part_edges.loc[idx]
                        edge_geom = edge.geometry
                        edge_coords = list(edge_geom.coords)

                        for j in range(i + 1, len(local_indices)):
                            idx2 = local_indices[j]
                            edge2 = part_edges.loc[idx2]
                            edge2_geom = edge2.geometry
                            edge2_coords = list(edge2_geom.coords)

                            shared_point = None
                            p1, p2 = None, None

                            if edge_coords[-1] == edge2_coords[0]:
                                shared_point = edge_coords[-1]
                                p1 = edge_coords[-2] if len(edge_coords) >= 2 else edge_coords[0]
                                p2 = edge2_coords[1] if len(edge2_coords) >= 2 else edge2_coords[-1]
                            elif edge_coords[0] == edge2_coords[-1]:
                                shared_point = edge_coords[0]
                                p1 = edge2_coords[-2] if len(edge2_coords) >= 2 else edge2_coords[0]
                                p2 = edge_coords[1] if len(edge_coords) >= 2 else edge_coords[-1]

                            if shared_point and p1 and p2:
                                v1 = np.array([p1[0] - shared_point[0], p1[1] - shared_point[1]])
                                v2 = np.array([p2[0] - shared_point[0], p2[1] - shared_point[1]])
                                v1_norm = v1 / (np.linalg.norm(v1) + 1e-10)
                                v2_norm = v2 / (np.linalg.norm(v2) + 1e-10)
                                dot_product = np.clip(np.dot(v1_norm, v2_norm), -1.0, 1.0)
                                angle = np.degrees(np.arccos(dot_product))

                                # Consider as one edge if angle > 100 degrees
                                # or one edge is too small
                                if (
                                    angle > 100
                                    or edge_geom.length < min_edge_length
                                    or edge2_geom.length < min_edge_length
                                ):
                                    edges_to_merge.add((min(idx, idx2), max(idx, idx2)))

                    # Adjust local count by subtracting merged edge pairs
                    local_count -= len(edges_to_merge)
                part_type = "loc_loc" if local_count >= 2 else "loc"

            # TODO:
            #  Remove: Save it to layer
            for _idx, edge in part_edges.iterrows():
                edges_for_blocks.append(edge)
                edges_for_parts.append(edge)

            try:
                color = ColorTypes[part_type]
            except KeyError:
                color = "rgb(255,255,255)"

            all_parts.append(
                {
                    "geometry": part_geom,
                    "class": "part",
                    "block_id": block_id,
                    "block_type": block_row.get("block_type", ""),
                    "site": block_row.get("site", ""),
                    "type": part_type,
                    "color": color,
                }
            )
    # Create GeoDataFrame
    gdf_out = gpd.GeoDataFrame(edges_for_blocks, crs=blocks_layer.crs)
    gdf_out.to_file(output_path, layer=output_layer_name + "-edges", driver="GPKG")
    gdf_out = gpd.GeoDataFrame(edges_for_parts, crs=blocks_layer.crs)
    gdf_out.to_file(output_path, layer=output_layer_name + "-edges-parts", driver="GPKG")

    # Create GeoDataFrame
    gdf_out = gpd.GeoDataFrame(all_parts, crs=blocks_layer.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")
    return output_layer_name
