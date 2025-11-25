# src/rue_lib/cluster/parts_generator.py
"""Generate block parts (corner, side, off-grid) for all blocks."""

import geopandas as gpd
from shapely.geometry import Polygon

from rue_lib.cluster.block_parts import create_block_parts_from_off_grid
from rue_lib.cluster.off_grid import create_off_grid_area


def generate_parts_with_off_grid(
    blocks: gpd.GeoDataFrame,
    roads: gpd.GeoDataFrame,
    part_art_d: float = 40.0,
    part_sec_d: float = 30.0,
    part_loc_d: float = 20.0,
    part_og_d: float = 140.0,
    part_og_w: float = 140.0,
    angle_threshold: float = 155.0,
    corner_distance: float = 1000.0,
) -> gpd.GeoDataFrame:
    """
    Generate "parts" (sub-blocks) for all blocks using the off-grid area.

    Logic for each block:
    1. Build extended road polylines along its perimeter
    2. Offset those roads inward to create an "off-grid" polygon
    3. If the off-grid area is too small, keep block as-is with class="block"
    4. Otherwise, split the block into parts based on that off-grid polygon:
       - corner parts (at sharp corners)
       - side parts (strips along edges)
       - off-grid center
       and tag them as class="part" with inherited attributes

    Args:
        blocks: GeoDataFrame with block polygons
        (should have 'type', 'site', 'block_id', 'block_type')
        roads: GeoDataFrame with roads (must have 'road_type' column)
        part_art_d: Depth of partition along arterial roads
        part_sec_d: Depth of partition along secondary roads
        part_loc_d: Depth of partition along local roads
        part_og_d: Off-grid depth (for minimum area threshold)
        part_og_w: Off-grid width (for minimum area threshold)
        angle_threshold: Minimum angle (degrees) to consider a corner "sharp"
        corner_distance: Distance to extend corner fan perpendiculars

    Returns:
        GeoDataFrame with all parts (original blocks and generated parts)

    Example:
        >>> blocks = gpd.read_file("blocks.geojson")
        >>> roads = gpd.read_file("roads.geojson")
        >>> parts = generate_parts_with_off_grid(blocks, roads)
    """
    all_parts = []
    min_useful_area = part_og_d * part_og_w * 0.5

    print(f"Processing {len(blocks)} blocks...")
    print(f"  Minimum off-grid area threshold: {min_useful_area:.2f} m²")

    for idx, block_row in blocks.iterrows():
        block = block_row.geometry

        if not isinstance(block, Polygon):
            print(f"  Warning: Block {idx} is not a Polygon, skipping")
            continue

        # Get block attributes
        site_id = block_row.get("site", None)
        block_id = block_row.get("block_id", idx)
        block_type = block_row.get("block_type", "unknown")

        # Step 1 & 2: Create off-grid area
        off_grid = create_off_grid_area(block, roads, part_art_d, part_sec_d, part_loc_d)

        # If no off-grid created, keep block as-is
        if off_grid is None:
            all_parts.append(
                {
                    "geometry": block,
                    "class": "block",
                    "part_type": "block",
                    "site": site_id,
                    "block_id": block_id,
                    "block_type": block_type,
                }
            )
            continue

        # Step 3: Check if off-grid is big enough
        if off_grid.area < min_useful_area:
            all_parts.append(
                {
                    "geometry": block,
                    "class": "block",
                    "part_type": "block",
                    "site": site_id,
                    "block_id": block_id,
                    "block_type": block_type,
                }
            )
            continue

        # Step 4: Split into parts
        try:
            corner_parts, side_parts, off_grid_final = create_block_parts_from_off_grid(
                block, off_grid, angle_threshold, corner_distance
            )

            # Add corner parts
            for i, corner in enumerate(corner_parts):
                all_parts.append(
                    {
                        "geometry": corner,
                        "class": "part",
                        "part_type": "corner",
                        "part_index": i,
                        "site": site_id,
                        "block_id": block_id,
                        "block_type": block_type,
                    }
                )

            # Add side parts
            for i, side in enumerate(side_parts):
                all_parts.append(
                    {
                        "geometry": side,
                        "class": "part",
                        "part_type": "side",
                        "part_index": i,
                        "site": site_id,
                        "block_id": block_id,
                        "block_type": block_type,
                    }
                )

            # Add off-grid center
            all_parts.append(
                {
                    "geometry": off_grid_final,
                    "class": "part",
                    "part_type": "off_grid",
                    "site": site_id,
                    "block_id": block_id,
                    "block_type": block_type,
                }
            )

            print(
                f"  Block {block_id}: {len(corner_parts)} "
                f"corners, {len(side_parts)} sides, 1 off-grid"
            )

        except Exception as e:
            print(f"  Warning: Failed to create parts for block {block_id}: {e}")
            # Keep block as-is
            all_parts.append(
                {
                    "geometry": block,
                    "class": "block",
                    "part_type": "block",
                    "site": site_id,
                    "block_id": block_id,
                    "block_type": block_type,
                }
            )

    # Convert to GeoDataFrame
    if not all_parts:
        return gpd.GeoDataFrame(
            columns=["geometry", "class", "part_type", "site", "block_id", "block_type"]
        )

    parts_gdf = gpd.GeoDataFrame(all_parts, crs=blocks.crs)

    print(f"\nGenerated {len(parts_gdf)} total parts:")
    print(f"  Corner parts: {len(parts_gdf[parts_gdf['part_type'] == 'corner'])}")
    print(f"  Side parts: {len(parts_gdf[parts_gdf['part_type'] == 'side'])}")
    print(f"  Off-grid parts: {len(parts_gdf[parts_gdf['part_type'] == 'off_grid'])}")
    print(f"  Blocks (no subdivision): {len(parts_gdf[parts_gdf['part_type'] == 'block'])}")

    return parts_gdf
