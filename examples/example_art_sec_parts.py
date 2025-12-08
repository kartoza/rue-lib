#!/usr/bin/env python3
"""Example: Generate arterial/secondary block parts without off-grid subdivision.

This example shows how to:
1. Load blocks and roads
2. Extract block edges and assign road types
3. Generate parts using the art_sec_no_offgrid module
"""

from pathlib import Path

import geopandas as gpd

from rue_lib.cluster.art_sec_no_offgrid import generate_art_sec_parts_no_offgrid
from rue_lib.cluster.block_edges import extract_block_edges


def main():
    """Run example for arterial/secondary parts generation."""
    # Input files (adjust paths as needed)
    blocks_path = "outputs/step2_streets/all_grids_merged.geojson"
    roads_path = "outputs/real_data_test/input_roads.geojson"
    output_path = "outputs/art_sec_parts.geojson"

    print("Loading input data...")
    blocks_gdf = gpd.read_file(blocks_path)
    roads_gdf = gpd.read_file(roads_path)

    print(f"  Loaded {len(blocks_gdf)} blocks")
    print(f"  Loaded {len(roads_gdf)} roads")

    # Add block_id if not present
    if 'block_id' not in blocks_gdf.columns:
        blocks_gdf['block_id'] = range(len(blocks_gdf))

    print("\nExtracting block edges and assigning road types...")
    # Extract edges from blocks and assign road types based on proximity
    block_edges_gdf = extract_block_edges(
        blocks_gdf=blocks_gdf,
        roads_gdf=roads_gdf,
        tolerance=5.0  # 5 meters tolerance for matching edges to roads
    )

    print(f"  Extracted {len(block_edges_gdf)} edges")

    # Show road type distribution
    if 'road_type' in block_edges_gdf.columns:
        road_type_counts = block_edges_gdf['road_type'].value_counts()
        print("\n  Edge road type distribution:")
        for road_type, count in road_type_counts.items():
            print(f"    {road_type}: {count}")

    print("\nGenerating arterial/secondary parts...")
    # Generate parts
    parts_gdf = generate_art_sec_parts_no_offgrid(
        blocks_gdf=blocks_gdf,
        block_edges_gdf=block_edges_gdf,
        part_art_d=40.0,  # Arterial road offset depth (meters)
        part_sec_d=30.0,  # Secondary road offset depth (meters)
        part_loc_d=20.0,  # Local road offset depth (meters)
    )

    print(f"  Generated {len(parts_gdf)} parts")

    # Show part type distribution
    if 'type' in parts_gdf.columns:
        part_type_counts = parts_gdf['type'].value_counts()
        print("\n  Part type distribution:")
        for part_type, count in part_type_counts.items():
            print(f"    {part_type}: {count}")

    # Save output
    print(f"\nSaving results to {output_path}...")
    output_dir = Path(output_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    parts_gdf.to_file(output_path, driver='GeoJSON')

    print("Done!")


if __name__ == '__main__':
    main()