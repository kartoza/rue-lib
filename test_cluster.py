#!/usr/bin/env python3
"""Test the parts generation workflow with real block and road data."""

import sys
from pathlib import Path

# Add src to path for direct imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

import geopandas as gpd

from rue_lib.cluster import (
    classify_plot_by_area,
    create_block_parts_from_off_grid,
    create_off_grid_area,
    subdivide_off_grid,
)


def test_with_real_data():
    """Test with the provided block_test.geojson and road_test.geojson."""
    print("=" * 80)
    print("TESTING WITH REAL BLOCK AND ROAD DATA")
    print("=" * 80)

    # Load data
    print("\nStep 1: Loading data...")
    blocks = gpd.read_file("block_test.geojson")
    roads = gpd.read_file("road_test.geojson")

    print(f"  Loaded {len(blocks)} blocks")
    print(f"  Loaded {len(roads)} roads")
    print(f"  CRS: {blocks.crs}")

    # Check road types
    if "road_type" in roads.columns:
        print("\n  Road types:")
        for rt in roads["road_type"].unique():
            count = len(roads[roads["road_type"] == rt])
            print(f"    {rt}: {count} road(s)")
    else:
        print("  ⚠ Warning: No 'road_type' column in roads!")

    # Check blocks
    print("\n  Block info:")
    print(f"    Total area: {blocks.geometry.area.sum():.2f} m²")
    print(f"    Average area: {blocks.geometry.area.mean():.2f} m²")
    print(f"    Bounds: {blocks.total_bounds}")

    # Test off-grid creation for ALL blocks
    print("\n" + "=" * 80)
    print("Step 2: Testing off-grid creation on ALL blocks...")
    print("=" * 80)

    # Test configuration
    config = {"part_art_d": 10.0, "part_sec_d": 10.0, "part_loc_d": 15.0}
    print(
        f"\nUsing config: art={config['part_art_d']}m, sec={config['part_sec_d']}m, "
        f"loc={config['part_loc_d']}m\n"
    )

    output_dir = Path("outputs/real_data_test")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Test each block
    off_grids = []
    for idx, block_row in blocks.iterrows():
        block = block_row.geometry
        block_id = block_row.get("id", idx)

        print(f"Block {block_id}:")
        print(f"  Area: {block.area:.2f} m²")
        print(f"  Bounds: {block.bounds}")

        off_grid = create_off_grid_area(
            block,
            roads,
            part_art_d=config["part_art_d"],
            part_sec_d=config["part_sec_d"],
            part_loc_d=config["part_loc_d"],
        )

        if off_grid:
            reduction = ((block.area - off_grid.area) / block.area) * 100
            print(f"  ✓ Off-grid area: {off_grid.area:.2f} m² ({reduction:.1f}% reduction)")
            off_grids.append(
                {
                    "geometry": off_grid,
                    "block_id": block_id,
                    "original_area": block.area,
                    "off_grid_area": off_grid.area,
                    "reduction_pct": reduction,
                }
            )
        else:
            print("  ✗ No off-grid area created")
            off_grids.append(
                {
                    "geometry": block,  # Keep original
                    "block_id": block_id,
                    "original_area": block.area,
                    "off_grid_area": 0,
                    "reduction_pct": 0,
                }
            )
        print()

    # Create GeoDataFrame with results
    off_grids_gdf = gpd.GeoDataFrame(off_grids, crs=blocks.crs)

    # Save off-grid results
    output_file = output_dir / "off_grid_results.geojson"
    off_grids_gdf.to_file(output_file, driver="GeoJSON")
    print(f"✓ Off-grid results saved to: {output_file}")

    # Create frames (blocks with holes) by subtracting off-grid from original blocks
    print("\n" + "=" * 80)
    print("Step 3: Creating frames (blocks - off_grid)...")
    print("=" * 80)

    frames = []
    for idx, block_row in blocks.iterrows():
        block = block_row.geometry
        block_id = block_row.get("id", idx)

        # Find corresponding off-grid
        off_grid_row = off_grids_gdf[off_grids_gdf["block_id"] == block_id].iloc[0]
        off_grid = off_grid_row.geometry

        # Create frame by subtracting off-grid from block
        try:
            frame = block.difference(off_grid)

            if not frame.is_empty:
                frame_area = frame.area
                print(f"Block {block_id}: Frame area = {frame_area:.2f} m²")
                frames.append(
                    {
                        "geometry": frame,
                        "block_id": block_id,
                        "frame_area": frame_area,
                        "original_area": block.area,
                        "off_grid_area": off_grid.area,
                    }
                )
            else:
                print(f"Block {block_id}: ✗ Empty frame (off-grid = entire block)")
        except Exception as e:
            print(f"Block {block_id}: ✗ Error creating frame: {e}")

    # Create GeoDataFrame with frames
    if frames:
        frames_gdf = gpd.GeoDataFrame(frames, crs=blocks.crs)

        # Save frames
        frames_file = output_dir / "frames.geojson"
        frames_gdf.to_file(frames_file, driver="GeoJSON")
        print(f"\n✓ Frames saved to: {frames_file}")
    else:
        print("\n✗ No frames created")

    # Also save input data
    blocks.to_file(output_dir / "input_blocks.geojson", driver="GeoJSON")
    roads.to_file(output_dir / "input_roads.geojson", driver="GeoJSON")

    # Step 4: Test block parts generation on frames
    print("\n" + "=" * 80)
    print("Step 4: Creating block parts (corners + sides) from frames...")
    print("=" * 80)

    # Parameters for parts generation
    angle_threshold = 155.0  # Degrees - detect sharp corners
    corner_distance = 50.0  # Meters - smaller to avoid overlaps

    print(
        f"\nUsing params: angle_threshold={angle_threshold}°, corner_distance={corner_distance}m\n"
    )

    all_corner_parts = []
    all_side_parts = []
    all_off_grid_parts = []

    for _idx, frame_row in frames_gdf.iterrows():
        block_id = frame_row["block_id"]
        frame = frame_row.geometry  # The frame (block - off_grid) - for reference only

        # Get original block (without hole)
        original_block = blocks[blocks.get("id", blocks.index) == block_id].iloc[0].geometry

        # Find corresponding off-grid
        off_grid_row = off_grids_gdf[off_grids_gdf["block_id"] == block_id].iloc[0]
        off_grid = off_grid_row.geometry

        print(f"Block {block_id}:")
        print(f"  Original area: {original_block.area:.2f} m²")
        print(f"  Frame area: {frame_row['frame_area']:.2f} m²")
        print(f"  Off-grid area: {off_grid.area:.2f} m²")

        try:
            # Pass original block (not frame) - the function creates parts from block and off_grid
            corner_parts, side_parts, off_grid_final = create_block_parts_from_off_grid(
                original_block,  # Original block boundary
                off_grid,  # Off-grid center polygon
                angle_threshold=angle_threshold,
                corner_distance=corner_distance,
            )

            print(f"  ✓ Created {len(corner_parts)} corners, {len(side_parts)} sides")

            # Collect corner parts
            for i, corner in enumerate(corner_parts):
                all_corner_parts.append(
                    {
                        "geometry": corner,
                        "block_id": block_id,
                        "part_type": "corner",
                        "part_index": i,
                        "area": corner.area,
                    }
                )

            # Collect side parts
            for i, side in enumerate(side_parts):
                all_side_parts.append(
                    {
                        "geometry": side,
                        "block_id": block_id,
                        "part_type": "side",
                        "part_index": i,
                        "area": side.area,
                    }
                )

            # Collect off-grid part
            all_off_grid_parts.append(
                {
                    "geometry": off_grid_final,
                    "block_id": block_id,
                    "part_type": "off_grid",
                    "area": off_grid_final.area,
                }
            )

        except Exception as e:
            print(f"  ✗ Error creating parts: {e}")

        print()

    # Save parts to separate files
    if all_corner_parts:
        corners_gdf = gpd.GeoDataFrame(all_corner_parts, crs=blocks.crs)
        corners_file = output_dir / "corner_parts.geojson"
        corners_gdf.to_file(corners_file, driver="GeoJSON")
        print(f"✓ Corner parts saved to: {corners_file}")

    if all_side_parts:
        sides_gdf = gpd.GeoDataFrame(all_side_parts, crs=blocks.crs)
        sides_file = output_dir / "side_parts.geojson"
        sides_gdf.to_file(sides_file, driver="GeoJSON")
        print(f"✓ Side parts saved to: {sides_file}")

    if all_off_grid_parts:
        off_grid_parts_gdf = gpd.GeoDataFrame(all_off_grid_parts, crs=blocks.crs)
        off_grid_parts_file = output_dir / "off_grid_parts.geojson"
        off_grid_parts_gdf.to_file(off_grid_parts_file, driver="GeoJSON")
        print(f"✓ Off-grid parts saved to: {off_grid_parts_file}")

    # Create combined parts file
    all_parts = all_corner_parts + all_side_parts + all_off_grid_parts
    if all_parts:
        all_parts_gdf = gpd.GeoDataFrame(all_parts, crs=blocks.crs)
        all_parts_file = output_dir / "all_parts.geojson"
        all_parts_gdf.to_file(all_parts_file, driver="GeoJSON")
        print(f"✓ All parts saved to: {all_parts_file}")

        # Calculate coverage
        total_parts_area = all_parts_gdf.geometry.area.sum()
        total_blocks_area = blocks.geometry.area.sum()
        coverage = (total_parts_area / total_blocks_area) * 100
        print(
            f"\nParts coverage: {total_parts_area:.2f} / {total_blocks_area:.2f} m²"
            f"({coverage:.1f}%)"
        )

    # Step 5: Test off-grid subdivision into plots
    print("\n" + "=" * 80)
    print("Step 5: Subdividing off-grid areas into plots...")
    print("=" * 80)

    # Parameters for subdivision
    part_og_w = 30.0  # Plot width in meters
    part_og_d = 35.0  # Plot depth in meters

    print(f"\nUsing params: part_og_w={part_og_w}m, part_og_d={part_og_d}m\n")

    all_plots = []

    for _idx, off_grid_part in enumerate(all_off_grid_parts):
        block_id = off_grid_part["block_id"]
        off_grid_geom = off_grid_part["geometry"]

        print(f"Block {block_id}:")
        print(f"  Off-grid area: {off_grid_geom.area:.2f} m²")

        try:
            # Subdivide the off-grid area using oriented approach
            plots = subdivide_off_grid(
                off_grid_geom,
                part_og_w=part_og_w,
                part_og_d=part_og_d,
                min_plot_area=part_og_w * part_og_d * 0.3,
            )

            print(f"  ✓ Created {len(plots)} plots")

            # Collect plots
            for i, plot in enumerate(plots):
                all_plots.append(
                    {
                        "geometry": plot,
                        "block_id": block_id,
                        "plot_type": classify_plot_by_area(plot, part_og_w, part_og_d),
                        "plot_index": i,
                        "area": plot.area,
                        "parent_part": "off_grid",
                    }
                )

        except Exception as e:
            print(f"  ✗ Error subdividing off-grid: {e}")

        print()

    # Save plots to file
    if all_plots:
        plots_gdf = gpd.GeoDataFrame(all_plots, crs=blocks.crs)
        plots_file = output_dir / "off_grid_plots.geojson"
        plots_gdf.to_file(plots_file, driver="GeoJSON")
        print(f"✓ Off-grid plots saved to: {plots_file}")

        # Calculate statistics
        total_plot_area = plots_gdf.geometry.area.sum()
        total_off_grid_area = sum(p["area"] for p in all_off_grid_parts)
        plot_coverage = (
            (total_plot_area / total_off_grid_area) * 100 if total_off_grid_area > 0 else 0
        )

        print(f"  Total plots: {len(plots_gdf)}")
        print(f"  Total plot area: {total_plot_area:.2f} m²")
        print(f"  Off-grid area coverage: {plot_coverage:.1f}%")
        print(f"  Average plot size: {plots_gdf.geometry.area.mean():.2f} m²")
    else:
        print("✗ No plots created")

    # Summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)

    successful = len([g for g in off_grids if g["off_grid_area"] > 0])
    print(f"\nOff-grid creation: {successful}/{len(blocks)} blocks")

    if successful > 0:
        avg_reduction = (
            sum(g["reduction_pct"] for g in off_grids if g["off_grid_area"] > 0) / successful
        )
        print(f"  Average reduction: {avg_reduction:.1f}%")

    if frames:
        print(f"\nFrames created: {len(frames)} blocks")
        total_frame_area = sum(f["frame_area"] for f in frames)
        total_block_area = sum(f["original_area"] for f in frames)
        frame_pct = (total_frame_area / total_block_area) * 100
        print(f"  Total frame area: {total_frame_area:.2f} m² ({frame_pct:.1f}% of original)")

    if all_parts:
        print("\nBlock parts created:")
        print(f"  Corner parts: {len(all_corner_parts)}")
        print(f"  Side parts: {len(all_side_parts)}")
        print(f"  Off-grid parts: {len(all_off_grid_parts)}")
        print(f"  Total parts: {len(all_parts)}")

    if all_plots:
        print("\nOff-grid plots created:")
        print(f"  Total plots: {len(all_plots)}")
        print(f"  Average plot area: {sum(p['area'] for p in all_plots) / len(all_plots):.2f} m²")

    print("\n" + "=" * 80)
    print("TEST COMPLETE")
    print("=" * 80)

    return off_grids_gdf, frames_gdf if frames else None


if __name__ == "__main__":
    print("\n")
    parts = test_with_real_data()
    print("\n")
