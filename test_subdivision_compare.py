#!/usr/bin/env python3
"""Compare different subdivision approaches."""

import sys
from pathlib import Path

# Add src to path for direct imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

import geopandas as gpd
from rue_lib.cluster.off_grid_subdivision_simple import subdivide_off_grid_simple
from rue_lib.cluster.off_grid_subdivision_v2 import (
    subdivide_off_grid_oriented,
    subdivide_with_inward_buffer,
)
from shapely.geometry import Polygon


def test_subdivision_methods():
    """Compare different subdivision methods on a simple rectangle."""
    print("=" * 80)
    print("COMPARING SUBDIVISION METHODS")
    print("=" * 80)

    # Create a simple rectangular off-grid area
    off_grid = Polygon([(0, 0), (150, 0), (150, 100), (0, 100), (0, 0)])

    print("\nTest polygon: 150m x 100m rectangle")
    print(f"Area: {off_grid.area:.2f} m²")
    print(f"Target plot size: 30m x 35m = {30 * 35} m²")
    print()

    # Test parameters
    part_og_w = 30.0
    part_og_d = 35.0
    min_area = part_og_w * part_og_d * 0.3

    # Method 1: Simple axis-aligned grid
    print("-" * 80)
    print("Method 1: Axis-Aligned Grid")
    print("-" * 80)
    plots1 = subdivide_off_grid_simple(
        off_grid,
        part_og_w=part_og_w,
        part_og_d=part_og_d,
        rotation_angle=0.0,
        min_plot_area=min_area,
    )
    total_area1 = sum(p.area for p in plots1)
    coverage1 = (total_area1 / off_grid.area) * 100
    print(f"Plots created: {len(plots1)}")
    print(f"Total area: {total_area1:.2f} m²")
    print(f"Coverage: {coverage1:.1f}%")
    print(f"Average plot size: {total_area1 / len(plots1):.2f} m²")

    # Method 2: Oriented grid
    print()
    print("-" * 80)
    print("Method 2: Oriented Grid (auto-detect orientation)")
    print("-" * 80)
    plots2 = subdivide_off_grid_oriented(
        off_grid,
        part_og_w=part_og_w,
        part_og_d=part_og_d,
        auto_orient=True,
        min_plot_area=min_area,
    )
    total_area2 = sum(p.area for p in plots2)
    coverage2 = (total_area2 / off_grid.area) * 100
    print(f"Plots created: {len(plots2)}")
    print(f"Total area: {total_area2:.2f} m²")
    print(f"Coverage: {coverage2:.1f}%")
    print(f"Average plot size: {total_area2 / len(plots2):.2f} m²")

    # Method 3: Inward buffering
    print()
    print("-" * 80)
    print("Method 3: Inward Buffering")
    print("-" * 80)
    plots3 = subdivide_with_inward_buffer(
        off_grid,
        part_og_w=part_og_w,
        part_og_d=part_og_d,
        min_plot_area=min_area,
    )
    total_area3 = sum(p.area for p in plots3)
    coverage3 = (total_area3 / off_grid.area) * 100
    print(f"Plots created: {len(plots3)}")
    print(f"Total area: {total_area3:.2f} m²")
    print(f"Coverage: {coverage3:.1f}%")
    print(f"Average plot size: {total_area3 / len(plots3):.2f} m²")

    # Save outputs for visualization
    output_dir = Path("outputs/subdivision_compare")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save each method's output
    for method_num, plots in enumerate([plots1, plots2, plots3], 1):
        plots_data = [
            {"geometry": plot, "plot_id": i, "area": plot.area, "method": method_num}
            for i, plot in enumerate(plots)
        ]
        gdf = gpd.GeoDataFrame(plots_data, crs="EPSG:32633")
        gdf.to_file(output_dir / f"method{method_num}_plots.geojson", driver="GeoJSON")

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\nBest coverage: Method 2 (Oriented) with {coverage2:.1f}%")
    print(f"Most plots: Method 1 (Axis-Aligned) with {len(plots1)} plots")
    print(f"\nOutputs saved to: {output_dir}")


if __name__ == "__main__":
    test_subdivision_methods()
