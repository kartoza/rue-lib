#!/usr/bin/env python3
"""Test script for fix_area_too_large function."""

from shapely.geometry import Polygon

from src.rue_lib.streets.cell import fix_area_too_large
from src.rue_lib.streets.grids import is_good_cell


def test_fix_oversized_cell():
    """Test splitting an oversized cell into good cells."""
    print("=" * 80)
    print("TEST: Fix Oversized Cell")
    print("=" * 80)

    # Create an oversized cell (200m x 100m = 20,000 m²)
    # Target area is 10,000 m² (100m x 100m)
    oversized_cell = Polygon([[0, 0, 0], [200, 0, 0], [200, 100, 0], [0, 100, 0], [0, 0, 0]])

    target_area = 100.0 * 100.0  # 10,000 m²

    # Check quality of original cell
    quality = is_good_cell(oversized_cell, target_area)

    print("\nOriginal cell:")
    print(f"  Area: {oversized_cell.area:.1f} m²")
    print(f"  Target area: {target_area:.1f} m²")
    print(f"  Area ratio: {quality['area_ratio']:.3f}")
    print(f"  Is good: {quality['is_good']}")
    print(f"  Reason: {quality['reason']}")

    # Try to fix it
    print("\nAttempting to split oversized cell...")
    fixed_cells = fix_area_too_large(oversized_cell, target_area, quality)

    print(f"\nResult: {len(fixed_cells)} cell(s) after split")

    # Check quality of each resulting cell
    good_count = 0
    for i, cell in enumerate(fixed_cells):
        cell_quality = is_good_cell(cell, target_area)
        status = "✓ GOOD" if cell_quality["is_good"] else "✗ BAD"

        print(f"\n  Cell {i}: {status}")
        print(f"    Area: {cell.area:.1f} m²")
        print(f"    Area ratio: {cell_quality['area_ratio']:.3f}")
        print(f"    Vertices: {cell_quality['num_vertices']}")
        print(f"    Right angles: {cell_quality['right_angles']}")
        print(f"    Reason: {cell_quality['reason']}")

        if cell_quality["is_good"]:
            good_count += 1

    print(f"\n{'=' * 80}")
    print(f"Summary: {good_count}/{len(fixed_cells)} cells are good")

    if good_count > 0 and len(fixed_cells) > 1:
        print("✓ SUCCESS: Oversized cell was successfully split into good cells")
    elif len(fixed_cells) == 1:
        print("✗ UNABLE TO FIX: Cell could not be split")
    else:
        print("✗ PARTIAL: Split created some bad cells")

    print("=" * 80)


def test_fix_very_oversized_cell():
    """Test splitting a very oversized cell (3x target area)."""
    print("\n" + "=" * 80)
    print("TEST: Fix Very Oversized Cell (3x target)")
    print("=" * 80)

    # Create a very oversized cell (300m x 100m = 30,000 m²)
    # Target area is 10,000 m² (100m x 100m)
    oversized_cell = Polygon([[0, 0, 0], [300, 0, 0], [300, 100, 0], [0, 100, 0], [0, 0, 0]])

    target_area = 100.0 * 100.0  # 10,000 m²

    quality = is_good_cell(oversized_cell, target_area)

    print("\nOriginal cell:")
    print(f"  Area: {oversized_cell.area:.1f} m²")
    print(f"  Target area: {target_area:.1f} m²")
    print(f"  Area ratio: {quality['area_ratio']:.3f}")

    # Try to fix it
    print("\nAttempting to split very oversized cell...")
    fixed_cells = fix_area_too_large(oversized_cell, target_area, quality)

    print(f"\nResult: {len(fixed_cells)} cell(s) after split")

    good_count = 0
    for i, cell in enumerate(fixed_cells):
        cell_quality = is_good_cell(cell, target_area)
        status = "✓ GOOD" if cell_quality["is_good"] else "✗ BAD"

        print(
            f"  Cell {i}: {status} | area={cell.area:.1f}m² | "
            f"ratio={cell_quality['area_ratio']:.3f} | {cell_quality['reason']}"
        )

        if cell_quality["is_good"]:
            good_count += 1

    print(f"\n{'=' * 80}")
    print(f"Summary: {good_count}/{len(fixed_cells)} cells are good")
    print("=" * 80)


def test_fix_non_oversized_cell():
    """Test that non-oversized cells are returned unchanged."""
    print("\n" + "=" * 80)
    print("TEST: Non-Oversized Cell (should return unchanged)")
    print("=" * 80)

    # Create a good cell (100m x 100m = 10,000 m²)
    good_cell = Polygon([[0, 0, 0], [100, 0, 0], [100, 100, 0], [0, 100, 0], [0, 0, 0]])

    target_area = 100.0 * 100.0  # 10,000 m²

    quality = is_good_cell(good_cell, target_area)

    print("\nOriginal cell:")
    print(f"  Is good: {quality['is_good']}")
    print(f"  Reason: {quality['reason']}")

    # Try to "fix" it
    result = fix_area_too_large(good_cell, target_area, quality)

    print(f"\nResult: {len(result)} cell(s)")

    if len(result) == 1 and result[0] == good_cell:
        print("✓ SUCCESS: Good cell returned unchanged")
    else:
        print("✗ FAILED: Good cell was modified")

    print("=" * 80)


if __name__ == "__main__":
    test_fix_oversized_cell()
    test_fix_very_oversized_cell()
    test_fix_non_oversized_cell()

    print("\n" + "=" * 80)
    print("ALL TESTS COMPLETED")
    print("=" * 80)
