#!/usr/bin/env python3
"""Standalone test script for inspect_and_fix_cells with embedded sample data."""

import json
from pathlib import Path

from pyproj import Transformer
from shapely.geometry import shape
from shapely.ops import transform

# Import from the library
from src.rue_lib.streets.cell import Cell, inspect_and_fix_cells
from src.rue_lib.streets.grids import is_good_cell

# Create transformer from WGS84 (EPSG:4326) to UTM Zone 51N (EPSG:32651)
# Philippines region around 123°E, 10°N
transformer = Transformer.from_crs("EPSG:4326", "EPSG:32651", always_xy=True)


def create_sample_grid_geojson():
    geojson = {
        "type": "FeatureCollection",
        "name": "all_grids_merged",
        "features": [
            {
                "type": "Feature",
                "id": 29,
                "properties": {"id": "1", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0566845, 10.8751022, 0.0],
                            [123.055952, 10.8742709, 0.0],
                            [123.0559947, 10.8736584, 0.0],
                            [123.0562152, 10.8740303, 0.0],
                            [123.0562152, 10.8740303, 0.0],
                            [123.0562152, 10.8740303, 0.0],
                            [123.0566984, 10.8742492, 0.0],
                            [123.0571459, 10.8742493, 0.0],
                            [123.057146, 10.8742493, 0.0],
                            [123.0573147, 10.8743676, 0.0],
                            [123.0574063, 10.874432, 0.0],
                            [123.0570693, 10.8746148, 0.0],
                            [123.0567117, 10.8748764, 0.0],
                            [123.0566187, 10.8749445, 0.0],
                            [123.0566187, 10.8749445, 0.0],
                            [123.0566187, 10.8749445, 0.0],
                            [123.056903, 10.8750952, 0.0],
                            [123.0568951, 10.8750954, 0.0],
                            [123.0566845, 10.8751022, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 30,
                "properties": {"id": "2", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0548142, 10.8728316, 0.0],
                            [123.0548594, 10.8728128, 0.0],
                            [123.0547791, 10.8728819, 0.0],
                            [123.0548142, 10.8728316, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 31,
                "properties": {"id": "3", "grid_type": "off_grid", "is_good": False},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0545159, 10.8731085, 0.0],
                            [123.0540601, 10.8735009, 0.0],
                            [123.0536918, 10.873083, 0.0],
                            [123.0542806, 10.8726948, 0.0],
                            [123.054219, 10.8729039, 0.0],
                            [123.0545159, 10.8731085, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 32,
                "properties": {"id": "4", "grid_type": "off_grid", "is_good": False},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0546611, 10.8741829, 0.0],
                            [123.0540601, 10.8735009, 0.0],
                            [123.0545159, 10.8731085, 0.0],
                            [123.0545869, 10.8731575, 0.0],
                            [123.0547791, 10.8728819, 0.0],
                            [123.0548594, 10.8728128, 0.0],
                            [123.0549563, 10.8727726, 0.0],
                            [123.0550755, 10.8727231, 0.0],
                            [123.0551121, 10.8726984, 0.0],
                            [123.0554496, 10.8724703, 0.0],
                            [123.0553721, 10.8727593, 0.0],
                            [123.055368, 10.8731213, 0.0],
                            [123.0556275, 10.8731939, 0.0],
                            [123.0557125, 10.8732777, 0.0],
                            [123.0546611, 10.8741829, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 33,
                "properties": {"id": "5", "grid_type": "off_grid", "is_good": False},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0559947, 10.8736584, 0.0],
                            [123.055952, 10.8742709, 0.0],
                            [123.055262, 10.8748649, 0.0],
                            [123.0546611, 10.8741829, 0.0],
                            [123.0557125, 10.8732777, 0.0],
                            [123.0558484, 10.8734116, 0.0],
                            [123.0559947, 10.8736584, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 34,
                "properties": {"id": "6", "grid_type": "off_grid", "is_good": False},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.055262, 10.8748649, 0.0],
                            [123.055952, 10.8742709, 0.0],
                            [123.0566845, 10.8751022, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.0563344, 10.8751135, 0.0],
                            [123.055862, 10.8755457, 0.0],
                            [123.055262, 10.8748649, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 35,
                "properties": {"id": "7", "grid_type": "off_grid", "is_good": False},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0540601, 10.8735009, 0.0],
                            [123.0533701, 10.8740949, 0.0],
                            [123.0529243, 10.873589, 0.0],
                            [123.0536918, 10.873083, 0.0],
                            [123.0540601, 10.8735009, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 36,
                "properties": {"id": "8", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0539711, 10.8747769, 0.0],
                            [123.0533701, 10.8740949, 0.0],
                            [123.0540601, 10.8735009, 0.0],
                            [123.0546611, 10.8741829, 0.0],
                            [123.0539711, 10.8747769, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 37,
                "properties": {"id": "9", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.055262, 10.8748649, 0.0],
                            [123.054572, 10.8754589, 0.0],
                            [123.0539711, 10.8747769, 0.0],
                            [123.0546611, 10.8741829, 0.0],
                            [123.055262, 10.8748649, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 38,
                "properties": {"id": "10", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.054572, 10.8754589, 0.0],
                            [123.055262, 10.8748649, 0.0],
                            [123.055862, 10.8755457, 0.0],
                            [123.0551902, 10.8761604, 0.0],
                            [123.054572, 10.8754589, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 39,
                "properties": {"id": "11", "grid_type": "off_grid", "is_good": False},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0533701, 10.8740949, 0.0],
                            [123.0526801, 10.8746889, 0.0],
                            [123.0521567, 10.8740949, 0.0],
                            [123.0529243, 10.873589, 0.0],
                            [123.0533701, 10.8740949, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 40,
                "properties": {"id": "12", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0526801, 10.8746889, 0.0],
                            [123.0533701, 10.8740949, 0.0],
                            [123.0539711, 10.8747769, 0.0],
                            [123.0532811, 10.8753709, 0.0],
                            [123.0526801, 10.8746889, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 41,
                "properties": {"id": "13", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0532811, 10.8753709, 0.0],
                            [123.0539711, 10.8747769, 0.0],
                            [123.054572, 10.8754589, 0.0],
                            [123.0538821, 10.8760529, 0.0],
                            [123.0532811, 10.8753709, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 42,
                "properties": {"id": "14", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0538821, 10.8760529, 0.0],
                            [123.054572, 10.8754589, 0.0],
                            [123.0551902, 10.8761604, 0.0],
                            [123.0545184, 10.876775, 0.0],
                            [123.0538821, 10.8760529, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 43,
                "properties": {"id": "15", "grid_type": "off_grid", "is_good": False},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0521567, 10.8740949, 0.0],
                            [123.0526801, 10.8746889, 0.0],
                            [123.0519902, 10.8752829, 0.0],
                            [123.0513892, 10.8746009, 0.0],
                            [123.0521567, 10.8740949, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 44,
                "properties": {"id": "16", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0526801, 10.8746889, 0.0],
                            [123.0532811, 10.8753709, 0.0],
                            [123.0525911, 10.8759649, 0.0],
                            [123.0523493, 10.8756905, 0.0],
                            [123.0519902, 10.8752829, 0.0],
                            [123.0526801, 10.8746889, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 45,
                "properties": {"id": "17", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0532811, 10.8753709, 0.0],
                            [123.0538821, 10.8760529, 0.0],
                            [123.0531921, 10.8766469, 0.0],
                            [123.0526441, 10.8760251, 0.0],
                            [123.0525911, 10.8759649, 0.0],
                            [123.0532811, 10.8753709, 0.0],
                        ]
                    ],
                },
            },
            {
                "type": "Feature",
                "id": 46,
                "properties": {"id": "18", "grid_type": "off_grid", "is_good": True},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [
                        [
                            [123.0531921, 10.8766469, 0.0],
                            [123.0538821, 10.8760529, 0.0],
                            [123.0545184, 10.876775, 0.0],
                            [123.0538465, 10.8773896, 0.0],
                            [123.0537015, 10.8772251, 0.0],
                            [123.0531921, 10.8766469, 0.0],
                        ]
                    ],
                },
            },
        ],
    }

    return geojson


def test_single_bad_cell():
    """Test with real-world geographic data from Philippines."""
    print("=" * 80)
    print("TEST: Real-world Grid Data (Philippines)")
    print("=" * 80)

    # Create sample grid with geographic coordinates
    grid_data = create_sample_grid_geojson()

    all_cells = []

    # Save to file
    output_path = Path("outputs/sample_grid_philippines.geojson")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(grid_data, f, indent=2)
    print(f"\n✓ Saved grid to: {output_path}")

    # Load geometries, reproject to UTM, and compute quality
    preferred_width = 100
    preferred_height = 100
    target_area = preferred_width * preferred_height

    cells_utm = []
    qualities = []
    areas_m2 = []

    print("\nReprojecting from EPSG:4326 (lat/lon) to EPSG:32651 (UTM Zone 51N)...")

    for _i, feature in enumerate(grid_data["features"]):
        geom_latlon = shape(feature["geometry"])

        geom_utm = transform(transformer.transform, geom_latlon)
        cells_utm.append(geom_utm)

        area_m2 = geom_utm.area
        areas_m2.append(area_m2)

    for i, (geom_utm, _area_m2) in enumerate(zip(cells_utm, areas_m2)):
        quality = is_good_cell(geom_utm, target_area)
        qualities.append(quality)

        _cell = Cell(i, geom_utm, quality)
        all_cells.append(_cell)

    # Print quality information
    print("\nCell Quality Analysis:")
    print("-" * 80)
    for i, (quality, area_m2) in enumerate(zip(qualities, areas_m2)):
        feature_id = grid_data["features"][i]["id"]
        expected_good = grid_data["features"][i]["properties"]["is_good"]
        status = "✓ GOOD" if quality["is_good"] else "✗ BAD"
        expected = "GOOD" if expected_good else "BAD"
        match = "✓" if (quality["is_good"] == expected_good) else "✗"

        print(
            f"{match} Cell {i:2d} (id={feature_id}): {status:6s} (expected {expected:4s}) | "
            f"area={area_m2:8.1f}m² | "
            f"ratio={quality['area_ratio']:.3f} | "
            f"v={quality['num_vertices']} | "
            f"ra={quality['right_angles']} | "
            f"{quality['reason']}"
        )

    # Run inspection with fixing enabled
    print("\n" + "-" * 80)
    print("Running inspect_and_fix_cells (with auto-fix enabled)...")
    print("-" * 80)

    fixed_cells = inspect_and_fix_cells(all_cells, target_area, fix_oversized=True)

    # Summary - compare before and after
    bad_count_before = sum(1 for q in qualities if not q["is_good"])
    bad_count_after = sum(1 for q in fixed_cells if not q.quality["is_good"])

    print("\n✓ Test completed:")
    print(f"  Before fixing: {bad_count_before}/{len(qualities)} bad cells")
    print(f"  After fixing:  {bad_count_after}/{len(fixed_cells)} bad cells")
    print(f"  Cells added:   {len(fixed_cells) - len(qualities)} (from splitting)")
    print()

    # Save fixed cells back to GeoJSON
    print("Saving fixed cells to GeoJSON...")

    # Create inverse transformer to convert back from UTM to lat/lon
    inverse_transformer = Transformer.from_crs("EPSG:32651", "EPSG:4326", always_xy=True)

    # Build GeoJSON features for fixed cells
    fixed_features = []
    for i, cell in enumerate(fixed_cells):
        cell_utm = cell.geom
        quality = cell.quality
        # Reproject back to lat/lon
        cell_latlon = transform(inverse_transformer.transform, cell_utm)

        # Create feature
        feature = {
            "type": "Feature",
            "id": i,
            "properties": {
                "id": str(i),
                "is_good": quality["is_good"],
                "reason": quality["reason"],
                "area": cell_utm.area,
                "area_ratio": quality["area_ratio"],
                "num_vertices": quality["num_vertices"],
                "right_angles": quality["right_angles"],
            },
            "geometry": {
                "type": "Polygon",
                "coordinates": [
                    [[coord[0], coord[1], 0.0] for coord in cell_latlon.exterior.coords]
                ],
            },
        }
        fixed_features.append(feature)

    # Create GeoJSON FeatureCollection
    fixed_geojson = {"type": "FeatureCollection", "name": "fixed_grids", "features": fixed_features}

    # Save to file
    fixed_output_path = Path("outputs/sample_grid_philippines_fixed.geo.json")
    with open(fixed_output_path, "w") as f:
        json.dump(fixed_geojson, f, indent=2)

    print(f"✓ Saved fixed grid to: {fixed_output_path}")


if __name__ == "__main__":
    print("\n" + "=" * 80)
    print(" CELL INSPECTION AND FIX TESTS")
    print("=" * 80 + "\n")

    # Run both tests
    test_single_bad_cell()

    print("=" * 80)
    print("TEST COMPLETED")
    print("=" * 80)
    print("\nGenerated GeoJSON files:")
    print("  - outputs/sample_grid_philippines.geojson (original)")
    print("  - outputs/sample_grid_philippines_fixed.geojson (after fixing)")
