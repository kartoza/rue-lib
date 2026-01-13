#!/usr/bin/env python3
"""Test GDAL buffering implementation."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent / "src"))

import geopandas as gpd
from shapely.geometry import LineString

from rue_lib.site.roads import buffer_geometry_gdal, buffer_roads


def test_buffer_geometry_gdal():
    """Test the GDAL buffer function directly."""
    print("=" * 80)
    print("TEST 1: Testing buffer_geometry_gdal()")
    print("=" * 80)

    # Create a simple line
    line = LineString([(0, 0), (100, 0)])
    print(f"\nInput line: {line.wkt}")
    print(f"Length: {line.length:.2f}")

    # Buffer with GDAL
    buffered_wkt = buffer_geometry_gdal(line.wkt, distance=10.0, cap_style=2)
    print("\nBuffered (GDAL, flat caps, distance=10):")
    print(f"WKT length: {len(buffered_wkt)} chars")

    from shapely import wkt

    buffered_geom = wkt.loads(buffered_wkt)
    print(f"Buffered area: {buffered_geom.area:.2f}")
    # length * 2 * radius for flat caps
    print(f"Expected area (approx): {line.length * 2 * 10:.2f}")

    # Test with round caps
    buffered_round_wkt = buffer_geometry_gdal(line.wkt, distance=10.0, cap_style=1)
    buffered_round = wkt.loads(buffered_round_wkt)
    print("\nBuffered (GDAL, round caps, distance=10):")
    print(f"Buffered area: {buffered_round.area:.2f}")
    # includes end caps
    print(f"Expected area (approx): {line.length * 2 * 10 + 3.14159 * 10**2:.2f}")

    # Test with invalid geometry
    print("\nTesting with invalid WKT:")
    result = buffer_geometry_gdal("INVALID WKT", distance=10.0)
    print(f"Result: {result} (should return original)")

    print("\n✓ Test 1 completed\n")


def test_buffer_roads():
    """Test buffer_roads with real road data."""
    print("=" * 80)
    print("TEST 2: Testing buffer_roads() with GDAL")
    print("=" * 80)

    # Create sample roads GeoDataFrame
    roads = gpd.GeoDataFrame(
        {
            "type": ["road_art", "road_sec", "road_sec_new"],
            "geometry": [
                LineString([(0, 0), (100, 0)]),
                LineString([(0, 50), (100, 50)]),
                LineString([(0, 100), (100, 100)]),
            ],
        },
        crs="EPSG:32633",  # UTM zone 33N (metric)
    )

    print(f"\nInput roads: {len(roads)} roads")
    for _idx, row in roads.iterrows():
        print(f"  {row['type']}: length={row.geometry.length:.2f}m")

    # Buffer roads
    art_width = 16.0
    sec_width = 10.0
    buffered = buffer_roads(
        roads, road_arterial_width_m=art_width, road_secondary_width_m=sec_width
    )

    print(f"\nBuffered roads: {len(buffered)} polygons")
    for _idx, row in buffered.iterrows():
        print(f"  {row['type']}: area={row.geometry.area:.2f}m²")

    # Check arterial road
    art = buffered[buffered["type"] == "road_art"].iloc[0]
    art_expected_area = 100 * art_width  # length * width for flat caps
    print("\nArterial road check:")
    print(f"  Expected area: ~{art_expected_area:.2f}m²")
    print(f"  Actual area: {art.geometry.area:.2f}m²")
    print(f"  Difference: {abs(art.geometry.area - art_expected_area):.2f}m²")

    # Check secondary road
    sec = buffered[buffered["type"] == "road_sec"].iloc[0]
    sec_expected_area = 100 * sec_width
    print("\nSecondary road check:")
    print(f"  Expected area: ~{sec_expected_area:.2f}m²")
    print(f"  Actual area: {sec.geometry.area:.2f}m²")
    print(f"  Difference: {abs(sec.geometry.area - sec_expected_area):.2f}m²")

    print("\n✓ Test 2 completed\n")


def test_with_real_data():
    """Test with real block and road data if available."""
    print("=" * 80)
    print("TEST 3: Testing with real data (if available)")
    print("=" * 80)

    block_file = Path("block_test.geojson")
    road_file = Path("road_test.geojson")

    if not block_file.exists() or not road_file.exists():
        print("\n⚠ Real data files not found, skipping test")
        print(f"  Looking for: {block_file}, {road_file}")
        return

    print("\nLoading real data...")
    roads = gpd.read_file(road_file)
    print(f"  Loaded {len(roads)} roads")

    if "road_type" in roads.columns:
        # Rename column to match expected format
        roads["type"] = roads["road_type"]

    # Buffer roads
    buffered = buffer_roads(roads, road_arterial_width_m=16.0, road_secondary_width_m=10.0)

    print("\nBuffered results:")
    print(f"  Total buffered roads: {len(buffered)}")
    print(f"  Total area: {buffered.geometry.area.sum():.2f}m²")

    if not buffered.empty:
        print(f"  Average area: {buffered.geometry.area.mean():.2f}m²")
        print(f"  Min area: {buffered.geometry.area.min():.2f}m²")
        print(f"  Max area: {buffered.geometry.area.max():.2f}m²")

    print("\n✓ Test 3 completed\n")


if __name__ == "__main__":
    print("\n")
    test_buffer_geometry_gdal()
    test_buffer_roads()
    test_with_real_data()
    print("=" * 80)
    print("ALL TESTS COMPLETED")
    print("=" * 80)
    print("\n")
