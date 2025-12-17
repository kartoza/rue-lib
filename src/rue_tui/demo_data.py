"""Generate demo data for RUE TUI if input files don't exist."""

import json
from pathlib import Path


def create_demo_site(output_path: str) -> None:
    """Create a demo site polygon."""
    # Simple rectangular site with some irregularity
    coordinates = [
        [
            [28.0, -25.7],  # SW corner
            [28.1, -25.7],  # SE corner
            [28.1, -25.6],  # NE corner
            [28.05, -25.55],  # North indent
            [28.0, -25.55],  # NW corner
            [28.0, -25.7],  # Close polygon
        ]
    ]

    geojson = {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {
                    "name": "Demo Site",
                    "area_hectares": 150.0,
                    "zone": "urban_development",
                },
                "geometry": {"type": "Polygon", "coordinates": coordinates},
            }
        ],
    }

    with open(output_path, "w") as f:
        json.dump(geojson, f, indent=2)


def create_demo_roads(output_path: str) -> None:
    """Create a demo road network."""
    features = []

    # Main arterial road (horizontal)
    features.append(
        {
            "type": "Feature",
            "properties": {"road_type": "road_art", "name": "Main Avenue", "width": 20},
            "geometry": {"type": "LineString", "coordinates": [[27.95, -25.65], [28.15, -25.65]]},
        }
    )

    # Secondary road (vertical)
    features.append(
        {
            "type": "Feature",
            "properties": {"road_type": "road_sec", "name": "Central Street", "width": 15},
            "geometry": {"type": "LineString", "coordinates": [[28.05, -25.75], [28.05, -25.5]]},
        }
    )

    # Local roads
    local_roads = [
        # Horizontal local roads
        ([27.98, -25.68], [28.12, -25.68]),
        ([27.98, -25.62], [28.12, -25.62]),
        ([27.98, -25.58], [28.12, -25.58]),
        # Vertical local roads
        ([28.02, -25.72], [28.02, -25.52]),
        ([28.08, -25.72], [28.08, -25.52]),
    ]

    for i, (start, end) in enumerate(local_roads):
        features.append(
            {
                "type": "Feature",
                "properties": {"road_type": "road_loc", "name": f"Local Road {i + 1}", "width": 10},
                "geometry": {"type": "LineString", "coordinates": [start, end]},
            }
        )

    geojson = {"type": "FeatureCollection", "features": features}

    with open(output_path, "w") as f:
        json.dump(geojson, f, indent=2)


def ensure_demo_data(data_dir: str = "data") -> tuple[bool, list[str]]:
    """
    Ensure demo data files exist, creating them if necessary.

    Returns:
        Tuple of (created_files, file_paths)
    """
    data_path = Path(data_dir)
    data_path.mkdir(exist_ok=True)

    site_path = data_path / "site.geojson"
    roads_path = data_path / "roads.geojson"

    created = []

    if not site_path.exists():
        create_demo_site(str(site_path))
        created.append("site.geojson")

    if not roads_path.exists():
        create_demo_roads(str(roads_path))
        created.append("roads.geojson")

    return len(created) > 0, [str(site_path), str(roads_path)]


if __name__ == "__main__":
    # Create demo data in current directory
    created, paths = ensure_demo_data()

    if created:
        print("✅ Created demo data files:")
        for path in paths:
            if Path(path).exists():
                print(f"   📁 {path}")
    else:
        print("ℹ️ Demo data files already exist:")
        for path in paths:
            print(f"   📁 {path}")

    print("\n🚀 You can now run: python scripts/run_tui.py")
