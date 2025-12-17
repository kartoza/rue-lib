"""
Step 1: Generate Parcels
-------------------------
This example demonstrates how to generate ownership parcels from a site
polygon using a grid-based approach.

This is the first step in the RUE workflow and produces the base parcels
that will be subdivided in subsequent steps.

Usage:
    python examples/step1_generate_parcels.py [--output-dir DIR]
"""

import argparse

from rue_lib.config import MainConfig
from rue_lib.site.runner import SiteConfig, generate_parcels


def main():
    """
    Step 1: Generate parcels using grid overlay.

    This creates initial ownership parcels by:
    1. Reading site polygon and roads
    2. Creating a rectangular grid over the site
    3. Intersecting grid with site to create parcels
    4. Optionally subtracting road corridors
    """
    parser = argparse.ArgumentParser(description="Generate parcels from site polygon")
    parser.add_argument(
        "--output-dir",
        default="outputs/step1_parcels",
        help="Output directory for generated parcels",
    )
    parser.add_argument(
        "--geopackage", default="outputs/output.gpkg", help="Path to output geopackage file"
    )
    args = parser.parse_args()

    # Configure grid generation
    geopackage_path = args.geopackage

    config = SiteConfig(
        site_path="data/site.geojson",
        roads_path="data/roads.geojson",
        output_dir=args.output_dir,
        geopackage_path=geopackage_path,
        road_arterial_width_m=MainConfig.road_arterial_width_m,
        road_secondary_width_m=MainConfig.road_secondary_width_m,
    )

    # Generate parcels
    print("=" * 60)
    print("STEP 1: Generating Parcels")
    print("=" * 60)
    print(f"Site: {config.site_path}")
    print(f"Roads: {config.roads_path}")
    print(f"Output: {config.output_dir}")
    print()

    output_path = generate_parcels(config)

    print()
    print("✔ Step 1 completed successfully!")
    print(f"  Output: {output_path}")
    print()
    print("Next step: Run step2_subdivide_site.py to subdivide parcels")
    print("           into zones, subsites, and blocks")


if __name__ == "__main__":
    main()
