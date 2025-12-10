"""
Step 1: Generate Parcels
-------------------------
This example demonstrates how to generate ownership parcels from a site
polygon using a grid-based approach.

This is the first step in the RUE workflow and produces the base parcels
that will be subdivided in subsequent steps.

Usage:
    python examples/step1_generate_parcels.py
"""

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

    # Configure grid generation
    config = SiteConfig(
        site_path="data/site.geojson",
        roads_path="data/roads.geojson",
        output_dir="outputs/step1_parcels",
        geopackage_path="outputs/output.gpkg",
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
