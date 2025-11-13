"""CLI for Step 1."""

import argparse
import sys
from pathlib import Path

from rue_lib.site.runner import SiteConfig, generate_parcels


def main():
    """Phase 1 CLI entry point."""
    script_dir = Path.cwd()
    default_site = script_dir / "site.geojson"
    default_roads = script_dir / "roads.geojson"

    parser = argparse.ArgumentParser(
        description="Step 1 (site_mob) – produce final landowners GeoJSON"
    )
    parser.add_argument(
        "--site",
        type=str,
        dest="site_path",
        default=str(default_site),
        help="Vector file for site polygon",
    )
    parser.add_argument(
        "--roads",
        type=str,
        dest="roads_path",
        default=str(default_roads),
        help="Vector file for roads",
    )
    parser.add_argument("--rows", type=int, default=4, help="Grid rows")
    parser.add_argument("--cols", type=int, default=4, help="Grid cols")
    parser.add_argument("--pad", type=float, default=10.0, help="Grid padding (meters)")
    parser.add_argument(
        "--trim-by-roads", action="store_true", help="Carve road corridors out of landowners"
    )
    parser.add_argument(
        "--min_owner_area", type=float, default=5.0, help="Min owner polygon area (m²)"
    )
    parser.add_argument(
        "--out", type=str, default="outputs/site_mob_phase1", help="Output directory"
    )

    args = parser.parse_args()

    if not Path(args.site_path).exists():
        print(f"✖ Site file not found: {args.site_path}", file=sys.stderr)
        sys.exit(1)

    if not Path(args.roads_path).exists():
        print(f"✖ Roads file not found: {args.roads_path}", file=sys.stderr)
        sys.exit(1)

    cfg = SiteConfig(
        site_path=args.site_path,
        roads_path=args.roads_path,
        output_dir=args.out,
        rows=args.rows,
        cols=args.cols,
        pad_m=args.pad,
    )

    out = generate_parcels(cfg)
    print(f"✔ Final landowners GeoJSON: {out}")


if __name__ == "__main__":
    main()
