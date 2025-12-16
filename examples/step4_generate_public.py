# examples/step4_generate_public.py

from rue_lib.config import MainConfig
from rue_lib.public import PublicConfig, generate_public


def main():
    """
    Step 4: Generate public spaces from clusters.

    This step takes the output from step 3 (clusters) and creates public spaces.
    """

    config = PublicConfig(
        site_path="outputs/step1_parcels/parcels.geojson",
        input_path="outputs/step3_clusters/outputs.geojson",
        output_dir="outputs/step4_public",
        open_percent=MainConfig.open_percent,
        amen_percent=MainConfig.amen_percent,
    )

    print("=" * 60)
    print("STEP 4: Generating Public Spaces")
    print("=" * 60)
    print(f"Input: {config.input_path}")
    print(f"Output: {config.output_dir}")
    print()

    output_path = generate_public(config)

    print()
    print("✔ Step 4 completed successfully!")
    print(f"  Output: {output_path}")


if __name__ == "__main__":
    main()
