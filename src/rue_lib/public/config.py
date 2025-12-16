# src/rue_lib/public/config.py
from dataclasses import dataclass


@dataclass
class PublicConfig:
    """Configuration for public spaces generation."""

    input_path: str
    site_path: str = "outputs/step1_parcels/parcels.geojson"
    output_dir: str = "outputs/public"

    open_percent: float = 4.0  # Minimum percentage of open space required
    amen_percent: float = 10.0  # Minimum percentage of amenities required
