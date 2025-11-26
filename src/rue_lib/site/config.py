from dataclasses import dataclass


@dataclass
class SiteConfig:
    """Configuration for site generation."""

    site_path: str
    roads_path: str
    output_dir: str = "outputs/site"
    subtract_roads: bool = True

    road_arterial_width_m: float = 20.0
    road_secondary_width_m: float = 15.0
