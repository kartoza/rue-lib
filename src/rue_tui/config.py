"""Configuration for the RUE TUI application."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class TuiConfig:
    """Configuration for TUI application."""

    # Input files
    site_path: str = "data/site.geojson"
    roads_path: str = "data/roads.geojson"

    # Output settings
    output_dir: str = "outputs"
    temp_dir: Optional[str] = None  # Will use system temp if None

    # Visualization settings
    image_width: int = 800
    image_height: int = 600
    dpi: int = 150

    # Map rendering
    show_labels: bool = True
    color_scheme: str = "viridis"
    background_color: str = "white"

    # Processing settings
    auto_advance: bool = False  # Auto-advance to next step on completion
    save_intermediate: bool = True  # Save intermediate results

    # File paths (computed)
    step1_output: str = field(init=False)
    step2_output: str = field(init=False)
    step3_output: str = field(init=False)
    geopackage_path: str = field(init=False)

    def __post_init__(self):
        """Initialize computed file paths."""
        output_path = Path(self.output_dir)
        self.step1_output = str(output_path / "step1_parcels")
        self.step2_output = str(output_path / "step2_streets")
        self.step3_output = str(output_path / "step3_clusters")
        self.geopackage_path = str(output_path / "output.gpkg")

        # Create output directory
        output_path.mkdir(parents=True, exist_ok=True)
