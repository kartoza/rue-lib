"""Map visualization utilities for the TUI."""

import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


class MapVisualizer:
    """Handles map rendering and visualization for the TUI."""

    def __init__(self, width: int = 800, height: int = 600, dpi: int = 150):
        self.width = width
        self.height = height
        self.dpi = dpi

    def render_step_results(
        self, step: int, geopackage_path: str, geojson_path: Optional[str] = None
    ) -> str:
        """
        Render results for a specific step to an image file.

        Args:
            step: Step number (1, 2, or 3)
            geopackage_path: Path to the main geopackage
            geojson_path: Optional path to step-specific geojson

        Returns:
            Path to the generated image file
        """
        fig, ax = plt.subplots(figsize=(self.width / 100, self.height / 100), dpi=self.dpi)

        try:
            if step == 1:
                return self._render_step1(ax, geopackage_path, geojson_path)
            elif step == 2:
                return self._render_step2(ax, geopackage_path, geojson_path)
            elif step == 3:
                return self._render_step3(ax, geopackage_path, geojson_path)
            else:
                raise ValueError(f"Invalid step number: {step}")

        except Exception as e:
            # Fallback: render error message
            ax.text(
                0.5,
                0.5,
                f"Error rendering step {step}:\n{str(e)}",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=12,
                color="#CC0403",
            )
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            return self._save_figure(fig, f"step{step}_error")

    def _render_step1(self, ax, geopackage_path: str, geojson_path: str) -> str:
        """Render Step 1: Parcels."""
        ax.set_title("Step 1: Generated Parcels", fontsize=16, fontweight="bold")

        # Load and display site and roads
        try:
            site = gpd.read_file(geopackage_path, layer="site")
            roads = gpd.read_file(geopackage_path, layer="roads")

            # Plot site boundary
            site.plot(
                ax=ax,
                color="lightblue",
                alpha=0.3,
                edgecolor="#569FC6",
                linewidth=2,
                label="Site Boundary",
            )

            # Plot roads
            roads.plot(ax=ax, color="#8A8B8B", linewidth=2, label="Roads")

        except Exception:
            # Fallback to basic plot
            pass

        # Load and display parcels
        if geojson_path and Path(geojson_path).exists():
            parcels = gpd.read_file(geojson_path)
            parcels.plot(
                ax=ax,
                color="lightgreen",
                alpha=0.7,
                edgecolor="#06969A",
                linewidth=1,
                label="Parcels",
            )

        self._style_map(ax)
        ax.legend(loc="upper right")
        return self._save_figure(plt.gcf(), "step1_parcels")

    def _render_step2(self, ax, geopackage_path: str, geojson_path: str) -> str:
        """Render Step 2: Streets and Blocks."""
        ax.set_title("Step 2: Street Blocks", fontsize=16, fontweight="bold")

        # Kartoza color scheme for different grid types
        colors = {
            "off_grid": "#CC0403",  # Kartoza red
            "on_grid_arterial": "#06969A",  # Kartoza teal
            "on_grid_secondary": "#569FC6",  # Kartoza blue
            "cold_boundary": "#DF9E2F",  # Kartoza yellow/orange
            "arterial_roads": "#569FC6",  # Kartoza blue
            "secondary_roads": "#8A8B8B",  # Kartoza grey
            "local_roads": "#8A8B8B",  # Kartoza grey (lighter)
        }

        try:
            # Try to load the merged grids layer
            if Path(geopackage_path).exists():
                # Load different layers from geopackage
                layers_to_plot = [
                    ("04_arterial_roads", colors["arterial_roads"], "Arterial Roads"),
                    ("05_secondary_roads", colors["secondary_roads"], "Secondary Roads"),
                    ("local_roads", colors["local_roads"], "Local Roads"),
                    ("17_all_grids_merged", None, "Grid Cells"),
                ]

                legend_elements = []

                for layer_name, color, label in layers_to_plot:
                    try:
                        gdf = gpd.read_file(geopackage_path, layer=layer_name)
                        if not gdf.empty:
                            if layer_name == "17_all_grids_merged":
                                # Color by grid_type if available
                                if "grid_type" in gdf.columns:
                                    for grid_type in gdf["grid_type"].unique():
                                        subset = gdf[gdf["grid_type"] == grid_type]
                                        color = colors.get(grid_type, "#CCCCCC")
                                        subset.plot(
                                            ax=ax,
                                            color=color,
                                            alpha=0.7,
                                            edgecolor="black",
                                            linewidth=0.5,
                                        )
                                        legend_elements.append(
                                            mpatches.Patch(
                                                color=color,
                                                label=grid_type.replace("_", " ").title(),
                                            )
                                        )
                                else:
                                    gdf.plot(
                                        ax=ax,
                                        color="lightblue",
                                        alpha=0.7,
                                        edgecolor="#569FC6",
                                        linewidth=0.5,
                                    )
                            else:
                                gdf.plot(ax=ax, color=color, linewidth=2, label=label)
                                legend_elements.append(mpatches.Patch(color=color, label=label))
                    except Exception:
                        continue

                if legend_elements:
                    ax.legend(handles=legend_elements, loc="upper right", fontsize=8)

        except Exception:
            # Fallback: try to load from geojson
            if geojson_path and Path(geojson_path).exists():
                streets = gpd.read_file(geojson_path)
                streets.plot(ax=ax, color="#DF9E2F", alpha=0.7, edgecolor="#CC0403", linewidth=1)

        self._style_map(ax)
        return self._save_figure(plt.gcf(), "step2_streets")

    def _render_step3(self, ax, geopackage_path: str, geojson_path: str) -> str:
        """Render Step 3: Clusters and Final Layout."""
        ax.set_title("Step 3: Urban Clusters", fontsize=16, fontweight="bold")

        try:
            # Load final clusters layer
            if Path(geopackage_path).exists():
                try:
                    final = gpd.read_file(geopackage_path, layer="300_final")

                    # Create colormap for different feature types
                    if "type" in final.columns:
                        unique_types = final["type"].unique()
                        cmap = plt.cm.Set3(np.linspace(0, 1, len(unique_types)))
                        type_colors = dict(zip(unique_types, cmap))

                        for ftype in unique_types:
                            subset = final[final["type"] == ftype]
                            color = type_colors[ftype]
                            subset.plot(
                                ax=ax,
                                color=color,
                                alpha=0.8,
                                edgecolor="black",
                                linewidth=0.5,
                                label=ftype,
                            )
                    else:
                        final.plot(
                            ax=ax, color="#06969A", alpha=0.7, edgecolor="#569FC6", linewidth=0.5
                        )

                    ax.legend(loc="upper right", fontsize=8)

                except Exception:
                    # Try other final layers
                    for layer in ["111_warm_block_final", "17_all_grids_merged"]:
                        try:
                            gdf = gpd.read_file(geopackage_path, layer=layer)
                            gdf.plot(
                                ax=ax,
                                color="#06969A",
                                alpha=0.7,
                                edgecolor="#569FC6",
                                linewidth=0.5,
                            )
                            break
                        except Exception:
                            continue

        except Exception:
            # Fallback text
            ax.text(
                0.5,
                0.5,
                "Step 3 Results\n(Visualization pending)",
                ha="center",
                va="center",
                transform=ax.transAxes,
                fontsize=14,
                color="#06969A",
            )

        self._style_map(ax)
        return self._save_figure(plt.gcf(), "step3_clusters")

    def _style_map(self, ax):
        """Apply consistent styling to map."""
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("Longitude", fontsize=10)
        ax.set_ylabel("Latitude", fontsize=10)

        # Remove ticks for cleaner look
        ax.tick_params(axis="both", which="major", labelsize=8)

    def _save_figure(self, fig, name: str) -> str:
        """Save figure to temporary file."""
        # Create temp file
        temp_file = tempfile.NamedTemporaryFile(
            delete=False, suffix=".png", prefix=f"rue_tui_{name}_"
        )
        temp_path = temp_file.name
        temp_file.close()

        # Save with high quality
        fig.savefig(
            temp_path, dpi=self.dpi, bbox_inches="tight", facecolor="white", edgecolor="none"
        )
        plt.close(fig)

        return temp_path

    def display_image_in_terminal(self, image_path: str) -> bool:
        """
        Display image in terminal using fim or similar tool.

        Returns:
            True if successful, False otherwise
        """
        if not Path(image_path).exists():
            return False

        # Try different image viewers in order of preference
        viewers = [
            ["fim", "-a", image_path],  # Frame buffer image viewer
            ["feh", "--auto-zoom", image_path],  # X11 image viewer
            ["display", image_path],  # ImageMagick
            ["eog", image_path],  # GNOME image viewer
        ]

        for viewer in viewers:
            try:
                # Check if command exists
                subprocess.run(["which", viewer[0]], capture_output=True, check=True)

                # Run viewer in background
                subprocess.Popen(viewer, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                return True

            except (subprocess.CalledProcessError, FileNotFoundError):
                continue

        return False

    def get_layer_info(self, geopackage_path: str) -> list[tuple[str, int]]:
        """
        Get information about layers in the geopackage.

        Returns:
            List of (layer_name, feature_count) tuples
        """
        info = []
        try:
            import fiona

            with fiona.open(geopackage_path):
                pass  # Just to check if file is valid

            # Try to read common layer names
            common_layers = [
                "site",
                "roads",
                "parcels",
                "sites",
                "04_arterial_roads",
                "05_secondary_roads",
                "local_roads",
                "17_all_grids_merged",
                "111_warm_block_final",
                "300_final",
            ]

            for layer in common_layers:
                try:
                    gdf = gpd.read_file(geopackage_path, layer=layer)
                    info.append((layer, len(gdf)))
                except Exception:
                    continue

        except Exception:
            pass

        return info
