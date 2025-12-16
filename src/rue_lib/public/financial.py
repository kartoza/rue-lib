# src/rue_lib/public/financial.py

from pathlib import Path

import geopandas as gpd

from rue_lib.financial import FinancialModel
from rue_lib.public import PublicConfig


class FinancialPublic(FinancialModel):
    """Financial attributes of public."""

    # Open space arteries
    open_art_area: float
    open_art_art_area: float
    open_art_sec_area: float
    open_art_loc_area: float

    # Open space secondaries
    open_sec_area: float
    open_sec_sec_area: float
    open_sec_loc_area: float

    # Open space local
    open_loc_area: float
    open_loc_loc_area: float

    # Open space off grid
    open_og_clus0_on_art_area: float
    open_og_clus0_on_sec_area: float
    open_og_clus0_on_loc_area: float
    open_og_clus1_on_loc_area: float
    open_og_clus2_on_loc_area: float
    open_total_area: float

    # Amenities arteries
    amen_art_area: float
    amen_art_art_area: float
    amen_art_sec_area: float
    amen_art_loc_area: float

    # Amenities secondary
    amen_sec_area: float
    amen_sec_sec_area: float
    amen_sec_loc_area: float

    # Amenities local
    amen_loc_area: float
    amen_loc_loc_area: float

    # Amenities off grid
    amen_og_clus0_on_art_area: float
    amen_og_clus0_on_sec_area: float
    amen_og_clus0_on_loc_area: float
    amen_og_clus1_on_loc_area: float
    amen_og_clus2_on_loc_area: float
    amen_total_area: float

    def __init__(self, config: PublicConfig):
        """Initialize a Site object."""
        out_dir = Path(config.output_dir)
        gpkg_path = out_dir / "outputs.gpkg"
        gpkg_path = str(gpkg_path)

        cluster_layer = gpd.read_file(gpkg_path, layer="04_final")

        # ------------------------
        # Open space
        # ------------------------
        # Open space arteries
        self.open_art_area = cluster_layer[cluster_layer["type"] == "art_os"].area.sum()
        self.open_art_art_area = cluster_layer[cluster_layer["type"] == "art_art_os"].area.sum()
        self.open_art_sec_area = cluster_layer[cluster_layer["type"] == "art_sec_os"].area.sum()
        self.open_art_loc_area = cluster_layer[cluster_layer["type"] == "art_loc_os"].area.sum()

        # Open space secondaries
        self.open_sec_area = cluster_layer[cluster_layer["type"] == "sec_os"].area.sum()
        self.open_sec_sec_area = cluster_layer[cluster_layer["type"] == "sec_sec_os"].area.sum()
        self.open_sec_loc_area = cluster_layer[cluster_layer["type"] == "sec_loc_os"].area.sum()

        # Open space local
        self.open_loc_area = cluster_layer[cluster_layer["type"] == "loc_os"].area.sum()
        self.open_loc_loc_area = cluster_layer[cluster_layer["type"] == "loc_loc_os"].area.sum()

        # Open space off grid
        open_space_off_grid = cluster_layer[cluster_layer["type"] == "off_grid0_os"]
        self.open_og_clus0_on_art_area = open_space_off_grid[cluster_layer["cluster_type"] == "art_os"].area.sum()
        self.open_og_clus0_on_sec_area = open_space_off_grid[cluster_layer["cluster_type"] == "sec_os"].area.sum()
        self.open_og_clus0_on_loc_area = open_space_off_grid[cluster_layer["cluster_type"] == "loc_os"].area.sum()
        self.open_og_clus1_on_loc_area = open_space_off_grid[cluster_layer["type"] == "off_grid1_os"].area.sum()
        self.open_og_clus2_on_loc_area = open_space_off_grid[cluster_layer["type"] == "off_grid2_os"].area.sum()

        self.open_total_area = (
            self.open_art_area +
            self.open_art_art_area +
            self.open_art_sec_area +
            self.open_art_loc_area +
            self.open_sec_area +
            self.open_sec_sec_area +
            self.open_sec_loc_area +
            self.open_loc_area +
            self.open_loc_loc_area +
            self.open_og_clus0_on_art_area +
            self.open_og_clus0_on_sec_area +
            self.open_og_clus0_on_loc_area +
            self.open_og_clus1_on_loc_area +
            self.open_og_clus2_on_loc_area
        )

        # ------------------------
        # Amenities
        # ------------------------
        # Amenities arteries
        self.amen_art_area = cluster_layer[cluster_layer["type"] == "art_am"].area.sum()
        self.amen_art_art_area = cluster_layer[cluster_layer["type"] == "art_art_am"].area.sum()
        self.amen_art_sec_area = cluster_layer[cluster_layer["type"] == "art_sec_am"].area.sum()
        self.amen_art_loc_area = cluster_layer[cluster_layer["type"] == "art_loc_am"].area.sum()

        # Amenities secondaries
        self.amen_sec_area = cluster_layer[cluster_layer["type"] == "sec_am"].area.sum()
        self.amen_sec_sec_area = cluster_layer[cluster_layer["type"] == "sec_sec_am"].area.sum()
        self.amen_sec_loc_area = cluster_layer[cluster_layer["type"] == "sec_loc_am"].area.sum()

        # Amenities local
        self.amen_loc_area = cluster_layer[cluster_layer["type"] == "loc_am"].area.sum()
        self.amen_loc_loc_area = cluster_layer[cluster_layer["type"] == "loc_loc_am"].area.sum()

        # Amenities off grid
        amenities_off_grid = cluster_layer[cluster_layer["type"] == "off_grid0_am"]
        self.amen_og_clus0_on_art_area = amenities_off_grid[cluster_layer["cluster_type"] == "art_am"].area.sum()
        self.amen_og_clus0_on_sec_area = amenities_off_grid[cluster_layer["cluster_type"] == "sec_am"].area.sum()
        self.amen_og_clus0_on_loc_area = amenities_off_grid[cluster_layer["cluster_type"] == "loc_am"].area.sum()
        self.amen_og_clus1_on_loc_area = amenities_off_grid[cluster_layer["type"] == "off_grid1_am"].area.sum()
        self.amen_og_clus2_on_loc_area = amenities_off_grid[cluster_layer["type"] == "off_grid2_am"].area.sum()

        self.amen_total_area = (
            self.amen_art_area +
            self.amen_art_art_area +
            self.amen_art_sec_area +
            self.amen_art_loc_area +
            self.amen_sec_area +
            self.amen_sec_sec_area +
            self.amen_sec_loc_area +
            self.amen_loc_area +
            self.amen_loc_loc_area +
            self.amen_og_clus0_on_art_area +
            self.amen_og_clus0_on_sec_area +
            self.amen_og_clus0_on_loc_area +
            self.amen_og_clus1_on_loc_area +
            self.amen_og_clus2_on_loc_area
        )

        self.save(config.output_dir)
