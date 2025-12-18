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

    # Lot area
    lot_art_area: float
    lot_art_art_area: float
    lot_art_sec_area: float
    lot_art_loc_area: float

    lot_sec_area: float
    lot_sec_sec_area: float
    lot_sec_loc_area: float

    lot_loc_area: float
    lot_loc_loc_area: float

    # Lot area off grid
    og_clus0_on_art_area: float
    og_clus0_on_sec_area: float
    og_clus0_on_loc_area: float
    og_clus1_on_loc_area: float
    og_clus2_on_loc_area: float

    # Lot number
    lot_art_num: int
    lot_art_art_num: int
    lot_art_sec_num: int
    lot_art_loc_num: int

    lot_sec_num: int
    lot_sec_sec_num: int
    lot_sec_loc_num: int

    lot_loc_num: int
    lot_loc_loc_num: int

    # Lot number off grid
    og_clus0_on_art_num: int
    og_clus0_on_sec_num: int
    og_clus0_on_loc_num: int
    og_clus1_on_loc_num: int
    og_clus2_on_loc_num: int

    og_entr0_on_art_area: float
    og_entr0_on_sec_area: float
    og_entr0_on_loc_area: float
    og_entr1_on_loc_area: float

    # Other
    site_total_area: float

    def type_sum(self, layer: gpd.GeoDataFrame, type: str):
        """Sum by type."""
        return layer[layer["type"] == type].area.sum()

    def cluster_type_sum(self, layer: gpd.GeoDataFrame, type: str):
        """Sum by type."""
        return layer[layer["cluster_type"] == type].area.sum()

    def type_num(self, layer: gpd.GeoDataFrame, type: str):
        """Numby type."""
        return len(layer[layer["type"] == type])

    def cluster_type_num(self, layer: gpd.GeoDataFrame, type: str):
        """Num by type."""
        return len(layer[layer["cluster_type"] == type])

    def __init__(self, config: PublicConfig):
        """Initialize a Site object."""
        out_dir = Path(config.output_dir)
        gpkg_path = out_dir / "outputs.gpkg"
        gpkg_path = str(gpkg_path)

        layer = gpd.read_file(gpkg_path, layer="04_final")

        # Reproject to appropriate projected CRS for accurate area calculations
        if layer.crs and layer.crs.is_geographic:
            layer = layer.to_crs(layer.estimate_utm_crs())

        # ------------------------
        # Open space
        # ------------------------
        # Open space arteries
        self.open_art_area = self.type_sum(layer, "art_os")
        self.open_art_art_area = self.type_sum(layer, "art_art_os")
        self.open_art_sec_area = self.type_sum(layer, "art_sec_os")
        self.open_art_loc_area = self.type_sum(layer, "art_loc_os")

        # Open space secondaries
        self.open_sec_area = self.type_sum(layer, "sec_os")
        self.open_sec_sec_area = self.type_sum(layer, "sec_sec_os")
        self.open_sec_loc_area = self.type_sum(layer, "sec_loc_os")

        # Open space local
        self.open_loc_area = self.type_sum(layer, "loc_os")
        self.open_loc_loc_area = self.type_sum(layer, "loc_loc_os")

        # Open space off grid
        off_grid = layer[layer["type"] == "off_grid0_os"]
        self.open_og_clus0_on_art_area = self.cluster_type_sum(off_grid, "art_os")
        self.open_og_clus0_on_sec_area = self.cluster_type_sum(off_grid, "sec_os")
        self.open_og_clus0_on_loc_area = self.cluster_type_sum(off_grid, "loc_os")

        self.open_og_clus1_on_loc_area = self.type_sum(layer, "off_grid1_os")
        self.open_og_clus2_on_loc_area = self.type_sum(layer, "off_grid2_os")

        self.open_total_area = (
            self.open_art_area
            + self.open_art_art_area
            + self.open_art_sec_area
            + self.open_art_loc_area
            + self.open_sec_area
            + self.open_sec_sec_area
            + self.open_sec_loc_area
            + self.open_loc_area
            + self.open_loc_loc_area
            + self.open_og_clus0_on_art_area
            + self.open_og_clus0_on_sec_area
            + self.open_og_clus0_on_loc_area
            + self.open_og_clus1_on_loc_area
            + self.open_og_clus2_on_loc_area
        )

        # ------------------------
        # Amenities
        # ------------------------
        # Amenities arteries
        self.amen_art_area = self.type_sum(layer, "art_am")
        self.amen_art_art_area = self.type_sum(layer, "art_art_am")
        self.amen_art_sec_area = self.type_sum(layer, "art_sec_am")
        self.amen_art_loc_area = self.type_sum(layer, "art_loc_am")

        # Amenities secondaries
        self.amen_sec_area = self.type_sum(layer, "sec_am")
        self.amen_sec_sec_area = self.type_sum(layer, "sec_sec_am")
        self.amen_sec_loc_area = self.type_sum(layer, "sec_loc_am")

        # Amenities local
        self.amen_loc_area = self.type_sum(layer, "loc_am")
        self.amen_loc_loc_area = self.type_sum(layer, "loc_loc_am")

        # Amenities off grid
        off_grid = layer[layer["type"] == "off_grid0_am"]
        self.amen_og_clus0_on_art_area = self.cluster_type_sum(off_grid, "art_am")
        self.amen_og_clus0_on_sec_area = self.cluster_type_sum(off_grid, "sec_am")
        self.amen_og_clus0_on_loc_area = self.cluster_type_sum(off_grid, "loc_am")

        self.amen_og_clus1_on_loc_area = self.type_sum(layer, "off_grid1_am")
        self.amen_og_clus2_on_loc_area = self.type_sum(layer, "off_grid2_am")

        self.amen_total_area = (
            self.amen_art_area
            + self.amen_art_art_area
            + self.amen_art_sec_area
            + self.amen_art_loc_area
            + self.amen_sec_area
            + self.amen_sec_sec_area
            + self.amen_sec_loc_area
            + self.amen_loc_area
            + self.amen_loc_loc_area
            + self.amen_og_clus0_on_art_area
            + self.amen_og_clus0_on_sec_area
            + self.amen_og_clus0_on_loc_area
            + self.amen_og_clus1_on_loc_area
            + self.amen_og_clus2_on_loc_area
        )

        # ------------------------
        # Lot Area
        # ------------------------
        self.lot_art_area = self.type_sum(layer, "art")
        self.lot_art_art_area = self.type_sum(layer, "art_art")
        self.lot_art_sec_area = self.type_sum(layer, "art_sec")
        self.lot_art_loc_area = self.type_sum(layer, "art_loc")

        self.lot_sec_area = self.type_sum(layer, "sec")
        self.lot_sec_sec_area = self.type_sum(layer, "sec_sec")
        self.lot_sec_loc_area = self.type_sum(layer, "sec_loc")

        self.lot_loc_area = self.type_sum(layer, "loc")
        self.lot_loc_loc_area = self.type_sum(layer, "loc_loc")

        # ------------------------
        # Lot Area Off Grid
        # ------------------------
        off_grid = layer[layer["type"] == "off_grid0"]
        self.og_clus0_on_art_area = self.cluster_type_sum(off_grid, "art")
        self.og_clus0_on_sec_area = self.cluster_type_sum(off_grid, "sec")
        self.og_clus0_on_loc_area = self.cluster_type_sum(off_grid, "loc")
        self.og_clus1_on_loc_area = self.type_sum(layer, "off_grid1")
        self.og_clus2_on_loc_area = self.type_sum(layer, "off_grid2")

        # ------------------------
        # Lot number
        # ------------------------
        self.lot_art_num = self.type_num(layer, "art")
        self.lot_art_art_num = self.type_num(layer, "art_art")
        self.lot_art_sec_num = self.type_num(layer, "art_sec")
        self.lot_art_loc_num = self.type_num(layer, "art_loc")

        self.lot_sec_num = self.type_num(layer, "sec")
        self.lot_sec_sec_num = self.type_num(layer, "sec_sec")
        self.lot_sec_loc_num = self.type_num(layer, "sec_loc")

        self.lot_loc_num = self.type_num(layer, "loc")
        self.lot_loc_loc_num = self.type_num(layer, "loc_loc")

        # ------------------------
        # Lot number off grid
        # ------------------------
        off_grid = layer[layer["type"] == "off_grid0"]
        self.og_clus0_on_art_num = self.type_num(off_grid, "art")
        self.og_clus0_on_sec_num = self.type_num(off_grid, "sec")
        self.og_clus0_on_loc_num = self.type_num(off_grid, "loc")
        self.og_clus1_on_loc_num = self.type_num(layer, "off_grid1")
        self.og_clus2_on_loc_num = self.type_num(layer, "off_grid2")

        self.save(config.output_dir)
