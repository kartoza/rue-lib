# src/rue_lib/streets/financial.py

import geopandas as gpd

from rue_lib.financial import FinancialModel
from rue_lib.site.config import SiteConfig


class FinancialSite(FinancialModel):
    """Financial attributes of site."""

    site_area_total: float
    site_roads_area: float

    def __init__(
            self, config: SiteConfig,
            site: gpd.GeoDataFrame,
            roads: gpd.GeoDataFrame
    ):
        """Initialize a Site object."""
        self.site_area_total = site.area.sum()
        self.site_roads_area = roads.area.sum()

        self.save(config.output_dir)
