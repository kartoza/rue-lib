# src/rue_lib/streets/financial.py


import geopandas as gpd

from rue_lib.financial import FinancialModel
from rue_lib.site.config import SiteConfig


class FinancialSite(FinancialModel):
    """Financial attributes of site."""

    site_area_total: float
    site_roads_area: float

    def __init__(self, config: SiteConfig):
        """Initialize a Site object."""
        gpkg_path = config.geopackage_path

        site_gdf = gpd.read_file(gpkg_path, layer="site")
        roads_buffer_gdf = gpd.read_file(gpkg_path, layer="roads_buffered")
        self.site_area_total = site_gdf.area.sum()
        self.site_roads_area = roads_buffer_gdf.area.sum()

        site_gdf = None
        roads_buffer_gdf = None

        self.save(config.output_dir)
