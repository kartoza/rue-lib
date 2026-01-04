# src/rue_lib/streets/financial.py


import geopandas as gpd

from rue_lib.financial import FinancialModel
from rue_lib.streets import StreetConfig


class FinancialStreet(FinancialModel):
    """Financial attributes of a street."""

    road_len_art_100: float
    road_len_sec_100: float
    road_len_loc_100: float

    road_len_art_50: float
    road_len_sec_50: float
    road_len_loc_50: float

    # Area of roads
    road_area_art: float
    road_area_sec: float
    road_area_loc: float

    def __init__(
        self,
        config: StreetConfig,
        local_roads_layer: str = "18_local_roads",
        local_roads_buffer_layer: str = "17_local_roads_buffer",
    ):
        """Initialize a FinancialStreet object."""
        output_path = config.geopackage_path

        roads = gpd.read_file(output_path, layer="00_roads")
        roads_local = gpd.read_file(output_path, layer=local_roads_layer)

        # Filter for 100% roads (road_pcent = 100 or road_pcent doesn't exist)
        art_100 = roads[
            (roads["road_type"] == "road_art")
            & ((roads["road_pcent"] == 100) | (roads["road_pcent"].isna()))
        ]
        sec_100 = roads[
            (roads["road_type"] == "road_sec")
            & ((roads["road_pcent"] == 100) | (roads["road_pcent"].isna()))
        ]

        # Filter for 50% roads
        art_50 = roads[(roads["road_type"] == "road_art") & (roads["road_pcent"] == 50)]
        sec_50 = roads[(roads["road_type"] == "road_sec") & (roads["road_pcent"] == 50)]

        # Calculate lengths
        self.road_len_art_100 = art_100.length.sum()
        self.road_len_sec_100 = sec_100.length.sum()
        self.road_len_loc_100 = roads_local.length.sum()

        self.road_len_art_50 = art_50.length.sum()
        self.road_len_sec_50 = sec_50.length.sum()
        self.road_len_loc_50 = 0

        # Calculate areas
        road_buffer = gpd.read_file(output_path, layer="00_roads_buffer")
        local_roads_buffer = gpd.read_file(output_path, layer=local_roads_buffer_layer)
        self.road_area_art = road_buffer[(road_buffer["road_type"] == "road_art")].area.sum()
        self.road_area_sec = road_buffer[(road_buffer["road_type"] == "road_sec")].area.sum()
        self.road_area_loc = local_roads_buffer.area.sum()

        roads = None
        roads_local = None
        road_buffer = None
        local_roads_buffer = None

        self.save(output_dir=config.output_dir)
