# src/rue_lib/streets/financial.py

from pathlib import Path

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

    def __init__(self, config: StreetConfig):
        """Initialize a FinancialStreet object."""
        out_dir = Path(config.output_dir)
        output_path = out_dir / "outputs.gpkg"
        output_path = str(output_path)

        roads = gpd.read_file(output_path, layer="00_roads")
        roads_local = gpd.read_file(output_path, layer="18_local_roads")

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

        roads = None
        roads_local = None

        self.save(output_dir=config.output_dir)
