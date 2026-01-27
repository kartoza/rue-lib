# src/rue_lib/streets/financial.py

from pathlib import Path

import geopandas as gpd
from shapely.ops import unary_union

from rue_lib.financial import FinancialModel
from rue_lib.streets import StreetConfig
from rue_lib.streets.geometry_utils import trim_line


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
        out_dir = Path(config.output_dir)
        output_path = out_dir / "outputs.gpkg"
        output_path = str(output_path)
        roads = gpd.read_file(output_path, layer="00_roads")
        roads_local = gpd.read_file(output_path, layer=local_roads_layer)

        site = gpd.read_file(output_path, layer="00_site")
        site_union = unary_union(site.geometry)

        # Get all arterial and secondary roads
        art_all = roads[roads["road_type"] == "road_art"]
        sec_all = roads[roads["road_type"] == "road_sec"]

        art_half_width = config.road_arterial_width_m / 2.0
        site_buffered_art = site_union.buffer(art_half_width, cap_style=2, join_style=2)

        art_intersecting = art_all[art_all.intersects(site_buffered_art)]

        if not art_intersecting.empty and not site_buffered_art.is_empty:
            art_clipped = art_intersecting.copy()
            art_clipped["geometry"] = art_clipped.geometry.intersection(site_buffered_art)
            art_clipped = art_clipped[~art_clipped.geometry.is_empty]

            art_clipped["geometry"] = art_clipped.geometry.apply(
                lambda geom: trim_line(geom, art_half_width)
            )
            art_clipped = art_clipped[art_clipped.geometry.notna()]
            art_clipped = art_clipped[~art_clipped.geometry.is_empty]

            art_100_clipped = art_clipped[
                (art_clipped["road_pcent"] == 100) | (art_clipped["road_pcent"].isna())
            ]
            art_50_clipped = art_clipped[art_clipped["road_pcent"] == 50]

            art_100_length = sum(
                s.geometry.length + config.road_locals_width_m
                for _, s in art_100_clipped.iterrows()
            )
            art_50_length = sum(
                s.geometry.length + config.road_locals_width_m for _, s in art_50_clipped.iterrows()
            )

            art_clipped.to_file(output_path, layer="000_art_clipped", driver="GPKG")

            art_clipped_buffered = art_clipped.copy()
            art_clipped_buffered["geometry"] = art_clipped.geometry.buffer(
                art_half_width, cap_style="flat", join_style="mitre"
            )
            art_clipped_buffered.to_file(
                output_path, layer="001_art_clipped_buffered", driver="GPKG"
            )
            art_area = art_clipped_buffered.area.sum()
        else:
            art_100_length = 0.0
            art_50_length = 0.0
            art_area = 0.0

        sec_half_width = config.road_secondary_width_m / 2.0
        site_buffered_sec = site_union.buffer(sec_half_width, cap_style=2, join_style=2)

        sec_intersecting = sec_all[sec_all.intersects(site_buffered_sec)]

        if not sec_intersecting.empty and not site_buffered_sec.is_empty:
            sec_clipped = sec_intersecting.copy()
            sec_clipped["geometry"] = sec_clipped.geometry.intersection(site_buffered_sec)
            sec_clipped = sec_clipped[~sec_clipped.geometry.is_empty]

            sec_clipped["geometry"] = sec_clipped.geometry
            sec_clipped = sec_clipped[sec_clipped.geometry.notna()]
            sec_clipped = sec_clipped[~sec_clipped.geometry.is_empty]

            sec_100_clipped = sec_clipped[
                (sec_clipped["road_pcent"] == 100) | (sec_clipped["road_pcent"].isna())
            ]
            sec_50_clipped = sec_clipped[sec_clipped["road_pcent"] == 50]

            sec_100_length = sum(
                s.geometry.length - sec_half_width + config.road_locals_width_m
                for _, s in sec_100_clipped.iterrows()
            )
            sec_50_length = sum(
                s.geometry.length - sec_half_width + config.road_locals_width_m
                for _, s in sec_50_clipped.iterrows()
            )

            sec_clipped.to_file(output_path, layer="002_sec_clipped", driver="GPKG")

            sec_clipped_buffered = sec_clipped.copy()
            sec_clipped_buffered["geometry"] = sec_clipped.geometry.buffer(
                sec_half_width, cap_style="flat", join_style="mitre"
            )
            sec_clipped_buffered.to_file(
                output_path, layer="003_sec_clipped_buffered", driver="GPKG"
            )
            sec_area = sum(
                s.geometry.area - sec_half_width * (2 * sec_half_width)
                for _, s in sec_clipped_buffered.iterrows()
            )
        else:
            sec_100_length = 0.0
            sec_50_length = 0.0
            sec_area = 0.0

        self.road_len_art_100 = art_100_length
        self.road_len_sec_100 = sec_100_length
        self.road_len_loc_100 = roads_local.length.sum()

        self.road_len_art_50 = art_50_length
        self.road_len_sec_50 = sec_50_length
        self.road_len_loc_50 = 0

        # Calculate areas
        local_roads_buffer = gpd.read_file(output_path, layer=local_roads_buffer_layer)
        self.road_area_art = art_area
        self.road_area_sec = sec_area
        self.road_area_loc = local_roads_buffer.area.sum()

        roads = None
        roads_local = None
        local_roads_buffer = None

        self.save(output_dir=config.output_dir)
