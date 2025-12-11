# src/rue_lib/core/roads.py
"""Functions for processing road geometries."""

import geopandas as gpd

from rue_lib.core.definitions import RoadTypes


def extract_roads_buffer(
    input_path: str,
    input_layer_name: str,
    output_path: str,
    output_layer_name: str,
    road_arterial_width_m: float,
    road_secondary_width_m: float,
    road_local_width_m: float,
):
    """
    Create polygon buffers from road lines based on road type.

    This function reads road line geometries and creates buffered polygons
    around them, with buffer width determined by the road type.

    Args:
        input_path: Path to input GeoPackage containing roads
        input_layer_name: Name of the layer with road line geometries
        output_path: Path to output GeoPackage
        output_layer_name: Name for the output buffered polygons layer
        road_arterial_width_m: Buffer width for arterial roads in meters
        road_secondary_width_m: Buffer width for secondary roads in meters
        road_local_width_m: Buffer width for local roads in meters

    Returns:
        Name of the created output layer
    """
    # Read roads layer
    gdf_roads = gpd.read_file(input_path, layer=input_layer_name)

    # Define road type to width mapping
    width_mapping = {
        RoadTypes.Artery: road_arterial_width_m,
        RoadTypes.Secondary: road_secondary_width_m,
        RoadTypes.Local: road_local_width_m,
    }

    output = []

    # Process each road feature
    for _idx, row in gdf_roads.iterrows():
        road_type = row.get("type", None)
        geometry = row.geometry

        if road_type is None or geometry is None:
            continue

        # Get buffer width based on road type
        buffer_width = width_mapping.get(road_type, road_local_width_m)

        # Create buffered polygon with sharp corners (mitered joins)
        # join_style=2 creates sharp corners instead of rounded
        buffered_geom = geometry.buffer(buffer_width / 2.0, join_style=2)

        # Keep all properties and update geometry
        row_dict = row.to_dict()
        row_dict["geometry"] = buffered_geom
        row_dict["type"] = road_type
        row_dict["buffer_width"] = buffer_width

        output.append(row_dict)

    # Create output GeoDataFrame
    gdf_out = gpd.GeoDataFrame(output, crs=gdf_roads.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")
    return output_layer_name
