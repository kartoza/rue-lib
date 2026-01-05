import os

import geopandas as gpd
from shapely import difference
from shapely.geometry import CAP_STYLE, JOIN_STYLE, LinearRing, Point


def create_local_streets_zone(
    input_path,
    input_layer_name,
    output_path,
    output_layer_name,
    sidewalk_width_m,
    road_width_m,
):
    """Create local streets zone with sidewalks from grid blocks.

    This creates a zone for local streets by:
    1. Extracting vertices from the geometry
    2. Removing duplicate vertices within 0.01m distance
    3. Creating a closed linestring (LinearRing) from cleaned vertices
    4. Creating an inner buffer from the linestring using half of road_width
    5. Keeping only the buffer that's inside the original geometry
    6. Creating an outer buffer from the inner result using sidewalk_width

    The result represents the area where local streets and sidewalks will be placed.
    Both inner and outer buffer zones are saved as separate layers.

    Args:
        input_path (str): Path to the input dataset containing grid blocks.
        input_layer_name (str): Name of the layer with grid blocks.
        output_path (str): Path to the output GeoPackage (.gpkg). Created if missing.
        output_layer_name (str): Base name for the output layers.
        sidewalk_width_m (float): Width of sidewalk in meters.
        road_width_m (float): Width of local road in meters.

    Returns:
        tuple[str, str]: Names of (outer_layer, inner_layer) created.

    Raises:
        Exception: Propagated GDAL/OGR errors.
    """
    inner_buffer_distance = -(road_width_m / 2.0)
    outer_buffer_distance = sidewalk_width_m

    input_gdf = gpd.read_file(input_path, layer=input_layer_name)
    if input_gdf.empty:
        raise RuntimeError(f"Could not find or read layer: {input_layer_name}")

    crs = input_gdf.crs

    inner_geoms = []
    outer_geoms = []
    block_types = []

    for _idx, row in input_gdf.iterrows():
        geom = row.geometry
        if geom is None or geom.is_empty:
            continue

        if not geom.is_valid:
            geom = geom.buffer(0)

        if hasattr(geom, "exterior"):
            coords = list(geom.exterior.coords)
        elif hasattr(geom, "geoms"):
            coords = list(geom.geoms[0].exterior.coords) if geom.geoms else []
        else:
            coords = list(geom.coords) if hasattr(geom, "coords") else []

        if len(coords) < 2:
            continue

        cleaned_coords = [coords[0]]
        for coord in coords[1:]:
            last_point = Point(cleaned_coords[-1])
            current_point = Point(coord)
            if last_point.distance(current_point) > 0.01:
                cleaned_coords.append(coord)

        if len(cleaned_coords) < 3:
            continue

        if cleaned_coords[0] != cleaned_coords[-1]:
            cleaned_coords.append(cleaned_coords[0])

        linestring = LinearRing(cleaned_coords)

        inner_buffered = linestring.buffer(
            abs(inner_buffer_distance),
            join_style=JOIN_STYLE.mitre,
            cap_style=CAP_STYLE.square,
        )

        inner_buffered = difference(geom, inner_buffered).normalize()

        if inner_buffered.is_empty:
            continue

        if not inner_buffered.is_valid:
            inner_buffered = inner_buffered.buffer(0)
        if inner_buffered.is_empty:
            continue

        inner_buffered = inner_buffered.simplify(0.01, preserve_topology=True)
        outer_buffered = inner_buffered.buffer(
            outer_buffer_distance,
            join_style=JOIN_STYLE.mitre,
            cap_style=CAP_STYLE.square,
        )

        if outer_buffered.is_empty:
            continue

        if not outer_buffered.is_valid:
            outer_buffered = outer_buffered.buffer(0)
        if outer_buffered.is_empty:
            continue

        outer_buffered = outer_buffered.simplify(0.01, preserve_topology=True)

        inner_geoms.append(inner_buffered)
        outer_geoms.append(outer_buffered)

        block_type = row.get("block_type", "")
        if block_type is None:
            block_type = ""
        block_types.append(block_type)

    if not inner_geoms:
        print(f"No valid geometries after inner buffer from {input_layer_name}")
        return (None, None)

    inner_layer_name = f"{output_layer_name}_inner"
    outer_layer_name = f"{output_layer_name}_outer"

    inner_data = []
    outer_data = []

    for inner_geom, outer_geom, block_type in zip(inner_geoms, outer_geoms, block_types):
        inner_area = inner_geom.area
        outer_area = outer_geom.area
        sidewalk_area = outer_area - inner_area

        inner_data.append(
            {
                "geometry": inner_geom,
                "area_m2": inner_area,
                "sidewalk_area": sidewalk_area,
                "buffer_dist": abs(inner_buffer_distance),
                "sidewalk_w": sidewalk_width_m,
                "road_w": road_width_m,
                "zone_type": "buildable",
                "block_type": block_type,
            }
        )

        outer_data.append(
            {
                "geometry": outer_geom,
                "area_m2": outer_area,
                "sidewalk_area": sidewalk_area,
                "buffer_dist": abs(outer_buffer_distance),
                "sidewalk_w": sidewalk_width_m,
                "road_w": road_width_m,
                "zone_type": "street_sidewalk",
            }
        )

    inner_gdf = gpd.GeoDataFrame(inner_data, crs=crs)
    outer_gdf = gpd.GeoDataFrame(outer_data, crs=crs)

    if os.path.exists(output_path):
        from osgeo import ogr

        ds = ogr.Open(output_path, 1)
        if ds:
            layers_to_delete = []
            for i in range(ds.GetLayerCount()):
                layer = ds.GetLayerByIndex(i)
                if layer.GetName() in [inner_layer_name, outer_layer_name]:
                    layers_to_delete.append(i)

            for i in reversed(layers_to_delete):
                ds.DeleteLayer(i)
            ds = None

    if os.path.exists(output_path):
        inner_gdf.to_file(output_path, layer=inner_layer_name, driver="GPKG", mode="a")
        outer_gdf.to_file(output_path, layer=outer_layer_name, driver="GPKG", mode="a")
    else:
        inner_gdf.to_file(output_path, layer=inner_layer_name, driver="GPKG")
        outer_gdf.to_file(output_path, layer=outer_layer_name, driver="GPKG", mode="a")

    return (outer_layer_name, inner_layer_name)
