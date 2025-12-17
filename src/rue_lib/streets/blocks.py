# src/rue_lib/streets/blocks.py
from pathlib import Path

from osgeo import ogr

from rue_lib.streets.config import StreetConfig
from rue_lib.streets.operations import (
    break_multipart_features,
    cleanup_grid_blocks,
    clip_layer,
    create_grid_from_on_grid,
    erase_layer,
)

from .geometry_utils import break_linestring_by_angle


def explode_multipolygon(geom):
    """Convert a multipolygon to individual polygons."""
    polygons = []
    geom_type = ogr.GT_Flatten(geom.GetGeometryType())

    if geom_type == ogr.wkbPolygon:
        polygons.append(geom.Clone())
    elif geom_type == ogr.wkbMultiPolygon:
        for i in range(geom.GetGeometryCount()):
            sub_geom = geom.GetGeometryRef(i)
            sub_geom_type = ogr.GT_Flatten(sub_geom.GetGeometryType())
            if sub_geom_type == ogr.wkbPolygon:
                polygons.append(sub_geom.Clone())
    return polygons


def merge_lines(
    gpkg_path,
    perpendicular_lines_layer,
    arterial_edge_lines_layer,
    street_blocks_layer,
    secondary_road_setback_layer,
    arterial_road_setback_layer,
    output_layer_name,
):
    """
    Merge multiple line sources into one layer, tagging setback provenance.

    Writes a unified line layer from:
      * `perpendicular_lines_layer`
      * `arterial_edge_lines_layer`
      * polygon boundaries of `street_blocks_layer`
      * polygon boundaries of `secondary_road_setback_layer`
      * polygon boundaries of `arterial_road_setback_layer`
      * (optional) site boundary

    Each output feature has:
      - source (str): producer of the line (e.g., 'perpendicular', 'arterial_edge',
        'street_blocks', 'secondary_setback', 'arterial_setback', 'site_boundary')
      - setback (str|NULL): 'secondary' or 'arterial' for setback boundaries; NULL otherwise
      - poly_id (int|NULL): FID of the source polygon when applicable
      - poly_area (float|NULL): Area of the source polygon when applicable

    Args:
        gpkg_path (str): GeoPackage path to read/write.
        perpendicular_lines_layer (str): Name of perpendicular lines layer.
        arterial_edge_lines_layer (str): Name of arterial edge lines layer.
        street_blocks_layer (str): Name of street blocks polygon layer.
        secondary_road_setback_layer (str): Name of secondary setback polygon layer.
        arterial_road_setback_layer (str): Name of arterial setback polygon layer.
        output_layer_name (str): Name of merged output line layer.

    Returns:
        None
    """
    ds = ogr.Open(gpkg_path, 1)
    line_id = 1

    perp_layer = ds.GetLayerByName(perpendicular_lines_layer)
    blocks_layer = ds.GetLayerByName(street_blocks_layer)
    secondary_layer = ds.GetLayerByName(secondary_road_setback_layer)
    arterial_setback_layer = ds.GetLayerByName(arterial_road_setback_layer)

    srs = perp_layer.GetSpatialRef()

    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbLineString25D)

    source_field = ogr.FieldDefn("source", ogr.OFTString)
    output_layer.CreateField(source_field)

    setback_field = ogr.FieldDefn("setback", ogr.OFTString)
    output_layer.CreateField(setback_field)

    poly_id_field = ogr.FieldDefn("line_id", ogr.OFTInteger)
    output_layer.CreateField(poly_id_field)

    poly_area_field = ogr.FieldDefn("poly_area", ogr.OFTReal)
    output_layer.CreateField(poly_area_field)

    line_count = 0

    for feature in perp_layer:
        geom = feature.GetGeometryRef()
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(geom)
        out_feature.SetField("source", "perpendicular")
        out_feature.SetField("setback", None)
        out_feature.SetField("line_id", line_id)
        output_layer.CreateFeature(out_feature)
        out_feature = None
        line_count += 1
        line_id += 1

    # Add arterial_road_setback boundaries
    arterial_lines = []
    for feature in arterial_setback_layer:
        poly = feature.GetGeometryRef()
        boundary = poly.Boundary()
        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(boundary)
        out_feature.SetField("source", "arterial_setback")
        out_feature.SetField("setback", "arterial")
        out_feature.SetField("line_id", line_id)
        output_layer.CreateFeature(out_feature)
        arterial_lines.append(boundary)
        out_feature = None
        line_count += 1
        line_id += 1

    for feature in blocks_layer:
        poly = feature.GetGeometryRef()
        boundary = poly.Boundary()
        for i in range(boundary.GetGeometryCount()):
            line = boundary.GetGeometryRef(i)
            segments = break_linestring_by_angle(line)
            for segment in segments:
                print(segment)
                out_feature = ogr.Feature(output_layer.GetLayerDefn())
                out_feature.SetGeometry(segment)
                out_feature.SetField("source", "street_blocks")
                out_feature.SetField("setback", None)
                out_feature.SetField("line_id", line_id)
                output_layer.CreateFeature(out_feature)
                out_feature = None
                line_count += 1
                line_id += 1

    for feature in secondary_layer:
        poly = feature.GetGeometryRef()
        boundary = poly.Boundary()

        out_feature = ogr.Feature(output_layer.GetLayerDefn())
        out_feature.SetGeometry(boundary)
        out_feature.SetField("source", "secondary_setback")
        out_feature.SetField("setback", "secondary")
        out_feature.SetField("line_id", feature.GetFID())
        output_layer.CreateFeature(out_feature)
        out_feature = None
        line_count += 1
        line_id += 1

    ds = None
    print(f"Merged {line_count} lines into '{output_layer_name}'")


def polygonize_and_classify_blocks(
    gpkg_path,
    lines_layer_name,
    arterial_setback_layer_name,
    secondary_setback_layer_name,
    output_layer_name,
):
    """
    Create polygon blocks from merged lines and classify them by setback type.

    A block is classified based on spatial containment:
    - 'arterial_setback' if block centroid is inside arterial setback zone
    - 'secondary_setback' if block centroid is inside secondary setback zone (but not arterial)
    - 'off_grid' if block centroid is outside both setback zones

    Args:
        gpkg_path (str): Path to the GeoPackage
        lines_layer_name (str): Name of the merged lines layer
        arterial_setback_layer_name (str): Name of arterial setback polygon layer
        secondary_setback_layer_name (str): Name of secondary setback polygon layer
        output_layer_name (str): Name of the output classified blocks layer

    Returns:
        None
    """
    ds = ogr.Open(gpkg_path, 1)

    lines_layer = ds.GetLayerByName(lines_layer_name)
    arterial_layer = ds.GetLayerByName(arterial_setback_layer_name)
    secondary_layer = ds.GetLayerByName(secondary_setback_layer_name)
    srs = lines_layer.GetSpatialRef()

    # Create union of arterial setback polygons
    arterial_geoms = []
    for feature in arterial_layer:
        arterial_geoms.append(feature.GetGeometryRef().Clone())

    arterial_union = arterial_geoms[0]
    for geom in arterial_geoms[1:]:
        arterial_union = arterial_union.Union(geom)

    # Create union of secondary setback polygons
    secondary_geoms = []
    for feature in secondary_layer:
        secondary_geoms.append(feature.GetGeometryRef().Clone())

    secondary_union = secondary_geoms[0]
    for geom in secondary_geoms[1:]:
        secondary_union = secondary_union.Union(geom)

    # Collect all lines
    all_line_geoms = []
    for feature in lines_layer:
        geom = feature.GetGeometryRef().Clone()
        all_line_geoms.append(geom)

    # Create a union of all lines
    union_geom = all_line_geoms[0]
    for geom in all_line_geoms[1:]:
        union_geom = union_geom.Union(geom)

    # Polygonize
    polygons = []
    poly_geom = union_geom.Polygonize()

    if poly_geom is not None:
        geom_type = ogr.GT_Flatten(poly_geom.GetGeometryType())

        if geom_type == ogr.wkbPolygon:
            polygons.append(poly_geom.Clone())
        elif geom_type == ogr.wkbMultiPolygon or geom_type == ogr.wkbGeometryCollection:
            for i in range(poly_geom.GetGeometryCount()):
                sub_geom = poly_geom.GetGeometryRef(i)
                if ogr.GT_Flatten(sub_geom.GetGeometryType()) == ogr.wkbPolygon:
                    polygons.append(sub_geom.Clone())

    print(f"Created {len(polygons)} polygons from lines")

    # Remove existing output layer if it exists
    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    # Create output layer
    output_layer = ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon25D)

    # Add classification fields
    block_id_field = ogr.FieldDefn("block_id", ogr.OFTInteger)
    output_layer.CreateField(block_id_field)

    classification_field = ogr.FieldDefn("block_type", ogr.OFTString)
    output_layer.CreateField(classification_field)

    area_field = ogr.FieldDefn("area", ogr.OFTReal)
    output_layer.CreateField(area_field)

    # Classify each polygon
    arterial_blocks = 0
    secondary_blocks = 0
    off_grid_blocks = 0

    for block_id, poly in enumerate(polygons):
        # Get centroid for classification
        centroid = poly.Centroid()

        # Classify based on spatial containment
        if arterial_union.Contains(centroid):
            block_type = "arterial_setback"
            arterial_blocks += 1
        elif secondary_union.Contains(centroid):
            block_type = "secondary_setback"
            secondary_blocks += 1
        else:
            block_type = "off_grid"
            off_grid_blocks += 1

        # Create output features - explode multipolygons into individual polygons
        polygons = explode_multipolygon(poly)
        for polygon in polygons:
            out_feature = ogr.Feature(output_layer.GetLayerDefn())
            out_feature.SetGeometry(polygon)
            out_feature.SetField("block_id", block_id)
            out_feature.SetField("block_type", block_type)
            out_feature.SetField("area", polygon.GetArea())

            output_layer.CreateFeature(out_feature)
            out_feature = None
            block_id += 1

    ds = None

    print("Block classification complete:")
    print(f"  - Arterial setback blocks: {arterial_blocks}")
    print(f"  - Secondary setback blocks: {secondary_blocks}")
    print(f"  - Off-grid blocks: {off_grid_blocks}")


def filter_classified_blocks(
    gpkg_path: str,
    input_layer_name: str,
    output_layer_name: str,
    allowed_types=("arterial_setback", "secondary_setback"),
):
    """
    Copy only features whose block_type is in allowed_types into output_layer_name.
    """
    ds = ogr.Open(gpkg_path, 1)
    in_lyr = ds.GetLayerByName(input_layer_name)
    if in_lyr is None:
        raise RuntimeError(f"Layer not found: {input_layer_name}")

    srs = in_lyr.GetSpatialRef()
    # drop output if exists
    for i in range(ds.GetLayerCount()):
        if ds.GetLayerByIndex(i).GetName() == output_layer_name:
            ds.DeleteLayer(i)
            break

    out_lyr = ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon25D)

    # Copy all fields from input
    in_defn = in_lyr.GetLayerDefn()
    for i in range(in_defn.GetFieldCount()):
        out_lyr.CreateField(in_defn.GetFieldDefn(i).Clone())

    kept = 0
    total = 0
    in_lyr.ResetReading()
    for feat in in_lyr:
        total += 1
        if feat.GetField("block_type") in allowed_types:
            out_f = ogr.Feature(out_lyr.GetLayerDefn())
            out_f.SetGeometry(feat.GetGeometryRef().Clone())
            # copy fields
            for i in range(in_defn.GetFieldCount()):
                name = in_defn.GetFieldDefn(i).GetNameRef()
                out_f.SetField(name, feat.GetField(name))
            out_lyr.CreateFeature(out_f)
            out_f = None
            kept += 1

    ds = None
    print(f"Filtered {kept}/{total} classified blocks -> '{output_layer_name}'")


def generate_on_grid_blocks(output_path: Path, site_layer_name: str, cfg: StreetConfig) -> Path:
    print("Clip arterial setback to site boundary...")
    clip_layer(
        output_path,
        "06_arterial_setback",
        output_path,
        site_layer_name,
        output_path,
        "10a_arterial_setback_clipped",
    )

    print("Clip secondary setback to site boundary...")
    clip_layer(
        output_path,
        "07_secondary_setback",
        output_path,
        site_layer_name,
        output_path,
        "10b_secondary_setback_clipped",
    )

    print("Intersect arterial and secondary setbacks...")
    clip_layer(
        output_path,
        "10a_arterial_setback_clipped",
        output_path,
        "10b_secondary_setback_clipped",
        output_path,
        "10_intersected_setbacks",
    )

    print("Arterial setback without intersection...")
    erase_layer(
        output_path,
        "10a_arterial_setback_clipped",
        output_path,
        "10_intersected_setbacks",
        output_path,
        "11_arterial_setback_final",
    )

    print("Secondary setback without intersection...")
    erase_layer(
        output_path,
        "10b_secondary_setback_clipped",
        output_path,
        "10_intersected_setbacks",
        output_path,
        "12_secondary_setback_final",
    )

    print("Breaking arterial setback multipart features...")
    break_multipart_features(
        output_path,
        "11_arterial_setback_final",
        output_path,
        "11_arterial_setback_final",
    )

    print("Breaking secondary setback multipart features...")
    break_multipart_features(
        output_path,
        "12_secondary_setback_final",
        output_path,
        "12_secondary_setback_final",
    )

    print("Create grid from on-grid arterial setback...")
    create_grid_from_on_grid(
        output_path,
        "11_arterial_setback_final",
        "04_arterial_roads",
        cfg.off_grid_partitions_preferred_width,
        output_path,
        "11a_arterial_setback_grid",
        road_buffer_distance=cfg.road_arterial_width_m,
    )

    print("Clean up grid blocks...")
    cleanup_grid_blocks(
        output_path,
        "11a_arterial_setback_grid",
        output_path,
        "11_arterial_setback_grid_cleaned",
        cfg.on_grid_partition_depth_arterial_roads * cfg.off_grid_partitions_preferred_width * 0.5,
    )

    print("Create grid from on-grid secondary setback...")
    create_grid_from_on_grid(
        output_path,
        "12_secondary_setback_final",
        "05_secondary_roads",
        cfg.off_grid_partitions_preferred_width,
        output_path,
        "12a_secondary_setback_grid",
        intersected_setbacks_layer_name="10_intersected_setbacks",
        road_buffer_distance=cfg.road_secondary_width_m,
    )

    print("Clean up grid blocks...")
    cleanup_grid_blocks(
        output_path,
        "12a_secondary_setback_grid",
        output_path,
        "12_secondary_setback_grid_cleaned",
        (
            cfg.on_grid_partition_depth_secondary_roads
            * cfg.off_grid_partitions_preferred_width
            * 0.5
        ),
    )

    return output_path
