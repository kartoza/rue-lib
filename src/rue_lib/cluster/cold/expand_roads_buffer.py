from osgeo import ogr
from shapely import wkt

from rue_lib.cluster.config import ClusterConfig


def create_buffered_lines_from_boundary_lines(
    input_gpkg: str,
    boundary_lines_layer_name: str,
    _erased_grid_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
    cfg: ClusterConfig,
) -> str:
    """
    Buffer road-adjacent boundary lines to create road strips.

    Uses the boundary lines layer (e.g. 202_cold_boundary_lines) instead of points
    to generate rectangular buffers per road type.
    """
    # Read first, then close before writing to avoid overlapping SQLite transactions
    read_ds = ogr.Open(input_gpkg, 0)
    if read_ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    lines_layer = read_ds.GetLayerByName(boundary_lines_layer_name)
    if lines_layer is None:
        raise ValueError(f"Layer {boundary_lines_layer_name} not found")

    srs = lines_layer.GetSpatialRef()

    buffer_depths = {
        "road_local": cfg.on_grid_partition_depth_local_roads * 2,
        "road_secondary": cfg.on_grid_partition_depth_secondary_roads * 2,
        "road_arterial": cfg.on_grid_partition_depth_arterial_roads * 2,
    }

    buffered_lines = []
    for feat in lines_layer:
        road_type = feat.GetField("road_type")
        block_id = feat.GetField("block_id")
        geom = feat.GetGeometryRef()
        if geom is None:
            continue

        depth = buffer_depths.get(road_type, 0)
        if depth <= 0:
            continue

        shapely_geom = wkt.loads(str(geom))
        shapely_buffered = shapely_geom.buffer(depth / 2.0, join_style=2, cap_style=3)
        buffered_geom = ogr.CreateGeometryFromWkb(shapely_buffered.wkb)

        if buffered_geom and not buffered_geom.IsEmpty():
            buffered_lines.append(
                {
                    "geometry": buffered_geom.Clone(),
                    "block_id": block_id,
                    "road_type": road_type,
                }
            )

    lines_layer = None
    read_ds = None

    out_ds = ogr.Open(output_gpkg, 1)
    if out_ds is None:
        raise ValueError(f"Could not open {output_gpkg} for writing")

    for i in range(out_ds.GetLayerCount()):
        layer = out_ds.GetLayerByIndex(i)
        if layer.GetName() == output_layer_name:
            out_ds.DeleteLayer(i)
            break

    out_layer = out_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)
    out_layer.CreateField(ogr.FieldDefn("block_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("road_type", ogr.OFTString))

    for line_data in buffered_lines:
        out_feat = ogr.Feature(out_layer.GetLayerDefn())
        out_feat.SetGeometry(line_data["geometry"])
        out_feat.SetField("block_id", line_data["block_id"])
        out_feat.SetField("road_type", line_data["road_type"])

        out_layer.CreateFeature(out_feat)
        out_feat = None

    out_layer = None
    out_ds = None

    print(f"  Created buffered lines layer: {output_layer_name}")
    return output_layer_name


def clip_buffered_lines_to_cold_grid(
    input_gpkg: str,
    buffered_lines_layer_name: str,
    cold_grid_layer_name: str,
    output_gpkg: str,
    output_layer_name: str,
) -> str:
    """
    Clip buffered lines to only keep geometries inside cold grid blocks.

    Args:
        input_gpkg: Path to input GeoPackage
        buffered_lines_layer_name: Name of buffered lines layer (e.g., 206_buffered_lines)
        cold_grid_layer_name: Name of cold grid erased layer (e.g., 201_cold_grid_erased)
        output_gpkg: Path to output GeoPackage
        output_layer_name: Name for output clipped layer

    Returns:
        Name of the output layer
    """
    # Read input data
    read_ds = ogr.Open(input_gpkg, 0)
    if read_ds is None:
        raise ValueError(f"Could not open {input_gpkg}")

    buffered_layer = read_ds.GetLayerByName(buffered_lines_layer_name)
    if buffered_layer is None:
        raise ValueError(f"Layer {buffered_lines_layer_name} not found")

    cold_grid_layer = read_ds.GetLayerByName(cold_grid_layer_name)
    if cold_grid_layer is None:
        raise ValueError(f"Layer {cold_grid_layer_name} not found")

    srs = buffered_layer.GetSpatialRef()

    # Collect cold grid geometries by block_id
    cold_grid_geoms = {}
    for feat in cold_grid_layer:
        geom = feat.GetGeometryRef()
        if geom:
            # Try to get block_id or use fid
            block_id = (
                feat.GetField("block_id") if feat.GetFieldIndex("block_id") >= 0 else feat.GetFID()
            )
            cold_grid_geoms[block_id] = geom.Clone()

    # Process buffered lines
    clipped_geometries = []
    for feat in buffered_layer:
        buffered_geom = feat.GetGeometryRef()
        block_id = feat.GetField("block_id")
        road_type = feat.GetField("road_type")

        if buffered_geom is None:
            continue

        # Find the corresponding cold grid block
        cold_grid_geom = cold_grid_geoms.get(block_id)
        if cold_grid_geom is None:
            continue

        # Intersect buffered line with cold grid block
        try:
            clipped_geom = buffered_geom.Intersection(cold_grid_geom)
            if clipped_geom and not clipped_geom.IsEmpty():
                clipped_geometries.append(
                    {
                        "geometry": clipped_geom.Clone(),
                        "block_id": block_id,
                        "road_type": road_type,
                    }
                )
        except RuntimeError as e:
            print(f"    Warning: Failed to clip geometry for block {block_id}: {e}")
            continue

    read_ds = None

    # Write output
    out_ds = ogr.Open(output_gpkg, 1)
    if out_ds is None:
        raise ValueError(f"Could not open {output_gpkg} for writing")

    # Delete existing layer if it exists
    for i in range(out_ds.GetLayerCount()):
        layer = out_ds.GetLayerByIndex(i)
        if layer.GetName() == output_layer_name:
            out_ds.DeleteLayer(i)
            break

    # Create output layer
    out_layer = out_ds.CreateLayer(output_layer_name, srs, ogr.wkbPolygon)
    out_layer.CreateField(ogr.FieldDefn("block_id", ogr.OFTInteger))
    out_layer.CreateField(ogr.FieldDefn("road_type", ogr.OFTString))

    for geom_data in clipped_geometries:
        out_feat = ogr.Feature(out_layer.GetLayerDefn())
        out_feat.SetGeometry(geom_data["geometry"])
        out_feat.SetField("block_id", geom_data["block_id"])
        out_feat.SetField("road_type", geom_data["road_type"])

        out_layer.CreateFeature(out_feat)
        out_feat = None

    out_layer = None
    out_ds = None

    print(f"  Created clipped layer: {output_layer_name} ({len(clipped_geometries)} features)")
    return output_layer_name
