from shapely import wkt


def feature_geom_to_shapely(feat):
    """
    Convert a feature geometry (QGIS or OGR) to Shapely.
    """
    if hasattr(feat, "geometry"):
        geom = feat.geometry()
        wkt_str = geom.asWkt() if hasattr(geom, "asWkt") else geom.ExportToWkt()
    else:
        geom = feat.GetGeometryRef()
        wkt_str = geom.ExportToWkt()
    return wkt.loads(wkt_str)
