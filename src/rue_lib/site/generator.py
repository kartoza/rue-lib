"""Grid generation and parcel creation."""

import geopandas as gpd
from shapely.geometry import box


def make_grid_over_bbox(
    bbox: tuple[float, float, float, float], crs, rows: int = 4, cols: int = 4, pad: float = 10.0
) -> gpd.GeoDataFrame:
    """
    Create a rectangular grid over the bounding box.

    Args:
        bbox: Bounding box (minx, miny, maxx, maxy)
        crs: Coordinate reference system
        rows: Number of rows in the grid
        cols: Number of columns in the grid
        pad: Padding to add around bbox in meters

    Returns:
        GeoDataFrame with grid cells
    """
    minx, miny, maxx, maxy = bbox
    minx -= pad
    miny -= pad
    maxx += pad
    maxy += pad

    dx = (maxx - minx) / cols
    dy = (maxy - miny) / rows

    cells = []
    cell_ids = []
    idx = 0

    for r in range(rows):
        for c in range(cols):
            x0 = minx + c * dx
            y0 = miny + r * dy
            x1 = x0 + dx
            y1 = y0 + dy
            cells.append(box(x0, y0, x1, y1))
            cell_ids.append(idx)
            idx += 1

    return gpd.GeoDataFrame({"cell_id": cell_ids}, geometry=cells, crs=crs)


def create_parcel_grid(
    site_gdf: gpd.GeoDataFrame, rows: int = 4, cols: int = 4, pad: float = 10.0
) -> gpd.GeoDataFrame:
    """
    Create parcels by intersecting grid with site polygons.

    Args:
        site_gdf: Site geometry in metric CRS
        rows: Number of grid rows
        cols: Number of grid columns
        pad: Grid padding in meters

    Returns:
        GeoDataFrame with parcel polygons
    """
    bbox = site_gdf.total_bounds
    grid = make_grid_over_bbox(tuple(bbox), crs=site_gdf.crs, rows=rows, cols=cols, pad=pad)

    parcels = gpd.overlay(grid, site_gdf, how="intersection")

    if parcels.empty:
        return parcels

    parcels["geometry"] = parcels.buffer(0)

    # Dissolve by cell to maintain contiguous parcels
    parcels = parcels.dissolve(by="cell_id", as_index=False, aggfunc="sum")
    parcels = parcels.explode(index_parts=False, ignore_index=True)

    parcels["area_m2"] = parcels.geometry.area

    # Remove tiny slivers only
    parcels = parcels[parcels["area_m2"] >= 5.0].copy()
    parcels = parcels[["cell_id", "area_m2", "geometry"]]

    return parcels
