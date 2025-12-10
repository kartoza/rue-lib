# src/rue_lib/cluster/frame.py
"""Create frame geometries from blocks by subtracting off-grid areas."""

from pathlib import Path

import geopandas as gpd


def extract_frame(
    output_path: Path,
    off_grid_layer_name: str,
    off_grids_inside_layer_name: str,
    output_layer_name: str,
):
    """Extract frame geometries by subtracting off-grid areas from blocks.

    Creates frame polygons representing the buildable perimeter area around
    off-grid spaces within blocks. The frame is computed as the geometric
    difference between the original block and its interior off-grid area.

    Args:
        output_path: Path to the GeoPackage file containing input layers and
            where output will be saved.
        off_grid_layer_name: Name of the layer containing block geometries.
        off_grids_inside_layer_name: Name of the layer containing off-grid
            geometries with their associated block_id values.
        output_layer_name: Name for the output layer that will contain the
            extracted frames.

    Returns:
        str: The name of the output layer.

    Notes:
        - Each frame includes metadata: block_id, frame_area, original_area,
          and off_grid_area.
        - Empty frames (where off-grid equals entire block) are skipped.
        - Errors during frame creation are caught and logged.
    """
    off_grid_layer = gpd.read_file(output_path, layer=off_grid_layer_name)
    off_grids_inside_layer = gpd.read_file(output_path, layer=off_grids_inside_layer_name)
    frames = []
    for idx, block_row in off_grid_layer.iterrows():
        block = block_row.geometry
        block_id = idx

        # Find corresponding off-grid
        try:
            off_grid_row = off_grids_inside_layer[
                off_grids_inside_layer["block_id"] == block_id
            ].iloc[0]
        except IndexError:
            continue

        off_grid = off_grid_row.geometry

        # Create frame by subtracting off-grid from block
        try:
            frame = block.difference(off_grid)

            if not frame.is_empty:
                frame_area = frame.area
                print(f"    Block {block_id}: Frame area = {frame_area:.2f} m²")
                frames.append(
                    {
                        "geometry": frame,
                        "block_id": block_id,
                        "frame_area": frame_area,
                        "original_area": block.area,
                        "off_grid_area": off_grid.area,
                    }
                )
            else:
                print(f"    {block_id}: ✗ Empty frame (off-grid = entire block)")
        except Exception as e:
            print(f"    Block {block_id}: ✗ Error creating frame: {e}")

    gdf_out = gpd.GeoDataFrame(frames, crs=off_grid_layer.crs)
    gdf_out.to_file(output_path, layer=output_layer_name, driver="GPKG")
    return output_layer_name
