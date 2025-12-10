import geopandas as gpd
from shapely.geometry import Polygon

from rue_lib.cluster.helpers import get_edge_angle
from rue_lib.core.definitions import RoadTypes


def classify_part_type(
        part: Polygon,
        part_edges_gdf: gpd.GeoDataFrame,
        part_art_d: float,
        part_sec_d: float,
        part_loc_d: float,
) -> str:
    """
    Classify part type based on adjacent road types and geometry.

    This function analyzes the edges of a block part to determine its type based on
    which road types it borders. It identifies corner parts (adjacent to multiple
    road types) and side parts (adjacent to one road type). Corner parts that are
    too large are automatically reclassified as side parts.

    The classification process:
    1. Groups consecutive edges by road type to identify distinct road segments
    2. Counts road segments using angle analysis at vertices
    3. Determines if part is a corner (multiple road types) or side (single type)
    4. Validates corner parts against expected area thresholds
    5. Simplifies oversized corners to side parts

    Args:
        part: Part polygon to classify
        part_edges_gdf: GeoDataFrame with part edges and road_type attributes
        part_art_d: Depth/offset for arterial roads in meters
        part_sec_d: Depth/offset for secondary roads in meters
        part_loc_d: Depth/offset for local roads in meters

    Returns:
        Part type string indicating road adjacency:
            - Single road types: 'art', 'sec', 'loc'
            - Corner types: 'art_art', 'art_sec', 'art_loc', 'sec_sec', 'sec_loc', 'loc_loc'
            - Off-grid: 'off_grid' (no adjacent roads)

    Example:
        >>> part_edges = extract_block_edges(part_gdf, roads_gdf)
        >>> part_type = classify_part_type(part, part_edges, 40.0, 30.0, 20.0)
        >>> print(part_type)  # e.g., 'art_sec' for arterial-secondary corner
    """
    areas_dict = {
        'art_art': part_art_d * part_art_d,
        'art_sec': part_art_d * part_sec_d,
        'art_loc': part_art_d * part_loc_d,
        'sec_sec': part_sec_d * part_sec_d,
        'sec_loc': part_sec_d * part_loc_d,
    }
    road_types = []

    # Check each road type
    for road_type_check in [
        RoadTypes.Artery, RoadTypes.Secondary, RoadTypes.Local
    ]:
        # Get edges matching this road type
        matching_edges = part_edges_gdf[
            part_edges_gdf.get('road_type', '') == road_type_check
            ]
        if matching_edges.empty:
            continue

        # Group consecutive edges
        edge_groups = [[]]
        prev_found = False

        for _, edge in part_edges_gdf.iterrows():
            if edge.get('road_type') == road_type_check:
                edge_groups[-1].append(edge)
                prev_found = True
            else:
                if prev_found and len(edge_groups[-1]) > 0:
                    edge_groups.append([])
                prev_found = False

        # Remove empty groups
        edge_groups = [g for g in edge_groups if len(g) > 0]

        if not edge_groups:
            continue

        # Merge first and last if they connect
        if len(edge_groups) == 2:
            edges = edge_groups[1] + edge_groups[0]
        else:
            edges = edge_groups[0] if edge_groups else []

        if len(edges) == 0:
            continue

        if len(edges) == 1:
            road_types.append(road_type_check)
            continue

        # Check angles at vertices to count distinct road segments
        coords = list(part.exterior.coords)

        # For multiple edges, check angles to see if they're separate segments
        segment_count = 1
        for i in range(1, len(edges)):
            # Get vertex angle
            try:
                angle = get_edge_angle(coords, i)
                if abs(angle - 180) > 45:
                    segment_count += 1
            except Exception:
                pass

        road_types.extend([road_type_check] * segment_count)

    road_type = "off_grid"

    # Map to block type
    block_types_dict = {
        RoadTypes.Artery: 'art',
        RoadTypes.Secondary: 'sec',
        RoadTypes.Tertiary: 'ter',
        RoadTypes.Local: 'loc'
    }
    if len(road_types) == 0:
        road_type = 'off_grid'
    elif len(road_types) == 1:
        road_type = block_types_dict.get(road_types[0], 'off_grid')
    else:
        # Multiple road types
        type0 = block_types_dict.get(road_types[0], '')
        type1 = block_types_dict.get(road_types[1], '')
        road_type = f"{type0}_{type1}" if type0 and type1 else type0

    part_type = road_type

    # Check if corner part should be simplified to side part
    corner_expected_area = areas_dict.get(road_type)
    if corner_expected_area is not None:
        if part.area > (corner_expected_area * 3):
            # It's too large for a corner, make it a side part
            part_type = road_type[:3]

        # Check for local road length
    if len(road_type) == 7:  # e.g., 'art_loc'
        # Check local road edge lengths
        loc_edges = part_edges_gdf[
            part_edges_gdf.get('road_type') == RoadTypes.Local]
        if not loc_edges.empty:
            max_length = loc_edges.geometry.length.max()
            if max_length > (part_loc_d * 2):
                part_type = road_type[:3]

        # Check secondary road edge lengths
        sec_edges = part_edges_gdf[
            part_edges_gdf.get('road_type') == RoadTypes.Secondary]
        if not sec_edges.empty:
            max_length = sec_edges.geometry.length.max()
            if max_length > (part_sec_d * 2):
                part_type = road_type[:3]

    return part_type


def classify_plot_by_area(
        plot: Polygon,
        part_og_w: float,
        part_og_d: float,
        threshold: float = 0.3,
) -> str:
    """
    Classify a plot as 'plot' or 'park' based on area comparison.

    Determines whether a plot is large enough to be a buildable plot or
    should be classified as park/open space based on a percentage of the
    target plot area.

    Args:
        plot: Plot polygon to classify
        part_og_w: Target plot width (meters)
        part_og_d: Target plot depth (meters)
        threshold: Area threshold as fraction of target (0.3 = 30% of target)

    Returns:
        'plot' if area >= threshold * target_area, otherwise 'park'
    """
    target_area = part_og_w * part_og_d

    if plot.area >= target_area * threshold:
        return "plot"
    else:
        return "park"
