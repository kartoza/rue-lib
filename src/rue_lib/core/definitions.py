"""Core definitions for rue-lib.

This module contains enumeration types and constants used throughout the rue-lib library
for defining properties and types related to street and grid generation.

Classes:
    Properties: Enum for property field names used in GIS layers.
    RoadTypes: Enum for road type classifications.
"""

from enum import EnumType


class PropertyKeys(EnumType):
    """Property field names used in GIS layers.

    Attributes:
        RoadType: Field name for road type classification.
        GridType: Field name for grid type classification.
    """

    RoadType = "type"
    GridType = "type"
    BlockType = "type"


class RoadTypes(EnumType):
    """Road type classifications used in the street network.

    Attributes:
        Artery: Arterial road type identifier.
        Secondary: Secondary road type identifier.
        Local: Local road type identifier.
    """

    Artery = "road_arterial"
    Secondary = "road_secondary"
    Tertiary = "road_tertiary"
    Local = "road_local"


class BlockTypes(EnumType):
    """Block type classifications used in grid clustering.

    Attributes:
        ON_GRID_ART: On-grid arterial block type.
        ON_GRID_SEC: On-grid secondary block type.
        OFF_GRID: Off-grid block type.
        COLD_GRID: Cold boundary block type.
    """

    ON_GRID_ART = "on_grid_art"
    ON_GRID_SEC = "on_grid_sec"
    OFF_GRID = "off_grid"
    COLD_GRID = "cold_boundaries"


class ClusterTypes(EnumType):
    """Cluster type classifications used in grid clustering.

    Attributes:
        Single road types - blocks adjacent to one road type:
            ON_GRID_ART: Arterial road block
            ON_GRID_SEC: Secondary road block
            ON_GRID_LOC: Local road block

        Corner types - blocks at corners between two road types:
            ON_GRID_ART_ART: Corner with arterial and arterial
            ON_GRID_ART_SEC: Corner with arterial and secondary
            ON_GRID_ART_LOC: Corner with arterial and local
            ON_GRID_SEC_SEC: Corner with secondary and secondary
            ON_GRID_SEC_LOC: Corner with secondary and local
            ON_GRID_LOC_LOC: Corner with local and local

        Off-grid types - blocks not aligned with grid:
            OFF_GRID_WARM: Warm off-grid block (interior)
            OFF_GRID_COLD: Cold off-grid block (boundary)
            OFF_GRID_COLD2: Secondary cold off-grid block
            CONCAVE_CORNER: Concave corner block
    """

    # Single road types
    ON_GRID_ART = "art"
    ON_GRID_SEC = "sec"
    ON_GRID_LOC = "loc"

    # Corner types - arterial combinations
    ON_GRID_ART_ART = "art_art"
    ON_GRID_ART_SEC = "art_sec"
    ON_GRID_ART_LOC = "art_loc"

    # Corner types - secondary combinations
    ON_GRID_SEC_SEC = "sec_sec"
    ON_GRID_SEC_LOC = "sec_loc"

    # Corner types - local combinations
    ON_GRID_LOC_LOC = "loc_loc"

    # Off-grid types
    OFF_GRID_WARM = "off_grid0"
    OFF_GRID_COLD = "off_grid1"
    OFF_GRID_COLD2 = "off_grid2"
    CONCAVE_CORNER = "concave_corner"


ColorTypes = {
    # Single road types
    ClusterTypes.ON_GRID_ART: "rgb(235,201,199)",
    ClusterTypes.ON_GRID_SEC: "rgb(250,232,219)",
    ClusterTypes.ON_GRID_LOC: "rgb(250,242,212)",
    # Corner types - arterial combinations
    ClusterTypes.ON_GRID_ART_ART: "rgb(235,201,199)",
    ClusterTypes.ON_GRID_ART_SEC: "rgb(217,158,153)",
    ClusterTypes.ON_GRID_ART_LOC: "rgb(230,186,184)",
    # Corner types - secondary combinations
    ClusterTypes.ON_GRID_SEC_SEC: "rgb(250,232,219)",
    ClusterTypes.ON_GRID_SEC_LOC: "rgb(250,219,199)",
    # Corner types - local combinations
    ClusterTypes.ON_GRID_LOC_LOC: "rgb(250,235,186)",
    # Off-grid types
    ClusterTypes.OFF_GRID_WARM: "rgb(207,230,242)",
    ClusterTypes.OFF_GRID_COLD: "rgb(207,230,242)",
    ClusterTypes.OFF_GRID_COLD2: "rgb(173,212,230)",
    ClusterTypes.CONCAVE_CORNER: "rgb(207,230,242)",
}
