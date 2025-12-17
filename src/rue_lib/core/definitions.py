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

        Amenity types - amenity allocations:
            AMENITY: Generic amenity
            SEC_AMENITY: Secondary amenity
            SEC_SEC_AMENITY: Secondary-secondary amenity
            SEC_LOC_AMENITY: Secondary-local amenity
            LOC_AMENITY: Local amenity
            LOC_LOC_AMENITY: Local-local amenity
            OFF_GRID0_AMENITY: Off-grid 0 amenity
            OFF_GRID1_AMENITY: Off-grid 1 amenity
            OFF_GRID2_AMENITY: Off-grid 2 amenity

        Open space types:
            OPEN_SPACE: Generic open space
            GREEN0: Green space type 0
            GREEN1: Green space type 1
            GREEN2: Green space type 2
            SEC_OPEN_SPACE: Secondary open space
            SEC_SEC_OPEN_SPACE: Secondary-secondary open space
            SEC_LOC_OPEN_SPACE: Secondary-local open space
            LOC_OPEN_SPACE: Local open space
            LOC_LOC_OPEN_SPACE: Local-local open space
            OFF_GRID0_OPEN_SPACE: Off-grid 0 open space
            OFF_GRID1_OPEN_SPACE: Off-grid 1 open space
            OFF_GRID2_OPEN_SPACE: Off-grid 2 open space
            CORNER_PARK: Corner park

        Path and entrance types:
            INTERNAL_PATH: Internal path
            PATH0: Path type 0
            PATH1: Path type 1
            PATH2: Path type 2
            ACCESS_PATH: Cluster access path
            ENTRANCE0: Entrance type 0
            ENTRANCE1: Entrance type 1

        Other types:
            OFF_GRID_CLUSTER_2ND: Off-grid cluster 2nd
            PUBLIC_STREETS: Public streets
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
    OFF_GRID = "off_grid"
    OFF_GRID_WARM = "off_grid0"
    OFF_GRID_COLD = "off_grid1"
    OFF_GRID_COLD2 = "off_grid2"
    CONCAVE_CORNER = "concave_corner"

    # Amenity types
    AMENITY = "am"
    ART_AMENITY = "art_am"
    SEC_AMENITY = "sec_am"
    SEC_SEC_AMENITY = "sec_sec_am"
    SEC_LOC_AMENITY = "sec_loc_am"
    LOC_AMENITY = "loc_am"
    LOC_LOC_AMENITY = "loc_loc_am"
    OFF_GRID0_AMENITY = "off_grid0_am"
    OFF_GRID1_AMENITY = "off_grid1_am"
    OFF_GRID2_AMENITY = "off_grid2_am"
    CONCAVE_CORNER_AMENITY = "concave_corner_am"

    # Green/Open space types
    GREEN0 = "green0"
    GREEN1 = "green1"
    GREEN2 = "green2"
    OPEN_SPACE = "os"
    SEC_OPEN_SPACE = "sec_os"
    SEC_SEC_OPEN_SPACE = "sec_sec_os"
    SEC_LOC_OPEN_SPACE = "sec_loc_os"
    LOC_OPEN_SPACE = "loc_os"
    LOC_LOC_OPEN_SPACE = "loc_loc_os"
    OFF_GRID0_OPEN_SPACE = "off_grid0_os"
    OFF_GRID1_OPEN_SPACE = "off_grid1_os"
    OFF_GRID2_OPEN_SPACE = "off_grid2_os"
    CONCAVE_CORNER_OPEN_SPACE = "concave_corner_os"

    # Path and entrance types
    INTERNAL_PATH = "internal path, any cluster"
    PATH0 = "path0"
    PATH1 = "path1"
    PATH2 = "path2"
    ACCESS_PATH = "Cluster, access path"
    ENTRANCE0 = "entr0"
    ENTRANCE1 = "entr1"

    # Other types
    OFF_GRID_CLUSTER_2ND = "og cluster, 2nd"
    PUBLIC_STREETS = "Public streets"


ColorTypes = {
    # Single road types
    ClusterTypes.ON_GRID_ART: "rgb(235,201,199)",
    ClusterTypes.ON_GRID_SEC: "rgb(250,232,219)",
    ClusterTypes.ON_GRID_LOC: "rgb(250,242,212)",
    # Corner types - arterial combinations
    ClusterTypes.ON_GRID_ART_ART: "rgb(204,130,122)",
    ClusterTypes.ON_GRID_ART_SEC: "rgb(217,158,153)",
    ClusterTypes.ON_GRID_ART_LOC: "rgb(230,186,184)",
    # Corner types - secondary combinations
    ClusterTypes.ON_GRID_SEC_SEC: "rgb(250,207,176)",
    ClusterTypes.ON_GRID_SEC_LOC: "rgb(250,219,199)",
    # Corner types - local combinations
    ClusterTypes.ON_GRID_LOC_LOC: "rgb(250,235,186)",
    # Off-grid types
    ClusterTypes.OFF_GRID: "rgb(207,230,242)",
    ClusterTypes.OFF_GRID_WARM: "rgb(207,230,242)",
    ClusterTypes.OFF_GRID_COLD: "rgb(207,230,242)",
    ClusterTypes.OFF_GRID_COLD2: "rgb(173,212,230)",
    ClusterTypes.CONCAVE_CORNER: "rgb(207,230,242)",
    # Amenity types
    ClusterTypes.AMENITY: "rgb(217,153,168)",
    ClusterTypes.ART_AMENITY: "rgb(217,153,168)",
    ClusterTypes.SEC_AMENITY: "rgb(217,153,168)",
    ClusterTypes.SEC_SEC_AMENITY: "rgb(217,153,168)",
    ClusterTypes.SEC_LOC_AMENITY: "rgb(217,153,168)",
    ClusterTypes.LOC_AMENITY: "rgb(217,153,168)",
    ClusterTypes.LOC_LOC_AMENITY: "rgb(217,153,168)",
    ClusterTypes.OFF_GRID0_AMENITY: "rgb(217,153,168)",
    ClusterTypes.OFF_GRID1_AMENITY: "rgb(217,153,168)",
    ClusterTypes.OFF_GRID2_AMENITY: "rgb(217,153,168)",
    ClusterTypes.CONCAVE_CORNER_AMENITY: "rgb(217,153,168)",
    # Green/Open space types
    ClusterTypes.GREEN0: "rgb(219,242,207)",
    ClusterTypes.GREEN1: "rgb(219,242,207)",
    ClusterTypes.GREEN2: "rgb(219,242,207)",
    ClusterTypes.OPEN_SPACE: "rgb(191,230,171)",
    ClusterTypes.SEC_OPEN_SPACE: "rgb(191,230,171)",
    ClusterTypes.SEC_SEC_OPEN_SPACE: "rgb(191,230,171)",
    ClusterTypes.SEC_LOC_OPEN_SPACE: "rgb(191,230,171)",
    ClusterTypes.LOC_OPEN_SPACE: "rgb(191,230,171)",
    ClusterTypes.LOC_LOC_OPEN_SPACE: "rgb(191,230,171)",
    ClusterTypes.OFF_GRID0_OPEN_SPACE: "rgb(191,230,171)",
    ClusterTypes.OFF_GRID1_OPEN_SPACE: "rgb(191,230,171)",
    ClusterTypes.OFF_GRID2_OPEN_SPACE: "rgb(191,230,171)",
    ClusterTypes.CONCAVE_CORNER_OPEN_SPACE: "rgb(191,230,171)",
    # Path and entrance types
    ClusterTypes.INTERNAL_PATH: "rgb(237,255,230)",
    ClusterTypes.PATH0: "rgb(237,255,230)",
    ClusterTypes.PATH1: "rgb(237,255,230)",
    ClusterTypes.PATH2: "rgb(237,255,230)",
    ClusterTypes.ACCESS_PATH: "rgb(255,250,235)",
    ClusterTypes.ENTRANCE0: "rgb(255,250,235)",
    ClusterTypes.ENTRANCE1: "rgb(255,250,235)",
    # Other types
    ClusterTypes.OFF_GRID_CLUSTER_2ND: "rgb(173,212,230)",
    ClusterTypes.PUBLIC_STREETS: "rgb(255,255,255)",
}
