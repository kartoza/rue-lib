"""Core definitions for rue-lib.

This module contains enumeration types and constants used throughout the rue-lib library
for defining properties and types related to street and grid generation.

Classes:
    Properties: Enum for property field names used in GIS layers.
    RoadTypes: Enum for road type classifications.
"""
from enum import EnumType


class PropertyNames(EnumType):
    """Property field names used in GIS layers.

    Attributes:
        RoadType: Field name for road type classification.
        GridType: Field name for grid type classification.
    """
    RoadType = "road_type"
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
