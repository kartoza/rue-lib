"""Core functionality for RUE."""

from rue_lib.core.exceptions import GeometryError, ProjectionError, RUEError
from rue_lib.core.geometry import bounce, clean_angles, clean_edges
from rue_lib.core.helpers import merge_gpkg_layers, remove_layer_from_gpkg
from rue_lib.core.projector import Projector

__all__ = [
    "Projector",
    "clean_edges",
    "clean_angles",
    "bounce",
    "RUEError",
    "ProjectionError",
    "GeometryError",
    "merge_gpkg_layers",
    "remove_layer_from_gpkg",
]
