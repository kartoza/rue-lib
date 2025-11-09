"""Projection utilities for converting between WGS84 and local metric CRS."""

from __future__ import annotations

from dataclasses import dataclass

from pyproj import CRS, Transformer
from shapely.geometry.base import BaseGeometry
from shapely.ops import transform as shp_transform

from rue_lib.core.exceptions import ProjectionError


@dataclass
class Projector:
    """Handles coordinate transformations between WGS84 and local UTM CRS."""

    fwd: Transformer
    rev: Transformer

    @staticmethod
    def from_geom(geom: BaseGeometry) -> Projector:
        """
        Create a Projector from a geometry's centroid.

        Args:
            geom: A Shapely geometry in WGS84 (EPSG:4326)

        Returns:
            A Projector configured for the appropriate UTM zone

        Raises:
            ProjectionError: If projection setup fails
        """
        try:
            cx, cy = geom.centroid.x, geom.centroid.y
            zone = int((cx + 180.0) // 6.0) + 1
            epsg = 32600 + zone if cy >= 0 else 32700 + zone

            wgs84 = CRS.from_epsg(4326)
            local = CRS.from_epsg(epsg)

            return Projector(
                Transformer.from_crs(wgs84, local, always_xy=True),
                Transformer.from_crs(local, wgs84, always_xy=True),
            )
        except Exception as e:
            raise ProjectionError(f"Failed to create projector: {e}") from e

    def to_local(self, geom: BaseGeometry) -> BaseGeometry:
        """Transform geometry from WGS84 to local metric CRS."""
        return shp_transform(self.fwd.transform, geom)

    def to_wgs84(self, geom: BaseGeometry) -> BaseGeometry:
        """Transform geometry from local metric CRS to WGS84."""
        return shp_transform(self.rev.transform, geom)
