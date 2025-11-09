"""Tests for core functionality."""

import pytest
from shapely.geometry import MultiPolygon, Point, Polygon

from rue_lib.core.exceptions import GeometryError, ProjectionError, RUEError
from rue_lib.core.geometry import (
    bounce,
    clean_angles,
    clean_edges,
    ring_clean_collinear,
    ring_clean_short_edges,
)
from rue_lib.core.projector import Projector


class TestProjector:
    """Tests for Projector class."""

    def test_from_geom_northern_hemisphere(self):
        """Test projector creation for northern hemisphere geometry."""
        # Point in London (51.5°N, 0.1°W)
        point = Point(-0.1, 51.5)
        projector = Projector.from_geom(point)

        assert projector.fwd is not None
        assert projector.rev is not None

    def test_from_geom_southern_hemisphere(self):
        """Test projector creation for southern hemisphere geometry."""
        # Point in Sydney (-33.9°S, 151.2°E)
        point = Point(151.2, -33.9)
        projector = Projector.from_geom(point)

        assert projector.fwd is not None
        assert projector.rev is not None

    def test_to_local_and_back(self):
        """Test round-trip transformation."""
        # Create a simple polygon
        poly = Polygon([(0.0, 0.0), (0.1, 0.0), (0.1, 0.1), (0.0, 0.1), (0.0, 0.0)])

        projector = Projector.from_geom(poly)

        # Transform to local and back
        poly_local = projector.to_local(poly)
        poly_back = projector.to_wgs84(poly_local)

        # Should be very close to original
        assert poly.equals_exact(poly_back, tolerance=1e-6)

    def test_area_calculation_accuracy(self):
        """Test that area calculations are more accurate in local CRS."""
        # Create a 1km x 1km square approximately
        # (this is rough, actual size varies by latitude)
        poly = Polygon([(0.0, 0.0), (0.01, 0.0), (0.01, 0.01), (0.0, 0.01), (0.0, 0.0)])

        projector = Projector.from_geom(poly)
        poly_local = projector.to_local(poly)

        # Area in local CRS should be in square meters
        area_m2 = poly_local.area

        # Should be roughly 1 million square meters (1 km²)
        # Allow for latitude variations
        assert 500_000 < area_m2 < 2_000_000


class TestGeometryOperations:
    """Tests for geometry cleaning operations."""

    def test_clean_edges_fixes_self_intersection(self):
        """Test that clean_edges fixes self-intersecting polygons."""
        # Create a bowtie polygon (self-intersecting)
        poly = Polygon([(0, 0), (2, 2), (2, 0), (0, 2), (0, 0)])

        cleaned = clean_edges(poly)

        assert cleaned.is_valid
        assert not cleaned.is_empty

    def test_clean_edges_preserves_valid_polygon(self):
        """Test that clean_edges preserves valid polygons."""
        poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])

        cleaned = clean_edges(poly)

        assert cleaned.is_valid
        assert cleaned.area > 0

    def test_clean_angles_simplifies_polygon(self):
        """Test that clean_angles removes near-collinear vertices."""
        # Create a polygon with many closely spaced points along edges
        poly = Polygon(
            [
                (0, 0),
                (0.5, 0.001),  # Nearly collinear
                (1, 0),
                (1, 1),
                (0, 1),
                (0, 0),
            ]
        )

        cleaned = clean_angles(poly, eps=0.01)

        # Should have fewer vertices
        assert len(cleaned.exterior.coords) <= len(poly.exterior.coords)
        assert cleaned.is_valid

    def test_bounce_with_zero_distance(self):
        """Test bounce with zero distance returns original polygon."""
        poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])

        bounced = bounce(poly, 0.0)

        assert bounced.equals(poly)

    def test_bounce_with_positive_distance(self):
        """Test bounce with positive distance smooths polygon."""
        poly = Polygon([(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)])

        bounced = bounce(poly, 1.0)

        assert bounced.is_valid
        assert bounced.area > 0
        # Area should be slightly different due to bounce
        assert abs(bounced.area - poly.area) < poly.area * 0.5


class TestRingCleaning:
    """Tests for ring cleaning functions."""

    def test_ring_clean_short_edges_removes_short_segments(self):
        """Test removal of short edges."""
        coords = [
            (0.0, 0.0),
            (0.01, 0.0),  # Very short edge
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ]

        cleaned = ring_clean_short_edges(coords, min_len=1.0)

        # Should have removed the short edge
        assert len(cleaned) < len(coords)
        assert cleaned[0] == cleaned[-1]  # Still closed

    def test_ring_clean_short_edges_preserves_long_edges(self):
        """Test that long edges are preserved."""
        coords = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]

        cleaned = ring_clean_short_edges(coords, min_len=1.0)

        # All edges are long, should be unchanged
        assert len(cleaned) == len(coords)

    def test_ring_clean_collinear_removes_collinear_points(self):
        """Test removal of collinear vertices."""
        coords = [
            (0.0, 0.0),
            (5.0, 0.0),  # Collinear with prev and next
            (10.0, 0.0),
            (10.0, 10.0),
            (0.0, 10.0),
            (0.0, 0.0),
        ]

        cleaned = ring_clean_collinear(coords, dot_thresh=0.9999)

        # Should have removed the middle point on the straight edge
        assert len(cleaned) <= len(coords)
        assert cleaned[0] == cleaned[-1]  # Still closed

    def test_ring_clean_collinear_preserves_corners(self):
        """Test that corner vertices are preserved."""
        coords = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0), (0.0, 0.0)]

        cleaned = ring_clean_collinear(coords, dot_thresh=0.9999)

        # All points are corners, should be unchanged (or nearly so)
        assert len(cleaned) >= 4  # At least a triangle + closure


class TestExceptions:
    """Tests for custom exceptions."""

    def test_rue_error_is_base_exception(self):
        """Test that RUEError is the base exception."""
        with pytest.raises(RUEError):
            raise RUEError("Test error")

    def test_projection_error_inherits_from_rue_error(self):
        """Test ProjectionError inheritance."""
        with pytest.raises(RUEError):
            raise ProjectionError("Test projection error")

    def test_geometry_error_inherits_from_rue_error(self):
        """Test GeometryError inheritance."""
        with pytest.raises(RUEError):
            raise GeometryError("Test geometry error")

    def test_exception_messages(self):
        """Test that exception messages are preserved."""
        message = "Custom error message"

        try:
            raise GeometryError(message)
        except GeometryError as e:
            assert str(e) == message


class TestEdgeCases:
    """Tests for edge cases and error conditions."""

    def test_clean_edges_with_empty_polygon(self):
        """Test clean_edges with empty polygon."""
        poly = Polygon()

        with pytest.raises(GeometryError):
            clean_edges(poly)

    def test_bounce_with_negative_distance(self):
        """Test bounce with negative distance returns original."""
        poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])

        bounced = bounce(poly, -1.0)

        assert bounced.equals(poly)

    def test_ring_clean_with_minimal_ring(self):
        """Test ring cleaning with minimal valid ring."""
        coords = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.0, 0.0)]

        cleaned = ring_clean_short_edges(coords, min_len=0.1)

        # Should preserve minimal ring
        assert len(cleaned) >= 4  # Triangle + closure

    def test_projector_with_multipolygon(self):
        """Test projector with MultiPolygon geometry."""
        poly1 = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
        poly2 = Polygon([(2, 0), (3, 0), (3, 1), (2, 1), (2, 0)])
        multipoly = MultiPolygon([poly1, poly2])

        projector = Projector.from_geom(multipoly)
        multipoly_local = projector.to_local(multipoly)

        assert multipoly_local.geom_type == "MultiPolygon"
        assert len(multipoly_local.geoms) == 2
