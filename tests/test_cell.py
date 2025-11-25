#!/usr/bin/env python3
"""Comprehensive tests for cell.py module."""

import pytest
from shapely.geometry import Polygon

from rue_lib.streets.cell import (
    Cell,
    find_neighbor_corner_split_lines,
    fix_area_too_large,
    get_neighbors,
    inspect_and_fix_cells,
)
from rue_lib.streets.grids import is_good_cell


class TestCell:
    """Tests for the Cell class."""

    def test_cell_creation(self):
        """Test creating a Cell instance."""
        geom = Polygon([[0, 0], [100, 0], [100, 100], [0, 100], [0, 0]])
        quality = {"is_good": True, "reason": "good"}
        cell = Cell(id=1, geom=geom, quality=quality)

        assert cell.id == 1
        assert cell.geom == geom
        assert cell.quality == quality

    def test_cell_attributes(self):
        """Test Cell attributes are correctly set."""
        geom = Polygon([[0, 0], [50, 0], [50, 50], [0, 50], [0, 0]])
        quality = {
            "is_good": False,
            "reason": "area_too_small",
            "area_ratio": 0.25,
            "num_vertices": 4,
            "right_angles": 4,
        }
        cell = Cell(id=42, geom=geom, quality=quality)

        assert cell.id == 42
        assert cell.geom.area == 2500.0
        assert cell.quality["is_good"] is False
        assert cell.quality["reason"] == "area_too_small"


class TestGetNeighbors:
    """Tests for get_neighbors function."""

    def test_get_neighbors_3x3_grid(self):
        """Test finding neighbors in a 3x3 grid."""
        target_area = 10000.0  # 100m × 100m

        # Create 3x3 grid
        cells = []
        for row in range(3):
            for col in range(3):
                x = col * 100
                y = row * 100
                geom = Polygon([[x, y], [x + 100, y], [x + 100, y + 100], [x, y + 100], [x, y]])
                quality = is_good_cell(geom, target_area)
                cell = Cell(id=len(cells), geom=geom, quality=quality)
                cells.append(cell)

        # Test center cell (id=4) - should have 8 neighbors
        center_cell = cells[4]
        neighbors = get_neighbors(center_cell, cells)
        assert len(neighbors) == 8

        # Test corner cell (id=0) - should have 3 neighbors
        corner_cell = cells[0]
        neighbors = get_neighbors(corner_cell, cells)
        assert len(neighbors) == 3

        # Test edge cell (id=1) - should have 5 neighbors
        edge_cell = cells[1]
        neighbors = get_neighbors(edge_cell, cells)
        assert len(neighbors) == 5

    def test_get_neighbors_excludes_self(self):
        """Test that get_neighbors excludes the cell itself."""
        target_area = 10000.0
        geom = Polygon([[0, 0], [100, 0], [100, 100], [0, 100], [0, 0]])
        quality = is_good_cell(geom, target_area)
        cell = Cell(id=0, geom=geom, quality=quality)

        neighbors = get_neighbors(cell, [cell])
        assert len(neighbors) == 0

    def test_get_neighbors_empty_list(self):
        """Test get_neighbors with empty cell list."""
        target_area = 10000.0
        geom = Polygon([[0, 0], [100, 0], [100, 100], [0, 100], [0, 0]])
        quality = is_good_cell(geom, target_area)
        cell = Cell(id=0, geom=geom, quality=quality)

        neighbors = get_neighbors(cell, [])
        assert len(neighbors) == 0


class TestFindNeighborCornerSplitLines:
    """Tests for find_neighbor_corner_split_lines function."""

    def test_find_split_lines_with_good_neighbors(self):
        """Test finding split lines with good neighboring cells."""
        target_area = 10000.0

        # Create grid: 3 good cells + 1 oversized center cell
        cells = []

        # Top good cell
        geom1 = Polygon([[100, 200], [200, 200], [200, 300], [100, 300], [100, 200]])
        quality1 = is_good_cell(geom1, target_area)
        cells.append(Cell(id=0, geom=geom1, quality=quality1))

        # Left good cell
        geom2 = Polygon([[0, 100], [100, 100], [100, 200], [0, 200], [0, 100]])
        quality2 = is_good_cell(geom2, target_area)
        cells.append(Cell(id=1, geom=geom2, quality=quality2))

        # Right good cell
        geom3 = Polygon([[200, 100], [300, 100], [300, 200], [200, 200], [200, 100]])
        quality3 = is_good_cell(geom3, target_area)
        cells.append(Cell(id=2, geom=geom3, quality=quality3))

        # Center oversized cell (200m × 100m)
        bad_geom = Polygon([[100, 100], [300, 100], [300, 200], [100, 200], [100, 100]])
        bad_quality = is_good_cell(bad_geom, target_area)
        bad_cell = Cell(id=3, geom=bad_geom, quality=bad_quality)
        cells.append(bad_cell)

        # Find split lines
        split_lines = find_neighbor_corner_split_lines(bad_cell, cells)

        # Should find multiple potential split lines
        assert len(split_lines) > 0

    def test_find_split_lines_no_neighbors(self):
        """Test finding split lines with no neighbors."""
        target_area = 10000.0
        geom = Polygon([[0, 0], [200, 0], [200, 100], [0, 100], [0, 0]])
        quality = is_good_cell(geom, target_area)
        cell = Cell(id=0, geom=geom, quality=quality)

        split_lines = find_neighbor_corner_split_lines(cell, [cell])
        assert len(split_lines) == 0

    def test_find_split_lines_insufficient_corners(self):
        """Test finding split lines with only 1 neighbor corner."""
        target_area = 10000.0

        # Bad cell
        bad_geom = Polygon([[0, 0], [200, 0], [200, 100], [0, 100], [0, 0]])
        bad_quality = is_good_cell(bad_geom, target_area)
        bad_cell = Cell(id=0, geom=bad_geom, quality=bad_quality)

        # Neighbor that only touches at one point
        neighbor_geom = Polygon([[200, 0], [300, 0], [300, 100], [200, 100], [200, 0]])
        neighbor_quality = is_good_cell(neighbor_geom, target_area)
        neighbor = Cell(id=1, geom=neighbor_geom, quality=neighbor_quality)

        cells = [bad_cell, neighbor]
        split_lines = find_neighbor_corner_split_lines(bad_cell, cells)

        # Should still find lines even with few neighbors
        assert isinstance(split_lines, list)


class TestFixAreaTooLarge:
    """Tests for fix_area_too_large function."""

    def test_fix_oversized_cell_2x(self):
        """Test fixing a cell that is 2x target area."""
        target_area = 10000.0

        # Create 3x3 grid with center cell oversized
        cells = []
        for row in range(3):
            for col in range(3):
                if row == 1 and col == 1:
                    # Center cell: 200m × 100m (oversized)
                    geom = Polygon([[100, 100], [300, 100], [300, 200], [100, 200], [100, 100]])
                else:
                    # Regular 100m × 100m cells
                    x = col * 100
                    y = row * 100
                    geom = Polygon([[x, y], [x + 100, y], [x + 100, y + 100], [x, y + 100], [x, y]])
                quality = is_good_cell(geom, target_area)
                cell = Cell(id=len(cells), geom=geom, quality=quality)
                cells.append(cell)

        # Fix the oversized center cell
        bad_cell = cells[4]
        assert bad_cell.quality["reason"] == "area_too_large"

        fixed = fix_area_too_large(bad_cell, target_area, all_cells=cells)

        # Should split into multiple cells
        assert len(fixed) > 1

    def test_fix_returns_unchanged_if_not_oversized(self):
        """Test that fix returns unchanged cell if not oversized."""
        target_area = 10000.0
        geom = Polygon([[0, 0], [100, 0], [100, 100], [0, 100], [0, 0]])
        quality = is_good_cell(geom, target_area)
        cell = Cell(id=0, geom=geom, quality=quality)

        fixed = fix_area_too_large(cell, target_area)

        assert len(fixed) == 1
        assert fixed[0] == cell

    def test_fix_without_neighbors(self):
        """Test that fix returns unchanged cell without neighbor info."""
        target_area = 10000.0
        geom = Polygon([[0, 0], [200, 0], [200, 100], [0, 100], [0, 0]])
        quality = is_good_cell(geom, target_area)
        cell = Cell(id=0, geom=geom, quality=quality)

        fixed = fix_area_too_large(cell, target_area, all_cells=None)

        # Without neighbors, should return unchanged
        assert len(fixed) == 1
        assert fixed[0] == cell


class TestInspectAndFixCells:
    """Tests for inspect_and_fix_cells function."""

    def test_inspect_without_fixing(self):
        """Test inspection without fixing enabled."""
        target_area = 10000.0

        # Create 3 cells: 2 good, 1 bad
        cells = []

        # Good cell
        geom1 = Polygon([[0, 0], [100, 0], [100, 100], [0, 100], [0, 0]])
        quality1 = is_good_cell(geom1, target_area)
        cells.append(Cell(id=0, geom=geom1, quality=quality1))

        # Good cell
        geom2 = Polygon([[100, 0], [200, 0], [200, 100], [100, 100], [100, 0]])
        quality2 = is_good_cell(geom2, target_area)
        cells.append(Cell(id=1, geom=geom2, quality=quality2))

        # Oversized cell
        geom3 = Polygon([[0, 100], [200, 100], [200, 200], [0, 200], [0, 100]])
        quality3 = is_good_cell(geom3, target_area)
        cells.append(Cell(id=2, geom=geom3, quality=quality3))

        # Inspect without fixing
        result = inspect_and_fix_cells(cells, target_area, fix_oversized=False)

        # Should return same cells unchanged
        assert len(result) == len(cells)

    def test_inspect_with_fixing(self):
        """Test inspection with fixing enabled."""
        target_area = 10000.0

        # Create 3x3 grid with one oversized center cell
        cells = []
        for row in range(3):
            for col in range(3):
                if row == 1 and col == 1:
                    # Center cell: 200m × 100m (oversized)
                    geom = Polygon([[100, 100], [300, 100], [300, 200], [100, 200], [100, 100]])
                else:
                    # Regular 100m × 100m cells
                    x = col * 100
                    y = row * 100
                    geom = Polygon([[x, y], [x + 100, y], [x + 100, y + 100], [x, y + 100], [x, y]])
                quality = is_good_cell(geom, target_area)
                cell = Cell(id=len(cells), geom=geom, quality=quality)
                cells.append(cell)

        # Count bad cells before
        bad_before = sum(1 for c in cells if not c.quality["is_good"])

        # Inspect with fixing
        result = inspect_and_fix_cells(cells, target_area, fix_oversized=True)

        # Should have more cells (from splitting)
        assert len(result) >= len(cells)

        # Should have fewer or same bad cells
        bad_after = sum(1 for c in result if not c.quality["is_good"])
        assert bad_after <= bad_before

    def test_inspect_all_good_cells(self):
        """Test inspection when all cells are good."""
        target_area = 10000.0

        # Create 4 good cells
        cells = []
        for row in range(2):
            for col in range(2):
                x = col * 100
                y = row * 100
                geom = Polygon([[x, y], [x + 100, y], [x + 100, y + 100], [x, y + 100], [x, y]])
                quality = is_good_cell(geom, target_area)
                cell = Cell(id=len(cells), geom=geom, quality=quality)
                cells.append(cell)

        # Inspect with fixing
        result = inspect_and_fix_cells(cells, target_area, fix_oversized=True)

        # Should return same cells unchanged
        assert len(result) == len(cells)

    def test_sequential_fixing(self):
        """Test that sequential fixing allows later fixes to benefit from earlier ones."""
        target_area = 10000.0

        # Create grid with 2 oversized cells
        cells = []

        # Row 1 - good cells
        for col in range(3):
            x = col * 100
            geom = Polygon([[x, 200], [x + 100, 200], [x + 100, 300], [x, 300], [x, 200]])
            quality = is_good_cell(geom, target_area)
            cells.append(Cell(id=len(cells), geom=geom, quality=quality))

        # Row 2 - left good, center and right oversized
        geom_left = Polygon([[0, 100], [100, 100], [100, 200], [0, 200], [0, 100]])
        quality_left = is_good_cell(geom_left, target_area)
        cells.append(Cell(id=len(cells), geom=geom_left, quality=quality_left))

        geom_center = Polygon([[100, 100], [200, 100], [200, 200], [100, 200], [100, 100]])
        quality_center = is_good_cell(geom_center, target_area)
        # Make it oversized
        geom_center = Polygon([[100, 100], [300, 100], [300, 200], [100, 200], [100, 100]])
        quality_center = is_good_cell(geom_center, target_area)
        cells.append(Cell(id=len(cells), geom=geom_center, quality=quality_center))

        geom_right = Polygon([[200, 100], [400, 100], [400, 200], [200, 200], [200, 100]])
        quality_right = is_good_cell(geom_right, target_area)
        cells.append(Cell(id=len(cells), geom=geom_right, quality=quality_right))

        # Row 3 - good cells
        for col in range(3):
            x = col * 100
            geom = Polygon([[x, 0], [x + 100, 0], [x + 100, 100], [x, 100], [x, 0]])
            quality = is_good_cell(geom, target_area)
            cells.append(Cell(id=len(cells), geom=geom, quality=quality))

        # Fix with sequential fixing
        result = inspect_and_fix_cells(cells, target_area, fix_oversized=True)

        # Should have more cells than original
        assert len(result) > len(cells)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
