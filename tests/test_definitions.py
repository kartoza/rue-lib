"""Tests for core definitions module."""

from rue_lib.core.definitions import BlockTypes, ClusterTypes, ColorTypes, PropertyKeys, RoadTypes


class TestPropertyKeys:
    """Test PropertyKeys enumeration."""

    def test_property_keys_values(self):
        """Test PropertyKeys have expected values."""
        assert PropertyKeys.RoadType == "type"
        assert PropertyKeys.GridType == "type"
        assert PropertyKeys.BlockType == "type"

    def test_property_keys_exist(self):
        """Test all property keys are accessible."""
        assert hasattr(PropertyKeys, "RoadType")
        assert hasattr(PropertyKeys, "GridType")
        assert hasattr(PropertyKeys, "BlockType")


class TestRoadTypes:
    """Test RoadTypes enumeration."""

    def test_road_types_values(self):
        """Test RoadTypes have expected values."""
        assert RoadTypes.Artery == "road_arterial"
        assert RoadTypes.Secondary == "road_secondary"
        assert RoadTypes.Tertiary == "road_tertiary"
        assert RoadTypes.Local == "road_local"

    def test_road_types_exist(self):
        """Test all road types are accessible."""
        assert hasattr(RoadTypes, "Artery")
        assert hasattr(RoadTypes, "Secondary")
        assert hasattr(RoadTypes, "Tertiary")
        assert hasattr(RoadTypes, "Local")


class TestBlockTypes:
    """Test BlockTypes enumeration."""

    def test_block_types_values(self):
        """Test BlockTypes have expected values."""
        assert BlockTypes.ON_GRID_ART == "on_grid_art"
        assert BlockTypes.ON_GRID_SEC == "on_grid_sec"
        assert BlockTypes.OFF_GRID == "off_grid"
        assert BlockTypes.COLD_GRID == "cold_boundaries"

    def test_block_types_exist(self):
        """Test all block types are accessible."""
        assert hasattr(BlockTypes, "ON_GRID_ART")
        assert hasattr(BlockTypes, "ON_GRID_SEC")
        assert hasattr(BlockTypes, "OFF_GRID")
        assert hasattr(BlockTypes, "COLD_GRID")


class TestClusterTypes:
    """Test ClusterTypes enumeration."""

    def test_single_road_types(self):
        """Test single road type cluster definitions."""
        assert ClusterTypes.ON_GRID_ART == "art"
        assert ClusterTypes.ON_GRID_SEC == "sec"
        assert ClusterTypes.ON_GRID_LOC == "loc"

    def test_corner_types_arterial(self):
        """Test arterial corner type cluster definitions."""
        assert ClusterTypes.ON_GRID_ART_ART == "art_art"
        assert ClusterTypes.ON_GRID_ART_SEC == "art_sec"
        assert ClusterTypes.ON_GRID_ART_LOC == "art_loc"

    def test_corner_types_secondary(self):
        """Test secondary corner type cluster definitions."""
        assert ClusterTypes.ON_GRID_SEC_SEC == "sec_sec"
        assert ClusterTypes.ON_GRID_SEC_LOC == "sec_loc"

    def test_corner_types_local(self):
        """Test local corner type cluster definitions."""
        assert ClusterTypes.ON_GRID_LOC_LOC == "loc_loc"

    def test_off_grid_types(self):
        """Test off-grid cluster type definitions."""
        assert ClusterTypes.OFF_GRID == "off_grid"
        assert ClusterTypes.OFF_GRID_WARM == "off_grid0"
        assert ClusterTypes.OFF_GRID_COLD == "off_grid1"
        assert ClusterTypes.OFF_GRID_COLD2 == "off_grid2"
        assert ClusterTypes.CONCAVE_CORNER == "concave_corner"

    def test_amenity_types(self):
        """Test amenity cluster type definitions."""
        assert ClusterTypes.AMENITY == "am"
        assert ClusterTypes.ART_AMENITY == "art_am"
        assert ClusterTypes.SEC_AMENITY == "sec_am"
        assert ClusterTypes.SEC_SEC_AMENITY == "sec_sec_am"
        assert ClusterTypes.SEC_LOC_AMENITY == "sec_loc_am"
        assert ClusterTypes.LOC_AMENITY == "loc_am"
        assert ClusterTypes.LOC_LOC_AMENITY == "loc_loc_am"
        assert ClusterTypes.OFF_GRID0_AMENITY == "off_grid0_am"
        assert ClusterTypes.OFF_GRID1_AMENITY == "off_grid1_am"
        assert ClusterTypes.OFF_GRID2_AMENITY == "off_grid2_am"
        assert ClusterTypes.CONCAVE_CORNER_AMENITY == "concave_corner_am"

    def test_open_space_types(self):
        """Test open space cluster type definitions."""
        assert ClusterTypes.GREEN0 == "green0"
        assert ClusterTypes.GREEN1 == "green1"
        assert ClusterTypes.GREEN2 == "green2"
        assert ClusterTypes.OPEN_SPACE == "os"
        assert ClusterTypes.SEC_OPEN_SPACE == "sec_os"
        assert ClusterTypes.SEC_SEC_OPEN_SPACE == "sec_sec_os"
        assert ClusterTypes.SEC_LOC_OPEN_SPACE == "sec_loc_os"
        assert ClusterTypes.LOC_OPEN_SPACE == "loc_os"
        assert ClusterTypes.LOC_LOC_OPEN_SPACE == "loc_loc_os"
        assert ClusterTypes.OFF_GRID0_OPEN_SPACE == "off_grid0_os"
        assert ClusterTypes.OFF_GRID1_OPEN_SPACE == "off_grid1_os"
        assert ClusterTypes.OFF_GRID2_OPEN_SPACE == "off_grid2_os"
        assert ClusterTypes.CONCAVE_CORNER_OPEN_SPACE == "concave_corner_os"

    def test_path_and_entrance_types(self):
        """Test path and entrance cluster type definitions."""
        assert ClusterTypes.INTERNAL_PATH == "internal path, any cluster"
        assert ClusterTypes.PATH0 == "path0"
        assert ClusterTypes.PATH1 == "path1"
        assert ClusterTypes.PATH2 == "path2"
        assert ClusterTypes.ACCESS_PATH == "Cluster, access path"
        assert ClusterTypes.ENTRANCE0 == "entr0"
        assert ClusterTypes.ENTRANCE1 == "entr1"

    def test_other_types(self):
        """Test other cluster type definitions."""
        assert ClusterTypes.OFF_GRID_CLUSTER_2ND == "og cluster, 2nd"
        assert ClusterTypes.PUBLIC_STREETS == "Public streets"


class TestColorTypes:
    """Test ColorTypes mapping."""

    def test_color_types_is_dict(self):
        """Test ColorTypes is a dictionary."""
        assert isinstance(ColorTypes, dict)

    def test_single_road_colors(self):
        """Test single road type colors are defined."""
        assert ClusterTypes.ON_GRID_ART in ColorTypes
        assert ClusterTypes.ON_GRID_SEC in ColorTypes
        assert ClusterTypes.ON_GRID_LOC in ColorTypes

    def test_corner_type_colors(self):
        """Test corner type colors are defined."""
        assert ClusterTypes.ON_GRID_ART_ART in ColorTypes
        assert ClusterTypes.ON_GRID_ART_SEC in ColorTypes

    def test_color_format(self):
        """Test colors are in expected RGB format."""
        for color in ColorTypes.values():
            assert isinstance(color, str)
            assert color.startswith("rgb(")
            assert color.endswith(")")

    def test_color_values(self):
        """Test specific color values."""
        assert ColorTypes[ClusterTypes.ON_GRID_ART] == "rgb(235,201,199)"
        assert ColorTypes[ClusterTypes.ON_GRID_SEC] == "rgb(250,232,219)"
        assert ColorTypes[ClusterTypes.ON_GRID_LOC] == "rgb(250,242,212)"

    def test_all_cluster_types_have_colors(self):
        """Test that key cluster types have color mappings."""
        expected_types = [
            ClusterTypes.ON_GRID_ART,
            ClusterTypes.ON_GRID_SEC,
            ClusterTypes.ON_GRID_LOC,
            ClusterTypes.ON_GRID_ART_ART,
            ClusterTypes.ON_GRID_ART_SEC,
        ]
        for cluster_type in expected_types:
            assert cluster_type in ColorTypes, f"Missing color for {cluster_type}"
