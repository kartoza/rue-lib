"""Tests for core functionality."""

from rue_lib.core import display_info, greet


def test_greet_default():
    """Test greet with default name."""
    result = greet()
    assert result == "Hello, World!"


def test_greet_with_name():
    """Test greet with custom name."""
    result = greet("Alice")
    assert result == "Hello, Alice!"


def test_greet_with_none():
    """Test greet with None explicitly."""
    result = greet(None)
    assert result == "Hello, World!"


def test_display_info():
    """Test display_info runs without error."""
    # This should not raise an exception
    display_info()
