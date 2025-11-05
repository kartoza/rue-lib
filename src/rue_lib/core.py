"""Core functionality for rue-lib."""

from typing import Optional

from rich.console import Console
from rich.table import Table

console = Console()


def greet(name: Optional[str] = None) -> str:
    """
    Return a greeting message.

    Args:
        name: Optional name to greet. Defaults to "World".

    Returns:
        A greeting string.
    """
    if name is None:
        name = "World"
    return f"Hello, {name}!"


def display_info() -> None:
    """Display information about the library using rich formatting."""
    table = Table(title="rue-lib Information")

    table.add_column("Property", style="cyan", no_wrap=True)
    table.add_column("Value", style="magenta")

    table.add_row("Name", "rue-lib")
    table.add_row("Purpose", "Rapid Urbanisation Explorer")
    table.add_row("Version", "0.1.0")

    console.print(table)
