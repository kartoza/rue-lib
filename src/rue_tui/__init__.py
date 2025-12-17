"""
RUE TUI - Terminal User Interface for RUE-lib
===========================================

A beautiful terminal interface for the Rapid Urbanisation Explorer library.
Guides users through parcel generation, street planning, and cluster analysis
with real-time visualization of results.

Features:
- Step-by-step guided workflow
- Real-time progress tracking
- Interactive result browser
- Map visualization with fim
- Rich terminal graphics
"""

__version__ = "0.1.0"
__author__ = "RUE-lib Team"

from .app import RueTuiApp
from .config import TuiConfig

__all__ = ["RueTuiApp", "TuiConfig"]
