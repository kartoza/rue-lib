#!/usr/bin/env python3
"""
Launcher script for RUE TUI application.

This script can be used to launch the TUI application with various options
and ensures the proper Python environment is configured.
"""

import sys
from pathlib import Path

# Add the rue_tui directory to Python path
script_dir = Path(__file__).parent
repo_root = script_dir.parent
tui_dir = repo_root / "rue_tui"
src_dir = repo_root / "src"

sys.path.insert(0, str(tui_dir))
sys.path.insert(0, str(src_dir))


def main():
    """Launch the RUE TUI application."""
    try:
        from rue_tui.app import main as run_tui

        print("🏙️ Launching RUE TUI...")
        print(f"📁 Repository root: {repo_root}")
        print(f"🐍 Python: {sys.executable}")
        print("=" * 50)

        run_tui()

    except ImportError as e:
        print(f"❌ Import error: {e}")
        print("\n💡 Make sure you have installed the required dependencies:")
        print("   pip install textual rich geopandas matplotlib")
        sys.exit(1)

    except KeyboardInterrupt:
        print("\n👋 Goodbye!")
        sys.exit(0)

    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
