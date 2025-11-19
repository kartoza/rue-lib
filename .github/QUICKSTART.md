# Quick Start Guide

Get started with rue-lib in minutes!

## Option 1: Quick Start with Nix (Recommended)

If you have Nix installed with flakes enabled:

```bash
# Clone and enter the repository
git clone https://github.com/kartoza/rue-lib.git
cd rue-lib

# Enter development environment (downloads everything automatically)
nix develop

# Run tests
nix run .#test

# Try the examples
python examples/basic_usage.py
```

That's it! You now have a complete development environment with Python, GDAL, GeoPandas, Neovim with timvim, and all development tools.

## Option 2: Quick Start with pip

```bash
# Clone the repository
git clone https://github.com/kartoza/rue-lib.git
cd rue-lib

# Create and activate virtual environment
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install the package
pip install -e .

# Try the basic example (works without geo dependencies)
python examples/basic_usage.py
```

### For Geospatial Features

To use geospatial features, you need GDAL and GeoPandas:

```bash
# Install geo dependencies (requires system GDAL libraries)
pip install -e ".[geo]"

# Try the geo example
python examples/geo_usage.py
```

**Note:** Installing GDAL can be tricky. If you encounter issues, use the Nix approach or see the [CONTRIBUTING.md](CONTRIBUTING.md) for detailed instructions.

## Using the Library

### Basic Usage

```python
from rue_lib.core import greet, display_info

# Simple greeting
print(greet("World"))

# Display library info with rich formatting
display_info()
```

### Geospatial Usage (requires geo dependencies)

```python
from rue_lib.geo import create_sample_geodataframe

# Create a sample GeoDataFrame
gdf = create_sample_geodataframe()
print(gdf)
```

## Running Tests

```bash
# With Nix
nix run .#test

# With pytest directly
pytest tests/

# With Make
make test
```

## Development Tasks

```bash
# Install pre-commit hooks
pre-commit install

# Run linter
ruff check src/ tests/

# Run security checks
bandit -c pyproject.toml -r src/

# Format code
ruff format src/ tests/

# Or use Make for any of these:
make lint
make security
make format
make precommit
```

## Next Steps

- Read [README.md](README.md) for comprehensive documentation
- Check [CONTRIBUTING.md](CONTRIBUTING.md) for development guidelines
- See [NIX_SETUP.md](NIX_SETUP.md) for detailed Nix information
- Review [PUBLISHING.md](PUBLISHING.md) if you want to publish to PyPI
- Explore [examples/](examples/) for more usage examples

## Getting Help

- Check the documentation files listed above
- Open an issue on GitHub if you encounter problems
- All tests should pass; if they don't, something is wrong with the setup

## What You Get with Nix

The Nix development environment provides:
- ✅ Python 3.11 with all dependencies
- ✅ GDAL, GeoPandas, Rich pre-installed
- ✅ Neovim with timvim configuration
- ✅ All development tools (pytest, ruff, bandit, cspell)
- ✅ Pre-commit hooks automatically installed
- ✅ No system pollution - everything is isolated
- ✅ Reproducible across all machines

## Common Commands

```bash
# Nix
nix develop              # Enter dev environment
nix run .#test          # Run tests
nix run .#test-cov      # Run tests with coverage

# Make
make help               # Show all available commands
make test               # Run tests
make lint               # Run linter
make format             # Format code
make security           # Security checks

# Direct Python
pytest tests/           # Run tests
ruff check src/         # Lint
ruff format src/        # Format
bandit -c pyproject.toml -r src/  # Security
```

Happy coding! 🚀
