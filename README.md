# rue-lib

Python library for the Rapid Urbanisation Explorer project

## Quick Start

**New to the project?** See [QUICKSTART.md](QUICKSTART.md) for the fastest way to get started!

## Features

- Geospatial data processing with GDAL/OGR and GeoPandas
- Rich terminal output formatting
- Ready for PyPI publishing
- Nix Flake for reproducible development environment
- Comprehensive pre-commit QA checks

## Installation

### System Requirements

This library requires GDAL to be installed on your system. Install it before installing rue-lib:

**macOS (using Homebrew):**
```bash
brew install gdal
```

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install -y gdal-bin libgdal-dev
```

**Windows:**
Download from [OSGeo4W](https://trac.osgeo.org/osgeo4w/) or use conda.

### From PyPI (when published)

```bash
# Install GDAL Python bindings matching your system GDAL version
pip install GDAL=="$(gdal-config --version)"

# Then install rue-lib
pip install rue-lib
```

### From Source

```bash
git clone https://github.com/kartoza/rue-lib.git
cd rue-lib

# Install GDAL Python bindings first
pip install GDAL=="$(gdal-config --version)"

# Then install rue-lib
pip install -e .
```

### Using Nix

This project includes a Nix flake for reproducible development:

```bash
# Enter development environment
nix develop

# Run tests
nix run .#test

# Run tests with coverage
nix run .#test-cov
```

## Development Environment

The Nix flake provides:
- Python 3.11 with all dependencies (GDAL, GeoPandas, Rich)
- Neovim with timvim configuration
- Pre-commit hooks (Ruff, Bandit, CSpell)
- All necessary build tools

### Setting up Pre-commit Hooks

```bash
# In the nix develop shell
pre-commit install

# Run checks manually
pre-commit run --all-files
```

## Usage

```python
from rue_lib.core import greet, display_info
from rue_lib.geo import create_sample_geodataframe

# Basic greeting
print(greet("World"))

# Display library info with rich formatting
display_info()

# Work with geospatial data
gdf = create_sample_geodataframe()
print(gdf.head())
```

## Testing

```bash
# Using pytest directly
pytest tests/

# With coverage
pytest tests/ --cov=rue_lib --cov-report=html

# Using Nix
nix run .#test

# Using Make
make test
```

## Dependencies

Core dependencies:
- Rich >= 13.0.0
- GeoPandas >= 0.14.0
- Shapely >= 2.0.0
- PyProj >= 3.6.0
- Pandas >= 2.0.0
- GDAL >= 3.4.0 (must be installed separately to match system library)

Development dependencies:
- pytest >= 7.0.0
- pytest-cov >= 4.0.0
- ruff >= 0.1.0
- bandit >= 1.7.0
- pre-commit >= 3.0.0

**Note:** GDAL must be installed separately and match your system's GDAL library version. See the Installation section above for details.

## Documentation

- [CONTRIBUTING.md](CONTRIBUTING.md) - Development setup and contribution guidelines
- [NIX_SETUP.md](NIX_SETUP.md) - Detailed Nix flake setup guide
- [PUBLISHING.md](PUBLISHING.md) - How to publish to PyPI
- [examples/](examples/) - Usage examples

## License

MIT License - see LICENSE file for details

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for detailed instructions on:
- Setting up the development environment
- Running tests and checks
- Code style guidelines
- Submitting pull requests
