# rue-lib

Python library for the Rapid Urbanisation Explorer project

## Features

- Geospatial data processing with GDAL/OGR and GeoPandas
- Rich terminal output formatting
- Ready for PyPI publishing
- Nix Flake for reproducible development environment
- Comprehensive pre-commit QA checks

## Installation

### From PyPI (when published)

```bash
pip install rue-lib
```

### From Source

```bash
git clone https://github.com/kartoza/rue-lib.git
cd rue-lib
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
```

## Dependencies

Core dependencies:
- GDAL >= 3.0.0
- GeoPandas >= 0.10.0
- Rich >= 13.0.0

## License

MIT License - see LICENSE file for details

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run pre-commit checks: `pre-commit run --all-files`
5. Submit a pull request
