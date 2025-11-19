# Contributing to rue-lib

Thank you for your interest in contributing to rue-lib! This guide will help you get started with the development environment.

## Development Environment Setup

### Option 1: Using Nix (Recommended)

This project uses Nix flakes for a reproducible development environment. If you have Nix installed with flakes enabled, you can get started immediately:

```bash
# Enter the development environment
nix develop

# The environment includes:
# - Python 3.11 with all dependencies
# - Neovim with timvim configuration
# - GDAL, GeoPandas, Rich
# - Development tools (pytest, ruff, bandit, pre-commit)
# - All required system libraries
```

#### Running Tests with Nix

```bash
# Run all tests
nix run .#test

# Run tests with coverage
nix run .#test-cov
```

### Option 2: Using Python Virtual Environment

If you don't have Nix installed:

```bash
# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -e ".[dev]"

# Install GDAL (system-dependent, you may need to install GDAL libraries first)
# On Ubuntu/Debian:
# sudo apt-get install gdal-bin libgdal-dev
# On macOS:
# brew install gdal

pip install gdal geopandas
```

## Pre-commit Hooks

We use pre-commit hooks to ensure code quality. These are automatically installed when you enter the Nix development environment.

### Manual Installation

```bash
# Install pre-commit hooks
pre-commit install

# Run all checks manually
pre-commit run --all-files
```

### Included Checks

- **Ruff**: Fast Python linter and formatter
- **Bandit**: Security vulnerability scanner
- **CSpell**: Spell checker for code and documentation
- **Standard hooks**: Trailing whitespace, end-of-file fixer, YAML/JSON/TOML validation

## Code Style

We use Ruff for both linting and formatting. The configuration is in `pyproject.toml`.

```bash
# Check code style
ruff check src/ tests/

# Auto-fix issues
ruff check --fix src/ tests/

# Format code
ruff format src/ tests/
```

## Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_core.py

# Run with coverage
pytest --cov=rue_lib --cov-report=html

# Run with verbose output
pytest -v
```

## Security Checks

```bash
# Run Bandit security scanner
bandit -c pyproject.toml -r src/
```

## Project Structure

```
rue-lib/
├── src/
│   └── rue_lib/          # Main package
│       ├── __init__.py   # Package initialization
│       ├── core.py       # Core functionality
│       └── geo.py        # Geospatial functionality
├── tests/                # Test files
│   ├── test_core.py
│   └── test_geo.py
├── pyproject.toml        # Project configuration and dependencies
├── flake.nix             # Nix flake for reproducible environment
├── .pre-commit-config.yaml  # Pre-commit hook configuration
├── .cspell.json          # Spell checker configuration
└── README.md             # Project documentation
```

## Making Changes

1. **Create a branch** for your feature or bugfix
2. **Make your changes** following the code style guidelines
3. **Write tests** for new functionality
4. **Run the test suite** to ensure everything works
5. **Run pre-commit checks** to validate your changes
6. **Commit your changes** with a clear message
7. **Push your branch** and create a pull request

## Commit Messages

Follow conventional commit format:

- `feat: add new feature`
- `fix: resolve bug in X`
- `docs: update documentation`
- `test: add tests for Y`
- `refactor: restructure code`
- `chore: update dependencies`

## Neovim with timvim

The Nix development environment includes Neovim configured with timvim. To use it:

```bash
# In the nix develop shell
nvim
```

The timvim configuration is automatically loaded, providing a powerful development environment for Python.

## Questions?

If you have questions or run into issues, please:
- Check the README.md for basic usage
- Open an issue on GitHub
- Review existing issues and discussions

Thank you for contributing!
