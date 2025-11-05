# Nix Flake Setup Guide

This guide explains how to use the Nix flake for rue-lib development.

## What is Nix?

Nix is a powerful package manager that enables reproducible development environments. With Nix, everyone working on the project gets the exact same dependencies and tools, regardless of their operating system.

## Installing Nix

### Linux and macOS

```bash
# Install Nix with flakes enabled
sh <(curl -L https://nixos.org/nix/install) --daemon

# Enable flakes (if not already enabled)
mkdir -p ~/.config/nix
echo "experimental-features = nix-command flakes" >> ~/.config/nix/nix.conf
```

### Windows

Use WSL2 (Windows Subsystem for Linux) and follow the Linux instructions.

## Using the Development Environment

### Enter the development shell

```bash
# Clone the repository
git clone https://github.com/kartoza/rue-lib.git
cd rue-lib

# Enter the Nix development environment
nix develop
```

This gives you:
- Python 3.11 with all dependencies (GDAL, GeoPandas, Rich)
- Neovim with timvim configuration from github.com/timlinux/timvim
- Development tools (pytest, ruff, bandit, pre-commit, cspell)
- GDAL system libraries and tools
- Pre-commit hooks automatically installed

### What's included in the environment?

When you run `nix develop`, you get access to:

**Python Packages:**
- gdal
- geopandas
- rich
- shapely
- pytest
- pytest-cov
- ruff
- bandit
- pip, setuptools, wheel
- build, twine

**System Tools:**
- GDAL command-line tools
- GEOS libraries
- PROJ libraries
- git
- pre-commit
- cspell
- Neovim with timvim configuration

**Environment Variables:**
- `PYTHONPATH` configured to include `src/`
- `GDAL_DATA` pointing to GDAL data files
- `PROJ_LIB` pointing to PROJ library

## Running Tests

The flake provides convenient commands to run tests:

```bash
# Run all tests
nix run .#test

# Run tests with coverage report
nix run .#test-cov
```

These commands:
- Set up the correct Python environment
- Configure GDAL and PROJ paths
- Run pytest with the appropriate options

## Building the Package

```bash
# Build the Python package using Nix
nix build

# The result will be in ./result/
ls -l result/
```

## Neovim with timvim

The development environment includes Neovim configured with [timvim](https://github.com/timlinux/timvim), a comprehensive Neovim configuration.

```bash
# In the nix develop shell
nvim
```

Features:
- LSP (Language Server Protocol) support for Python
- Treesitter for syntax highlighting
- timvim plugins and configurations
- Nix-aware syntax highlighting

## Customizing the Flake

The flake configuration is in `flake.nix`. You can customize it for your needs:

### Adding Python packages

Edit the `pythonEnv` section:

```nix
pythonEnv = pkgs.python311.withPackages (ps: with ps; [
  # Add your packages here
  gdal
  geopandas
  rich
  your-package-here
]);
```

### Adding system tools

Edit the `buildInputs` in the `devShells.default` section:

```nix
buildInputs = [
  pythonEnv
  neovimWithTimvim
  pkgs.your-tool-here
];
```

## Continuous Integration

The project includes a GitHub Actions workflow that uses Nix for testing. See `.github/workflows/test.yml`.

## Flake Lock File

The `flake.lock` file pins all dependencies to specific versions, ensuring reproducibility. To update dependencies:

```bash
# Update all inputs
nix flake update

# Update a specific input
nix flake lock --update-input nixpkgs
```

## Advantages of Using Nix

1. **Reproducibility**: Everyone gets the same environment
2. **No conflicts**: Nix packages don't interfere with system packages
3. **Easy cleanup**: Just exit the shell; nothing is left behind
4. **Multiple versions**: Work on different projects with different dependency versions
5. **Cross-platform**: Works on Linux, macOS, and WSL2

## Troubleshooting

### Flake evaluation errors

```bash
# Check flake syntax
nix flake check

# Show detailed error messages
nix develop --show-trace
```

### Build failures

```bash
# Clean build cache
nix-collect-garbage

# Rebuild
nix develop --refresh
```

### Slow first build

The first time you run `nix develop`, Nix downloads and builds all dependencies. This can take several minutes but is cached for future use.

## Alternative: direnv

You can use [direnv](https://direnv.net/) to automatically enter the Nix environment when you `cd` into the project:

```bash
# Install direnv
nix-env -i direnv

# Create .envrc file
echo "use flake" > .envrc

# Allow direnv for this directory
direnv allow

# Now the environment activates automatically when you cd into the directory
```

## Resources

- [Nix Manual](https://nixos.org/manual/nix/stable/)
- [Nix Flakes Documentation](https://nixos.wiki/wiki/Flakes)
- [nix.dev](https://nix.dev/) - Nix tutorials
- [timvim Repository](https://github.com/timlinux/timvim)

## Getting Help

If you encounter issues with Nix:
1. Check the Nix manual and documentation
2. Ask in the project's GitHub issues
3. Visit the [NixOS Discourse](https://discourse.nixos.org/)
