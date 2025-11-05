# rue-lib Examples

This directory contains example scripts demonstrating the usage of rue-lib.

## Running Examples

### Using Nix (Recommended)

```bash
# Enter the development environment
nix develop

# Run the basic example
python examples/basic_usage.py

# Run the geospatial example
python examples/geo_usage.py
```

### Using Virtual Environment

```bash
# Activate your virtual environment
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Run the basic example
python examples/basic_usage.py

# Run the geospatial example (requires GDAL and GeoPandas)
python examples/geo_usage.py
```

## Examples

### basic_usage.py

Demonstrates core functionality:
- Simple greeting function
- Rich-formatted library information display

### geo_usage.py

Demonstrates geospatial functionality:
- GDAL version information
- OGR driver enumeration
- Creating and working with GeoDataFrames

**Note:** The geospatial example requires GDAL and GeoPandas to be installed. Use the Nix development environment for the easiest setup.
