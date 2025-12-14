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

### Basic Examples

#### basic_usage.py

Demonstrates core functionality:
- Simple greeting function
- Rich-formatted library information display

#### geo_usage.py

Demonstrates geospatial functionality:
- GDAL version information
- OGR driver enumeration
- Creating and working with GeoDataFrames

**Note:** The geospatial example requires GDAL and GeoPandas to be installed. Use the Nix development environment for the easiest setup.

### Workflow Examples

The workflow examples demonstrate the complete process of generating street blocks and clusters from site boundaries and roads.

#### step1_generate_parcels.py

Generates parcels from a site boundary.

**Input:**
- Site boundary GeoJSON

**Output:**
- `outputs/step1_parcels/parcels.geojson`

```bash
python examples/step1_generate_parcels.py
```

#### step2_generate_streets.py

Generates street blocks from parcels and roads.

**Input:**
- Parcels from step 1
- Roads GeoJSON

**Output:**
- `outputs/step2_streets/outputs.gpkg` - Full GeoPackage with all intermediate layers
- `outputs/step2_streets/all_grids_merged.geojson` - Merged grid cells with grid_type classification

The output includes:
- On-grid blocks (along arterial, secondary, and intersected setback areas)
- Off-grid blocks (away from road setbacks)

```bash
python examples/step2_generate_streets.py
```

#### step3_generate_clusters.py

Generates clusters/partitions from street blocks.

**Input:**
- Site boundary (from step 1)
- Roads GeoJSON
- Street blocks (from step 2)

**Output:**
- `outputs/step3_clusters/clusters.gpkg` - GeoPackage with three layers:
  - `off_grid_corner_parts` - Corner fan parts at sharp angles
  - `off_grid_side_parts` - Side strip parts along edges
  - `off_grid_center_parts` - Center parts of off-grid blocks

```bash
python examples/step3_generate_clusters.py
```

#### step4_generate_public.py

Generates public spaces from clusters.

**Input:**
- Clusters (from step 3)

**Output:**
- `outputs/step4_public/outputs.gpkg` - GeoPackage with public spaces

```bash
python examples/step4_generate_public.py
```

### Running the Full Workflow

To run the complete workflow from start to finish:

```bash
# Step 1: Generate parcels
python examples/step1_generate_parcels.py

# Step 2: Generate street blocks
python examples/step2_generate_streets.py

# Step 3: Generate clusters
python examples/step3_generate_clusters.py

# Step 4: Generate public spaces
python examples/step4_generate_public.py
```
