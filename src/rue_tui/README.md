# RUE TUI - Terminal User Interface for RUE-lib

A beautiful, interactive terminal application for the Rapid Urbanisation Explorer (RUE) library. Provides a guided workflow for urban planning analysis with real-time visualization and progress tracking.

![RUE TUI Demo](demo.gif)

## Features

🏗️ **Guided Workflow**: Step-by-step process through parcel generation, street planning, and cluster analysis

📊 **Real-time Progress**: Live progress tracking with detailed status updates

🖼️ **Map Visualization**: Automatic rendering of analysis results to images with terminal display via `fim`

📁 **Results Browser**: Interactive browser for exploring geopackage layers and outputs

⚙️ **Configuration**: Easy configuration of input files, output settings, and visualization options

📝 **Comprehensive Logging**: Detailed logging with different levels and timestamps

## Installation

### Prerequisites

RUE TUI requires the RUE-lib library and some additional dependencies:

```bash
# Install RUE-lib first
pip install rue-lib

# Install RUE TUI
pip install rue-tui
```

### Optional: Image Viewing

For terminal-based image viewing, install one of these image viewers:

```bash
# Frame buffer image viewer (recommended for terminals)
sudo apt-get install fim

# Or X11 image viewer
sudo apt-get install feh

# Or ImageMagick
sudo apt-get install imagemagick
```

## Usage

### Basic Usage

Launch the TUI application:

```bash
rue-tui
```

Or run as a Python module:

```bash
python -m rue_tui
```

### Workflow

1. **Configure**: Set up input files and output settings in the Config tab
2. **Step 1 - Generate Parcels**: Create ownership parcels from site boundary and roads
3. **Step 2 - Generate Streets**: Design street networks and block structures
4. **Step 3 - Generate Clusters**: Create final urban clusters and layout
5. **Browse Results**: Explore generated layers and visualizations

### Input Files

The TUI expects two input files:

- **Site GeoJSON** (`data/site.geojson`): Polygon defining the development site boundary
- **Roads GeoJSON** (`data/roads.geojson`): LineString features representing the road network

### Output

Results are saved to the configured output directory (default: `outputs/`) including:

- **GeoPackage** (`output.gpkg`): All analysis layers in a single file
- **Step Results**: Individual GeoJSON files for each step
- **Visualizations**: PNG images of map renderings

## Interface

### Tabs

- **Workflow**: Main interface with step cards, progress tracking, and action buttons
- **Config**: Configuration panel for input files, output settings, and visualization options
- **Results**: Layer browser, log viewer, and result visualization controls

### Keyboard Shortcuts

- `Tab`: Navigate between interface elements
- `Enter`: Activate buttons and controls
- `Escape`: Close modal dialogs
- `Ctrl+C`: Exit application

## Configuration Options

### File Paths
- `site_path`: Path to site boundary GeoJSON
- `roads_path`: Path to roads network GeoJSON
- `output_dir`: Directory for output files

### Visualization
- `image_width/height`: Size of rendered map images
- `dpi`: Resolution for map rendering
- `color_scheme`: Color scheme for visualizations

### Workflow
- `auto_advance`: Automatically proceed to next step on completion
- `save_intermediate`: Save intermediate processing results

## API Integration

RUE TUI is built as a separate module that imports and uses RUE-lib:

```python
from rue_tui import RueTuiApp, TuiConfig

# Custom configuration
config = TuiConfig(
    site_path="my_site.geojson",
    roads_path="my_roads.geojson",
    output_dir="my_outputs",
    auto_advance=True
)

# Run with custom config
app = RueTuiApp()
app.config = config
app.run()
```

## Architecture

RUE TUI is built using:

- **[Textual](https://textual.textualize.io/)**: Modern TUI framework with rich widgets
- **[Rich](https://rich.readthedocs.io/)**: Beautiful terminal formatting and display
- **[Matplotlib](https://matplotlib.org/)**: Map rendering and visualization
- **[GeoPandas](https://geopandas.org/)**: Geospatial data processing

### Components

- `app.py`: Main application and workflow orchestration
- `widgets.py`: Custom UI components (step cards, progress displays, etc.)
- `visualizer.py`: Map rendering and image generation
- `config.py`: Configuration management

## Development

### Setup Development Environment

```bash
git clone https://github.com/timlinux/rue-lib.git
cd rue-lib/rue_tui

# Install in development mode
pip install -e .[dev]
```

### Running Tests

```bash
pytest tests/
```

### Code Formatting

```bash
black rue_tui/
isort rue_tui/
```

### Type Checking

```bash
mypy rue_tui/
```

## Troubleshooting

### Image Viewing Issues

If images don't display in the terminal:

1. Install an image viewer: `sudo apt-get install fim`
2. Check if running in a compatible terminal (framebuffer support)
3. Images are saved to temp files that can be opened manually

### Performance

For large datasets:

- Reduce image resolution in config
- Disable auto-advance for manual control
- Monitor memory usage during processing

### Dependencies

Ensure all required packages are installed:

```bash
pip install textual rich geopandas matplotlib fiona
```

## Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## License

MIT License - see LICENSE file for details.

## Related Projects

- [RUE-lib](https://github.com/timlinux/rue-lib): Core urban planning analysis library
- [Textual](https://textual.textualize.io/): TUI framework
- [Rich](https://rich.readthedocs.io/): Terminal formatting library
