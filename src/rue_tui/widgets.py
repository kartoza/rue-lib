"""Custom widgets for the RUE TUI application."""

import subprocess
from pathlib import Path
from typing import Optional

from rich import box
from rich.align import Align
from rich.console import Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text
from textual.app import ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.message import Message
from textual.widgets import (
    Button,
    Checkbox,
    DirectoryTree,
    Input,
    Label,
    Static,
    Tree,
)


class StepCard(Static):
    """A card widget representing a workflow step with beautiful Rich styling."""

    def __init__(
        self,
        step_num: int,
        title: str,
        description: str,
        status: str = "pending",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.step_num = step_num
        self.title = title
        self.description = description
        self._status = status

    def update_status(self, status: str):
        """Update the step status."""
        self._status = status
        self.update_display()

    def update_display(self):
        """Update the visual display of the card with Rich styling."""
        # Kartoza color scheme configuration
        status_config = {
            "pending": {
                "icon": "⏳",
                "emoji": ":hourglass_flowing_sand:",
                "color": "#569FC6",
                "bg_color": "#E8F4FD",
                "border": "blue",
                "style": "dim",
            },
            "running": {
                "icon": "⚡",
                "emoji": ":zap:",
                "color": "#DF9E2F",
                "bg_color": "#FFF8DC",
                "border": "yellow",
                "style": "bold blink",
            },
            "completed": {
                "icon": "✅",
                "emoji": ":white_check_mark:",
                "color": "#06969A",
                "bg_color": "#F6FFED",
                "border": "green",
                "style": "bold",
            },
            "error": {
                "icon": "❌",
                "emoji": ":cross_mark:",
                "color": "#CC0403",
                "bg_color": "#FFF2F0",
                "border": "red",
                "style": "bold",
            },
        }

        config = status_config.get(self._status, status_config["pending"])

        # Simple compact display
        compact_content = Text()
        compact_content.append(f"{config['icon']} ", style=f"{config['color']} bold")
        compact_content.append(
            f"Step {self.step_num}: {self.title} ", style=f"bold {config['color']}"
        )
        compact_content.append(f"({self._status})", style=f"dim {config['color']}")

        self.update(compact_content)


class ProgressDisplay(Static):
    """Beautiful progress display using Rich components."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.progress = None
        self.task_id = None
        self.spinner = None

    def start_progress(self, description: str, total: Optional[int] = None):
        """Start compact progress display."""
        # Compact progress display
        progress_text = Text()
        progress_text.append("🔄 ", style="bold #569FC6")
        progress_text.append(description, style="bold white")
        progress_text.append(" - Working...", style="dim #8A8B8B")

        self.update(progress_text)

    def complete_progress(self, message: str = "Completed"):
        """Show compact completion message."""
        # Compact success display
        success_text = Text()
        success_text.append("✅ ", style="bold #06969A")
        success_text.append(message, style="bold white")
        success_text.append(" - Done!", style="dim #06969A")

        self.update(success_text)

    def error_progress(self, message: str = "Error occurred"):
        """Show compact error message."""
        # Compact error display
        error_text = Text()
        error_text.append("❌ ", style="bold #CC0403")
        error_text.append(message, style="bold white")
        error_text.append(" - Failed!", style="dim #CC0403")

        self.update(error_text)


class LayerBrowser(Static):
    """Beautiful widget for browsing geopackage layers with Rich styling."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._layers = []

    def update_layers(self, layer_info: list[tuple]):
        """Update the layer information."""
        self._layers = layer_info
        self._render_beautiful_table()

    def _render_beautiful_table(self):
        """Render a beautiful layers table with Rich styling."""
        if not self._layers:
            empty_content = Group(
                "",
                Align.center(Text("📋 No layers available yet", style="bold yellow")),
                "",
                Align.center(
                    Text(
                        "Layers will appear here after processing steps",
                        style="dim italic",
                    )
                ),
                "",
            )

            panel = Panel(
                empty_content,
                title="[bold blue]🗂️ GEOPACKAGE LAYERS 🗂️[/]",
                border_style="yellow",
                box=box.ROUNDED,
                padding=(1, 2),
            )
            self.update(panel)
            return

        # Create clean, simple table
        table = Table(
            title="Geopackage Layers",
            box=box.SIMPLE,
            title_style="bold white",
            header_style="bold white",
            border_style="dim white",
            show_lines=False,
        )

        table.add_column("Layer Name", style="white", no_wrap=True)
        table.add_column("Features", justify="right", style="cyan")
        table.add_column("Status", justify="center")
        table.add_column("Description", style="dim white")

        # Layer descriptions for better UX
        layer_descriptions = {
            "site": "Development site boundary",
            "roads": "Road network infrastructure",
            "parcels": "Generated ownership parcels",
            "04_arterial_roads": "Main arterial routes",
            "05_secondary_roads": "Secondary street network",
            "local_roads": "Local access roads",
            "17_all_grids_merged": "Combined street blocks",
            "111_warm_block_final": "Warm development zones",
            "300_final": "Complete urban layout",
        }

        # Sort layers by importance/step order
        layer_order = [
            "site",
            "roads",
            "parcels",
            "04_arterial_roads",
            "05_secondary_roads",
            "local_roads",
            "17_all_grids_merged",
            "111_warm_block_final",
            "300_final",
        ]

        def sort_key(item):
            layer_name = item[0]
            try:
                return layer_order.index(layer_name)
            except ValueError:
                return len(layer_order)

        sorted_layers = sorted(self._layers, key=sort_key)

        for layer_name, feature_count in sorted_layers:
            if feature_count > 0:
                status = Text("✅ Ready", style="bold green")
                count_style = "bold green"
            else:
                status = Text("⚠️ Empty", style="bold yellow")
                count_style = "yellow"

            # Get description
            description = layer_descriptions.get(layer_name, "Analysis layer")

            # Format feature count with nice styling
            if feature_count > 1000:
                count_display = f"{feature_count:,}"
            else:
                count_display = str(feature_count)

            table.add_row(
                layer_name,
                count_display,
                status,
                description,
                style=count_style if feature_count > 0 else "dim",
            )

        # Create summary statistics
        total_features = sum(count for _, count in self._layers)
        total_layers = len([1 for _, count in self._layers if count > 0])

        summary = Group(
            table,
            "",
            Text(
                f"Summary: {total_layers} active layers • {total_features:,} total features",
                style="dim white",
                justify="center",
            ),
        )

        panel = Panel(
            summary,
            title="Layers",
            border_style="dim white",
            box=box.ROUNDED,
            padding=(0, 1),
        )
        self.update(panel)


class LogViewer(Static):
    """Beautiful widget for displaying log messages with Rich styling."""

    def __init__(self, max_lines: int = 100, **kwargs):
        super().__init__(**kwargs)
        self.max_lines = max_lines
        self.lines = []

    def add_log(self, message: str, level: str = "info"):
        """Add a beautifully formatted log message."""
        # Kartoza color scheme for log levels
        level_config = {
            "info": {
                "icon": "💡",
                "emoji": ":information:",
                "color": "#569FC6",
                "bg": "#E8F4FD",
            },
            "success": {
                "icon": "✨",
                "emoji": ":sparkles:",
                "color": "#06969A",
                "bg": "#F0FFF0",
            },
            "warning": {
                "icon": "⚠️",
                "emoji": ":warning:",
                "color": "#DF9E2F",
                "bg": "#FFFACD",
            },
            "error": {
                "icon": "💥",
                "emoji": ":collision:",
                "color": "#CC0403",
                "bg": "#FFE4E1",
            },
            "debug": {
                "icon": "🔍",
                "emoji": ":mag:",
                "color": "#8A8B8B",
                "bg": "#F3E5F5",
            },
            "step": {
                "icon": "🚀",
                "emoji": ":rocket:",
                "color": "#06969A",
                "bg": "#E0F7FF",
            },
        }

        config = level_config.get(level, level_config["info"])

        # Add detailed timestamp
        from datetime import datetime

        timestamp = datetime.now().strftime("%H:%M:%S.%f")[:-3]  # Include milliseconds

        # Create rich formatted line
        formatted_line = {
            "timestamp": timestamp,
            "level": level.upper(),
            "message": message,
            "config": config,
        }

        self.lines.append(formatted_line)

        # Trim if too many lines
        if len(self.lines) > self.max_lines:
            self.lines = self.lines[-self.max_lines :]

        self._render_beautiful_log()

    def _render_beautiful_log(self):
        """Render compact log content."""
        if not self.lines:
            content = Text("No activity yet", style="dim italic")
            self.update(content)
            return

        # Show more lines for detailed output
        recent_lines = self.lines[-15:]  # Show last 15 lines for better visibility
        log_content = []

        for log_entry in recent_lines:
            config = log_entry["config"]

            # Compact single line format: [time] icon message
            line = Text()
            line.append(f"[{log_entry['timestamp'][-8:]}] ", style="dim")  # Show only seconds
            line.append(f"{config['icon']} ", style=f"bold {config['color']}")
            line.append(log_entry["message"], style="white")

            log_content.append(line)

        # Simple group without panel for compactness
        from rich.console import Group

        content = Group(*log_content)
        self.update(content)

    def clear(self):
        """Clear all log messages."""
        self.lines.clear()
        content = Text("🧹 Log cleared - Ready for new activity", style="dim yellow")
        self.update(content)


class LayerTreeView(Tree):
    """Tree view widget for browsing geopackage layers."""

    class LayerSelected(Message):
        """Message sent when a layer is selected."""

        def __init__(self, layer_name: str) -> None:
            super().__init__()
            self.layer_name = layer_name

    def __init__(self, **kwargs):
        super().__init__("📦 GeoPackage Layers", **kwargs)
        self._layers = []
        self.geopackage_path = None

    def update_layers(self, geopackage_path: str, layer_info: list[tuple]):
        """Update the layers in the tree view."""
        self.geopackage_path = geopackage_path
        self._layers = layer_info
        self._rebuild_tree()

    def _rebuild_tree(self):
        """Rebuild the tree with current layer information."""
        self.clear()

        if not self._layers:
            self.root.add_leaf("📋 No layers available", data=None)
            return

        # Group layers by category
        categories = {
            "Input Data": ["site", "roads"],
            "Step 1 - Parcels": ["parcels", "sites", "roads_buffered"],
            "Step 2 - Streets": [
                "04_arterial_roads",
                "05_secondary_roads",
                "local_roads",
                "17_all_grids_merged",
            ],
            "Step 3 - Clusters": ["111_warm_block_final", "300_final", "cold", "warm"],
        }

        # Create categorized tree structure
        category_nodes = {}
        uncategorized_layers = {layer[0] for layer in self._layers}

        for category, layer_names in categories.items():
            category_node = self.root.add(f"📁 {category}")
            category_nodes[category] = category_node

            for layer_name in layer_names:
                matching_layers = [layer for layer in self._layers if layer[0] == layer_name]
                if matching_layers:
                    layer_name, feature_count = matching_layers[0]
                    if feature_count > 0:
                        icon = "✅"
                    else:
                        icon = "⚠️"

                    category_node.add_leaf(
                        f"{icon} {layer_name} ({feature_count})", data=layer_name
                    )
                    uncategorized_layers.discard(layer_name)

        # Add any remaining uncategorized layers
        if uncategorized_layers:
            other_node = self.root.add("📁 Other Layers")
            for layer_name in sorted(uncategorized_layers):
                matching_layers = [layer for layer in self._layers if layer[0] == layer_name]
                if matching_layers:
                    _, feature_count = matching_layers[0]
                    icon = "✅" if feature_count > 0 else "⚠️"
                    other_node.add_leaf(f"{icon} {layer_name} ({feature_count})", data=layer_name)

    def on_tree_node_selected(self, event: Tree.NodeSelected) -> None:
        """Handle tree node selection."""
        if event.node.data:  # Only handle leaf nodes with layer data
            self.post_message(self.LayerSelected(event.node.data))


class LayerDetailPanel(Static):
    """Panel showing detailed information about a selected layer."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.current_layer = None
        self.geopackage_path = None
        self._show_welcome()

    def _show_welcome(self):
        """Show welcome message when no layer is selected."""
        welcome_text = Group(
            "",
            Align.center(Text("🗂️ Layer Information", style="bold blue")),
            "",
            Align.center(Text("Select a layer from the tree to view details", style="dim italic")),
            "",
        )

        panel = Panel(
            welcome_text,
            title="Layer Details",
            border_style="blue",
            box=box.ROUNDED,
        )
        self.update(panel)

    def update_layer_details(self, geopackage_path: str, layer_name: str):
        """Update the panel with details for the selected layer."""
        self.geopackage_path = geopackage_path
        self.current_layer = layer_name

        try:
            # Get layer information using GDAL/OGR
            layer_info = self._get_layer_info(geopackage_path, layer_name)
            preview_path = self._generate_preview(geopackage_path, layer_name)

            # Create detailed information display
            details_table = Table(show_header=False, box=box.SIMPLE, expand=True)
            details_table.add_column("Property", style="bold cyan", width=15)
            details_table.add_column("Value", style="white")

            details_table.add_row("Layer Name", layer_name)
            details_table.add_row("Feature Type", layer_info.get("geometry_type", "Unknown"))
            details_table.add_row("Feature Count", str(layer_info.get("feature_count", 0)))
            details_table.add_row("CRS", layer_info.get("crs", "Unknown"))
            details_table.add_row("Extent", layer_info.get("extent", "Unknown"))

            # Add field information if available
            if "fields" in layer_info and layer_info["fields"]:
                fields_text = ", ".join(layer_info["fields"][:5])  # Show first 5 fields
                if len(layer_info["fields"]) > 5:
                    fields_text += f" (+{len(layer_info['fields']) - 5} more)"
                details_table.add_row("Fields", fields_text)

            content_parts = [details_table]

            # Add preview if available
            if preview_path and Path(preview_path).exists():
                content_parts.extend(
                    [
                        "",
                        Text("🖼️ Preview:", style="bold green"),
                        Text(f"Image saved to: {preview_path}", style="dim"),
                        Text(
                            "Use 'feh', 'fim' or image viewer to open",
                            style="dim italic",
                        ),
                    ]
                )

                # Try to display preview in terminal using chafa
                terminal_preview = self._get_terminal_preview(preview_path)
                if terminal_preview:
                    content_parts.extend(["", terminal_preview])

            content = Group(*content_parts)

            panel = Panel(
                content,
                title=f"[bold green]📋 {layer_name}[/]",
                border_style="green",
                box=box.ROUNDED,
            )
            self.update(panel)

        except Exception as e:
            error_content = Group(
                "",
                Align.center(Text("❌ Error Loading Layer", style="bold red")),
                "",
                Align.center(Text(str(e), style="red")),
                "",
            )

            panel = Panel(
                error_content,
                title="Error",
                border_style="red",
                box=box.ROUNDED,
            )
            self.update(panel)

    def _get_layer_info(self, geopackage_path: str, layer_name: str) -> dict:
        """Get detailed information about a layer using GDAL/OGR."""
        try:
            from osgeo import gdal, ogr

            # Enable GDAL exceptions
            gdal.UseExceptions()

            # Open the geopackage
            ds = ogr.Open(geopackage_path, 0)  # 0 = read-only
            if not ds:
                return {"error": "Could not open geopackage"}

            # Get the layer
            layer = ds.GetLayerByName(layer_name)
            if not layer:
                return {"error": f"Layer '{layer_name}' not found"}

            # Get basic information
            feature_count = layer.GetFeatureCount()
            layer_defn = layer.GetLayerDefn()
            geometry_type = ogr.GeometryTypeToName(layer_defn.GetGeomType())

            # Get spatial reference system
            srs = layer.GetSpatialRef()
            crs_info = "Unknown"
            if srs:
                try:
                    crs_info = srs.GetAuthorityCode("PROJCS") or srs.GetAuthorityCode("GEOGCS")
                    if crs_info:
                        crs_info = f"EPSG:{crs_info}"
                    else:
                        crs_info = srs.GetName() or "Custom CRS"
                except Exception:
                    crs_info = "Unknown CRS"

            # Get extent
            extent = layer.GetExtent()
            extent_str = f"({extent[0]:.2f}, {extent[2]:.2f}) to ({extent[1]:.2f}, {extent[3]:.2f})"

            # Get field information
            fields = []
            for i in range(layer_defn.GetFieldCount()):
                field_defn = layer_defn.GetFieldDefn(i)
                field_name = field_defn.GetName()
                field_type = field_defn.GetTypeName()
                fields.append(f"{field_name} ({field_type})")

            return {
                "feature_count": feature_count,
                "geometry_type": geometry_type,
                "crs": crs_info,
                "extent": extent_str,
                "fields": fields,
            }

        except Exception as e:
            return {"error": str(e)}

    def _generate_preview(self, geopackage_path: str, layer_name: str) -> str:
        """Generate a preview image of the layer using GDAL."""
        try:
            import tempfile

            # Create temporary output file
            temp_dir = Path(tempfile.gettempdir()) / "rue_tui_previews"
            temp_dir.mkdir(exist_ok=True)

            preview_path = temp_dir / f"{layer_name}_preview.png"

            # Use ogr2ogr to create a simple preview with secure temp file
            import tempfile

            with tempfile.NamedTemporaryFile(suffix=".geojson", delete=False) as temp_file:
                temp_geojson_path = temp_file.name

            cmd = [
                "ogr2ogr",
                "-f",
                "GeoJSON",
                temp_geojson_path,
                geopackage_path,
                layer_name,
            ]

            result = subprocess.run(cmd, capture_output=True, text=True)

            try:
                if result.returncode == 0:
                    # If ogr2ogr succeeded, we could use a mapping library to create a preview
                    # For now, just return the path where a preview would be saved
                    pass
            finally:
                # Clean up temporary file
                try:
                    Path(temp_geojson_path).unlink()
                except OSError:
                    pass

            return str(preview_path) if result.returncode == 0 else None

        except Exception:
            return None

    def _get_terminal_preview(self, image_path: str) -> Optional[Text]:
        """Try to get a terminal preview using chafa or similar tool."""
        try:
            # Try chafa first
            if subprocess.run(["which", "chafa"], capture_output=True).returncode == 0:
                result = subprocess.run(
                    ["chafa", "--size", "40x20", image_path],
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                if result.returncode == 0 and result.stdout:
                    return Text(result.stdout)

            # Try fim as fallback (though it won't work in terminal)
            if subprocess.run(["which", "fim"], capture_output=True).returncode == 0:
                return Text("Preview available - use 'fim' to view", style="dim italic")

            return None

        except Exception:
            return None


class SetupPanel(Container):
    """Panel for setup and configuration with file selection and preview."""

    class FileSelected(Message):
        """Message sent when a file is selected."""

        def __init__(self, file_type: str, file_path: str) -> None:
            super().__init__()
            self.file_type = file_type
            self.file_path = file_path

    class PreviewRequested(Message):
        """Message sent when a preview is requested."""

        def __init__(self, file_path: str) -> None:
            super().__init__()
            self.file_path = file_path

    class BrowseRequested(Message):
        """Message sent when browse button is pressed."""

        def __init__(self, browse_type: str) -> None:
            super().__init__()
            self.browse_type = browse_type

    def __init__(self, config, **kwargs):
        super().__init__(**kwargs)
        self.config = config

    def compose(self) -> ComposeResult:
        """Compose the setup panel in Qt-style form layout."""

        with Vertical(classes="form-container"):
            # Inputs section
            yield Label("Inputs", classes="section-header")

            # Site Boundary
            with Horizontal(classes="form-row"):
                yield Label("Site Boundary:", classes="form-label")
                yield Input(
                    value=self.config.site_path,
                    id="site_path",
                    placeholder="data/site.geojson",
                    classes="form-path",
                )
                yield Button("Preview", id="preview-site", variant="default")

            # Roads
            with Horizontal(classes="form-row"):
                yield Label("Roads:", classes="form-label")
                yield Input(
                    value=self.config.roads_path,
                    id="roads_path",
                    placeholder="data/roads.geojson",
                    classes="form-path",
                )
                yield Button("Preview", id="preview-roads", variant="default")

            # Outputs section
            yield Label("Outputs", classes="section-header")

            # Output folder
            with Horizontal(classes="form-row"):
                yield Label("Folder:", classes="form-label")
                yield Input(
                    value=self.config.output_dir,
                    id="output_dir",
                    placeholder="outputs",
                    classes="form-path",
                )
                yield Button("Browse", id="browse-output", variant="default")

            # Options section
            yield Label("Options", classes="section-header")

            # Store intermediate results
            with Horizontal(classes="form-row"):
                yield Label("", classes="form-label")
                yield Checkbox(
                    "Store intermediate results",
                    value=self.config.save_intermediate,
                    id="save_intermediate",
                )
                yield Static("")  # Spacer

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses."""
        button_id = event.button.id

        if button_id == "preview-site":
            site_path = self.query_one("#site_path", Input).value
            if site_path and Path(site_path).exists():
                self.post_message(self.PreviewRequested(site_path))
        elif button_id == "preview-roads":
            roads_path = self.query_one("#roads_path", Input).value
            if roads_path and Path(roads_path).exists():
                self.post_message(self.PreviewRequested(roads_path))
        elif button_id == "browse-output":
            # Request file browser for output directory
            self.post_message(self.BrowseRequested("output"))

    def on_input_changed(self, event: Input.Changed) -> None:
        """Handle input changes."""
        if event.input.id == "output_dir":
            # Update config when output directory is manually changed
            new_output_dir = event.value.strip()
            if new_output_dir and new_output_dir != self.config.output_dir:
                self.config.output_dir = new_output_dir
                self.config.update_output_paths()

    def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
        """Handle checkbox changes."""
        pass


class FileBrowserWidget(Container):
    """File browser widget for selecting directories and files."""

    class DirectorySelected(Message):
        """Message sent when a directory is selected."""

        def __init__(self, directory_path: str) -> None:
            super().__init__()
            self.directory_path = directory_path

    class FileSelected(Message):
        """Message sent when a file is selected."""

        def __init__(self, file_path: str) -> None:
            super().__init__()
            self.file_path = file_path

    def __init__(self, initial_path: str = ".", **kwargs):
        super().__init__(**kwargs)
        self.current_path = Path(initial_path).resolve()

    def compose(self) -> ComposeResult:
        """Compose the file browser interface."""
        with Vertical():
            # Current path display and navigation
            with Horizontal(classes="form-row"):
                yield Label("Current Path:", classes="form-label")
                yield Input(value=str(self.current_path), id="current-path", disabled=True)
                yield Button("↑", id="up-dir", variant="default")
                yield Button("🏠", id="home-dir", variant="default")

            # Directory tree browser
            yield DirectoryTree(str(self.current_path), id="directory-tree")

            # Action buttons
            with Horizontal(classes="form-row"):
                yield Button("Select Directory", id="select-dir", variant="primary")
                yield Button("Create New Dir", id="create-dir", variant="default")
                yield Input(placeholder="New directory name", id="new-dir-name")

    def on_mount(self) -> None:
        """Initialize the file browser."""
        self._update_current_path_display()

    def on_directory_tree_directory_selected(self, event: DirectoryTree.DirectorySelected) -> None:
        """Handle directory selection in the tree."""
        self.current_path = Path(event.path)
        self._update_current_path_display()

    def on_directory_tree_file_selected(self, event: DirectoryTree.FileSelected) -> None:
        """Handle file selection in the tree."""
        self.post_message(self.FileSelected(event.path))

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses."""
        button_id = event.button.id

        if button_id == "up-dir":
            if self.current_path.parent != self.current_path:  # Not at root
                self.current_path = self.current_path.parent
                self._update_directory_tree()
                self._update_current_path_display()

        elif button_id == "home-dir":
            self.current_path = Path.home()
            self._update_directory_tree()
            self._update_current_path_display()

        elif button_id == "select-dir":
            self.post_message(self.DirectorySelected(str(self.current_path)))

        elif button_id == "create-dir":
            new_dir_input = self.query_one("#new-dir-name", Input)
            new_dir_name = new_dir_input.value.strip()
            if new_dir_name:
                try:
                    new_dir_path = self.current_path / new_dir_name
                    new_dir_path.mkdir(exist_ok=True)
                    new_dir_input.value = ""
                    self._update_directory_tree()
                except Exception:
                    # Could show an error message here
                    pass

    def _update_current_path_display(self) -> None:
        """Update the current path display."""
        path_input = self.query_one("#current-path", Input)
        path_input.value = str(self.current_path)

    def _update_directory_tree(self) -> None:
        """Update the directory tree with the current path."""
        try:
            directory_tree = self.query_one("#directory-tree", DirectoryTree)
            directory_tree.path = str(self.current_path)
        except Exception:
            # Tree might not be mounted yet
            pass


class GeojsonPreviewWidget(Static):
    """Widget for displaying GeoJSON file previews using matplotlib."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.current_file = None

    def show_preview(self, geojson_path: str):
        """Generate and display a preview of the GeoJSON file."""
        self.current_file = geojson_path

        try:
            preview_image = self._generate_matplotlib_preview(geojson_path)
            if preview_image:
                self._display_preview(preview_image)
            else:
                self._show_error("Failed to generate preview")

        except Exception as e:
            self._show_error(f"Preview error: {str(e)}")

    def _generate_matplotlib_preview(self, geojson_path: str) -> Optional[str]:
        """Generate a preview image using matplotlib."""
        try:
            import tempfile

            import geopandas as gpd
            import matplotlib.pyplot as plt

            # Read the GeoJSON file
            gdf = gpd.read_file(geojson_path)

            if gdf.empty:
                return None

            # Create the plot
            fig, ax = plt.subplots(figsize=(10, 8))

            # Plot the data
            if gdf.geom_type.iloc[0] in ["Point", "MultiPoint"]:
                gdf.plot(ax=ax, markersize=50, alpha=0.7, color="red")
            elif gdf.geom_type.iloc[0] in ["LineString", "MultiLineString"]:
                gdf.plot(ax=ax, linewidth=2, alpha=0.8, color="blue")
            else:  # Polygon, MultiPolygon
                gdf.plot(ax=ax, alpha=0.7, facecolor="lightblue", edgecolor="blue")

            # Style the plot
            ax.set_title(f"Preview: {Path(geojson_path).name}", fontsize=14, fontweight="bold")
            ax.set_aspect("equal")

            # Remove axes for cleaner look
            ax.set_xticks([])
            ax.set_yticks([])

            # Add feature count
            feature_count = len(gdf)
            geom_type = gdf.geom_type.iloc[0] if not gdf.empty else "Unknown"
            ax.text(
                0.02,
                0.98,
                f"Features: {feature_count}\nType: {geom_type}",
                transform=ax.transAxes,
                fontsize=10,
                verticalalignment="top",
                bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.8},
            )

            # Save to temporary file
            temp_dir = Path(tempfile.gettempdir()) / "rue_tui_previews"
            temp_dir.mkdir(exist_ok=True)

            preview_path = temp_dir / f"{Path(geojson_path).stem}_preview.png"
            plt.savefig(
                preview_path,
                dpi=150,
                bbox_inches="tight",
                facecolor="white",
                edgecolor="none",
            )
            plt.close()

            return str(preview_path)

        except Exception as e:
            print(f"Preview generation error: {e}")
            return None

    def _display_preview(self, image_path: str):
        """Display the preview image using terminal tools."""
        content_parts = [
            Text("🖼️ GeoJSON Preview Generated", style="bold green"),
            "",
            Text(f"File: {Path(self.current_file).name}", style="bold white"),
            Text(f"Preview: {image_path}", style="dim"),
            "",
        ]

        # Try to show in terminal
        terminal_preview = self._get_terminal_preview(image_path)
        if terminal_preview:
            content_parts.append(terminal_preview)
        else:
            content_parts.extend(
                [
                    Text("Use one of these commands to view:", style="yellow"),
                    Text(f"fim {image_path}", style="cyan"),
                    Text(f"feh {image_path}", style="cyan"),
                    Text(f"eog {image_path}", style="cyan"),
                ]
            )

        content = Group(*content_parts)
        panel = Panel(
            content,
            title="Preview",
            border_style="green",
            box=box.ROUNDED,
        )
        self.update(panel)

    def _show_error(self, error_msg: str):
        """Show an error message."""
        content = Group(
            "",
            Align.center(Text("❌ Preview Error", style="bold red")),
            "",
            Align.center(Text(error_msg, style="red")),
            "",
        )
        panel = Panel(
            content,
            title="Error",
            border_style="red",
            box=box.ROUNDED,
        )
        self.update(panel)

    def _get_terminal_preview(self, image_path: str) -> Optional[Text]:
        """Try to get a terminal preview using chafa."""
        try:
            if subprocess.run(["which", "chafa"], capture_output=True).returncode == 0:
                result = subprocess.run(
                    ["chafa", "--size", "60x30", image_path],
                    capture_output=True,
                    text=True,
                    timeout=10,
                )
                if result.returncode == 0 and result.stdout:
                    return Text(result.stdout)
            return None
        except Exception:
            return None
