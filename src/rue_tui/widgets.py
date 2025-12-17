"""Custom widgets for the RUE TUI application."""

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
    Input,
    Label,
    Static,
)


class StepCard(Static):
    """A card widget representing a workflow step with beautiful Rich styling."""

    def __init__(
        self, step_num: int, title: str, description: str, status: str = "pending", **kwargs
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

    def update_progress(self, completed: int, total: int, description: str = ""):
        """Update with beautiful animated progress bar."""
        percentage = (completed / total) * 100 if total > 0 else 0

        # Create animated progress bar
        filled = int(percentage // 3.33)  # 30 char bar
        bar = "█" * filled + "▓" * (1 if filled < 30 else 0) + "░" * (29 - filled)

        # Color the progress bar based on completion
        if percentage < 30:
            bar_style = "red on dark_red"
        elif percentage < 70:
            bar_style = "yellow on dark_blue"
        else:
            bar_style = "green on dark_green"

        # Beautiful percentage display
        percent_text = Text(f"{percentage:.1f}%", style="bold white")

        # Status messages
        if percentage < 10:
            status = "🌱 Just getting started..."
        elif percentage < 50:
            status = "🏗️ Making good progress..."
        elif percentage < 90:
            status = "🚀 Almost there..."
        else:
            status = "✨ Finishing up..."

        progress_content = Group(
            Align.center(Text(description, style="bold white")),
            "",
            Align.center(Text(bar, style=bar_style)),
            "",
            Align.center(percent_text),
            "",
            Align.center(Text(f"({completed}/{total})", style="dim")),
            "",
            Align.center(Text(status, style="italic cyan")),
        )

        content = Panel(
            progress_content,
            title="[bold blue]🚀 PROCESSING 🚀[/]",
            border_style="bold blue",
            box=box.DOUBLE,
            padding=(1, 2),
        )
        self.update(content)

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
                    Text("Layers will appear here after processing steps", style="dim italic")
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
            summary, title="Layers", border_style="dim white", box=box.ROUNDED, padding=(0, 1)
        )
        self.update(panel)


class ConfigPanel(Container):
    """Panel for configuration options."""

    class ConfigChanged(Message):
        """Message sent when configuration changes."""

        def __init__(self, key: str, value: str) -> None:
            super().__init__()
            self.key = key
            self.value = value

    def __init__(self, config, **kwargs):
        super().__init__(**kwargs)
        self.config = config

    def compose(self) -> ComposeResult:
        """Compose the config panel."""
        yield Label("Configuration", classes="config-title")

        with Vertical(classes="config-form"):
            # File paths
            yield Label("Input Files:", classes="config-section")
            yield Input(
                placeholder="Site GeoJSON path", value=self.config.site_path, id="site_path"
            )
            yield Input(
                placeholder="Roads GeoJSON path", value=self.config.roads_path, id="roads_path"
            )

            # Output settings
            yield Label("Output Settings:", classes="config-section")
            yield Input(
                placeholder="Output directory", value=self.config.output_dir, id="output_dir"
            )

            # Visualization settings
            yield Label("Visualization:", classes="config-section")
            yield Horizontal(
                Input(placeholder="Width", value=str(self.config.image_width), id="image_width"),
                Input(placeholder="Height", value=str(self.config.image_height), id="image_height"),
                classes="config-row",
            )

            # Options
            yield Checkbox("Auto-advance steps", value=self.config.auto_advance, id="auto_advance")
            yield Checkbox(
                "Save intermediate results",
                value=self.config.save_intermediate,
                id="save_intermediate",
            )

    def on_input_changed(self, event: Input.Changed) -> None:
        """Handle input changes."""
        self.post_message(self.ConfigChanged(event.input.id, event.value))

    def on_checkbox_changed(self, event: Checkbox.Changed) -> None:
        """Handle checkbox changes."""
        self.post_message(self.ConfigChanged(event.checkbox.id, str(event.value)))


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
            "info": {"icon": "💡", "emoji": ":information:", "color": "#569FC6", "bg": "#E8F4FD"},
            "success": {"icon": "✨", "emoji": ":sparkles:", "color": "#06969A", "bg": "#F0FFF0"},
            "warning": {"icon": "⚠️", "emoji": ":warning:", "color": "#DF9E2F", "bg": "#FFFACD"},
            "error": {"icon": "💥", "emoji": ":collision:", "color": "#CC0403", "bg": "#FFE4E1"},
            "debug": {"icon": "🔍", "emoji": ":mag:", "color": "#8A8B8B", "bg": "#F3E5F5"},
            "step": {"icon": "🚀", "emoji": ":rocket:", "color": "#06969A", "bg": "#E0F7FF"},
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


class ActionButton(Button):
    """Enhanced button with status indication."""

    def __init__(self, label: str, action: str, enabled: bool = True, **kwargs):
        super().__init__(label, **kwargs)
        self.action = action
        self._enabled = enabled
        self._success = False
        self.update_appearance()

    def update_appearance(self):
        """Update button appearance based on state."""
        if self._success:
            self.variant = "success"
            self.disabled = False
        elif self._enabled:
            self.variant = "primary"
            self.disabled = False
        else:
            self.variant = "default"
            self.disabled = True

    def enable(self):
        """Enable the button."""
        self._enabled = True
        self.update_appearance()

    def disable(self):
        """Disable the button."""
        self._enabled = False
        self.update_appearance()

    def set_success(self):
        """Set button to success state (green)."""
        self._success = True
        self._enabled = True
        self.update_appearance()

    async def _on_click(self, event) -> None:
        """Handle button click and add debug logging."""
        # Add debug logging to see if clicks are detected
        print(f"DEBUG: ActionButton {self.id} clicked!")
        # Let the default button behavior handle the rest
        await super()._on_click(event)
