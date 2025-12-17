"""Main TUI application for RUE-lib."""

import asyncio
import shutil
from pathlib import Path

from rich.text import Text
from textual import events
from textual.app import App, ComposeResult
from textual.containers import Container, Horizontal, Vertical
from textual.reactive import reactive
from textual.screen import ModalScreen
from textual.widgets import Button, Footer, Header, Input, Label, Static, TabbedContent, TabPane

from rue_lib.cluster.runner import ClusterConfig, generate_clusters

# RUE-lib imports
from rue_lib.config import MainConfig
from rue_lib.site.runner import SiteConfig, generate_parcels
from rue_lib.streets.runner import StreetConfig, generate_streets

# TUI imports
from .config import TuiConfig
from .visualizer import MapVisualizer
from .widgets import (
    GeojsonPreviewWidget,
    LayerBrowser,
    LayerDetailPanel,
    LayerTreeView,
    LogViewer,
    ProgressDisplay,
    SetupPanel,
    StepCard,
)


class ImageViewerScreen(ModalScreen):
    """Modal screen for viewing rendered maps."""

    def __init__(self, image_path: str, title: str = "Map View"):
        super().__init__()
        self.image_path = image_path
        self.title = title

    def compose(self) -> ComposeResult:
        with Container(id="image-viewer"):
            yield Label(self.title, classes="viewer-title")
            yield Label(f"Image: {self.image_path}", classes="viewer-path")
            yield Label("Press 'q' or ESC to close", classes="viewer-help")
            yield Button("Close", variant="primary", id="close-viewer")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "close-viewer":
            self.dismiss()

    def on_key(self, event: events.Key) -> None:
        if event.key in ("q", "escape"):
            self.dismiss()


class RueTuiApp(App):
    """Main RUE TUI application."""

    CSS = """
    /* Compact layout optimized for small screens */
    .title {
        text-align: center;
        padding: 0 1;
        text-style: bold;
        height: 2;
    }

    .subtitle {
        text-align: center;
        color: $accent;
        text-style: italic;
        padding: 0;
        height: 1;
    }

    .welcome-banner {
        padding: 1;
        margin: 0;
        height: 4;
    }

    /* Compact step cards */
    .step-container {
        padding: 0;
        margin: 0;
        height: auto;
    }

    .step-card {
        margin: 0 0 1 0;
        height: 3;
    }

    /* Compact config panel */
    .config-title {
        text-style: bold;
        color: $accent;
    }

    .config-section {
        margin: 0;
        padding: 0;
        text-style: bold;
        color: $primary;
    }

    .config-form {
        padding: 0;
    }

    .config-row {
        height: 2;
        margin: 0;
        padding: 0;
    }

    /* Compact progress display */
    .progress-container {
        padding: 0;
        height: 4;
    }

    /* Compact action buttons */
    .action-bar {
        padding: 0;
        margin: 0;
        height: 3;
    }

    /* Expanded log viewer for detailed output */
    .log-container {
        padding: 1;
        height: 1fr;
    }

    /* Compact layer browser */
    .layer-container {
        padding: 0;
        height: 6;
    }

    /* Compact image viewer modal */
    #image-viewer {
        background: $surface;
        border: round $primary;
        width: 90%;
        height: 80%;
        margin: 1;
        padding: 1;
    }

    .viewer-title {
        text-style: bold;
        text-align: center;
        color: $primary;
    }

    .viewer-path {
        color: $text-muted;
        text-align: center;
    }

    .viewer-help {
        text-align: center;
        color: $accent;
        margin: 1 0;
    }

    /* Rounded right edge for final button (viewer button) */
    #btn-view {
        border: round $primary;
    }

    /* Results tab layout */
    .results-main {
        height: 1fr;
        padding: 0;
    }

    .tree-container {
        width: 1fr;
        min-width: 30;
        max-width: 40%;
        border-right: solid $primary;
        padding: 0 1 0 0;
    }

    .detail-container {
        width: 2fr;
        padding: 0 0 0 1;
    }

    .hidden {
        display: none;
    }

    /* Qt-style form layout */
    .form-container {
        padding: 0;
        height: auto;
    }

    .section-header {
        text-style: bold;
        color: $primary;
        margin: 0 0 0 0;
        height: 1;
    }

    .form-row {
        height: 2;
        margin: 0;
        padding: 0;
        align: left middle;
    }

    .form-label {
        width: 16;
        text-style: bold;
        color: $accent;
        text-align: right;
        padding: 0 1 0 0;
    }

    /* Setup content area */
    .setup-content {
        height: 12;
        max-height: 12;
        margin: 0;
        padding: 0;
        overflow-y: auto;
    }

    .setup-content.hidden {
        display: none;
    }
    """

    TITLE = "🏙️ RUE TUI - Rapid Urbanisation Explorer"
    SUB_TITLE = "Terminal Interface for Urban Planning Analysis"

    # Reactive state
    current_step: reactive[int] = reactive(0)

    def __init__(self):
        super().__init__()
        self.config = TuiConfig()
        self.visualizer = MapVisualizer(
            width=self.config.image_width, height=self.config.image_height, dpi=self.config.dpi
        )

        # State tracking
        self.step_status = {0: "pending", 1: "pending", 2: "pending", 3: "pending"}
        self.results = {}  # Store step results

        # UI components
        self.step_cards = {}
        self.progress_display = None
        self.log_viewer = None
        self.layer_browser = None
        self.layer_tree = None
        self.layer_detail_panel = None
        self.setup_panel = None
        self.preview_widget = None

    def compose(self) -> ComposeResult:
        """Compose the main UI."""
        yield Header()

        with TabbedContent(initial="workflow"):
            # Main workflow tab
            with TabPane("Workflow", id="workflow"):
                # Compact welcome banner
                yield self._create_welcome_banner()

                # Compact step cards
                with Vertical(classes="step-container"):
                    self.step_cards[0] = StepCard(
                        0, "Setup Project", "Configure inputs and settings", status="pending"
                    )
                    yield self.step_cards[0]

                    self.step_cards[1] = StepCard(
                        1, "Generate Parcels", "Create ownership parcels", status="pending"
                    )
                    yield self.step_cards[1]

                    self.step_cards[2] = StepCard(
                        2, "Generate Streets", "Create street networks", status="pending"
                    )
                    yield self.step_cards[2]

                    self.step_cards[3] = StepCard(
                        3, "Generate Clusters", "Design urban clusters", status="pending"
                    )
                    yield self.step_cards[3]

                # Step 0 setup content area
                with Container(id="setup-content", classes="setup-content"):
                    self.setup_panel = SetupPanel(self.config)
                    yield self.setup_panel

                    # Preview widget for GeoJSON files
                    self.preview_widget = GeojsonPreviewWidget()
                    yield self.preview_widget

                # Compact progress display
                with Container(classes="progress-container"):
                    self.progress_display = ProgressDisplay()
                    yield self.progress_display

                # Compact action buttons
                with Horizontal(classes="action-bar"):
                    yield Button("▶️ Step 0", variant="primary", id="btn-step0")
                    yield Button("▶️ Step 1", variant="default", disabled=True, id="btn-step1")
                    yield Button("▶️ Step 2", variant="default", disabled=True, id="btn-step2")
                    yield Button("▶️ Step 3", variant="default", disabled=True, id="btn-step3")
                    yield Button("👁️ View", variant="default", disabled=True, id="btn-view")

            # Results tab with split layout
            with TabPane("Results", id="results"):
                # Main content area with tree view and layer details
                with Horizontal(classes="results-main"):
                    # Left 1/3: Tree view
                    with Container(classes="tree-container"):
                        self.layer_tree = LayerTreeView()
                        yield self.layer_tree

                    # Right 2/3: Layer details panel
                    with Container(classes="detail-container"):
                        self.layer_detail_panel = LayerDetailPanel()
                        yield self.layer_detail_panel

                # Compact results actions
                with Horizontal(classes="action-bar"):
                    yield Button("📊 S1", id="view-step1", disabled=True)
                    yield Button("📊 S2", id="view-step2", disabled=True)
                    yield Button("📊 S3", id="view-step3", disabled=True)
                    yield Button("📁 Folder", id="open-folder")
                    yield Button("🧹 Clear", id="clear-results")

                # Keep the old layer browser for compatibility (hidden)
                with Container(classes="hidden"):
                    self.layer_browser = LayerBrowser()
                    yield self.layer_browser

            # Logs tab
            with TabPane("Logs", id="logs"):
                with Container(classes="log-container"):
                    self.log_viewer = LogViewer()
                    yield self.log_viewer

        yield Footer()

    def _create_welcome_banner(self) -> Static:
        """Create a beautiful welcome banner using Rich."""
        from rich.align import Align
        from rich.console import Group

        # Create gradient title
        title = Text()
        title.append("🏙️ ", style="bold #DF9E2F")
        title.append("RUE", style="bold #CC0403")  # Kartoza red
        title.append("-", style="bold white")
        title.append("lib", style="bold #06969A")  # Kartoza teal
        title.append(" Urban Planning Workflow", style="bold #569FC6")  # Kartoza blue
        title.append(" 🏙️", style="bold #DF9E2F")

        # Create subtitle with emojis
        subtitle = Text()
        subtitle.append("✨ ", style="yellow")
        subtitle.append("Transform your urban development ideas into reality", style="italic cyan")
        subtitle.append(" ✨", style="yellow")

        # Create feature highlights
        features = Text()
        features.append("🏗️ Smart Parcel Generation  ", style="green")
        features.append("🛣️ Intelligent Street Networks  ", style="blue")
        features.append("🏘️ Optimized Urban Clusters", style="purple")

        # Compact title for small screens
        title = Text()
        title.append("🏙️ ", style="bold #DF9E2F")
        title.append("RUE", style="bold #CC0403")
        title.append("-", style="bold white")
        title.append("lib", style="bold #06969A")
        title.append(" Planning", style="bold #569FC6")

        # Simple subtitle
        subtitle = Text("Parcels → Streets → Clusters", style="italic #8A8B8B")

        # Minimal layout
        content = Group(Align.center(title), Align.center(subtitle))

        return Static(content, classes="welcome-banner")

    def on_mount(self) -> None:
        """Initialize the app with beautiful welcome messages."""
        self.log_viewer.add_log("Welcome to RUE TUI! 🎉", "step")
        self.log_viewer.add_log("Urban planning made beautiful and interactive", "info")
        self.log_viewer.add_log(f"Output directory configured: {self.config.output_dir}", "info")

        # Show setup content initially for Step 0
        self._update_setup_visibility()

        # Check if input files exist with detailed messages
        site_exists = Path(self.config.site_path).exists()
        roads_exists = Path(self.config.roads_path).exists()

        if site_exists and roads_exists:
            self.log_viewer.add_log("All input files found and ready! 🎯", "success")
            self.log_viewer.add_log("You can start with Step 1: Generate Parcels", "info")
        else:
            missing = []
            if not site_exists:
                missing.append(f"site boundary: {self.config.site_path}")
            if not roads_exists:
                missing.append(f"roads network: {self.config.roads_path}")
            self.log_viewer.add_log(f"Missing input files: {', '.join(missing)}", "warning")
            self.log_viewer.add_log("Demo data will be created automatically when needed", "info")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        """Handle button presses."""
        button_id = event.button.id
        print(f"DEBUG: on_button_pressed called with button_id: {button_id}")

        # Debug popup to confirm click detection
        self.notify(f"Button clicked: {button_id}")
        print(f"DEBUG: Notification sent for button: {button_id}")

        if button_id == "btn-step0":
            event.button.disabled = True
            event.button.label = "🔄 Running..."
            self.log_viewer.add_log("🚀 Step 0 started - Setup and configuration", "step")
            self.run_worker(self.run_step(0))
        elif button_id == "btn-step1":
            event.button.disabled = True
            event.button.label = "🔄 Running..."
            self.log_viewer.add_log("🚀 Step 1 started - Generating parcels", "step")
            self.run_worker(self.run_step(1))
        elif button_id == "btn-step2":
            event.button.disabled = True
            event.button.label = "🔄 Running..."
            self.log_viewer.add_log("🚀 Step 2 started - Creating streets", "step")
            self.run_worker(self.run_step(2))
        elif button_id == "btn-step3":
            event.button.disabled = True
            event.button.label = "🔄 Running..."
            self.log_viewer.add_log("🚀 Step 3 started - Designing clusters", "step")
            self.run_worker(self.run_step(3))
        elif button_id == "btn-view":
            self.view_current_results()
        elif button_id == "btn-browse":
            self.browse_layers()
        elif button_id.startswith("view-step"):
            step_num = int(button_id.split("-")[1][-1])
            self.view_step_results(step_num)
        elif button_id == "open-folder":
            self.open_output_folder()
        elif button_id == "clear-results":
            self.clear_results()

    def on_layer_tree_view_layer_selected(self, event: LayerTreeView.LayerSelected) -> None:
        """Handle layer selection from tree view."""
        if Path(self.config.geopackage_path).exists():
            self.layer_detail_panel.update_layer_details(
                self.config.geopackage_path, event.layer_name
            )
            self.log_viewer.add_log(f"🗂️ Selected layer: {event.layer_name}", "info")

    def on_setup_panel_file_selected(self, event: SetupPanel.FileSelected) -> None:
        """Handle file selection from setup panel."""
        self.log_viewer.add_log(f"📁 Selected {event.file_type} file: {event.file_path}", "info")

        # Update config with new file path
        if event.file_type == "site":
            self.config.site_path = event.file_path
        elif event.file_type == "roads":
            self.config.roads_path = event.file_path

    def on_setup_panel_preview_requested(self, event: SetupPanel.PreviewRequested) -> None:
        """Handle preview request from setup panel."""
        self.log_viewer.add_log(f"🖼️ Generating preview for: {event.file_path}", "info")

        if self.preview_widget:
            self.preview_widget.show_preview(event.file_path)
            # The preview will show in the setup content area
            self.log_viewer.add_log("Preview generated in setup area", "info")

    async def run_step(self, step: int):
        """Run a specific step of the analysis."""
        try:
            self.step_status[step] = "running"
            self.step_cards[step].update_status("running")
            self.progress_display.start_progress(f"Running Step {step}")
            self.log_viewer.add_log(f"▶️ Starting Step {step}", "info")

            if step == 0:
                result = await self._run_step0()
            elif step == 1:
                result = await self._run_step1()
            elif step == 2:
                result = await self._run_step2()
            elif step == 3:
                result = await self._run_step3()
            else:
                raise ValueError(f"Invalid step: {step}")

            # Success
            self.step_status[step] = "completed"
            self.step_cards[step].update_status("completed")
            self.results[step] = result
            self.progress_display.complete_progress(f"Step {step} completed")
            self.log_viewer.add_log(f"✅ Step {step} completed successfully", "success")

            # Restore button state and change to green
            current_btn = self.query_one(f"#btn-step{step}")
            current_btn.label = f"✅ Step {step}"
            current_btn.variant = "success"
            current_btn.disabled = False

            # Enable next step and focus it
            if step < 3:
                next_btn = self.query_one(f"#btn-step{step + 1}")
                next_btn.disabled = False
                next_btn.variant = "primary"
                # Auto-focus the next step button
                next_btn.focus()
            else:
                # After step 3, focus the viewer button
                view_btn = self.query_one("#btn-view")
                view_btn.focus()

            # Enable view buttons (only for steps 1, 2, 3)
            if step > 0:
                view_btn = self.query_one("#btn-view")
                view_btn.disabled = False
                self.query_one(f"#view-step{step}").disabled = False

            # Update layer browser
            self._update_layer_browser()

            # Update setup visibility
            self._update_setup_visibility()

            # Auto-advance if enabled
            if self.config.auto_advance and step < 3:
                await asyncio.sleep(2)  # Brief pause
                await self.run_step(step + 1)

        except Exception as e:
            self.step_status[step] = "error"
            self.step_cards[step].update_status("error")
            self.progress_display.error_progress(f"Step {step} failed: {str(e)}")
            self.log_viewer.add_log(f"❌ Step {step} failed: {str(e)}", "error")

            # Restore button state on error
            current_btn = self.query_one(f"#btn-step{step}")
            current_btn.disabled = False
            current_btn.variant = "error"
            current_btn.label = f"❌ Step {step} (retry)"

    async def _run_step0(self) -> str:
        """Run Step 0: Setup and Configuration."""
        self.log_viewer.add_log("📋 Validating setup and configuration", "info")

        # Simulate setup process
        await asyncio.sleep(1)  # Brief delay for UI feedback

        # Check if required files are set and exist
        site_path = (
            self.setup_panel.query_one("#site_path", Input).value
            if self.setup_panel
            else self.config.site_path
        )
        roads_path = (
            self.setup_panel.query_one("#roads_path", Input).value
            if self.setup_panel
            else self.config.roads_path
        )

        missing_files = []
        if not site_path or not Path(site_path).exists():
            missing_files.append("site boundary file")
        if not roads_path or not Path(roads_path).exists():
            missing_files.append("roads network file")

        if missing_files:
            self.log_viewer.add_log(
                f"⚠️ Missing required files: {', '.join(missing_files)}", "warning"
            )
            self.log_viewer.add_log("Demo data will be created automatically when needed", "info")
        else:
            self.log_viewer.add_log("✅ All required files are available", "success")

        self.log_viewer.add_log("🔧 Configuration validated successfully", "success")
        return "setup_complete"

    async def _run_step1(self) -> str:
        """Run Step 1: Generate Parcels."""
        self.log_viewer.add_log("📋 Configuring Step 1 parameters", "info")

        config = SiteConfig(
            site_path=self.config.site_path,
            roads_path=self.config.roads_path,
            output_dir=self.config.step1_output,
            geopackage_path=self.config.geopackage_path,
            road_arterial_width_m=MainConfig.road_arterial_width_m,
            road_secondary_width_m=MainConfig.road_secondary_width_m,
        )

        self.log_viewer.add_log(f"📁 Site path: {self.config.site_path}", "info")
        self.log_viewer.add_log(f"🛣️ Roads path: {self.config.roads_path}", "info")
        self.log_viewer.add_log(f"📤 Output directory: {self.config.step1_output}", "info")

        # Update the ui to show we are processing
        self.progress_display.update("Generating parcels for step1 ...")
        self.log_viewer.add_log("🔄 Starting parcel generation process...", "info")

        # Capture detailed progress with a wrapper function
        def generate_parcels_with_logging(config):
            try:
                # Log major processing steps
                self.call_from_thread(self.log_viewer.add_log, "📂 Reading input files...", "info")

                # Import here to capture any import-time logging
                from rue_lib.geo import to_metric_crs
                from rue_lib.site.io import read_roads, read_site

                self.call_from_thread(self.log_viewer.add_log, "🗺️ Loading site boundary", "info")
                site = read_site(config.site_path)

                self.call_from_thread(self.log_viewer.add_log, "🛣️ Loading roads network", "info")
                roads = read_roads(config.roads_path)

                self.call_from_thread(
                    self.log_viewer.add_log, "🔧 Converting to metric CRS if needed", "info"
                )
                if site.crs and site.crs.is_projected:
                    site_m = site
                else:
                    site_m = to_metric_crs(site)

                if roads.crs and roads.crs.is_projected:
                    roads_m = roads
                else:
                    roads_m = to_metric_crs(roads)

                self.call_from_thread(
                    self.log_viewer.add_log, f"📊 Site features: {len(site_m)}", "info"
                )
                self.call_from_thread(
                    self.log_viewer.add_log, f"🛣️ Road features: {len(roads_m)}", "info"
                )

                self.call_from_thread(
                    self.log_viewer.add_log, "⚡ Processing road buffers...", "info"
                )

                # Call the actual generate_parcels function
                result = generate_parcels(config)

                self.call_from_thread(
                    self.log_viewer.add_log, f"✅ Generated parcels saved to: {result}", "success"
                )
                return result

            except Exception as e:
                self.call_from_thread(
                    self.log_viewer.add_log, f"❌ Error in parcel generation: {str(e)}", "error"
                )
                raise

        # Run in thread to avoid blocking UI
        result = await asyncio.get_event_loop().run_in_executor(
            None, generate_parcels_with_logging, config
        )

        self.log_viewer.add_log("🎉 Step 1 completed successfully!", "success")
        return str(result)

    async def _run_step2(self) -> str:
        """Run Step 2: Generate Streets."""
        # Get parcels from step 1
        if 1 not in self.results:
            raise RuntimeError("Step 1 must be completed first")

        config = StreetConfig(
            parcel_path=self.results[1],
            roads_path=self.config.roads_path,
            output_dir=self.config.step2_output,
            geopackage_path=self.config.geopackage_path,
            road_arterial_width_m=MainConfig.road_arterial_width_m,
            road_secondary_width_m=MainConfig.road_secondary_width_m,
            road_locals_width_m=MainConfig.road_local_width_m,
            part_art_d=MainConfig.on_grid_partition_depth_arterial_roads,
            part_sec_d=MainConfig.on_grid_partition_depth_secondary_roads,
            part_loc_d=MainConfig.on_grid_partition_depth_local_roads,
            on_grid_partition_depth_arterial_roads=MainConfig.on_grid_partition_depth_arterial_roads,
            on_grid_partition_depth_secondary_roads=MainConfig.on_grid_partition_depth_secondary_roads,
            off_grid_cluster_depth=MainConfig.off_grid_cluster_depth,
            off_grid_cluster_width=MainConfig.off_grid_cluster_width,
            off_grid_arterial_clusters_depth=MainConfig.off_grid_arterial_clusters_depth,
            off_grid_secondary_clusters_depth=MainConfig.off_grid_secondary_clusters_depth,
            off_grid_local_clusters_depth=MainConfig.off_grid_local_clusters_depth,
            off_grid_local_clusters_width=MainConfig.off_grid_local_clusters_width,
            sidewalk_width_m=MainConfig.sidewalk_width_m,
        )

        result = await asyncio.get_event_loop().run_in_executor(None, generate_streets, config)

        return str(result)

    async def _run_step3(self) -> str:
        """Run Step 3: Generate Clusters."""
        # Check step 2 results
        if 2 not in self.results:
            raise RuntimeError("Step 2 must be completed first")

        # Look for the all_grids_merged.geojson file
        step2_output_dir = Path(self.config.step2_output)
        grids_file = step2_output_dir / "all_grids_merged.geojson"

        if not grids_file.exists():
            raise FileNotFoundError(f"Required file not found: {grids_file}")

        config = ClusterConfig(
            roads_path=self.config.roads_path,
            input_path=str(grids_file),
            output_dir=self.config.step3_output,
            geopackage_path=self.config.geopackage_path,
            road_arterial_width_m=MainConfig.road_arterial_width_m,
            road_secondary_width_m=MainConfig.road_secondary_width_m,
            road_local_width_m=MainConfig.road_local_width_m,
            on_grid_partition_depth_arterial_roads=MainConfig.on_grid_partition_depth_arterial_roads,
            on_grid_partition_depth_secondary_roads=MainConfig.on_grid_partition_depth_secondary_roads,
            off_grid_cluster_depth=MainConfig.off_grid_cluster_depth,
            off_grid_cluster_width=MainConfig.off_grid_cluster_width,
            sidewalk_width_m=MainConfig.sidewalk_width_m,
        )

        result = await asyncio.get_event_loop().run_in_executor(None, generate_clusters, config)

        return str(result)

    def view_current_results(self):
        """View results of the latest completed step."""
        latest_step = max(
            [step for step, status in self.step_status.items() if status == "completed"], default=0
        )
        if latest_step > 0:
            self.view_step_results(latest_step)

    def view_step_results(self, step: int):
        """View results for a specific step."""
        try:
            if step not in self.results:
                self.log_viewer.add_log(f"⚠️ No results for step {step}", "warning")
                return

            # Render the step results
            geojson_path = self.results.get(step)
            image_path = self.visualizer.render_step_results(
                step, self.config.geopackage_path, geojson_path
            )

            self.log_viewer.add_log(f"🖼️ Rendered step {step} visualization", "info")

            # Try to display in terminal
            if self.visualizer.display_image_in_terminal(image_path):
                self.log_viewer.add_log(f"📺 Opened step {step} visualization", "success")
            else:
                # Show modal with path info
                self.push_screen(ImageViewerScreen(image_path, f"Step {step} Results"))
                self.log_viewer.add_log(f"💾 Image saved: {image_path}", "info")

        except Exception as e:
            self.log_viewer.add_log(f"❌ Failed to view step {step}: {str(e)}", "error")

    def browse_layers(self):
        """Browse layers in the geopackage."""
        self._update_layer_browser()
        # Switch to results tab
        tabbed_content = self.query_one(TabbedContent)
        tabbed_content.active = "results"

    def _update_layer_browser(self):
        """Update the layer browser and tree view with current data."""
        try:
            if Path(self.config.geopackage_path).exists():
                layer_info = self.visualizer.get_layer_info(self.config.geopackage_path)

                # Update both old layer browser and new tree view
                self.layer_browser.update_layers(layer_info)
                self.layer_tree.update_layers(self.config.geopackage_path, layer_info)

                self.log_viewer.add_log(f"📊 Found {len(layer_info)} layers", "info")
            else:
                self.layer_browser.update_layers([])
                self.layer_tree.update_layers("", [])

        except Exception as e:
            self.log_viewer.add_log(f"❌ Failed to browse layers: {str(e)}", "error")

    def _update_setup_visibility(self):
        """Show/hide setup content based on current step."""
        try:
            setup_content = self.query_one("#setup-content")

            # Show setup content only if Step 0 is not completed
            if self.step_status.get(0) == "completed":
                setup_content.add_class("hidden")
            else:
                setup_content.remove_class("hidden")

        except Exception:
            # Setup content might not be mounted yet
            pass

    def open_output_folder(self):
        """Open the output folder in file manager."""
        import platform
        import subprocess

        output_path = Path(self.config.output_dir)
        if not output_path.exists():
            self.log_viewer.add_log("⚠️ Output directory doesn't exist yet", "warning")
            return

        try:
            system = platform.system()
            if system == "Darwin":  # macOS
                subprocess.run(["open", str(output_path)])
            elif system == "Windows":
                subprocess.run(["explorer", str(output_path)])
            else:  # Linux
                subprocess.run(["xdg-open", str(output_path)])

            self.log_viewer.add_log(f"📁 Opened {output_path}", "success")

        except Exception as e:
            self.log_viewer.add_log(f"❌ Failed to open folder: {str(e)}", "error")

    def clear_results(self):
        """Clear all results and reset the workflow."""
        try:
            # Reset state
            self.step_status = {0: "pending", 1: "pending", 2: "pending", 3: "pending"}
            self.results.clear()

            # Reset UI
            for _step, card in self.step_cards.items():
                card.update_status("pending")

            # Reset Step 0 button
            self.query_one("#btn-step0").disabled = False
            self.query_one("#btn-step0").variant = "primary"
            self.query_one("#btn-step0").label = "▶️ Step 0"

            # Disable buttons (except step 0)
            for step in [1, 2, 3]:
                self.query_one(f"#btn-step{step}").disable()
                if step < 3:  # Only steps 1 and 2 have view buttons
                    self.query_one(f"#view-step{step}").disabled = True

            self.query_one("#btn-view").disable()
            self.query_one("#btn-browse").disable()

            # Clear displays
            self.progress_display.update("Ready to start")
            self.layer_browser.update_layers([])
            self.layer_tree.update_layers("", [])
            self.layer_detail_panel._show_welcome()

            # Clear output directory
            output_path = Path(self.config.output_dir)
            if output_path.exists():
                shutil.rmtree(output_path)
                output_path.mkdir(parents=True, exist_ok=True)

            # Update setup visibility
            self._update_setup_visibility()

            self.log_viewer.add_log("🧹 Results cleared, ready to start over", "info")

        except Exception as e:
            self.log_viewer.add_log(f"❌ Failed to clear results: {str(e)}", "error")


def main():
    """Main entry point for the TUI application."""
    app = RueTuiApp()
    app.run()


if __name__ == "__main__":
    main()
