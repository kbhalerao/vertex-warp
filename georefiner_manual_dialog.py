# -*- coding: utf-8 -*-
"""
GeoRefiner Manual - Main Dialog

Side-by-side synchronized map views for interactive georeferencing refinement.

Layout:
+-----------------------------------------------------------------------+
| [Input Panel]                                                          |
| Target: [file picker]  Boundary: [layer combo]  Buffer: [spinbox] m   |
| [Load]                                                                 |
+-----------------------------------------------------------------------+
| [Left Canvas: Reference]        |  [Right Canvas: Target]              |
| Base XYZ + Original Boundary    |  Clipped Target + Editable Vertices  |
+-----------------------------------------------------------------------+
| [Vertex Navigation]  [< Prev] Vertex 3/15 - Moved: Yes [Next >]       |
| [Warp Options]       Algorithm: [dropdown]   [Refine]   [Export]      |
+-----------------------------------------------------------------------+
| [Status Bar]         Quality metrics, progress                         |
+-----------------------------------------------------------------------+
"""

from pathlib import Path
from typing import Optional, List

from qgis.PyQt.QtCore import Qt, QSettings
from qgis.PyQt.QtWidgets import (
    QApplication,
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QLabel,
    QPushButton,
    QSpinBox,
    QDoubleSpinBox,
    QGroupBox,
    QFileDialog,
    QMessageBox,
    QSplitter,
    QComboBox,
    QCheckBox,
    QSlider,
    QFrame,
    QSizePolicy,
    QStatusBar,
)
from qgis.core import (
    QgsProject,
    QgsMapLayerProxyModel,
    QgsRasterLayer,
    QgsVectorLayer,
    QgsCoordinateReferenceSystem,
)
from qgis.gui import QgsMapLayerComboBox, QgsMapCanvas, QgsMapToolPan, QgsMapToolZoom


class GeoRefinerManualDialog(QDialog):
    """
    Main dialog for GeoRefiner Manual plugin.

    Provides side-by-side synchronized map views for interactive
    georeferencing refinement.
    """

    SETTINGS_PREFIX = "GeoRefinerManual/"
    TARGET_CRS = QgsCoordinateReferenceSystem("EPSG:3857")

    def __init__(self, iface, parent=None):
        """
        Initialize the dialog.

        Args:
            iface: QgisInterface instance.
            parent: Parent widget.
        """
        super().__init__(parent)
        self.iface = iface
        self.settings = QSettings()

        # State
        self.target_path: Optional[Path] = None
        self.clipped_path: Optional[Path] = None
        self.boundary_layer: Optional[QgsVectorLayer] = None
        self.vertices: List[dict] = []  # Original vertex positions
        self.current_vertex_idx: int = 0
        self.vertex_manager = None  # VertexLayerManager instance
        self.highlight_layer = None  # Current vertex highlight

        # Layer references
        self.target_layer = None
        self.boundary_display = None
        self.boundary_target = None

        # Map tools
        self.left_pan_tool = None
        self.left_zoom_tool = None
        self.right_pan_tool = None
        self.right_zoom_tool = None
        self.vertex_drag_tool = None  # Custom tool for dragging vertices

        # Synchronization flag
        self._syncing = False

        # Setup UI
        self._setup_ui()
        self._setup_connections()
        self._load_settings()

    def _setup_ui(self):
        """Create the UI layout and widgets."""
        self.setWindowTitle("VertexWarp - Interactive Georeferencing Refinement")
        self.setMinimumSize(1200, 800)
        self.resize(1400, 900)

        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(8, 8, 8, 8)
        main_layout.setSpacing(8)

        # =================================================================
        # Input Panel
        # =================================================================
        input_group = QGroupBox("Input")
        input_layout = QGridLayout(input_group)
        input_layout.setSpacing(6)

        # Target raster file
        input_layout.addWidget(QLabel("Target Image:"), 0, 0)
        self.target_path_edit = QPushButton("Select target raster...")
        self.target_path_edit.clicked.connect(self._browse_target)
        input_layout.addWidget(self.target_path_edit, 0, 1)

        # Boundary layer selection
        input_layout.addWidget(QLabel("Boundary:"), 0, 2)
        self.boundary_combo = QgsMapLayerComboBox()
        self.boundary_combo.setFilters(QgsMapLayerProxyModel.PolygonLayer)
        self.boundary_combo.setAllowEmptyLayer(True)
        input_layout.addWidget(self.boundary_combo, 0, 3)

        # Buffer distance
        input_layout.addWidget(QLabel("Buffer:"), 0, 4)
        self.buffer_spin = QSpinBox()
        self.buffer_spin.setRange(10, 500)
        self.buffer_spin.setValue(100)
        self.buffer_spin.setSuffix(" m")
        input_layout.addWidget(self.buffer_spin, 0, 5)

        # Load button
        self.load_btn = QPushButton("Load")
        self.load_btn.setMinimumWidth(80)
        self.load_btn.clicked.connect(self._load_data)
        input_layout.addWidget(self.load_btn, 0, 6)

        main_layout.addWidget(input_group)

        # =================================================================
        # Dual Map Canvas (Side-by-Side)
        # =================================================================
        canvas_splitter = QSplitter(Qt.Horizontal)
        canvas_splitter.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        # Left canvas frame (Reference)
        left_frame = QFrame()
        left_frame.setFrameStyle(QFrame.StyledPanel)
        left_layout = QVBoxLayout(left_frame)
        left_layout.setContentsMargins(2, 2, 2, 2)
        left_layout.setSpacing(2)

        left_label = QLabel("Reference (Base Layer + Original Boundary)")
        left_label.setStyleSheet("font-weight: bold; color: #2196F3;")
        left_layout.addWidget(left_label)

        self.left_canvas = QgsMapCanvas()
        self.left_canvas.setCanvasColor(Qt.white)
        self.left_canvas.setDestinationCrs(self.TARGET_CRS)
        self.left_canvas.enableAntiAliasing(True)
        left_layout.addWidget(self.left_canvas)

        # Right canvas frame (Target)
        right_frame = QFrame()
        right_frame.setFrameStyle(QFrame.StyledPanel)
        right_layout = QVBoxLayout(right_frame)
        right_layout.setContentsMargins(2, 2, 2, 2)
        right_layout.setSpacing(2)

        right_label = QLabel("Target (Clipped Image + Editable Vertices)")
        right_label.setStyleSheet("font-weight: bold; color: #4CAF50;")
        right_layout.addWidget(right_label)

        self.right_canvas = QgsMapCanvas()
        self.right_canvas.setCanvasColor(Qt.white)
        self.right_canvas.setDestinationCrs(self.TARGET_CRS)
        self.right_canvas.enableAntiAliasing(True)
        right_layout.addWidget(self.right_canvas)

        canvas_splitter.addWidget(left_frame)
        canvas_splitter.addWidget(right_frame)
        canvas_splitter.setSizes([600, 600])

        main_layout.addWidget(canvas_splitter, stretch=1)

        # =================================================================
        # Control Panel
        # =================================================================
        control_group = QGroupBox("Controls")
        control_layout = QHBoxLayout(control_group)
        control_layout.setSpacing(12)

        # Vertex navigation
        nav_frame = QFrame()
        nav_layout = QHBoxLayout(nav_frame)
        nav_layout.setContentsMargins(0, 0, 0, 0)
        nav_layout.setSpacing(4)

        self.prev_btn = QPushButton("< Prev")
        self.prev_btn.setMinimumWidth(60)
        self.prev_btn.clicked.connect(self._prev_vertex)
        nav_layout.addWidget(self.prev_btn)

        self.vertex_label = QLabel("Vertex 0/0")
        self.vertex_label.setMinimumWidth(150)
        self.vertex_label.setAlignment(Qt.AlignCenter)
        nav_layout.addWidget(self.vertex_label)

        self.next_btn = QPushButton("Next >")
        self.next_btn.setMinimumWidth(60)
        self.next_btn.clicked.connect(self._next_vertex)
        nav_layout.addWidget(self.next_btn)

        self.start_over_btn = QPushButton("Start Over")
        self.start_over_btn.setMinimumWidth(80)
        self.start_over_btn.setStyleSheet("background-color: #FF9800; color: white;")
        self.start_over_btn.clicked.connect(self._start_over)
        nav_layout.addWidget(self.start_over_btn)

        control_layout.addWidget(nav_frame)

        # Separator
        sep1 = QFrame()
        sep1.setFrameShape(QFrame.VLine)
        control_layout.addWidget(sep1)

        # Topology-preserving drag
        topo_frame = QFrame()
        topo_layout = QHBoxLayout(topo_frame)
        topo_layout.setContentsMargins(0, 0, 0, 0)
        topo_layout.setSpacing(4)

        self.topo_check = QCheckBox("Topology-preserving drag")
        self.topo_check.setChecked(True)
        topo_layout.addWidget(self.topo_check)

        topo_layout.addWidget(QLabel("Falloff:"))
        self.falloff_slider = QSlider(Qt.Horizontal)
        self.falloff_slider.setRange(1, 50)  # 1-50m range, default 2m
        self.falloff_slider.setValue(2)
        self.falloff_slider.setMaximumWidth(100)
        topo_layout.addWidget(self.falloff_slider)

        self.falloff_label = QLabel("2m")
        self.falloff_label.setMinimumWidth(40)
        self.falloff_slider.valueChanged.connect(
            lambda v: self.falloff_label.setText(f"{v}m")
        )
        topo_layout.addWidget(self.falloff_label)

        control_layout.addWidget(topo_frame)

        # Separator
        sep2 = QFrame()
        sep2.setFrameShape(QFrame.VLine)
        control_layout.addWidget(sep2)

        # Warp algorithm
        warp_frame = QFrame()
        warp_layout = QHBoxLayout(warp_frame)
        warp_layout.setContentsMargins(0, 0, 0, 0)
        warp_layout.setSpacing(4)

        warp_layout.addWidget(QLabel("Algorithm:"))
        self.warp_combo = QComboBox()
        self.warp_combo.addItems([
            "Polynomial 1 (Affine)",
            "Polynomial 2",
            "Polynomial 3",
            "Thin-Plate Spline (TPS)",
        ])
        self.warp_combo.setCurrentIndex(0)
        warp_layout.addWidget(self.warp_combo)

        control_layout.addWidget(warp_frame)

        # Spacer
        control_layout.addStretch()

        # Action buttons
        self.refine_btn = QPushButton("Refine")
        self.refine_btn.setMinimumWidth(80)
        self.refine_btn.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold;")
        self.refine_btn.clicked.connect(self._refine)
        control_layout.addWidget(self.refine_btn)

        self.export_btn = QPushButton("Export")
        self.export_btn.setMinimumWidth(80)
        self.export_btn.clicked.connect(self._export)
        control_layout.addWidget(self.export_btn)

        main_layout.addWidget(control_group)

        # =================================================================
        # Status Bar
        # =================================================================
        self.status_bar = QStatusBar()
        self.status_bar.showMessage("Ready. Load a target image and boundary to begin.")
        main_layout.addWidget(self.status_bar)

        # Setup map tools
        self._setup_map_tools()

    def _setup_map_tools(self):
        """Setup pan and zoom tools for both canvases."""
        # Left canvas tools
        self.left_pan_tool = QgsMapToolPan(self.left_canvas)
        self.left_canvas.setMapTool(self.left_pan_tool)

        # Right canvas tools
        self.right_pan_tool = QgsMapToolPan(self.right_canvas)
        self.right_canvas.setMapTool(self.right_pan_tool)

    def _setup_connections(self):
        """Setup signal/slot connections."""
        # Synchronize canvas extents
        self.left_canvas.extentsChanged.connect(self._sync_left_to_right)
        self.right_canvas.extentsChanged.connect(self._sync_right_to_left)

    def _sync_left_to_right(self):
        """Sync right canvas extent to match left canvas."""
        if not self._syncing:
            self._syncing = True
            self.right_canvas.setExtent(self.left_canvas.extent())
            self.right_canvas.refresh()
            self._syncing = False

    def _sync_right_to_left(self):
        """Sync left canvas extent to match right canvas."""
        if not self._syncing:
            self._syncing = True
            self.left_canvas.setExtent(self.right_canvas.extent())
            self.left_canvas.refresh()
            self._syncing = False

    def _load_settings(self):
        """Load saved settings."""
        buffer = self.settings.value(
            self.SETTINGS_PREFIX + "buffer_m", 100, type=int
        )
        self.buffer_spin.setValue(buffer)

        warp_idx = self.settings.value(
            self.SETTINGS_PREFIX + "warp_algorithm", 0, type=int
        )
        self.warp_combo.setCurrentIndex(warp_idx)

        topo_enabled = self.settings.value(
            self.SETTINGS_PREFIX + "topo_drag", True, type=bool
        )
        self.topo_check.setChecked(topo_enabled)

        falloff = self.settings.value(
            self.SETTINGS_PREFIX + "falloff_m", 2, type=int
        )
        self.falloff_slider.setValue(falloff)

    def _save_settings(self):
        """Save current settings."""
        self.settings.setValue(
            self.SETTINGS_PREFIX + "buffer_m", self.buffer_spin.value()
        )
        self.settings.setValue(
            self.SETTINGS_PREFIX + "warp_algorithm", self.warp_combo.currentIndex()
        )
        self.settings.setValue(
            self.SETTINGS_PREFIX + "topo_drag", self.topo_check.isChecked()
        )
        self.settings.setValue(
            self.SETTINGS_PREFIX + "falloff_m", self.falloff_slider.value()
        )

    def _browse_target(self):
        """Open file dialog to select target raster."""
        last_dir = self.settings.value(
            self.SETTINGS_PREFIX + "last_target_dir", ""
        )
        filepath, _ = QFileDialog.getOpenFileName(
            self,
            "Select Target Raster",
            last_dir,
            "GeoTIFF (*.tif *.tiff);;All files (*.*)"
        )
        if filepath:
            self.target_path = Path(filepath)
            self.target_path_edit.setText(self.target_path.name)
            self.settings.setValue(
                self.SETTINGS_PREFIX + "last_target_dir",
                str(self.target_path.parent)
            )

    def _load_data(self):
        """Load target image and boundary, clip and display."""
        from .clipper import clip_raster_to_boundary, extract_vertices_from_boundary

        # Validate inputs
        if self.target_path is None:
            QMessageBox.warning(self, "Missing Input", "Please select a target raster.")
            return

        boundary_layer = self.boundary_combo.currentLayer()
        if boundary_layer is None:
            QMessageBox.warning(self, "Missing Input", "Please select a boundary layer.")
            return

        self.boundary_layer = boundary_layer
        buffer_m = self.buffer_spin.value()

        self.status_bar.showMessage(f"Step 1/4: Preparing boundary...")
        self.load_btn.setEnabled(False)
        QApplication.processEvents()  # Force UI update

        try:
            # Step 1: Extract vertices from boundary
            self.status_bar.showMessage(f"Step 1/4: Extracting {boundary_layer.featureCount()} boundary vertices...")
            QApplication.processEvents()
            self.vertices = extract_vertices_from_boundary(boundary_layer)
            self.current_vertex_idx = 0

            # Step 2: Clip the target raster
            self.status_bar.showMessage(f"Step 2/4: Clipping {self.target_path.name} (this may take a moment)...")
            QApplication.processEvents()
            self.clipped_path, buffered_boundary_path = clip_raster_to_boundary(
                self.target_path,
                boundary_layer,
                buffer_m,
            )

            # Step 3: Load layers into canvases
            self.status_bar.showMessage(f"Step 3/4: Loading layers into canvas...")
            QApplication.processEvents()
            self._setup_canvas_layers(buffered_boundary_path)

            # Step 4: Final setup
            self.status_bar.showMessage(f"Step 4/4: Setting up vertex editing...")
            QApplication.processEvents()
            self._update_vertex_label()

            # Zoom to first vertex
            self.current_vertex_idx = 0
            self._zoom_to_current_vertex()

            self.status_bar.showMessage(
                f"Ready! Loaded {len(self.vertices)} vertices. Drag vertices on right panel, then click Refine."
            )

        except Exception as e:
            QMessageBox.critical(
                self, "Error",
                f"Failed to load data:\n{str(e)}"
            )
            self.status_bar.showMessage("Error loading data.")

        finally:
            self.load_btn.setEnabled(True)

    def _setup_canvas_layers(self, buffered_boundary_path: Path):
        """Setup layers in both canvases after loading."""
        from .vertex_layer import VertexLayerManager, create_vertex_highlight_layer
        from qgis.core import QgsFillSymbol

        # Get all visible XYZ/WMS layers from main QGIS as reference base layers
        base_layers = []
        for layer in QgsProject.instance().mapLayers().values():
            if isinstance(layer, QgsRasterLayer):
                # Check if it's an XYZ or WMS tile layer
                provider = layer.providerType()
                if provider in ("wms", "xyz"):
                    base_layers.append(layer)

        # Load the original boundary into left canvas
        boundary_display = QgsVectorLayer(
            str(buffered_boundary_path), "Boundary (Reference)", "ogr"
        )
        if boundary_display.isValid():
            # Style the boundary
            symbol = QgsFillSymbol.createSimple({
                'color': '0,0,0,0',  # Transparent fill
                'outline_color': '#2196F3',  # Blue outline
                'outline_width': '2',
            })
            boundary_display.renderer().setSymbol(symbol)

        # Load the clipped target into right canvas
        target_layer = QgsRasterLayer(str(self.clipped_path), "Target (Clipped)")
        if not target_layer.isValid():
            raise RuntimeError(f"Failed to load clipped raster: {self.clipped_path}")

        # Create vertex layer manager
        self.vertex_manager = VertexLayerManager(
            self.vertices,
            on_vertex_moved=self._on_vertex_moved
        )
        original_vertices_layer, current_vertices_layer = self.vertex_manager.get_layers()

        # Setup left canvas (reference): base layers + boundary + original vertices
        left_layers = []
        if boundary_display.isValid():
            left_layers.append(boundary_display)
        left_layers.append(original_vertices_layer)
        left_layers.extend(base_layers)
        self.left_canvas.setLayers(left_layers)

        # Setup right canvas (target): clipped target + boundary + current vertices + original vertices
        boundary_target = QgsVectorLayer(
            str(buffered_boundary_path), "Boundary (Target)", "ogr"
        )
        if boundary_target.isValid():
            symbol = QgsFillSymbol.createSimple({
                'color': '0,0,0,0',
                'outline_color': '#4CAF50',  # Green outline
                'outline_width': '2',
            })
            boundary_target.renderer().setSymbol(symbol)

        # Include original vertices on right canvas too (so user can see target positions)
        right_layers = [current_vertices_layer, original_vertices_layer]
        if boundary_target.isValid():
            right_layers.append(boundary_target)
        right_layers.append(target_layer)
        self.right_canvas.setLayers(right_layers)

        # Store references
        self.target_layer = target_layer
        self.boundary_display = boundary_display
        self.boundary_target = boundary_target

        # Zoom to extent
        if target_layer.isValid():
            extent = target_layer.extent()
            # Sync both canvases to this extent
            self._syncing = True
            self.left_canvas.setExtent(extent)
            self.right_canvas.setExtent(extent)
            self.left_canvas.refresh()
            self.right_canvas.refresh()
            self._syncing = False

        # Setup vertex drag tool on right canvas
        self._setup_vertex_drag_tool(current_vertices_layer)

    def _setup_vertex_drag_tool(self, vertex_layer):
        """Setup the custom vertex drag tool on the right canvas."""
        from .vertex_drag_tool import VertexDragTool

        # Create drag tool
        self.vertex_drag_tool = VertexDragTool(
            self.right_canvas,
            vertex_layer,
            snap_tolerance=15.0
        )

        # Connect signals
        self.vertex_drag_tool.vertexMoved.connect(self._handle_vertex_drag)
        self.vertex_drag_tool.dragStarted.connect(self._on_drag_started)
        self.vertex_drag_tool.dragEnded.connect(self._on_drag_ended)

        # Set as active tool on right canvas
        self.right_canvas.setMapTool(self.vertex_drag_tool)

    def _handle_vertex_drag(self, vertex_id: int, new_x: float, new_y: float):
        """Handle vertex drag with optional topology preservation."""
        from .topo_drag import apply_topology_drag

        topo_enabled = self.topo_check.isChecked()
        sigma = float(self.falloff_slider.value())

        if topo_enabled and len(self.vertices) >= 3:
            # Apply topology-preserving drag
            updated = apply_topology_drag(
                self.vertices,
                vertex_id,
                new_x,
                new_y,
                sigma=sigma,
                enabled=True
            )

            # Update all affected vertices
            for v_update in updated:
                vid = v_update["id"]
                vx = v_update["x"]
                vy = v_update["y"]
                moved = v_update["moved"]

                # Update internal state
                for v in self.vertices:
                    if v["id"] == vid:
                        v["current_x"] = vx
                        v["current_y"] = vy
                        v["moved"] = moved
                        break

                # Update vertex manager (use _update_layer_geometry to avoid recursion)
                if self.vertex_manager and moved:
                    self.vertex_manager._update_layer_geometry(vid, vx, vy)
        else:
            # Simple drag - only move one vertex
            self._on_vertex_moved(vertex_id, new_x, new_y)

        # Refresh canvas
        self.right_canvas.refresh()
        self._update_vertex_label()

    def _on_drag_started(self, vertex_id: int):
        """Called when vertex drag starts."""
        self.status_bar.showMessage(f"Dragging vertex {vertex_id}...")

    def _on_drag_ended(self, vertex_id: int):
        """Called when vertex drag ends."""
        moved_count = sum(1 for v in self.vertices if v.get("moved", False))
        self.status_bar.showMessage(
            f"Vertex {vertex_id} moved. Total moved: {moved_count}/{len(self.vertices)}"
        )

    def _on_vertex_moved(self, vertex_id: int, new_x: float, new_y: float):
        """Callback when a vertex is moved (single vertex, no topology)."""
        # Update internal state
        for v in self.vertices:
            if v["id"] == vertex_id:
                v["current_x"] = new_x
                v["current_y"] = new_y
                v["moved"] = True
                break

        # Update vertex manager (updates the layer geometry)
        # Note: Don't use the callback version to avoid recursion
        if self.vertex_manager:
            self.vertex_manager._update_layer_geometry(vertex_id, new_x, new_y)

        # Refresh canvas
        self.right_canvas.refresh()

    def _update_vertex_label(self):
        """Update the vertex navigation label."""
        if not self.vertices:
            self.vertex_label.setText("No vertices")
            return

        vertex = self.vertices[self.current_vertex_idx]
        moved = vertex.get("moved", False)
        self.vertex_label.setText(
            f"Vertex {self.current_vertex_idx + 1}/{len(self.vertices)} - "
            f"Moved: {'Yes' if moved else 'No'}"
        )

    def _prev_vertex(self):
        """Navigate to previous vertex."""
        if self.vertices and self.current_vertex_idx > 0:
            self.current_vertex_idx -= 1
            self._zoom_to_current_vertex()
            self._update_vertex_label()

    def _next_vertex(self):
        """Navigate to next vertex."""
        if self.vertices and self.current_vertex_idx < len(self.vertices) - 1:
            self.current_vertex_idx += 1
            self._zoom_to_current_vertex()
            self._update_vertex_label()

    def _reset_vertices_to_original(self):
        """Reset all vertices to original positions (used after warp preview)."""
        if not self.vertices:
            return

        # Reset internal vertex state
        for v in self.vertices:
            v["current_x"] = v["orig_x"]
            v["current_y"] = v["orig_y"]
            v["moved"] = False

        # Reset vertex manager layers
        if self.vertex_manager:
            self.vertex_manager.reset_all_vertices()

            # Update right canvas with new current layer + original layer
            original_layer, current_layer = self.vertex_manager.get_layers()
            right_layers = [current_layer, original_layer]
            if self.boundary_target and self.boundary_target.isValid():
                right_layers.append(self.boundary_target)
            if self.target_layer and self.target_layer.isValid():
                right_layers.append(self.target_layer)
            self.right_canvas.setLayers(right_layers)

            # Update drag tool
            if self.vertex_drag_tool:
                self.vertex_drag_tool.update_vertex_layer(current_layer)

        self._update_vertex_label()
        self.right_canvas.refresh()

    def _start_over(self):
        """Reset all vertices to original positions and reload."""
        if not self.vertices:
            return

        reply = QMessageBox.question(
            self, "Start Over",
            "Reset all vertices to their original positions?\n\n"
            "This will undo all vertex movements.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply != QMessageBox.Yes:
            return

        self._reset_vertices_to_original()

        # Reset to first vertex and zoom
        self.current_vertex_idx = 0
        self._zoom_to_current_vertex()

        self.status_bar.showMessage("Reset all vertices to original positions. Start fresh!")

    def _zoom_to_current_vertex(self):
        """Zoom both canvases to current vertex with padding."""
        from qgis.core import QgsRectangle, QgsPointXY

        if not self.vertices:
            return

        vertex = self.vertices[self.current_vertex_idx]
        x = vertex["current_x"]
        y = vertex["current_y"]

        # Create extent centered on vertex with padding (50m on each side)
        padding = 50.0
        extent = QgsRectangle(
            x - padding, y - padding,
            x + padding, y + padding
        )

        # Sync both canvases to this extent
        self._syncing = True
        self.left_canvas.setExtent(extent)
        self.right_canvas.setExtent(extent)
        self.left_canvas.refresh()
        self.right_canvas.refresh()
        self._syncing = False

    def _refine(self):
        """Execute refinement warp."""
        from .warper import refine_image, get_transform_type, get_min_gcps

        if self.clipped_path is None or not self.vertices:
            QMessageBox.warning(self, "Not Ready", "Please load data first.")
            return

        algorithm_idx = self.warp_combo.currentIndex()
        transform_type = get_transform_type(algorithm_idx)
        min_gcps = get_min_gcps(transform_type)

        if len(self.vertices) < min_gcps:
            QMessageBox.warning(
                self, "Not Enough GCPs",
                f"{transform_type} requires at least {min_gcps} GCPs.\n"
                f"You have {len(self.vertices)} vertices."
            )
            return

        self.status_bar.showMessage(f"Refining with {transform_type}...")
        self.refine_btn.setEnabled(False)

        try:
            # Get output directory (same as clipped file)
            output_dir = self.clipped_path.parent

            # Get vertices from manager
            vertices = self.vertex_manager.get_gcps() if self.vertex_manager else self.vertices

            # Perform refinement
            result = refine_image(
                self.clipped_path,
                vertices,
                output_dir,
                algorithm_idx,
            )

            if result.success:
                # Load the refined result into the right canvas
                self._load_refined_result(result.output_path)

                # Update clipped_path to the refined version for next iteration
                self.clipped_path = result.output_path

                # Reset vertices to original positions (image is now warped to align)
                # This lets user see if the warp worked - features should align with blue markers
                self._reset_vertices_to_original()

                self.status_bar.showMessage(
                    f"Preview: GCPs: {result.gcp_count}, "
                    f"Mean shift: {result.mean_shift_m:.2f}m. "
                    f"Check alignment - continue refining or Export when satisfied."
                )

                QMessageBox.information(
                    self, "Preview Ready",
                    f"Warp applied - check alignment!\n\n"
                    f"Transform: {result.transform_type}\n"
                    f"GCPs used: {result.gcp_count}\n"
                    f"Mean shift: {result.mean_shift_m:.2f}m\n\n"
                    "Blue markers show target positions.\n"
                    "If image features align with blue markers, the warp worked!\n\n"
                    "Continue refining if needed, or Export when satisfied."
                )
            else:
                QMessageBox.critical(
                    self, "Refinement Failed",
                    f"Failed to refine image:\n{result.error_message}"
                )
                self.status_bar.showMessage("Refinement failed.")

        except Exception as e:
            QMessageBox.critical(
                self, "Error",
                f"Error during refinement:\n{str(e)}"
            )
            self.status_bar.showMessage("Error during refinement.")

        finally:
            self.refine_btn.setEnabled(True)

    def _load_refined_result(self, refined_path: Path):
        """Load refined result into right canvas, replacing previous target."""
        # Create new raster layer
        refined_layer = QgsRasterLayer(str(refined_path), "Refined Target")
        if not refined_layer.isValid():
            raise RuntimeError(f"Failed to load refined raster: {refined_path}")

        # Get vertex layers
        original_layer = None
        current_layer = None
        if self.vertex_manager:
            original_layer, current_layer = self.vertex_manager.get_layers()

        # Update right canvas layers (include original vertices for alignment check)
        right_layers = []
        if current_layer:
            right_layers.append(current_layer)
        if original_layer:
            right_layers.append(original_layer)
        if self.boundary_target and self.boundary_target.isValid():
            right_layers.append(self.boundary_target)
        right_layers.append(refined_layer)
        self.right_canvas.setLayers(right_layers)

        # Update reference
        self.target_layer = refined_layer

        # Refresh canvas
        self.right_canvas.refresh()

    def _export(self):
        """Export the refined result, optionally load into QGIS and close."""
        if self.clipped_path is None:
            QMessageBox.warning(self, "Nothing to Export", "Please refine an image first.")
            return

        # Get export path
        default_name = f"{self.clipped_path.stem}_final.tif"
        filepath, _ = QFileDialog.getSaveFileName(
            self,
            "Export Refined Image",
            str(self.clipped_path.parent / default_name),
            "GeoTIFF (*.tif *.tiff);;All files (*.*)"
        )

        if not filepath:
            return

        try:
            import shutil
            shutil.copy2(self.clipped_path, filepath)
            self.status_bar.showMessage(f"Exported to: {filepath}")

            # Ask if user wants to load into QGIS and close
            reply = QMessageBox.question(
                self, "Export Complete",
                f"Successfully exported to:\n{Path(filepath).name}\n\n"
                "Load into QGIS and close this dialog?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.Yes
            )

            if reply == QMessageBox.Yes:
                # Load into main QGIS window
                layer_name = Path(filepath).stem
                layer = QgsRasterLayer(filepath, layer_name)
                if layer.isValid():
                    QgsProject.instance().addMapLayer(layer)
                    self.iface.messageBar().pushSuccess(
                        "VertexWarp",
                        f"Loaded refined image: {layer_name}"
                    )
                # Close dialog
                self.accept()

        except Exception as e:
            QMessageBox.critical(
                self, "Export Failed",
                f"Failed to export:\n{str(e)}"
            )

    def keyPressEvent(self, event):
        """Handle keyboard shortcuts."""
        from qgis.PyQt.QtCore import Qt

        key = event.key()

        # Left/Right arrows for vertex navigation
        if key == Qt.Key_Left:
            self._prev_vertex()
            event.accept()
        elif key == Qt.Key_Right:
            self._next_vertex()
            event.accept()
        # Home/End for first/last vertex
        elif key == Qt.Key_Home:
            if self.vertices:
                self.current_vertex_idx = 0
                self._zoom_to_current_vertex()
                self._update_vertex_label()
            event.accept()
        elif key == Qt.Key_End:
            if self.vertices:
                self.current_vertex_idx = len(self.vertices) - 1
                self._zoom_to_current_vertex()
                self._update_vertex_label()
            event.accept()
        # R for refine
        elif key == Qt.Key_R and event.modifiers() == Qt.ControlModifier:
            self._refine()
            event.accept()
        # E for export
        elif key == Qt.Key_E and event.modifiers() == Qt.ControlModifier:
            self._export()
            event.accept()
        else:
            super().keyPressEvent(event)

    def closeEvent(self, event):
        """Handle dialog close."""
        self._save_settings()
        event.accept()
