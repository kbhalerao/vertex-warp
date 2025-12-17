# -*- coding: utf-8 -*-
"""
GeoRefiner Manual - Vertex Drag Map Tool

Custom QgsMapTool that allows dragging vertex points on the map canvas.
Supports topology-preserving drag when enabled.
"""

from typing import Optional, Callable, List, Dict

from qgis.PyQt.QtCore import Qt, pyqtSignal
from qgis.PyQt.QtGui import QCursor, QColor
from qgis.PyQt.QtWidgets import QApplication
from qgis.core import (
    QgsPointXY,
    QgsRectangle,
    QgsGeometry,
    QgsFeature,
    QgsVectorLayer,
    QgsWkbTypes,
)
from qgis.gui import QgsMapTool, QgsRubberBand, QgsMapCanvas


class VertexDragTool(QgsMapTool):
    """
    Map tool for dragging vertices interactively.

    Emits signals when vertices are moved, allowing the dialog
    to update the vertex layer and optionally apply topology-preserving drag.
    """

    # Signal emitted when a vertex is dragged
    # Arguments: vertex_id (int), new_x (float), new_y (float)
    vertexMoved = pyqtSignal(int, float, float)

    # Signal emitted when drag starts
    # Argument: vertex_id (int)
    dragStarted = pyqtSignal(int)

    # Signal emitted when drag ends
    # Argument: vertex_id (int)
    dragEnded = pyqtSignal(int)

    def __init__(
        self,
        canvas: QgsMapCanvas,
        vertex_layer: QgsVectorLayer,
        snap_tolerance: float = 15.0,
    ):
        """
        Initialize the vertex drag tool.

        Args:
            canvas: The QgsMapCanvas to attach to.
            vertex_layer: The vertex point layer to edit.
            snap_tolerance: Pixel tolerance for clicking on vertices.
        """
        super().__init__(canvas)
        self.canvas = canvas
        self.vertex_layer = vertex_layer
        self.snap_tolerance = snap_tolerance

        # Dragging state
        self.dragging = False
        self.dragged_feature_id: Optional[int] = None
        self.dragged_vertex_id: Optional[int] = None
        self.start_point: Optional[QgsPointXY] = None

        # Rubber band for visual feedback during drag
        self.rubber_band = QgsRubberBand(canvas, QgsWkbTypes.PointGeometry)
        self.rubber_band.setColor(QColor(255, 0, 0, 200))
        self.rubber_band.setWidth(3)
        self.rubber_band.setIcon(QgsRubberBand.ICON_CIRCLE)
        self.rubber_band.setIconSize(12)

        # Set cursor
        self.setCursor(Qt.OpenHandCursor)

    def canvasPressEvent(self, event):
        """Handle mouse press - start drag if on a vertex."""
        if event.button() != Qt.LeftButton:
            return

        # Get click point in map coordinates
        click_point = self.toMapCoordinates(event.pos())

        # Find nearest vertex within tolerance
        feature_id, vertex_id = self._find_nearest_vertex(click_point)

        if feature_id is not None:
            self.dragging = True
            self.dragged_feature_id = feature_id
            self.dragged_vertex_id = vertex_id
            self.start_point = click_point

            # Set up rubber band at current position
            self.rubber_band.reset(QgsWkbTypes.PointGeometry)
            self.rubber_band.addPoint(click_point)

            # Change cursor
            self.setCursor(Qt.ClosedHandCursor)

            # Emit signal
            self.dragStarted.emit(vertex_id)

    def canvasMoveEvent(self, event):
        """Handle mouse move - update rubber band during drag."""
        if not self.dragging:
            # Update cursor based on whether we're hovering over a vertex
            click_point = self.toMapCoordinates(event.pos())
            feature_id, _ = self._find_nearest_vertex(click_point)
            if feature_id is not None:
                self.setCursor(Qt.OpenHandCursor)
            else:
                self.setCursor(Qt.ArrowCursor)
            return

        # Update rubber band position
        current_point = self.toMapCoordinates(event.pos())
        self.rubber_band.reset(QgsWkbTypes.PointGeometry)
        self.rubber_band.addPoint(current_point)

    def canvasReleaseEvent(self, event):
        """Handle mouse release - complete drag."""
        if event.button() != Qt.LeftButton or not self.dragging:
            return

        # Get final position
        final_point = self.toMapCoordinates(event.pos())

        # Emit the move signal
        if self.dragged_vertex_id is not None:
            self.vertexMoved.emit(
                self.dragged_vertex_id,
                final_point.x(),
                final_point.y()
            )
            self.dragEnded.emit(self.dragged_vertex_id)

        # Reset state
        self.dragging = False
        self.dragged_feature_id = None
        self.dragged_vertex_id = None
        self.start_point = None

        # Clear rubber band
        self.rubber_band.reset(QgsWkbTypes.PointGeometry)

        # Reset cursor
        self.setCursor(Qt.OpenHandCursor)

    def _find_nearest_vertex(self, point: QgsPointXY) -> tuple:
        """
        Find the nearest vertex to a point within snap tolerance.

        Args:
            point: Map coordinates to search from.

        Returns:
            Tuple of (feature_id, vertex_id) or (None, None) if not found.
        """
        if self.vertex_layer is None or not self.vertex_layer.isValid():
            return None, None

        # Convert snap tolerance from pixels to map units
        map_tolerance = self.snap_tolerance * self.canvas.mapUnitsPerPixel()

        # Create search rectangle
        search_rect = QgsRectangle(
            point.x() - map_tolerance,
            point.y() - map_tolerance,
            point.x() + map_tolerance,
            point.y() + map_tolerance
        )

        # Search for features in rectangle
        best_feature_id = None
        best_vertex_id = None
        best_distance = float('inf')

        for feature in self.vertex_layer.getFeatures():
            geom = feature.geometry()
            if geom.isEmpty():
                continue

            feature_point = geom.asPoint()
            distance = point.distance(feature_point)

            if distance < map_tolerance and distance < best_distance:
                best_distance = distance
                best_feature_id = feature.id()
                best_vertex_id = feature["vertex_id"]

        return best_feature_id, best_vertex_id

    def deactivate(self):
        """Clean up when tool is deactivated."""
        self.rubber_band.reset(QgsWkbTypes.PointGeometry)
        super().deactivate()

    def update_vertex_layer(self, new_layer: QgsVectorLayer):
        """Update the vertex layer reference."""
        self.vertex_layer = new_layer
