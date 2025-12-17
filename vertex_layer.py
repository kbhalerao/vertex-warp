# -*- coding: utf-8 -*-
"""
GeoRefiner Manual - Vertex Layer Management

Creates and manages vertex point layers for:
- Original positions (non-editable, hollow markers)
- Current positions (editable, solid markers)

Tracks vertex movements and provides visual feedback.
"""

from typing import List, Dict, Optional, Callable

from qgis.PyQt.QtCore import QVariant
from qgis.PyQt.QtGui import QColor
from qgis.core import (
    QgsVectorLayer,
    QgsFeature,
    QgsGeometry,
    QgsPointXY,
    QgsField,
    QgsFields,
    QgsCoordinateReferenceSystem,
    QgsMarkerSymbol,
    QgsSingleSymbolRenderer,
    QgsCategorizedSymbolRenderer,
    QgsRendererCategory,
)


TARGET_CRS = QgsCoordinateReferenceSystem("EPSG:3857")


class VertexLayerManager:
    """
    Manages vertex point layers for original and current positions.

    Creates two memory layers:
    - Original markers: hollow circles showing where vertices should be
    - Current markers: solid circles that can be dragged

    Tracks movements and provides callbacks for vertex drag events.
    """

    def __init__(self, vertices: List[Dict], on_vertex_moved: Optional[Callable] = None):
        """
        Initialize vertex layer manager.

        Args:
            vertices: List of vertex dictionaries from clipper.extract_vertices_from_boundary().
            on_vertex_moved: Callback function(vertex_id, new_x, new_y) called when vertex moves.
        """
        self.vertices = vertices
        self.on_vertex_moved = on_vertex_moved

        # Create layers
        self.original_layer = self._create_original_layer()
        self.current_layer = self._create_current_layer()

        # Apply styling
        self._style_original_layer()
        self._style_current_layer()

    def _create_original_layer(self) -> QgsVectorLayer:
        """Create the original positions layer (non-editable markers)."""
        layer = QgsVectorLayer(
            f"Point?crs={TARGET_CRS.authid()}",
            "Original Vertices",
            "memory"
        )

        # Add fields
        provider = layer.dataProvider()
        fields = QgsFields()
        fields.append(QgsField("vertex_id", QVariant.Int))
        fields.append(QgsField("x", QVariant.Double))
        fields.append(QgsField("y", QVariant.Double))
        provider.addAttributes(fields)
        layer.updateFields()

        # Add features
        features = []
        for v in self.vertices:
            feat = QgsFeature()
            feat.setGeometry(QgsGeometry.fromPointXY(
                QgsPointXY(v["orig_x"], v["orig_y"])
            ))
            feat.setAttributes([v["id"], v["orig_x"], v["orig_y"]])
            features.append(feat)

        provider.addFeatures(features)
        layer.updateExtents()

        return layer

    def _create_current_layer(self) -> QgsVectorLayer:
        """Create the current positions layer (editable markers)."""
        layer = QgsVectorLayer(
            f"Point?crs={TARGET_CRS.authid()}",
            "Current Vertices",
            "memory"
        )

        # Add fields
        provider = layer.dataProvider()
        fields = QgsFields()
        fields.append(QgsField("vertex_id", QVariant.Int))
        fields.append(QgsField("orig_x", QVariant.Double))
        fields.append(QgsField("orig_y", QVariant.Double))
        fields.append(QgsField("current_x", QVariant.Double))
        fields.append(QgsField("current_y", QVariant.Double))
        fields.append(QgsField("moved", QVariant.Int))  # 0 or 1
        provider.addAttributes(fields)
        layer.updateFields()

        # Add features
        features = []
        for v in self.vertices:
            feat = QgsFeature()
            feat.setGeometry(QgsGeometry.fromPointXY(
                QgsPointXY(v["current_x"], v["current_y"])
            ))
            feat.setAttributes([
                v["id"],
                v["orig_x"],
                v["orig_y"],
                v["current_x"],
                v["current_y"],
                1 if v["moved"] else 0,
            ])
            features.append(feat)

        provider.addFeatures(features)
        layer.updateExtents()

        return layer

    def _style_original_layer(self):
        """Apply hollow circle styling to original layer."""
        # Hollow cyan circles (target positions)
        symbol = QgsMarkerSymbol.createSimple({
            'name': 'circle',
            'color': '0,0,0,0',  # Transparent fill
            'outline_color': '#00FFFF',  # Cyan
            'outline_width': '1.5',
            'size': '8',
        })
        self.original_layer.setRenderer(QgsSingleSymbolRenderer(symbol))

    def _style_current_layer(self):
        """Apply categorized styling to current layer (moved vs unmoved)."""
        # Unmoved: hollow blue circles (matches original)
        symbol_unmoved = QgsMarkerSymbol.createSimple({
            'name': 'circle',
            'color': '0,0,0,0',  # Transparent fill
            'outline_color': '#2196F3',  # Blue
            'outline_width': '2',
            'size': '10',
        })

        # Moved: hollow royal blue circles (visible on warm NDVI)
        symbol_moved = QgsMarkerSymbol.createSimple({
            'name': 'circle',
            'color': '0,0,0,0',  # Transparent fill
            'outline_color': '#4169E1',  # Royal Blue
            'outline_width': '2',
            'size': '10',
        })

        # Create categorized renderer
        categories = [
            QgsRendererCategory(0, symbol_unmoved, "Not moved"),
            QgsRendererCategory(1, symbol_moved, "Moved"),
        ]
        renderer = QgsCategorizedSymbolRenderer("moved", categories)
        self.current_layer.setRenderer(renderer)

    def get_layers(self) -> tuple:
        """
        Get both vertex layers.

        Returns:
            Tuple of (original_layer, current_layer).
        """
        return self.original_layer, self.current_layer

    def update_vertex_position(self, vertex_id: int, new_x: float, new_y: float):
        """
        Update a vertex's current position (with callback).

        Args:
            vertex_id: ID of the vertex to update.
            new_x: New X coordinate in EPSG:3857.
            new_y: New Y coordinate in EPSG:3857.
        """
        self._update_layer_geometry(vertex_id, new_x, new_y)

        # Callback (don't use this if already handling updates externally)
        if self.on_vertex_moved:
            self.on_vertex_moved(vertex_id, new_x, new_y)

    def _update_layer_geometry(self, vertex_id: int, new_x: float, new_y: float):
        """
        Update layer geometry without triggering callback.
        Use this to avoid recursion when called from the callback handler.

        Args:
            vertex_id: ID of the vertex to update.
            new_x: New X coordinate in EPSG:3857.
            new_y: New Y coordinate in EPSG:3857.
        """
        # Update internal vertices list
        for v in self.vertices:
            if v["id"] == vertex_id:
                v["current_x"] = new_x
                v["current_y"] = new_y
                v["moved"] = True
                break

        # Update the current layer
        self.current_layer.startEditing()
        for feature in self.current_layer.getFeatures():
            if feature["vertex_id"] == vertex_id:
                # Update geometry
                self.current_layer.changeGeometry(
                    feature.id(),
                    QgsGeometry.fromPointXY(QgsPointXY(new_x, new_y))
                )
                # Update attributes
                self.current_layer.changeAttributeValue(
                    feature.id(),
                    self.current_layer.fields().indexOf("current_x"),
                    new_x
                )
                self.current_layer.changeAttributeValue(
                    feature.id(),
                    self.current_layer.fields().indexOf("current_y"),
                    new_y
                )
                self.current_layer.changeAttributeValue(
                    feature.id(),
                    self.current_layer.fields().indexOf("moved"),
                    1
                )
                break
        self.current_layer.commitChanges()

    def reset_vertex(self, vertex_id: int):
        """
        Reset a vertex to its original position.

        Args:
            vertex_id: ID of the vertex to reset.
        """
        for v in self.vertices:
            if v["id"] == vertex_id:
                self.update_vertex_position(vertex_id, v["orig_x"], v["orig_y"])
                v["moved"] = False
                # Update moved flag in layer
                self.current_layer.startEditing()
                for feature in self.current_layer.getFeatures():
                    if feature["vertex_id"] == vertex_id:
                        self.current_layer.changeAttributeValue(
                            feature.id(),
                            self.current_layer.fields().indexOf("moved"),
                            0
                        )
                        break
                self.current_layer.commitChanges()
                break

    def reset_all_vertices(self):
        """Reset all vertices to their original positions."""
        for v in self.vertices:
            v["current_x"] = v["orig_x"]
            v["current_y"] = v["orig_y"]
            v["moved"] = False

        # Recreate current layer
        self.current_layer = self._create_current_layer()
        self._style_current_layer()

    def get_gcps(self) -> List[Dict]:
        """
        Get GCPs from vertex movements.

        Returns:
            List of GCP dictionaries:
            [
                {
                    "vertex_id": 0,
                    "orig_x": ..., "orig_y": ...,  # Target position (where it should end up)
                    "current_x": ..., "current_y": ...,  # Source position (where it is now)
                    "moved": True/False,
                },
                ...
            ]

        Note: For warping, we want to warp FROM current position TO original position,
        so the target raster pixels at current_x/y get mapped to orig_x/y.
        """
        return [
            {
                "vertex_id": v["id"],
                "orig_x": v["orig_x"],
                "orig_y": v["orig_y"],
                "current_x": v["current_x"],
                "current_y": v["current_y"],
                "moved": v["moved"],
            }
            for v in self.vertices
        ]

    def get_moved_count(self) -> int:
        """Get count of moved vertices."""
        return sum(1 for v in self.vertices if v["moved"])

    def get_vertex_by_id(self, vertex_id: int) -> Optional[Dict]:
        """Get vertex by ID."""
        for v in self.vertices:
            if v["id"] == vertex_id:
                return v
        return None


def create_vertex_highlight_layer(x: float, y: float, crs: QgsCoordinateReferenceSystem) -> QgsVectorLayer:
    """
    Create a temporary layer to highlight the current vertex.

    Args:
        x: X coordinate.
        y: Y coordinate.
        crs: Coordinate reference system.

    Returns:
        QgsVectorLayer with a single highlighted point.
    """
    layer = QgsVectorLayer(
        f"Point?crs={crs.authid()}",
        "Current Vertex Highlight",
        "memory"
    )

    # Add a single feature
    provider = layer.dataProvider()
    feat = QgsFeature()
    feat.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
    provider.addFeatures([feat])
    layer.updateExtents()

    # Style: large red circle
    symbol = QgsMarkerSymbol.createSimple({
        'name': 'circle',
        'color': '255,0,0,100',  # Semi-transparent red
        'outline_color': '#FF0000',
        'outline_width': '2',
        'size': '16',
    })
    layer.setRenderer(QgsSingleSymbolRenderer(symbol))

    return layer
