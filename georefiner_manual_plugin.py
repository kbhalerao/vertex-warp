# -*- coding: utf-8 -*-
"""
VertexWarp - Main Plugin Class

This module contains the main plugin class that handles QGIS integration:
- Menu and toolbar creation
- Plugin lifecycle (initGui, unload)
- Dialog management

Copyright (C) 2025 Kaustubh Bhalerao, AgSci LLC
License: GPLv2
"""

from pathlib import Path

from qgis.PyQt.QtCore import QSettings
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction
from qgis.core import QgsApplication


class GeoRefinerManualPlugin:
    """
    Main QGIS Plugin class for VertexWarp.

    Handles plugin initialization, menu/toolbar setup, and dialog management.
    """

    def __init__(self, iface):
        """
        Initialize the plugin.

        Args:
            iface: QgisInterface instance providing access to QGIS application.
        """
        self.iface = iface
        self.plugin_dir = Path(__file__).parent
        self.actions = []
        self.menu_name = "&VertexWarp"
        self.toolbar = None
        self.dialog = None

        # Settings key prefix
        self.settings_prefix = "VertexWarp/"

    def initGui(self):
        """
        Create the menu entries, toolbar, and action.

        Called by QGIS when plugin is activated.
        """
        # Create toolbar
        self.toolbar = self.iface.addToolBar("VertexWarp")
        self.toolbar.setObjectName("VertexWarpToolbar")

        # Get icon path
        icon_path = self.plugin_dir / "resources" / "icon.png"
        if icon_path.exists():
            icon = QIcon(str(icon_path))
        else:
            # Use a default QGIS icon if custom icon not available
            icon = QgsApplication.getThemeIcon("/mActionShowPluginManager.svg")

        # Create main action
        self.action_open = QAction(
            icon,
            "VertexWarp",
            self.iface.mainWindow()
        )
        self.action_open.setObjectName("vertexwarp_open")
        self.action_open.setStatusTip("Open VertexWarp for interactive georeferencing refinement")
        self.action_open.triggered.connect(self.open_dialog)

        # Add to toolbar and menu
        self.toolbar.addAction(self.action_open)
        self.iface.addPluginToRasterMenu(self.menu_name, self.action_open)
        self.actions.append(self.action_open)

    def open_dialog(self):
        """Open the main VertexWarp dialog."""
        from .georefiner_manual_dialog import GeoRefinerManualDialog

        # Create dialog if needed
        if self.dialog is None:
            self.dialog = GeoRefinerManualDialog(self.iface, self.iface.mainWindow())

        # Show dialog (non-modal to allow interaction with main QGIS)
        self.dialog.show()
        self.dialog.raise_()
        self.dialog.activateWindow()

    def unload(self):
        """
        Remove the plugin menu items and icons from QGIS GUI.

        Called by QGIS when plugin is deactivated.
        """
        # Close dialog if open
        if self.dialog is not None:
            self.dialog.close()
            self.dialog.deleteLater()
            self.dialog = None

        # Remove actions from menu
        for action in self.actions:
            self.iface.removePluginRasterMenu(self.menu_name, action)

        # Remove toolbar
        if self.toolbar:
            del self.toolbar
            self.toolbar = None

        self.actions = []
