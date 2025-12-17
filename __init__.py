# -*- coding: utf-8 -*-
"""
GeoRefiner Manual - QGIS Plugin

Interactive georeferencing refinement with side-by-side synchronized views.
Provides manual vertex adjustment with topology-preserving drag and
iterative refinement using Polynomial or TPS warping.
"""


def classFactory(iface):
    """
    Load the GeoRefiner Manual plugin class.

    This function is called by QGIS when the plugin is loaded.

    Args:
        iface: QgisInterface instance providing access to QGIS application.

    Returns:
        GeoRefinerManualPlugin instance.
    """
    from .georefiner_manual_plugin import GeoRefinerManualPlugin
    return GeoRefinerManualPlugin(iface)
