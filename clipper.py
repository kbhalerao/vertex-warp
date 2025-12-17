# -*- coding: utf-8 -*-
"""
GeoRefiner Manual - Image Clipping Module

Clips target raster to field boundary with buffer and reprojects to EPSG:3857.
Uses GDAL for efficient processing.
"""

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Optional, Tuple

from qgis.core import (
    QgsVectorLayer,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsProject,
    QgsGeometry,
    QgsFeature,
)


TARGET_CRS = QgsCoordinateReferenceSystem("EPSG:3857")


def clip_raster_to_boundary(
    target_path: Path,
    boundary_layer: QgsVectorLayer,
    buffer_m: float,
    output_dir: Optional[Path] = None,
) -> Tuple[Path, Path]:
    """
    Clip a target raster to a buffered boundary polygon and reproject to EPSG:3857.

    Args:
        target_path: Path to the input raster.
        boundary_layer: QgsVectorLayer with the boundary polygon.
        buffer_m: Buffer distance in meters to apply around the boundary.
        output_dir: Directory for output files. If None, uses temp directory.

    Returns:
        Tuple of (clipped_raster_path, buffered_boundary_path).

    Raises:
        ValueError: If inputs are invalid.
        RuntimeError: If GDAL processing fails.
    """
    if not target_path.exists():
        raise ValueError(f"Target raster not found: {target_path}")

    if boundary_layer is None or not boundary_layer.isValid():
        raise ValueError("Invalid boundary layer")

    # Create output directory
    if output_dir is None:
        output_dir = Path(tempfile.mkdtemp(prefix="georefiner_"))
    output_dir.mkdir(parents=True, exist_ok=True)

    # Extract and transform boundary to EPSG:3857
    boundary_3857_path = output_dir / "boundary_3857.geojson"
    _export_boundary_to_3857(boundary_layer, boundary_3857_path)

    # Buffer the boundary
    buffered_path = output_dir / "boundary_buffered.geojson"
    _buffer_boundary(boundary_3857_path, buffered_path, buffer_m)

    # Clip and reproject target raster
    clipped_path = output_dir / f"{target_path.stem}_clipped.tif"
    _clip_raster(target_path, buffered_path, clipped_path)

    return clipped_path, buffered_path


def _export_boundary_to_3857(
    boundary_layer: QgsVectorLayer,
    output_path: Path,
) -> None:
    """
    Export boundary layer to GeoJSON in EPSG:3857.

    Args:
        boundary_layer: Source boundary layer.
        output_path: Path for output GeoJSON.
    """
    source_crs = boundary_layer.crs()

    # Create coordinate transform if needed
    transform = None
    if source_crs != TARGET_CRS:
        transform = QgsCoordinateTransform(
            source_crs, TARGET_CRS, QgsProject.instance()
        )

    # Extract all geometries and merge
    features = []
    for feature in boundary_layer.getFeatures():
        geom = QgsGeometry(feature.geometry())
        if transform:
            geom.transform(transform)
        features.append(geom)

    if not features:
        raise ValueError("No features found in boundary layer")

    # Merge geometries if multiple
    if len(features) == 1:
        merged_geom = features[0]
    else:
        merged_geom = QgsGeometry.unaryUnion(features)

    # Convert to GeoJSON
    geojson = {
        "type": "FeatureCollection",
        "crs": {
            "type": "name",
            "properties": {"name": "EPSG:3857"}
        },
        "features": [{
            "type": "Feature",
            "properties": {},
            "geometry": json.loads(merged_geom.asJson())
        }]
    }

    with open(output_path, "w") as f:
        json.dump(geojson, f)


def _buffer_boundary(
    input_path: Path,
    output_path: Path,
    buffer_m: float,
) -> None:
    """
    Buffer a GeoJSON polygon using ogr2ogr.

    Args:
        input_path: Input GeoJSON path.
        output_path: Output GeoJSON path.
        buffer_m: Buffer distance in meters (EPSG:3857 units).
    """
    # Use ogr2ogr with SQLite dialect for buffer operation
    cmd = [
        "ogr2ogr",
        "-f", "GeoJSON",
        str(output_path),
        str(input_path),
        "-dialect", "SQLite",
        "-sql", f"SELECT ST_Buffer(geometry, {buffer_m}) AS geometry FROM boundary_3857"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        # Fallback: try with QGIS buffer if ogr2ogr fails
        _buffer_boundary_qgis(input_path, output_path, buffer_m)


def _buffer_boundary_qgis(
    input_path: Path,
    output_path: Path,
    buffer_m: float,
) -> None:
    """
    Fallback buffer using QGIS geometry operations.

    Args:
        input_path: Input GeoJSON path.
        output_path: Output GeoJSON path.
        buffer_m: Buffer distance in meters.
    """
    # Load GeoJSON
    layer = QgsVectorLayer(str(input_path), "boundary", "ogr")
    if not layer.isValid():
        raise RuntimeError(f"Failed to load boundary: {input_path}")

    # Buffer the geometry
    for feature in layer.getFeatures():
        geom = feature.geometry()
        buffered = geom.buffer(buffer_m, 8)  # 8 segments per quarter circle

        # Write buffered GeoJSON
        geojson = {
            "type": "FeatureCollection",
            "crs": {
                "type": "name",
                "properties": {"name": "EPSG:3857"}
            },
            "features": [{
                "type": "Feature",
                "properties": {},
                "geometry": json.loads(buffered.asJson())
            }]
        }

        with open(output_path, "w") as f:
            json.dump(geojson, f)
        break  # Only need first feature


def _clip_raster(
    input_path: Path,
    cutline_path: Path,
    output_path: Path,
) -> None:
    """
    Clip and reproject raster using gdalwarp.

    Args:
        input_path: Input raster path.
        cutline_path: GeoJSON cutline polygon path.
        output_path: Output raster path.
    """
    cmd = [
        "gdalwarp",
        "-cutline", str(cutline_path),
        "-crop_to_cutline",
        "-t_srs", "EPSG:3857",
        "-r", "bilinear",
        "-co", "COMPRESS=LZW",
        "-co", "TILED=YES",
        "-co", "BIGTIFF=IF_SAFER",
        str(input_path),
        str(output_path),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            f"gdalwarp failed: {result.stderr}\nCommand: {' '.join(cmd)}"
        )


def extract_vertices_from_boundary(
    boundary_layer: QgsVectorLayer,
) -> list:
    """
    Extract vertices from boundary polygon and transform to EPSG:3857.

    Args:
        boundary_layer: Source boundary layer.

    Returns:
        List of dicts with vertex info:
        [
            {
                "id": 0,
                "orig_x": 1234567.89,
                "orig_y": 4567890.12,
                "current_x": 1234567.89,
                "current_y": 4567890.12,
                "moved": False
            },
            ...
        ]
    """
    source_crs = boundary_layer.crs()

    # Create coordinate transform if needed
    transform = None
    if source_crs != TARGET_CRS:
        transform = QgsCoordinateTransform(
            source_crs, TARGET_CRS, QgsProject.instance()
        )

    vertices = []
    vertex_id = 0

    for feature in boundary_layer.getFeatures():
        geom = feature.geometry()

        # Handle multi-part geometries
        if geom.isMultipart():
            parts = geom.asMultiPolygon()
        else:
            parts = [geom.asPolygon()]

        for polygon in parts:
            if not polygon:
                continue

            # Get exterior ring (first ring)
            exterior = polygon[0]

            # Skip the last point (it's a duplicate of the first for closed polygons)
            for i, point in enumerate(exterior[:-1]):
                x, y = point.x(), point.y()

                # Transform to EPSG:3857
                if transform:
                    transformed = transform.transform(x, y)
                    x, y = transformed.x(), transformed.y()

                vertices.append({
                    "id": vertex_id,
                    "orig_x": x,
                    "orig_y": y,
                    "current_x": x,
                    "current_y": y,
                    "moved": False,
                })
                vertex_id += 1

    return vertices
