# -*- coding: utf-8 -*-
"""
GeoRefiner Manual - Warping Module

Handles GCP extraction and GDAL-based image warping with multiple transform types:
- Polynomial 1 (affine) - requires 3+ GCPs
- Polynomial 2 - requires 6+ GCPs
- Polynomial 3 - requires 10+ GCPs
- Thin-Plate Spline (TPS) - requires 3+ GCPs
"""

import subprocess
import tempfile
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

from osgeo import gdal


@dataclass
class WarpResult:
    """Result of a warp operation."""
    success: bool
    output_path: Optional[Path]
    error_message: Optional[str] = None
    gcp_count: int = 0
    transform_type: str = ""
    mean_shift_m: float = 0.0


def get_transform_type(algorithm_idx: int) -> str:
    """
    Get GDAL transform type string from algorithm index.

    Args:
        algorithm_idx: Index from UI dropdown:
            0 = Polynomial 1 (Affine)
            1 = Polynomial 2
            2 = Polynomial 3
            3 = Thin-Plate Spline (TPS)

    Returns:
        Transform type string for gdalwarp.
    """
    transforms = {
        0: "polynomial1",
        1: "polynomial2",
        2: "polynomial3",
        3: "tps",
    }
    return transforms.get(algorithm_idx, "polynomial1")


def get_min_gcps(transform_type: str) -> int:
    """
    Get minimum number of GCPs required for a transform type.

    Args:
        transform_type: GDAL transform type string.

    Returns:
        Minimum GCP count.
    """
    minimums = {
        "polynomial1": 3,
        "polynomial2": 6,
        "polynomial3": 10,
        "tps": 3,
    }
    return minimums.get(transform_type, 3)


def extract_gcps_for_warp(
    vertices: List[Dict],
    raster_path: Path,
) -> Tuple[List[Tuple], float]:
    """
    Extract GCPs from vertex movements for GDAL warping.

    For refinement, we want to warp the target image so that:
    - Pixels at "current" positions move to "original" positions

    GDAL GCPs specify: (pixel_x, pixel_y, geo_x, geo_y)
    where pixel coords are SOURCE and geo coords are DESTINATION.

    Args:
        vertices: List of vertex dictionaries with orig_x/y and current_x/y.
        raster_path: Path to the raster to get pixel coordinates.

    Returns:
        Tuple of (list of GCP tuples, mean shift in meters).
    """
    # Open raster to get geotransform
    ds = gdal.Open(str(raster_path))
    if ds is None:
        raise RuntimeError(f"Failed to open raster: {raster_path}")

    gt = ds.GetGeoTransform()
    ds = None

    # Inverse geotransform for geo -> pixel conversion
    # gt[0] = x origin, gt[1] = pixel width, gt[2] = rotation
    # gt[3] = y origin, gt[4] = rotation, gt[5] = pixel height (negative)
    def geo_to_pixel(geo_x, geo_y):
        """Convert geographic coordinates to pixel coordinates."""
        px = (geo_x - gt[0]) / gt[1]
        py = (geo_y - gt[3]) / gt[5]
        return px, py

    gcps = []
    total_shift = 0.0
    moved_count = 0

    for v in vertices:
        # Only use vertices that have been moved (or use all for initial warp)
        # For refinement, we use ALL vertices as GCPs
        current_x = v["current_x"]
        current_y = v["current_y"]
        orig_x = v["orig_x"]
        orig_y = v["orig_y"]

        # Convert current position (where vertex appears in target) to pixel coords
        px, py = geo_to_pixel(current_x, current_y)

        # GCP: pixel coords in source -> geo coords in destination
        # We want pixels at current_x/y to move to orig_x/y
        gcp = (px, py, orig_x, orig_y)
        gcps.append(gcp)

        # Calculate shift
        shift = ((current_x - orig_x) ** 2 + (current_y - orig_y) ** 2) ** 0.5
        total_shift += shift
        if v.get("moved", False):
            moved_count += 1

    mean_shift = total_shift / len(vertices) if vertices else 0.0

    return gcps, mean_shift


def warp_with_gcps(
    input_path: Path,
    gcps: List[Tuple],
    output_path: Path,
    transform_type: str = "polynomial1",
    resampling: str = "bilinear",
) -> WarpResult:
    """
    Warp a raster using GCPs.

    Args:
        input_path: Input raster path.
        gcps: List of GCP tuples (pixel_x, pixel_y, geo_x, geo_y).
        output_path: Output raster path.
        transform_type: GDAL transform type (polynomial1, polynomial2, polynomial3, tps).
        resampling: Resampling method (bilinear, nearest, cubic, lanczos).

    Returns:
        WarpResult with success status and details.
    """
    min_gcps = get_min_gcps(transform_type)
    if len(gcps) < min_gcps:
        return WarpResult(
            success=False,
            output_path=None,
            error_message=f"Not enough GCPs. {transform_type} requires at least {min_gcps}, got {len(gcps)}.",
            gcp_count=len(gcps),
            transform_type=transform_type,
        )

    try:
        # Step 1: Create a temporary file with GCPs embedded
        temp_with_gcps = output_path.parent / f"{output_path.stem}_gcps.tif"

        # Build gdal_translate command to add GCPs
        cmd_translate = [
            "gdal_translate",
            "-of", "GTiff",
        ]

        # Add each GCP
        for gcp in gcps:
            px, py, geo_x, geo_y = gcp
            cmd_translate.extend(["-gcp", str(px), str(py), str(geo_x), str(geo_y)])

        cmd_translate.extend([str(input_path), str(temp_with_gcps)])

        result = subprocess.run(cmd_translate, capture_output=True, text=True)
        if result.returncode != 0:
            return WarpResult(
                success=False,
                output_path=None,
                error_message=f"gdal_translate failed: {result.stderr}",
                gcp_count=len(gcps),
                transform_type=transform_type,
            )

        # Step 2: Warp using the GCPs
        cmd_warp = [
            "gdalwarp",
            "-r", resampling,
            "-order" if transform_type.startswith("polynomial") else "-tps",
        ]

        # Add polynomial order if applicable
        if transform_type == "polynomial1":
            cmd_warp.extend(["1"])
        elif transform_type == "polynomial2":
            cmd_warp.extend(["2"])
        elif transform_type == "polynomial3":
            cmd_warp.extend(["3"])
        # For TPS, the -tps flag was already added, no additional param needed
        elif transform_type == "tps":
            # Remove the empty element added by the conditional
            cmd_warp = [
                "gdalwarp",
                "-r", resampling,
                "-tps",
            ]

        cmd_warp.extend([
            "-co", "COMPRESS=LZW",
            "-co", "TILED=YES",
            "-co", "BIGTIFF=IF_SAFER",
            str(temp_with_gcps),
            str(output_path),
        ])

        result = subprocess.run(cmd_warp, capture_output=True, text=True)

        # Cleanup temp file
        if temp_with_gcps.exists():
            temp_with_gcps.unlink()

        if result.returncode != 0:
            return WarpResult(
                success=False,
                output_path=None,
                error_message=f"gdalwarp failed: {result.stderr}",
                gcp_count=len(gcps),
                transform_type=transform_type,
            )

        return WarpResult(
            success=True,
            output_path=output_path,
            gcp_count=len(gcps),
            transform_type=transform_type,
        )

    except Exception as e:
        return WarpResult(
            success=False,
            output_path=None,
            error_message=str(e),
            gcp_count=len(gcps),
            transform_type=transform_type,
        )


def refine_image(
    input_path: Path,
    vertices: List[Dict],
    output_dir: Path,
    algorithm_idx: int = 0,
) -> WarpResult:
    """
    Perform a refinement warp on the input image.

    Args:
        input_path: Path to the clipped target raster.
        vertices: List of vertex dictionaries from VertexLayerManager.
        output_dir: Directory for output files.
        algorithm_idx: Warp algorithm index from UI dropdown.

    Returns:
        WarpResult with refined image path and statistics.
    """
    transform_type = get_transform_type(algorithm_idx)

    # Extract GCPs
    gcps, mean_shift = extract_gcps_for_warp(vertices, input_path)

    # Generate output filename
    iteration = 1
    output_path = output_dir / f"{input_path.stem}_refined_{iteration:02d}.tif"
    while output_path.exists():
        iteration += 1
        output_path = output_dir / f"{input_path.stem}_refined_{iteration:02d}.tif"

    # Perform warp
    result = warp_with_gcps(input_path, gcps, output_path, transform_type)
    result.mean_shift_m = mean_shift

    return result
