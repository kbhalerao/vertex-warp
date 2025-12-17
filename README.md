# VertexWarp - QGIS Plugin

Interactive georeferencing refinement for aerial imagery using manual vertex adjustment.

## Purpose

Aerial imagery from planes (2-5k ft AGL) often has approximate georeferencing with distortions that prevent sub-meter accuracy. This plugin provides a **manual, iterative refinement workflow** to achieve 0.5m or better accuracy needed for precision agriculture applications like zonal statistics on 3m x 10m grids.

## How It Works

1. You have a **target image** (aerial/drone imagery) that's approximately georeferenced
2. You have a **known-good boundary polygon** (e.g., field boundary traced from NAIP)
3. The plugin shows side-by-side views: reference (base map + boundary) | target (your image)
4. You drag vertices on the target to match where features actually appear
5. Click Refine - the image warps to align vertices with their target positions
6. Check alignment, refine again if needed, then export

## Installation

1. Copy the `georefiner_manual` folder to your QGIS plugins directory:
   - Linux: `~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/`
   - Windows: `%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\`
   - macOS: `~/Library/Application Support/QGIS/QGIS3/profiles/default/python/plugins/`

2. Restart QGIS

3. Enable the plugin: Plugins - Manage and Install Plugins - search "VertexWarp" - Enable

4. Access via: Plugins menu - VertexWarp, or the toolbar button

## Requirements

- QGIS 3.x
- No external dependencies (everything ships with QGIS)

## Usage

### Step 1: Prepare Your Data

- Load your base reference layer (XYZ tiles) in QGIS
- Have your field boundary polygon layer loaded
- Have your target raster file ready

### Step 2: Load Data

1. Click **Select target raster** and choose your aerial image
2. Select your **Boundary** layer from the dropdown
3. Set **Buffer** distance (default 100m) - how much area around the boundary to include
4. Click **Load**

The plugin will:
- Clip your target image to the boundary + buffer
- Reproject everything to EPSG:3857
- Display side-by-side views with vertex markers

### Step 3: Adjust Vertices

- **Left panel**: Reference view (base map + original boundary)
- **Right panel**: Target view (your image + draggable vertices)

**Cyan hollow circles** = Original/target positions (where vertices should end up)
**Blue hollow circles** = Current positions, not yet moved
**Royal blue hollow circles** = Current positions, moved

**Drag vertices** on the right panel to match where features actually appear in your image.

**Topology-preserving drag** (enabled by default): When you drag a vertex, nearby vertices move together proportionally. Adjust the **Falloff** slider (default 2m) to control the influence radius.

### Step 4: Preview & Refine

1. Select warp **Algorithm**:
   - Polynomial 1 (Affine) - simple, requires 3+ points
   - Polynomial 2 - more flexible, requires 6+ points
   - Polynomial 3 - most flexible, requires 10+ points
   - Thin-Plate Spline (TPS) - rubber-sheet, requires 3+ points

2. Click **Refine**

The image warps and vertices reset. Now check:
- Do image features align with the cyan markers?
- If yes - ready to export
- If no - drag vertices and refine again

### Step 5: Export

Click **Export** to save the refined image. You'll be asked if you want to load it into QGIS and close the dialog.

## Controls

### Buttons

| Button | Action |
|--------|--------|
| Load | Clip target and setup workspace |
| < Prev / Next > | Navigate between vertices |
| Start Over | Reset all vertices to original positions |
| Refine | Apply warp transformation (preview) |
| Export | Save final result |

### Keyboard Shortcuts

| Key | Action |
|-----|--------|
| Left / Right | Previous / Next vertex |
| Home / End | First / Last vertex |
| Ctrl+R | Refine |
| Ctrl+E | Export |

### Options

| Option | Description |
|--------|-------------|
| Topology-preserving drag | When enabled, nearby vertices move together |
| Falloff | Influence radius for topology drag (default 2m) |
| Algorithm | Warp transformation type |

## Tips

1. **Start with corners**: Corners are usually the most misaligned. Fix them first.

2. **Use the right algorithm**:
   - Few adjustments? Use Polynomial 1 or TPS
   - Many adjustments? Use Polynomial 2 or 3

3. **Iterate**: Don't try to get it perfect in one pass. Refine, check, refine again.

4. **Check alignment**: After refining, the cyan markers show where vertices should be. If your image features align with them, you're done.

5. **Falloff radius**: Keep it small (2-5m) for precise control. Increase if you want larger areas to move together.

## Technical Details

### Coordinate System

All processing is done in **EPSG:3857** (Web Mercator) for consistency with common base maps.

### Warp Process

1. Vertices define Ground Control Points (GCPs)
2. GCPs map current positions to original positions
3. GDAL `gdalwarp` applies the transformation
4. Output is LZW-compressed GeoTIFF

### Dependencies

**Zero external dependencies** - everything ships with QGIS:
- PyQt5 (Qt bindings)
- QGIS Python API
- GDAL/OGR (raster/vector processing)
- NumPy (numerical operations)

### File Outputs

Files are saved in the same directory as your target image:
- `{name}_clipped.tif` - Initial clipped image
- `{name}_clipped_refined_01.tif` - After first refinement
- `{name}_clipped_refined_02.tif` - After second refinement
- etc.

## Troubleshooting

**"No base layers found"**: Load an XYZ tile layer (like Google Satellite) in QGIS before opening the plugin.

**Vertices not dragging**: Make sure you're clicking directly on a vertex circle. The snap tolerance is 15 pixels.

**Warp looks wrong**: Try a different algorithm. TPS is more forgiving, polynomials need good GCP distribution.

## License

GPLv2

## Author

Kaustubh Bhalerao
AgSci LLC
kaustubh@agsci.com

https://github.com/kbhalerao/vertex-warp
