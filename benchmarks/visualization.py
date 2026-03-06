"""
Unified Benchmark Visualization using Napari

Interactive napari-based visualization for benchmarks and raw microscopy tiles.
Streams data from remote server and saves screenshots back.

Visualization modes:
- segmentation: Compare Cellpose, Cellpose4, StarDist segmentation methods
- spots: Compare Standard vs Spotiflow spot detection
- phenotype: Generate 4 grayscale channel images (DAPI, Tubulin, yH2AX, Actin) - zoomed
- sbs: Generate GTAC composite overlay image - zoomed
- illum: Illumination correction function + example DAPI - full field with corner markers

Usage:
    # Benchmark visualizations
    python visualization.py --mode segmentation --sample P-1_W-B3_T-169
    python visualization.py --mode spots --sample P-1_W-B1_T-2

    # Raw tile visualizations (random sample if not specified)
    python visualization.py --mode phenotype
    python visualization.py --mode sbs --sample P-1_W-A1_T-0
    python visualization.py --mode illum

    # List available samples
    python visualization.py --list-samples

    # Clear local cache
    python visualization.py --clear-cache

Remote Configuration:
    Benchmarks: fry.wi.mit.edu:/lab/ops_analysis/cheeseman/aconcagua-analysis-new/benchmarks/results
    Phenotype:  fry.wi.mit.edu:/lab/ops_analysis/cheeseman/aconcagua-analysis-new/analysis/brieflow_output/phenotype/images
    SBS:        fry.wi.mit.edu:/lab/ops_analysis/cheeseman/aconcagua-analysis-new/analysis/brieflow_output/sbs/images
    Illum:      fry.wi.mit.edu:/lab/ops_analysis/cheeseman/aconcagua-analysis-new/analysis/brieflow_output/phenotype/illum
    Screenshots saved to: remote:results/visualizations/screenshots/
"""

import argparse
import tempfile
import subprocess
from pathlib import Path
import numpy as np

# ============================================================================
# CONFIGURATION
# ============================================================================

# Remote server settings
REMOTE_HOST = "fry.wi.mit.edu"
REMOTE_BASE = (
    "/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis-new/benchmarks/results"
)

# Local cache for downloaded files
LOCAL_CACHE = Path(tempfile.gettempdir()) / "napari_viz_cache"
LOCAL_CACHE.mkdir(exist_ok=True)

# Visualization parameters
SEGMENTATION_CONFIG = {
    "window_size": 800,
    "center_position": [1800, 1700],  # [y, x]
    "view_width": 1000,
    "nuclei_opacity": 0.6,
    "cell_opacity": 0.6,
    "contour_thickness": 2,
    "pixels_per_micron": 1 / 0.3035714285714286,
}

SPOT_CONFIG = {
    "window_size": 1200,
    "center_position": [800, 550],  # [y, x]
    "view_width": 200,
    "spot_size": 6,
    "spot_edge_width": 0.5,
    "standard_color": "red",
    "spotiflow_color": "yellow",
    "standard_threshold": 400,
    "pixels_per_micron": 1 / 1.214286,
    "cycle_index": 0,
    "contrast_limits": {
        "G": [4000, 30000],
        "T": [4000, 30000],
        "A": [12000, 30000],
        "C": [12000, 30000],
    },
}

PHENOTYPE_CONFIG = {
    "window_size": 800,
    "view_width": 800,  # Zoomed view (less zoom = larger view_width)
    "margin": 250,  # Margin from edges for random center
    "pixels_per_micron": 1 / 0.3035714285714286,
    "channel_names": ["DAPI", "Tubulin", "yH2AX", "Actin"],
    "contrast_limits": {
        "DAPI": [254.0, 10000.0],
        "Tubulin": [600, 1500.0],
        "yH2AX": [0.0, 4000.0],
        "Actin": [350, 1500.0],
    },
}

SBS_CONFIG = {
    "window_size": 800,
    "view_width": 400,  # Zoomed view
    "margin": 250,  # Margin from edges for random center
    "pixels_per_micron": 1 / 1.214286,
    "cycle_index": 0,
    "dapi_opacity": 0.5,
    "dapi_contrast": 0.10,
    "channel_names": ["G", "T", "A", "C"],
    "channel_colors": ["green", "red", "magenta", "cyan"],
    "contrast_limits": {
        "G": [4000, 18000],
        "T": [4000, 18000],
        "A": [12000, 20000],
        "C": [12000, 20000],
    },
}

ILLUM_CONFIG = {
    "window_size": 800,
    "pixels_per_micron": 1 / 0.3035714285714286,
}

# Remote paths for raw pipeline output
REMOTE_PHENOTYPE_DIR = "/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis-new/analysis/brieflow_output/phenotype/images"
REMOTE_SBS_DIR = "/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis-new/analysis/brieflow_output/sbs/images"
REMOTE_IC_FIELDS_DIR = "/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis-new/analysis/brieflow_output/preprocess/ic_fields/phenotype"
REMOTE_PREPROCESS_DIR = "/archive/cheeseman/ops_analysis/aconcagua-analysis/analysis/brieflow_output/preprocess/images/phenotype"

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================


def get_random_center(height, width, margin):
    """Get a random center position within image bounds.

    Args:
        height: Image height
        width: Image width
        margin: Minimum distance from edges

    Returns:
        [y, x] center position
    """
    import random

    y = random.randint(margin, height - margin)
    x = random.randint(margin, width - margin)
    return [y, x]


# ============================================================================
# REMOTE DATA UTILITIES
# ============================================================================


def run_ssh_command(cmd):
    """Run a command on the remote server and return output."""
    result = subprocess.run(
        ["ssh", REMOTE_HOST, cmd],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"SSH command failed: {result.stderr}")
    return result.stdout.strip()


def download_file(remote_path, local_path=None):
    """Download a file from remote server using scp."""
    if local_path is None:
        local_path = LOCAL_CACHE / Path(remote_path).name

    local_path = Path(local_path)
    local_path.parent.mkdir(parents=True, exist_ok=True)

    # Skip if already cached
    if local_path.exists():
        print(f"  Using cached: {local_path.name}")
        return local_path

    print(f"  Downloading: {Path(remote_path).name}...")
    result = subprocess.run(
        ["scp", f"{REMOTE_HOST}:{remote_path}", str(local_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Download failed: {result.stderr}")

    return local_path


def upload_file(local_path, remote_path):
    """Upload a file to remote server using scp."""
    print(f"  Uploading: {Path(local_path).name}...")

    # Ensure remote directory exists
    remote_dir = str(Path(remote_path).parent)
    run_ssh_command(f"mkdir -p {remote_dir}")

    result = subprocess.run(
        ["scp", str(local_path), f"{REMOTE_HOST}:{remote_path}"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Upload failed: {result.stderr}")


def list_remote_files(remote_dir, pattern="*"):
    """List files matching pattern in remote directory."""
    try:
        output = run_ssh_command(f"ls {remote_dir}/{pattern} 2>/dev/null")
        return [Path(f).name for f in output.split("\n") if f]
    except RuntimeError:
        return []


def list_segmentation_samples():
    """List available segmentation samples on remote."""
    remote_dir = f"{REMOTE_BASE}/segmentation/aligned"
    files = list_remote_files(remote_dir, "*__aligned.tiff")
    return [f.replace("__aligned.tiff", "") for f in files]


def list_spot_samples():
    """List available spot calling samples on remote."""
    remote_dir = f"{REMOTE_BASE}/spot_calling/aligned"
    files = list_remote_files(remote_dir, "*__aligned.tiff")
    return [f.replace("__aligned.tiff", "") for f in files]


# ============================================================================
# SEGMENTATION VISUALIZATION
# ============================================================================


def load_segmentation_data(sample_id):
    """Download and load segmentation data for a sample."""
    print(f"\nLoading segmentation data for {sample_id}...")

    # Download aligned image
    aligned_remote = f"{REMOTE_BASE}/segmentation/aligned/{sample_id}__aligned.tiff"
    aligned_local = download_file(aligned_remote)

    # Download segmentation results for each method
    methods = ["cellpose", "cellpose4", "stardist"]
    data = {"aligned": aligned_local, "nuclei": {}, "cells": {}}

    for method in methods:
        nuclei_remote = f"{REMOTE_BASE}/segmentation/{method}/{sample_id}_nuclei.tiff"
        cells_remote = f"{REMOTE_BASE}/segmentation/{method}/{sample_id}_cells.tiff"

        # Use unique local names to avoid cache collision between methods
        nuclei_local = LOCAL_CACHE / f"{sample_id}_{method}_nuclei.tiff"
        cells_local = LOCAL_CACHE / f"{sample_id}_{method}_cells.tiff"

        try:
            data["nuclei"][method] = download_file(nuclei_remote, nuclei_local)
            data["cells"][method] = download_file(cells_remote, cells_local)
        except RuntimeError:
            print(f"  Warning: {method} data not found")
            data["nuclei"][method] = None
            data["cells"][method] = None

    return data


def run_segmentation_visualization(sample_id, output_dir=None):
    """Run interactive segmentation visualization."""
    import napari
    import skimage.io
    from qtpy.QtWidgets import QApplication
    from qtpy.QtCore import QTimer
    import time

    # Load data
    data = load_segmentation_data(sample_id)

    # Read images
    print("\nReading image files...")
    image = skimage.io.imread(data["aligned"])
    print(f"  Image shape: {image.shape}")

    nuclei = {}
    cells = {}
    for method in ["cellpose", "cellpose4", "stardist"]:
        if data["nuclei"][method]:
            nuclei[method] = skimage.io.imread(data["nuclei"][method])
            cells[method] = skimage.io.imread(data["cells"][method])
            print(
                f"  Loaded {method}: nuclei={nuclei[method].shape}, cells={cells[method].shape}"
            )

    # Create viewer
    print("\nStarting napari viewer...")
    viewer = napari.Viewer()

    # Extract channels (assuming shape is H x W x C or similar)
    if image.ndim == 3:
        if image.shape[2] <= 5:  # Channels last
            dapi = image[:, :, 0]
            phalloidin = image[:, :, 3] if image.shape[2] > 3 else image[:, :, -1]
        else:  # Channels first
            dapi = image[0, :, :]
            phalloidin = image[3, :, :] if image.shape[0] > 3 else image[-1, :, :]
    else:
        dapi = image
        phalloidin = image

    # Add DAPI channel
    viewer.add_image(
        dapi,
        name="DAPI",
        colormap="bop blue",
        contrast_limits=[dapi.min(), dapi.max() * 0.3],
    )

    # Add Phalloidin channel
    viewer.add_image(
        phalloidin,
        name="Phalloidin",
        colormap="gray",
        contrast_limits=[phalloidin.min(), phalloidin.max() * 0.45],
        opacity=0.6,
    )

    # Add segmentation layers (initially hidden)
    nuclei_layers = {}
    cell_layers = {}
    cfg = SEGMENTATION_CONFIG

    for method in ["cellpose", "cellpose4", "stardist"]:
        if method in nuclei:
            nuclei_layers[method] = viewer.add_labels(
                nuclei[method],
                name=f"Nuclei_{method}",
                opacity=cfg["nuclei_opacity"],
                visible=False,
            )
            cell_layers[method] = viewer.add_labels(
                cells[method],
                name=f"Cells_{method}",
                opacity=cfg["cell_opacity"],
                visible=False,
            )
            cell_layers[method].contour = cfg["contour_thickness"]

    # Configure viewer
    viewer.scale_bar.unit = "μm"
    viewer.scale_bar.position = "bottom_right"
    viewer.scale_bar.ticks = False
    viewer.window.resize(cfg["window_size"], cfg["window_size"])

    zoom_level = cfg["window_size"] / cfg["view_width"]
    viewer.camera.center = cfg["center_position"]
    viewer.camera.zoom = zoom_level

    # Setup output directory
    if output_dir is None:
        output_dir = LOCAL_CACHE / "screenshots" / "segmentation"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    def take_screenshots():
        """Capture screenshots for all segmentation methods."""
        print("\nCapturing screenshots...")

        # Hide UI panels
        viewer.window.qt_viewer.dockLayerList.setVisible(False)
        viewer.window.qt_viewer.dockLayerControls.setVisible(False)
        canvas = viewer.window._qt_viewer.canvas.native

        scale_bar_modes = [
            ("with_text", True, 20),
            ("no_text", True, 0),
            ("no_scalebar", False, 0),
        ]

        for mode_name, visible, font_size in scale_bar_modes:
            viewer.scale_bar.visible = visible
            if visible:
                viewer.scale_bar.font_size = font_size

            # Base image (no segmentation)
            path = output_dir / f"base_image_{mode_name}.png"
            QApplication.processEvents()
            time.sleep(0.3)
            canvas.grab().save(str(path))
            print(f"  Saved: {path.name}")

            # Each segmentation method
            for method in nuclei_layers:
                nuclei_layers[method].visible = True
                cell_layers[method].visible = True

                QApplication.processEvents()
                time.sleep(0.3)

                path = output_dir / f"segmentation_{method}_{mode_name}.png"
                canvas.grab().save(str(path))
                print(f"  Saved: {path.name}")

                nuclei_layers[method].visible = False
                cell_layers[method].visible = False

        print(f"\nScreenshots saved to: {output_dir}")

        # Upload to remote
        print("\nUploading screenshots to remote...")
        remote_output = f"{REMOTE_BASE}/visualizations/screenshots/segmentation"
        for png_file in output_dir.glob("*.png"):
            upload_file(png_file, f"{remote_output}/{png_file.name}")
        print("Upload complete!")

        # Close the viewer
        viewer.close()

    # Auto-screenshot after viewer initializes
    timer = QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(take_screenshots)
    timer.start(2000)

    print("Screenshots will be captured automatically in 2 seconds...")
    print("You can also manually adjust the view before screenshots are taken.")
    napari.run()


# ============================================================================
# SPOT CALLING VISUALIZATION
# ============================================================================


def load_spot_data(sample_id):
    """Download and load spot calling data for a sample."""
    print(f"\nLoading spot calling data for {sample_id}...")

    data = {
        "aligned": None,
        "max_filtered": None,
        "standard_peaks": None,
        "spotiflow_peaks": None,
    }

    # Download aligned image from spot_calling/aligned/
    aligned_remote = f"{REMOTE_BASE}/spot_calling/aligned/{sample_id}__aligned.tiff"
    try:
        data["aligned"] = download_file(aligned_remote)
    except RuntimeError:
        print("  Warning: Aligned image not found")

    # Download max_filtered image (optional, for better visualization)
    max_filtered_remote = (
        f"{REMOTE_BASE}/spot_calling/aligned/{sample_id}__max_filtered.tiff"
    )
    try:
        data["max_filtered"] = download_file(max_filtered_remote)
    except RuntimeError:
        print("  Warning: Max filtered image not found (will use aligned)")

    # Download peaks TIFF files (use unique local names to avoid cache collision)
    standard_peaks_remote = (
        f"{REMOTE_BASE}/spot_calling/standard/{sample_id}_peaks.tiff"
    )
    try:
        standard_local = LOCAL_CACHE / f"{sample_id}_standard_peaks.tiff"
        data["standard_peaks"] = download_file(standard_peaks_remote, standard_local)
    except RuntimeError:
        print("  Warning: Standard peaks not found")

    spotiflow_peaks_remote = (
        f"{REMOTE_BASE}/spot_calling/spotiflow/{sample_id}_peaks.tiff"
    )
    try:
        spotiflow_local = LOCAL_CACHE / f"{sample_id}_spotiflow_peaks.tiff"
        data["spotiflow_peaks"] = download_file(spotiflow_peaks_remote, spotiflow_local)
    except RuntimeError:
        print("  Warning: Spotiflow peaks not found")

    return data


def binary_to_points(binary_mask, threshold=0, is_standard=False):
    """Convert binary peak mask to (y, x) point coordinates.

    Args:
        binary_mask: 2D or 3D array of peak locations
        threshold: Only include peaks > threshold (for standard method)
        is_standard: If True, apply intensity threshold filter

    Returns:
        Array of [y, x] coordinates
    """
    if binary_mask.ndim == 3:
        # Handle multi-channel mask (height x width x channels) or (channels x height x width)
        # Determine axis order
        if binary_mask.shape[-1] <= 5:
            # (H, W, C) format
            n_channels = binary_mask.shape[-1]
            points = []
            for c in range(n_channels):
                if is_standard:
                    channel_points = np.where(binary_mask[:, :, c] > threshold)
                else:
                    channel_points = np.where(binary_mask[:, :, c] > 0)
                for y, x in zip(channel_points[0], channel_points[1]):
                    points.append([y, x])
            return np.array(points) if points else np.empty((0, 2))
        else:
            # (C, H, W) format
            n_channels = binary_mask.shape[0]
            points = []
            for c in range(n_channels):
                if is_standard:
                    channel_points = np.where(binary_mask[c, :, :] > threshold)
                else:
                    channel_points = np.where(binary_mask[c, :, :] > 0)
                for y, x in zip(channel_points[0], channel_points[1]):
                    points.append([y, x])
            return np.array(points) if points else np.empty((0, 2))
    else:
        # 2D mask
        if is_standard:
            coords = np.where(binary_mask > threshold)
        else:
            coords = np.where(binary_mask > 0)
        return (
            np.array(list(zip(coords[0], coords[1])))
            if len(coords[0]) > 0
            else np.empty((0, 2))
        )


def run_spot_visualization(sample_id, output_dir=None):
    """Run interactive spot calling visualization."""
    import napari
    import tifffile
    from qtpy.QtWidgets import QApplication
    from qtpy.QtCore import QTimer
    import time

    # Load data
    data = load_spot_data(sample_id)

    if not data.get("aligned"):
        print("Error: No aligned image found for this sample")
        return

    if not data.get("standard_peaks") and not data.get("spotiflow_peaks"):
        print("Error: No peaks data found for this sample")
        return

    # Read images
    print("\nReading image files...")
    cfg = SPOT_CONFIG

    aligned_image = tifffile.imread(data["aligned"])
    print(f"  Aligned image shape: {aligned_image.shape}")

    # Load max_filtered if available (better for spot visualization)
    if data.get("max_filtered"):
        max_image = tifffile.imread(data["max_filtered"])
        print(f"  Max filtered image shape: {max_image.shape}")
    else:
        max_image = None

    # Convert peaks TIFF to point coordinates
    standard_points = np.empty((0, 2))
    spotiflow_points = np.empty((0, 2))

    if data.get("standard_peaks"):
        standard_peaks = tifffile.imread(data["standard_peaks"])
        standard_points = binary_to_points(
            standard_peaks, threshold=cfg["standard_threshold"], is_standard=True
        )
        print(
            f"  Standard spots (threshold>{cfg['standard_threshold']}): {len(standard_points)}"
        )

    if data.get("spotiflow_peaks"):
        spotiflow_peaks = tifffile.imread(data["spotiflow_peaks"])
        spotiflow_points = binary_to_points(
            spotiflow_peaks, threshold=0, is_standard=False
        )
        print(f"  Spotiflow spots: {len(spotiflow_points)}")

    # Create viewer
    print("\nStarting napari viewer...")
    viewer = napari.Viewer()

    # Extract channels from aligned image
    # Shape is typically (cycles, channels, H, W) or (channels, H, W)
    cycle_idx = cfg["cycle_index"]

    if aligned_image.ndim == 4:
        # (cycles, channels, H, W)
        dapi = aligned_image[cycle_idx, 0, :, :]
        channel_data = aligned_image[cycle_idx, 1:5, :, :]  # G, T, A, C
    elif aligned_image.ndim == 3:
        if aligned_image.shape[0] <= 5:
            # (channels, H, W)
            dapi = aligned_image[0, :, :]
            channel_data = aligned_image[1:5, :, :]
        else:
            # (H, W, channels)
            dapi = aligned_image[:, :, 0]
            channel_data = np.moveaxis(aligned_image[:, :, 1:5], -1, 0)
    else:
        dapi = aligned_image
        channel_data = None

    # Add DAPI channel
    viewer.add_image(
        dapi,
        name="DAPI",
        colormap="gray",
        blending="additive",
        contrast_limits=[dapi.min(), dapi.max() * 0.15],
        opacity=0.5,
    )

    # Add GTAC channels
    if channel_data is not None:
        channel_names = ["G", "T", "A", "C"]
        cmaps = ["green", "red", "magenta", "cyan"]

        for i, (name, cmap) in enumerate(zip(channel_names, cmaps)):
            if i < channel_data.shape[0]:
                limits = cfg["contrast_limits"].get(name, [0, 30000])
                viewer.add_image(
                    channel_data[i],
                    name=name,
                    colormap=cmap,
                    blending="additive",
                    contrast_limits=limits,
                )

    # Add spot layers
    standard_layer = None
    spotiflow_layer = None

    if len(standard_points) > 0:
        standard_layer = viewer.add_points(
            standard_points,
            name="Standard Spots",
            size=cfg["spot_size"],
            face_color=[0, 0, 0, 0],  # Transparent face (hollow circles)
            symbol="o",
            border_color=cfg["standard_color"],
            border_width=cfg["spot_edge_width"],
            border_width_is_relative=False,
            visible=False,
        )

    if len(spotiflow_points) > 0:
        spotiflow_layer = viewer.add_points(
            spotiflow_points,
            name="Spotiflow Spots",
            size=cfg["spot_size"],
            face_color=[0, 0, 0, 0],  # Transparent face (hollow circles)
            symbol="o",
            border_color=cfg["spotiflow_color"],
            border_width=cfg["spot_edge_width"],
            border_width_is_relative=False,
            visible=False,
        )

    # Configure viewer
    viewer.scale_bar.unit = "μm"
    viewer.scale_bar.position = "bottom_right"
    viewer.scale_bar.ticks = False
    viewer.window.resize(cfg["window_size"], cfg["window_size"])

    zoom_level = cfg["window_size"] / cfg["view_width"]
    viewer.camera.center = cfg["center_position"]
    viewer.camera.zoom = zoom_level

    # Setup output directory
    if output_dir is None:
        output_dir = LOCAL_CACHE / "screenshots" / "spots"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    def take_screenshots():
        """Capture screenshots for spot visualization."""
        print("\nCapturing screenshots...")

        viewer.window.qt_viewer.dockLayerList.setVisible(False)
        viewer.window.qt_viewer.dockLayerControls.setVisible(False)
        canvas = viewer.window._qt_viewer.canvas.native

        spot_modes = [
            ("no_spots", False, False),
            ("standard_spots", True, False),
            ("spotiflow_spots", False, True),
            ("both_spots", True, True),
        ]

        scale_bar_modes = [
            ("with_text", True, 20),
            ("no_text", True, 0),
            ("no_scalebar", False, 0),
        ]

        for scale_name, visible, font_size in scale_bar_modes:
            viewer.scale_bar.visible = visible
            if visible:
                viewer.scale_bar.font_size = font_size

            for spot_name, show_standard, show_spotiflow in spot_modes:
                if standard_layer:
                    standard_layer.visible = show_standard
                if spotiflow_layer:
                    spotiflow_layer.visible = show_spotiflow

                QApplication.processEvents()
                time.sleep(0.3)

                path = output_dir / f"{spot_name}_{scale_name}.png"
                canvas.grab().save(str(path))
                print(f"  Saved: {path.name}")

        print(f"\nScreenshots saved to: {output_dir}")

        # Upload to remote
        print("\nUploading screenshots to remote...")
        remote_output = f"{REMOTE_BASE}/visualizations/screenshots/spots"
        for png_file in output_dir.glob("*.png"):
            upload_file(png_file, f"{remote_output}/{png_file.name}")
        print("Upload complete!")

        # Close the viewer
        viewer.close()

    timer = QTimer()
    timer.setSingleShot(True)
    timer.timeout.connect(take_screenshots)
    timer.start(2000)

    print("Screenshots will be captured automatically in 2 seconds...")
    napari.run()


# ============================================================================
# PHENOTYPE TILE VISUALIZATION
# ============================================================================


# Default plate/well to use for phenotype visualization (for consistent contrast)
PHENOTYPE_DEFAULT_PLATE_WELL = "P-1_W-B3"


def list_phenotype_samples():
    """List available phenotype samples on remote (random selection from fixed plate/well)."""
    try:
        # Filter to specific plate/well for consistent contrast limits
        output = run_ssh_command(
            f"ls {REMOTE_PHENOTYPE_DIR}/ | grep '{PHENOTYPE_DEFAULT_PLATE_WELL}' | grep '__aligned.tiff' | shuf | head -20"
        )
        files = [f for f in output.split("\n") if f]
        return [f.replace("__aligned.tiff", "") for f in files]
    except RuntimeError:
        return []


def load_phenotype_data(sample_id):
    """Download phenotype aligned image."""
    print(f"\nLoading phenotype data for {sample_id}...")

    aligned_remote = f"{REMOTE_PHENOTYPE_DIR}/{sample_id}__aligned.tiff"
    try:
        aligned_local = download_file(aligned_remote)
        return {"aligned": aligned_local}
    except RuntimeError as e:
        print(f"  Error: {e}")
        return {"aligned": None}


def run_phenotype_visualization(sample_id, output_dir=None):
    """Run phenotype tile visualization - 4 separate grayscale channel images."""
    import napari
    import tifffile
    from qtpy.QtWidgets import QApplication
    import time

    # Load data
    data = load_phenotype_data(sample_id)

    if not data.get("aligned"):
        print("Error: No aligned image found for this sample")
        return

    # Read image
    print("\nReading image file...")
    cfg = PHENOTYPE_CONFIG

    image = tifffile.imread(data["aligned"])
    print(f"  Image shape: {image.shape}")

    # Extract channels (assuming shape is H x W x C or C x H x W)
    if image.ndim == 3:
        if image.shape[2] <= 5:  # Channels last (H, W, C)
            channels = [image[:, :, i] for i in range(min(4, image.shape[2]))]
            height, width = image.shape[0], image.shape[1]
        else:  # Channels first (C, H, W)
            channels = [image[i, :, :] for i in range(min(4, image.shape[0]))]
            height, width = image.shape[1], image.shape[2]
    else:
        channels = [image]
        height, width = image.shape[0], image.shape[1]

    # Get random center position (same for all channels)
    center_position = get_random_center(height, width, cfg["margin"])
    print(f"  Random center position: {center_position}")

    # Setup output directory
    if output_dir is None:
        output_dir = LOCAL_CACHE / "screenshots" / "phenotype"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create separate viewer for each channel
    for ch_idx, (ch_data, ch_name) in enumerate(zip(channels, cfg["channel_names"])):
        print(f"\nProcessing channel: {ch_name}")

        viewer = napari.Viewer()

        # Add grayscale image with manual contrast limits from config
        viewer.add_image(
            ch_data,
            name=ch_name,
            colormap="gray",
            contrast_limits=cfg["contrast_limits"].get(
                ch_name, [ch_data.min(), ch_data.max()]
            ),
        )

        # Configure viewer
        viewer.scale_bar.unit = "μm"
        viewer.scale_bar.position = "bottom_right"
        viewer.scale_bar.ticks = False
        viewer.window.resize(cfg["window_size"], cfg["window_size"])

        # Set zoom and position (random center, same for all channels)
        zoom_level = cfg["window_size"] / cfg["view_width"]
        viewer.camera.center = center_position
        viewer.camera.zoom = zoom_level

        # Hide UI panels
        viewer.window.qt_viewer.dockLayerList.setVisible(False)
        viewer.window.qt_viewer.dockLayerControls.setVisible(False)

        # Take screenshots with different scale bar modes
        scale_bar_modes = [
            ("with_text", True, 20),
            ("no_text", True, 0),
            ("no_scalebar", False, 0),
        ]

        canvas = viewer.window._qt_viewer.canvas.native

        for mode_name, sb_visible, sb_font_size in scale_bar_modes:
            viewer.scale_bar.visible = sb_visible
            if sb_visible:
                viewer.scale_bar.font_size = sb_font_size

            QApplication.processEvents()
            time.sleep(0.3)

            path = output_dir / f"{sample_id}_{ch_name}_{mode_name}.png"
            canvas.grab().save(str(path))
            print(f"  Saved: {path.name}")

        viewer.close()

    print(f"\nScreenshots saved to: {output_dir}")

    # Upload to remote
    print("\nUploading screenshots to remote...")
    remote_output = f"{REMOTE_BASE}/visualizations/screenshots/phenotype"
    for png_file in output_dir.glob(f"{sample_id}_*.png"):
        upload_file(png_file, f"{remote_output}/{png_file.name}")
    print("Upload complete!")


# ============================================================================
# SBS TILE VISUALIZATION
# ============================================================================


def list_sbs_samples():
    """List available SBS samples on remote (random selection)."""
    try:
        output = run_ssh_command(
            f"ls {REMOTE_SBS_DIR}/*__aligned.tiff 2>/dev/null | shuf | head -20"
        )
        files = [Path(f).name for f in output.split("\n") if f]
        return [f.replace("__aligned.tiff", "") for f in files]
    except RuntimeError:
        return []


def load_sbs_data(sample_id):
    """Download SBS aligned image."""
    print(f"\nLoading SBS data for {sample_id}...")

    aligned_remote = f"{REMOTE_SBS_DIR}/{sample_id}__aligned.tiff"
    try:
        aligned_local = download_file(aligned_remote)
        return {"aligned": aligned_local}
    except RuntimeError as e:
        print(f"  Error: {e}")
        return {"aligned": None}


def run_sbs_visualization(sample_id, output_dir=None):
    """Run SBS tile visualization - DAPI + GTAC composite overlay."""
    import napari
    import tifffile
    from qtpy.QtWidgets import QApplication
    import time

    # Load data
    data = load_sbs_data(sample_id)

    if not data.get("aligned"):
        print("Error: No aligned image found for this sample")
        return

    # Read image
    print("\nReading image file...")
    cfg = SBS_CONFIG

    aligned_image = tifffile.imread(data["aligned"])
    print(f"  Image shape: {aligned_image.shape}")

    # Extract channels from aligned image
    # Shape is typically (cycles, channels, H, W)
    cycle_idx = cfg["cycle_index"]

    if aligned_image.ndim == 4:
        # (cycles, channels, H, W)
        dapi = aligned_image[cycle_idx, 0, :, :]
        channel_data = [aligned_image[cycle_idx, i + 1, :, :] for i in range(4)]
        height, width = aligned_image.shape[2], aligned_image.shape[3]
    elif aligned_image.ndim == 3:
        if aligned_image.shape[0] <= 5:
            # (channels, H, W)
            dapi = aligned_image[0, :, :]
            channel_data = [
                aligned_image[i + 1, :, :]
                for i in range(min(4, aligned_image.shape[0] - 1))
            ]
            height, width = aligned_image.shape[1], aligned_image.shape[2]
        else:
            # (H, W, channels)
            dapi = aligned_image[:, :, 0]
            channel_data = [
                aligned_image[:, :, i + 1]
                for i in range(min(4, aligned_image.shape[2] - 1))
            ]
            height, width = aligned_image.shape[0], aligned_image.shape[1]
    else:
        print("Error: Unexpected image dimensions")
        return

    # Get random center position
    center_position = get_random_center(height, width, cfg["margin"])
    print(f"  Random center position: {center_position}")

    # Setup output directory
    if output_dir is None:
        output_dir = LOCAL_CACHE / "screenshots" / "sbs"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create viewer
    print("\nStarting napari viewer...")
    viewer = napari.Viewer()

    # Add DAPI channel (gray)
    viewer.add_image(
        dapi,
        name="DAPI",
        colormap="gray",
        blending="additive",
        contrast_limits=[dapi.min(), dapi.max() * cfg["dapi_contrast"]],
        opacity=cfg["dapi_opacity"],
    )

    # Add GTAC channels with colors
    for i, (ch_name, ch_color) in enumerate(
        zip(cfg["channel_names"], cfg["channel_colors"])
    ):
        if i < len(channel_data):
            limits = cfg["contrast_limits"].get(ch_name, [0, 30000])
            viewer.add_image(
                channel_data[i],
                name=ch_name,
                colormap=ch_color,
                blending="additive",
                contrast_limits=limits,
            )

    # Configure viewer
    viewer.scale_bar.unit = "μm"
    viewer.scale_bar.position = "bottom_right"
    viewer.scale_bar.ticks = False
    viewer.window.resize(cfg["window_size"], cfg["window_size"])

    # Set zoom and position (random center)
    zoom_level = cfg["window_size"] / cfg["view_width"]
    viewer.camera.center = center_position
    viewer.camera.zoom = zoom_level

    # Hide UI panels
    viewer.window.qt_viewer.dockLayerList.setVisible(False)
    viewer.window.qt_viewer.dockLayerControls.setVisible(False)

    # Take screenshots with different scale bar modes
    scale_bar_modes = [
        ("with_text", True, 20),
        ("no_text", True, 0),
        ("no_scalebar", False, 0),
    ]

    canvas = viewer.window._qt_viewer.canvas.native

    for mode_name, sb_visible, sb_font_size in scale_bar_modes:
        viewer.scale_bar.visible = sb_visible
        if sb_visible:
            viewer.scale_bar.font_size = sb_font_size

        QApplication.processEvents()
        time.sleep(0.3)

        path = output_dir / f"{sample_id}_GTAC_{mode_name}.png"
        canvas.grab().save(str(path))
        print(f"  Saved: {path.name}")

    viewer.close()

    print(f"\nScreenshots saved to: {output_dir}")

    # Upload to remote
    print("\nUploading screenshots to remote...")
    remote_output = f"{REMOTE_BASE}/visualizations/screenshots/sbs"
    for png_file in output_dir.glob(f"{sample_id}_*.png"):
        upload_file(png_file, f"{remote_output}/{png_file.name}")
    print("Upload complete!")


# ============================================================================
# ILLUMINATION CORRECTION VISUALIZATION
# ============================================================================


def list_ic_fields():
    """List available IC field files on remote."""
    try:
        output = run_ssh_command(
            f"ls {REMOTE_IC_FIELDS_DIR}/ | grep '__ic_field.tiff' | shuf | head -20"
        )
        files = [f for f in output.split("\n") if f]
        return files
    except RuntimeError:
        return []


def load_illum_data(sample_id=None):
    """Download illumination correction function and a corresponding raw image.

    IC fields are per plate-well (e.g., P-1_W-A1__ic_field.tiff).
    Raw images are per tile (e.g., P-1_W-A1_T-0__image.tiff).
    """
    print("\nLoading illumination correction data...")

    data = {"illum_function": None, "raw_image": None, "corrected_image": None}

    # Get list of IC fields
    ic_files = list_ic_fields()
    if not ic_files:
        print("  Warning: No IC field files found")
        return data

    # Pick a random IC field or use specified sample
    import random

    if sample_id:
        # sample_id might be full tile ID like P-1_W-A1_T-0, extract plate-well
        parts = sample_id.split("_T-")
        plate_well = parts[0]  # e.g., P-1_W-A1
        ic_file = f"{plate_well}__ic_field.tiff"
    else:
        ic_file = random.choice(ic_files)
        # Extract plate-well from IC file name (e.g., P-1_W-A1__ic_field.tiff -> P-1_W-A1)
        plate_well = ic_file.replace("__ic_field.tiff", "")

    # Download IC field
    ic_remote = f"{REMOTE_IC_FIELDS_DIR}/{ic_file}"
    try:
        data["illum_function"] = download_file(ic_remote)
        data["ic_name"] = ic_file
        print(f"  Downloaded IC field: {ic_file}")
    except RuntimeError as e:
        print(f"  Error downloading IC field: {e}")
        return data

    # Find a random tile for this plate-well from preprocess images (raw)
    try:
        output = run_ssh_command(
            f"ls {REMOTE_PREPROCESS_DIR}/ | grep '{plate_well}_T-' | grep '__image.tiff' | shuf | head -1"
        )
        if output.strip():
            raw_image_file = output.strip()
            raw_remote = f"{REMOTE_PREPROCESS_DIR}/{raw_image_file}"
            data["raw_image"] = download_file(raw_remote)
            # Extract tile sample ID (e.g., P-1_W-A1_T-0)
            tile_sample = raw_image_file.replace("__image.tiff", "")
            data["tile_sample"] = tile_sample
            print(f"  Downloaded raw image: {raw_image_file}")

            # Download the matching corrected (aligned) image from phenotype output
            corrected_remote = f"{REMOTE_PHENOTYPE_DIR}/{tile_sample}__aligned.tiff"
            try:
                data["corrected_image"] = download_file(corrected_remote)
                print(f"  Downloaded corrected image: {tile_sample}__aligned.tiff")
            except RuntimeError as e:
                print(f"  Warning: Could not download corrected image: {e}")
    except RuntimeError as e:
        print(f"  Error downloading raw image: {e}")

    return data


def run_illum_visualization(sample_id=None, output_dir=None):
    """Run illumination correction visualization - illum function + example DAPI image."""
    import napari
    import tifffile
    from qtpy.QtWidgets import QApplication
    import time

    # Load data
    data = load_illum_data(sample_id)

    if not data.get("illum_function") and not data.get("raw_image"):
        print("Error: No illumination data found")
        return

    cfg = ILLUM_CONFIG

    # Setup output directory
    if output_dir is None:
        output_dir = LOCAL_CACHE / "screenshots" / "illum"
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Helper to extract DAPI channel from multi-channel image
    def extract_dapi(img):
        if img.ndim == 3:
            return img[0] if img.shape[0] <= 5 else img[:, :, 0]
        return img

    # Images to process: (image_array, title, filename_prefix)
    images_to_process = []

    # Load and prepare illumination function
    if data.get("illum_function"):
        illum_img = tifffile.imread(data["illum_function"])
        print(f"  Illum function shape: {illum_img.shape}")
        illum_dapi = extract_dapi(illum_img)
        images_to_process.append(
            (
                illum_dapi,
                "Illumination correction function, DAPI channel",
                "illum_function",
            )
        )

    # Load raw DAPI (before IC)
    if data.get("raw_image"):
        raw_img = tifffile.imread(data["raw_image"])
        print(f"  Raw image shape: {raw_img.shape}")
        raw_dapi = extract_dapi(raw_img)
        images_to_process.append(
            (raw_dapi, "Raw DAPI (before illumination correction)", "raw_dapi")
        )

    # Load corrected DAPI (after IC) — the aligned phenotype image
    if data.get("corrected_image"):
        corrected_img = tifffile.imread(data["corrected_image"])
        print(f"  Corrected image shape: {corrected_img.shape}")
        corrected_dapi = extract_dapi(corrected_img)
        images_to_process.append(
            (
                corrected_dapi,
                "Corrected DAPI (after illumination correction)",
                "corrected_dapi",
            )
        )

    # Scale bar modes
    scale_bar_modes = [
        ("with_text", True, 20),
        ("no_text", True, 0),
        ("no_scalebar", False, 0),
    ]

    # Process each image
    for img_data, title, prefix in images_to_process:
        print(f"\nProcessing: {title}")
        height, width = img_data.shape

        for mode_name, sb_visible, sb_font_size in scale_bar_modes:
            viewer = napari.Viewer()

            # Add grayscale image (full field, no zoom)
            viewer.add_image(
                img_data,
                name="image",
                colormap="gray",
                contrast_limits=[img_data.min(), img_data.max()],
            )

            # Configure viewer - full field view matching image aspect ratio
            viewer.scale_bar.unit = "μm"
            viewer.scale_bar.position = "bottom_right"
            viewer.scale_bar.ticks = False
            viewer.scale_bar.visible = sb_visible
            if sb_visible:
                viewer.scale_bar.font_size = sb_font_size

            # Size window to match image (scaled to fit within max size)
            max_dim = cfg["window_size"]
            if height >= width:
                win_h = max_dim
                win_w = int(max_dim * width / height)
            else:
                win_w = max_dim
                win_h = int(max_dim * height / width)
            viewer.window.resize(win_w, win_h)

            # Set camera to show exactly the image bounds (no extra padding)
            viewer.camera.center = (height / 2, width / 2)
            viewer.camera.zoom = min(win_h / height, win_w / width)

            # Hide UI panels
            viewer.window.qt_viewer.dockLayerList.setVisible(False)
            viewer.window.qt_viewer.dockLayerControls.setVisible(False)

            # Take screenshot
            QApplication.processEvents()
            time.sleep(0.5)

            canvas = viewer.window._qt_viewer.canvas.native
            path = output_dir / f"{prefix}_{mode_name}.png"
            canvas.grab().save(str(path))
            print(f"  Saved: {path.name}")

            viewer.close()

    print(f"\nScreenshots saved to: {output_dir}")

    # Upload to remote
    print("\nUploading screenshots to remote...")
    remote_output = f"{REMOTE_BASE}/visualizations/screenshots/illum"
    for png_file in output_dir.glob("*.png"):
        upload_file(png_file, f"{remote_output}/{png_file.name}")
    print("Upload complete!")


# ============================================================================
# MAIN
# ============================================================================


def parse_args():
    parser = argparse.ArgumentParser(
        description="Unified benchmark visualization with remote data streaming",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Benchmark visualizations
    python visualization.py --mode segmentation --sample P-1_W-B3_T-169
    python visualization.py --mode spots --sample P-1_W-B1_T-2

    # Raw tile visualizations (random sample if not specified)
    python visualization.py --mode phenotype
    python visualization.py --mode sbs --sample P-1_W-A1_T-0
    python visualization.py --mode illum

    # Utilities
    python visualization.py --list-samples
    python visualization.py --clear-cache
        """,
    )
    parser.add_argument(
        "--mode",
        choices=["segmentation", "spots", "phenotype", "sbs", "illum", "all"],
        help="Visualization mode ('all' runs each mode sequentially)",
    )
    parser.add_argument(
        "--sample",
        type=str,
        help="Sample ID (e.g., P-1_W-B3_T-169)",
    )
    parser.add_argument(
        "--list-samples",
        action="store_true",
        help="List available samples on remote",
    )
    parser.add_argument(
        "--clear-cache",
        action="store_true",
        help="Clear local file cache",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        help="Local output directory for screenshots",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print("=" * 60)
    print("BENCHMARK VISUALIZATION (Napari + Remote Data)")
    print("=" * 60)
    print(f"Remote: {REMOTE_HOST}:{REMOTE_BASE}")
    print(f"Cache:  {LOCAL_CACHE}")
    print()

    if args.clear_cache:
        import shutil

        if LOCAL_CACHE.exists():
            shutil.rmtree(LOCAL_CACHE)
            print("Cache cleared.")
        return

    if args.list_samples:
        print("Segmentation samples (benchmark):")
        for s in list_segmentation_samples():
            print(f"  {s}")
        print()
        print("Spot calling samples (benchmark):")
        for s in list_spot_samples():
            print(f"  {s}")
        print()
        print("Phenotype samples (random 20):")
        for s in list_phenotype_samples()[:5]:
            print(f"  {s}")
        print("  ...")
        print()
        print("SBS samples (random 20):")
        for s in list_sbs_samples()[:5]:
            print(f"  {s}")
        print("  ...")
        return

    if not args.mode:
        print(
            "Error: --mode is required (segmentation, spots, phenotype, sbs, illum, or all)"
        )
        print("Use --list-samples to see available samples")
        return

    modes = (
        ["segmentation", "spots", "phenotype", "sbs", "illum"]
        if args.mode == "all"
        else [args.mode]
    )

    for mode in modes:
        print(f"\n{'=' * 60}")
        print(f"  MODE: {mode}")
        print(f"{'=' * 60}")

        if mode == "segmentation":
            samples = list_segmentation_samples()
            sample = args.sample or (samples[0] if samples else None)
            if not sample:
                print("Error: No segmentation samples found")
                continue
            if sample not in samples:
                print(f"Warning: {sample} not in available samples, trying anyway...")
            run_segmentation_visualization(sample, args.output_dir)

        elif mode == "spots":
            samples = list_spot_samples()
            sample = args.sample or (samples[0] if samples else None)
            if not sample:
                print("Error: No spot calling samples found")
                continue
            if sample not in samples:
                print(f"Warning: {sample} not in available samples, trying anyway...")
            run_spot_visualization(sample, args.output_dir)

        elif mode == "phenotype":
            if args.sample:
                sample = args.sample
            else:
                samples = list_phenotype_samples()
                if not samples:
                    print("Error: No phenotype samples found")
                    continue
                import random

                sample = random.choice(samples)
            run_phenotype_visualization(sample, args.output_dir)

        elif mode == "sbs":
            if args.sample:
                sample = args.sample
            else:
                samples = list_sbs_samples()
                if not samples:
                    print("Error: No SBS samples found")
                    continue
                import random

                sample = random.choice(samples)
            run_sbs_visualization(sample, args.output_dir)

        elif mode == "illum":
            run_illum_visualization(args.sample, args.output_dir)


if __name__ == "__main__":
    main()
