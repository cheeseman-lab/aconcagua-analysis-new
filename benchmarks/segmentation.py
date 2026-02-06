"""
Segmentation Benchmark

Compares cell/nuclei segmentation methods:
- Cellpose 3 (cyto3 model)
- Cellpose 4 (CPSAM model - requires separate environment)
- StarDist (2D_versatile_fluo)

Metrics: runtime, memory, nuclei/cell counts, retention rates, pipeline efficiency.

Usage:
    python segmentation.py                     # Run all available methods
    python segmentation.py --method cellpose   # Run Cellpose (v3, cyto3)
    python segmentation.py --method cellpose4  # Run Cellpose 4 (CPSAM)
    python segmentation.py --method stardist   # Run only StarDist
"""

import os
import sys
import time
import random
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tifffile import imread, imwrite
from pathlib import Path
from datetime import datetime
import psutil
import gc
from tqdm import tqdm

# Add brieflow to path (relative to benchmarks directory)
BENCHMARKS_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BENCHMARKS_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT / "brieflow" / "workflow"))


# Parse command line arguments first to determine which imports are needed
def parse_args():
    parser = argparse.ArgumentParser(description="Run segmentation benchmarks")
    parser.add_argument(
        "--method",
        type=str,
        choices=["cellpose", "cellpose4", "stardist", "all"],
        default="all",
        help="Segmentation method to benchmark (default: all available)",
    )
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Only regenerate plots from existing results",
    )
    return parser.parse_args()


# Parse args at module level so imports can be conditional
args = parse_args()
SELECTED_METHOD = args.method

# Conditionally import segmentation methods based on --method argument
CELLPOSE_AVAILABLE = False
CELLPOSE4_AVAILABLE = False
STARDIST_AVAILABLE = False

segment_cellpose = None
segment_stardist = None

if SELECTED_METHOD in ["cellpose", "all"]:
    try:
        from lib.shared.segment_cellpose import segment_cellpose

        CELLPOSE_AVAILABLE = True
        print("Cellpose 3 loaded successfully")
    except ImportError as e:
        print(f"Warning: Cellpose not available ({e}). Skipping Cellpose benchmark.")

if SELECTED_METHOD in ["cellpose4", "all"]:
    try:
        from lib.shared.segment_cellpose import segment_cellpose

        CELLPOSE4_AVAILABLE = True
        print("Cellpose 4 (CPSAM) loaded successfully")
    except ImportError as e:
        print(
            f"Warning: Cellpose 4 not available ({e}). Skipping Cellpose 4 benchmark."
        )

if SELECTED_METHOD in ["stardist", "all"]:
    try:
        from lib.shared.segment_stardist import segment_stardist

        STARDIST_AVAILABLE = True
        print("StarDist loaded successfully")
    except ImportError as e:
        print(f"Warning: StarDist not available ({e}). Skipping StarDist benchmark.")

random.seed(42)

# Configure paths (relative to project root)
SOURCE_DIR = PROJECT_ROOT / "analysis" / "brieflow_output" / "phenotype" / "images"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "segmentation"
ALIGNED_DIR = OUTPUT_DIR / "aligned"

# Create directories
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)
ALIGNED_DIR.mkdir(exist_ok=True, parents=True)

# Create subdirectories for each segmentation method
for method in ["cellpose", "cellpose4", "stardist"]:
    (OUTPUT_DIR / method).mkdir(exist_ok=True)

# Parameters for segmentation methods (from config.yml)
CELLPOSE_PARAMS = {
    "dapi_index": 0,
    "cyto_index": 3,
    "nuclei_diameter": 42.684959745981,
    "cell_diameter": 68.05913211132304,
    "cellpose_model": "cyto3",
    "nuclei_flow_threshold": 0.4,
    "nuclei_cellprob_threshold": 0.0,
    "cell_flow_threshold": 1,
    "cell_cellprob_threshold": 0,
    "reconcile": "contained_in_cells",
    "return_counts": True,
    "gpu": False,
    "segment_cells": True,
}

# Cellpose 4 with CPSAM model (requires cellpose 4.x)
CELLPOSE4_PARAMS = {
    "dapi_index": 0,
    "cyto_index": 3,
    "nuclei_diameter": 42.684959745981,
    "cell_diameter": 68.05913211132304,
    "cellpose_model": "cpsam",  # CPSAM model in Cellpose 4.x
    "nuclei_flow_threshold": 0.4,
    "nuclei_cellprob_threshold": 0.0,
    "cell_flow_threshold": 1,
    "cell_cellprob_threshold": 0,
    "reconcile": "contained_in_cells",
    "return_counts": True,
    "gpu": False,
    "segment_cells": True,
    "helper_index": 1,  # TUBULIN channel to help with CPSAM segmentation
}

STARDIST_PARAMS = {
    "dapi_index": 0,
    "cyto_index": 3,
    "stardist_model": "2D_versatile_fluo",
    "nuclei_prob_threshold": 0.479071,
    "nuclei_nms_threshold": 0.3,
    "cell_prob_threshold": 0.479071,
    "cell_nms_threshold": 0.3,
    "reconcile": "contained_in_cells",
    "return_counts": True,
    "gpu": False,
    "segment_cells": True,
}


def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024  # Convert to MB


def prepare_test_images():
    """Find one aligned image per plate randomly and copy to benchmark directory."""
    # Pattern to match aligned image files
    pattern = str(SOURCE_DIR / "*__aligned.tiff")
    all_images = glob.glob(pattern)
    print(f"Found {len(all_images)} total aligned images")

    # Extract plate info from filenames
    images_by_plate = {}
    for img_path in all_images:
        img_name = os.path.basename(img_path)
        if img_name.startswith("P-"):
            plate = img_name.split("_")[0].replace("P-", "")
            if plate not in images_by_plate:
                images_by_plate[plate] = []
            images_by_plate[plate].append(img_path)

    # Select a random image from each plate
    selected_images = []
    for plate, images in images_by_plate.items():
        if images:
            # Choose a random image from this plate
            random_image = random.choice(images)
            dst_path = ALIGNED_DIR / os.path.basename(random_image)
            print(
                f"Copying random tile from plate {plate}: {os.path.basename(random_image)}"
            )
            if not dst_path.exists():
                os.system(f"cp {random_image} {dst_path}")
            selected_images.append(dst_path)
        else:
            print(f"No images found for plate {plate}")

    print(f"\nCopied {len(selected_images)} images to {ALIGNED_DIR}")
    return selected_images


def collect_segmentation_stats():
    """Collect all segmentation stats into a single benchmark results file."""
    all_stats = []

    for method in ["cellpose", "cellpose4", "stardist"]:
        stat_files = list((OUTPUT_DIR / method).glob("*_segmentation_stats.tsv"))
        for file in stat_files:
            try:
                stats = pd.read_csv(file, sep="\t")
                all_stats.append(stats)
            except Exception as e:
                print(f"Error reading {file}: {e}")

    if all_stats:
        # Combine all stats into a single DataFrame
        benchmark_results = pd.concat(all_stats, ignore_index=True)
        # Save to CSV
        benchmark_results.to_csv(OUTPUT_DIR / "benchmark_results.csv", index=False)
        print(f"Saved benchmark results to {OUTPUT_DIR / 'benchmark_results.csv'}")
        return benchmark_results
    else:
        print("No segmentation stats found.")
        return None


def benchmark_segmentation_methods():
    """Run benchmark tests on all three segmentation methods."""
    # Get test images
    image_paths = list(ALIGNED_DIR.glob("*.tiff"))
    if not image_paths:
        print("No images found. Preparing test images...")
        image_paths = prepare_test_images()

    print(f"Found {len(image_paths)} aligned images for benchmarking")

    # Process each image
    for img_path in tqdm(image_paths, desc="Processing images"):
        img_name = img_path.name

        # Extract plate, well, tile info from filename
        if "P-" in img_name and "W-" in img_name and "T-" in img_name:
            plate = img_name.split("_")[0].replace("P-", "")
            well = img_name.split("_")[1].replace("W-", "")
            tile = img_name.split("_")[2].replace("T-", "")
        else:
            # Default values if pattern doesn't match
            plate = "unknown"
            well = "unknown"
            tile = "unknown"

        print(f"\nProcessing plate {plate}, well {well}, tile {tile}")

        # Load image - these are already aligned
        try:
            image = imread(img_path)
            print(f"  Image shape: {image.shape}")
        except Exception as e:
            print(f"  Error loading image: {e}")
            continue

        # Benchmark Cellpose
        if CELLPOSE_AVAILABLE:
            print("  Running Cellpose...")
            try:
                # Measure memory before GC to establish true baseline
                gc.collect()
                start_mem = get_memory_usage()
                start_time = time.time()

                # Get segmentation with reconciliation
                nuclei_cp, cells_cp, counts_df = segment_cellpose(
                    data=image,
                    dapi_index=CELLPOSE_PARAMS["dapi_index"],
                    cyto_index=CELLPOSE_PARAMS["cyto_index"],
                    nuclei_diameter=CELLPOSE_PARAMS["nuclei_diameter"],
                    cell_diameter=CELLPOSE_PARAMS["cell_diameter"],
                    cyto_model=CELLPOSE_PARAMS["cellpose_model"],
                    cellpose_kwargs=dict(
                        nuclei_flow_threshold=CELLPOSE_PARAMS["nuclei_flow_threshold"],
                        nuclei_cellprob_threshold=CELLPOSE_PARAMS[
                            "nuclei_cellprob_threshold"
                        ],
                        cell_flow_threshold=CELLPOSE_PARAMS["cell_flow_threshold"],
                        cell_cellprob_threshold=CELLPOSE_PARAMS[
                            "cell_cellprob_threshold"
                        ],
                    ),
                    reconcile=CELLPOSE_PARAMS["reconcile"],
                    return_counts=CELLPOSE_PARAMS["return_counts"],
                    gpu=CELLPOSE_PARAMS["gpu"],
                    cells=CELLPOSE_PARAMS["segment_cells"],
                )

                cellpose_time = time.time() - start_time
                end_mem = get_memory_usage()
                mem_usage = max(0, end_mem - start_mem)

                # Save segmentation results
                nuclei_out = (
                    OUTPUT_DIR / "cellpose" / f"P-{plate}_W-{well}_T-{tile}_nuclei.tiff"
                )
                cells_out = (
                    OUTPUT_DIR / "cellpose" / f"P-{plate}_W-{well}_T-{tile}_cells.tiff"
                )
                segmentation_stats_out = (
                    OUTPUT_DIR
                    / "cellpose"
                    / f"P-{plate}_W-{well}_T-{tile}_segmentation_stats.tsv"
                )

                # Add info to the beginning of the counts_df (metadata)
                counts_df.insert(0, "method", "cellpose")
                counts_df.insert(1, "plate", plate)
                counts_df.insert(2, "well", well)
                counts_df.insert(3, "tile", tile)

                # Add metrics to the end of the counts_df (performance metrics)
                counts_df["runtime_seconds"] = cellpose_time
                counts_df["memory_mb"] = mem_usage
                counts_df["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                imwrite(nuclei_out, nuclei_cp.astype(np.uint16))
                imwrite(cells_out, cells_cp.astype(np.uint16))
                counts_df.to_csv(segmentation_stats_out, index=False, sep="\t")

                # Clean up to reduce memory usage
                del nuclei_cp, cells_cp, counts_df
                gc.collect()

            except Exception as e:
                print(f"  Error running Cellpose: {e}")

        # Benchmark StarDist
        if STARDIST_AVAILABLE:
            print("  Running StarDist...")
            try:
                # Measure memory before GC to establish true baseline
                gc.collect()
                start_mem = get_memory_usage()
                start_time = time.time()

                # Get segmentation with reconciliation
                nuclei_sd, cells_sd, counts_df = segment_stardist(
                    data=image,
                    dapi_index=STARDIST_PARAMS["dapi_index"],
                    cyto_index=STARDIST_PARAMS["cyto_index"],
                    model_type=STARDIST_PARAMS["stardist_model"],
                    stardist_kwargs=dict(
                        nuclei_prob_threshold=STARDIST_PARAMS["nuclei_prob_threshold"],
                        nuclei_nms_threshold=STARDIST_PARAMS["nuclei_nms_threshold"],
                        cell_prob_threshold=STARDIST_PARAMS["cell_prob_threshold"],
                        cell_nms_threshold=STARDIST_PARAMS["cell_nms_threshold"],
                    ),
                    reconcile=STARDIST_PARAMS["reconcile"],
                    return_counts=STARDIST_PARAMS["return_counts"],
                    gpu=STARDIST_PARAMS["gpu"],
                    cells=STARDIST_PARAMS["segment_cells"],
                )

                stardist_time = time.time() - start_time
                end_mem = get_memory_usage()
                mem_usage = max(0, end_mem - start_mem)

                # Save segmentation results
                nuclei_out = (
                    OUTPUT_DIR / "stardist" / f"P-{plate}_W-{well}_T-{tile}_nuclei.tiff"
                )
                cells_out = (
                    OUTPUT_DIR / "stardist" / f"P-{plate}_W-{well}_T-{tile}_cells.tiff"
                )
                segmentation_stats_out = (
                    OUTPUT_DIR
                    / "stardist"
                    / f"P-{plate}_W-{well}_T-{tile}_segmentation_stats.tsv"
                )

                # Add info to the beginning of the counts_df (metadata)
                counts_df.insert(0, "method", "stardist")
                counts_df.insert(1, "plate", plate)
                counts_df.insert(2, "well", well)
                counts_df.insert(3, "tile", tile)

                # Add metrics to the end of the counts_df (performance metrics)
                counts_df["runtime_seconds"] = stardist_time
                counts_df["memory_mb"] = mem_usage
                counts_df["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                imwrite(nuclei_out, nuclei_sd.astype(np.uint16))
                imwrite(cells_out, cells_sd.astype(np.uint16))
                counts_df.to_csv(segmentation_stats_out, index=False, sep="\t")

                # Clean up to reduce memory usage
                del nuclei_sd, cells_sd, counts_df
                gc.collect()

            except Exception as e:
                print(f"  Error running StarDist: {e}")

        # Benchmark Cellpose 4 (CPSAM)
        if CELLPOSE4_AVAILABLE:
            print("  Running Cellpose 4 (CPSAM)...")
            try:
                # Measure memory before GC to establish true baseline
                gc.collect()
                start_mem = get_memory_usage()
                start_time = time.time()

                # Get segmentation with reconciliation using CPSAM model
                nuclei_cp4, cells_cp4, counts_df = segment_cellpose(
                    data=image,
                    dapi_index=CELLPOSE4_PARAMS["dapi_index"],
                    cyto_index=CELLPOSE4_PARAMS["cyto_index"],
                    nuclei_diameter=CELLPOSE4_PARAMS["nuclei_diameter"],
                    cell_diameter=CELLPOSE4_PARAMS["cell_diameter"],
                    cyto_model=CELLPOSE4_PARAMS["cellpose_model"],
                    cellpose_kwargs=dict(
                        nuclei_flow_threshold=CELLPOSE4_PARAMS["nuclei_flow_threshold"],
                        nuclei_cellprob_threshold=CELLPOSE4_PARAMS[
                            "nuclei_cellprob_threshold"
                        ],
                        cell_flow_threshold=CELLPOSE4_PARAMS["cell_flow_threshold"],
                        cell_cellprob_threshold=CELLPOSE4_PARAMS[
                            "cell_cellprob_threshold"
                        ],
                    ),
                    helper_index=CELLPOSE4_PARAMS["helper_index"],
                    reconcile=CELLPOSE4_PARAMS["reconcile"],
                    return_counts=CELLPOSE4_PARAMS["return_counts"],
                    gpu=CELLPOSE4_PARAMS["gpu"],
                    cells=CELLPOSE4_PARAMS["segment_cells"],
                )

                cellpose4_time = time.time() - start_time
                end_mem = get_memory_usage()
                mem_usage = max(0, end_mem - start_mem)

                # Save segmentation results
                nuclei_out = (
                    OUTPUT_DIR
                    / "cellpose4"
                    / f"P-{plate}_W-{well}_T-{tile}_nuclei.tiff"
                )
                cells_out = (
                    OUTPUT_DIR / "cellpose4" / f"P-{plate}_W-{well}_T-{tile}_cells.tiff"
                )
                segmentation_stats_out = (
                    OUTPUT_DIR
                    / "cellpose4"
                    / f"P-{plate}_W-{well}_T-{tile}_segmentation_stats.tsv"
                )

                # Add info to the beginning of the counts_df (metadata)
                counts_df.insert(0, "method", "cellpose4")
                counts_df.insert(1, "plate", plate)
                counts_df.insert(2, "well", well)
                counts_df.insert(3, "tile", tile)

                # Add metrics to the end of the counts_df (performance metrics)
                counts_df["runtime_seconds"] = cellpose4_time
                counts_df["memory_mb"] = mem_usage
                counts_df["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

                imwrite(nuclei_out, nuclei_cp4.astype(np.uint16))
                imwrite(cells_out, cells_cp4.astype(np.uint16))
                counts_df.to_csv(segmentation_stats_out, index=False, sep="\t")

                # Clean up to reduce memory usage
                del nuclei_cp4, cells_cp4, counts_df
                gc.collect()

            except Exception as e:
                print(f"  Error running Cellpose 4 (CPSAM): {e}")

        # Clean up image to reduce memory usage
        del image
        gc.collect()

    # Collect all segmentation stats into benchmark results
    results_df = collect_segmentation_stats()

    # Generate summary statistics
    if results_df is not None and len(results_df) > 0:
        # Generate summary visualizations
        generate_benchmark_summary(results_df)
    else:
        print("No results generated. Check for errors.")

    return results_df


def generate_benchmark_summary(results_df):
    """Generate summary statistics and plots from benchmark results"""
    import numpy as np
    from plot_style import setup_plot_style

    # Apply consistent plot styling
    setup_plot_style()

    # Use viridis color map for bars
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(results_df["method"].unique())))

    # Sort methods alphabetically
    results_df = results_df.sort_values("method")

    # Calculate retention rates if they don't exist already
    if "nuclei_retention" not in results_df.columns:
        # Calculate retention rates, avoiding division by zero
        results_df["nuclei_retention"] = np.where(
            results_df["initial_nuclei"] > 0,
            results_df["final_nuclei"] / results_df["initial_nuclei"] * 100,
            0,
        )

    if "cell_retention" not in results_df.columns:
        results_df["cell_retention"] = np.where(
            results_df["initial_cells"] > 0,
            results_df["final_cells"] / results_df["initial_cells"] * 100,
            0,
        )

    # Calculate summary statistics by method
    summary = (
        results_df.groupby("method")
        .agg(
            {
                "initial_nuclei": ["mean", "std", "min", "max"],
                "initial_cells": ["mean", "std", "min", "max"],
                "after_edge_removal_nuclei": ["mean", "std", "min", "max"],
                "after_edge_removal_cells": ["mean", "std", "min", "max"],
                "final_nuclei": ["mean", "std", "min", "max"],
                "final_cells": ["mean", "std", "min", "max"],
                "nuclei_retention": ["mean", "std", "min", "max"],
                "cell_retention": ["mean", "std", "min", "max"],
                "runtime_seconds": ["mean", "std", "min", "max"],
                "memory_mb": ["mean", "std", "min", "max"],
            }
        )
        .reset_index()
    )

    # Save summary to file
    summary.to_csv(OUTPUT_DIR / "benchmark_summary.csv")

    # Define a function for adding value labels to bars
    def add_value_labels(ax, fmt):
        for container in ax.containers:
            ax.bar_label(container, fmt=fmt)

    # Create runtime comparison plot
    plt.figure(figsize=(8, 5))
    # Create a dictionary that maps methods to colors
    method_colors = dict(zip(sorted(results_df["method"].unique()), colors))
    # Plot bars with specific colors based on method
    ax = sns.barplot(
        x="method",
        y="runtime_seconds",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        palette=method_colors,
    )
    add_value_labels(ax, "%.2f")
    plt.title("Segmentation Method Runtime Comparison", fontweight="bold")
    plt.ylabel("Runtime (seconds)")
    plt.xlabel("Segmentation Method")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "runtime_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Create memory usage comparison plot
    plt.figure(figsize=(8, 5))
    ax = sns.barplot(
        x="method",
        y="memory_mb",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        palette=method_colors,
    )
    add_value_labels(ax, "%.1f")
    plt.title("Segmentation Method Memory Usage", fontweight="bold")
    plt.ylabel("Memory Usage (MB)")
    plt.xlabel("Segmentation Method")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "memory_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Create count comparison plots
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)
    fig.suptitle("Cell and Nuclei Counts by Method", fontweight="bold", y=0.98)

    ax1 = axes[0, 0]
    sns.barplot(
        x="method",
        y="initial_nuclei",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        ax=ax1,
        palette=method_colors,
    )
    add_value_labels(ax1, "%.0f")
    ax1.set_title("Initial Nuclei Count", fontweight="bold")
    ax1.set_ylabel("Count")
    ax1.set_xlabel("")

    ax2 = axes[0, 1]
    sns.barplot(
        x="method",
        y="initial_cells",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        ax=ax2,
        palette=method_colors,
    )
    add_value_labels(ax2, "%.0f")
    ax2.set_title("Initial Cell Count", fontweight="bold")
    ax2.set_ylabel("Count")
    ax2.set_xlabel("")

    ax3 = axes[1, 0]
    sns.barplot(
        x="method",
        y="final_cells",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        ax=ax3,
        palette=method_colors,
    )
    add_value_labels(ax3, "%.0f")
    ax3.set_title("Final Cell Count", fontweight="bold")
    ax3.set_ylabel("Count")
    ax3.set_xlabel("Method")

    ax4 = axes[1, 1]
    sns.barplot(
        x="method",
        y="runtime_seconds",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        ax=ax4,
        palette=method_colors,
    )
    add_value_labels(ax4, "%.0f")
    ax4.set_title("Segmentation Method Runtime Comparison", fontweight="bold")
    ax4.set_ylabel("Runtime (seconds)")
    ax4.set_xlabel("Method")

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "count_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Create retention rate plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Retention Rates by Method", fontweight="bold", y=0.98)

    sns.barplot(
        x="method",
        y="nuclei_retention",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        ax=ax1,
        palette=method_colors,
    )
    add_value_labels(ax1, "%.1f%%")
    ax1.set_title("Nuclei Retention Rate")
    ax1.set_ylabel("Retention Rate (%)")
    ax1.set_xlabel("Method")
    ax1.set_ylim(0, 105)

    sns.barplot(
        x="method",
        y="cell_retention",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        ax=ax2,
        palette=method_colors,
    )
    add_value_labels(ax2, "%.1f%%")
    ax2.set_title("Cell Retention Rate")
    ax2.set_ylabel("Retention Rate (%)")
    ax2.set_xlabel("Method")
    ax2.set_ylim(0, 105)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "retention_rates.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Create scatter plot comparing initial nuclei vs. cell counts
    plt.figure(figsize=(8, 6))
    scatter = sns.scatterplot(
        data=results_df,
        x="initial_nuclei",
        y="initial_cells",
        hue="method",
        style="method",
        s=80,
        alpha=0.8,
        palette=method_colors,
    )

    # Add reference line (y=x)
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    lims = [max(xmin, ymin), min(xmax, ymax)]
    plt.plot(lims, lims, "k--", alpha=0.5, zorder=0)

    # Improve legend
    plt.legend(title="Method", frameon=True, framealpha=0.9, edgecolor="lightgray")

    plt.title("Initial Nuclei vs. Cell Counts", fontweight="bold")
    plt.xlabel("Initial Nuclei Count")
    plt.ylabel("Initial Cell Count")
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "initial_count_correlation.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    # Create retention by stage plots
    # Calculate the differences for the stacked bars
    results_df["edge_removed_nuclei"] = (
        results_df["initial_nuclei"] - results_df["after_edge_removal_nuclei"]
    )
    results_df["edge_removed_cells"] = (
        results_df["initial_cells"] - results_df["after_edge_removal_cells"]
    )
    results_df["reconciled_nuclei"] = (
        results_df["after_edge_removal_nuclei"] - results_df["final_nuclei"]
    )
    results_df["reconciled_cells"] = (
        results_df["after_edge_removal_cells"] - results_df["final_cells"]
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Processing Pipeline by Method", fontweight="bold", y=0.98)

    # Use viridis colormap for stacked bars
    stack_colors = plt.cm.viridis(np.linspace(0.1, 0.9, 3))

    # Nuclei pipeline with cleaner appearance
    data = results_df.groupby("method")[
        ["final_nuclei", "reconciled_nuclei", "edge_removed_nuclei"]
    ].mean()
    data = data.reindex(sorted(data.index))  # Sort alphabetically
    data.plot(kind="bar", stacked=True, ax=ax1, alpha=0.8, color=stack_colors)
    ax1.set_title("Nuclei Processing Pipeline")
    ax1.set_ylabel("Count")
    ax1.set_xlabel("Method")
    ax1.legend(
        ["Final Nuclei", "Lost in Reconciliation", "Edge Removed"],
        frameon=True,
        framealpha=0.9,
        edgecolor="lightgray",
    )

    # Cells pipeline with cleaner appearance
    data = results_df.groupby("method")[
        ["final_cells", "reconciled_cells", "edge_removed_cells"]
    ].mean()
    data = data.reindex(sorted(data.index))  # Sort alphabetically
    data.plot(kind="bar", stacked=True, ax=ax2, alpha=0.8, color=stack_colors)
    ax2.set_title("Cell Processing Pipeline")
    ax2.set_ylabel("Count")
    ax2.set_xlabel("Method")
    ax2.legend(
        ["Final Cells", "Lost in Reconciliation", "Edge Removed"],
        frameon=True,
        framealpha=0.9,
        edgecolor="lightgray",
    )

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "processing_pipeline.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Calculate edge removal and reconciliation efficiency by method
    # Edge removal efficiency (percentage of objects kept after edge removal)
    results_df["edge_removal_nuclei_efficiency"] = np.where(
        results_df["initial_nuclei"] > 0,
        (results_df["after_edge_removal_nuclei"] / results_df["initial_nuclei"]) * 100,
        0,
    )

    results_df["edge_removal_cells_efficiency"] = np.where(
        results_df["initial_cells"] > 0,
        (results_df["after_edge_removal_cells"] / results_df["initial_cells"]) * 100,
        0,
    )

    # Reconciliation efficiency (percentage of edge-removed objects kept after reconciliation)
    results_df["reconciliation_nuclei_efficiency"] = np.where(
        results_df["after_edge_removal_nuclei"] > 0,
        (results_df["final_nuclei"] / results_df["after_edge_removal_nuclei"]) * 100,
        0,
    )

    results_df["reconciliation_cells_efficiency"] = np.where(
        results_df["after_edge_removal_cells"] > 0,
        (results_df["final_cells"] / results_df["after_edge_removal_cells"]) * 100,
        0,
    )

    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    fig.suptitle("Pipeline Efficiency by Method", fontweight="bold", y=0.98)

    ax1 = axes[0, 0]
    sns.barplot(
        x="method",
        y="edge_removal_nuclei_efficiency",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        ax=ax1,
        palette=method_colors,
    )
    add_value_labels(ax1, "%.1f%%")
    ax1.set_title("Edge Removal Efficiency - Nuclei", fontweight="bold")
    ax1.set_ylabel("Efficiency (%)")
    ax1.set_xlabel("")
    ax1.set_ylim(0, 105)

    ax2 = axes[0, 1]
    sns.barplot(
        x="method",
        y="edge_removal_cells_efficiency",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        ax=ax2,
        palette=method_colors,
    )
    add_value_labels(ax2, "%.1f%%")
    ax2.set_title("Edge Removal Efficiency - Cells", fontweight="bold")
    ax2.set_ylabel("Efficiency (%)")
    ax2.set_xlabel("")
    ax2.set_ylim(0, 105)

    ax3 = axes[1, 0]
    sns.barplot(
        x="method",
        y="reconciliation_nuclei_efficiency",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        ax=ax3,
        palette=method_colors,
    )
    add_value_labels(ax3, "%.1f%%")
    ax3.set_title("Reconciliation Efficiency - Nuclei", fontweight="bold")
    ax3.set_ylabel("Efficiency (%)")
    ax3.set_xlabel("Method")
    ax3.set_ylim(0, 105)

    ax4 = axes[1, 1]
    sns.barplot(
        x="method",
        y="reconciliation_cells_efficiency",
        data=results_df,
        errorbar=None,
        alpha=0.8,
        ax=ax4,
        palette=method_colors,
    )
    add_value_labels(ax4, "%.1f%%")
    ax4.set_title("Reconciliation Efficiency - Cells", fontweight="bold")
    ax4.set_ylabel("Efficiency (%)")
    ax4.set_xlabel("Method")
    ax4.set_ylim(0, 105)

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "pipeline_efficiency.png", dpi=300, bbox_inches="tight")
    plt.close()

    # Print key findings to console
    print("\n=== SEGMENTATION BENCHMARK SUMMARY ===")
    for method in results_df["method"].unique():
        method_df = results_df[results_df["method"] == method]
        print(f"\n--- {method.upper()} ---")
        print(f"Average runtime: {method_df['runtime_seconds'].mean():.2f} seconds")
        print(f"Average memory usage: {method_df['memory_mb'].mean():.1f} MB")
        print(
            f"Initial counts: {method_df['initial_nuclei'].mean():.1f} nuclei, {method_df['initial_cells'].mean():.1f} cells"
        )
        print(
            f"After edge removal: {method_df['after_edge_removal_nuclei'].mean():.1f} nuclei, {method_df['after_edge_removal_cells'].mean():.1f} cells"
        )
        print(
            f"Final counts: {method_df['final_nuclei'].mean():.1f} nuclei, {method_df['final_cells'].mean():.1f} cells"
        )
        print(
            f"Retention rates: {method_df['nuclei_retention'].mean():.1f}% nuclei, {method_df['cell_retention'].mean():.1f}% cells"
        )
        print(
            f"Edge removal efficiency: {method_df['edge_removal_nuclei_efficiency'].mean():.1f}% nuclei, {method_df['edge_removal_cells_efficiency'].mean():.1f}% cells"
        )
        print(
            f"Reconciliation efficiency: {method_df['reconciliation_nuclei_efficiency'].mean():.1f}% nuclei, {method_df['reconciliation_cells_efficiency'].mean():.1f}% cells"
        )


if __name__ == "__main__":
    if args.plots_only:
        # Load existing results and regenerate plots
        results_path = OUTPUT_DIR / "benchmark_results.csv"
        if results_path.exists():
            print("Regenerating plots from existing results...")
            results_df = pd.read_csv(results_path)
            generate_benchmark_summary(results_df)
        else:
            print(f"No results found at {results_path}. Run benchmark first.")
    else:
        benchmark_segmentation_methods()
