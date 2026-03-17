"""
Spot Calling Benchmark

Compares two spot detection methods for SBS (sequencing-by-synthesis):
- Standard: Classical peak detection (LOG filter -> std dev -> find peaks -> threshold)
- Spotiflow: Deep learning-based spot detection

Metrics: runtime, memory, peaks detected, read mapping rates, cell-level assignment.
See BENCHMARKING_SUMMARY.md for detailed documentation.
"""

import os
import sys
import time
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tifffile import imread, imwrite
import shutil
from pathlib import Path
from datetime import datetime
import psutil
import gc
from tqdm import tqdm

# Add brieflow to path (relative to benchmarks directory)
BENCHMARKS_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BENCHMARKS_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT / "brieflow" / "workflow"))

# Import spot calling functions and utilities
from lib.sbs.compute_standard_deviation import compute_standard_deviation
from lib.sbs.find_peaks import find_peaks, find_peaks_spotiflow
from lib.sbs.extract_bases import extract_bases
from lib.sbs.call_reads import call_reads
from lib.sbs.call_cells import call_cells
from lib.shared.extract_phenotype_minimal import extract_phenotype_minimal

# Import evaluation functions
from lib.sbs.eval_mapping import plot_read_mapping_heatmap, mapping_overview

random.seed(42)

# Configure paths (relative to project root)
PIPELINE_OUTPUT_DIR = PROJECT_ROOT / "analysis" / "brieflow_output" / "sbs" / "images"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "spot_calling"
BARCODE_LIBRARY_FILE = PROJECT_ROOT / "analysis" / "config" / "barcode_library.tsv"

# Create directories
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Create subdirectories for each spot calling method
for method in ["standard", "spotiflow"]:
    (OUTPUT_DIR / method).mkdir(exist_ok=True)

# Create directories for visualization
(OUTPUT_DIR / "aligned").mkdir(exist_ok=True)
(OUTPUT_DIR / "max_filtered").mkdir(exist_ok=True)

# Parameters for spot calling methods (from config.yml)
STANDARD_PARAMS = {
    "name": "standard",
    "dapi_index": 0,
    "bases": ["G", "T", "A", "C"],
    "threshold_peaks": 400,
    "call_reads_method": "median",
    "error_correct": False,
    "q_min": 0,
    "prefix_col": "prefix",
    "sort_calls": "count",
}

SPOTIFLOW_PARAMS = {
    "name": "spotiflow",
    "dapi_index": 0,
    "spotiflow_model": "general",
    "spotiflow_threshold": 0.3,
    "spotiflow_cycle_index": 0,
    "spotiflow_min_distance": 1,
    "bases": ["G", "T", "A", "C"],
    "threshold_peaks": 0,
    "call_reads_method": "median",
    "error_correct": False,
    "q_min": 0,
    "prefix_col": "prefix",
    "sort_calls": "count",
}


def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024


def find_sample_data(one_per_plate=True):
    """Find sample SBS data with complete permanent intermediate files.

    Args:
        one_per_plate: If True, select one random sample per plate (default).
                      This gives 8 total samples for 8 plates.
    """
    # Look for samples with complete permanent outputs
    aligned_files = list(PIPELINE_OUTPUT_DIR.glob("P-*_W-*_T-*__aligned.tiff"))

    if not aligned_files:
        print(
            "No aligned files found! Have you run the pipeline with temp→None changes?"
        )
        return []

    # Group samples by plate
    samples_by_plate = {}

    for aligned_file in aligned_files:
        base_filename = aligned_file.name.replace("__aligned.tiff", "")

        # Parse plate, well, tile from filename
        parts = base_filename.split("_")
        plate = parts[0].replace("P-", "")
        well = parts[1].replace("W-", "")
        tile = parts[2].replace("T-", "")

        # Check for all required permanent files
        log_filtered_file = PIPELINE_OUTPUT_DIR / f"{base_filename}__log_filtered.tiff"
        max_filtered_file = PIPELINE_OUTPUT_DIR / f"{base_filename}__max_filtered.tiff"
        nuclei_file = PIPELINE_OUTPUT_DIR / f"{base_filename}__nuclei.tiff"
        cells_file = PIPELINE_OUTPUT_DIR / f"{base_filename}__cells.tiff"

        # Only include if ALL permanent files exist
        if all(
            [
                log_filtered_file.exists(),
                max_filtered_file.exists(),
                nuclei_file.exists(),
                cells_file.exists(),
            ]
        ):
            sample = {
                "plate": plate,
                "well": well,
                "tile": tile,
                "base_name": base_filename,
                "aligned_file": str(aligned_file),
                "log_filtered_file": str(log_filtered_file),
                "max_filtered_file": str(max_filtered_file),
                "nuclei_file": str(nuclei_file),
                "cells_file": str(cells_file),
            }

            if plate not in samples_by_plate:
                samples_by_plate[plate] = []
            samples_by_plate[plate].append(sample)

    # Select one random sample per plate
    selected_samples = []
    for plate in sorted(samples_by_plate.keys()):
        samples = samples_by_plate[plate]
        if samples:
            selected = random.choice(samples)
            selected_samples.append(selected)
            print(
                f"Selected sample: plate {plate}, well {selected['well']}, tile {selected['tile']}"
            )

    if not selected_samples:
        print("Could not find samples with complete permanent files!")
        print("Make sure pipeline has been run with temp→None changes.")
    else:
        print(
            f"\nFound {len(selected_samples)} samples (one per plate) with complete permanent files"
        )

    return selected_samples


def process_standard_method(sample, df_barcode_library=None):
    """Process a sample with the standard spot detection method."""
    print(
        f"Processing plate {sample['plate']}, well {sample['well']}, tile {sample['tile']} with standard method..."
    )

    # Track metrics
    metrics = {
        "method": "standard",
        "plate": sample["plate"],
        "well": sample["well"],
        "tile": sample["tile"],
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    try:
        # Load pre-existing permanent files
        print("  Loading intermediate files from pipeline...")

        aligned_images = imread(sample["aligned_file"])
        log_filtered_data = imread(sample["log_filtered_file"])
        max_filtered_data = imread(sample["max_filtered_file"])
        cells_data = imread(sample["cells_file"])
        nuclei_data = imread(sample["nuclei_file"])

        print(f"  Loaded aligned images: {aligned_images.shape}")
        print(f"  Loaded LOG filtered data: {log_filtered_data.shape}")
        print(f"  Loaded max filtered data: {max_filtered_data.shape}")

        # Start timing
        gc.collect()
        start_mem = get_memory_usage()
        start_time = time.time()

        # Compute standard deviation (method-specific for standard peak detection)
        standard_deviation_data = compute_standard_deviation(
            log_filtered_data=log_filtered_data,
            remove_index=STANDARD_PARAMS["dapi_index"],
        )

        # Find peaks using standard method
        peaks_data = find_peaks(standard_deviation_data=standard_deviation_data)

        # Save peaks TIFF for visualization
        peaks_path = OUTPUT_DIR / "standard" / f"{sample['base_name']}_peaks.tiff"
        imwrite(peaks_path, peaks_data)

        # Compute runtime and memory usage
        runtime = time.time() - start_time
        end_mem = get_memory_usage()
        mem_usage = max(0, end_mem - start_mem)

        print(f"  Standard peak detection completed in {runtime:.2f} seconds")

        # Continue with downstream analysis
        # Extract bases at peak locations
        bases_data = extract_bases(
            peaks_data=peaks_data,
            max_filtered_data=max_filtered_data,
            cells_data=cells_data,
            threshold_peaks=STANDARD_PARAMS["threshold_peaks"],
            bases=STANDARD_PARAMS["bases"],
            wildcards={
                "plate": sample["plate"],
                "well": sample["well"],
                "tile": sample["tile"],
            },
        )

        # Save bases data
        bases_path = OUTPUT_DIR / "standard" / f"{sample['base_name']}_bases.tsv"
        bases_data.to_csv(bases_path, sep="\t", index=False)

        # Call reads from bases
        reads_data = call_reads(
            bases_data=bases_data,
            peaks_data=peaks_data,
            method=STANDARD_PARAMS["call_reads_method"],
        )

        # Save reads data
        reads_path = OUTPUT_DIR / "standard" / f"{sample['base_name']}_reads.tsv"
        reads_data.to_csv(reads_path, sep="\t", index=False)

        # Add metrics
        metrics["runtime_seconds"] = runtime
        metrics["memory_mb"] = mem_usage
        metrics["total_peaks"] = len(peaks_data.nonzero()[0])
        metrics["total_bases"] = len(bases_data)
        metrics["total_reads"] = len(reads_data)

        # If barcode library is provided, evaluate mapping quality
        if df_barcode_library is not None:
            # Get barcodes for mapping evaluation
            barcodes = df_barcode_library["prefix"].tolist()

            # Calculate read mapping rate
            read_mapping = plot_read_mapping_heatmap(
                reads_data,
                barcodes,
                shape="6W_sbs",
                return_plot=False,
                return_summary=True,
            )

            if read_mapping is not None and len(read_mapping) > 0:
                metrics["read_mapping_fraction"] = read_mapping[
                    "fraction of reads mapping"
                ].mean()

            # Call cells (assign barcodes to cells)
            cells_data_df = call_cells(
                reads_data=reads_data,
                df_barcode_library=df_barcode_library,
                q_min=STANDARD_PARAMS["q_min"],
                prefix_col=STANDARD_PARAMS["prefix_col"],
                sort_calls=STANDARD_PARAMS["sort_calls"],
                error_correct=STANDARD_PARAMS["error_correct"],
            )

            # Save cells data
            cells_path = OUTPUT_DIR / "standard" / f"{sample['base_name']}_cells.tsv"
            cells_data_df.to_csv(cells_path, sep="\t", index=False)

            # Extract minimal phenotype info
            sbs_info = extract_phenotype_minimal(
                phenotype_data=nuclei_data,
                nuclei_data=nuclei_data,
                wildcards={
                    "plate": sample["plate"],
                    "well": sample["well"],
                    "tile": sample["tile"],
                },
            )

            # Save SBS info
            sbs_info_path = (
                OUTPUT_DIR / "standard" / f"{sample['base_name']}_sbs_info.tsv"
            )
            sbs_info.to_csv(sbs_info_path, sep="\t", index=False)

            # Calculate cell mapping statistics
            mapping_stats = mapping_overview(sbs_info, cells_data_df)
            for col in mapping_stats.columns:
                if col != "well":
                    metrics[col] = mapping_stats.iloc[0][col]

        print(f"  Successfully processed standard method in {runtime:.2f} seconds")
        return {"metrics": metrics}

    except Exception as e:
        print(f"  Error in standard method: {e}")
        import traceback

        traceback.print_exc()
        metrics["error"] = str(e)
        return {"metrics": metrics}


def process_spotiflow_method(sample, df_barcode_library=None):
    """Process a sample with the Spotiflow spot detection method."""
    print(
        f"Processing plate {sample['plate']}, well {sample['well']}, tile {sample['tile']} with spotiflow method..."
    )

    # Track metrics
    metrics = {
        "method": "spotiflow",
        "plate": sample["plate"],
        "well": sample["well"],
        "tile": sample["tile"],
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }

    try:
        # Load pre-existing permanent files
        print("  Loading intermediate files from pipeline...")

        aligned_images = imread(sample["aligned_file"])
        max_filtered_data = imread(sample["max_filtered_file"])
        cells_data = imread(sample["cells_file"])
        nuclei_data = imread(sample["nuclei_file"])

        print(f"  Loaded aligned images: {aligned_images.shape}")
        print(f"  Loaded max filtered data: {max_filtered_data.shape}")

        # Start timing
        gc.collect()
        start_mem = get_memory_usage()
        start_time = time.time()

        # Find peaks using Spotiflow
        peaks_data, _ = find_peaks_spotiflow(
            aligned_images=aligned_images,
            cycle_idx=SPOTIFLOW_PARAMS["spotiflow_cycle_index"],
            model=SPOTIFLOW_PARAMS["spotiflow_model"],
            prob_thresh=SPOTIFLOW_PARAMS["spotiflow_threshold"],
            min_distance=SPOTIFLOW_PARAMS["spotiflow_min_distance"],
            remove_index=SPOTIFLOW_PARAMS["dapi_index"],
            verbose=True,
        )

        # Save peaks TIFF for visualization
        peaks_path = OUTPUT_DIR / "spotiflow" / f"{sample['base_name']}_peaks.tiff"
        imwrite(peaks_path, peaks_data)

        # Compute runtime and memory usage
        runtime = time.time() - start_time
        end_mem = get_memory_usage()
        mem_usage = max(0, end_mem - start_mem)

        print(f"  Spotiflow peak detection completed in {runtime:.2f} seconds")

        # Continue with downstream analysis
        # Extract bases at peak locations
        bases_data = extract_bases(
            peaks_data=peaks_data,
            max_filtered_data=max_filtered_data,
            cells_data=cells_data,
            threshold_peaks=SPOTIFLOW_PARAMS["threshold_peaks"],
            bases=SPOTIFLOW_PARAMS["bases"],
            wildcards={
                "plate": sample["plate"],
                "well": sample["well"],
                "tile": sample["tile"],
            },
        )

        # Save bases data
        bases_path = OUTPUT_DIR / "spotiflow" / f"{sample['base_name']}_bases.tsv"
        bases_data.to_csv(bases_path, sep="\t", index=False)

        # Call reads from bases
        reads_data = call_reads(
            bases_data=bases_data,
            peaks_data=peaks_data,
            method=SPOTIFLOW_PARAMS["call_reads_method"],
        )

        # Save reads data
        reads_path = OUTPUT_DIR / "spotiflow" / f"{sample['base_name']}_reads.tsv"
        reads_data.to_csv(reads_path, sep="\t", index=False)

        # Add metrics
        metrics["runtime_seconds"] = runtime
        metrics["memory_mb"] = mem_usage
        metrics["total_peaks"] = len(peaks_data.nonzero()[0])
        metrics["total_bases"] = len(bases_data)
        metrics["total_reads"] = len(reads_data)

        # If barcode library is provided, evaluate mapping quality
        if df_barcode_library is not None:
            # Get barcodes for mapping evaluation
            barcodes = df_barcode_library["prefix"].tolist()

            # Calculate read mapping rate
            read_mapping = plot_read_mapping_heatmap(
                reads_data,
                barcodes,
                shape="6W_sbs",
                return_plot=False,
                return_summary=True,
            )

            if read_mapping is not None and len(read_mapping) > 0:
                metrics["read_mapping_fraction"] = read_mapping[
                    "fraction of reads mapping"
                ].mean()

            # Call cells (assign barcodes to cells)
            cells_data_df = call_cells(
                reads_data=reads_data,
                df_barcode_library=df_barcode_library,
                q_min=SPOTIFLOW_PARAMS["q_min"],
                prefix_col=SPOTIFLOW_PARAMS["prefix_col"],
                sort_calls=SPOTIFLOW_PARAMS["sort_calls"],
                error_correct=SPOTIFLOW_PARAMS["error_correct"],
            )

            # Save cells data
            cells_path = OUTPUT_DIR / "spotiflow" / f"{sample['base_name']}_cells.tsv"
            cells_data_df.to_csv(cells_path, sep="\t", index=False)

            # Extract minimal phenotype info
            sbs_info = extract_phenotype_minimal(
                phenotype_data=nuclei_data,
                nuclei_data=nuclei_data,
                wildcards={
                    "plate": sample["plate"],
                    "well": sample["well"],
                    "tile": sample["tile"],
                },
            )

            # Save SBS info
            sbs_info_path = (
                OUTPUT_DIR / "spotiflow" / f"{sample['base_name']}_sbs_info.tsv"
            )
            sbs_info.to_csv(sbs_info_path, sep="\t", index=False)

            # Calculate cell mapping statistics
            mapping_stats = mapping_overview(sbs_info, cells_data_df)
            for col in mapping_stats.columns:
                if col != "well":
                    metrics[col] = mapping_stats.iloc[0][col]

        print(f"  Successfully processed spotiflow method in {runtime:.2f} seconds")
        return {"metrics": metrics}

    except Exception as e:
        print(f"  Error in spotiflow method: {e}")
        import traceback

        traceback.print_exc()
        metrics["error"] = str(e)
        return {"metrics": metrics}


def run_benchmark():
    """Main benchmark workflow comparing standard vs spotiflow spot detection."""
    print("=" * 70)
    print("SPOT CALLING BENCHMARK")
    print("=" * 70)

    # Find sample data with complete permanent files (one per plate = 8 samples)
    samples = find_sample_data(one_per_plate=True)

    if not samples:
        print("\nERROR: Could not find suitable sample data for benchmarking")
        print("Make sure the pipeline has been run with temp→None changes.")
        return None

    print(f"\nFound {len(samples)} samples for benchmarking\n")

    # Initialize results collection
    all_metrics = []

    # Load barcode library if available
    df_barcode_library = None
    try:
        if BARCODE_LIBRARY_FILE.exists():
            df_barcode_library = pd.read_csv(BARCODE_LIBRARY_FILE, sep="\t")
            print(f"Loaded barcode library with {len(df_barcode_library)} entries\n")
    except Exception as e:
        print(f"Warning: Error loading barcode library: {e}")
        print("Will continue without evaluating mapping quality\n")

    # Process each sample with each method
    for sample in tqdm(samples, desc="Processing samples"):
        sample_id = f"P-{sample['plate']}_W-{sample['well']}_T-{sample['tile']}"

        # Copy aligned image to results directory for visualization
        aligned_dst = OUTPUT_DIR / "aligned" / f"{sample['base_name']}__aligned.tiff"
        if not aligned_dst.exists():
            shutil.copy(sample["aligned_file"], aligned_dst)

        # Copy max_filtered image to results directory for visualization
        max_filtered_dst = (
            OUTPUT_DIR / "max_filtered" / f"{sample['base_name']}__max_filtered.tiff"
        )
        if not max_filtered_dst.exists():
            shutil.copy(sample["max_filtered_file"], max_filtered_dst)

        # Process with standard method
        print(f"\n{'=' * 60}")
        print(f"Sample: {sample_id} - Standard Method")
        print(f"{'=' * 60}")
        standard_results = process_standard_method(sample, df_barcode_library)
        all_metrics.append(standard_results["metrics"])

        # Process with spotiflow method
        print(f"\n{'=' * 60}")
        print(f"Sample: {sample_id} - Spotiflow Method")
        print(f"{'=' * 60}")
        spotiflow_results = process_spotiflow_method(sample, df_barcode_library)
        all_metrics.append(spotiflow_results["metrics"])

    # Combine all metrics into DataFrame
    results_df = pd.DataFrame(all_metrics)

    # Save results
    results_df.to_csv(OUTPUT_DIR / "benchmark_results.csv", index=False)
    print(f"\n{'=' * 60}")
    print(f"Results saved to {OUTPUT_DIR / 'benchmark_results.csv'}")
    print(f"{'=' * 60}")

    # Print summary
    print("\n=== BENCHMARK SUMMARY ===")
    for method in results_df["method"].unique():
        method_data = results_df[results_df["method"] == method]
        print(f"\n--- {method.upper()} ---")
        print(
            f"Runtime: {method_data['runtime_seconds'].mean():.2f} sec (±{method_data['runtime_seconds'].std():.2f})"
        )
        print(
            f"Memory: {method_data['memory_mb'].mean():.1f} MB (±{method_data['memory_mb'].std():.1f})"
        )
        print(
            f"Peaks: {method_data['total_peaks'].mean():.0f} (±{method_data['total_peaks'].std():.0f})"
        )
        print(
            f"Reads: {method_data['total_reads'].mean():.0f} (±{method_data['total_reads'].std():.0f})"
        )

    # Generate visualization
    print("\n=== GENERATING VISUALIZATION ===")
    generate_visualization(results_df)

    return results_df


def generate_visualization(results_df=None):
    """Generate publication-quality visualization comparing Spotiflow vs Standard methods.

    Args:
        results_df: DataFrame with benchmark results. If None, loads from benchmark_results.csv
    """
    from plot_style import setup_plot_style, box_strip, FIGSIZE, COLORS, save_figure, print_summary_table

    # Apply consistent plot styling
    setup_plot_style()

    # Load raw per-tile results
    if results_df is None:
        raw_path = OUTPUT_DIR / "benchmark_results.csv"
        if raw_path.exists():
            results_df = pd.read_csv(raw_path)
        else:
            print("No benchmark_results.csv found. Run benchmark first.")
            return

    # Also save the summary for reference
    numeric_cols = results_df.select_dtypes(include=[np.number]).columns.tolist()
    method_summary = results_df.groupby("method")[numeric_cols].agg(["mean", "std"])
    method_summary.to_csv(OUTPUT_DIR / "method_summary.csv")
    # Display labels and palette
    label_map = {"spotiflow": "Spotiflow", "standard": "Standard"}
    results_df["Method"] = results_df["method"].map(label_map)
    method_order = ["Spotiflow", "Standard"]
    palette = {
        "Spotiflow": COLORS.get("spotiflow", "#9467bd"),
        "Standard": COLORS.get("standard", "#d62728"),
    }

    # Runtime & Memory comparison
    fig, (ax_rt, ax_mem) = plt.subplots(1, 2, figsize=FIGSIZE["double"])
    box_strip(ax_rt, results_df, "Method", "runtime_seconds", palette, method_order,
              ylabel="Runtime (seconds)", title="Runtime per Tile")
    box_strip(ax_mem, results_df, "Method", "memory_mb", palette, method_order,
              ylabel="Memory (MB)", title="Peak Memory per Tile")
    plt.tight_layout()
    save_figure(fig, OUTPUT_DIR / "runtime_memory_comparison.png")
    plt.close()
    print(f"Saved: {OUTPUT_DIR / 'runtime_memory_comparison.png'}")

    # Spots + mapping rate comparison (square layout, replaces old bar chart)
    if "total_reads" in results_df.columns:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 6))
        box_strip(ax1, results_df, "Method", "total_reads", palette, method_order,
                  ylabel="Number of Spots", title="Filtered Spots",
                  fmt="int")
        if "1_or_more_genes__percent" in results_df.columns:
            box_strip(ax2, results_df, "Method", "1_or_more_genes__percent",
                      palette, method_order,
                      ylabel="Percentage of Cells (%)", title="Mapped Cells",
                      fmt="pct")
        plt.tight_layout()
        save_figure(fig, OUTPUT_DIR / "spots_mapping_comparison.png")
        plt.close()
        print(f"Saved: {OUTPUT_DIR / 'spots_mapping_comparison.png'}")

        # Square version for publication (larger labels to match original)
        fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=FIGSIZE["square"])
        box_strip(ax3, results_df, "Method", "total_reads", palette, method_order,
                  ylabel="Number of Spots", title="Filtered Spots",
                  fmt="int")
        if "1_or_more_genes__percent" in results_df.columns:
            box_strip(ax4, results_df, "Method", "1_or_more_genes__percent",
                      palette, method_order,
                      ylabel="Percentage of Cells (%)", title="Mapped Cells",
                      fmt="pct")
        for ax in (ax3, ax4):
            ax.set_title(ax.get_title(), fontsize=20, fontweight="bold")
            ax.set_ylabel(ax.get_ylabel(), fontsize=16)
            ax.tick_params(labelsize=14)
        plt.tight_layout()
        save_figure(fig2, OUTPUT_DIR / "spotiflow_vs_standard_comparison_square.png",
                    transparent=True)
        plt.close()
        print(f"Saved: {OUTPUT_DIR / 'spotiflow_vs_standard_comparison_square.png'}")

    # Print publication-ready summary tables
    pct_fmt = lambda m, lo, hi: f"{m:.1f}% [{lo:.1f}–{hi:.1f}%]"
    spot_metrics = [
        ("runtime_seconds", "Runtime (s)"),
        ("memory_mb", "Memory (MB)"),
        ("total_reads", "Filtered Spots"),
    ]
    spot_fmt = {}
    if "1_or_more_genes__percent" in results_df.columns:
        spot_metrics.append(("1_or_more_genes__percent", "Mapped Cells (%)"))
        spot_fmt["1_or_more_genes__percent"] = pct_fmt
    print_summary_table(
        "Spot Calling — Spotiflow vs Standard",
        results_df, "Method", spot_metrics, fmt_map=spot_fmt,
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Spot calling benchmark")
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Only regenerate plots from existing results",
    )
    args = parser.parse_args()

    if args.plots_only:
        # Load existing results and regenerate plots
        results_path = OUTPUT_DIR / "benchmark_results.csv"
        if results_path.exists():
            print("Regenerating plots from existing results...")
            results_df = pd.read_csv(results_path)
            generate_visualization(results_df)
        else:
            print(f"No results found at {results_path}. Run benchmark first.")
    else:
        # Run full benchmark
        run_benchmark()
