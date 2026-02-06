"""
Merge Benchmark

Compares merge approaches from actual Snakemake pipeline runs:
- Fast: Tile-by-tile alignment using initial_sites configuration (brieflow_output/merge/)
- Stitch: Full well stitching with triangle hash alignment (brieflow_output/merge_stitch/)

Extracts metrics from slurm logs and parquet outputs for well P-1_W-A1.

Directory structure:
    brieflow_output/
        merge/           # Fast approach (default, all wells)
            parquets/    # P-{plate}_W-{well}__fast_merge.parquet, __merge_final.parquet
            eval/        # Evaluation outputs
        merge_stitch/    # Stitch approach (test on P-1_W-A1)
            parquets/    # P-{plate}_W-{well}__merged_cells.parquet, __merge_final.parquet
            eval/        # Evaluation outputs

Usage:
    python merge.py              # Generate all outputs
    python merge.py --plots-only # Regenerate plots only
"""

import re
import argparse
import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

# Import shared plot styling
from plot_style import setup_plot_style, FIGSIZE, COLORS, add_value_labels, save_figure

# Configure paths
BENCHMARKS_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BENCHMARKS_DIR.parent
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
SLURM_DIR = ANALYSIS_DIR / "slurm" / "slurm_output" / "rule"
SLURM_CONFIG = ANALYSIS_DIR / "slurm" / "config.yaml"
BRIEFLOW_OUTPUT = ANALYSIS_DIR / "brieflow_output"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "merge"

# Output directories for each approach
FAST_OUTPUT = BRIEFLOW_OUTPUT / "merge"
STITCH_OUTPUT = BRIEFLOW_OUTPUT / "merge_stitch"

# Well to analyze
PLATE = 1
WELL = "A1"

# Slurm job IDs for each approach (P-1_W-A1)
FAST_JOBS = {
    "fast_alignment": "5762209",
    "fast_merge": "5763763",
    "format_merge": "5763817",
    "deduplicate_merge": "5763917",
    "final_merge": "5764197",
}

# Stitch jobs from Jan 7 run (merge_stitch output)
STITCH_JOBS = {
    "estimate_stitch_sbs": "5798715",
    "estimate_stitch_phenotype": "5798716",
    "stitch_sbs": "5798717",
    "stitch_phenotype": "5798718",
    "stitch_alignment": "5812454",
    "stitch_merge": "5817135",
    "format_merge": "5841527",
    "deduplicate_merge": "5841552",
    "final_merge": "5842342",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Merge benchmark from snakemake run")
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Only regenerate plots from existing results",
    )
    return parser.parse_args()


def parse_slurm_log(log_path):
    """Extract timing and resource info from a slurm log file."""
    if not log_path.exists():
        return None

    with open(log_path, "r") as f:
        content = f.read()

    # Extract timestamps
    timestamps = re.findall(
        r"\[([A-Za-z]{3} [A-Za-z]{3} +\d+ \d+:\d+:\d+ \d+)\]", content
    )

    # Extract memory requested
    mem_match = re.search(r"mem_mb=(\d+)", content)
    mem_requested = int(mem_match.group(1)) if mem_match else None

    # Extract rule name
    rule_match = re.search(r"(?:local)?rule (\w+):", content)
    rule_name = rule_match.group(1) if rule_match else None

    # Calculate runtime from first and last timestamp
    runtime_seconds = None
    if len(timestamps) >= 2:
        try:
            start = datetime.strptime(timestamps[0], "%a %b %d %H:%M:%S %Y")
            end = datetime.strptime(timestamps[-1], "%a %b %d %H:%M:%S %Y")
            runtime_seconds = (end - start).total_seconds()
        except ValueError:
            pass

    # Check for errors
    has_error = "error" in content.lower() or "failed" in content.lower()
    identity_fallback = "identity transformation" in content.lower()

    return {
        "rule": rule_name,
        "runtime_seconds": runtime_seconds,
        "memory_requested_mb": mem_requested,
        "has_error": has_error,
        "identity_fallback": identity_fallback,
    }


def extract_timing_from_logs():
    """Extract timing data from all slurm logs."""
    results = []

    # Fast approach
    for step, job_id in FAST_JOBS.items():
        log_path = SLURM_DIR / f"{job_id}.out"
        info = parse_slurm_log(log_path)
        if info:
            results.append(
                {
                    "approach": "fast",
                    "step": step,
                    "runtime_seconds": info["runtime_seconds"],
                    "memory_requested_mb": info["memory_requested_mb"],
                    "status": "error" if info["has_error"] else "success",
                }
            )

    # Stitch approach
    for step, job_id in STITCH_JOBS.items():
        log_path = SLURM_DIR / f"{job_id}.out"
        info = parse_slurm_log(log_path)
        if info:
            status = "success"
            notes = ""
            if info["identity_fallback"]:
                status = "alignment_degraded"
                notes = "identity transform fallback due to memory"
            elif info["has_error"]:
                status = "error"

            results.append(
                {
                    "approach": "stitch",
                    "step": step,
                    "runtime_seconds": info["runtime_seconds"],
                    "memory_requested_mb": info["memory_requested_mb"],
                    "status": status,
                    "notes": notes,
                }
            )

    return pd.DataFrame(results)


def extract_quality_metrics():
    """Extract quality metrics from parquet files."""
    results = []

    # Get input cell counts (shared between approaches)
    ph_info_path = (
        BRIEFLOW_OUTPUT
        / "phenotype"
        / "parquets"
        / f"P-{PLATE}_W-{WELL}__phenotype_info.parquet"
    )
    sbs_info_path = (
        BRIEFLOW_OUTPUT / "sbs" / "parquets" / f"P-{PLATE}_W-{WELL}__sbs_info.parquet"
    )

    phenotype_cells = (
        len(pd.read_parquet(ph_info_path)) if ph_info_path.exists() else None
    )
    sbs_cells = len(pd.read_parquet(sbs_info_path)) if sbs_info_path.exists() else None

    # Fast approach
    fast_merge_path = (
        FAST_OUTPUT / "parquets" / f"P-{PLATE}_W-{WELL}__fast_merge.parquet"
    )
    fast_final_path = (
        FAST_OUTPUT / "parquets" / f"P-{PLATE}_W-{WELL}__merge_final.parquet"
    )

    if fast_merge_path.exists():
        fast_merge = pd.read_parquet(fast_merge_path)
        fast_final = (
            pd.read_parquet(fast_final_path) if fast_final_path.exists() else None
        )

        # Barcode metrics
        cells_with_barcode = (
            fast_final["gene_symbol_0"].notna().sum()
            if fast_final is not None and "gene_symbol_0" in fast_final.columns
            else None
        )
        cells_single_gene = (
            fast_final["mapped_single_gene"].sum()
            if fast_final is not None and "mapped_single_gene" in fast_final.columns
            else None
        )

        results.append(
            {
                "approach": "fast",
                "plate": PLATE,
                "well": WELL,
                "phenotype_cells": phenotype_cells,
                "sbs_cells": sbs_cells,
                "matched_cells": len(fast_merge),
                "final_cells": len(fast_final) if fast_final is not None else None,
                "cells_with_barcode": cells_with_barcode,
                "cells_single_gene": cells_single_gene,
                "mean_distance": fast_merge["distance"].mean(),
                "median_distance": fast_merge["distance"].median(),
                "std_distance": fast_merge["distance"].std(),
                "max_distance": fast_merge["distance"].max(),
                "matches_under_1px": (fast_merge["distance"] < 1).sum(),
                "matches_under_2px": (fast_merge["distance"] < 2).sum(),
                "matches_under_5px": (fast_merge["distance"] < 5).sum(),
                "matches_under_10px": (fast_merge["distance"] < 10).sum(),
                "pct_under_1px": (fast_merge["distance"] < 1).mean() * 100,
                "pct_under_2px": (fast_merge["distance"] < 2).mean() * 100,
                "barcode_rate": cells_with_barcode / len(fast_final) * 100
                if fast_final is not None and cells_with_barcode
                else None,
                "single_gene_rate": cells_single_gene / len(fast_final) * 100
                if fast_final is not None and cells_single_gene
                else None,
                "phenotype_retention": len(fast_final) / phenotype_cells * 100
                if phenotype_cells and fast_final is not None
                else None,
                "sbs_retention": len(fast_final) / sbs_cells * 100
                if sbs_cells and fast_final is not None
                else None,
            }
        )

    # Stitch approach
    stitch_merge_path = (
        STITCH_OUTPUT / "parquets" / f"P-{PLATE}_W-{WELL}__merged_cells.parquet"
    )
    stitch_final_path = (
        STITCH_OUTPUT / "parquets" / f"P-{PLATE}_W-{WELL}__merge_final.parquet"
    )

    if stitch_merge_path.exists():
        stitch_merge = pd.read_parquet(stitch_merge_path)
        stitch_final = (
            pd.read_parquet(stitch_final_path) if stitch_final_path.exists() else None
        )

        # Barcode metrics
        cells_with_barcode = (
            stitch_final["gene_symbol_0"].notna().sum()
            if stitch_final is not None and "gene_symbol_0" in stitch_final.columns
            else None
        )
        cells_single_gene = (
            stitch_final["mapped_single_gene"].sum()
            if stitch_final is not None and "mapped_single_gene" in stitch_final.columns
            else None
        )

        results.append(
            {
                "approach": "stitch",
                "plate": PLATE,
                "well": WELL,
                "phenotype_cells": phenotype_cells,
                "sbs_cells": sbs_cells,
                "matched_cells": len(stitch_merge),
                "final_cells": len(stitch_final) if stitch_final is not None else None,
                "cells_with_barcode": cells_with_barcode,
                "cells_single_gene": cells_single_gene,
                "mean_distance": stitch_merge["distance"].mean(),
                "median_distance": stitch_merge["distance"].median(),
                "std_distance": stitch_merge["distance"].std(),
                "max_distance": stitch_merge["distance"].max(),
                "matches_under_1px": (stitch_merge["distance"] < 1).sum(),
                "matches_under_2px": (stitch_merge["distance"] < 2).sum(),
                "matches_under_5px": (stitch_merge["distance"] < 5).sum(),
                "matches_under_10px": (stitch_merge["distance"] < 10).sum(),
                "pct_under_1px": (stitch_merge["distance"] < 1).mean() * 100,
                "pct_under_2px": (stitch_merge["distance"] < 2).mean() * 100,
                "barcode_rate": cells_with_barcode / len(stitch_final) * 100
                if stitch_final is not None and cells_with_barcode
                else None,
                "single_gene_rate": cells_single_gene / len(stitch_final) * 100
                if stitch_final is not None and cells_single_gene
                else None,
                "phenotype_retention": len(stitch_final) / phenotype_cells * 100
                if phenotype_cells and stitch_final is not None
                else None,
                "sbs_retention": len(stitch_final) / sbs_cells * 100
                if sbs_cells and stitch_final is not None
                else None,
            }
        )

    return pd.DataFrame(results)


def load_distance_data():
    """Load raw distance data for distribution plots."""
    distances = {}

    fast_path = FAST_OUTPUT / "parquets" / f"P-{PLATE}_W-{WELL}__fast_merge.parquet"
    if fast_path.exists():
        fast_df = pd.read_parquet(fast_path, columns=["distance"])
        distances["fast"] = fast_df["distance"].values

    stitch_path = (
        STITCH_OUTPUT / "parquets" / f"P-{PLATE}_W-{WELL}__merged_cells.parquet"
    )
    if stitch_path.exists():
        stitch_df = pd.read_parquet(stitch_path, columns=["distance"])
        distances["stitch"] = stitch_df["distance"].values

    return distances


def plot_runtime_comparison(timing_df, output_path):
    """Create horizontal stacked bar chart comparing runtime by approach."""
    setup_plot_style()

    fig, ax = plt.subplots(figsize=(12, 4))

    approaches = ["fast", "stitch"]
    y_positions = [1, 0]  # Fast on top, Stitch below
    bar_height = 0.6

    for approach, y_pos in zip(approaches, y_positions):
        approach_data = timing_df[timing_df["approach"] == approach].copy()
        if approach_data.empty:
            continue

        # Sort by execution order (alphabetical approximates this)
        approach_data = approach_data.sort_values("step")

        left = 0
        color = COLORS.get(approach, "#7f7f7f")
        total_time = 0

        # Only label these major stitch steps
        label_steps = {"stitch_merge", "stitch_phenotype", "stitch_sbs"}

        for _, row in approach_data.iterrows():
            runtime_min = row["runtime_seconds"] / 60 if row["runtime_seconds"] else 0
            total_time += runtime_min

            # Draw bar segment
            ax.barh(
                y_pos,
                runtime_min,
                left=left,
                height=bar_height,
                color=color,
                alpha=0.8,
                edgecolor="white",
                linewidth=1,
            )

            # Only label major stitch steps
            if row["step"] in label_steps and runtime_min > 30:
                # Format step name for display (replace underscore with space)
                step_label = row["step"].replace("_", " ")
                ax.text(
                    left + runtime_min / 2,
                    y_pos,
                    f"{step_label}\n({runtime_min:.0f}m)",
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white",
                    fontweight="bold",
                )

            left += runtime_min

        # Total time label at end
        ax.text(
            left + 5,
            y_pos,
            f"Total: {total_time:.0f} min",
            va="center",
            fontsize=10,
            fontweight="bold",
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(["Fast", "Stitch"], fontsize=11, fontweight="bold")
    ax.set_xlabel("Runtime (minutes)")
    ax.set_title(f"Merge Runtime Comparison - P-{PLATE}_W-{WELL}")
    ax.set_ylim(-0.5, 2)
    ax.set_xlim(0, ax.get_xlim()[1] * 1.15)

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_distance_distribution(distances, output_path):
    """Create violin plot comparing distance distributions."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE["double"])

    # Violin plot
    ax = axes[0]
    data = []
    labels = []
    colors = []
    for approach, dist in distances.items():
        # Sample for visualization (full data too large)
        sample = np.random.choice(dist, min(50000, len(dist)), replace=False)
        data.append(sample)
        labels.append(f"{approach.capitalize()}\n(n={len(dist):,})")
        colors.append(COLORS.get(approach, "#7f7f7f"))

    parts = ax.violinplot(data, showmeans=False, showmedians=False, showextrema=False)
    # Color the violin bodies
    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(colors[i])
        pc.set_alpha(0.7)

    ax.set_xticks([1, 2])
    ax.set_xticklabels(labels)
    ax.set_ylabel("Distance (pixels)")
    ax.set_title("Match Distance Distribution")

    # Histogram comparison
    ax = axes[1]
    bins = np.linspace(0, 12, 50)
    for approach, dist in distances.items():
        ax.hist(
            dist,
            bins=bins,
            alpha=0.6,
            label=f"{approach.capitalize()}",
            density=True,
            color=COLORS.get(approach, "#7f7f7f"),
            edgecolor="white",
            linewidth=0.5,
        )

    ax.set_xlabel("Distance (pixels)")
    ax.set_ylabel("Density")
    ax.set_title("Distance Histogram")
    ax.legend()

    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_cell_funnel(quality_df, output_path):
    """Create cell funnel plot showing retention through pipeline stages."""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left plot: Absolute counts through pipeline
    ax = axes[0]
    stages = [
        "Phenotype\nInput",
        "SBS\nInput",
        "Matched",
        "Final",
        "With\nBarcode",
        "Single\nGene",
    ]
    x = np.arange(len(stages))
    width = 0.35

    for i, approach in enumerate(quality_df["approach"].tolist()):
        row = quality_df[quality_df["approach"] == approach].iloc[0]
        values = [
            row["phenotype_cells"],
            row["sbs_cells"],
            row["matched_cells"],
            row["final_cells"],
            row["cells_with_barcode"],
            row["cells_single_gene"],
        ]
        color = COLORS.get(approach, "#7f7f7f")
        bars = ax.bar(
            x + i * width - width / 2,
            values,
            width,
            label=approach.capitalize(),
            color=color,
            edgecolor="white",
            linewidth=0.5,
        )
        # Add value labels
        for bar, val in zip(bars, values):
            if val:
                ax.annotate(
                    f"{val / 1e6:.2f}M",
                    xy=(bar.get_x() + bar.get_width() / 2, val),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                    fontweight="bold",
                )

    ax.set_ylabel("Cell Count")
    ax.set_title("Cell Counts Through Pipeline")
    ax.set_xticks(x)
    ax.set_xticklabels(stages, fontsize=9)
    ax.legend(loc="upper right")

    # Right plot: Retention rates all relative to phenotype input
    ax = axes[1]
    stages = ["Matched", "Final", "With Barcode", "Single Gene"]
    x = np.arange(len(stages))

    for i, approach in enumerate(quality_df["approach"].tolist()):
        row = quality_df[quality_df["approach"] == approach].iloc[0]
        ph_cells = row["phenotype_cells"]
        values = [
            row["matched_cells"] / ph_cells * 100,
            row["final_cells"] / ph_cells * 100,
            row["cells_with_barcode"] / ph_cells * 100,
            row["cells_single_gene"] / ph_cells * 100,
        ]
        color = COLORS.get(approach, "#7f7f7f")
        bars = ax.bar(
            x + i * width - width / 2,
            values,
            width,
            label=approach.capitalize(),
            color=color,
            edgecolor="white",
            linewidth=0.5,
        )
        # Add value labels
        for bar, val in zip(bars, values):
            ax.annotate(
                f"{val:.0f}%",
                xy=(bar.get_x() + bar.get_width() / 2, val),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
            )

    ax.set_ylabel("% of Phenotype Input")
    ax.set_title("Retention Through Pipeline")
    ax.set_xticks(x)
    ax.set_xticklabels(stages, fontsize=9)
    ax.legend(loc="upper right")
    ax.set_ylim(0, 130)

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_tile_match_visualization(output_path):
    """Create visualization comparing merged vs unmerged cells for fast and stitch approaches."""
    setup_plot_style()

    # Load tile metadata for positions
    meta = pd.read_parquet(
        BRIEFLOW_OUTPUT
        / "preprocess"
        / "metadata"
        / "phenotype"
        / f"P-{PLATE}_W-{WELL}__combined_metadata.parquet"
    )
    tile_pos = meta[
        ["tile", "x_pos", "y_pos", "pixel_size_x", "pixel_size_y"]
    ].drop_duplicates()

    # Load phenotype cells
    ph = pd.read_parquet(
        BRIEFLOW_OUTPUT
        / "phenotype"
        / "parquets"
        / f"P-{PLATE}_W-{WELL}__phenotype_info.parquet"
    )

    # Merge to get tile positions
    ph = ph.merge(tile_pos, on="tile")

    # Compute global coordinates (microns)
    ph["global_j"] = ph["x_pos"] + ph["j"] * ph["pixel_size_x"]
    ph["global_i"] = ph["y_pos"] + ph["i"] * ph["pixel_size_y"]

    # Load FINAL merge results (after deduplication)
    # tile = phenotype tile, site = SBS tile, cell_0 = phenotype cell, cell_1 = SBS cell
    fast = pd.read_parquet(
        FAST_OUTPUT / "parquets" / f"P-{PLATE}_W-{WELL}__merge_final.parquet",
        columns=["tile", "site", "cell_0", "cell_1"],
    )

    # Try to load stitch final results
    stitch_path = (
        STITCH_OUTPUT / "parquets" / f"P-{PLATE}_W-{WELL}__merge_final.parquet"
    )
    has_stitch = stitch_path.exists()
    if has_stitch:
        stitch = pd.read_parquet(
            stitch_path, columns=["tile", "site", "cell_0", "cell_1"]
        )

    # Mark merged cells - need to match on tile + cell
    fast_matched = set(zip(fast["tile"], fast["cell_0"]))
    ph["fast_merged"] = list(zip(ph["tile"], ph["cell"]))
    ph["fast_merged"] = ph["fast_merged"].apply(lambda x: x in fast_matched)

    if has_stitch:
        stitch_matched = set(zip(stitch["tile"], stitch["cell_0"]))
        ph["stitch_merged"] = list(zip(ph["tile"], ph["cell"]))
        ph["stitch_merged"] = ph["stitch_merged"].apply(lambda x: x in stitch_matched)

    # Select center region for visualization (sample 10% of cells for performance)
    j_center = ph["global_j"].median()
    i_center = ph["global_i"].median()
    j_range = (ph["global_j"].max() - ph["global_j"].min()) * 0.15
    i_range = (ph["global_i"].max() - ph["global_i"].min()) * 0.15

    region = ph[
        (ph["global_j"] > j_center - j_range)
        & (ph["global_j"] < j_center + j_range)
        & (ph["global_i"] > i_center - i_range)
        & (ph["global_i"] < i_center + i_range)
    ].copy()

    # Get axis limits from cell data
    x_min, x_max = region["global_j"].min(), region["global_j"].max()
    y_min, y_max = region["global_i"].min(), region["global_i"].max()
    x_pad = (x_max - x_min) * 0.02
    y_pad = (y_max - y_min) * 0.02

    # Compute phenotype tile outlines
    ph_tile_size = 2960 * tile_pos["pixel_size_x"].iloc[0]  # ~898 µm
    ph_tiles = []
    for _, row in tile_pos.iterrows():
        t_x, t_y = row["x_pos"], row["y_pos"]
        # Check if tile overlaps with region
        if (
            t_x + ph_tile_size > x_min
            and t_x < x_max
            and t_y + ph_tile_size > y_min
            and t_y < y_max
        ):
            ph_tiles.append((t_x, t_y, ph_tile_size))

    # Load SBS cells and metadata
    sbs_meta = pd.read_parquet(
        BRIEFLOW_OUTPUT
        / "preprocess"
        / "metadata"
        / "sbs"
        / f"P-{PLATE}_W-{WELL}__combined_metadata.parquet"
    )
    sbs_tile_pos = sbs_meta[
        ["tile", "x_pos", "y_pos", "pixel_size_x", "pixel_size_y"]
    ].drop_duplicates()

    sbs_cells = pd.read_parquet(
        BRIEFLOW_OUTPUT / "sbs" / "parquets" / f"P-{PLATE}_W-{WELL}__sbs_info.parquet"
    )
    sbs_cells = sbs_cells.merge(sbs_tile_pos, on="tile")
    sbs_cells["global_j"] = (
        sbs_cells["x_pos"] + sbs_cells["j"] * sbs_cells["pixel_size_x"]
    )
    sbs_cells["global_i"] = (
        sbs_cells["y_pos"] + sbs_cells["i"] * sbs_cells["pixel_size_y"]
    )

    # Mark SBS cells as matched - use 'site' (SBS tile) not 'tile' (phenotype tile)
    fast_sbs_matched = set(zip(fast["site"], fast["cell_1"]))
    sbs_cells["fast_merged"] = list(zip(sbs_cells["tile"], sbs_cells["cell"]))
    sbs_cells["fast_merged"] = sbs_cells["fast_merged"].apply(
        lambda x: x in fast_sbs_matched
    )

    if has_stitch:
        stitch_sbs_matched = set(zip(stitch["site"], stitch["cell_1"]))
        sbs_cells["stitch_merged"] = list(zip(sbs_cells["tile"], sbs_cells["cell"]))
        sbs_cells["stitch_merged"] = sbs_cells["stitch_merged"].apply(
            lambda x: x in stitch_sbs_matched
        )

    # Select SBS cells in same region
    sbs_region = sbs_cells[
        (sbs_cells["global_j"] > x_min)
        & (sbs_cells["global_j"] < x_max)
        & (sbs_cells["global_i"] > y_min)
        & (sbs_cells["global_i"] < y_max)
    ].copy()

    # Compute SBS tile outlines
    sbs_tile_size = 1476 * sbs_tile_pos["pixel_size_x"].iloc[0]  # ~1787 µm
    sbs_tiles = []
    for _, row in sbs_tile_pos.iterrows():
        t_x, t_y = row["x_pos"], row["y_pos"]
        if (
            t_x + sbs_tile_size > x_min
            and t_x < x_max
            and t_y + sbs_tile_size > y_min
            and t_y < y_max
        ):
            sbs_tiles.append((t_x, t_y, sbs_tile_size))

    import matplotlib.patches as mpatches

    def draw_tiles(ax, tiles, color="#888888"):
        for t_x, t_y, size in tiles:
            rect = mpatches.Rectangle(
                (t_x, t_y),
                size,
                size,
                linewidth=1,
                edgecolor=color,
                facecolor="none",
                alpha=0.5,
                clip_on=True,
            )
            ax.add_patch(rect)

    # ========== PLOT 1: Phenotype cells with phenotype tile outlines ==========
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))

    # Top-left: Fast matched
    ax = axes[0, 0]
    matched = region[region["fast_merged"]]
    ax.scatter(
        matched["global_j"],
        matched["global_i"],
        c="#2ca02c",
        s=2,
        alpha=0.5,
        rasterized=True,
    )
    draw_tiles(ax, ph_tiles, "#1f77b4")
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_title(f"Fast - Matched ({len(matched):,})")
    ax.set_aspect("equal")

    # Top-right: Fast unmatched
    ax = axes[0, 1]
    unmatched = region[~region["fast_merged"]]
    ax.scatter(
        unmatched["global_j"],
        unmatched["global_i"],
        c="#d62728",
        s=2,
        alpha=0.5,
        rasterized=True,
    )
    draw_tiles(ax, ph_tiles, "#1f77b4")
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_title(f"Fast - Unmatched ({len(unmatched):,})")
    ax.set_aspect("equal")

    # Bottom-left: Stitch matched
    ax = axes[1, 0]
    matched = (
        region[region["stitch_merged"]] if has_stitch else region[region["fast_merged"]]
    )
    ax.scatter(
        matched["global_j"],
        matched["global_i"],
        c="#2ca02c",
        s=2,
        alpha=0.5,
        rasterized=True,
    )
    draw_tiles(ax, ph_tiles, "#1f77b4")
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_title(f"Stitch - Matched ({len(matched):,})")
    ax.set_aspect("equal")

    # Bottom-right: Stitch unmatched
    ax = axes[1, 1]
    unmatched = (
        region[~region["stitch_merged"]]
        if has_stitch
        else region[~region["fast_merged"]]
    )
    ax.scatter(
        unmatched["global_j"],
        unmatched["global_i"],
        c="#d62728",
        s=2,
        alpha=0.5,
        rasterized=True,
    )
    draw_tiles(ax, ph_tiles, "#1f77b4")
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_title(f"Stitch - Unmatched ({len(unmatched):,})")
    ax.set_aspect("equal")

    plt.suptitle(
        f"Phenotype Cells with Phenotype Tile Boundaries - P-{PLATE}_W-{WELL}",
        fontsize=14,
        fontweight="bold",
    )
    plt.tight_layout()
    save_figure(fig, output_path.parent / "tiles_phenotype.png")
    plt.close()
    print(f"Saved: {output_path.parent / 'tiles_phenotype.png'}")

    # ========== PLOT 2: SBS cells with SBS tile outlines ==========
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))

    # Top-left: Fast matched
    ax = axes[0, 0]
    matched = sbs_region[sbs_region["fast_merged"]]
    ax.scatter(
        matched["global_j"],
        matched["global_i"],
        c="#2ca02c",
        s=2,
        alpha=0.5,
        rasterized=True,
    )
    draw_tiles(ax, sbs_tiles, "#ff7f0e")
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_title(f"Fast - Matched ({len(matched):,})")
    ax.set_aspect("equal")

    # Top-right: Fast unmatched
    ax = axes[0, 1]
    unmatched = sbs_region[~sbs_region["fast_merged"]]
    ax.scatter(
        unmatched["global_j"],
        unmatched["global_i"],
        c="#d62728",
        s=2,
        alpha=0.5,
        rasterized=True,
    )
    draw_tiles(ax, sbs_tiles, "#ff7f0e")
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_title(f"Fast - Unmatched ({len(unmatched):,})")
    ax.set_aspect("equal")

    # Bottom-left: Stitch matched
    ax = axes[1, 0]
    matched = (
        sbs_region[sbs_region["stitch_merged"]]
        if has_stitch
        else sbs_region[sbs_region["fast_merged"]]
    )
    ax.scatter(
        matched["global_j"],
        matched["global_i"],
        c="#2ca02c",
        s=2,
        alpha=0.5,
        rasterized=True,
    )
    draw_tiles(ax, sbs_tiles, "#ff7f0e")
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_title(f"Stitch - Matched ({len(matched):,})")
    ax.set_aspect("equal")

    # Bottom-right: Stitch unmatched
    ax = axes[1, 1]
    unmatched = (
        sbs_region[~sbs_region["stitch_merged"]]
        if has_stitch
        else sbs_region[~sbs_region["fast_merged"]]
    )
    ax.scatter(
        unmatched["global_j"],
        unmatched["global_i"],
        c="#d62728",
        s=2,
        alpha=0.5,
        rasterized=True,
    )
    draw_tiles(ax, sbs_tiles, "#ff7f0e")
    ax.set_xlim(x_min - x_pad, x_max + x_pad)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)
    ax.set_title(f"Stitch - Unmatched ({len(unmatched):,})")
    ax.set_aspect("equal")

    plt.suptitle(
        f"SBS Cells with SBS Tile Boundaries - P-{PLATE}_W-{WELL}",
        fontsize=14,
        fontweight="bold",
    )
    plt.tight_layout()
    save_figure(fig, output_path.parent / "tiles_sbs.png")
    plt.close()
    print(f"Saved: {output_path.parent / 'tiles_sbs.png'}")

    # ========== OVERLAP ANALYSIS ==========
    # For each cell in the plotted region, determine if it's in a tile overlap
    print("\n=== Overlap Analysis (plotted region) ===")

    def count_covering_tiles(global_j, global_i, tiles_df, tile_size):
        """Count how many tiles cover a given position."""
        count = 0
        for _, t in tiles_df.iterrows():
            if (
                t["x_pos"] <= global_j <= t["x_pos"] + tile_size
                and t["y_pos"] <= global_i <= t["y_pos"] + tile_size
            ):
                count += 1
        return count

    # Sample cells for efficiency (use all if < 5000)
    if len(region) > 5000:
        analysis_sample = region.sample(n=5000, random_state=42).copy()
    else:
        analysis_sample = region.copy()

    # Compute overlap status for phenotype tiles
    ph_tile_size_val = 2960 * tile_pos["pixel_size_x"].iloc[0]
    analysis_sample["n_tiles"] = analysis_sample.apply(
        lambda row: count_covering_tiles(
            row["global_j"], row["global_i"], tile_pos, ph_tile_size_val
        ),
        axis=1,
    )

    no_overlap = analysis_sample[analysis_sample["n_tiles"] == 1]
    in_overlap = analysis_sample[analysis_sample["n_tiles"] >= 2]

    print(f"\nPhenotype cells in plotted region: {len(analysis_sample):,}")
    print(
        f"  Outside overlap (1 tile):  {len(no_overlap):,} ({len(no_overlap) / len(analysis_sample) * 100:.1f}%)"
    )
    print(
        f"  Inside overlap (2+ tiles): {len(in_overlap):,} ({len(in_overlap) / len(analysis_sample) * 100:.1f}%)"
    )

    # Fast approach breakdown
    print("\n--- Fast Approach ---")
    fast_no_overlap_matched = no_overlap["fast_merged"].sum()
    fast_no_overlap_unmatched = len(no_overlap) - fast_no_overlap_matched
    fast_in_overlap_matched = in_overlap["fast_merged"].sum()
    fast_in_overlap_unmatched = len(in_overlap) - fast_in_overlap_matched

    print(
        f"  Outside overlap: {fast_no_overlap_matched:,} matched ({fast_no_overlap_matched / len(no_overlap) * 100:.1f}%), "
        f"{fast_no_overlap_unmatched:,} unmatched ({fast_no_overlap_unmatched / len(no_overlap) * 100:.1f}%)"
    )
    print(
        f"  Inside overlap:  {fast_in_overlap_matched:,} matched ({fast_in_overlap_matched / len(in_overlap) * 100:.1f}%), "
        f"{fast_in_overlap_unmatched:,} unmatched ({fast_in_overlap_unmatched / len(in_overlap) * 100:.1f}%)"
    )

    # Stitch approach breakdown
    if has_stitch:
        print("\n--- Stitch Approach ---")
        stitch_no_overlap_matched = no_overlap["stitch_merged"].sum()
        stitch_no_overlap_unmatched = len(no_overlap) - stitch_no_overlap_matched
        stitch_in_overlap_matched = in_overlap["stitch_merged"].sum()
        stitch_in_overlap_unmatched = len(in_overlap) - stitch_in_overlap_matched

        print(
            f"  Outside overlap: {stitch_no_overlap_matched:,} matched ({stitch_no_overlap_matched / len(no_overlap) * 100:.1f}%), "
            f"{stitch_no_overlap_unmatched:,} unmatched ({stitch_no_overlap_unmatched / len(no_overlap) * 100:.1f}%)"
        )
        print(
            f"  Inside overlap:  {stitch_in_overlap_matched:,} matched ({stitch_in_overlap_matched / len(in_overlap) * 100:.1f}%), "
            f"{stitch_in_overlap_unmatched:,} unmatched ({stitch_in_overlap_unmatched / len(in_overlap) * 100:.1f}%)"
        )

    # Save overlap analysis to CSV
    overlap_data = {
        "region": ["outside_overlap", "inside_overlap"],
        "n_cells": [len(no_overlap), len(in_overlap)],
        "fast_matched": [fast_no_overlap_matched, fast_in_overlap_matched],
        "fast_unmatched": [fast_no_overlap_unmatched, fast_in_overlap_unmatched],
        "fast_match_rate": [
            fast_no_overlap_matched / len(no_overlap) * 100,
            fast_in_overlap_matched / len(in_overlap) * 100,
        ],
    }
    if has_stitch:
        overlap_data["stitch_matched"] = [
            stitch_no_overlap_matched,
            stitch_in_overlap_matched,
        ]
        overlap_data["stitch_unmatched"] = [
            stitch_no_overlap_unmatched,
            stitch_in_overlap_unmatched,
        ]
        overlap_data["stitch_match_rate"] = [
            stitch_no_overlap_matched / len(no_overlap) * 100,
            stitch_in_overlap_matched / len(in_overlap) * 100,
        ]

    overlap_df = pd.DataFrame(overlap_data)
    overlap_df.to_csv(output_path.parent / "overlap_analysis.csv", index=False)
    print(f"\nSaved: {output_path.parent / 'overlap_analysis.csv'}")


def plot_memory_requirements(output_path):
    """Plot memory requirements from slurm config for fast vs stitch approaches."""
    setup_plot_style()

    # Load slurm config
    with open(SLURM_CONFIG, "r") as f:
        config = yaml.safe_load(f)

    resources = config.get("set-resources", {})

    # Define steps for each approach (matching FAST_JOBS and STITCH_JOBS)
    fast_steps = list(FAST_JOBS.keys())
    stitch_steps = list(STITCH_JOBS.keys())

    # Extract memory for each step
    data = []
    for step in fast_steps:
        mem = resources.get(step, {}).get(
            "mem_mb", config.get("default-resources", {}).get("mem_mb", 3000)
        )
        data.append({"approach": "Fast", "step": step, "memory_gb": mem / 1000})

    for step in stitch_steps:
        mem = resources.get(step, {}).get(
            "mem_mb", config.get("default-resources", {}).get("mem_mb", 3000)
        )
        data.append({"approach": "Stitch", "step": step, "memory_gb": mem / 1000})

    df = pd.DataFrame(data)

    # Calculate totals
    fast_total = df[df["approach"] == "Fast"]["memory_gb"].sum()
    stitch_total = df[df["approach"] == "Stitch"]["memory_gb"].sum()

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=FIGSIZE["double"])

    # Use consistent colors from plot_style
    fast_color = COLORS["fast"]
    stitch_color = COLORS["stitch"]
    fast_df = df[df["approach"] == "Fast"]
    stitch_df = df[df["approach"] == "Stitch"]

    width = 0.7

    # Plot fast steps
    ax1.barh(
        np.arange(len(fast_df)),
        fast_df["memory_gb"].values,
        width,
        label="Fast",
        color=fast_color,
        alpha=0.8,
    )

    # Plot stitch steps (offset)
    ax1.barh(
        np.arange(len(stitch_df)) + len(fast_df) + 1,
        stitch_df["memory_gb"].values,
        width,
        label="Stitch",
        color=stitch_color,
        alpha=0.8,
    )

    # Labels
    all_steps = list(fast_df["step"].values) + [""] + list(stitch_df["step"].values)
    ax1.set_yticks(range(len(all_steps)))
    ax1.set_yticklabels(all_steps)
    ax1.set_xlabel("Memory Requested (GB)")
    ax1.set_title("Memory Requirements by Step")
    ax1.axvline(x=1000, color="red", linestyle="--", alpha=0.5, label="1 TB limit")
    ax1.legend()

    # Right: Total comparison
    totals = [fast_total, stitch_total]
    approaches = ["Fast", "Stitch"]
    bars = ax2.bar(approaches, totals, color=[fast_color, stitch_color], alpha=0.8)

    ax2.axhline(
        y=1000, color="red", linestyle="--", alpha=0.7, label="1 TB cluster limit"
    )
    ax2.set_ylabel("Total Memory Requested (GB)")
    ax2.set_title("Total Memory Requirements")

    # Add value labels using helper
    add_value_labels(ax2, bars, fmt="{:.0f} GB")

    ax2.legend()

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def main():
    args = parse_args()

    # Ensure output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not args.plots_only:
        print("=" * 60)
        print("MERGE SNAKEMAKE BENCHMARK")
        print(f"Well: P-{PLATE}_W-{WELL}")
        print("=" * 60)

        # Extract timing data
        print("\nExtracting timing from slurm logs...")
        timing_df = extract_timing_from_logs()
        print(f"  Found {len(timing_df)} timing records")

        # Extract quality metrics
        print("\nExtracting quality metrics from parquets...")
        quality_df = extract_quality_metrics()
        print(f"  Found {len(quality_df)} quality records")

        # Save CSV
        results_df = quality_df.merge(
            timing_df.groupby("approach")
            .agg(
                {
                    "runtime_seconds": "sum",
                    "memory_requested_mb": "max",
                    "status": lambda x: "success"
                    if all(s == "success" for s in x)
                    else "degraded",
                }
            )
            .reset_index()
            .rename(columns={"runtime_seconds": "total_runtime_seconds"}),
            on="approach",
        )
        results_df["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        csv_path = OUTPUT_DIR / "snakemake_benchmark_results.csv"
        results_df.to_csv(csv_path, index=False)
        print(f"\nSaved: {csv_path}")

        # Save timing detail
        timing_csv = OUTPUT_DIR / "snakemake_timing_detail.csv"
        timing_df.to_csv(timing_csv, index=False)
        print(f"Saved: {timing_csv}")
    else:
        # Load existing data
        csv_path = OUTPUT_DIR / "snakemake_benchmark_results.csv"
        timing_csv = OUTPUT_DIR / "snakemake_timing_detail.csv"
        quality_df = pd.read_csv(csv_path)
        timing_df = pd.read_csv(timing_csv)

    # Generate plots
    print("\nGenerating plots...")

    plot_runtime_comparison(timing_df, OUTPUT_DIR / "runtime.png")

    distances = load_distance_data()
    if distances:
        plot_distance_distribution(distances, OUTPUT_DIR / "alignment.png")

    plot_cell_funnel(quality_df, OUTPUT_DIR / "retention.png")
    plot_tile_match_visualization(OUTPUT_DIR / "tiles.png")
    plot_memory_requirements(OUTPUT_DIR / "memory.png")

    print("\n" + "=" * 60)
    print("BENCHMARK COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
