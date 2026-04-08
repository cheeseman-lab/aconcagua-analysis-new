"""
Clustering Benchmark - Brieflow vs Funk et al. 2022

Rigorous head-to-head comparison addressing reviewer feedback:
- Quantitative benchmarking with STRING, CORUM, KEGG enrichment
- Precision-Recall curves comparing clustering methods
- Resolution sweep analysis (k=5-15) for both methods

Metrics compared:
- STRING F1 score (protein-protein interaction enrichment)
- CORUM enrichment (% clusters with significant complex enrichment)
- KEGG enrichment (% clusters with significant pathway enrichment)
- Precision-Recall curves across clustering resolutions

Usage:
    python cluster_enrichment.py              # Generate all outputs
    python cluster_enrichment.py --plots-only # Regenerate plots only
"""

import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Import shared plot styling
from plot_style import setup_plot_style, save_figure, COLORS

# Configure paths
BENCHMARKS_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BENCHMARKS_DIR.parent
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
BRIEFLOW_OUTPUT = ANALYSIS_DIR / "brieflow_output" / "cluster"
FUNK_RESULTS = BENCHMARKS_DIR / "external" / "results" / "cluster"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "cluster" / "enrichment"

# Clustering configurations to analyze (ideal resolutions)
IDEAL_CONFIGS = [
    {
        "name": "Interphase",
        "cell_class": "Interphase",
        "brieflow_k": 12,
        "funk_k": 10,  # Funk's published ideal resolution
        "display_name": "Interphase",
    },
    {
        "name": "Mitotic",
        "cell_class": "Mitotic",
        "brieflow_k": 5,
        "funk_k": 9,  # Funk's published ideal resolution
        "display_name": "Mitotic",
    },
]

# Resolution sweep range
K_RANGE = range(5, 16)  # k=5 to k=15


def parse_args():
    parser = argparse.ArgumentParser(description="Clustering benchmark")
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Only regenerate plots from existing results",
    )
    return parser.parse_args()


def load_funk_enrichment(cell_class, k):
    """Load Funk clustering enrichment results for a specific k value."""
    results_dir = FUNK_RESULTS / f"Funk_{cell_class}_k{k}"
    metrics_file = results_dir / "CB-Real__global_metrics.json"

    if not metrics_file.exists():
        return None

    with open(metrics_file, "r") as f:
        metrics = json.load(f)

    # Load clustering file to get cluster count
    clustering_file = results_dir / "phate_leiden_clustering.tsv"
    if clustering_file.exists():
        clustering_df = pd.read_table(clustering_file)
        num_clusters = clustering_df["cluster"].nunique()
    else:
        num_clusters = None

    return {
        "metrics": metrics,
        "num_clusters": num_clusters,
        "k": k,
        "cell_class": cell_class,
    }


def load_brieflow_enrichment(cell_class, k):
    """Load Brieflow clustering enrichment results for a specific k value."""
    # Load from resolution sweep data
    results_dir = (
        BRIEFLOW_OUTPUT / "DAPI_TUBULIN_GH2AX_PHALLOIDIN" / cell_class / str(k)
    )
    metrics_file = results_dir / "CB-Real__global_metrics.json"

    if not metrics_file.exists():
        return None

    with open(metrics_file, "r") as f:
        metrics = json.load(f)

    clustering_file = results_dir / "phate_leiden_clustering.tsv"
    if clustering_file.exists():
        clustering_df = pd.read_table(clustering_file)
        num_clusters = clustering_df["cluster"].nunique()
    else:
        num_clusters = None

    return {
        "metrics": metrics,
        "num_clusters": num_clusters,
        "k": k,
        "cell_class": cell_class,
    }


def build_resolution_sweep_data():
    """Build data for resolution sweep comparison."""
    results = []

    for cell_class in ["Interphase", "Mitotic"]:
        for k in K_RANGE:
            # Load Funk results
            funk = load_funk_enrichment(cell_class, k)
            if funk:
                row = {
                    "method": "Funk et al. 2022",
                    "cell_class": cell_class,
                    "k": k,
                    "num_clusters": funk["num_clusters"],
                    "string_f1": funk["metrics"]
                    .get("STRING", {})
                    .get("f1_score", np.nan),
                    "string_precision": funk["metrics"]
                    .get("STRING", {})
                    .get("precision", np.nan),
                    "string_recall": funk["metrics"]
                    .get("STRING", {})
                    .get("recall", np.nan),
                    "corum_enrich": funk["metrics"]
                    .get("CORUM", {})
                    .get("proportion_enriched", np.nan),
                    "kegg_enrich": funk["metrics"]
                    .get("KEGG", {})
                    .get("proportion_enriched", np.nan),
                }
                results.append(row)

            # Load Brieflow results (only available at ideal k for now)
            brieflow = load_brieflow_enrichment(cell_class, k)
            if brieflow:
                row = {
                    "method": "Brieflow",
                    "cell_class": cell_class,
                    "k": k,
                    "num_clusters": brieflow["num_clusters"],
                    "string_f1": brieflow["metrics"]
                    .get("STRING", {})
                    .get("f1_score", np.nan),
                    "string_precision": brieflow["metrics"]
                    .get("STRING", {})
                    .get("precision", np.nan),
                    "string_recall": brieflow["metrics"]
                    .get("STRING", {})
                    .get("recall", np.nan),
                    "corum_enrich": brieflow["metrics"]
                    .get("CORUM", {})
                    .get("proportion_enriched", np.nan),
                    "kegg_enrich": brieflow["metrics"]
                    .get("KEGG", {})
                    .get("proportion_enriched", np.nan),
                }
                results.append(row)

    return pd.DataFrame(results)


def _extract_metrics(data):
    """Extract standard metrics dict from loaded enrichment data."""
    m = data["metrics"]
    return {
        "string_f1": m.get("STRING", {}).get("f1_score", np.nan),
        "string_precision": m.get("STRING", {}).get("precision", np.nan),
        "string_recall": m.get("STRING", {}).get("recall", np.nan),
        "corum_enrich": m.get("CORUM", {}).get("proportion_enriched", np.nan),
        "kegg_enrich": m.get("KEGG", {}).get("proportion_enriched", np.nan),
        "num_clusters": data["num_clusters"],
    }


def load_shuffled_enrichment(cell_class, k):
    """Load shuffled control enrichment results."""
    results_dir = (
        BRIEFLOW_OUTPUT / "DAPI_TUBULIN_GH2AX_PHALLOIDIN" / cell_class / str(k)
    )
    metrics_file = results_dir / "CB-Shuffled__global_metrics.json"

    if not metrics_file.exists():
        return None

    with open(metrics_file, "r") as f:
        metrics = json.load(f)

    # Shuffled uses same clustering file for cluster count
    clustering_file = results_dir / "phate_leiden_clustering_shuffled.tsv"
    if clustering_file.exists():
        clustering_df = pd.read_table(clustering_file)
        num_clusters = clustering_df["cluster"].nunique()
    else:
        # Fall back to the real clustering count
        clustering_file = results_dir / "phate_leiden_clustering.tsv"
        if clustering_file.exists():
            clustering_df = pd.read_table(clustering_file)
            num_clusters = clustering_df["cluster"].nunique()
        else:
            num_clusters = None

    return {
        "metrics": metrics,
        "num_clusters": num_clusters,
        "k": k,
        "cell_class": cell_class,
    }


def build_ideal_comparison_table():
    """Build comparison table at ideal resolutions (Brieflow, Funk, Shuffled)."""
    results = []

    for config in IDEAL_CONFIGS:
        cell_class = config["cell_class"]

        funk = load_funk_enrichment(cell_class, config["funk_k"])
        brieflow = load_brieflow_enrichment(cell_class, config["brieflow_k"])
        shuffled = load_shuffled_enrichment(cell_class, config["brieflow_k"])

        if funk and brieflow:
            row = {
                "cell_class": cell_class,
                "display_name": config["display_name"],
            }
            # Add prefixed metrics for each pipeline
            for prefix, data in [
                ("funk", funk),
                ("brieflow", brieflow),
            ]:
                m = _extract_metrics(data)
                row[f"{prefix}_k"] = data["k"]
                for key, val in m.items():
                    row[f"{prefix}_{key}"] = val

            # Shuffled (may not exist)
            if shuffled:
                m = _extract_metrics(shuffled)
                row["shuffled_k"] = shuffled["k"]
                for key, val in m.items():
                    row[f"shuffled_{key}"] = val
            else:
                row["shuffled_k"] = np.nan
                for key in ["string_f1", "corum_enrich", "kegg_enrich", "num_clusters"]:
                    row[f"shuffled_{key}"] = np.nan

            results.append(row)

    df = pd.DataFrame(results)

    # Calculate fold changes
    df["f1_improvement"] = (df["brieflow_string_f1"] - df["funk_string_f1"]) / df[
        "funk_string_f1"
    ]
    df["corum_improvement"] = (
        df["brieflow_corum_enrich"] - df["funk_corum_enrich"]
    ) / df["funk_corum_enrich"]
    df["kegg_improvement"] = (df["brieflow_kegg_enrich"] - df["funk_kegg_enrich"]) / df[
        "funk_kegg_enrich"
    ]

    return df


def plot_ideal_comparison(comparison_df, output_path):
    """Plot head-to-head comparison at ideal resolutions (3-panel: STRING, CORUM, KEGG).

    Each panel shows grouped bars: Brieflow, Funk, Shuffled per cell class.
    """
    setup_plot_style()

    has_shuffled = "shuffled_string_f1" in comparison_df.columns and comparison_df[
        "shuffled_string_f1"
    ].notna().any()

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    metrics = [
        ("string_f1", "STRING F1 Score", axes[0]),
        ("corum_enrich", "CORUM Enrichment (%)", axes[1]),
        ("kegg_enrich", "KEGG Enrichment (%)", axes[2]),
    ]

    n_groups = len(comparison_df)
    n_bars = 3 if has_shuffled else 2
    width = 0.25
    offsets = np.linspace(-(n_bars - 1) * width / 2, (n_bars - 1) * width / 2, n_bars)

    bar_configs = [
        ("brieflow", "Brieflow", COLORS["brieflow"]),
        ("funk", "Funk", COLORS["funk"]),
    ]
    if has_shuffled:
        bar_configs.append(("shuffled", "Shuffled", "#999999"))

    for metric, title, ax in metrics:
        x = np.arange(n_groups)

        all_values = []
        for i, (prefix, label, color) in enumerate(bar_configs):
            values = comparison_df[f"{prefix}_{metric}"].values.copy()
            if "enrich" in metric:
                values = values * 100
            all_values.append(values)

            ax.bar(
                x + offsets[i],
                values,
                width,
                label=label,
                color=color,
                alpha=0.85,
            )

            # Value labels on bars
            for j, val in enumerate(values):
                if np.isnan(val):
                    continue
                if "f1" in metric:
                    txt = f"{val:.3f}"
                else:
                    txt = f"{val:.1f}%"
                ax.annotate(
                    txt,
                    xy=(x[j] + offsets[i], max(val, 0)),
                    xytext=(0, 4),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    fontweight="bold",
                )

        ax.set_ylabel(title, fontsize=13, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(comparison_df["display_name"], fontsize=12)
        ax.tick_params(axis="y", labelsize=11)
        ax.grid(True, alpha=0.3, axis="y")

        # Pad y-axis
        ymax = max(np.nanmax(v) for v in all_values)
        ax.set_ylim(0, ymax * 1.3)

    # Single shared legend in the rightmost panel
    axes[2].legend(loc="upper right", fontsize=11, framealpha=0.9)

    fig.suptitle(
        "Clustering Method Comparison at Ideal Resolutions",
        fontsize=16,
        fontweight="bold",
    )
    fig.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_resolution_sweep(sweep_df, output_path):
    """Plot resolution sweep comparison across k=5-15."""
    setup_plot_style()

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    cell_classes = ["Interphase", "Mitotic"]
    colors = {"Funk et al. 2022": COLORS["funk"], "Brieflow": COLORS["brieflow"]}
    markers = {"Funk et al. 2022": "o", "Brieflow": "s"}

    # Plot all subplots but only add labels on the first one
    plot_idx = 0
    for row in range(2):
        for col in range(2):
            ax = axes[row, col]

            # Determine which metric to plot
            if plot_idx == 0:
                metric_key = "string_f1"
                ylabel = "STRING F1 Score"
                title = "STRING Protein Interaction Enrichment"
                scale = 1
            elif plot_idx == 1:
                metric_key = "corum_enrich"
                ylabel = "CORUM Enrichment (%)"
                title = "Protein Complex Enrichment"
                scale = 100
            elif plot_idx == 2:
                metric_key = "kegg_enrich"
                ylabel = "KEGG Enrichment (%)"
                title = "Pathway Enrichment"
                scale = 100
            else:
                metric_key = "num_clusters"
                ylabel = "Number of Clusters"
                title = "Clustering Granularity"
                scale = 1

            for cell_class in cell_classes:
                for method in sweep_df["method"].unique():
                    data = sweep_df[
                        (sweep_df["cell_class"] == cell_class)
                        & (sweep_df["method"] == method)
                    ]
                    if len(data) > 0:
                        linestyle = "-" if cell_class == "Interphase" else "--"
                        # Only add label on first subplot
                        label = f"{method} ({cell_class})" if plot_idx == 0 else None
                        ax.plot(
                            data["k"],
                            data[metric_key] * scale,
                            marker=markers[method],
                            color=colors[method],
                            linestyle=linestyle,
                            linewidth=2,
                            markersize=6,
                            label=label,
                            alpha=0.8,
                        )

            ax.set_xlabel("Clustering Resolution (k)")
            ax.set_ylabel(ylabel)
            ax.set_title(title)
            ax.grid(True, alpha=0.3)

            plot_idx += 1

    # Add legend only to first subplot
    axes[0, 0].legend(loc="best", fontsize=9)

    plt.suptitle(
        "Resolution Sweep: Brieflow vs Funk et al. 2022 (k=5-15)",
        fontsize=14,
        fontweight="bold",
        y=0.995,
    )
    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_pr_curves(sweep_df, output_path):
    """Plot Precision-Recall curves for both methods across resolutions."""
    setup_plot_style()

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    cell_classes = ["Interphase", "Mitotic"]
    colors = {"Funk et al. 2022": COLORS["funk"], "Brieflow": COLORS["brieflow"]}

    # First pass: determine global axis limits
    all_recall = []
    all_precision = []

    for cell_class in cell_classes:
        for method in sweep_df["method"].unique():
            data = sweep_df[
                (sweep_df["cell_class"] == cell_class) & (sweep_df["method"] == method)
            ].sort_values("k")
            if len(data) > 0:
                recall = data["string_recall"].values
                precision = data["string_precision"].values
                valid = ~(np.isnan(recall) | np.isnan(precision))
                all_recall.extend(recall[valid])
                all_precision.extend(precision[valid])

    # Set matched axis limits with some padding
    max_recall = max(all_recall) * 1.1
    max_precision = max(all_precision) * 1.1

    # Second pass: plot data
    for idx, cell_class in enumerate(cell_classes):
        ax = axes[idx]

        for method in sweep_df["method"].unique():
            data = sweep_df[
                (sweep_df["cell_class"] == cell_class) & (sweep_df["method"] == method)
            ].sort_values("k")

            if len(data) > 0:
                # Plot STRING precision-recall curve across resolutions
                recall = data["string_recall"].values
                precision = data["string_precision"].values

                # Remove NaN values
                valid = ~(np.isnan(recall) | np.isnan(precision))
                recall = recall[valid]
                precision = precision[valid]

                if len(recall) > 0:
                    # Only add label for the first plot to avoid duplicate legends
                    label = method if idx == 0 else None
                    ax.plot(
                        recall,
                        precision,
                        "-o",
                        color=colors[method],
                        linewidth=2,
                        markersize=6,
                        label=label,
                        alpha=0.8,
                    )

                    # Annotate k values at key points
                    k_values = data["k"].values[valid]
                    for i in [0, len(k_values) // 2, -1]:  # Start, middle, end
                        if i < len(k_values):
                            ax.annotate(
                                f"k={k_values[i]}",
                                xy=(recall[i], precision[i]),
                                xytext=(5, 5),
                                textcoords="offset points",
                                fontsize=8,
                                alpha=0.7,
                            )

        # Add diagonal reference line (only add to legend on first plot)
        diagonal_label = "Random" if idx == 0 else None
        ax.plot(
            [0, max_recall],
            [0, max_precision],
            "k--",
            alpha=0.3,
            linewidth=1,
            label=diagonal_label,
        )

        ax.set_xlabel("Recall (STRING)", fontsize=11)
        ax.set_ylabel("Precision (STRING)", fontsize=11)
        ax.set_title(f"{cell_class}", fontsize=12, fontweight="bold")

        # Apply matched axis limits
        ax.set_xlim([0, max_recall])
        ax.set_ylim([0, max_precision])

        # Only show legend on first plot
        if idx == 0:
            ax.legend(loc="upper left", fontsize=10)

        ax.grid(True, alpha=0.3)

    plt.suptitle(
        "Precision-Recall Curves: STRING Enrichment Across Resolutions",
        fontsize=14,
        fontweight="bold",
    )
    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def generate_summary_table(comparison_df, output_path):
    """Generate summary table comparing methods at ideal resolutions."""

    display_df = pd.DataFrame(
        {
            "Cell Class": comparison_df["display_name"],
            "Funk k": comparison_df["funk_k"],
            "Funk STRING F1": comparison_df["funk_string_f1"].apply(
                lambda x: f"{x:.4f}"
            ),
            "Funk CORUM %": (comparison_df["funk_corum_enrich"] * 100).apply(
                lambda x: f"{x:.1f}%"
            ),
            "Funk KEGG %": (comparison_df["funk_kegg_enrich"] * 100).apply(
                lambda x: f"{x:.1f}%"
            ),
            "Brieflow k": comparison_df["brieflow_k"],
            "Brieflow STRING F1": comparison_df["brieflow_string_f1"].apply(
                lambda x: f"{x:.4f}"
            ),
            "Brieflow CORUM %": (comparison_df["brieflow_corum_enrich"] * 100).apply(
                lambda x: f"{x:.1f}%"
            ),
            "Brieflow KEGG %": (comparison_df["brieflow_kegg_enrich"] * 100).apply(
                lambda x: f"{x:.1f}%"
            ),
            "STRING F1 Δ": comparison_df["f1_improvement"].apply(
                lambda x: f"{x * 100:+.1f}%"
            ),
            "CORUM Δ": comparison_df["corum_improvement"].apply(
                lambda x: f"{x * 100:+.1f}%"
            ),
            "KEGG Δ": comparison_df["kegg_improvement"].apply(
                lambda x: f"{x * 100:+.1f}%"
            ),
        }
    )

    display_df.to_csv(output_path, index=False)
    print(f"Saved: {output_path}")

    # Print to console
    print("\n" + "=" * 100)
    print("CLUSTERING METHOD COMPARISON AT IDEAL RESOLUTIONS")
    print("=" * 100)
    print(display_df.to_string(index=False))
    print("=" * 100)
    print("\nKey Findings:")
    for _, row in comparison_df.iterrows():
        print(f"\n{row['display_name']}:")
        print(f"  STRING F1: {row['f1_improvement'] * 100:+.1f}% improvement")
        print(f"  CORUM enrichment: {row['corum_improvement'] * 100:+.1f}% improvement")
        print(f"  KEGG enrichment: {row['kegg_improvement'] * 100:+.1f}% improvement")


def main():
    args = parse_args()

    # Ensure output directory exists
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("CLUSTERING BENCHMARK - Brieflow vs Funk et al. 2022")
    print("=" * 80)

    if not args.plots_only:
        # Build resolution sweep data
        print("\nBuilding resolution sweep data...")
        sweep_df = build_resolution_sweep_data()
        sweep_df.to_csv(OUTPUT_DIR / "resolution_sweep.csv", index=False)
        print(f"Saved: {OUTPUT_DIR / 'resolution_sweep.csv'}")

        # Build ideal comparison table
        print("\nBuilding ideal resolution comparison...")
        comparison_df = build_ideal_comparison_table()
        comparison_df.to_csv(OUTPUT_DIR / "ideal_comparison.csv", index=False)
        print(f"Saved: {OUTPUT_DIR / 'ideal_comparison.csv'}")

        # Generate summary table
        generate_summary_table(comparison_df, OUTPUT_DIR / "summary_table.csv")
    else:
        # Load existing data
        sweep_df = pd.read_csv(OUTPUT_DIR / "resolution_sweep.csv")
        comparison_df = pd.read_csv(OUTPUT_DIR / "ideal_comparison.csv")

    # Generate plots
    print("\nGenerating plots...")

    plot_ideal_comparison(comparison_df, OUTPUT_DIR / "ideal_comparison.png")
    plot_resolution_sweep(sweep_df, OUTPUT_DIR / "resolution_sweep.png")
    plot_pr_curves(sweep_df, OUTPUT_DIR / "pr_curves.png")

    print("\n" + "=" * 80)
    print("BENCHMARK COMPLETE")
    print("=" * 80)
    print(f"\nResults saved to: {OUTPUT_DIR}")
    print("\nOutputs:")
    print("- ideal_comparison.png - Head-to-head at ideal k values")
    print("- resolution_sweep.png - Performance across k=5-15")
    print("- pr_curves.png - Precision-recall curves")
    print("- summary_table.csv - Detailed comparison metrics")


if __name__ == "__main__":
    main()
