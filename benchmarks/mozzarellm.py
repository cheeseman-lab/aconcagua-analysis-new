"""
MozzareLLM Manuscript Figures A–C

Generates publication figures from mozzarellm summaries:
  A: High Confidence Cluster Comparison (Brieflow vs Funk, Interphase & Mitotic)
  B: Gene Discovery Potential — stacked bar (Established / Novel / Uncharacterized)
  C: Mean Genes per High Confidence Cluster — grouped bar

Reads summaries from:
  Brieflow:  .../Interphase/12/mozzarellm/claude-opus-4-6_summaries.tsv
             .../Mitotic/5/mozzarellm/claude-opus-4-6_summaries.tsv
  Funk:      benchmarks/external/results/cluster/Funk_Interphase_k10/mozzarellm/claude-opus-4-6_summaries.tsv
             benchmarks/external/results/cluster/Funk_Mitotic_k9/mozzarellm/claude-opus-4-6_summaries.tsv
  Shuffled:  .../Interphase/12/mozzarellm_shuffled/claude-opus-4-6_summaries.tsv (optional)

Usage:
    python mozzarellm.py
    python mozzarellm.py --include-shuffled
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

from plot_style import setup_plot_style, save_figure, COLORS

# Paths
SCRIPT_DIR = Path(__file__).parent
ANALYSIS_DIR = SCRIPT_DIR.parent / "analysis"
CLUSTER_DIR = ANALYSIS_DIR / "brieflow_output" / "cluster" / "DAPI_TUBULIN_GH2AX_PHALLOIDIN"
FUNK_DIR = SCRIPT_DIR / "external" / "results" / "cluster"
OUTPUT_DIR = SCRIPT_DIR / "results" / "mozzarellm"

# Pipeline → cell class → path
PIPELINES = {
    "Brieflow": {
        "Interphase": CLUSTER_DIR / "Interphase" / "12" / "mozzarellm" / "claude-opus-4-6_summaries.tsv",
        "Mitotic": CLUSTER_DIR / "Mitotic" / "5" / "mozzarellm" / "claude-opus-4-6_summaries.tsv",
    },
    "Funk": {
        "Interphase": FUNK_DIR / "Funk_Interphase_k10" / "mozzarellm" / "claude-opus-4-6_summaries.tsv",
        "Mitotic": FUNK_DIR / "Funk_Mitotic_k9" / "mozzarellm" / "claude-opus-4-6_summaries.tsv",
    },
}

SHUFFLED_PATH = CLUSTER_DIR / "Interphase" / "12" / "mozzarellm_shuffled" / "claude-opus-4-6_summaries.tsv"

PIPELINE_COLORS = {
    "Brieflow": COLORS["brieflow"],  # Green
    "Funk": COLORS["funk"],          # Red
    "Shuffled": "#999999",           # Gray
}

CELL_CLASSES = ["Interphase", "Mitotic"]


def load_summaries(include_shuffled=False):
    """Load mozzarellm summaries for all pipelines and cell classes.

    Returns:
        dict of {label: DataFrame} where label is e.g. "Brieflow Interphase"
    """
    data = {}
    for pipeline, classes in PIPELINES.items():
        for cell_class, path in classes.items():
            label = f"{pipeline} {cell_class}"
            if not path.exists():
                print(f"WARNING: {label} summaries not found at {path}")
                continue
            df = pd.read_csv(path, sep="\t")
            data[label] = df
            print(f"Loaded {label}: {len(df)} clusters")

    if include_shuffled:
        if SHUFFLED_PATH.exists():
            df = pd.read_csv(SHUFFLED_PATH, sep="\t")
            data["Shuffled"] = df
            print(f"Loaded Shuffled: {len(df)} clusters")
        else:
            print(f"WARNING: Shuffled summaries not found at {SHUFFLED_PATH}")

    return data


def get_confidence_stats(df):
    """Compute confidence-level statistics for a summaries DataFrame."""
    total = len(df)
    high = (df["pathway_confidence"] == "High").sum()
    medium = (df["pathway_confidence"] == "Medium").sum()
    low = (df["pathway_confidence"] == "Low").sum()

    high_df = df[df["pathway_confidence"] == "High"]
    n_established = high_df["num_established"].sum()
    n_novel = high_df["num_novel"].sum()
    n_uncharacterized = high_df["num_uncharacterized"].sum()
    n_high = len(high_df)

    return {
        "total": total,
        "high": high,
        "medium": medium,
        "low": low,
        "pct_high": 100 * high / total if total else 0,
        "pct_medium": 100 * medium / total if total else 0,
        "pct_low": 100 * low / total if total else 0,
        "genes_established": n_established,
        "genes_novel": n_novel,
        "genes_uncharacterized": n_uncharacterized,
        "genes_total_highconf": n_established + n_novel + n_uncharacterized,
        "mean_established": n_established / n_high if n_high else 0,
        "mean_novel": n_novel / n_high if n_high else 0,
        "mean_uncharacterized": n_uncharacterized / n_high if n_high else 0,
    }


def figure_a(data, output_dir):
    """Figure A: High Confidence Cluster Comparison — grouped by cell class, Brieflow vs Funk."""
    setup_plot_style()

    stats = {name: get_confidence_stats(df) for name, df in data.items()}

    # Build groups: each cell class gets Brieflow + Funk bars; Shuffled is its own group
    groups = []
    for cc in CELL_CLASSES:
        group_labels = []
        for pipeline in PIPELINES:
            label = f"{pipeline} {cc}"
            if label in stats:
                group_labels.append(label)
        if group_labels:
            groups.append((cc, group_labels))
    if "Shuffled" in stats:
        groups.append(("Shuffled", ["Shuffled"]))

    # Compute bar positions with gaps between groups
    bar_positions = []
    bar_labels = []
    bar_colors = []
    group_centers = []
    group_names = []
    pos = 0
    width = 0.6
    group_gap = 0.4

    for group_name, labels in groups:
        center_positions = []
        for label in labels:
            bar_positions.append(pos)
            bar_labels.append(label)
            if label == "Shuffled":
                bar_colors.append(PIPELINE_COLORS["Shuffled"])
            else:
                pipeline = label.split(" ")[0]
                bar_colors.append(PIPELINE_COLORS[pipeline])
            center_positions.append(pos)
            pos += width + 0.1
        group_centers.append(np.mean(center_positions))
        group_names.append(group_name)
        pos += group_gap

    fig, ax = plt.subplots(figsize=(7, 5))

    pcts = [stats[l]["pct_high"] for l in bar_labels]
    bars = ax.bar(bar_positions, pcts, width, color=bar_colors, edgecolor="white", linewidth=1.5)

    for bar, label in zip(bars, bar_labels):
        h = bar.get_height()
        s = stats[label]
        ax.text(
            bar.get_x() + bar.get_width() / 2, h + 0.8,
            f"{h:.1f}%\n({s['high']}/{s['total']})",
            ha="center", va="bottom", fontsize=9, fontweight="bold",
        )

    ax.set_ylabel("% High Confidence Clusters")
    ax.set_xticks(group_centers)
    ax.set_xticklabels(group_names)
    ax.set_ylim(0, max(pcts) * 1.35)
    ax.set_title("Cluster Pathway Confidence")

    # Legend for pipelines
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=PIPELINE_COLORS[p], label=p) for p in PIPELINES]
    if "Shuffled" in stats:
        handles.append(Patch(facecolor=PIPELINE_COLORS["Shuffled"], label="Shuffled"))
    ax.legend(handles=handles, loc="upper right")

    save_figure(fig, output_dir / "figure_a_high_confidence.png")
    plt.close(fig)
    print("Saved: figure_a_high_confidence.png")


def figure_b(data, output_dir):
    """Figure B: Gene Discovery Potential (stacked bar, high-conf only, grouped by cell class)."""
    setup_plot_style()

    # Exclude Shuffled
    labels = [l for l in data if l != "Shuffled"]
    stats = {l: get_confidence_stats(data[l]) for l in labels}

    # Order: Brieflow Interphase, Funk Interphase, Brieflow Mitotic, Funk Mitotic
    ordered = []
    for cc in CELL_CLASSES:
        for pipeline in PIPELINES:
            label = f"{pipeline} {cc}"
            if label in stats:
                ordered.append(label)

    fig, ax = plt.subplots(figsize=(7, 5))

    # Positions with gap between cell classes
    positions = []
    width = 0.6
    pos = 0
    group_centers = []
    cc_idx = 0
    for i, label in enumerate(ordered):
        positions.append(pos)
        # Add gap after every 2 bars (between cell classes)
        if i > 0 and i % len(PIPELINES) == 0:
            pass  # gap already added below
        pos += width + 0.1
        if (i + 1) % len(PIPELINES) == 0:
            group_centers.append(np.mean(positions[-len(PIPELINES):]))
            pos += 0.3  # gap between groups

    x = np.array(positions)

    established = [stats[l]["genes_established"] for l in ordered]
    novel = [stats[l]["genes_novel"] for l in ordered]
    uncharacterized = [stats[l]["genes_uncharacterized"] for l in ordered]

    # Stack colors: consistent across all bars
    stack_colors = {"Established": "#2ca02c", "Novel role": "#ff7f0e", "Uncharacterized": "#d62728"}

    # Use hatching to distinguish pipelines
    hatches = [None if l.startswith("Brieflow") else "///" for l in ordered]

    for i, label in enumerate(ordered):
        h = hatches[i]
        ax.bar(x[i], established[i], width, color=stack_colors["Established"],
               edgecolor="white", hatch=h)
        ax.bar(x[i], novel[i], width, bottom=established[i], color=stack_colors["Novel role"],
               edgecolor="white", hatch=h)
        ax.bar(x[i], uncharacterized[i], width, bottom=established[i] + novel[i],
               color=stack_colors["Uncharacterized"], edgecolor="white", hatch=h)

    # Total labels
    totals = [stats[l]["genes_total_highconf"] for l in ordered]
    for i, total in enumerate(totals):
        ax.text(x[i], total + 10, f"{total}", ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_ylabel("Genes in High Confidence Clusters")
    ax.set_xticks(group_centers)
    ax.set_xticklabels(CELL_CLASSES)
    ax.set_ylim(0, max(totals) * 1.2)
    ax.set_title("Gene Discovery Potential")

    # Legend: gene categories + pipeline distinction
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=stack_colors["Established"], label="Established"),
        Patch(facecolor=stack_colors["Novel role"], label="Novel role"),
        Patch(facecolor=stack_colors["Uncharacterized"], label="Uncharacterized"),
        Patch(facecolor="#cccccc", edgecolor="black", label="Brieflow"),
        Patch(facecolor="#cccccc", edgecolor="black", hatch="///", label="Funk"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=9)

    save_figure(fig, output_dir / "figure_b_gene_discovery.png")
    plt.close(fig)
    print("Saved: figure_b_gene_discovery.png")


def figure_c(data, output_dir):
    """Figure C: Mean Genes per High Confidence Cluster (grouped bar by pipeline, per gene category)."""
    setup_plot_style()

    # Exclude Shuffled
    labels = [l for l in data if l != "Shuffled"]
    stats = {l: get_confidence_stats(data[l]) for l in labels}

    # Order labels
    ordered = []
    for cc in CELL_CLASSES:
        for pipeline in PIPELINES:
            label = f"{pipeline} {cc}"
            if label in stats:
                ordered.append(label)

    categories = ["Established", "Novel role", "Uncharacterized"]
    cat_keys = ["mean_established", "mean_novel", "mean_uncharacterized"]

    fig, ax = plt.subplots(figsize=(8, 5))

    x = np.arange(len(categories))
    n_groups = len(ordered)
    width = 0.18
    offsets = np.linspace(-(n_groups - 1) * width / 2, (n_groups - 1) * width / 2, n_groups)

    for i, label in enumerate(ordered):
        pipeline = label.split(" ")[0]
        cell_class = " ".join(label.split(" ")[1:])
        vals = [stats[label][k] for k in cat_keys]
        # Use pipeline color, with hatching for Mitotic to distinguish
        hatch = "//" if cell_class == "Mitotic" else None
        bars = ax.bar(
            x + offsets[i], vals, width,
            label=label, color=PIPELINE_COLORS[pipeline],
            edgecolor="white", linewidth=1.5, hatch=hatch, alpha=0.85,
        )
        for bar in bars:
            h = bar.get_height()
            if h > 0:
                ax.text(
                    bar.get_x() + bar.get_width() / 2, h + 0.2,
                    f"{h:.1f}", ha="center", va="bottom", fontsize=8,
                )

    ax.set_ylabel("Mean Genes per High Confidence Cluster")
    ax.set_xticks(x)
    ax.set_xticklabels(categories)
    ax.set_title("Cluster Composition")
    ax.legend(fontsize=9)

    save_figure(fig, output_dir / "figure_c_mean_genes.png")
    plt.close(fig)
    print("Saved: figure_c_mean_genes.png")


def save_summary_table(data, output_dir):
    """Save a summary table as TSV."""
    rows = []
    for name, df in data.items():
        s = get_confidence_stats(df)
        rows.append({
            "condition": name,
            "total_clusters": s["total"],
            "high_confidence": s["high"],
            "pct_high": round(s["pct_high"], 1),
            "medium_confidence": s["medium"],
            "pct_medium": round(s["pct_medium"], 1),
            "low_confidence": s["low"],
            "pct_low": round(s["pct_low"], 1),
            "genes_in_high_conf": s["genes_total_highconf"],
            "genes_established": s["genes_established"],
            "genes_novel": s["genes_novel"],
            "genes_uncharacterized": s["genes_uncharacterized"],
            "mean_established": round(s["mean_established"], 1),
            "mean_novel": round(s["mean_novel"], 1),
            "mean_uncharacterized": round(s["mean_uncharacterized"], 1),
        })
    table = pd.DataFrame(rows)
    path = output_dir / "summary_table.tsv"
    table.to_csv(path, sep="\t", index=False)
    print(f"Saved: summary_table.tsv")
    print(table.to_string(index=False))
    return table


def main():
    parser = argparse.ArgumentParser(description="Generate MozzareLLM manuscript figures A-C")
    parser.add_argument("--include-shuffled", action="store_true", help="Include shuffled negative control in Figure A")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading MozzareLLM summaries...")
    data = load_summaries(include_shuffled=args.include_shuffled)

    if len(data) == 0:
        print("ERROR: No summaries found. Are the mozzarellm runs complete?")
        return

    save_summary_table(data, OUTPUT_DIR)

    print("\nGenerating figures...")
    figure_a(data, OUTPUT_DIR)
    figure_b(data, OUTPUT_DIR)
    figure_c(data, OUTPUT_DIR)

    print(f"\nAll figures saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
