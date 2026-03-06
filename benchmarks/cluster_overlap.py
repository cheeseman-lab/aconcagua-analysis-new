"""
Cluster Overlap Metrics: Brieflow vs Funk

Computes systematic Jaccard similarity and overlap metrics between all
high-confidence clusters in Brieflow and Funk pipelines, in both directions.

Outputs:
  - TSV files with per-cluster overlap metrics
  - Markdown summary tables
  - Heatmap visualizations of top cluster matches

Usage:
    python cluster_overlap.py
    python cluster_overlap.py --min-confidence Medium  # include Medium+High
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path

from plot_style import setup_plot_style, save_figure, COLORS

# Paths
SCRIPT_DIR = Path(__file__).parent
ANALYSIS_DIR = SCRIPT_DIR.parent / "analysis"
CLUSTER_DIR = ANALYSIS_DIR / "brieflow_output" / "cluster" / "DAPI_TUBULIN_GH2AX_PHALLOIDIN"
FUNK_DIR = SCRIPT_DIR / "external" / "results" / "cluster"
OUTPUT_DIR = SCRIPT_DIR / "results" / "cluster" / "overlap"

# Configurations
CONFIGS = {
    "Interphase": {
        "brieflow": {
            "combined": CLUSTER_DIR / "Interphase" / "12" / "CB-Real__combined_table.tsv",
            "summaries": CLUSTER_DIR / "Interphase" / "12" / "mozzarellm" / "claude-opus-4-6_summaries.tsv",
            "flagged": CLUSTER_DIR / "Interphase" / "12" / "mozzarellm" / "claude-opus-4-6_flagged_genes.tsv",
        },
        "funk": {
            "combined": FUNK_DIR / "Funk_Interphase_k10" / "CB-Real__combined_table.tsv",
            "summaries": FUNK_DIR / "Funk_Interphase_k10" / "mozzarellm" / "claude-opus-4-6_summaries.tsv",
            "flagged": FUNK_DIR / "Funk_Interphase_k10" / "mozzarellm" / "claude-opus-4-6_flagged_genes.tsv",
        },
    },
    "Mitotic": {
        "brieflow": {
            "combined": CLUSTER_DIR / "Mitotic" / "5" / "CB-Real__combined_table.tsv",
            "summaries": CLUSTER_DIR / "Mitotic" / "5" / "mozzarellm" / "claude-opus-4-6_summaries.tsv",
            "flagged": CLUSTER_DIR / "Mitotic" / "5" / "mozzarellm" / "claude-opus-4-6_flagged_genes.tsv",
        },
        "funk": {
            "combined": FUNK_DIR / "Funk_Mitotic_k9" / "CB-Real__combined_table.tsv",
            "summaries": FUNK_DIR / "Funk_Mitotic_k9" / "mozzarellm" / "claude-opus-4-6_summaries.tsv",
            "flagged": FUNK_DIR / "Funk_Mitotic_k9" / "mozzarellm" / "claude-opus-4-6_flagged_genes.tsv",
        },
    },
}


def load_pipeline_data(config):
    """Load cluster gene sets, MozzareLLM summaries, and gene categories."""
    combined = pd.read_csv(config["combined"], sep="\t")
    summaries = pd.read_csv(config["summaries"], sep="\t")

    # Parse gene sets from semicolon-separated strings
    cluster_genes = {}
    for _, row in combined.iterrows():
        cl_id = int(row["cluster"])
        genes = {g.strip() for g in str(row["genes"]).split(";") if g.strip()}
        # Remove nontargeting controls
        genes = {g for g in genes if not g.startswith("nontargeting")}
        cluster_genes[cl_id] = genes

    # Build confidence lookup
    confidence = {}
    process = {}
    for _, row in summaries.iterrows():
        cl_id = int(row["cluster_id"])
        confidence[cl_id] = row["pathway_confidence"]
        process[cl_id] = row["dominant_process"]

    # Build gene category lookup: {cluster_id: {"novel_role": set, "uncharacterized": set}}
    gene_categories = {}
    flagged = pd.read_csv(config["flagged"], sep="\t")
    for _, row in flagged.iterrows():
        cl_id = int(row["cluster_id"])
        if cl_id not in gene_categories:
            gene_categories[cl_id] = {"novel_role": set(), "uncharacterized": set()}
        cat = row["category"]
        if cat in ("novel_role", "uncharacterized"):
            gene_categories[cl_id][cat].add(row["gene"])

    return cluster_genes, confidence, process, gene_categories


def compute_overlap(source_genes, target_genes):
    """Compute Jaccard and overlap metrics between two gene sets."""
    intersection = source_genes & target_genes
    union = source_genes | target_genes
    jaccard = len(intersection) / len(union) if union else 0
    # Overlap coefficient: fraction of source genes found in target
    overlap_frac = len(intersection) / len(source_genes) if source_genes else 0
    return {
        "jaccard": jaccard,
        "overlap_frac": overlap_frac,
        "shared_genes": len(intersection),
        "source_size": len(source_genes),
        "target_size": len(target_genes),
        "union_size": len(union),
    }


def find_best_matches(source_clusters, source_conf, source_proc,
                      target_clusters, target_conf, target_proc,
                      source_gene_categories=None,
                      min_confidence="High"):
    """For each source cluster, find best matching target cluster by Jaccard.

    Returns DataFrame with one row per source cluster.
    """
    conf_levels = ["High"] if min_confidence == "High" else ["High", "Medium"]
    if source_gene_categories is None:
        source_gene_categories = {}
    rows = []

    for src_id, src_genes in sorted(source_clusters.items()):
        if source_conf.get(src_id) not in conf_levels:
            continue

        best_jaccard = 0
        best_target = None
        best_metrics = None
        all_targets = []

        for tgt_id, tgt_genes in target_clusters.items():
            metrics = compute_overlap(src_genes, tgt_genes)
            all_targets.append((tgt_id, metrics))
            if metrics["jaccard"] > best_jaccard:
                best_jaccard = metrics["jaccard"]
                best_target = tgt_id
                best_metrics = metrics

        if best_target is None:
            continue

        # Count how many target clusters contain at least 1 source gene
        fragmentation = sum(
            1 for _, m in all_targets if m["shared_genes"] > 0
        )

        # Count total source genes found across ALL target clusters
        all_target_genes = set()
        for tgt_id, tgt_genes in target_clusters.items():
            all_target_genes |= (src_genes & tgt_genes)
        genes_found = len(all_target_genes)
        genes_missing = len(src_genes) - genes_found

        # Classify match quality
        j = best_jaccard
        if j >= 0.6:
            category = "Strong concordance"
        elif j >= 0.3:
            category = "Good concordance"
        elif j >= 0.15:
            category = "Partial match"
        elif j >= 0.05:
            category = "Weak match"
        else:
            category = "Unique to source"

        # Gene categories for this source cluster
        cats = source_gene_categories.get(src_id, {"novel_role": set(), "uncharacterized": set()})
        novel_genes = sorted(cats.get("novel_role", set()))
        unchar_genes = sorted(cats.get("uncharacterized", set()))
        # Established = all genes minus novel minus uncharacterized
        flagged = cats.get("novel_role", set()) | cats.get("uncharacterized", set())
        established_genes = sorted(src_genes - flagged)

        rows.append({
            "source_cluster": src_id,
            "source_process": source_proc.get(src_id, ""),
            "source_confidence": source_conf.get(src_id, ""),
            "source_size": len(src_genes),
            "best_target_cluster": best_target,
            "target_process": target_proc.get(best_target, ""),
            "target_confidence": target_conf.get(best_target, ""),
            "target_size": len(target_clusters[best_target]),
            "jaccard": best_jaccard,
            "overlap_frac": best_metrics["overlap_frac"],
            "shared_genes": best_metrics["shared_genes"],
            "fragmentation": fragmentation,
            "genes_found_total": genes_found,
            "genes_missing": genes_missing,
            "category": category,
            "genes_established": "; ".join(established_genes),
            "genes_novel_role": "; ".join(novel_genes),
            "genes_uncharacterized": "; ".join(unchar_genes),
        })

    return pd.DataFrame(rows)


def compute_gene_level_overlap(source_clusters, source_conf, source_proc,
                               target_clusters, min_confidence="High"):
    """Compute gene-level statistics: how many genes from high-conf source
    clusters appear in ANY target cluster, and vice versa."""
    conf_levels = ["High"] if min_confidence == "High" else ["High", "Medium"]

    # Collect all genes in high-conf source clusters
    source_genes_hc = set()
    for cl_id, genes in source_clusters.items():
        if source_conf.get(cl_id) in conf_levels:
            source_genes_hc |= genes

    # Collect all genes in ANY target cluster
    all_target_genes = set()
    for genes in target_clusters.values():
        all_target_genes |= genes

    shared = source_genes_hc & all_target_genes
    source_unique = source_genes_hc - all_target_genes

    return {
        "source_hc_genes": len(source_genes_hc),
        "in_target": len(shared),
        "source_unique": len(source_unique),
        "pct_found": 100 * len(shared) / len(source_genes_hc) if source_genes_hc else 0,
    }


def plot_jaccard_distribution(results_dict, output_path):
    """Plot Jaccard similarity distribution for each condition."""
    setup_plot_style()
    n = len(results_dict)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.5), squeeze=False)

    for idx, (label, df) in enumerate(results_dict.items()):
        ax = axes[0, idx]
        jaccards = df["jaccard"].values

        bins = np.arange(0, 1.05, 0.05)
        ax.hist(jaccards, bins=bins, color=COLORS["brieflow"], edgecolor="white",
                linewidth=0.5, alpha=0.85)

        ax.set_xlabel("Best-match Jaccard similarity")
        ax.set_ylabel("Number of clusters")
        ax.set_title(label, fontweight="bold")
        ax.set_xlim(0, 1)

        # Add median line
        med = np.median(jaccards)
        ax.axvline(med, color="#d62728", linestyle="--", linewidth=1.5, alpha=0.8)
        ax.text(med + 0.02, ax.get_ylim()[1] * 0.9, f"median={med:.2f}",
                color="#d62728", fontsize=10)

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close(fig)
    print(f"Saved: {output_path}")


def plot_category_summary(results_dict, output_path):
    """Stacked bar chart of match categories across conditions."""
    setup_plot_style()
    categories = [
        "Strong concordance", "Good concordance",
        "Partial match", "Weak match", "Unique to source",
    ]
    cat_colors = ["#2ca02c", "#98df8a", "#ffbb78", "#ff7f0e", "#d62728"]

    labels = list(results_dict.keys())
    n = len(labels)
    fig, ax = plt.subplots(figsize=(max(6, 2.5 * n), 5))

    x = np.arange(n)
    bottom = np.zeros(n)

    for cat, color in zip(categories, cat_colors):
        counts = []
        for label in labels:
            df = results_dict[label]
            counts.append((df["category"] == cat).sum())
        counts = np.array(counts, dtype=float)
        ax.bar(x, counts, bottom=bottom, color=color, edgecolor="white",
               linewidth=0.5, label=cat, width=0.6)
        # Label non-zero segments
        for i, c in enumerate(counts):
            if c > 0:
                ax.text(x[i], bottom[i] + c / 2, str(int(c)),
                        ha="center", va="center", fontsize=10, fontweight="bold")
        bottom += counts

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Number of high-conf clusters")
    ax.set_title("Cluster Match Quality: Brieflow → Funk", fontweight="bold")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=9)
    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close(fig)
    print(f"Saved: {output_path}")


def plot_overlap_heatmap(results_dict, output_path, top_n=20):
    """Heatmap of top Jaccard matches for each condition."""
    setup_plot_style()
    n = len(results_dict)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, max(6, 0.4 * top_n)),
                             squeeze=False)

    cmap = LinearSegmentedColormap.from_list(
        "jaccard", ["#ffffff", "#c7e9c0", "#41ab5d", "#006d2c"])

    for idx, (label, df) in enumerate(results_dict.items()):
        ax = axes[0, idx]
        df_sorted = df.sort_values("jaccard", ascending=False).head(top_n)

        # Build labels
        y_labels = []
        for _, row in df_sorted.iterrows():
            proc = row["source_process"]
            if len(proc) > 40:
                proc = proc[:37] + "..."
            y_labels.append(f"cl{int(row['source_cluster'])} ({proc})")

        data = df_sorted[["jaccard", "overlap_frac"]].values
        im = ax.imshow(data, aspect="auto", cmap=cmap, vmin=0, vmax=1)

        ax.set_yticks(range(len(y_labels)))
        ax.set_yticklabels(y_labels, fontsize=8)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Jaccard", "Overlap\nfraction"], fontsize=10)
        ax.set_title(label, fontweight="bold")

        # Annotate cells
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                val = data[i, j]
                color = "white" if val > 0.6 else "black"
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=8, color=color)

    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, label="Score")
    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Cluster overlap metrics: Brieflow vs Funk")
    parser.add_argument("--min-confidence", default="High",
                        choices=["High", "Medium"],
                        help="Minimum confidence level to include")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    all_bf_to_fk = {}
    all_fk_to_bf = {}
    gene_stats = {}

    for cell_class, pipelines in CONFIGS.items():
        print(f"\n{'='*60}")
        print(f"Processing {cell_class}")
        print(f"{'='*60}")

        # Load data
        bf_genes, bf_conf, bf_proc, bf_cats = load_pipeline_data(pipelines["brieflow"])
        fk_genes, fk_conf, fk_proc, fk_cats = load_pipeline_data(pipelines["funk"])

        print(f"  Brieflow: {len(bf_genes)} clusters")
        print(f"  Funk: {len(fk_genes)} clusters")

        # Brieflow -> Funk
        bf_to_fk = find_best_matches(
            bf_genes, bf_conf, bf_proc,
            fk_genes, fk_conf, fk_proc,
            source_gene_categories=bf_cats,
            min_confidence=args.min_confidence,
        )
        label_bf = f"BF→FK {cell_class}"
        all_bf_to_fk[label_bf] = bf_to_fk
        print(f"  BF→FK: {len(bf_to_fk)} high-conf clusters matched")

        # Funk -> Brieflow
        fk_to_bf = find_best_matches(
            fk_genes, fk_conf, fk_proc,
            bf_genes, bf_conf, bf_proc,
            source_gene_categories=fk_cats,
            min_confidence=args.min_confidence,
        )
        label_fk = f"FK→BF {cell_class}"
        all_fk_to_bf[label_fk] = fk_to_bf
        print(f"  FK→BF: {len(fk_to_bf)} high-conf clusters matched")

        # Save TSVs
        bf_to_fk.to_csv(OUTPUT_DIR / f"bf_to_fk_{cell_class.lower()}.tsv",
                         sep="\t", index=False)
        fk_to_bf.to_csv(OUTPUT_DIR / f"fk_to_bf_{cell_class.lower()}.tsv",
                         sep="\t", index=False)

        # Gene-level stats
        gene_stats[f"BF→FK {cell_class}"] = compute_gene_level_overlap(
            bf_genes, bf_conf, bf_proc, fk_genes,
            min_confidence=args.min_confidence,
        )
        gene_stats[f"FK→BF {cell_class}"] = compute_gene_level_overlap(
            fk_genes, fk_conf, fk_proc, bf_genes,
            min_confidence=args.min_confidence,
        )

    # Combine all results
    all_results = {**all_bf_to_fk, **all_fk_to_bf}

    # Print summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for label, df in all_results.items():
        print(f"\n{label}:")
        print(f"  Clusters: {len(df)}")
        print(f"  Median Jaccard: {df['jaccard'].median():.3f}")
        print(f"  Mean Jaccard: {df['jaccard'].mean():.3f}")
        cat_counts = df["category"].value_counts()
        for cat in ["Strong concordance", "Good concordance", "Partial match",
                     "Weak match", "Unique to source"]:
            c = cat_counts.get(cat, 0)
            print(f"  {cat}: {c}")

    print("\nGene-level:")
    for label, stats in gene_stats.items():
        print(f"  {label}: {stats['source_hc_genes']} HC genes, "
              f"{stats['in_target']} found ({stats['pct_found']:.1f}%), "
              f"{stats['source_unique']} unique")

    # Generate plots
    print("\nGenerating plots...")
    plot_jaccard_distribution(all_bf_to_fk, OUTPUT_DIR / "jaccard_dist_bf_to_fk.png")
    plot_jaccard_distribution(all_fk_to_bf, OUTPUT_DIR / "jaccard_dist_fk_to_bf.png")
    plot_category_summary(
        {**all_bf_to_fk, **all_fk_to_bf},
        OUTPUT_DIR / "category_summary.png",
    )
    plot_overlap_heatmap(all_bf_to_fk, OUTPUT_DIR / "heatmap_bf_to_fk.png")
    plot_overlap_heatmap(all_fk_to_bf, OUTPUT_DIR / "heatmap_fk_to_bf.png")

    print("\nDone! Results in:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
