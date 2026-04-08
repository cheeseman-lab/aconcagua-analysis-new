"""
Cluster Overlap: Brieflow vs Funk Pipeline-Unique Clusters

Computes bidirectional Jaccard similarity between all high-confidence clusters
in Brieflow and Funk pipelines. Identifies pipeline-unique clusters (Jaccard
< 0.15) and annotates them with novel-role and uncharacterized gene members.

Outputs:
    results/cluster/overlap/bf_to_fk_interphase.tsv   Brieflow→Funk interphase matches
    results/cluster/overlap/bf_to_fk_mitotic.tsv       Brieflow→Funk mitotic matches
    results/cluster/overlap/fk_to_bf_interphase.tsv    Funk→Brieflow interphase matches
    results/cluster/overlap/fk_to_bf_mitotic.tsv       Funk→Brieflow mitotic matches
    results/cluster/overlap/unique_clusters.png         Combined unique clusters figure

Usage:
    python cluster_overlap.py
    python cluster_overlap.py --min-confidence Medium  # include Medium+High
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


def _shorten_process(proc, maxlen=55):
    """Truncate process name for bar labels."""
    if len(proc) > maxlen:
        return proc[: maxlen - 3] + "..."
    return proc


def _format_gene_list(genes_str, max_genes=6):
    """Format semicolon-separated gene string, truncating if needed."""
    if not genes_str or pd.isna(genes_str) or genes_str.strip() == "":
        return ""
    genes = [g.strip() for g in str(genes_str).split(";") if g.strip()]
    if not genes:
        return ""
    if len(genes) <= max_genes:
        return ", ".join(genes)
    return ", ".join(genes[:max_genes]) + f" (+{len(genes) - max_genes} more)"


def _collect_low_jaccard(results_dict, jaccard_threshold=0.15):
    """Collect low-Jaccard clusters from results dict."""
    rows = []
    for label, df in results_dict.items():
        cc = label.split()[-1]
        low = df[df["jaccard"] < jaccard_threshold].copy()
        low["cell_class"] = cc
        rows.append(low)
    if not rows:
        return pd.DataFrame()
    combined = pd.concat(rows, ignore_index=True)
    return combined.sort_values("jaccard", ascending=True).reset_index(drop=True)


def plot_unique_combined(bf_to_fk_dict, fk_to_bf_dict, output_path,
                         jaccard_threshold=0.15):
    """Single combined figure: Brieflow-unique clusters on top, Funk-unique below."""
    setup_plot_style()

    bf_low = _collect_low_jaccard(bf_to_fk_dict, jaccard_threshold)
    fk_low = _collect_low_jaccard(fk_to_bf_dict, jaccard_threshold)
    n_bf = len(bf_low)
    n_fk = len(fk_low)
    gap = 1  # small space between groups

    fig, ax = plt.subplots(figsize=(9, 11))

    # Brieflow on top (higher y values), Funk on bottom
    # Funk: y = 0..n_fk-1, gap, Brieflow: y = n_fk+gap..n_fk+gap+n_bf-1
    fk_y = np.arange(n_fk)
    bf_y = np.arange(n_bf) + n_fk + gap

    # --- Draw bars ---
    for i, (_, row) in enumerate(fk_low.iterrows()):
        is_mitotic = row["cell_class"] == "Mitotic"
        ax.barh(fk_y[i], row["jaccard"], color=COLORS["funk"], alpha=0.85,
                height=0.65, hatch="//" if is_mitotic else None,
                edgecolor="white", linewidth=1.5)

    for i, (_, row) in enumerate(bf_low.iterrows()):
        is_mitotic = row["cell_class"] == "Mitotic"
        ax.barh(bf_y[i], row["jaccard"], color=COLORS["brieflow"], alpha=0.85,
                height=0.65, hatch="//" if is_mitotic else None,
                edgecolor="white", linewidth=1.5)


    # --- Labels ---
    def _add_labels(df, y_positions):
        for i, (_, row) in enumerate(df.iterrows()):
            proc = _shorten_process(row["source_process"])
            cc = row["cell_class"]
            suffix = " (mitotic)" if cc == "Mitotic" else ""
            n_genes = int(row["source_size"])
            j = row["jaccard"]

            label = f"{proc}{suffix} ({n_genes})"
            ax.text(j + 0.002, y_positions[i], label, va="center",
                    ha="left", fontsize=8.5, fontweight="bold")

    _add_labels(bf_low, bf_y)
    _add_labels(fk_low, fk_y)

    ax.set_yticks([])
    ax.set_xlabel("Best-match Jaccard similarity", fontsize=10, fontfamily="Arial")
    ax.set_xlim(0, 0.225)
    ax.set_ylim(-0.8, n_fk + gap + n_bf - 0.2)
    ax.tick_params(axis='x', labelsize=9)
    ax.spines["left"].set_visible(False)

    plt.tight_layout()
    save_figure(fig, output_path, dpi=400)
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Cluster overlap metrics: Brieflow vs Funk")
    parser.add_argument("--min-confidence", default="High",
                        choices=["High", "Medium"],
                        help="Minimum confidence level to include")
    parser.add_argument("--plots-only", action="store_true",
                        help="Regenerate figures from cached TSVs without recomputing")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if args.plots_only:
        print("Loading cached overlap TSVs...")
        all_bf_to_fk = {
            "BF→FK Interphase": pd.read_csv(OUTPUT_DIR / "bf_to_fk_interphase.tsv", sep="\t"),
            "BF→FK Mitotic":    pd.read_csv(OUTPUT_DIR / "bf_to_fk_mitotic.tsv", sep="\t"),
        }
        all_fk_to_bf = {
            "FK→BF Interphase": pd.read_csv(OUTPUT_DIR / "fk_to_bf_interphase.tsv", sep="\t"),
            "FK→BF Mitotic":    pd.read_csv(OUTPUT_DIR / "fk_to_bf_mitotic.tsv", sep="\t"),
        }
        plot_unique_combined(all_bf_to_fk, all_fk_to_bf, OUTPUT_DIR / "unique_clusters.png")
        print("Done. Results in:", OUTPUT_DIR)
        return

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
    plot_unique_combined(
        all_bf_to_fk, all_fk_to_bf, OUTPUT_DIR / "unique_clusters.png",
    )

    print("\nDone! Results in:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
