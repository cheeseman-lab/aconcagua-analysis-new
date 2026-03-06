"""
Cluster PHATE Embedding Visualizations

Generates publication-quality PHATE scatter plots colored by cluster assignment
for both Brieflow and Funk et al. 2022, plus zoomed-in Pearson correlation
heatmaps for user-specified cluster subsets (e.g. mitochondrial sub-modules).

Outputs:
  - PHATE scatter (all clusters): categorical color per cluster
  - PHATE scatter (confidence): High/Medium/Low MozzareLLM confidence
  - PHATE scatter (highlighted): specific cluster subsets overlaid
  - Zoomed heatmaps: correlation structure within highlighted clusters

Usage:
    python cluster_phate.py                # All outputs
    python cluster_phate.py --phate-only   # PHATE scatter plots only
    python cluster_phate.py --heatmap-only # Zoomed heatmaps only
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
from sklearn.metrics.pairwise import cosine_similarity

from plot_style import setup_plot_style, save_figure

# ============================================================================
# Paths
# ============================================================================

BENCHMARKS_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BENCHMARKS_DIR.parent
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
AGGREGATE_DIR = ANALYSIS_DIR / "brieflow_output" / "aggregate" / "tsvs"
CLUSTER_DIR = (
    ANALYSIS_DIR / "brieflow_output" / "cluster" / "DAPI_TUBULIN_GH2AX_PHALLOIDIN"
)
EXTERNAL_DIR = BENCHMARKS_DIR / "external"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "cluster" / "phate"

BRIEFLOW_METADATA_COLS = {
    "gene_symbol_0",
    "cell_count",
    "cell_stage",
    "cell_stage_confidence",
    "perturbation_auc",
}

IDEAL_CONFIGS = {
    "Interphase": {"brieflow_k": 12, "funk_k": 10},
    "Mitotic": {"brieflow_k": 5, "funk_k": 9},
}

CONFIDENCE_COLORS = {
    "High": "#0077BB",    # Teal blue
    "Medium": "#EE7733",  # Coral orange
    "Low": "#cccccc",     # Light grey
}

# ============================================================================
# Cluster highlight definitions
# ============================================================================

HIGHLIGHT_SETS = {
    "Interphase": {
        "mitochondrial": {
            "label": "Mitochondrial",
            "brieflow_clusters": [22, 69, 121, 129, 160],
            "funk_clusters": [135, 147, 149, 185],
            "description": (
                "Brieflow resolves mitochondrial biology into 5 distinct "
                "sub-modules (OXPHOS, mitoribosomes, ATP synthase, translation, "
                "membrane organization) vs Funk's 4 mixed clusters"
            ),
        },
    },
}

# ============================================================================
# Data loading
# ============================================================================


def load_phate_data(cell_class):
    """Load PHATE coordinates and cluster assignments for both pipelines."""
    cfg = IDEAL_CONFIGS[cell_class]

    # Brieflow
    brieflow_path = (
        CLUSTER_DIR
        / cell_class
        / str(cfg["brieflow_k"])
        / "phate_leiden_clustering.tsv"
    )
    brieflow_df = pd.read_csv(brieflow_path, sep="\t")
    brieflow_df = brieflow_df.rename(columns={"gene_symbol_0": "gene_symbol"})

    # Funk
    funk_path = (
        EXTERNAL_DIR
        / "results"
        / "cluster"
        / f"Funk_{cell_class}_k{cfg['funk_k']}"
        / "phate_leiden_clustering.tsv"
    )
    funk_df = pd.read_csv(funk_path, sep="\t")

    return brieflow_df, funk_df, cfg


def _load_mozzarellm_confidence(cell_class):
    """Load MozzareLLM confidence per cluster for both pipelines.

    Returns:
        (brieflow_conf, funk_conf): dicts mapping cluster_id → confidence level
    """
    cfg = IDEAL_CONFIGS[cell_class]

    brieflow_path = (
        CLUSTER_DIR
        / cell_class
        / str(cfg["brieflow_k"])
        / "mozzarellm"
        / "claude-opus-4-6_summaries.tsv"
    )
    funk_path = (
        EXTERNAL_DIR
        / "results"
        / "cluster"
        / f"Funk_{cell_class}_k{cfg['funk_k']}"
        / "mozzarellm"
        / "claude-opus-4-6_summaries.tsv"
    )

    def _load(path):
        df = pd.read_csv(path, sep="\t")
        return dict(zip(df["cluster_id"], df["pathway_confidence"]))

    return _load(brieflow_path), _load(funk_path)


def load_features_for_heatmap(cell_class):
    """Load PCA features for zoomed heatmaps (shared genes only)."""
    cfg = IDEAL_CONFIGS[cell_class]

    # Brieflow aggregated (PCA features)
    agg_path = (
        AGGREGATE_DIR
        / f"CeCl-{cell_class}_ChCo-DAPI_TUBULIN_GH2AX_PHALLOIDIN__aggregated.tsv"
    )
    brieflow_agg = pd.read_csv(agg_path, sep="\t")
    brieflow_pc_cols = [c for c in brieflow_agg.columns if c.startswith("PC_")]

    # Funk PCA features
    funk_pca_path = EXTERNAL_DIR / f"Funk_{cell_class}_pca_features.tsv"
    funk_pca = pd.read_csv(funk_pca_path, sep="\t")
    funk_pc_cols = [c for c in funk_pca.columns if c.startswith("PC_")]

    n_pcs = len(funk_pc_cols)
    brieflow_pc_matched = brieflow_pc_cols[:n_pcs]

    # Brieflow cluster assignments
    brieflow_clust_path = (
        CLUSTER_DIR
        / cell_class
        / str(cfg["brieflow_k"])
        / "phate_leiden_clustering.tsv"
    )
    brieflow_clust = pd.read_csv(brieflow_clust_path, sep="\t")[
        ["gene_symbol_0", "cluster"]
    ].rename(columns={"gene_symbol_0": "gene_symbol", "cluster": "brieflow_cluster"})

    # Funk cluster assignments
    funk_clust_path = (
        EXTERNAL_DIR
        / "results"
        / "cluster"
        / f"Funk_{cell_class}_k{cfg['funk_k']}"
        / "phate_leiden_clustering.tsv"
    )
    funk_clust = pd.read_csv(funk_clust_path, sep="\t")[
        ["gene_symbol", "cluster"]
    ].rename(columns={"cluster": "funk_cluster"})

    # Shared genes
    shared_genes = (
        set(brieflow_agg["gene_symbol_0"])
        & set(funk_pca["gene_symbol"])
        & set(brieflow_clust["gene_symbol"])
        & set(funk_clust["gene_symbol"])
    )
    gene_list = sorted(shared_genes)

    def _subset(df, gene_col, cols):
        sub = df[df[gene_col].isin(shared_genes)].drop_duplicates(gene_col)
        sub = sub.set_index(gene_col).loc[gene_list]
        return sub[cols].values

    brieflow_pca_arr = _subset(brieflow_agg, "gene_symbol_0", brieflow_pc_matched)
    funk_pca_arr = _subset(funk_pca, "gene_symbol", funk_pc_cols)

    # Merge cluster assignments
    clust_df = brieflow_clust[brieflow_clust["gene_symbol"].isin(shared_genes)].merge(
        funk_clust[funk_clust["gene_symbol"].isin(shared_genes)],
        on="gene_symbol",
        how="inner",
    )
    clust_df = clust_df.drop_duplicates("gene_symbol").set_index("gene_symbol")
    clust_df = clust_df.loc[gene_list].reset_index()

    return {
        "gene_list": gene_list,
        "brieflow_pca": brieflow_pca_arr,
        "funk_pca": funk_pca_arr,
        "n_pcs": n_pcs,
        "brieflow_cluster": clust_df["brieflow_cluster"].values,
        "funk_cluster": clust_df["funk_cluster"].values,
        "cfg": cfg,
    }


# ============================================================================
# PHATE scatter plots
# ============================================================================


def _cluster_colormap(n_clusters):
    """Generate a categorical colormap for n clusters."""
    if n_clusters <= 20:
        base = plt.colormaps["tab20"]
        return [base(i) for i in range(n_clusters)]
    base_colors = (
        list(plt.colormaps["tab20"](np.linspace(0, 1, 20)))
        + list(plt.colormaps["tab20b"](np.linspace(0, 1, 20)))
        + list(plt.colormaps["tab20c"](np.linspace(0, 1, 20)))
    )
    return [base_colors[i % len(base_colors)] for i in range(n_clusters)]


def _square_scatter_axes(ax):
    """Force a scatter subplot to be square with equal data range on both axes."""
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])


def plot_phate_scatter(cell_class, output_path, highlight_set=None):
    """Plot side-by-side PHATE embeddings for Brieflow vs Funk."""
    setup_plot_style()
    brieflow_df, funk_df, cfg = load_phate_data(cell_class)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    point_size = 4

    for ax, df, pipeline, k in [
        (axes[0], brieflow_df, "Brieflow", cfg["brieflow_k"]),
        (axes[1], funk_df, "Funk", cfg["funk_k"]),
    ]:
        clusters = df["cluster"].values
        unique_clusters = sorted(df["cluster"].unique())
        n_clusters = len(unique_clusters)
        cluster_to_idx = {c: i for i, c in enumerate(unique_clusters)}
        colors = _cluster_colormap(n_clusters)

        if highlight_set is not None:
            key = "brieflow_clusters" if pipeline == "Brieflow" else "funk_clusters"
            hl_clusters = set(highlight_set[key])

            bg_mask = ~df["cluster"].isin(hl_clusters)
            ax.scatter(
                df.loc[bg_mask, "PHATE_0"],
                df.loc[bg_mask, "PHATE_1"],
                c="#cccccc",
                s=point_size,
                alpha=0.3,
                rasterized=True,
                linewidths=0,
            )

            legend_handles = []
            hl_colors = _cluster_colormap(len(hl_clusters))
            for i, cl in enumerate(sorted(hl_clusters)):
                mask = df["cluster"] == cl
                if mask.sum() == 0:
                    continue
                ax.scatter(
                    df.loc[mask, "PHATE_0"],
                    df.loc[mask, "PHATE_1"],
                    c=[hl_colors[i]],
                    s=point_size * 2.5,
                    alpha=0.9,
                    rasterized=True,
                    linewidths=0.3,
                    edgecolors="black",
                    zorder=3,
                )
                legend_handles.append(
                    mpatches.Patch(
                        color=hl_colors[i],
                        label=f"Cluster {cl} ({mask.sum()} genes)",
                    )
                )
            ax.legend(
                handles=legend_handles,
                loc="upper right",
                fontsize=8,
                framealpha=0.9,
                markerscale=0.8,
            )
        else:
            c_arr = [colors[cluster_to_idx[c]] for c in clusters]
            ax.scatter(
                df["PHATE_0"],
                df["PHATE_1"],
                c=c_arr,
                s=point_size,
                alpha=0.6,
                rasterized=True,
                linewidths=0,
            )

        n_genes = len(df)
        ax.set_xlabel("PHATE 1")
        ax.set_ylabel("PHATE 2")
        ax.set_title(
            f"{pipeline} (k={k}, {n_clusters} clusters, {n_genes} genes)",
            fontweight="bold",
        )
        _square_scatter_axes(ax)

    hl_suffix = ""
    if highlight_set is not None:
        hl_suffix = f" — {highlight_set['label']} clusters highlighted"
    fig.suptitle(
        f"{cell_class} PHATE Embeddings{hl_suffix}",
        fontsize=16,
        fontweight="bold",
    )
    fig.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_phate_confidence(cell_class, output_path):
    """Plot PHATE embeddings colored by MozzareLLM pathway confidence."""
    setup_plot_style()
    brieflow_df, funk_df, cfg = load_phate_data(cell_class)
    bl_conf, fk_conf = _load_mozzarellm_confidence(cell_class)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    point_size = 4

    for ax, df, pipeline, conf_map in [
        (axes[0], brieflow_df, "Brieflow", bl_conf),
        (axes[1], funk_df, "Funk", fk_conf),
    ]:
        # Map each gene to its cluster's confidence level
        conf_levels = df["cluster"].map(conf_map).fillna("Low").values

        # Draw in order: Low (back), Medium, High (front)
        for level in ["Low", "Medium", "High"]:
            mask = conf_levels == level
            if mask.sum() == 0:
                continue
            alpha = 0.2 if level == "Low" else (0.6 if level == "Medium" else 0.85)
            size = point_size if level == "Low" else (
                point_size * 1.2 if level == "Medium" else point_size * 1.5
            )
            ax.scatter(
                df.loc[mask, "PHATE_0"],
                df.loc[mask, "PHATE_1"],
                c=CONFIDENCE_COLORS[level],
                s=size,
                alpha=alpha,
                rasterized=True,
                linewidths=0,
                zorder=1 if level == "Low" else (2 if level == "Medium" else 3),
            )

        # Count clusters per level
        cluster_conf = pd.Series(conf_map)
        counts = cluster_conf.value_counts()
        legend_handles = []
        for level in ["High", "Medium", "Low"]:
            n = counts.get(level, 0)
            legend_handles.append(
                mpatches.Patch(
                    color=CONFIDENCE_COLORS[level],
                    label=f"{level} ({n} clusters)",
                )
            )
        ax.legend(
            handles=legend_handles,
            loc="upper right",
            fontsize=9,
            framealpha=0.9,
            title="Confidence",
            title_fontsize=9,
        )

        k = cfg["brieflow_k"] if pipeline == "Brieflow" else cfg["funk_k"]
        ax.set_xlabel("PHATE 1")
        ax.set_ylabel("PHATE 2")
        ax.set_title(
            f"{pipeline} (k={k})",
            fontweight="bold",
        )
        _square_scatter_axes(ax)

    fig.suptitle(
        f"{cell_class} PHATE — MozzareLLM Pathway Confidence",
        fontsize=16,
        fontweight="bold",
    )
    fig.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"  Saved: {output_path}")


# ============================================================================
# Zoomed-in heatmaps
# ============================================================================


def _center(X):
    """Column-center so cosine_similarity becomes Pearson correlation."""
    return X - X.mean(axis=0, keepdims=True)


def plot_zoomed_heatmap(cell_class, highlight_set, output_path):
    """Plot zoomed Pearson correlation heatmap for highlighted clusters.

    Shows a side-by-side comparison:
      Left:  Brieflow highlighted clusters (genes ordered by cluster)
      Right: Funk highlighted clusters (same gene universe, Funk ordering)

    Only genes belonging to the highlighted clusters (union) are shown.
    """
    setup_plot_style()
    data = load_features_for_heatmap(cell_class)

    bl_clusters = data["brieflow_cluster"]
    fk_clusters = data["funk_cluster"]
    brieflow_pca = data["brieflow_pca"]
    funk_pca = data["funk_pca"]

    hl_bl = set(highlight_set["brieflow_clusters"])
    hl_fk = set(highlight_set["funk_clusters"])

    bl_mask = np.isin(bl_clusters, list(hl_bl))
    fk_mask = np.isin(fk_clusters, list(hl_fk))
    union_mask = bl_mask | fk_mask

    n_union = union_mask.sum()
    print(
        f"  {highlight_set['label']}: {bl_mask.sum()} genes in Brieflow clusters, "
        f"{fk_mask.sum()} in Funk clusters, {n_union} union"
    )

    sub_bl_pca = _center(brieflow_pca[union_mask])
    sub_fk_pca = _center(funk_pca[union_mask])
    sub_bl_cl = bl_clusters[union_mask]
    sub_fk_cl = fk_clusters[union_mask]

    corr_bl = cosine_similarity(sub_bl_pca)
    corr_fk = cosine_similarity(sub_fk_pca)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    for ax, corr, clusters, hl_set, pipeline in [
        (axes[0], corr_bl, sub_bl_cl, hl_bl, "Brieflow"),
        (axes[1], corr_fk, sub_fk_cl, hl_fk, "Funk"),
    ]:
        # Order: highlighted clusters first (sorted), then rest
        hl_mask_local = np.isin(clusters, list(hl_set))
        hl_idx = np.where(hl_mask_local)[0]
        other_idx = np.where(~hl_mask_local)[0]

        hl_order = hl_idx[np.argsort(clusters[hl_idx], kind="stable")]
        other_order = other_idx[np.argsort(clusters[other_idx], kind="stable")]
        order = np.concatenate([hl_order, other_order])

        corr_ordered = corr[np.ix_(order, order)]
        clusters_ordered = clusters[order]
        n_hl = len(hl_order)
        n_total = len(order)

        im = ax.pcolormesh(
            corr_ordered,
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            rasterized=True,
        )
        ax.set_aspect("equal")
        ax.invert_yaxis()

        # Cluster boundary lines and labels within highlighted region
        unique_hl_clusters = []
        seen = set()
        for c in clusters_ordered[:n_hl]:
            if c not in seen:
                unique_hl_clusters.append(c)
                seen.add(c)

        pos = 0
        for cl in unique_hl_clusters:
            cl_size = (clusters_ordered[:n_hl] == cl).sum()
            ax.axhline(
                y=pos, xmin=0, xmax=n_hl / n_total, color="black", linewidth=0.8
            )
            ax.axvline(
                x=pos, ymin=1 - n_hl / n_total, ymax=1, color="black", linewidth=0.8
            )
            mid = pos + cl_size / 2
            ax.text(
                mid,
                mid,
                str(cl),
                ha="center",
                va="center",
                fontsize=7,
                fontweight="bold",
                color="black",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.7, lw=0),
            )
            pos += cl_size

        # Boundary between highlighted and context genes
        ax.axhline(y=n_hl, xmin=0, xmax=1, color="black", linewidth=1.2, alpha=0.5)
        ax.axvline(x=n_hl, ymin=0, ymax=1, color="black", linewidth=1.2, alpha=0.5)

        n_hl_clusters = len(unique_hl_clusters)
        ax.set_title(
            f"{pipeline} — {n_hl_clusters} {highlight_set['label'].lower()} "
            f"clusters ({n_hl} genes), {n_total - n_hl} context genes",
            fontsize=10,
            fontweight="bold",
        )
        ax.set_xticks([])
        ax.set_yticks([])

    fig.colorbar(im, ax=axes, label="Pearson r", shrink=0.6, pad=0.02)
    fig.suptitle(
        f"{cell_class} — {highlight_set['label']} Cluster Correlation (PCA features)",
        fontsize=14,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 0.92, 0.95])
    save_figure(fig, output_path)
    plt.close()
    print(f"  Saved: {output_path}")


# ============================================================================
# Main
# ============================================================================


def main():
    parser = argparse.ArgumentParser(description="Cluster PHATE visualizations")
    parser.add_argument(
        "--phate-only", action="store_true", help="Only generate PHATE scatter plots"
    )
    parser.add_argument(
        "--heatmap-only", action="store_true", help="Only generate zoomed heatmaps"
    )
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    do_phate = not args.heatmap_only
    do_heatmap = not args.phate_only

    print("=" * 70)
    print("CLUSTER PHATE VISUALIZATIONS")
    print("=" * 70)

    for cell_class in IDEAL_CONFIGS:
        print(f"\n--- {cell_class} ---")

        if do_phate:
            # Full PHATE scatter (all clusters, categorical colors)
            plot_phate_scatter(
                cell_class,
                OUTPUT_DIR / f"phate_scatter_{cell_class.lower()}.png",
            )

            # Confidence-colored PHATE scatter
            plot_phate_confidence(
                cell_class,
                OUTPUT_DIR / f"phate_confidence_{cell_class.lower()}.png",
            )

            # Highlighted PHATE scatter for each highlight set
            if cell_class in HIGHLIGHT_SETS:
                for key, hl_set in HIGHLIGHT_SETS[cell_class].items():
                    plot_phate_scatter(
                        cell_class,
                        OUTPUT_DIR / f"phate_scatter_{cell_class.lower()}_{key}.png",
                        highlight_set=hl_set,
                    )

        if do_heatmap:
            if cell_class in HIGHLIGHT_SETS:
                for key, hl_set in HIGHLIGHT_SETS[cell_class].items():
                    plot_zoomed_heatmap(
                        cell_class,
                        hl_set,
                        OUTPUT_DIR / f"zoomed_heatmap_{cell_class.lower()}_{key}.png",
                    )

    print(f"\nOutputs saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
