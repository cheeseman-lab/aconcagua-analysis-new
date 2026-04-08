"""
Cluster PHATE Embedding Visualizations

Generates publication-quality PHATE scatter plots colored by cluster assignment
for both Brieflow and Funk et al. 2022, plus split-diagonal Pearson correlation
heatmaps for user-specified cluster subsets (e.g. mitochondrial sub-modules).

Outputs:
  - PHATE scatter (all clusters): continuous color by spatial position
  - PHATE scatter (confidence): High/Medium/Low MozzareLLM confidence
  - PHATE scatter (highlighted): specific cluster subsets overlaid
  - Split-diagonal heatmap: upper=Brieflow, lower=Funk, MozzareLLM labels

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

IDEAL_CONFIGS = {
    "Interphase": {"brieflow_k": 12, "funk_k": 10},
    "Mitotic": {"brieflow_k": 5, "funk_k": 9},
}

CONFIDENCE_COLORS = {
    "High": "#0077BB",
    "Medium": "#EE7733",
    "Low": "#cccccc",
}

BRIEFLOW_METADATA_COLS = {
    "gene_symbol_0",
    "cell_count",
    "cell_stage",
    "cell_stage_confidence",
    "perturbation_auc",
}

# Unified cluster palettes — used by PHATE scatter highlights AND heatmaps
# Brieflow: warm tones (orange, magenta-pink, olive gold, soft red, amber)
BRIEFLOW_CLUSTER_PALETTE = ["#E0833A", "#C75B8A", "#B8A038", "#D45D5D", "#CF8B2E"]
# Funk: cool tones (steel blue, teal, slate indigo, mid blue, forest green)
FUNK_CLUSTER_PALETTE = ["#2E6B8A", "#3A8C7C", "#5B5EA6", "#4682B4", "#2B5E4A"]

# ============================================================================
# Cluster highlight definitions
# ============================================================================

HIGHLIGHT_SETS = {
    "Interphase": {
        "mitochondrial": {
            "label": "Mitochondrial",
            "brieflow_clusters": [22, 69, 121, 129, 160],
            "funk_clusters": [135, 147, 149, 185],
        },
    },
}

# ============================================================================
# Data loading
# ============================================================================


def load_phate_data(cell_class):
    """Load PHATE coordinates and cluster assignments for both pipelines."""
    cfg = IDEAL_CONFIGS[cell_class]

    brieflow_path = (
        CLUSTER_DIR
        / cell_class
        / str(cfg["brieflow_k"])
        / "phate_leiden_clustering.tsv"
    )
    brieflow_df = pd.read_csv(brieflow_path, sep="\t")
    brieflow_df = brieflow_df.rename(columns={"gene_symbol_0": "gene_symbol"})

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
    """Load MozzareLLM confidence per cluster for both pipelines."""
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


def _load_mozzarellm_for_clusters(cell_class, clusters, pipeline="brieflow"):
    """Load MozzareLLM labels and gene categories for specified clusters.

    Returns:
        labels: {cluster_id: dominant_process}
        categories: {gene_name: category} (absent = "established")
    """
    cfg = IDEAL_CONFIGS[cell_class]

    if pipeline == "brieflow":
        base = (
            CLUSTER_DIR / cell_class / str(cfg["brieflow_k"]) / "mozzarellm"
        )
    else:
        base = (
            EXTERNAL_DIR
            / "results"
            / "cluster"
            / f"Funk_{cell_class}_k{cfg['funk_k']}"
            / "mozzarellm"
        )

    summaries = pd.read_csv(base / "claude-opus-4-6_summaries.tsv", sep="\t")
    summaries = summaries[summaries["cluster_id"].isin(clusters)]
    labels = dict(zip(summaries["cluster_id"], summaries["dominant_process"]))

    flagged = pd.read_csv(base / "claude-opus-4-6_flagged_genes.tsv", sep="\t")
    flagged = flagged[flagged["cluster_id"].isin(clusters)]
    categories = dict(zip(flagged["gene"], flagged["category"]))

    return labels, categories


def load_features_for_gene_heatmap(cell_class, highlight_set, select_by="brieflow"):
    """Load PCA + raw features for gene-level heatmaps.

    Args:
        select_by: "brieflow" selects genes in Brieflow mito clusters (117),
                   "funk" selects genes in Funk mito clusters (~73).
                   Genes are ordered by the selecting pipeline's clusters.

    Returns dict with gene_list, feature arrays, cluster info, and MozzareLLM data.
    """
    cfg = IDEAL_CONFIGS[cell_class]
    hl_bl = highlight_set["brieflow_clusters"]
    hl_fk = highlight_set["funk_clusters"]

    # --- Load cluster assignments ---
    brieflow_clust_path = (
        CLUSTER_DIR
        / cell_class
        / str(cfg["brieflow_k"])
        / "phate_leiden_clustering.tsv"
    )
    brieflow_clust = pd.read_csv(brieflow_clust_path, sep="\t")
    brieflow_clust = brieflow_clust.rename(
        columns={"gene_symbol_0": "gene_symbol"}
    )
    bl_all_clusters = dict(
        zip(brieflow_clust["gene_symbol"], brieflow_clust["cluster"])
    )

    funk_clust_path = (
        EXTERNAL_DIR
        / "results"
        / "cluster"
        / f"Funk_{cell_class}_k{cfg['funk_k']}"
        / "phate_leiden_clustering.tsv"
    )
    funk_clust = pd.read_csv(funk_clust_path, sep="\t")
    fk_all_clusters = dict(
        zip(funk_clust["gene_symbol"], funk_clust["cluster"])
    )

    # --- Select gene set based on ordering pipeline ---
    if select_by == "brieflow":
        select_df = brieflow_clust[brieflow_clust["cluster"].isin(hl_bl)].copy()
    else:
        select_df = funk_clust[funk_clust["cluster"].isin(hl_fk)].copy()

    # Ensure column name is consistent
    if "gene_symbol_0" in select_df.columns:
        select_df = select_df.rename(columns={"gene_symbol_0": "gene_symbol"})

    select_df = select_df.sort_values(
        ["cluster", "gene_symbol"]
    ).drop_duplicates("gene_symbol")
    gene_list = select_df["gene_symbol"].tolist()

    # --- Load MozzareLLM data ---
    bl_labels, bl_categories = _load_mozzarellm_for_clusters(
        cell_class, hl_bl, "brieflow"
    )
    fk_labels, fk_categories = _load_mozzarellm_for_clusters(
        cell_class, hl_fk, "funk"
    )

    # --- Load PCA features ---
    agg_path = (
        AGGREGATE_DIR
        / f"CeCl-{cell_class}_ChCo-DAPI_TUBULIN_GH2AX_PHALLOIDIN__aggregated.tsv"
    )
    brieflow_agg = pd.read_csv(agg_path, sep="\t")
    brieflow_pc_cols = [c for c in brieflow_agg.columns if c.startswith("PC_")]

    funk_pca_path = EXTERNAL_DIR / f"Funk_{cell_class}_pca_features.tsv"
    funk_pca = pd.read_csv(funk_pca_path, sep="\t")
    funk_pc_cols = [c for c in funk_pca.columns if c.startswith("PC_")]
    n_pcs = len(funk_pc_cols)
    brieflow_pc_matched = brieflow_pc_cols[:n_pcs]

    # --- Load raw features ---
    brieflow_feat_path = (
        AGGREGATE_DIR
        / f"CeCl-{cell_class}_ChCo-DAPI_TUBULIN_GH2AX_PHALLOIDIN__features_genes.tsv"
    )
    brieflow_feat = pd.read_csv(brieflow_feat_path, sep="\t")
    brieflow_raw_cols = [
        c
        for c in brieflow_feat.columns
        if c not in BRIEFLOW_METADATA_COLS
        and pd.api.types.is_numeric_dtype(brieflow_feat[c])
    ]

    funk_feat_path = EXTERNAL_DIR / f"Funk_{cell_class}_cleaned_features.tsv"
    funk_feat = pd.read_csv(funk_feat_path, sep="\t")
    funk_feat = funk_feat.rename(columns={"gene_symbol": "gene_symbol_0"})
    funk_raw_cols = [c for c in funk_feat.columns if c != "gene_symbol_0"]

    # Case-insensitive feature matching
    bl_lower = {c.lower(): c for c in brieflow_raw_cols}
    fk_lower = {c.lower(): c for c in funk_raw_cols}
    shared_keys = sorted(set(bl_lower) & set(fk_lower))
    bl_matched_raw = [bl_lower[k] for k in shared_keys]
    fk_matched_raw = [fk_lower[k] for k in shared_keys]

    # --- Subset to genes present in both pipelines ---
    shared_genes = (
        set(gene_list)
        & set(brieflow_agg["gene_symbol_0"])
        & set(funk_pca["gene_symbol"])
        & set(brieflow_feat["gene_symbol_0"])
        & set(funk_feat["gene_symbol_0"])
    )
    gene_list = [g for g in gene_list if g in shared_genes]

    def _extract(df, gene_col, cols):
        sub = df[df[gene_col].isin(shared_genes)].drop_duplicates(gene_col)
        sub = sub.set_index(gene_col).loc[gene_list]
        return sub[cols].values

    bl_pca_arr = _extract(brieflow_agg, "gene_symbol_0", brieflow_pc_matched)
    fk_pca_arr = _extract(funk_pca, "gene_symbol", funk_pc_cols)
    bl_raw_arr = _extract(brieflow_feat, "gene_symbol_0", bl_matched_raw)
    fk_raw_arr = _extract(funk_feat, "gene_symbol_0", fk_matched_raw)

    # Cluster maps for all genes in the set
    bl_gene_clusters = {g: bl_all_clusters.get(g) for g in gene_list}
    fk_gene_clusters = {g: fk_all_clusters.get(g) for g in gene_list}

    return {
        "gene_list": gene_list,
        "select_by": select_by,
        "bl_gene_clusters": bl_gene_clusters,
        "fk_gene_clusters": fk_gene_clusters,
        "bl_labels": bl_labels,
        "fk_labels": fk_labels,
        "bl_categories": bl_categories,
        "fk_categories": fk_categories,
        "bl_pca": bl_pca_arr,
        "fk_pca": fk_pca_arr,
        "bl_raw": bl_raw_arr,
        "fk_raw": fk_raw_arr,
        "hl_bl_clusters": hl_bl,
        "hl_fk_clusters": hl_fk,
    }


# ============================================================================
# PHATE scatter plots
# ============================================================================


def _spatial_cluster_colors(df, cmap_name="Spectral"):
    """Assign colors based on each cluster's mean PHATE position (angle from center).

    Clusters that are spatially near each other in the PHATE embedding get
    similar hues, producing a smooth spatial gradient instead of random noise.
    """
    cmap = plt.colormaps[cmap_name]
    cluster_ids = df["cluster"].values
    unique_clusters = sorted(df["cluster"].unique())

    # Compute mean PHATE position per cluster
    cluster_centers = df.groupby("cluster")[["PHATE_0", "PHATE_1"]].mean()

    # Use angular position from global center for color assignment
    cx = cluster_centers["PHATE_0"].median()
    cy = cluster_centers["PHATE_1"].median()
    angles = np.arctan2(
        cluster_centers["PHATE_1"] - cy,
        cluster_centers["PHATE_0"] - cx,
    )
    # Normalize to [0, 1]
    angle_norm = (angles - angles.min()) / (angles.max() - angles.min() + 1e-10)
    cluster_color_map = {cl: cmap(angle_norm[cl]) for cl in unique_clusters}

    return [cluster_color_map[c] for c in cluster_ids]


def _highlight_cluster_colors(n_clusters, pipeline="Brieflow"):
    """Pipeline-specific colors for highlighted clusters.

    Brieflow: blue tones, Funk: orange tones — matching the global color scheme.
    """
    palette = (BRIEFLOW_CLUSTER_PALETTE if pipeline == "Brieflow"
               else FUNK_CLUSTER_PALETTE)
    return [palette[i % len(palette)] for i in range(n_clusters)]


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
    point_size = 20

    for ax, df, pipeline, k in [
        (axes[0], brieflow_df, "Brieflow", cfg["brieflow_k"]),
        (axes[1], funk_df, "Funk", cfg["funk_k"]),
    ]:
        n_clusters = df["cluster"].nunique()

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
            hl_colors = _highlight_cluster_colors(len(hl_clusters), pipeline=pipeline)
            # Load MozzareLLM labels for this pipeline's clusters
            pip_key = "brieflow" if pipeline == "Brieflow" else "funk"
            hl_annots, _ = _load_mozzarellm_for_clusters(
                cell_class, list(hl_clusters), pip_key
            )
            for i, cl in enumerate(sorted(hl_clusters)):
                mask = df["cluster"] == cl
                if mask.sum() == 0:
                    continue
                ax.scatter(
                    df.loc[mask, "PHATE_0"],
                    df.loc[mask, "PHATE_1"],
                    c=[hl_colors[i]],
                    s=point_size,
                    alpha=0.9,
                    rasterized=True,
                    linewidths=0.3,
                    edgecolors="black",
                    zorder=3,
                )
                annot = hl_annots.get(cl, "")
                lbl = f"Cluster {cl}: {annot}" if annot else f"Cluster {cl} ({mask.sum()} genes)"
                legend_handles.append(
                    mpatches.Patch(color=hl_colors[i], label=lbl)
                )
            # Legend info stored but not drawn on the plot
        else:
            # Spatial coloring: nearby clusters get similar hues
            c_arr = _spatial_cluster_colors(df)
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
        ax.set_xlabel("PHATE 1", fontsize=14)
        ax.set_ylabel("PHATE 2", fontsize=14)
        if highlight_set is None:
            ax.set_title(
                f"{pipeline} (k={k}, {n_clusters} clusters, {n_genes} genes)",
                fontweight="bold",
            )
        _square_scatter_axes(ax)

    if highlight_set is None:
        fig.suptitle(
            f"{cell_class} PHATE Embeddings",
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
        conf_levels = df["cluster"].map(conf_map).fillna("Low").values

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
        ax.set_title(f"{pipeline} (k={k})", fontweight="bold")
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
# Split-diagonal heatmap
# ============================================================================


def _center(X):
    """Column-center so cosine_similarity becomes Pearson correlation."""
    return X - X.mean(axis=0, keepdims=True)


def _make_cluster_colormap(hl_clusters, pipeline="Brieflow"):
    """Build {cluster_id: color} for highlighted clusters + 'other' in gray.

    Uses the unified palette shared with PHATE scatter highlights.
    """
    palette = (BRIEFLOW_CLUSTER_PALETTE if pipeline == "Brieflow"
               else FUNK_CLUSTER_PALETTE)
    cmap = {}
    for i, cl in enumerate(sorted(hl_clusters)):
        cmap[cl] = palette[i % len(palette)]
    return cmap


OTHER_COLOR = "#dddddd"


def _draw_cluster_lines(ax, gene_list, gene_clusters, cluster_cmap, x_pos,
                        linewidth=3):
    """Draw colored line segments at the heatmap edge for cluster assignment.

    One horizontal line per gene row, colored by cluster. Drawn in data coords
    so they sit exactly at the heatmap boundary with no overlap.
    """
    for i, g in enumerate(gene_list):
        cl = gene_clusters.get(g)
        color = cluster_cmap.get(cl, OTHER_COLOR)
        ax.plot([x_pos, x_pos], [i, i + 1],
                color=color, linewidth=linewidth, solid_capstyle="butt",
                clip_on=False)


def _cluster_legend_handles(cluster_cmap, labels, pipeline, include_other=False):
    """Build legend handles mapping cluster colors to dominant_process."""
    handles = []
    for cl in sorted(cluster_cmap):
        desc = labels.get(cl, f"Cluster {cl}")
        handles.append(mpatches.Patch(
            color=cluster_cmap[cl],
            label=f"{pipeline} {cl}: {desc}",
        ))
    if include_other:
        handles.append(mpatches.Patch(
            color=OTHER_COLOR, label="Other cluster",
        ))
    return handles



def _save_companion_legend(bl_cmap, fk_cmap, bl_labels, fk_labels, output_path):
    """Save a single companion legend image with colorbar, category key,
    and cluster color swatches for both pipelines."""
    setup_plot_style()
    fig, axes = plt.subplots(3, 1, figsize=(7, 7),
                             gridspec_kw={"height_ratios": [0.5, 0.8, 4]})

    # --- Colorbar ---
    ax_cb = axes[0]
    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=plt.Normalize(-1, 1))
    fig.colorbar(sm, cax=ax_cb, orientation="horizontal")
    ax_cb.set_xlabel("Pearson r (per-triangle normalized)", fontsize=13)
    ax_cb.tick_params(labelsize=12)

    # --- Category key ---
    ax_cat = axes[1]
    ax_cat.axis("off")
    ax_cat.text(0.05, 0.5, "Established", transform=ax_cat.transAxes,
                fontsize=14, va="center", fontweight="bold")
    ax_cat.text(0.38, 0.5, "Novel role / Uncharacterized", transform=ax_cat.transAxes,
                fontsize=14, va="center", fontweight="bold", fontstyle="italic")
    ax_cat.set_title("MozzareLLM gene category", fontsize=14, fontweight="bold",
                     loc="left")

    # --- Cluster colors ---
    ax_cl = axes[2]
    ax_cl.axis("off")
    y = 0.95
    dy = 0.10
    ax_cl.text(0.0, y, "Brieflow clusters", fontsize=13,
               fontweight="bold", transform=ax_cl.transAxes, va="top")
    y -= dy
    for cl in sorted(bl_cmap):
        desc = bl_labels.get(cl, "")
        ax_cl.add_patch(plt.Rectangle((0.0, y - 0.03), 0.04, 0.05,
                                      color=bl_cmap[cl], transform=ax_cl.transAxes,
                                      clip_on=False))
        ax_cl.text(0.06, y - 0.005, f"Cluster {cl}: {desc}",
                   fontsize=12, transform=ax_cl.transAxes, va="center")
        y -= dy
    y -= 0.03
    ax_cl.text(0.0, y, "Funk clusters", fontsize=13,
               fontweight="bold", transform=ax_cl.transAxes, va="top")
    y -= dy
    for cl in sorted(fk_cmap):
        desc = fk_labels.get(cl, "")
        ax_cl.add_patch(plt.Rectangle((0.0, y - 0.03), 0.04, 0.05,
                                      color=fk_cmap[cl], transform=ax_cl.transAxes,
                                      clip_on=False))
        ax_cl.text(0.06, y - 0.005, f"Cluster {cl}: {desc}",
                   fontsize=12, transform=ax_cl.transAxes, va="center")
        y -= dy

    fig.tight_layout()
    save_figure(fig, output_path, transparent=True)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_mito_heatmap(data, feature_type, output_path):
    """Split-diagonal heatmap with staggered gene names on the Y-axis.

    Upper triangle = Brieflow correlations, lower = Funk correlations.
    Gene names are staggered at two horizontal offsets with connecting lines.
    Cluster color bars sit flush against the heatmap (no tick gaps).
    Companion legend and metrics TSV are saved alongside.

    Args:
        data: dict from load_features_for_gene_heatmap()
        feature_type: "pca" or "raw"
        output_path: where to save
    """
    setup_plot_style()

    gene_list = data["gene_list"]
    n = len(gene_list)
    select_by = data["select_by"]

    bl_cmap = _make_cluster_colormap(data["hl_bl_clusters"], "Brieflow")
    fk_cmap = _make_cluster_colormap(data["hl_fk_clusters"], "Funk")

    if feature_type == "pca":
        corr_bl = cosine_similarity(_center(data["bl_pca"]))
        corr_fk = cosine_similarity(_center(data["fk_pca"]))
    else:
        corr_bl = cosine_similarity(_center(data["bl_raw"]))
        corr_fk = cosine_similarity(_center(data["fk_raw"]))

    # Per-triangle 99th-percentile normalization
    upper_idx = np.triu_indices(n, k=1)
    lower_idx = np.tril_indices(n, k=-1)

    def _norm(vals):
        vmax = np.percentile(np.abs(vals), 99)
        return np.clip(vals / vmax, -1, 1) if vmax > 0 else vals

    combined = np.zeros((n, n))
    combined[upper_idx] = _norm(corr_bl[upper_idx])
    combined[lower_idx] = _norm(corr_fk[lower_idx])

    bl_gcl = data["bl_gene_clusters"]
    fk_gcl = data["fk_gene_clusters"]

    # Determine gene category for the selecting pipeline
    sel_categories = (data["bl_categories"] if select_by == "brieflow"
                      else data["fk_categories"])

    # Helper: font style from MozzareLLM category
    # established = bold; novel_role & uncharacterized = bold italic
    def _gene_style(gene):
        cat = sel_categories.get(gene, "established")
        if cat in ("novel_role", "uncharacterized"):
            return gene, {"fontweight": "bold", "fontstyle": "italic"}
        return gene, {"fontweight": "bold"}

    # Figure sizing — gene names should dominate visually
    # Funk-ordered heatmaps (fewer genes) get slightly smaller font
    gene_fontsize = 20 if n <= 90 else 20
    fig_side = max(20, n * 0.28)
    fig_width = fig_side + 8

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(fig_width, fig_side))
    fig.subplots_adjust(left=0.15, right=0.85, top=0.95, bottom=0.03)

    ax.pcolormesh(combined, cmap="RdBu_r", vmin=-1, vmax=1, rasterized=True)
    ax.set_aspect("equal")
    ax.set_xlim(0, n)
    ax.set_ylim(n, 0)  # inverted
    ax.plot([0, n], [0, n], color="black", linewidth=0.5, alpha=0.4)

    # Thin black outlines around each cluster block on the diagonal
    sel_gcl = bl_gcl if select_by == "brieflow" else fk_gcl
    prev_cl = sel_gcl.get(gene_list[0])
    block_start = 0
    for i in range(1, n + 1):
        cur_cl = sel_gcl.get(gene_list[i]) if i < n else None
        if cur_cl != prev_cl:
            ax.add_patch(plt.Rectangle(
                (block_start, block_start), i - block_start, i - block_start,
                fill=False, edgecolor="black", linewidth=0.5, clip_on=True,
            ))
            block_start = i
            prev_cl = cur_cl

    # Remove all tick marks — color bars sit flush
    ax.tick_params(axis="both", length=0, pad=0)
    ax.set_xticks([])
    ax.set_yticks([])

    # Cluster color lines tight against heatmap edges (data coords)
    cl_gap = n * 0.005  # tiny gap between heatmap and color line
    _draw_cluster_lines(ax, gene_list, fk_gcl, fk_cmap, x_pos=-cl_gap, linewidth=6)
    _draw_cluster_lines(ax, gene_list, bl_gcl, bl_cmap, x_pos=n + cl_gap, linewidth=6)

    # Staggered gene names — all in data coordinates for consistent spacing
    # Offsets are proportional to n so everything scales together
    label_gap = n * 0.02    # gap between cluster line and nearest gene label
    label_close = cl_gap + label_gap           # level 1 (even-indexed genes)
    label_far = cl_gap + label_gap + n * 0.05  # level 2 (odd-indexed genes)
    line_alpha = 0.3

    for i, gene in enumerate(gene_list):
        y_pos = i + 0.5
        label_text, font_kw = _gene_style(gene)
        is_odd = i % 2 == 1

        # Left side — Funk cluster colors
        fk_color = fk_cmap.get(fk_gcl.get(gene), OTHER_COLOR)
        x_left = -(label_far if is_odd else label_close)
        ax.text(
            x_left, y_pos, label_text,
            fontsize=gene_fontsize, color=fk_color,
            ha="right", va="center", clip_on=False,
            **font_kw,
        )
        # Connecting line from label to cluster line (all data coords)
        ax.annotate(
            "", xy=(-cl_gap, y_pos), xytext=(x_left + n * 0.002, y_pos),
            xycoords="data", textcoords="data",
            arrowprops=dict(arrowstyle="-", color=fk_color,
                            lw=0.4, alpha=line_alpha),
        )

        # Right side — Brieflow cluster colors
        bl_color = bl_cmap.get(bl_gcl.get(gene), OTHER_COLOR)
        x_right = n + (label_far if is_odd else label_close)
        ax.text(
            x_right, y_pos, label_text,
            fontsize=gene_fontsize, color=bl_color,
            ha="left", va="center", clip_on=False,
            **font_kw,
        )
        ax.annotate(
            "", xy=(n + cl_gap, y_pos), xytext=(x_right - n * 0.002, y_pos),
            xycoords="data", textcoords="data",
            arrowprops=dict(arrowstyle="-", color=bl_color,
                            lw=0.4, alpha=line_alpha),
        )

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
            plot_phate_scatter(
                cell_class,
                OUTPUT_DIR / f"phate_scatter_{cell_class.lower()}.png",
            )
            plot_phate_confidence(
                cell_class,
                OUTPUT_DIR / f"phate_confidence_{cell_class.lower()}.png",
            )
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
                    cc = cell_class.lower()
                    # Generate single companion legend per highlight set
                    bl_cmap = _make_cluster_colormap(hl_set["brieflow_clusters"], "Brieflow")
                    fk_cmap = _make_cluster_colormap(hl_set["funk_clusters"], "Funk")
                    bl_labels, _ = _load_mozzarellm_for_clusters(
                        cell_class, hl_set["brieflow_clusters"], "brieflow")
                    fk_labels, _ = _load_mozzarellm_for_clusters(
                        cell_class, hl_set["funk_clusters"], "funk")
                    _save_companion_legend(
                        bl_cmap, fk_cmap, bl_labels, fk_labels,
                        OUTPUT_DIR / f"heatmap_{cc}_{key}_legend.pdf",
                    )
                    for select_by in ["brieflow", "funk"]:
                        data = load_features_for_gene_heatmap(
                            cell_class, hl_set, select_by=select_by
                        )
                        for feat in ["pca", "raw"]:
                            plot_mito_heatmap(
                                data, feat,
                                OUTPUT_DIR / f"heatmap_{cc}_{key}_{feat}_{select_by}.png",
                            )

    print(f"\nOutputs saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
