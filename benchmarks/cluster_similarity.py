"""
Cluster Similarity: Brieflow vs Funk et al. 2022

Compares the feature embeddings from each pipeline via pairwise Pearson
correlation heatmaps.  For each cell class (Interphase, Mitotic) a 2x2 grid
is produced:

    Top row    — PCA features (Brieflow truncated to match Funk PC count)
    Bottom row — Raw / cleaned CellProfiler features

    Left col   — genes ordered by Brieflow clusters
    Right col  — genes ordered by Funk clusters

Each panel is a split-diagonal heatmap: upper triangle = Brieflow Pearson
correlation, lower triangle = Funk Pearson correlation.

Note: Pearson correlation is equivalent to cosine similarity on column-centered
data.  Raw cosine similarity on non-centered data is dominated by the shared
mean component, producing artifactually high values.

Usage:
    python cluster_similarity.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.metrics.pairwise import cosine_similarity

from plot_style import setup_plot_style, save_figure


def _center(X):
    """Column-center a matrix so cosine_similarity becomes Pearson correlation."""
    return X - X.mean(axis=0, keepdims=True)


# Paths
BENCHMARKS_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BENCHMARKS_DIR.parent
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
AGGREGATE_DIR = ANALYSIS_DIR / "brieflow_output" / "aggregate" / "tsvs"
CLUSTER_DIR = (
    ANALYSIS_DIR / "brieflow_output" / "cluster" / "DAPI_TUBULIN_GH2AX_PHALLOIDIN"
)
EXTERNAL_DIR = BENCHMARKS_DIR / "external"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "cluster" / "similarity"

BRIEFLOW_METADATA_COLS = {
    "gene_symbol_0",
    "cell_count",
    "cell_stage",
    "cell_stage_confidence",
    "perturbation_auc",
}

CELL_CONFIGS = [
    {
        "name": "Interphase",
        "brieflow_agg": AGGREGATE_DIR
        / "CeCl-Interphase_ChCo-DAPI_TUBULIN_GH2AX_PHALLOIDIN__aggregated.tsv",
        "brieflow_features": AGGREGATE_DIR
        / "CeCl-Interphase_ChCo-DAPI_TUBULIN_GH2AX_PHALLOIDIN__features_genes.tsv",
        "brieflow_cluster": CLUSTER_DIR
        / "Interphase"
        / "12"
        / "phate_leiden_clustering.tsv",
        "funk_pca": EXTERNAL_DIR / "Funk_Interphase_pca_features.tsv",
        "funk_features": EXTERNAL_DIR / "Funk_Interphase_cleaned_features.tsv",
        "funk_cluster": EXTERNAL_DIR
        / "results"
        / "cluster"
        / "Funk_Interphase_k10"
        / "phate_leiden_clustering.tsv",
        "brieflow_k": 12,
        "funk_k": 10,
    },
    {
        "name": "Mitotic",
        "brieflow_agg": AGGREGATE_DIR
        / "CeCl-Mitotic_ChCo-DAPI_TUBULIN_GH2AX_PHALLOIDIN__aggregated.tsv",
        "brieflow_features": AGGREGATE_DIR
        / "CeCl-Mitotic_ChCo-DAPI_TUBULIN_GH2AX_PHALLOIDIN__features_genes.tsv",
        "brieflow_cluster": CLUSTER_DIR
        / "Mitotic"
        / "5"
        / "phate_leiden_clustering.tsv",
        "funk_pca": EXTERNAL_DIR / "Funk_Mitotic_pca_features.tsv",
        "funk_features": EXTERNAL_DIR / "Funk_Mitotic_cleaned_features.tsv",
        "funk_cluster": EXTERNAL_DIR
        / "results"
        / "cluster"
        / "Funk_Mitotic_k9"
        / "phate_leiden_clustering.tsv",
        "brieflow_k": 5,
        "funk_k": 9,
    },
]


def load_cell_data(cfg):
    """Load all data for one cell class. Returns dict with arrays and metadata."""
    name = cfg["name"]
    print(f"\n{'—' * 60}")
    print(f"Loading {name} data...")

    # --- Brieflow PCA (from aggregate) ---
    brieflow_agg = pd.read_csv(cfg["brieflow_agg"], sep="\t")
    brieflow_pc_cols = [c for c in brieflow_agg.columns if c.startswith("PC_")]
    print(f"  Brieflow PCA: {len(brieflow_agg)} genes, {len(brieflow_pc_cols)} PCs")

    # --- Brieflow raw features ---
    brieflow_feat_df = pd.read_csv(cfg["brieflow_features"], sep="\t")
    brieflow_raw_cols = [
        c for c in brieflow_feat_df.columns if c not in BRIEFLOW_METADATA_COLS
    ]
    # Keep only numeric columns
    brieflow_raw_cols = [
        c
        for c in brieflow_raw_cols
        if pd.api.types.is_numeric_dtype(brieflow_feat_df[c])
    ]
    print(
        f"  Brieflow features: {len(brieflow_feat_df)} genes, {len(brieflow_raw_cols)} features"
    )

    # --- Brieflow cluster assignments ---
    brieflow_clust = pd.read_csv(cfg["brieflow_cluster"], sep="\t")
    brieflow_clust = brieflow_clust[["gene_symbol_0", "cluster"]].rename(
        columns={"cluster": "brieflow_cluster"}
    )

    # --- Funk PCA ---
    funk_pca_df = pd.read_csv(cfg["funk_pca"], sep="\t")
    funk_pc_cols = [c for c in funk_pca_df.columns if c.startswith("PC_")]
    funk_pca_df = funk_pca_df.rename(columns={"gene_symbol": "gene_symbol_0"})
    n_funk_pcs = len(funk_pc_cols)
    print(f"  Funk PCA: {len(funk_pca_df)} genes, {n_funk_pcs} PCs")

    # --- Funk cleaned features ---
    funk_feat_df = pd.read_csv(cfg["funk_features"], sep="\t")
    funk_feat_df = funk_feat_df.rename(columns={"gene_symbol": "gene_symbol_0"})
    funk_raw_cols = [c for c in funk_feat_df.columns if c != "gene_symbol_0"]
    print(f"  Funk features: {len(funk_feat_df)} genes, {len(funk_raw_cols)} features")

    # --- Match Brieflow features to Funk's feature set (case-insensitive) ---
    # All Funk features are a subset of Brieflow's (just different casing)
    bl_lower_to_orig = {c.lower(): c for c in brieflow_raw_cols}
    fk_lower_to_orig = {c.lower(): c for c in funk_raw_cols}
    shared_feat_keys = sorted(set(bl_lower_to_orig) & set(fk_lower_to_orig))
    brieflow_matched_cols = [bl_lower_to_orig[k] for k in shared_feat_keys]
    funk_matched_cols = [fk_lower_to_orig[k] for k in shared_feat_keys]
    n_shared_feats = len(shared_feat_keys)
    print(f"  Overlapping features: {n_shared_feats}")

    # --- Funk cluster assignments ---
    funk_clust = pd.read_csv(cfg["funk_cluster"], sep="\t")
    funk_clust = funk_clust[["gene_symbol", "cluster"]].rename(
        columns={"gene_symbol": "gene_symbol_0", "cluster": "funk_cluster"}
    )

    # --- Shared genes across all six datasets ---
    shared_genes = (
        set(brieflow_agg["gene_symbol_0"])
        & set(brieflow_feat_df["gene_symbol_0"])
        & set(brieflow_clust["gene_symbol_0"])
        & set(funk_pca_df["gene_symbol_0"])
        & set(funk_feat_df["gene_symbol_0"])
        & set(funk_clust["gene_symbol_0"])
    )
    print(f"  Shared genes: {len(shared_genes)}")

    # --- Build merged dataframe (gene_symbol_0 + clusters only) ---
    gene_list = sorted(shared_genes)
    gene_order = {g: i for i, g in enumerate(gene_list)}

    # Collect everything keyed on shared genes
    def subset_and_sort(df, cols=None):
        sub = df[df["gene_symbol_0"].isin(shared_genes)].drop_duplicates(
            "gene_symbol_0"
        )
        sub = sub.assign(_order=sub["gene_symbol_0"].map(gene_order))
        sub = sub.sort_values("_order").drop(columns="_order").reset_index(drop=True)
        if cols is not None:
            return sub[cols].values
        return sub

    # Cluster assignments
    clust_df = subset_and_sort(brieflow_clust)
    clust_df = clust_df.merge(
        funk_clust[funk_clust["gene_symbol_0"].isin(shared_genes)].drop_duplicates(
            "gene_symbol_0"
        ),
        on="gene_symbol_0",
        how="inner",
    )
    clust_df = clust_df.assign(_order=clust_df["gene_symbol_0"].map(gene_order))
    clust_df = (
        clust_df.sort_values("_order").drop(columns="_order").reset_index(drop=True)
    )

    # Truncate Brieflow PCs to match Funk PC count
    brieflow_pc_matched = brieflow_pc_cols[:n_funk_pcs]

    # Extract arrays
    brieflow_pca_arr = subset_and_sort(brieflow_agg, brieflow_pc_matched)
    funk_pca_arr = subset_and_sort(funk_pca_df, funk_pc_cols)
    brieflow_feat_arr = subset_and_sort(brieflow_feat_df, brieflow_matched_cols)
    funk_feat_arr = subset_and_sort(funk_feat_df, funk_matched_cols)

    n = len(clust_df)
    print(
        f"  Final: {n} genes | "
        f"PCA: {brieflow_pca_arr.shape[1]} matched PCs | "
        f"Features: {n_shared_feats} shared features"
    )

    return {
        "name": name,
        "n": n,
        "brieflow_k": cfg["brieflow_k"],
        "funk_k": cfg["funk_k"],
        "brieflow_cluster": clust_df["brieflow_cluster"].values,
        "funk_cluster": clust_df["funk_cluster"].values,
        "brieflow_pca": brieflow_pca_arr,
        "funk_pca": funk_pca_arr,
        "n_pcs": n_funk_pcs,
        "brieflow_features": brieflow_feat_arr,
        "funk_features": funk_feat_arr,
        "n_shared_feats": n_shared_feats,
    }


def _build_split_diagonal(cos_upper, cos_lower, order):
    """Build a combined matrix with each triangle independently normalized to [-1, 1].

    This ensures fair visual comparison between triangles that may have different
    cosine similarity distributions (e.g. due to different feature dimensionalities).
    """
    cos_u = cos_upper[np.ix_(order, order)]
    cos_l = cos_lower[np.ix_(order, order)]
    n = len(order)

    upper_idx = np.triu_indices(n, k=1)
    lower_idx = np.tril_indices(n, k=-1)
    diag_idx = np.diag_indices(n)

    # Normalize each triangle to [-1, 1] using symmetric percentile clipping
    def normalize(vals):
        vmax = np.percentile(np.abs(vals), 99)
        if vmax == 0:
            return vals
        return np.clip(vals / vmax, -1, 1)

    combined = np.zeros((n, n))
    combined[upper_idx] = normalize(cos_u[upper_idx])
    combined[lower_idx] = normalize(cos_l[lower_idx])
    combined[diag_idx] = 0  # neutral diagonal
    return combined


def plot_2x2_grid(data, output_path):
    """Generate 2×2 split-diagonal heatmap grid for one cell class."""
    setup_plot_style()
    name = data["name"]
    n = data["n"]
    bk = data["brieflow_k"]
    fk = data["funk_k"]

    print(f"\n  Computing Pearson correlations for {name}...")
    cos_pca_brieflow = cosine_similarity(_center(data["brieflow_pca"]))
    cos_pca_funk = cosine_similarity(_center(data["funk_pca"]))
    cos_feat_brieflow = cosine_similarity(_center(data["brieflow_features"]))
    cos_feat_funk = cosine_similarity(_center(data["funk_features"]))

    brieflow_order = data["brieflow_cluster"].argsort(kind="stable")
    funk_order = data["funk_cluster"].argsort(kind="stable")

    fig, axes = plt.subplots(2, 2, figsize=(16, 16))

    panels = [
        # (row, col, cos_upper, cos_lower, order, title)
        (
            0,
            0,
            cos_pca_brieflow,
            cos_pca_funk,
            brieflow_order,
            f"PCA ({data['n_pcs']} PCs) — Brieflow k={bk} ordering",
        ),
        (
            0,
            1,
            cos_pca_brieflow,
            cos_pca_funk,
            funk_order,
            f"PCA ({data['n_pcs']} PCs) — Funk k={fk} ordering",
        ),
        (
            1,
            0,
            cos_feat_brieflow,
            cos_feat_funk,
            brieflow_order,
            f"Shared Features ({data['n_shared_feats']}) — Brieflow k={bk} ordering",
        ),
        (
            1,
            1,
            cos_feat_brieflow,
            cos_feat_funk,
            funk_order,
            f"Shared Features ({data['n_shared_feats']}) — Funk k={fk} ordering",
        ),
    ]

    for row, col, cos_u, cos_l, order, title in panels:
        ax = axes[row, col]
        combined = _build_split_diagonal(cos_u, cos_l, order)
        ax.pcolormesh(combined, cmap="RdBu_r", vmin=-1, vmax=1, rasterized=True)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.invert_yaxis()
        ax.plot([0, n], [0, n], color="black", linewidth=0.5, alpha=0.4)

        # Triangle labels
        if row == 0:
            bl_label = f"Brieflow\n({data['n_pcs']} PCs)"
            fk_label = f"Funk\n({data['n_pcs']} PCs)"
        else:
            bl_label = f"Brieflow\n({data['n_shared_feats']} feats)"
            fk_label = f"Funk\n({data['n_shared_feats']} feats)"
        ax.text(
            n * 0.75,
            n * 0.25,
            bl_label,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            color="black",
            alpha=0.6,
        )
        ax.text(
            n * 0.25,
            n * 0.75,
            fk_label,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            color="black",
            alpha=0.6,
        )
        ax.set_title(title, fontsize=11, fontweight="bold")

    fig.suptitle(
        f"{name} — Pairwise Gene Pearson Correlation (each triangle independently normalized)",
        fontsize=14,
        fontweight="bold",
        y=0.93,
    )
    fig.subplots_adjust(
        left=0.02, right=0.98, top=0.90, bottom=0.02, wspace=0.05, hspace=0.10
    )
    save_figure(fig, output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("FEATURE EMBEDDING PEARSON CORRELATION — Brieflow vs Funk")
    print("=" * 70)

    for cfg in CELL_CONFIGS:
        data = load_cell_data(cfg)
        plot_2x2_grid(
            data, OUTPUT_DIR / f"pearson_correlation_{cfg['name'].lower()}.png"
        )

    print("\n" + "=" * 70)
    print("ALL PLOTS COMPLETE")
    print("=" * 70)
    print(f"Outputs in: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
