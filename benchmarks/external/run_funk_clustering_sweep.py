"""
Run Funk et al. 2022 clustering across multiple resolutions (k=5-15).

This script runs the exact Funk clustering pipeline at different resolution
parameters and saves the results for downstream enrichment analysis.

It also extracts and saves intermediate feature representations (cleaned
features and PCA) for downstream similarity comparisons with Brieflow.

Usage:
    conda activate ops_clustering

    # Full clustering sweep (slow — PHATE + Leiden at each k):
    python run_funk_clustering_sweep.py

    # Extract and save features only (no clustering):
    python run_funk_clustering_sweep.py --features-only
"""

import argparse
from igraph import Graph
import leidenalg
import numpy as np
import pandas as pd
import phate
from sklearn.preprocessing import StandardScaler
from sklearn import decomposition
from pathlib import Path

# Paths
EXTERNAL_DIR = Path(__file__).parent
INTERPHASE_HDF = EXTERNAL_DIR / "interphase-reclassified_cp_phenotype_gene_medians.20210429.hdf"
MITOTIC_HDF = EXTERNAL_DIR / "mitotic-reclassified_cp_phenotype_gene_medians.20210419.hdf"
FEATURES_I_TXT = EXTERNAL_DIR / "features_i.txt"
FEATURES_M_TXT = EXTERNAL_DIR / "features_m.txt"
CLUSTER_OUTPUT_DIR = EXTERNAL_DIR / "results" / "cluster"

# Load feature lists
with open(FEATURES_I_TXT, 'r') as f:
    features_i = [line.strip() for line in f if line.strip()]
with open(FEATURES_M_TXT, 'r') as f:
    features_m = [line.strip() for line in f if line.strip()]

VARIANCE_EXPLAINED_THRESHOLD = 0.95


def parse_args():
    parser = argparse.ArgumentParser(description="Funk et al. 2022 clustering pipeline")
    parser.add_argument(
        "--features-only",
        action="store_true",
        help="Only extract and save cleaned features + PCA (no clustering)",
    )
    return parser.parse_args()


# PHATE and clustering functions (exact same as Funk et al.)
def phate_(df, return_p=False, **kwargs):
    p = phate.PHATE(random_state=42, **kwargs)
    X_phate = p.fit_transform(df.values)
    df_phate = pd.DataFrame(X_phate, index=df.index, columns=['PHATE_0', 'PHATE_1'])
    return df_phate, p


def phate_leiden(df, knn=10, metric='euclidean', resolution_parameter=1, **kwargs):
    df_phate, p = phate_(df, knn=knn, knn_dist=metric, n_jobs=4, **kwargs)

    weights = np.asarray(p.graph.diff_op.todense())

    g = Graph().Weighted_Adjacency(matrix=weights.tolist(), mode='undirected')

    leiden_clusters = leidenalg.find_partition(
        g,
        partition_type=leidenalg.RBConfigurationVertexPartition,
        weights=g.es['weight'],
        n_iterations=-1,
        seed=42,
        resolution_parameter=resolution_parameter
    )

    df_phate['cluster'] = leiden_clusters.membership

    return df_phate


def extract_features(df_genes, features, cell_class):
    """Extract cleaned features and PCA from Funk pipeline.

    Returns:
        gene_symbols: list of gene symbols
        df_cleaned: DataFrame of StandardScaler-normalized features
        df_pca: DataFrame of PCA features (95% variance)
        n_pca: number of PCA components retained
    """
    print(f"\nExtracting features for {cell_class}...")

    gene_symbols = [idx[0] for idx in df_genes.index]

    # Normalize features (fit on controls, transform all)
    scaler = StandardScaler()
    scaler.fit(df_genes.query('gene_id=="-1"')[features].values)
    cleaned = scaler.transform(df_genes[features].values)
    df_cleaned = pd.DataFrame(cleaned, columns=features, index=df_genes.index)
    print(f"  Cleaned features: {df_cleaned.shape}")

    # PCA (to 95% variance explained)
    pca = decomposition.PCA()
    pca_values = pca.fit_transform(df_cleaned.values)
    n_pca = np.argwhere(
        pca.explained_variance_ratio_.cumsum() > VARIANCE_EXPLAINED_THRESHOLD
    ).T[0][0] + 1
    pca_cols = [f"PC_{n}" for n in range(n_pca)]
    df_pca = pd.DataFrame(pca_values[:, :n_pca], columns=pca_cols, index=df_genes.index)
    print(f"  PCA features: {df_pca.shape} ({n_pca} components, 95% variance)")

    return gene_symbols, df_cleaned, df_pca, n_pca


def save_features(gene_symbols, df_cleaned, df_pca, cell_class):
    """Save cleaned features and PCA to TSV files in external dir."""
    prefix = f"Funk_{cell_class}"

    # Cleaned features
    out_cleaned = EXTERNAL_DIR / f"{prefix}_cleaned_features.tsv"
    df_out = df_cleaned.copy()
    df_out.insert(0, "gene_symbol", gene_symbols)
    df_out.to_csv(out_cleaned, sep="\t", index=False)
    print(f"  Saved: {out_cleaned}")

    # PCA features
    out_pca = EXTERNAL_DIR / f"{prefix}_pca_features.tsv"
    df_out = df_pca.copy()
    df_out.insert(0, "gene_symbol", gene_symbols)
    df_out.to_csv(out_pca, sep="\t", index=False)
    print(f"  Saved: {out_pca}")


def run_clustering_for_resolution(df_pca, df_genes, cell_class, resolution_k):
    """Run PHATE + Leiden clustering on pre-computed PCA features."""
    print(f"\nClustering {cell_class} at k={resolution_k}...")

    df_phate = phate_leiden(
        df_pca,
        knn=5,
        metric='euclidean',
        n_pca=None,
        resolution_parameter=resolution_k
    )
    print(f"  Generated {df_phate['cluster'].nunique()} clusters")

    df_result = pd.DataFrame({
        'gene_symbol': [idx[0] for idx in df_genes.index],
        'gene_id': [idx[1] for idx in df_genes.index],
        'PHATE_0': df_phate['PHATE_0'].values,
        'PHATE_1': df_phate['PHATE_1'].values,
        'cluster': df_phate['cluster'].values
    })

    return df_result


def main():
    args = parse_args()

    print("=" * 80)
    print("FUNK ET AL. 2022 — FEATURE EXTRACTION & CLUSTERING")
    print("=" * 80)

    # Load data
    print("\nLoading HDF data...")
    df_genes_i = pd.read_hdf(INTERPHASE_HDF)
    df_genes_m = pd.read_hdf(MITOTIC_HDF)
    print(f"  Interphase: {df_genes_i.shape}")
    print(f"  Mitotic: {df_genes_m.shape}")

    # Extract and save features for both cell classes
    for df_genes, features, cell_class in [
        (df_genes_i, features_i, "Interphase"),
        (df_genes_m, features_m, "Mitotic"),
    ]:
        gene_symbols, df_cleaned, df_pca, n_pca = extract_features(
            df_genes, features, cell_class
        )
        save_features(gene_symbols, df_cleaned, df_pca, cell_class)

        # Run clustering sweep unless --features-only
        if not args.features_only:
            CLUSTER_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
            for k in range(5, 16):
                df_result = run_clustering_for_resolution(
                    df_pca, df_genes, cell_class, k
                )
                output_dir = CLUSTER_OUTPUT_DIR / f"Funk_{cell_class}_k{k}"
                output_dir.mkdir(parents=True, exist_ok=True)
                output_file = output_dir / "phate_leiden_clustering.tsv"
                df_result.to_csv(output_file, sep="\t", index=False)
                print(f"  Saved: {output_file}")

    print("\n" + "=" * 80)
    print("COMPLETE")
    print("=" * 80)
    if args.features_only:
        print("Feature extraction only (--features-only)")
    else:
        print("Full clustering sweep + feature extraction")
    print(f"Features saved to: {EXTERNAL_DIR}")


if __name__ == "__main__":
    main()
