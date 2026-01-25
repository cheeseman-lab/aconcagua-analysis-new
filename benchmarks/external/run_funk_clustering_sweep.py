"""
Run Funk et al. 2022 clustering across multiple resolutions (k=5-15).

This script runs the exact Funk clustering pipeline at different resolution
parameters and saves the results for downstream enrichment analysis.
"""

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
OUTPUT_DIR = EXTERNAL_DIR.parent / "results" / "cluster"

# Load feature lists
with open(FEATURES_I_TXT, 'r') as f:
    features_i = [line.strip() for line in f if line.strip()]
with open(FEATURES_M_TXT, 'r') as f:
    features_m = [line.strip() for line in f if line.strip()]

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


def run_clustering_for_resolution(df_genes, features, cell_class, resolution_k):
    """
    Run Funk clustering pipeline at a specific resolution.

    Args:
        df_genes: Gene profiles DataFrame
        features: List of features to use
        cell_class: 'Interphase' or 'Mitotic'
        resolution_k: Resolution parameter for Leiden clustering

    Returns:
        df_result: DataFrame with gene_symbol, PHATE_0, PHATE_1, cluster
    """
    print(f"\nProcessing {cell_class} at k={resolution_k}...")

    # Normalize features (fit on controls)
    scaler = StandardScaler()
    scaler.fit(df_genes.query('gene_id=="-1"')[features].values)
    df_genes_norm = df_genes[features].copy()
    df_genes_norm[features] = scaler.transform(df_genes[features].values)

    # PCA (to 95% variance explained)
    variance_explained_threshold = 0.95
    pca = decomposition.PCA()
    df_genes_pca = pd.DataFrame(
        pca.fit_transform(df_genes_norm.values),
        columns=[f'pca_{n}' for n in range(len(features))],
        index=df_genes_norm.index
    )
    n_pca = np.argwhere(pca.explained_variance_ratio_.cumsum() > variance_explained_threshold).T[0][0] + 1
    print(f"  Using {n_pca} PCA components")

    # PHATE + Leiden clustering
    df_phate = phate_leiden(
        df_genes_pca.iloc[:, :n_pca],
        knn=5,
        metric='euclidean',
        n_pca=None,
        resolution_parameter=resolution_k
    )
    print(f"  Generated {df_phate['cluster'].nunique()} clusters")

    # Create output dataframe with gene symbols
    df_result = pd.DataFrame({
        'gene_symbol': [idx[0] for idx in df_genes.index],
        'gene_id': [idx[1] for idx in df_genes.index],
        'PHATE_0': df_phate['PHATE_0'].values,
        'PHATE_1': df_phate['PHATE_1'].values,
        'cluster': df_phate['cluster'].values
    })

    return df_result


def main():
    print("="*80)
    print("FUNK ET AL. 2022 CLUSTERING - RESOLUTION SWEEP")
    print("="*80)
    print("\nRunning clustering at resolutions k=5 to k=15")

    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    print("\nLoading HDF data...")
    df_genes_i = pd.read_hdf(INTERPHASE_HDF)
    df_genes_m = pd.read_hdf(MITOTIC_HDF)
    print(f"  Interphase: {df_genes_i.shape}")
    print(f"  Mitotic: {df_genes_m.shape}")

    # Run clustering sweep for Interphase
    print("\n" + "="*80)
    print("INTERPHASE CLUSTERING SWEEP")
    print("="*80)

    for k in range(5, 16):
        df_result = run_clustering_for_resolution(df_genes_i, features_i, 'Interphase', k)

        # Save results
        output_dir = OUTPUT_DIR / f"Funk_Interphase_k{k}"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / "phate_leiden_clustering.tsv"
        df_result.to_csv(output_file, sep='\t', index=False)
        print(f"  Saved to {output_file}")

    # Run clustering sweep for Mitotic
    print("\n" + "="*80)
    print("MITOTIC CLUSTERING SWEEP")
    print("="*80)

    for k in range(5, 16):
        df_result = run_clustering_for_resolution(df_genes_m, features_m, 'Mitotic', k)

        # Save results
        output_dir = OUTPUT_DIR / f"Funk_Mitotic_k{k}"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / "phate_leiden_clustering.tsv"
        df_result.to_csv(output_file, sep='\t', index=False)
        print(f"  Saved to {output_file}")

    print("\n" + "="*80)
    print("CLUSTERING SWEEP COMPLETE")
    print("="*80)
    print(f"\nGenerated clustering results for:")
    print(f"  Interphase: k=5 to k=15 (11 resolutions)")
    print(f"  Mitotic: k=5 to k=15 (11 resolutions)")
    print(f"\nResults saved to: {OUTPUT_DIR}")
    print("\nNext step: Run enrichment analysis on all clustering results")


if __name__ == "__main__":
    main()
