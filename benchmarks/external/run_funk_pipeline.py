#!/usr/bin/env python3
"""
Consolidated Funk et al. 2022 pipeline: feature extraction, clustering, and enrichment.

Three subcommands:

    features  — Extract cleaned (StandardScaler) features + PCA (95% variance).
                Requires: conda activate ops_clustering

    cluster   — Full PHATE + Leiden sweep at k=5-15 (or specific k values).
                Requires: conda activate ops_clustering

    enrich    — Run Brieflow enrichment benchmarks in parallel.
                Requires: conda activate brieflow_aconcagua

Usage:
    # Extract features only (ops_clustering env):
    python run_funk_pipeline.py features

    # Full clustering sweep (ops_clustering env):
    python run_funk_pipeline.py cluster

    # Cluster at specific k values only:
    python run_funk_pipeline.py cluster --k-values 10 12

    # Run enrichment on all results (brieflow_aconcagua env):
    python run_funk_pipeline.py enrich --cell-class Interphase --k-values 10

    # Run enrichment on all k values for both classes:
    python run_funk_pipeline.py enrich --max-workers 22
"""

import argparse
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn import decomposition

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
EXTERNAL_DIR = Path(__file__).parent
INTERPHASE_HDF = (
    EXTERNAL_DIR / "interphase-reclassified_cp_phenotype_gene_medians.20210429.hdf"
)
MITOTIC_HDF = (
    EXTERNAL_DIR / "mitotic-reclassified_cp_phenotype_gene_medians.20210419.hdf"
)
FEATURES_I_TXT = EXTERNAL_DIR / "features_i.txt"
FEATURES_M_TXT = EXTERNAL_DIR / "features_m.txt"
CLUSTER_OUTPUT_DIR = EXTERNAL_DIR / "results" / "cluster"

VARIANCE_EXPLAINED_THRESHOLD = 0.95

# ---------------------------------------------------------------------------
# Feature lists (loaded lazily to avoid file-not-found when running enrich)
# ---------------------------------------------------------------------------
_features_cache = {}


def _load_feature_list(path):
    if path not in _features_cache:
        with open(path, "r") as f:
            _features_cache[path] = [line.strip() for line in f if line.strip()]
    return _features_cache[path]


# ===================================================================
# Shared helpers for features / cluster subcommands
# ===================================================================


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
    n_pca = (
        np.argwhere(
            pca.explained_variance_ratio_.cumsum() > VARIANCE_EXPLAINED_THRESHOLD
        ).T[0][0]
        + 1
    )
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


def _load_hdf_data():
    """Load HDF data for both cell classes."""
    print("\nLoading HDF data...")
    df_genes_i = pd.read_hdf(INTERPHASE_HDF)
    df_genes_m = pd.read_hdf(MITOTIC_HDF)
    print(f"  Interphase: {df_genes_i.shape}")
    print(f"  Mitotic: {df_genes_m.shape}")
    return df_genes_i, df_genes_m


# ===================================================================
# Clustering helpers (lazy-import igraph / leidenalg / phate)
# ===================================================================


def phate_(df, return_p=False, **kwargs):
    import phate as phate_mod

    p = phate_mod.PHATE(random_state=42, **kwargs)
    X_phate = p.fit_transform(df.values)
    df_phate = pd.DataFrame(X_phate, index=df.index, columns=["PHATE_0", "PHATE_1"])
    return df_phate, p


def phate_leiden(df, knn=10, metric="euclidean", resolution_parameter=1, **kwargs):
    from igraph import Graph
    import leidenalg

    df_phate, p = phate_(df, knn=knn, knn_dist=metric, n_jobs=4, **kwargs)

    weights = np.asarray(p.graph.diff_op.todense())

    g = Graph().Weighted_Adjacency(matrix=weights.tolist(), mode="undirected")

    leiden_clusters = leidenalg.find_partition(
        g,
        partition_type=leidenalg.RBConfigurationVertexPartition,
        weights=g.es["weight"],
        n_iterations=-1,
        seed=42,
        resolution_parameter=resolution_parameter,
    )

    df_phate["cluster"] = leiden_clusters.membership

    return df_phate


def run_clustering_for_resolution(df_pca, df_genes, cell_class, resolution_k):
    """Run PHATE + Leiden clustering on pre-computed PCA features."""
    print(f"\nClustering {cell_class} at k={resolution_k}...")

    df_phate = phate_leiden(
        df_pca,
        knn=5,
        metric="euclidean",
        n_pca=None,
        resolution_parameter=resolution_k,
    )
    print(f"  Generated {df_phate['cluster'].nunique()} clusters")

    df_result = pd.DataFrame(
        {
            "gene_symbol": [idx[0] for idx in df_genes.index],
            "gene_id": [idx[1] for idx in df_genes.index],
            "PHATE_0": df_phate["PHATE_0"].values,
            "PHATE_1": df_phate["PHATE_1"].values,
            "cluster": df_phate["cluster"].values,
        }
    )

    return df_result


# ===================================================================
# Enrichment helpers
# ===================================================================


def load_benchmark_data(config_dir):
    """Load STRING, CORUM, KEGG, and UniProt benchmark datasets."""
    string_benchmark = pd.read_csv(config_dir / "string_pair_benchmark.tsv", sep="\t")
    corum_benchmark = pd.read_csv(config_dir / "corum_group_benchmark.tsv", sep="\t")
    kegg_benchmark = pd.read_csv(config_dir / "kegg_group_benchmark.tsv", sep="\t")
    uniprot_data = pd.read_csv(config_dir / "uniprot_data.tsv", sep="\t")

    # Create UniProt lookup by expanding gene_names column
    uniprot_rows = []
    for _, row in uniprot_data.iterrows():
        gene_names = str(row["gene_names"]).split()
        for i, gene_name in enumerate(gene_names):
            uniprot_rows.append(
                {
                    "gene_name": gene_name,
                    "uniprot_entry": row["entry"],
                    "uniprot_function": row["function"],
                    "uniprot_link": row["link"],
                    "position": i,
                }
            )

    uniprot_lookup = pd.DataFrame(uniprot_rows)
    uniprot_lookup = (
        uniprot_lookup.sort_values("position")
        .drop_duplicates(subset="gene_name", keep="first")
        .drop(columns="position")
    )

    return string_benchmark, corum_benchmark, kegg_benchmark, uniprot_lookup


def convert_funk_to_brieflow_format(funk_df, uniprot_lookup):
    """Convert Funk clustering format to Brieflow format."""
    brieflow_df = funk_df.rename(columns={"gene_symbol": "gene_symbol_0"})
    brieflow_df["cell_count"] = 1

    brieflow_df = brieflow_df.merge(
        uniprot_lookup,
        left_on="gene_symbol_0",
        right_on="gene_name",
        how="left",
    )

    if "gene_name" in brieflow_df.columns:
        brieflow_df = brieflow_df.drop(columns="gene_name")

    column_order = [
        "gene_symbol_0",
        "cell_count",
        "PHATE_0",
        "PHATE_1",
        "cluster",
        "uniprot_entry",
        "uniprot_function",
        "uniprot_link",
    ]
    other_cols = [c for c in brieflow_df.columns if c not in column_order]
    brieflow_df = brieflow_df[column_order + other_cols]

    return brieflow_df


def run_enrichment_for_clustering(
    clustering_file,
    output_dir,
    string_benchmark,
    corum_benchmark,
    kegg_benchmark,
    uniprot_lookup,
):
    """Run enrichment analysis on a single Funk clustering result."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sys.path.insert(
        0, str(Path(__file__).parent.parent.parent / "brieflow" / "workflow")
    )
    from lib.cluster.benchmark_clusters import run_benchmark_analysis, save_json_results

    # Load Funk clustering result
    funk_df = pd.read_csv(clustering_file, sep="\t")

    # Convert to Brieflow format
    brieflow_df = convert_funk_to_brieflow_format(funk_df, uniprot_lookup)

    # Prepare cluster datasets
    cluster_datasets = {"Real": brieflow_df}

    # Run enrichment analysis
    (
        integrated_results,
        combined_tables,
        global_metrics,
        enrichment_pie_charts,
        enrichment_bar_charts,
    ) = run_benchmark_analysis(
        cluster_datasets,
        string_benchmark,
        corum_benchmark,
        kegg_benchmark,
        perturbation_col_name="gene_symbol_0",
        control_key=None,
        max_clusters=None,
    )

    # Save outputs
    output_dir.mkdir(parents=True, exist_ok=True)

    save_json_results(
        integrated_results["Real"],
        output_dir / "CB-Real__integrated_results.json",
    )

    combined_tables["Real"].to_csv(
        output_dir / "CB-Real__combined_table.tsv",
        sep="\t",
        index=False,
    )

    save_json_results(
        global_metrics["Real"],
        output_dir / "CB-Real__global_metrics.json",
    )

    enrichment_pie_charts["Real"].savefig(
        output_dir / "CB-Real__pie_chart.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(enrichment_pie_charts["Real"])

    enrichment_bar_charts["Real"].savefig(
        output_dir / "CB-Real__enrichment_bar_chart.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(enrichment_bar_charts["Real"])

    return global_metrics["Real"]


def _run_single_enrichment(args_tuple):
    """Worker function for parallel enrichment (must be top-level for pickle)."""
    cell_class, k, results_dir, config_dir = args_tuple

    run_name = f"Funk_{cell_class}_k{k}"
    clustering_file = results_dir / run_name / "phate_leiden_clustering.tsv"

    if not clustering_file.exists():
        return run_name, None, f"clustering file not found: {clustering_file}"

    try:
        string_benchmark, corum_benchmark, kegg_benchmark, uniprot_lookup = (
            load_benchmark_data(config_dir)
        )
        metrics = run_enrichment_for_clustering(
            clustering_file,
            clustering_file.parent,
            string_benchmark,
            corum_benchmark,
            kegg_benchmark,
            uniprot_lookup,
        )
        return run_name, metrics, None
    except Exception:
        import traceback

        return run_name, None, traceback.format_exc()


# ===================================================================
# Subcommand implementations
# ===================================================================


def cmd_features(args):
    """Extract cleaned features + PCA only (no clustering)."""
    print("=" * 80)
    print("FUNK ET AL. 2022 — FEATURE EXTRACTION")
    print("=" * 80)

    df_genes_i, df_genes_m = _load_hdf_data()
    features_i = _load_feature_list(FEATURES_I_TXT)
    features_m = _load_feature_list(FEATURES_M_TXT)

    for df_genes, features, cell_class in [
        (df_genes_i, features_i, "Interphase"),
        (df_genes_m, features_m, "Mitotic"),
    ]:
        gene_symbols, df_cleaned, df_pca, n_pca = extract_features(
            df_genes, features, cell_class
        )
        save_features(gene_symbols, df_cleaned, df_pca, cell_class)

    print("\n" + "=" * 80)
    print("FEATURE EXTRACTION COMPLETE")
    print("=" * 80)
    print(f"Features saved to: {EXTERNAL_DIR}")


def cmd_cluster(args):
    """Full PHATE + Leiden clustering sweep."""
    print("=" * 80)
    print("FUNK ET AL. 2022 — FEATURE EXTRACTION & CLUSTERING")
    print("=" * 80)

    k_values = args.k_values if args.k_values else list(range(5, 16))

    df_genes_i, df_genes_m = _load_hdf_data()
    features_i = _load_feature_list(FEATURES_I_TXT)
    features_m = _load_feature_list(FEATURES_M_TXT)

    CLUSTER_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for df_genes, features, cell_class in [
        (df_genes_i, features_i, "Interphase"),
        (df_genes_m, features_m, "Mitotic"),
    ]:
        gene_symbols, df_cleaned, df_pca, n_pca = extract_features(
            df_genes, features, cell_class
        )
        save_features(gene_symbols, df_cleaned, df_pca, cell_class)

        for k in k_values:
            df_result = run_clustering_for_resolution(df_pca, df_genes, cell_class, k)
            output_dir = CLUSTER_OUTPUT_DIR / f"Funk_{cell_class}_k{k}"
            output_dir.mkdir(parents=True, exist_ok=True)
            output_file = output_dir / "phate_leiden_clustering.tsv"
            df_result.to_csv(output_file, sep="\t", index=False)
            print(f"  Saved: {output_file}")

    print("\n" + "=" * 80)
    print("CLUSTERING COMPLETE")
    print("=" * 80)
    print(f"k values: {k_values}")
    print(f"Features saved to: {EXTERNAL_DIR}")
    print(f"Clustering saved to: {CLUSTER_OUTPUT_DIR}")


def cmd_enrich(args):
    """Run Brieflow enrichment in parallel via ProcessPoolExecutor."""
    from concurrent.futures import ProcessPoolExecutor, as_completed

    print("=" * 80)
    print("FUNK ENRICHMENT ANALYSIS (parallel)")
    print("=" * 80)

    cell_classes = [args.cell_class] if args.cell_class else ["Interphase", "Mitotic"]
    k_values = args.k_values if args.k_values else list(range(5, 16))
    max_workers = args.max_workers

    script_dir = Path(__file__).parent
    results_dir = script_dir / "results" / "cluster"
    config_dir = script_dir.parent.parent / "analysis" / "config" / "benchmark_clusters"

    # Build job list
    jobs = []
    for cell_class in cell_classes:
        for k in k_values:
            jobs.append((cell_class, k, results_dir, config_dir))

    print(f"  Cell classes: {cell_classes}")
    print(f"  k values: {k_values}")
    print(f"  Total jobs: {len(jobs)}")
    print(f"  Max workers: {max_workers}")
    print()

    # Run in parallel
    results = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_run_single_enrichment, job): job for job in jobs}
        for future in as_completed(futures):
            run_name, metrics, error = future.result()
            results[run_name] = (metrics, error)
            if error:
                print(f"  FAILED: {run_name}")
            else:
                string_f1 = metrics.get("STRING", {}).get("f1_score", 0.0)
                corum_pct = (
                    metrics.get("CORUM", {}).get("proportion_enriched", 0.0) * 100
                )
                kegg_pct = metrics.get("KEGG", {}).get("proportion_enriched", 0.0) * 100
                print(
                    f"  {run_name}: STRING F1={string_f1:.4f}  "
                    f"CORUM={corum_pct:.1f}%  KEGG={kegg_pct:.1f}%"
                )

    # Summary table
    print("\n" + "=" * 80)
    print("RESULTS SUMMARY")
    print("=" * 80)
    print(f"\n{'Run':<30s} {'STRING F1':>12s} {'CORUM %':>10s} {'KEGG %':>10s}")
    print("-" * 65)

    for cell_class in cell_classes:
        for k in k_values:
            run_name = f"Funk_{cell_class}_k{k}"
            metrics, error = results.get(run_name, (None, "not run"))
            if metrics:
                string_f1 = metrics.get("STRING", {}).get("f1_score", 0.0)
                corum_pct = (
                    metrics.get("CORUM", {}).get("proportion_enriched", 0.0) * 100
                )
                kegg_pct = metrics.get("KEGG", {}).get("proportion_enriched", 0.0) * 100
                print(
                    f"{run_name:<30s} {string_f1:>12.4f} {corum_pct:>10.1f} {kegg_pct:>10.1f}"
                )
            else:
                print(f"{run_name:<30s} {'FAILED':>12s} {'-':>10s} {'-':>10s}")

    success_count = sum(1 for m, e in results.values() if m is not None)
    print(f"\nSuccessfully processed {success_count}/{len(jobs)} clustering results")
    print("=" * 80)


# ===================================================================
# CLI entry point
# ===================================================================


def main():
    parser = argparse.ArgumentParser(
        description="Funk et al. 2022 pipeline: features, clustering, enrichment"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- features ---
    subparsers.add_parser(
        "features",
        help="Extract cleaned features + PCA (requires ops_clustering env)",
    )

    # --- cluster ---
    p_cluster = subparsers.add_parser(
        "cluster",
        help="PHATE + Leiden clustering sweep (requires ops_clustering env)",
    )
    p_cluster.add_argument(
        "--k-values",
        type=int,
        nargs="+",
        default=None,
        help="Specific k values to cluster at (default: 5-15)",
    )

    # --- enrich ---
    p_enrich = subparsers.add_parser(
        "enrich",
        help="Run Brieflow enrichment in parallel (requires brieflow_aconcagua env)",
    )
    p_enrich.add_argument(
        "--cell-class",
        choices=["Interphase", "Mitotic"],
        default=None,
        help="Single cell class to process (default: both)",
    )
    p_enrich.add_argument(
        "--k-values",
        type=int,
        nargs="+",
        default=None,
        help="Specific k values to enrich (default: 5-15)",
    )
    p_enrich.add_argument(
        "--max-workers",
        type=int,
        default=22,
        help="Max parallel workers (default: 22)",
    )

    args = parser.parse_args()

    if args.command == "features":
        cmd_features(args)
    elif args.command == "cluster":
        cmd_cluster(args)
    elif args.command == "enrich":
        cmd_enrich(args)


if __name__ == "__main__":
    main()
