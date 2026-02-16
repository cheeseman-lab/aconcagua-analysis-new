#!/usr/bin/env python3
"""Generate phate_leiden_clustering_shuffled.tsv from CB-Shuffled eval output.

Takes the CB-Shuffled__combined_table.tsv (which has the updated shuffled cluster
assignments) and joins with the real clustering data to produce the full format
needed by mozzarellm.
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster-dir", required=True, help="Path to cluster dir (e.g. .../Interphase/12)")
    args = parser.parse_args()

    cluster_dir = Path(args.cluster_dir)

    shuffled_table = cluster_dir / "CB-Shuffled__combined_table.tsv"
    real_clustering = cluster_dir / "phate_leiden_clustering.tsv"
    output_file = cluster_dir / "phate_leiden_clustering_shuffled.tsv"

    print(f"Loading CB-Shuffled combined table: {shuffled_table}")
    shuffled_df = pd.read_csv(shuffled_table, sep="\t")

    # Explode semicolon-separated gene lists into one row per gene
    rows = []
    for _, row in shuffled_df.iterrows():
        cluster_id = row["cluster"]
        genes = [g.strip() for g in str(row["genes"]).split(";")]
        for gene in genes:
            rows.append({"gene_symbol_0": gene, "cluster": cluster_id})

    gene_cluster_df = pd.DataFrame(rows)
    print(f"Shuffled assignments: {len(gene_cluster_df)} genes across {gene_cluster_df['cluster'].nunique()} clusters")

    # Load real clustering for metadata
    print(f"Loading real clustering for metadata: {real_clustering}")
    real_df = pd.read_csv(real_clustering, sep="\t")

    # Keep metadata columns from real data (everything except cluster, PHATE coords)
    metadata_cols = [c for c in real_df.columns if c not in ["cluster", "PHATE_0", "PHATE_1",
                                                              "mean_potential_to_nontargeting",
                                                              "normalized_potential_to_nontargeting"]]
    metadata_df = real_df[metadata_cols].copy()

    # Merge shuffled cluster assignments with metadata
    result = gene_cluster_df.merge(metadata_df, on="gene_symbol_0", how="left")

    # Set PHATE coordinates to 0 (not meaningful for shuffled data)
    result["PHATE_0"] = 0.0
    result["PHATE_1"] = 0.0

    # Reorder columns to match expected format
    output_cols = ["gene_symbol_0", "cell_count", "perturbation_auc",
                   "PHATE_0", "PHATE_1", "cluster",
                   "uniprot_entry", "uniprot_function", "uniprot_link"]
    result = result[[c for c in output_cols if c in result.columns]]

    # Drop genes that didn't match (e.g. nontargeting controls not in real data)
    n_before = len(result)
    result = result.dropna(subset=["cell_count"])
    n_after = len(result)
    if n_before != n_after:
        print(f"Dropped {n_before - n_after} genes without metadata (likely not in real clustering)")

    print(f"Output: {len(result)} genes across {result['cluster'].nunique()} clusters")
    result.to_csv(output_file, sep="\t", index=False)
    print(f"Saved to: {output_file}")


if __name__ == "__main__":
    main()
