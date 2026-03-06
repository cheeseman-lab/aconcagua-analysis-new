#!/usr/bin/env python3
"""Run mozzarellm on Funk et al. 2022 clustering results.

Loads Funk clustering TSV, adds UniProt annotations, and runs mozzarellm.
Requires: conda activate mozzarellm

Usage:
    python run_mozzarellm_funk.py --cell-class Interphase --k 10
    python run_mozzarellm_funk.py --cell-class Mitotic --k 9
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

from mozzarellm import ClusterAnalyzer, reshape_to_clusters

# Load API keys
SCRIPT_DIR = Path(__file__).parent
ANALYSIS_DIR = SCRIPT_DIR.parent / "analysis"
load_dotenv(ANALYSIS_DIR / ".env")

# Paths
CLUSTER_DIR = SCRIPT_DIR / "external" / "results" / "cluster"
CONFIG_DIR = ANALYSIS_DIR / "config" / "benchmark_clusters"

MODEL = "claude-opus-4-6"
TEMPERATURE = 0.0

SCREEN_CONTEXT_TEMPLATE = (
    "These clusters are from an optical pooled screen (OPS) that measured morphological "
    "phenotypes in human cells. "
    "The clusters are derived from {cell_class_desc}. "
    "The screen involved perturbing genes using CRISPR knockout and imaging "
    "the resulting cellular morphology via fluorescence microscopy using a "
    "cell-painting-derived panel, specifically DAPI, Tubulin, \u03b3H2AX (DNA damage), "
    "and Phalloidin. Genes grouped within a cluster tend to exhibit similar phenotypes, "
    "suggesting they may participate in the same biological process or pathway."
)

CELL_CLASS_DESCRIPTIONS = {
    "Interphase": (
        "interphase cells, so phenotypes related to cell division "
        "and mitotic processes are not expected"
    ),
    "Mitotic": (
        "mitotic cells, so phenotypes related to cell division, spindle assembly, "
        "and chromosome segregation are expected to be prominent"
    ),
}


def load_uniprot_lookup():
    """Load UniProt data and create gene-name lookup table."""
    uniprot_data = pd.read_csv(CONFIG_DIR / "uniprot_data.tsv", sep="\t")

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

    lookup = pd.DataFrame(uniprot_rows)
    lookup = (
        lookup.sort_values("position")
        .drop_duplicates(subset="gene_name", keep="first")
        .drop(columns="position")
    )
    return lookup


def convert_funk_to_mozzarellm_format(funk_df, uniprot_lookup):
    """Convert Funk clustering to format expected by mozzarellm."""
    df = funk_df.rename(columns={"gene_symbol": "gene_symbol_0"})

    df = df.merge(
        uniprot_lookup,
        left_on="gene_symbol_0",
        right_on="gene_name",
        how="left",
    )
    if "gene_name" in df.columns:
        df = df.drop(columns="gene_name")

    return df


def main():
    parser = argparse.ArgumentParser(description="Run mozzarellm on Funk clustering")
    parser.add_argument(
        "--cell-class",
        required=True,
        choices=["Interphase", "Mitotic"],
    )
    parser.add_argument("--k", type=int, required=True, help="Leiden resolution k")
    args = parser.parse_args()

    cell_class = args.cell_class
    k = args.k

    # Paths
    run_name = f"Funk_{cell_class}_k{k}"
    cluster_file = CLUSTER_DIR / run_name / "phate_leiden_clustering.tsv"
    output_dir = CLUSTER_DIR / run_name / "mozzarellm"

    print(f"Mozzarellm Analysis — Funk {cell_class} k={k}")
    print(f"{'=' * 60}")
    print(f"Model: {MODEL}")
    print(f"Input: {cluster_file}")
    print(f"Output: {output_dir}")
    print(f"{'=' * 60}")
    print()

    if not cluster_file.exists():
        print(f"ERROR: Clustering file not found: {cluster_file}")
        sys.exit(1)

    # Load and convert
    print("Loading clustering data...")
    funk_df = pd.read_csv(cluster_file, sep="\t")
    print(f"Loaded {len(funk_df)} genes across {funk_df['cluster'].nunique()} clusters")

    print("Loading UniProt annotations...")
    uniprot_lookup = load_uniprot_lookup()

    print("Converting to mozzarellm format...")
    gene_df = convert_funk_to_mozzarellm_format(funk_df, uniprot_lookup)

    # Reshape to cluster format
    print("Reshaping data to cluster format...")
    cluster_df, gene_annotations = reshape_to_clusters(
        input_df=gene_df,
        gene_col="gene_symbol_0",
        cluster_col="cluster",
        uniprot_col="uniprot_function",
        verbose=True,
    )
    print(f"Reshaped to {len(cluster_df)} clusters")

    # Build screen context
    screen_context = SCREEN_CONTEXT_TEMPLATE.format(
        cell_class_desc=CELL_CLASS_DESCRIPTIONS[cell_class]
    )

    # Run analysis
    print("\nRunning LLM analysis...")
    analyzer = ClusterAnalyzer(model=MODEL, temperature=TEMPERATURE, show_progress=True)

    analyzer.analyze(
        cluster_df,
        gene_annotations=gene_annotations,
        screen_context=screen_context,
        output_dir=output_dir,
    )

    print("\nDone!")
    print(f"Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
