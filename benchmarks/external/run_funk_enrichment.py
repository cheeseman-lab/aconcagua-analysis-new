#!/usr/bin/env python3
"""Run Brieflow enrichment pipeline on a single Funk clustering result.

Usage: python run_funk_enrichment_single.py <cell_class> <k_value>
Example: python run_funk_enrichment_single.py Interphase 10
"""

import sys
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for HPC
import matplotlib.pyplot as plt

# Add Brieflow workflow to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "brieflow" / "workflow"))

from lib.cluster.benchmark_clusters import run_benchmark_analysis, save_json_results


def load_benchmark_data(config_dir):
    """Load STRING, CORUM, KEGG, and UniProt benchmark datasets."""
    string_benchmark = pd.read_csv(config_dir / "string_pair_benchmark.tsv", sep="\t")
    corum_benchmark = pd.read_csv(config_dir / "corum_group_benchmark.tsv", sep="\t")
    kegg_benchmark = pd.read_csv(config_dir / "kegg_group_benchmark.tsv", sep="\t")
    uniprot_data = pd.read_csv(config_dir / "uniprot_data.tsv", sep="\t")

    # Create UniProt lookup by expanding gene_names column
    uniprot_rows = []
    for _, row in uniprot_data.iterrows():
        gene_names = str(row['gene_names']).split()
        for i, gene_name in enumerate(gene_names):
            uniprot_rows.append({
                'gene_name': gene_name,
                'uniprot_entry': row['entry'],
                'uniprot_function': row['function'],
                'uniprot_link': row['link'],
                'position': i
            })

    uniprot_lookup = pd.DataFrame(uniprot_rows)
    uniprot_lookup = (
        uniprot_lookup
        .sort_values('position')
        .drop_duplicates(subset='gene_name', keep='first')
        .drop(columns='position')
    )

    return string_benchmark, corum_benchmark, kegg_benchmark, uniprot_lookup


def convert_funk_to_brieflow_format(funk_df, uniprot_lookup):
    """Convert Funk clustering format to Brieflow format."""
    brieflow_df = funk_df.rename(columns={'gene_symbol': 'gene_symbol_0'})
    brieflow_df['cell_count'] = 1

    brieflow_df = brieflow_df.merge(
        uniprot_lookup,
        left_on='gene_symbol_0',
        right_on='gene_name',
        how='left'
    )

    if 'gene_name' in brieflow_df.columns:
        brieflow_df = brieflow_df.drop(columns='gene_name')

    column_order = [
        'gene_symbol_0', 'cell_count', 'PHATE_0', 'PHATE_1', 'cluster',
        'uniprot_entry', 'uniprot_function', 'uniprot_link'
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
    uniprot_lookup
):
    """Run enrichment analysis on a single Funk clustering result."""
    # Load Funk clustering result
    funk_df = pd.read_csv(clustering_file, sep="\t")

    # Convert to Brieflow format
    brieflow_df = convert_funk_to_brieflow_format(funk_df, uniprot_lookup)

    # Prepare cluster datasets
    cluster_datasets = {"Real": brieflow_df}

    # Run enrichment analysis
    (integrated_results, combined_tables, global_metrics,
     enrichment_pie_charts, enrichment_bar_charts) = run_benchmark_analysis(
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
        output_dir / "CB-Real__integrated_results.json"
    )

    combined_tables["Real"].to_csv(
        output_dir / "CB-Real__combined_table.tsv",
        sep="\t",
        index=False
    )

    save_json_results(
        global_metrics["Real"],
        output_dir / "CB-Real__global_metrics.json"
    )

    enrichment_pie_charts["Real"].savefig(
        output_dir / "CB-Real__pie_chart.png",
        dpi=300,
        bbox_inches='tight'
    )
    plt.close(enrichment_pie_charts["Real"])

    enrichment_bar_charts["Real"].savefig(
        output_dir / "CB-Real__enrichment_bar_chart.png",
        dpi=300,
        bbox_inches='tight'
    )
    plt.close(enrichment_bar_charts["Real"])

    return global_metrics["Real"]


def main():
    """Main enrichment runner for a single clustering result."""
    if len(sys.argv) != 3:
        print("Usage: python run_funk_enrichment_single.py <cell_class> <k_value>")
        print("Example: python run_funk_enrichment_single.py Interphase 10")
        sys.exit(1)

    cell_class = sys.argv[1]
    k_value = int(sys.argv[2])

    if cell_class not in ["Interphase", "Mitotic"]:
        print(f"Error: cell_class must be 'Interphase' or 'Mitotic', got '{cell_class}'")
        sys.exit(1)

    if k_value < 5 or k_value > 15:
        print(f"Error: k_value must be between 5 and 15, got {k_value}")
        sys.exit(1)

    # Set up paths
    script_dir = Path(__file__).parent
    results_dir = script_dir / "results" / "cluster"
    config_dir = script_dir.parent.parent / "analysis" / "config" / "benchmark_clusters"

    run_name = f"Funk_{cell_class}_k{k_value}"
    clustering_file = results_dir / run_name / "phate_leiden_clustering.tsv"

    if not clustering_file.exists():
        print(f"Error: {clustering_file} not found")
        sys.exit(1)

    print(f"Processing {run_name}...")
    print("Loading benchmark data...")
    string_benchmark, corum_benchmark, kegg_benchmark, uniprot_lookup = load_benchmark_data(config_dir)

    print("Running enrichment analysis...")
    output_dir = clustering_file.parent

    try:
        metrics = run_enrichment_for_clustering(
            clustering_file,
            output_dir,
            string_benchmark,
            corum_benchmark,
            kegg_benchmark,
            uniprot_lookup
        )

        # Print results
        string_f1 = metrics.get("STRING", {}).get("f1_score", 0.0)
        corum_pct = metrics.get("CORUM", {}).get("proportion_enriched", 0.0) * 100
        kegg_pct = metrics.get("KEGG", {}).get("proportion_enriched", 0.0) * 100

        print(f"\n✓ {run_name} complete:")
        print(f"  STRING F1: {string_f1:.4f}")
        print(f"  CORUM enriched: {corum_pct:.1f}%")
        print(f"  KEGG enriched: {kegg_pct:.1f}%")

    except Exception as e:
        print(f"\nError processing {run_name}: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
