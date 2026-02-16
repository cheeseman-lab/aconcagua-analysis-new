#!/bin/bash
# =============================================================================
# 13.run_mozzarellm_shuffled.sh - Run mozzarellm on shuffled interphase clusters
# =============================================================================
#
# Runs mozzarellm on the CB-Shuffled clustering (negative control).
# Uses the same config as interphase but reads phate_leiden_clustering_shuffled.tsv.
#
# Usage:
#   bash 13.run_mozzarellm_shuffled.sh
#
# =============================================================================

set -e
cd "$(dirname "$0")"

echo "Running mozzarellm analysis (SHUFFLED negative control)..."
echo "Reading configuration from config/config.yml"
echo ""

python3 << 'MOZZARELLM_SCRIPT'
"""Mozzarellm analysis on shuffled clusters — negative control."""

import sys
from pathlib import Path

import pandas as pd
import yaml
from dotenv import load_dotenv

from mozzarellm import ClusterAnalyzer, reshape_to_clusters

# Load environment variables (.env file for API keys)
load_dotenv(Path(".env"))

# =============================================================================
# Load configuration (same as interphase)
# =============================================================================

CONFIG_PATH = Path("config/config.yml")
with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

if "mozzarellm" not in config:
    print("ERROR: mozzarellm section not found in config.yml")
    sys.exit(1)

mzlm_config = config["mozzarellm"]

# Extract parameters
ROOT_FP = Path(config["all"]["root_fp"])
CELL_CLASS = mzlm_config["cell_class"]
CHANNEL_COMBO = mzlm_config["channel_combo"]
RESOLUTION = mzlm_config["leiden_resolution"]
MODEL = mzlm_config.get("model", "claude-sonnet-4-5-20250929")
TEMPERATURE = mzlm_config.get("temperature", 0.0)
SCREEN_CONTEXT = mzlm_config.get("screen_context", "")
GENE_COL = config["aggregate"]["perturbation_name_col"]

# Build paths — use SHUFFLED input and separate output dir
cluster_path = ROOT_FP / "cluster" / CHANNEL_COMBO / CELL_CLASS / str(RESOLUTION)

cluster_file = cluster_path / "phate_leiden_clustering_shuffled.tsv"
output_dir = cluster_path / "mozzarellm_shuffled"

print(f"Mozzarellm Analysis — SHUFFLED (Negative Control)")
print(f"{'=' * 60}")
print(f"Model: {MODEL}")
print(f"Cell class: {CELL_CLASS}")
print(f"Channel combo: {CHANNEL_COMBO}")
print(f"Resolution: {RESOLUTION}")
print(f"Input: {cluster_file}")
print(f"Output: {output_dir}")
print(f"{'=' * 60}")
print()

if not cluster_file.exists():
    print(f"ERROR: Shuffled clustering file not found: {cluster_file}")
    print(f"Run generate_shuffled_clustering.py first.")
    sys.exit(1)

# =============================================================================
# Run analysis
# =============================================================================

print("Loading shuffled clustering data...")
gene_df = pd.read_csv(cluster_file, sep="\t")

if GENE_COL not in gene_df.columns:
    for alt in ["gene_symbol_0", "gene_symbol", "gene"]:
        if alt in gene_df.columns:
            gene_df = gene_df.rename(columns={alt: GENE_COL})
            break

print(f"Loaded {len(gene_df)} genes across {gene_df['cluster'].nunique()} clusters")

print("Reshaping data to cluster format...")
cluster_df, gene_annotations = reshape_to_clusters(
    input_df=gene_df,
    gene_col=GENE_COL,
    cluster_col="cluster",
    uniprot_col="uniprot_function",
    verbose=True,
)
print(f"Reshaped to {len(cluster_df)} clusters")

print("\nRunning LLM analysis...")
analyzer = ClusterAnalyzer(model=MODEL, temperature=TEMPERATURE, show_progress=True)

results = analyzer.analyze(
    cluster_df,
    gene_annotations=gene_annotations,
    screen_context=SCREEN_CONTEXT,
    output_dir=output_dir,
)

print(f"\nDone!")
print(f"Results saved to: {output_dir}")

MOZZARELLM_SCRIPT
