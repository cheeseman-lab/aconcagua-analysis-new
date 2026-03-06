#!/usr/bin/env python3
"""Generate comprehensive cluster and gene tables (S6-S9) for the Aconcagua analysis.

For each of the 4 datasets (Brieflow Interphase/Mitotic, Funk Interphase/Mitotic),
produces two CSVs:
  - Cluster-level: one row per cluster with summary, confidence, gene counts
  - Gene-level: one row per gene with cluster assignment and classification

Sources:
  - phate_leiden_clustering.tsv (gene -> cluster assignment)
  - claude-opus-4-6_results.json (MozzareLLM annotations, gene classifications)

Naming convention:
  BI = Brieflow Interphase (k=12)
  BM = Brieflow Mitotic (k=5)
  FI = Funk Interphase (k=10)
  FM = Funk Mitotic (k=9)
"""

import csv
import json
from pathlib import Path

# --- Paths ---
ROOT = Path(__file__).resolve().parent.parent
CLUSTER_BASE = ROOT / "analysis" / "brieflow_output" / "cluster" / "DAPI_TUBULIN_GH2AX_PHALLOIDIN"
FUNK_BASE = ROOT / "benchmarks" / "external" / "results" / "cluster"
OUTPUT_DIR = ROOT / "benchmarks" / "results" / "supplemental_tables"

DATASETS = {
    "S6_brieflow_interphase": {
        "clustering": CLUSTER_BASE / "Interphase" / "12" / "phate_leiden_clustering.tsv",
        "results_json": CLUSTER_BASE / "Interphase" / "12" / "mozzarellm" / "claude-opus-4-6_results.json",
        "prefix": "BI",
        "gene_col": "gene_symbol_0",
        "description": "Brieflow Interphase (k=12)",
    },
    "S7_brieflow_mitotic": {
        "clustering": CLUSTER_BASE / "Mitotic" / "5" / "phate_leiden_clustering.tsv",
        "results_json": CLUSTER_BASE / "Mitotic" / "5" / "mozzarellm" / "claude-opus-4-6_results.json",
        "prefix": "BM",
        "gene_col": "gene_symbol_0",
        "description": "Brieflow Mitotic (k=5)",
    },
    "S8_funk_interphase": {
        "clustering": FUNK_BASE / "Funk_Interphase_k10" / "phate_leiden_clustering.tsv",
        "results_json": FUNK_BASE / "Funk_Interphase_k10" / "mozzarellm" / "claude-opus-4-6_results.json",
        "prefix": "FI",
        "gene_col": "gene_symbol",
        "description": "Funk Interphase (k=10)",
    },
    "S9_funk_mitotic": {
        "clustering": FUNK_BASE / "Funk_Mitotic_k9" / "phate_leiden_clustering.tsv",
        "results_json": FUNK_BASE / "Funk_Mitotic_k9" / "mozzarellm" / "claude-opus-4-6_results.json",
        "prefix": "FM",
        "gene_col": "gene_symbol",
        "description": "Funk Mitotic (k=9)",
    },
}

CLUSTER_FIELDS = [
    "cluster_label",
    "cluster_id",
    "dominant_process",
    "pathway_confidence",
    "summary",
    "total_genes",
    "num_established",
    "num_novel",
    "num_uncharacterized",
    "num_unclassified",
    "pct_established",
    "classification_completeness",
]

GENE_FIELDS = [
    "cluster_label",
    "cluster_id",
    "gene_symbol",
    "gene_classification",
    "priority",
    "rationale",
]


def load_clustering(path: Path, gene_col: str) -> dict[str, int]:
    """Read clustering TSV and return gene_symbol -> cluster_id mapping."""
    gene_to_cluster = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row[gene_col].strip()
            cluster = int(row["cluster"].strip())
            gene_to_cluster[gene] = cluster
    return gene_to_cluster


def load_results_json(path: Path) -> dict[int, dict]:
    """Read results JSON and return cluster_id -> cluster info dict."""
    with open(path) as f:
        data = json.load(f)

    clusters = {}
    for cl in data["clusters"]:
        cid = int(cl["cluster_id"])

        # Build gene -> classification mapping
        gene_info = {}

        for gene in cl.get("established_genes", []):
            gene_info[gene] = {"classification": "established", "priority": "", "rationale": ""}

        for item in cl.get("novel_role_genes", []):
            gene_info[item["gene"]] = {
                "classification": "novel_role",
                "priority": item.get("priority", ""),
                "rationale": item.get("rationale", ""),
            }

        for item in cl.get("uncharacterized_genes", []):
            gene_info[item["gene"]] = {
                "classification": "uncharacterized",
                "priority": item.get("priority", ""),
                "rationale": item.get("rationale", ""),
            }

        for gene in cl.get("missed_genes", []):
            gene_info[gene] = {"classification": "unclassified", "priority": "", "rationale": ""}

        clusters[cid] = {
            "dominant_process": cl.get("dominant_process", ""),
            "pathway_confidence": cl.get("pathway_confidence", ""),
            "summary": cl.get("summary", ""),
            "gene_info": gene_info,
            "num_established": len(cl.get("established_genes", [])),
            "num_novel": len(cl.get("novel_role_genes", [])),
            "num_uncharacterized": len(cl.get("uncharacterized_genes", [])),
            "num_unclassified": len(cl.get("missed_genes", [])),
            "classification_completeness": cl.get("classification_completeness", 0),
        }

    return clusters


def generate_tables(config: dict, cluster_path: Path, gene_path: Path):
    """Generate cluster-level and gene-level CSVs."""
    prefix = config["prefix"]

    gene_to_cluster = load_clustering(config["clustering"], config["gene_col"])
    cluster_data = load_results_json(config["results_json"])

    all_cluster_ids = sorted(set(gene_to_cluster.values()))

    cluster_rows = []
    gene_rows = []

    for cid in all_cluster_ids:
        cl = cluster_data.get(cid, {})
        label = f"{prefix}{cid}"
        cluster_genes = sorted([g for g, c in gene_to_cluster.items() if c == cid])
        total = len(cluster_genes)

        n_est = cl.get("num_established", 0)
        n_novel = cl.get("num_novel", 0)
        n_unchar = cl.get("num_uncharacterized", 0)
        n_unclass = cl.get("num_unclassified", 0)
        pct_est = f"{100 * n_est / total:.1f}" if total > 0 else ""
        raw_comp = cl.get("classification_completeness", 0)
        completeness = f"{100 * raw_comp:.1f}%" if isinstance(raw_comp, (int, float)) and raw_comp > 0 else ""

        # Cluster-level row
        cluster_rows.append({
            "cluster_label": label,
            "cluster_id": cid,
            "dominant_process": cl.get("dominant_process", ""),
            "pathway_confidence": cl.get("pathway_confidence", ""),
            "summary": cl.get("summary", ""),
            "total_genes": total,
            "num_established": n_est,
            "num_novel": n_novel,
            "num_uncharacterized": n_unchar,
            "num_unclassified": n_unclass,
            "pct_established": pct_est,
            "classification_completeness": completeness,
        })

        # Gene-level rows
        gene_info = cl.get("gene_info", {})
        for gene in cluster_genes:
            info = gene_info.get(gene, {"classification": "unclassified", "priority": "", "rationale": ""})
            gene_rows.append({
                "cluster_label": label,
                "cluster_id": cid,
                "gene_symbol": gene,
                "gene_classification": info["classification"],
                "priority": info["priority"],
                "rationale": info["rationale"],
            })

    # Write cluster-level CSV
    with open(cluster_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CLUSTER_FIELDS)
        writer.writeheader()
        writer.writerows(cluster_rows)

    # Write gene-level CSV
    with open(gene_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=GENE_FIELDS)
        writer.writeheader()
        writer.writerows(gene_rows)

    return len(cluster_rows), len(gene_rows)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Generating cluster-level and gene-level supplemental tables...")
    print(f"Output directory: {OUTPUT_DIR}\n")

    for name, config in DATASETS.items():
        cluster_path = OUTPUT_DIR / f"{name}_clusters.csv"
        gene_path = OUTPUT_DIR / f"{name}_genes.csv"

        print(f"  {config['description']} ({name}):")

        if not config["clustering"].exists():
            print(f"    ERROR: Clustering file not found!")
            continue
        if not config["results_json"].exists():
            print(f"    ERROR: Results JSON not found!")
            continue

        n_clusters, n_genes = generate_tables(config, cluster_path, gene_path)
        print(f"    -> {cluster_path.name}: {n_clusters} clusters")
        print(f"    -> {gene_path.name}: {n_genes} genes\n")

    # Clean up old combined files
    for name in DATASETS:
        old = OUTPUT_DIR / f"{name}_all_clusters.csv"
        if old.exists():
            old.unlink()
            print(f"  Removed old file: {old.name}")

    print("\nDone.")


if __name__ == "__main__":
    main()
