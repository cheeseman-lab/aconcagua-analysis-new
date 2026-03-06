#!/usr/bin/env python3
"""Detailed investigation of gene discrepancies between findings and supplemental tables."""

import csv
import re
from collections import defaultdict
from pathlib import Path

BASE = Path("/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis-new/benchmarks")
FINDINGS = BASE / "findings_validation_base.md"
S6_GENES = BASE / "results/supplemental_tables/S6_brieflow_interphase_genes.csv"
S7_GENES = BASE / "results/supplemental_tables/S7_brieflow_mitotic_genes.csv"

# Also check MozzareLLM source data
MOZZARELLM_INTER = Path("/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis-new/analysis/brieflow_output/cluster/DAPI_TUBULIN_GH2AX_PHALLOIDIN/Interphase/12/mozzarellm/claude-opus-4-6_summaries.tsv")
MOZZARELLM_MITO = Path("/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis-new/analysis/brieflow_output/cluster/DAPI_TUBULIN_GH2AX_PHALLOIDIN/Mitotic/5/mozzarellm/claude-opus-4-6_summaries.tsv")
MOZZARELLM_FLAGGED_INTER = Path("/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis-new/analysis/brieflow_output/cluster/DAPI_TUBULIN_GH2AX_PHALLOIDIN/Interphase/12/mozzarellm/claude-opus-4-6_flagged_genes.tsv")
MOZZARELLM_FLAGGED_MITO = Path("/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis-new/analysis/brieflow_output/cluster/DAPI_TUBULIN_GH2AX_PHALLOIDIN/Mitotic/5/mozzarellm/claude-opus-4-6_flagged_genes.tsv")


def load_genes_table(path):
    """Load gene table -> dict of {cluster_label: {gene: classification}}"""
    clusters = defaultdict(dict)
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            label = row["cluster_label"]
            gene = row["gene_symbol"]
            classification = row["gene_classification"]
            clusters[label][gene] = classification
    return clusters


def build_gene_index(genes_tables):
    """Build reverse index: gene -> list of clusters it appears in."""
    index = defaultdict(list)
    for cluster_label, genes in genes_tables.items():
        for gene in genes:
            index[gene].append(cluster_label)
    return index


def load_mozzarellm_summaries(path):
    """Load MozzareLLM summaries -> dict of {cluster_id: row}"""
    if not path.exists():
        return {}
    clusters = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            clusters[row.get("cluster_id", row.get("cluster", ""))] = row
    return clusters


def load_mozzarellm_flagged(path):
    """Load MozzareLLM flagged genes -> dict of {cluster_id: {gene: classification}}"""
    if not path.exists():
        return {}
    clusters = defaultdict(dict)
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cid = row.get("cluster_id", row.get("cluster", ""))
            gene = row.get("gene_symbol", row.get("gene", ""))
            classification = row.get("gene_classification", row.get("classification", ""))
            if gene:
                clusters[cid][gene] = classification
    return clusters


def parse_findings(path):
    """Parse findings_validation_base.md to extract gene lists per cluster."""
    text = path.read_text()
    sections = re.split(r'(?=^### \d)', text, flags=re.MULTILINE)
    findings = {}
    for section in sections:
        header_match = re.search(r'### \d+\.\d+\s+.*?\(([A-Z]+\d+)\)', section)
        if not header_match:
            continue
        cluster_id = header_match.group(1)
        genes_by_category = {}
        for category in ["Established", "Novel role", "Uncharacterized"]:
            pattern = rf'\*\*{category}\s*\(\d+\):\*\*\s*(.+?)(?:\n|$)'
            match = re.search(pattern, section)
            if match:
                gene_text = match.group(1).strip()
                if gene_text in ('—', '-'):
                    genes_by_category[category] = []
                else:
                    genes = [g.strip() for g in gene_text.split(',')]
                    genes_by_category[category] = genes
            else:
                genes_by_category[category] = []
        count_match = re.search(r'\*\*(\d+)\s+genes\s*\|', section)
        claimed_count = int(count_match.group(1)) if count_match else None
        findings[cluster_id] = {
            "genes_by_category": genes_by_category,
            "claimed_count": claimed_count,
        }
    return findings


def main():
    # Load all data sources
    s6_genes = load_genes_table(S6_GENES)
    s7_genes = load_genes_table(S7_GENES)
    all_genes = {}
    all_genes.update(s6_genes)
    all_genes.update(s7_genes)

    gene_index = build_gene_index(all_genes)

    # Load mozzarellm data
    moz_inter_flagged = load_mozzarellm_flagged(MOZZARELLM_FLAGGED_INTER)
    moz_mito_flagged = load_mozzarellm_flagged(MOZZARELLM_FLAGGED_MITO)

    findings = parse_findings(FINDINGS)

    # Focus on Brieflow clusters (BI* and BM*) that have issues
    print("=" * 80)
    print("DETAILED INVESTIGATION: Where do 'missing' genes actually live?")
    print("=" * 80)

    brieflow_clusters_with_issues = []
    for cluster_id, info in sorted(findings.items()):
        if not cluster_id.startswith(("BI", "BM")):
            continue
        if cluster_id not in all_genes:
            continue

        actual_genes = all_genes[cluster_id]
        all_claimed = set()
        for genes in info["genes_by_category"].values():
            all_claimed.update(g.strip() for g in genes if g.strip())

        missing = all_claimed - set(actual_genes.keys())
        if missing:
            brieflow_clusters_with_issues.append((cluster_id, missing, info))

    for cluster_id, missing_genes, info in brieflow_clusters_with_issues:
        actual_genes = all_genes[cluster_id]
        print(f"\n{'='*60}")
        print(f"CLUSTER {cluster_id}: {len(missing_genes)} missing genes")
        print(f"  Claimed gene count: {info['claimed_count']}")
        print(f"  Actual gene count in table: {len(actual_genes)}")
        print(f"  Claimed classified genes: {sum(len(v) for v in info['genes_by_category'].values())}")

        # Check mozzarellm source for this cluster
        cid = cluster_id[2:]  # Remove BI/BM prefix to get numeric cluster id
        if cluster_id.startswith("BI"):
            moz_genes = moz_inter_flagged.get(cid, {})
        else:
            moz_genes = moz_mito_flagged.get(cid, {})

        print(f"  MozzareLLM flagged genes for cluster {cid}: {len(moz_genes)}")

        for gene in sorted(missing_genes):
            # Where is this gene in the supplemental tables?
            locations = gene_index.get(gene, [])
            # Is it in MozzareLLM source?
            in_moz = gene in moz_genes
            moz_class = moz_genes.get(gene, "N/A")

            if locations:
                print(f"  MISSING: {gene} -> actually in clusters: {', '.join(locations)} | mozzarellm: {moz_class}")
            else:
                print(f"  MISSING: {gene} -> NOT FOUND IN ANY CLUSTER | mozzarellm: {moz_class}")

        # Also show which genes ARE in the actual cluster but not claimed
        actual_classified = {g: c for g, c in actual_genes.items() if c != "unclassified"}
        claimed_all = set()
        for genes in info["genes_by_category"].values():
            claimed_all.update(g.strip() for g in genes if g.strip())
        in_table_not_claimed = set(actual_classified.keys()) - claimed_all
        if in_table_not_claimed:
            print(f"\n  In table but NOT in findings ({len(in_table_not_claimed)}):")
            for g in sorted(in_table_not_claimed):
                print(f"    + {g} ({actual_classified[g]})")


if __name__ == "__main__":
    main()
