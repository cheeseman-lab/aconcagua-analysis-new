#!/usr/bin/env python3
"""Validate gene claims in findings_validation_base.md against supplemental tables."""

import csv
import re
from collections import defaultdict
from pathlib import Path

BASE = Path("/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis-new/benchmarks")
FINDINGS = BASE / "findings_validation_base.md"
S6_CLUSTERS = BASE / "results/supplemental_tables/S6_brieflow_interphase_clusters.csv"
S6_GENES = BASE / "results/supplemental_tables/S6_brieflow_interphase_genes.csv"
S7_CLUSTERS = BASE / "results/supplemental_tables/S7_brieflow_mitotic_clusters.csv"
S7_GENES = BASE / "results/supplemental_tables/S7_brieflow_mitotic_genes.csv"


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


def load_clusters_table(path):
    """Load cluster table -> dict of {cluster_label: {field: value}}"""
    clusters = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            clusters[row["cluster_label"]] = row
    return clusters


def parse_findings(path):
    """Parse findings_validation_base.md to extract gene lists per cluster."""
    text = path.read_text()

    # Find all cluster sections: ### N.N Title (CLUSTER_ID)
    # Pattern matches lines like: ### 1.1 Autophagy and PI3K Membrane Trafficking (BI115)
    sections = re.split(r'(?=^### \d)', text, flags=re.MULTILINE)

    findings = {}
    for section in sections:
        # Extract cluster ID from header
        header_match = re.search(r'### \d+\.\d+\s+.*?\(([A-Z]+\d+)\)', section)
        if not header_match:
            continue
        cluster_id = header_match.group(1)

        # Extract gene lists by category
        genes_by_category = {}

        # Pattern: - **Established (N):** GENE1, GENE2, ...
        for category in ["Established", "Novel role", "Uncharacterized"]:
            pattern = rf'\*\*{category}\s*\(\d+\):\*\*\s*(.+?)(?:\n|$)'
            match = re.search(pattern, section)
            if match:
                gene_text = match.group(1).strip()
                if gene_text == '—' or gene_text == '-':
                    genes_by_category[category] = []
                else:
                    genes = [g.strip() for g in gene_text.split(',')]
                    genes_by_category[category] = genes
            else:
                genes_by_category[category] = []

        # Also extract the gene count from header line
        count_match = re.search(r'\*\*(\d+)\s+genes\s*\|', section)
        claimed_count = int(count_match.group(1)) if count_match else None

        findings[cluster_id] = {
            "genes_by_category": genes_by_category,
            "claimed_count": claimed_count,
            "section_text": section[:200],  # first 200 chars for debugging
        }

    return findings


def map_category(classification):
    """Map supplemental table classification to findings category."""
    mapping = {
        "established": "Established",
        "novel_role": "Novel role",
        "uncharacterized": "Uncharacterized",
        "unclassified": None,  # Not listed in findings
    }
    return mapping.get(classification)


def main():
    print("=" * 80)
    print("FINDINGS VALIDATION: Cross-referencing against supplemental tables")
    print("=" * 80)

    # Load supplemental tables
    s6_genes = load_genes_table(S6_GENES)
    s7_genes = load_genes_table(S7_GENES)
    s6_clusters = load_clusters_table(S6_CLUSTERS)
    s7_clusters = load_clusters_table(S7_CLUSTERS)

    # Combine interphase and mitotic
    all_genes = {}
    all_genes.update(s6_genes)
    all_genes.update(s7_genes)

    all_cluster_info = {}
    all_cluster_info.update(s6_clusters)
    all_cluster_info.update(s7_clusters)

    # Parse findings
    findings = parse_findings(FINDINGS)
    print(f"\nFound {len(findings)} clusters in findings document")
    print(f"Supplemental tables: {len(s6_genes)} interphase clusters, {len(s7_genes)} mitotic clusters")

    total_issues = 0
    total_genes_checked = 0
    issues_by_type = defaultdict(list)

    for cluster_id, info in sorted(findings.items()):
        cluster_issues = []

        # Check if cluster exists in supplemental tables
        if cluster_id not in all_genes:
            cluster_issues.append(f"CLUSTER NOT FOUND in supplemental tables")
            issues_by_type["cluster_not_found"].append(cluster_id)
            print(f"\n### {cluster_id}: CLUSTER NOT FOUND IN SUPPLEMENTAL TABLES ###")
            total_issues += 1
            continue

        actual_genes = all_genes[cluster_id]

        # Check gene count
        claimed_count = info["claimed_count"]
        actual_count = len(actual_genes)
        if claimed_count and claimed_count != actual_count:
            cluster_issues.append(f"Gene count mismatch: claimed {claimed_count}, actual {actual_count}")
            issues_by_type["count_mismatch"].append((cluster_id, claimed_count, actual_count))

        # Check each claimed gene
        for category, genes in info["genes_by_category"].items():
            for gene in genes:
                total_genes_checked += 1
                gene = gene.strip()
                if not gene:
                    continue

                if gene not in actual_genes:
                    cluster_issues.append(f"GENE NOT IN CLUSTER: {gene} (claimed {category})")
                    issues_by_type["gene_not_in_cluster"].append((cluster_id, gene, category))
                else:
                    # Check classification
                    actual_class = actual_genes[gene]
                    expected_class = {
                        "Established": "established",
                        "Novel role": "novel_role",
                        "Uncharacterized": "uncharacterized",
                    }.get(category)

                    if expected_class and actual_class != expected_class:
                        cluster_issues.append(
                            f"CLASSIFICATION MISMATCH: {gene} - "
                            f"findings says '{category}', table says '{actual_class}'"
                        )
                        issues_by_type["classification_mismatch"].append(
                            (cluster_id, gene, category, actual_class)
                        )

        # Check for genes in table but not mentioned in findings
        claimed_genes = set()
        for genes in info["genes_by_category"].values():
            claimed_genes.update(g.strip() for g in genes if g.strip())

        # Only flag classified genes that are missing from findings (not unclassified)
        for gene, classification in actual_genes.items():
            if classification != "unclassified" and gene not in claimed_genes:
                # This is expected - findings only lists selected clusters, not all genes
                pass

        if cluster_issues:
            total_issues += len(cluster_issues)
            print(f"\n{'='*60}")
            print(f"ISSUES IN {cluster_id}:")
            for issue in cluster_issues:
                print(f"  - {issue}")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total clusters checked: {len(findings)}")
    print(f"Total genes checked: {total_genes_checked}")
    print(f"Total issues found: {total_issues}")

    if issues_by_type["cluster_not_found"]:
        print(f"\nClusters not found in tables ({len(issues_by_type['cluster_not_found'])}):")
        for c in issues_by_type["cluster_not_found"]:
            print(f"  - {c}")

    if issues_by_type["count_mismatch"]:
        print(f"\nGene count mismatches ({len(issues_by_type['count_mismatch'])}):")
        for c, claimed, actual in issues_by_type["count_mismatch"]:
            print(f"  - {c}: claimed {claimed}, actual {actual}")

    if issues_by_type["gene_not_in_cluster"]:
        print(f"\nGenes not found in cluster ({len(issues_by_type['gene_not_in_cluster'])}):")
        for c, g, cat in issues_by_type["gene_not_in_cluster"]:
            print(f"  - {c}: {g} (claimed {cat})")

    if issues_by_type["classification_mismatch"]:
        print(f"\nClassification mismatches ({len(issues_by_type['classification_mismatch'])}):")
        for c, g, claimed_cat, actual_cat in issues_by_type["classification_mismatch"]:
            print(f"  - {c}: {g} - findings='{claimed_cat}', table='{actual_cat}'")

    if total_issues == 0:
        print("\nAll gene claims validated successfully!")


if __name__ == "__main__":
    main()
