#!/usr/bin/env python3
"""Generate all supplemental tables (S1-S9) and merged Excel workbook.

Combines:
  S1: All high-confidence clusters
  S2: Gene annotations for high-confidence clusters
  S3: Cross-pipeline Jaccard similarity
  S4: Funk cluster retention in Brieflow (computed inline)
  S5: Shuffled control comparison
  S6-S9: Comprehensive cluster and gene tables (Brieflow + Funk)

Final output: all_supplemental_tables.xlsx with one sheet per table.

Naming convention:
  BI = Brieflow Interphase (k=12)
  BM = Brieflow Mitotic (k=5)
  FI = Funk Interphase (k=10)
  FM = Funk Mitotic (k=9)

Usage:
    python cluster_tables.py
"""

import csv
import json
from collections import Counter
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
BENCHMARKS_DIR = ROOT / "benchmarks"
ANALYSIS_DIR = ROOT / "analysis"
CLUSTER_BASE = ANALYSIS_DIR / "brieflow_output" / "cluster" / "DAPI_TUBULIN_GH2AX_PHALLOIDIN"
EXTERNAL_DIR = BENCHMARKS_DIR / "external"
FUNK_BASE = EXTERNAL_DIR / "results" / "cluster"
OVERLAP_DIR = BENCHMARKS_DIR / "results" / "cluster" / "overlap"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "cluster" / "tables"

# Brieflow data paths
INTERPHASE_SUMMARIES = CLUSTER_BASE / "Interphase" / "12" / "mozzarellm" / "claude-opus-4-6_summaries.tsv"
MITOTIC_SUMMARIES = CLUSTER_BASE / "Mitotic" / "5" / "mozzarellm" / "claude-opus-4-6_summaries.tsv"
INTERPHASE_FLAGGED = CLUSTER_BASE / "Interphase" / "12" / "mozzarellm" / "claude-opus-4-6_flagged_genes.tsv"
MITOTIC_FLAGGED = CLUSTER_BASE / "Mitotic" / "5" / "mozzarellm" / "claude-opus-4-6_flagged_genes.tsv"
SHUFFLED_SUMMARIES = CLUSTER_BASE / "Interphase" / "12" / "mozzarellm_shuffled" / "claude-opus-4-6_summaries.tsv"

# Overlap TSVs (from cluster_overlap.py)
BF_TO_FK_INTER = OVERLAP_DIR / "bf_to_fk_interphase.tsv"
FK_TO_BF_INTER = OVERLAP_DIR / "fk_to_bf_interphase.tsv"
BF_TO_FK_MITO = OVERLAP_DIR / "bf_to_fk_mitotic.tsv"
FK_TO_BF_MITO = OVERLAP_DIR / "fk_to_bf_mitotic.tsv"

# Funk published data (for retention analysis)
FUNK_PUBLISHED_CSV = EXTERNAL_DIR / "1-s2.0-S0092867422013599-mmc4.csv"

# Ideal configurations
IDEAL_CONFIGS = {
    "Interphase": {"brieflow_k": 12, "funk_k": 10},
    "Mitotic": {"brieflow_k": 5, "funk_k": 9},
}

FUNK_MOZZARELLM_CONFIGS = {
    "Interphase": "Funk_Interphase_k10",
    "Mitotic": "Funk_Mitotic_k9",
}

# Comprehensive table datasets (S6-S9)
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

# Funk published cluster definitions (for S4 retention analysis)
FUNK_PAPER_CLUSTERS = {
    "Interphase": [
        (66, "40S ribosome subunits"),
        (23, "60S ribosome subunits"),
        (14, "tRNA ligases & eIF2 translation initiation"),
        (136, "40S ribosome biogenesis"),
        (15, "60S ribosome biogenesis"),
        (21, "Nucleolar proteins / RNA helicases"),
        (112, "Related nucleolar factors (I)"),
        (203, "Related nucleolar factors (II)"),
        (216, "Related nucleolar factors (III)"),
        (192, "TFIID complex"),
        (60, "Mediator / general transcription factors / mRNA export"),
        (8, "Mediator / general transcription factors (II)"),
        (199, "RNA polymerase (I)"),
        (45, "RNA polymerase (II)"),
        (155, "RNA polymerase (III)"),
        (121, "MYC/MAX transcription regulators"),
        (146, "Transcriptional regulation (I)"),
        (39, "Transcriptional regulation (II)"),
        (52, "RNA processing / splicing (I)"),
        (138, "RNA processing / splicing (II)"),
        (110, "RNA processing / splicing (III)"),
        (215, "RNA processing / splicing (IV)"),
        (157, "RNA processing / splicing (V)"),
        (145, "RNA processing / splicing (VI)"),
        (197, "m6A modification (METTL3/METTL14/HNRNPD)"),
        (37, "Integrator complex / mTOR / ER-Golgi"),
        (217, "Integrator subunits (related to cluster 37)"),
        (149, "KRAS/BRAF signaling & mitochondria"),
        (167, "20S proteasome core (AKIRIN2)"),
        (106, "Protein degradation / ubiquitin-proteasome (I)"),
        (213, "Protein degradation / ubiquitin-proteasome (II)"),
        (200, "Protein degradation / ubiquitin-proteasome (III)"),
        (26, "DNA replication & damage response (I)"),
        (13, "DNA replication & damage response (II)"),
        (195, "DNA replication & damage response (III)"),
        (179, "DNA replication & damage response (IV)"),
        (3, "DNA replication & damage response (V)"),
        (148, "Cytokinesis"),
        (218, "Cell cycle control (I)"),
        (95, "Cell cycle control (II)"),
        (46, "Cell cycle control (III)"),
        (29, "Cell adhesion / integrins"),
        (184, "Actin cytoskeleton & adhesion"),
        (104, "Nuclear transport / pore"),
        (201, "Vesicle trafficking (I)"),
        (54, "Vesicle trafficking (II)"),
        (140, "Vesicle trafficking (III)"),
    ],
    "Mitotic": [
        (109, "ZNF335 / spindle bipolarity / gamma-tubulin"),
        (205, "Tubulin"),
        (214, "Augmin complex"),
        (11, "Chromosome alignment"),
        (6, "mRNA splicing (mitotic)"),
        (34, "DNA replication (mitotic)"),
        (88, "Proteasome (mitotic)"),
    ],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def read_tsv(path: Path) -> list[dict]:
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def write_csv(rows: list[dict], path: Path, fieldnames: list[str]) -> None:
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  {path.name}: {len(rows)} rows")


def brieflow_label(cluster_id: str | int, cell_class: str) -> str:
    prefix = "BI" if cell_class.lower().startswith("inter") else "BM"
    return f"{prefix}{cluster_id}"


def funk_label(cluster_id: str | int, cell_class: str) -> str:
    prefix = "FI" if cell_class.lower().startswith("inter") else "FM"
    return f"{prefix}{cluster_id}"


# ---------------------------------------------------------------------------
# S1: All high-confidence clusters
# ---------------------------------------------------------------------------
def generate_s1() -> pd.DataFrame:
    print("S1: All high-confidence clusters")
    rows = []
    for path, cell_class in [(INTERPHASE_SUMMARIES, "Interphase"), (MITOTIC_SUMMARIES, "Mitotic")]:
        for row in read_tsv(path):
            if row["pathway_confidence"] == "High":
                total = int(row["num_established"]) + int(row["num_novel"]) + int(row["num_uncharacterized"])
                rows.append({
                    "cluster_label": brieflow_label(row["cluster_id"], cell_class),
                    "cell_class": cell_class,
                    "cluster_id": int(row["cluster_id"]),
                    "dominant_process": row["dominant_process"],
                    "pathway_confidence": row["pathway_confidence"],
                    "num_established": int(row["num_established"]),
                    "num_novel": int(row["num_novel"]),
                    "num_uncharacterized": int(row["num_uncharacterized"]),
                    "total_genes": total,
                    "pct_established": round(int(row["num_established"]) / max(1, total) * 100, 1),
                    "classification_completeness": row["classification_completeness"],
                })
    rows.sort(key=lambda r: (r["cell_class"], r["cluster_id"]))
    fieldnames = [
        "cluster_label", "cell_class", "cluster_id", "dominant_process",
        "pathway_confidence", "num_established", "num_novel", "num_uncharacterized",
        "total_genes", "pct_established", "classification_completeness",
    ]
    path = OUTPUT_DIR / "S1_all_high_confidence_clusters.csv"
    write_csv(rows, path, fieldnames)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# S2: Gene annotations for high-confidence clusters
# ---------------------------------------------------------------------------
def generate_s2() -> pd.DataFrame:
    print("S2: High-confidence cluster gene annotations")
    hc_clusters: dict[str, set[str]] = {"Interphase": set(), "Mitotic": set()}
    for path, cell_class in [(INTERPHASE_SUMMARIES, "Interphase"), (MITOTIC_SUMMARIES, "Mitotic")]:
        for row in read_tsv(path):
            if row["pathway_confidence"] == "High":
                hc_clusters[cell_class].add(row["cluster_id"])

    rows = []
    for path, cell_class in [(INTERPHASE_FLAGGED, "Interphase"), (MITOTIC_FLAGGED, "Mitotic")]:
        for row in read_tsv(path):
            if row["cluster_id"] in hc_clusters[cell_class]:
                rows.append({
                    "cluster_label": brieflow_label(row["cluster_id"], cell_class),
                    "cell_class": cell_class,
                    "cluster_id": int(row["cluster_id"]),
                    "gene": row["gene"],
                    "category": row["category"],
                    "priority": int(row["priority"]),
                    "rationale": row["rationale"],
                    "dominant_process": row["dominant_process"],
                })
    rows.sort(key=lambda r: (r["cell_class"], r["cluster_id"], -r["priority"], r["gene"]))
    fieldnames = [
        "cluster_label", "cell_class", "cluster_id", "gene", "category",
        "priority", "rationale", "dominant_process",
    ]
    path = OUTPUT_DIR / "S2_high_confidence_gene_annotations.csv"
    write_csv(rows, path, fieldnames)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# S3: Cross-pipeline Jaccard similarity
# ---------------------------------------------------------------------------
def generate_s3() -> pd.DataFrame:
    print("S3: Cross-pipeline Jaccard similarity")
    rows = []
    configs = [
        (BF_TO_FK_INTER, "Interphase", "BF→FK", "BI", "FI"),
        (FK_TO_BF_INTER, "Interphase", "FK→BF", "FI", "BI"),
        (BF_TO_FK_MITO, "Mitotic", "BF→FK", "BM", "FM"),
        (FK_TO_BF_MITO, "Mitotic", "FK→BF", "FM", "BM"),
    ]
    for path, cell_class, direction, src_prefix, tgt_prefix in configs:
        if not path.exists():
            print(f"  WARNING: {path} not found, skipping")
            continue
        for row in read_tsv(path):
            rows.append({
                "cell_class": cell_class,
                "direction": direction,
                "source_label": f"{src_prefix}{row['source_cluster']}",
                "source_process": row["source_process"],
                "source_confidence": row["source_confidence"],
                "source_size": int(row["source_size"]),
                "target_label": f"{tgt_prefix}{row['best_target_cluster']}",
                "target_process": row["target_process"],
                "target_confidence": row["target_confidence"],
                "target_size": int(row["target_size"]),
                "jaccard": round(float(row["jaccard"]), 4),
                "overlap_frac": round(float(row["overlap_frac"]), 4),
                "shared_genes": row["shared_genes"],
                "fragmentation": int(row["fragmentation"]),
                "genes_found_total": int(row["genes_found_total"]),
                "genes_missing": int(row["genes_missing"]) if row["genes_missing"] else 0,
                "category": row["category"],
            })
    rows.sort(key=lambda r: (r["cell_class"], r["direction"], -r["jaccard"]))
    fieldnames = [
        "cell_class", "direction", "source_label", "source_process",
        "source_confidence", "source_size", "target_label", "target_process",
        "target_confidence", "target_size", "jaccard", "overlap_frac",
        "shared_genes", "fragmentation", "genes_found_total", "genes_missing",
        "category",
    ]
    path = OUTPUT_DIR / "S3_cross_pipeline_jaccard.csv"
    write_csv(rows, path, fieldnames)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# S4: Funk cluster retention (folded in from cluster_retention.py)
# ---------------------------------------------------------------------------
def _load_funk_published() -> pd.DataFrame:
    df = pd.read_csv(FUNK_PUBLISHED_CSV)
    df = df.rename(columns={"\ufeffGene symbol": "Gene symbol"})
    if "Gene symbol" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Gene symbol"})
    return df


def _load_brieflow_clustered(cell_class: str) -> pd.DataFrame:
    k = IDEAL_CONFIGS[cell_class]["brieflow_k"]
    path = CLUSTER_BASE / cell_class / str(k) / "phate_leiden_clustering.tsv"
    return pd.read_csv(path, sep="\t")


def _load_mozzarellm_summaries(cell_class: str) -> pd.DataFrame:
    k = IDEAL_CONFIGS[cell_class]["brieflow_k"]
    path = CLUSTER_BASE / cell_class / str(k) / "mozzarellm" / "claude-opus-4-6_summaries.tsv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def _load_funk_mozzarellm_summaries(cell_class: str) -> pd.DataFrame:
    path = FUNK_BASE / FUNK_MOZZARELLM_CONFIGS[cell_class] / "mozzarellm" / "claude-opus-4-6_summaries.tsv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def _analyze_retention_cluster(funk_cluster_id, cell_class, funk_df, brieflow_df, mozzarellm_df, funk_mozzarellm_df):
    """Analyze how one Funk cluster maps onto Brieflow clustering."""
    col = "Interphase cluster" if cell_class == "Interphase" else "Mitotic cluster"
    mask = pd.to_numeric(funk_df[col], errors="coerce") == funk_cluster_id
    genes = funk_df.loc[mask, "Gene symbol"].dropna().tolist()

    bl_lookup = dict(zip(brieflow_df["gene_symbol_0"], brieflow_df["cluster"]))

    found = []
    for gene in genes:
        bl_cl = bl_lookup.get(gene)
        if bl_cl is not None:
            found.append(bl_cl)

    n_found = len(found)
    if n_found == 0:
        return None

    distribution = Counter(found)
    dominant_cluster, dominant_count = distribution.most_common(1)[0]
    preservation = dominant_count / n_found

    # Get mozzarellm annotations
    dom_ann, dom_conf = "", ""
    if not mozzarellm_df.empty:
        match = mozzarellm_df[mozzarellm_df["cluster_id"] == int(dominant_cluster)]
        if not match.empty:
            dom_ann = match.iloc[0].get("dominant_process", "")
            dom_conf = match.iloc[0].get("pathway_confidence", "")

    funk_ann, funk_conf = "", ""
    if not funk_mozzarellm_df.empty:
        match = funk_mozzarellm_df[funk_mozzarellm_df["cluster_id"] == funk_cluster_id]
        if not match.empty:
            funk_ann = match.iloc[0].get("dominant_process", "")
            funk_conf = match.iloc[0].get("pathway_confidence", "")

    dist_parts = [f"{int(cl)}({cnt})" for cl, cnt in distribution.most_common()]

    return {
        "n_genes": len(genes),
        "n_found": n_found,
        "preservation": preservation,
        "funk_annotation": funk_ann,
        "funk_confidence": funk_conf,
        "dominant_brieflow_cluster": int(dominant_cluster),
        "dominant_count": dominant_count,
        "dominant_annotation": dom_ann,
        "dominant_confidence": dom_conf,
        "n_brieflow_clusters": len(distribution),
        "distribution": ", ".join(dist_parts),
    }


def _assess(preservation):
    if preservation >= 0.8:
        return "well"
    elif preservation >= 0.5:
        return "partially"
    else:
        return "poorly"


def generate_s4() -> pd.DataFrame:
    """S4: Funk cluster retention in Brieflow (computed inline)."""
    print("S4: Funk cluster retention in Brieflow")
    funk_df = _load_funk_published()

    brieflow_data = {}
    mozzarellm_data = {}
    funk_mozzarellm_data = {}
    for cell_class in ["Interphase", "Mitotic"]:
        brieflow_data[cell_class] = _load_brieflow_clustered(cell_class)
        mozzarellm_data[cell_class] = _load_mozzarellm_summaries(cell_class)
        funk_mozzarellm_data[cell_class] = _load_funk_mozzarellm_summaries(cell_class)

    rows = []
    for cell_class, clusters in FUNK_PAPER_CLUSTERS.items():
        for funk_id, desc in clusters:
            summary = _analyze_retention_cluster(
                funk_id, cell_class,
                funk_df, brieflow_data[cell_class],
                mozzarellm_data[cell_class],
                funk_mozzarellm_data[cell_class],
            )

            funk_lbl = funk_label(funk_id, cell_class)
            if summary is None:
                rows.append({
                    "funk_label": funk_lbl,
                    "cell_class": cell_class,
                    "funk_cluster": funk_id,
                    "funk_description": desc,
                    "assessment": "no genes found",
                    "n_genes": 0,
                    "n_found": 0,
                    "preservation": 0.0,
                    "funk_annotation": "",
                    "funk_confidence": "",
                    "dominant_brieflow_label": "",
                    "dominant_count": 0,
                    "dominant_annotation": "",
                    "dominant_confidence": "",
                    "n_brieflow_clusters": 0,
                    "distribution": "",
                })
            else:
                bl_prefix = "BI" if cell_class.lower().startswith("inter") else "BM"
                dom_bl_label = f"{bl_prefix}{summary['dominant_brieflow_cluster']}"
                rows.append({
                    "funk_label": funk_lbl,
                    "cell_class": cell_class,
                    "funk_cluster": funk_id,
                    "funk_description": desc,
                    "assessment": _assess(summary["preservation"]),
                    "n_genes": summary["n_genes"],
                    "n_found": summary["n_found"],
                    "preservation": round(summary["preservation"], 4),
                    "funk_annotation": summary["funk_annotation"],
                    "funk_confidence": summary["funk_confidence"],
                    "dominant_brieflow_label": dom_bl_label,
                    "dominant_count": summary["dominant_count"],
                    "dominant_annotation": summary["dominant_annotation"],
                    "dominant_confidence": summary["dominant_confidence"],
                    "n_brieflow_clusters": summary["n_brieflow_clusters"],
                    "distribution": summary["distribution"],
                })

    rows.sort(key=lambda r: (r["cell_class"], -r["preservation"]))
    fieldnames = [
        "funk_label", "cell_class", "funk_cluster", "funk_description",
        "assessment", "n_genes", "n_found", "preservation",
        "funk_annotation", "funk_confidence", "dominant_brieflow_label",
        "dominant_count", "dominant_annotation", "dominant_confidence",
        "n_brieflow_clusters", "distribution",
    ]
    path = OUTPUT_DIR / "S4_funk_cluster_retention.csv"
    write_csv(rows, path, fieldnames)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# S5: Shuffled control comparison
# ---------------------------------------------------------------------------
def generate_s5() -> pd.DataFrame:
    print("S5: Shuffled control summary")
    datasets = {
        "Brieflow Interphase (real)": INTERPHASE_SUMMARIES,
        "Brieflow Mitotic (real)": MITOTIC_SUMMARIES,
        "Brieflow Interphase (shuffled)": SHUFFLED_SUMMARIES,
    }

    rows = []
    for dataset_name, path in datasets.items():
        data = read_tsv(path)
        total = len(data)
        counts = {"High": 0, "Medium": 0, "Low": 0}
        gene_counts = {"High": 0, "Medium": 0, "Low": 0}
        for row in data:
            conf = row["pathway_confidence"]
            counts[conf] = counts.get(conf, 0) + 1
            n_genes = int(row["num_established"]) + int(row["num_novel"]) + int(row["num_uncharacterized"])
            gene_counts[conf] = gene_counts.get(conf, 0) + n_genes

        for conf in ["High", "Medium", "Low"]:
            rows.append({
                "dataset": dataset_name,
                "confidence_level": conf,
                "n_clusters": counts[conf],
                "pct_clusters": round(counts[conf] / max(1, total) * 100, 1),
                "n_genes_in_clusters": gene_counts[conf],
                "total_clusters": total,
            })

    fieldnames = [
        "dataset", "confidence_level", "n_clusters", "pct_clusters",
        "n_genes_in_clusters", "total_clusters",
    ]
    path = OUTPUT_DIR / "S5_shuffled_control_summary.csv"
    write_csv(rows, path, fieldnames)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# S6-S9: Comprehensive cluster and gene tables
# ---------------------------------------------------------------------------
CLUSTER_FIELDS = [
    "cluster_label", "cluster_id", "dominant_process", "pathway_confidence",
    "summary", "total_genes", "num_established", "num_novel",
    "num_uncharacterized", "num_unclassified", "pct_established",
    "classification_completeness",
]

GENE_FIELDS = [
    "cluster_label", "cluster_id", "gene_symbol", "gene_classification",
    "priority", "rationale",
]


def _load_clustering(path: Path, gene_col: str) -> dict[str, int]:
    gene_to_cluster = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row[gene_col].strip()
            cluster = int(row["cluster"].strip())
            gene_to_cluster[gene] = cluster
    return gene_to_cluster


def _load_results_json(path: Path) -> dict[int, dict]:
    with open(path) as f:
        data = json.load(f)

    clusters = {}
    for cl in data["clusters"]:
        cid = int(cl["cluster_id"])
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


def generate_comprehensive_tables(name: str, config: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate cluster-level and gene-level tables for one dataset."""
    prefix = config["prefix"]
    gene_to_cluster = _load_clustering(config["clustering"], config["gene_col"])
    cluster_data = _load_results_json(config["results_json"])
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

        cluster_rows.append({
            "cluster_label": label, "cluster_id": cid,
            "dominant_process": cl.get("dominant_process", ""),
            "pathway_confidence": cl.get("pathway_confidence", ""),
            "summary": cl.get("summary", ""),
            "total_genes": total,
            "num_established": n_est, "num_novel": n_novel,
            "num_uncharacterized": n_unchar, "num_unclassified": n_unclass,
            "pct_established": pct_est, "classification_completeness": completeness,
        })

        gene_info = cl.get("gene_info", {})
        for gene in cluster_genes:
            info = gene_info.get(gene, {"classification": "unclassified", "priority": "", "rationale": ""})
            gene_rows.append({
                "cluster_label": label, "cluster_id": cid,
                "gene_symbol": gene,
                "gene_classification": info["classification"],
                "priority": info["priority"],
                "rationale": info["rationale"],
            })

    cluster_path = OUTPUT_DIR / f"{name}_clusters.csv"
    gene_path = OUTPUT_DIR / f"{name}_genes.csv"

    with open(cluster_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CLUSTER_FIELDS)
        writer.writeheader()
        writer.writerows(cluster_rows)

    with open(gene_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=GENE_FIELDS)
        writer.writeheader()
        writer.writerows(gene_rows)

    print(f"  {cluster_path.name}: {len(cluster_rows)} clusters")
    print(f"  {gene_path.name}: {len(gene_rows)} genes")

    return pd.DataFrame(cluster_rows), pd.DataFrame(gene_rows)


def generate_s6_to_s9() -> dict[str, pd.DataFrame]:
    """Generate S6-S9 comprehensive tables."""
    print("S6-S9: Comprehensive cluster and gene tables")
    all_dfs = {}
    for name, config in DATASETS.items():
        if not config["clustering"].exists():
            print(f"  WARNING: {config['clustering']} not found, skipping {name}")
            continue
        if not config["results_json"].exists():
            print(f"  WARNING: {config['results_json']} not found, skipping {name}")
            continue
        cluster_df, gene_df = generate_comprehensive_tables(name, config)
        all_dfs[f"{name}_clusters"] = cluster_df
        all_dfs[f"{name}_genes"] = gene_df
    return all_dfs


# ---------------------------------------------------------------------------
# Excel export
# ---------------------------------------------------------------------------
def export_excel(all_tables: dict[str, pd.DataFrame]) -> None:
    """Write all tables to a single Excel workbook with one sheet per table."""
    excel_path = OUTPUT_DIR / "all_supplemental_tables.xlsx"
    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        for sheet_name, df in all_tables.items():
            # Excel sheet names max 31 chars
            short_name = sheet_name[:31]
            df.to_excel(writer, sheet_name=short_name, index=False)
    print(f"\nExcel workbook: {excel_path.name} ({len(all_tables)} sheets)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}\n")

    all_tables: dict[str, pd.DataFrame] = {}

    all_tables["S1_high_conf_clusters"] = generate_s1()
    all_tables["S2_gene_annotations"] = generate_s2()
    all_tables["S3_jaccard_similarity"] = generate_s3()
    all_tables["S4_funk_retention"] = generate_s4()
    all_tables["S5_shuffled_control"] = generate_s5()

    comprehensive = generate_s6_to_s9()
    all_tables.update(comprehensive)

    export_excel(all_tables)

    print(f"\nDone. All tables written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
