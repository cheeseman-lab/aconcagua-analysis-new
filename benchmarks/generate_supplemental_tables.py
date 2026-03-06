#!/usr/bin/env python3
"""Generate supplemental tables (S1-S5) for the Aconcagua analysis.

Reads source TSVs from MozzareLLM output and cluster overlap analysis,
produces clean CSVs with consistent BI/BM/FI/FM cluster labeling.

Naming convention:
  BI = Brieflow Interphase (k=12)
  BM = Brieflow Mitotic (k=5)
  FI = Funk Interphase (k=10)
  FM = Funk Mitotic (k=9)
"""

import csv
import os
from pathlib import Path

# --- Paths ---
ROOT = Path(__file__).resolve().parent.parent
CLUSTER_BASE = ROOT / "analysis" / "brieflow_output" / "cluster" / "DAPI_TUBULIN_GH2AX_PHALLOIDIN"

INTERPHASE_SUMMARIES = CLUSTER_BASE / "Interphase" / "12" / "mozzarellm" / "claude-opus-4-6_summaries.tsv"
MITOTIC_SUMMARIES = CLUSTER_BASE / "Mitotic" / "5" / "mozzarellm" / "claude-opus-4-6_summaries.tsv"
INTERPHASE_FLAGGED = CLUSTER_BASE / "Interphase" / "12" / "mozzarellm" / "claude-opus-4-6_flagged_genes.tsv"
MITOTIC_FLAGGED = CLUSTER_BASE / "Mitotic" / "5" / "mozzarellm" / "claude-opus-4-6_flagged_genes.tsv"
SHUFFLED_SUMMARIES = CLUSTER_BASE / "Interphase" / "12" / "mozzarellm_shuffled" / "claude-opus-4-6_summaries.tsv"

OVERLAP_DIR = ROOT / "benchmarks" / "results" / "cluster_overlap"
BF_TO_FK_INTER = OVERLAP_DIR / "bf_to_fk_interphase.tsv"
FK_TO_BF_INTER = OVERLAP_DIR / "fk_to_bf_interphase.tsv"
BF_TO_FK_MITO = OVERLAP_DIR / "bf_to_fk_mitotic.tsv"
FK_TO_BF_MITO = OVERLAP_DIR / "fk_to_bf_mitotic.tsv"

RETENTION_TSV = ROOT / "benchmarks" / "results" / "cluster" / "retention" / "funk_cluster_retention.tsv"

OUTPUT_DIR = ROOT / "benchmarks" / "results" / "supplemental_tables"


def read_tsv(path: Path) -> list[dict]:
    """Read a TSV file and return list of dicts."""
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def write_csv(rows: list[dict], path: Path, fieldnames: list[str]) -> None:
    """Write rows to CSV with given fieldnames."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Wrote {len(rows)} rows -> {path.name}")


def brieflow_label(cluster_id: str | int, cell_class: str) -> str:
    """Convert cluster_id + cell_class to BI/BM label."""
    prefix = "BI" if cell_class.lower().startswith("inter") else "BM"
    return f"{prefix}{cluster_id}"


def funk_label(cluster_id: str | int, cell_class: str) -> str:
    """Convert cluster_id + cell_class to FI/FM label."""
    prefix = "FI" if cell_class.lower().startswith("inter") else "FM"
    return f"{prefix}{cluster_id}"


def generate_s1() -> None:
    """S1: All high-confidence clusters across both cell classes."""
    print("S1: All high-confidence clusters")
    rows = []

    for path, cell_class in [(INTERPHASE_SUMMARIES, "Interphase"), (MITOTIC_SUMMARIES, "Mitotic")]:
        for row in read_tsv(path):
            if row["pathway_confidence"] == "High":
                rows.append({
                    "cluster_label": brieflow_label(row["cluster_id"], cell_class),
                    "cell_class": cell_class,
                    "cluster_id": int(row["cluster_id"]),
                    "dominant_process": row["dominant_process"],
                    "pathway_confidence": row["pathway_confidence"],
                    "num_established": int(row["num_established"]),
                    "num_novel": int(row["num_novel"]),
                    "num_uncharacterized": int(row["num_uncharacterized"]),
                    "total_genes": int(row["num_established"]) + int(row["num_novel"]) + int(row["num_uncharacterized"]),
                    "pct_established": round(
                        int(row["num_established"]) / max(1, int(row["num_established"]) + int(row["num_novel"]) + int(row["num_uncharacterized"])) * 100, 1
                    ),
                    "classification_completeness": row["classification_completeness"],
                })

    rows.sort(key=lambda r: (r["cell_class"], r["cluster_id"]))
    fieldnames = [
        "cluster_label", "cell_class", "cluster_id", "dominant_process",
        "pathway_confidence", "num_established", "num_novel", "num_uncharacterized",
        "total_genes", "pct_established", "classification_completeness",
    ]
    write_csv(rows, OUTPUT_DIR / "S1_all_high_confidence_clusters.csv", fieldnames)


def generate_s2() -> None:
    """S2: Gene annotations for all high-confidence clusters."""
    print("S2: High-confidence cluster gene annotations")

    # First collect HC cluster IDs per cell class
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
    write_csv(rows, OUTPUT_DIR / "S2_high_confidence_gene_annotations.csv", fieldnames)


def generate_s3() -> None:
    """S3: Cross-pipeline Jaccard similarity for all HC cluster pairs."""
    print("S3: Cross-pipeline Jaccard similarity")
    rows = []

    configs = [
        (BF_TO_FK_INTER, "Interphase", "BF→FK", "BI", "FI"),
        (FK_TO_BF_INTER, "Interphase", "FK→BF", "FI", "BI"),
        (BF_TO_FK_MITO, "Mitotic", "BF→FK", "BM", "FM"),
        (FK_TO_BF_MITO, "Mitotic", "FK→BF", "FM", "BM"),
    ]

    for path, cell_class, direction, src_prefix, tgt_prefix in configs:
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
    write_csv(rows, OUTPUT_DIR / "S3_cross_pipeline_jaccard.csv", fieldnames)


def generate_s4() -> None:
    """S4: Funk cluster retention in Brieflow pipeline."""
    print("S4: Funk cluster retention")
    rows = []

    for row in read_tsv(RETENTION_TSV):
        cell_class = row["cell_class"]
        funk_id = row["funk_cluster"]
        funk_lbl = funk_label(funk_id, cell_class)

        # Parse dominant brieflow cluster
        dom_bl = row.get("dominant_brieflow_cluster", "")
        dom_bl_label = ""
        if dom_bl:
            bl_prefix = "BI" if cell_class.lower().startswith("inter") else "BM"
            dom_bl_label = f"{bl_prefix}{dom_bl}"

        rows.append({
            "funk_label": funk_lbl,
            "cell_class": cell_class,
            "funk_cluster": int(funk_id),
            "funk_description": row.get("description", ""),
            "assessment": row.get("assessment", ""),
            "n_genes": int(row["n_genes"]),
            "n_found": int(row["n_found"]),
            "preservation": round(float(row["preservation"]), 4),
            "funk_annotation": row.get("funk_annotation", ""),
            "funk_confidence": row.get("funk_confidence", ""),
            "dominant_brieflow_label": dom_bl_label,
            "dominant_count": int(row["dominant_count"]) if row.get("dominant_count") else 0,
            "dominant_annotation": row.get("dominant_annotation", ""),
            "dominant_confidence": row.get("dominant_confidence", ""),
            "n_brieflow_clusters": int(row["n_brieflow_clusters"]) if row.get("n_brieflow_clusters") else 0,
            "distribution": row.get("distribution", ""),
        })

    rows.sort(key=lambda r: (r["cell_class"], -r["preservation"]))
    fieldnames = [
        "funk_label", "cell_class", "funk_cluster", "funk_description",
        "assessment", "n_genes", "n_found", "preservation",
        "funk_annotation", "funk_confidence", "dominant_brieflow_label",
        "dominant_count", "dominant_annotation", "dominant_confidence",
        "n_brieflow_clusters", "distribution",
    ]
    write_csv(rows, OUTPUT_DIR / "S4_funk_cluster_retention.csv", fieldnames)


def generate_s5() -> None:
    """S5: Shuffled control comparison — confidence distributions."""
    print("S5: Shuffled control summary")

    # Count confidence levels for each dataset
    datasets = {
        "Brieflow Interphase (real)": INTERPHASE_SUMMARIES,
        "Brieflow Mitotic (real)": MITOTIC_SUMMARIES,
        "Brieflow Interphase (shuffled)": SHUFFLED_SUMMARIES,
    }

    summary_rows = []
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
            summary_rows.append({
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
    write_csv(summary_rows, OUTPUT_DIR / "S5_shuffled_control_summary.csv", fieldnames)


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}\n")

    generate_s1()
    generate_s2()
    generate_s3()
    generate_s4()
    generate_s5()

    print(f"\nDone. All tables written to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
