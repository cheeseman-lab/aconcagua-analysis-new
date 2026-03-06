"""
Funk Cluster Retention Analysis

For each functional cluster defined in Funk et al. 2022, assess how well
brieflow's clustering retains the published gene groupings.

Outputs:
    results/cluster/retention/funk_cluster_retention.tsv   — per-cluster metrics
    results/cluster/retention/funk_gene_retention.tsv      — per-gene mapping

Usage:
    python cluster_retention.py
"""

from collections import Counter
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BENCHMARKS_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BENCHMARKS_DIR.parent
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
EXTERNAL_DIR = BENCHMARKS_DIR / "external"
CLUSTER_DIR = (
    ANALYSIS_DIR / "brieflow_output" / "cluster" / "DAPI_TUBULIN_GH2AX_PHALLOIDIN"
)
FUNK_DIR = BENCHMARKS_DIR / "external" / "results" / "cluster"
TSV_DIR = BENCHMARKS_DIR / "results" / "cluster" / "retention"

FUNK_PUBLISHED_CSV = EXTERNAL_DIR / "1-s2.0-S0092867422013599-mmc4.csv"

IDEAL_CONFIGS = {
    "Interphase": {"brieflow_k": 12},
    "Mitotic": {"brieflow_k": 5},
}

FUNK_MOZZARELLM_CONFIGS = {
    "Interphase": "Funk_Interphase_k10",
    "Mitotic": "Funk_Mitotic_k9",
}

# ---------------------------------------------------------------------------
# Cluster definitions from Funk et al. 2022
# ---------------------------------------------------------------------------
FUNK_PAPER_CLUSTERS = {
    "Interphase": [
        # --- Translation machinery ---
        (66, "40S ribosome subunits"),
        (23, "60S ribosome subunits"),
        (14, "tRNA ligases & eIF2 translation initiation"),
        (136, "40S ribosome biogenesis"),
        (15, "60S ribosome biogenesis"),
        (21, "Nucleolar proteins / RNA helicases"),
        (112, "Related nucleolar factors (I)"),
        (203, "Related nucleolar factors (II)"),
        (216, "Related nucleolar factors (III)"),
        # --- Transcription ---
        (192, "TFIID complex"),
        (60, "Mediator / general transcription factors / mRNA export"),
        (8, "Mediator / general transcription factors (II)"),
        (199, "RNA polymerase (I)"),
        (45, "RNA polymerase (II)"),
        (155, "RNA polymerase (III)"),
        (121, "MYC/MAX transcription regulators"),
        (146, "Transcriptional regulation (I)"),
        (39, "Transcriptional regulation (II)"),
        # --- RNA processing ---
        (52, "RNA processing / splicing (I)"),
        (138, "RNA processing / splicing (II)"),
        (110, "RNA processing / splicing (III)"),
        (215, "RNA processing / splicing (IV)"),
        (157, "RNA processing / splicing (V)"),
        (145, "RNA processing / splicing (VI)"),
        (197, "m6A modification (METTL3/METTL14/HNRNPD)"),
        # --- Integrator / signaling ---
        (37, "Integrator complex / mTOR / ER-Golgi"),
        (217, "Integrator subunits (related to cluster 37)"),
        (149, "KRAS/BRAF signaling & mitochondria"),
        # --- Protein degradation ---
        (167, "20S proteasome core (AKIRIN2)"),
        (106, "Protein degradation / ubiquitin-proteasome (I)"),
        (213, "Protein degradation / ubiquitin-proteasome (II)"),
        (200, "Protein degradation / ubiquitin-proteasome (III)"),
        # --- DNA replication & damage ---
        (26, "DNA replication & damage response (I)"),
        (13, "DNA replication & damage response (II)"),
        (195, "DNA replication & damage response (III)"),
        (179, "DNA replication & damage response (IV)"),
        (3, "DNA replication & damage response (V)"),
        # --- Cell cycle ---
        (148, "Cytokinesis"),
        (218, "Cell cycle control (I)"),
        (95, "Cell cycle control (II)"),
        (46, "Cell cycle control (III)"),
        # --- Cytoskeleton & adhesion ---
        (29, "Cell adhesion / integrins"),
        (184, "Actin cytoskeleton & adhesion"),
        # --- Nuclear transport ---
        (104, "Nuclear transport / pore"),
        # --- Vesicle trafficking ---
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
# Data loading
# ---------------------------------------------------------------------------
def load_funk_published():
    df = pd.read_csv(FUNK_PUBLISHED_CSV)
    df = df.rename(columns={"\ufeffGene symbol": "Gene symbol"})
    if "Gene symbol" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Gene symbol"})
    return df


def load_brieflow_clustered(cell_class):
    k = IDEAL_CONFIGS[cell_class]["brieflow_k"]
    path = CLUSTER_DIR / cell_class / str(k) / "phate_leiden_clustering.tsv"
    return pd.read_csv(path, sep="\t")


def load_mozzarellm_summaries(cell_class):
    """Load mozzarellm summaries — returns full dataframe."""
    k = IDEAL_CONFIGS[cell_class]["brieflow_k"]
    path = (
        CLUSTER_DIR / cell_class / str(k) / "mozzarellm"
        / "claude-opus-4-6_summaries.tsv"
    )
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def load_funk_mozzarellm_summaries(cell_class):
    """Load Funk pipeline mozzarellm summaries."""
    path = (
        FUNK_DIR / FUNK_MOZZARELLM_CONFIGS[cell_class]
        / "mozzarellm" / "claude-opus-4-6_summaries.tsv"
    )
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
def analyze_cluster(funk_cluster_id, cell_class, funk_df, brieflow_df, mozzarellm_df, funk_mozzarellm_df):
    """Analyze how one Funk cluster maps onto brieflow clustering."""
    col = "Interphase cluster" if cell_class == "Interphase" else "Mitotic cluster"

    mask = pd.to_numeric(funk_df[col], errors="coerce") == funk_cluster_id
    genes = funk_df.loc[mask, "Gene symbol"].dropna().tolist()

    bl_lookup = dict(zip(brieflow_df["gene_symbol_0"], brieflow_df["cluster"]))

    # Build per-gene rows
    gene_rows = []
    for gene in sorted(genes):
        bl_cl = bl_lookup.get(gene)
        gene_rows.append({
            "gene": gene,
            "funk_cluster": funk_cluster_id,
            "brieflow_cluster": bl_cl,
        })

    found = [r["brieflow_cluster"] for r in gene_rows if r["brieflow_cluster"] is not None]
    n_found = len(found)
    distribution = Counter(found)

    if n_found == 0:
        return None, gene_rows

    dominant_cluster, dominant_count = distribution.most_common(1)[0]
    preservation = dominant_count / n_found

    # Get mozzarellm annotation for dominant cluster
    dom_ann = ""
    dom_conf = ""
    if not mozzarellm_df.empty:
        match = mozzarellm_df[mozzarellm_df["cluster_id"] == int(dominant_cluster)]
        if not match.empty:
            dom_ann = match.iloc[0].get("dominant_process", "")
            dom_conf = match.iloc[0].get("pathway_confidence", "")

    # Get Funk mozzarellm annotation for this Funk cluster
    funk_ann = ""
    funk_conf = ""
    if not funk_mozzarellm_df.empty:
        match = funk_mozzarellm_df[funk_mozzarellm_df["cluster_id"] == funk_cluster_id]
        if not match.empty:
            funk_ann = match.iloc[0].get("dominant_process", "")
            funk_conf = match.iloc[0].get("pathway_confidence", "")

    # Build distribution string: "42(28), 7(1)" format
    dist_parts = []
    for cl, cnt in distribution.most_common():
        dist_parts.append(f"{int(cl)}({cnt})")
    dist_str = ", ".join(dist_parts)

    # Annotate all non-dominant clusters
    other_annotations = []
    for cl, cnt in distribution.most_common():
        if cl == dominant_cluster:
            continue
        ann = ""
        if not mozzarellm_df.empty:
            m = mozzarellm_df[mozzarellm_df["cluster_id"] == int(cl)]
            if not m.empty:
                ann = m.iloc[0].get("dominant_process", "")
        other_annotations.append(f"BL{int(cl)}({cnt}): {ann}" if ann else f"BL{int(cl)}({cnt})")

    # Mark dominant in gene rows
    for row in gene_rows:
        row["in_dominant"] = (row["brieflow_cluster"] == dominant_cluster) if row["brieflow_cluster"] is not None else False
        # Add the mozzarellm annotation for this gene's brieflow cluster
        row["brieflow_annotation"] = ""
        if row["brieflow_cluster"] is not None and not mozzarellm_df.empty:
            m = mozzarellm_df[mozzarellm_df["cluster_id"] == int(row["brieflow_cluster"])]
            if not m.empty:
                row["brieflow_annotation"] = m.iloc[0].get("dominant_process", "")
        row["funk_annotation"] = funk_ann

    summary = {
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
        "distribution": dist_str,
        "other_clusters": "; ".join(other_annotations) if other_annotations else "",
    }

    return summary, gene_rows


def assess(preservation):
    if preservation >= 0.8:
        return "well"
    elif preservation >= 0.5:
        return "partially"
    else:
        return "poorly"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("FUNK CLUSTER RETENTION ANALYSIS")
    print("=" * 70)

    # Load data
    print("\nLoading data...")
    funk_df = load_funk_published()
    print(f"  Funk published: {len(funk_df)} genes")

    brieflow_data = {}
    mozzarellm_data = {}
    funk_mozzarellm_data = {}
    for cell_class in ["Interphase", "Mitotic"]:
        bl = load_brieflow_clustered(cell_class)
        brieflow_data[cell_class] = bl
        k = IDEAL_CONFIGS[cell_class]["brieflow_k"]
        print(f"  Brieflow {cell_class} k={k}: {len(bl)} genes")

        mozz = load_mozzarellm_summaries(cell_class)
        mozzarellm_data[cell_class] = mozz
        print(f"  Brieflow MozzareLLM annotations: {len(mozz)} clusters")

        funk_mozz = load_funk_mozzarellm_summaries(cell_class)
        funk_mozzarellm_data[cell_class] = funk_mozz
        print(f"  Funk MozzareLLM annotations: {len(funk_mozz)} clusters")

    # Analyze each cluster
    print("\nAnalyzing clusters...")
    cluster_rows = []
    all_gene_rows = []

    for cell_class, clusters in FUNK_PAPER_CLUSTERS.items():
        for funk_id, desc in clusters:
            label = f"M{funk_id}" if cell_class == "Mitotic" else f"Cluster {funk_id}"

            summary, gene_rows = analyze_cluster(
                funk_id, cell_class,
                funk_df, brieflow_data[cell_class],
                mozzarellm_data[cell_class],
                funk_mozzarellm_data[cell_class],
            )

            # Gene-level rows
            for gr in gene_rows:
                all_gene_rows.append({
                    "funk_cluster": funk_id,
                    "funk_description": desc,
                    "cell_class": cell_class,
                    "gene": gr["gene"],
                    "brieflow_cluster": int(gr["brieflow_cluster"]) if gr["brieflow_cluster"] is not None else None,
                    "in_dominant": gr.get("in_dominant", False),
                    "brieflow_annotation": gr.get("brieflow_annotation", ""),
                    "funk_annotation": gr.get("funk_annotation", ""),
                })

            if summary is None:
                cluster_rows.append({
                    "funk_cluster": funk_id,
                    "cell_class": cell_class,
                    "label": label,
                    "description": desc,
                    "n_genes": len(gene_rows),
                    "n_found": 0,
                    "preservation": 0.0,
                    "funk_annotation": "",
                    "funk_confidence": "",
                    "assessment": "no genes found",
                    "dominant_brieflow_cluster": None,
                    "dominant_count": 0,
                    "dominant_annotation": "",
                    "dominant_confidence": "",
                    "n_brieflow_clusters": 0,
                    "distribution": "",
                    "other_clusters": "",
                })
            else:
                cluster_rows.append({
                    "funk_cluster": funk_id,
                    "cell_class": cell_class,
                    "label": label,
                    "description": desc,
                    "assessment": assess(summary["preservation"]),
                    **summary,
                })

            pres_str = f"{summary['preservation']:.0%}" if summary else "N/A"
            n_bl = summary["n_brieflow_clusters"] if summary else 0
            n_genes = summary["n_genes"] if summary else len(gene_rows)
            n_found = summary["n_found"] if summary else 0
            print(f"  {label:>12s} | {n_genes:>3d} genes | "
                  f"{n_found:>3d} found | pres={pres_str:>4s} | "
                  f"{n_bl} BL clusters | {desc}")

    # Build DataFrames
    cluster_df = pd.DataFrame(cluster_rows)
    gene_df = pd.DataFrame(all_gene_rows)

    # Save TSVs
    TSV_DIR.mkdir(parents=True, exist_ok=True)

    cluster_tsv = TSV_DIR / "funk_cluster_retention.tsv"
    cluster_df.to_csv(cluster_tsv, sep="\t", index=False)
    print(f"\nSaved: {cluster_tsv}")

    gene_tsv = TSV_DIR / "funk_gene_retention.tsv"
    gene_df.to_csv(gene_tsv, sep="\t", index=False)
    print(f"Saved: {gene_tsv}")

    # Summary
    found = cluster_df[cluster_df["n_found"] > 0]
    mean_p = found["preservation"].mean()
    print(f"\nMean preservation: {mean_p:.0%} across {len(found)} clusters")


if __name__ == "__main__":
    main()
