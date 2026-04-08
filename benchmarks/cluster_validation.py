"""
Cluster Validation: Brieflow Preservation of Funk et al. 2022 Highlighted Clusters

Tracks all 57 individual Funk clusters highlighted in the paper (text discussion
and heatmap pathway groups) through the Brieflow pipeline. For each Funk cluster:

  1. Finds the dominant Brieflow cluster its genes land in
  2. Computes preservation (fraction of genes retained)
  3. Looks up MozzareLLM confidence on both the Funk and receiving Brieflow cluster

For fragmented clusters (<50% preservation), analyzes all receiving Brieflow
clusters to show that gene redistribution typically moves genes into more
specific, higher-confidence modules.

Outputs:
    results/cluster/validation/cluster_tracking.tsv          per-cluster metrics
    results/cluster/validation/redistribution_detail.tsv     redistribution flow detail
    results/cluster/validation/cluster_preservation.png      Fig A: all 57 clusters
    results/cluster/validation/cluster_preservation_text.png Fig B: 15 text-discussed clusters
    results/cluster/validation/redistribution_sankey.png     Fig C: Sankey diagram

Usage:
    python cluster_validation.py
"""

import argparse
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from collections import Counter
from pathlib import Path

from plot_style import setup_plot_style, save_figure

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
FUNK_CLUSTER_DIR = EXTERNAL_DIR / "results" / "cluster"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "cluster" / "validation"

# Funk supplementary CSV with published cluster assignments
FUNK_PUBLISHED_CSV = EXTERNAL_DIR / "1-s2.0-S0092867422013599-mmc4.csv"

# Ideal resolutions
IDEAL_CONFIGS = {
    "Interphase": {"brieflow_k": 12, "funk_k": 10},
    "Mitotic": {"brieflow_k": 5, "funk_k": 9},
}

# ---------------------------------------------------------------------------
# All Funk clusters highlighted in the paper (text + heatmap groups)
# ---------------------------------------------------------------------------

# Clusters discussed in the paper text (for filtering)
FUNK_TEXT_CLUSTERS = {66, 23, 14, 136, 15, 21, 167, 197, 37, 149, 121, 148, 104, 29}
FUNK_TEXT_CLUSTERS_MITOTIC = {109}

FUNK_HIGHLIGHTED_CLUSTERS = {
    "Interphase": {
        # Discussed in text
        66: "40S ribosome subunits",
        23: "60S ribosome subunits",
        14: "tRNA ligases & eIF2 initiation",
        136: "40S ribosome biogenesis",
        15: "60S ribosome biogenesis",
        21: "Nucleolar proteins / RNA helicases",
        167: "20S proteasome core",
        197: "m6A modification (METTL3/14)",
        37: "Integrator complex / mTOR / ER-Golgi",
        149: "KRAS/BRAF signaling & mitochondria",
        121: "MYC/MAX transcription regulators",
        148: "Cytokinesis factors",
        104: "Nuclear transport / pore",
        29: "Cell adhesion / integrins",
        # From heatmap pathway groups
        112: "Translation factors",
        216: "Translation regulation",
        203: "Translation-associated",
        155: "Transcriptional regulators",
        8: "Transcription co-factors",
        60: "Transcription / RNA pol",
        199: "Transcription-associated",
        192: "Transcription regulation",
        45: "Transcription factors",
        52: "Splicing factors",
        138: "RNA modification",
        110: "Spliceosome components",
        215: "RNA processing",
        157: "mRNA processing",
        145: "RNA metabolism",
        217: "Exosome / RNA decay",
        213: "Ubiquitin ligases",
        106: "Proteasome regulators",
        200: "Ubiquitin-proteasome",
        95: "Cell cycle regulators",
        46: "Cell cycle control",
        218: "Cell cycle-associated",
        26: "DNA replication",
        13: "DNA repair",
        195: "DNA damage response",
        179: "Replication-associated",
        3: "DNA damage / replication",
        212: "Chaperonin / tubulin folding",
        184: "Actin cytoskeleton",
        201: "Vesicle transport",
        54: "ER-Golgi trafficking",
        140: "Membrane trafficking",
    },
    "Mitotic": {
        109: "ZNF335 / spindle / gamma-tubulin",
        205: "Tubulin / spindle",
        214: "Augmin complex",
        11: "Chromosome alignment",
        17: "Translation (mitotic)",
        21: "Ribosome (mitotic)",
        0: "Translation factors (mitotic)",
        33: "Splicing (mitotic)",
        6: "mRNA processing (mitotic)",
        88: "Proteasome (mitotic)",
        34: "DNA replication (mitotic)",
    },
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_funk_published():
    """Load Funk's published cluster assignments from supplementary CSV."""
    df = pd.read_csv(FUNK_PUBLISHED_CSV)
    df = df.rename(columns={"\ufeffGene symbol": "Gene symbol"})
    if "Gene symbol" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Gene symbol"})
    return df


def load_mozzarellm_summaries(cell_class, pipeline):
    """Load MozzareLLM summaries. Returns {cluster_id: confidence}, {cluster_id: process}."""
    if pipeline == "brieflow":
        k = IDEAL_CONFIGS[cell_class]["brieflow_k"]
        path = (
            CLUSTER_DIR
            / cell_class
            / str(k)
            / "mozzarellm"
            / "claude-opus-4-6_summaries.tsv"
        )
    else:
        k = IDEAL_CONFIGS[cell_class]["funk_k"]
        path = (
            FUNK_CLUSTER_DIR
            / f"Funk_{cell_class}_k{k}"
            / "mozzarellm"
            / "claude-opus-4-6_summaries.tsv"
        )
    df = pd.read_csv(path, sep="\t")
    confidence = dict(zip(df["cluster_id"].astype(int), df["pathway_confidence"]))
    process = dict(zip(df["cluster_id"].astype(int), df["dominant_process"]))
    return confidence, process


def load_brieflow_clustered(cell_class):
    """Load Brieflow clustering at ideal k."""
    k = IDEAL_CONFIGS[cell_class]["brieflow_k"]
    path = CLUSTER_DIR / cell_class / str(k) / "phate_leiden_clustering.tsv"
    return pd.read_csv(path, sep="\t")


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------


def get_funk_genes(funk_published, cell_class, cluster_id):
    """Get gene list for a Funk published cluster."""
    if cell_class == "Interphase":
        genes = funk_published[
            funk_published["Interphase cluster"] == cluster_id
        ]["Gene symbol"].dropna().tolist()
    else:
        mitotic_cl = pd.to_numeric(
            funk_published["Mitotic cluster"], errors="coerce"
        )
        genes = funk_published[mitotic_cl == cluster_id][
            "Gene symbol"
        ].dropna().tolist()
    return genes


def track_cluster(genes, brieflow_df):
    """Track a set of genes into Brieflow clustering.

    Returns:
        preservation: fraction of found genes in dominant cluster
        dominant_cluster: Brieflow cluster receiving the most genes
        n_found: genes found in Brieflow
        distribution: Counter of Brieflow cluster assignments
    """
    matched = brieflow_df[brieflow_df["gene_symbol_0"].isin(genes)]
    n_found = len(matched)

    if n_found == 0:
        return 0.0, None, 0, Counter()

    distribution = Counter(matched["cluster"].values)
    dominant_cluster, dominant_count = distribution.most_common(1)[0]
    preservation = dominant_count / n_found

    return preservation, dominant_cluster, n_found, distribution


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def run_validation():
    """Run the full validation analysis."""
    print("=" * 80)
    print("CLUSTER VALIDATION: Brieflow Preservation of Funk Highlighted Clusters")
    print("=" * 80)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    print("\nLoading data...")
    funk_published = load_funk_published()
    print(f"  Funk published: {len(funk_published)} genes")

    brieflow_data = {}
    for cell_class in ["Interphase", "Mitotic"]:
        brieflow_data[cell_class] = load_brieflow_clustered(cell_class)
        print(f"  Brieflow {cell_class}: {len(brieflow_data[cell_class])} genes")

    # Load MozzareLLM confidence (keyed by cell class to avoid ID collisions)
    funk_conf = {}
    funk_proc = {}
    brieflow_conf = {}
    brieflow_proc = {}
    for cell_class in ["Interphase", "Mitotic"]:
        fc, fp = load_mozzarellm_summaries(cell_class, "funk")
        funk_conf[cell_class] = fc
        funk_proc[cell_class] = fp
        bc, bp = load_mozzarellm_summaries(cell_class, "brieflow")
        brieflow_conf[cell_class] = bc
        brieflow_proc[cell_class] = bp

    # Count total highlighted clusters
    n_total_clusters = sum(
        len(clusters) for clusters in FUNK_HIGHLIGHTED_CLUSTERS.values()
    )
    print(f"\n  Funk highlighted clusters: {n_total_clusters}")

    # -----------------------------------------------------------------------
    # Track each cluster through Brieflow
    # -----------------------------------------------------------------------
    print("\nTracking clusters through Brieflow...")
    rows = []

    for cell_class, clusters in FUNK_HIGHLIGHTED_CLUSTERS.items():
        for funk_cl_id, description in sorted(clusters.items()):
            genes = get_funk_genes(funk_published, cell_class, funk_cl_id)
            if not genes:
                continue

            preservation, bl_dom_cl, n_found, distribution = track_cluster(
                genes, brieflow_data[cell_class]
            )

            # Classify preservation
            if preservation >= 0.80:
                tier = "Well-preserved"
            elif preservation >= 0.50:
                tier = "Partially preserved"
            else:
                tier = "Fragmented"

            # MozzareLLM confidence (cell_class-aware to avoid ID collision)
            fk_confidence = funk_conf[cell_class].get(funk_cl_id, "N/A")
            fk_process = funk_proc[cell_class].get(funk_cl_id, "")
            bl_confidence = (
                brieflow_conf[cell_class].get(bl_dom_cl, "N/A")
                if bl_dom_cl is not None
                else "N/A"
            )
            bl_process = (
                brieflow_proc[cell_class].get(bl_dom_cl, "")
                if bl_dom_cl is not None
                else ""
            )

            # Track whether discussed in text vs heatmap only
            in_text = (
                (cell_class == "Interphase" and funk_cl_id in FUNK_TEXT_CLUSTERS)
                or (cell_class == "Mitotic" and funk_cl_id in FUNK_TEXT_CLUSTERS_MITOTIC)
            )

            rows.append({
                "cell_class": cell_class,
                "funk_cluster": funk_cl_id,
                "description": description,
                "highlighted_in": "text" if in_text else "heatmap",
                "n_genes": len(genes),
                "n_found_in_brieflow": n_found,
                "preservation": preservation,
                "tier": tier,
                "brieflow_dominant_cluster": bl_dom_cl,
                "n_brieflow_clusters": len(distribution),
                "funk_confidence": fk_confidence,
                "funk_process": fk_process,
                "brieflow_confidence": bl_confidence,
                "brieflow_process": bl_process,
            })

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_DIR / "cluster_tracking.tsv", sep="\t", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'cluster_tracking.tsv'}")

    # -----------------------------------------------------------------------
    # Summary statistics
    # -----------------------------------------------------------------------
    n_total = len(df)
    n_well = (df["tier"] == "Well-preserved").sum()
    n_partial = (df["tier"] == "Partially preserved").sum()
    n_frag = (df["tier"] == "Fragmented").sum()
    n_preserved = n_well + n_partial

    funk_hc = (df["funk_confidence"] == "High").sum()
    bl_hc = (df["brieflow_confidence"] == "High").sum()
    mean_well_pres = df.loc[df["tier"] == "Well-preserved", "preservation"].mean()

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"\n  Total Funk highlighted clusters: {n_total}")
    print(f"  Well-preserved (>=80%):          {n_well} ({100*n_well/n_total:.1f}%)")
    print(f"  Partially preserved (>=50%):     {n_partial} ({100*n_partial/n_total:.1f}%)")
    print(f"  Fragmented (<50%):               {n_frag} ({100*n_frag/n_total:.1f}%)")
    print(f"  Total preserved (>=50%):         {n_preserved} ({100*n_preserved/n_total:.1f}%)")
    if n_well > 0:
        print(f"  Mean preservation (well-pres):   {mean_well_pres:.1%}")
    print(f"\n  Funk high-confidence:            {funk_hc}/{n_total} ({100*funk_hc/n_total:.1f}%)")
    print(f"  Brieflow receiving high-conf:    {bl_hc}/{n_total} ({100*bl_hc/n_total:.1f}%)")

    # Detailed table
    print(f"\n{'CL':>4s} {'Class':<11s} {'N':>3s} {'Pres':>6s} {'Tier':<20s} "
          f"{'FK Conf':<8s} {'BL Cl':>5s} {'BL Conf':<8s} {'Description'}")
    print("-" * 110)
    for _, row in df.sort_values(
        ["cell_class", "preservation"], ascending=[True, False]
    ).iterrows():
        bl_cl_str = (
            str(int(row["brieflow_dominant_cluster"]))
            if pd.notna(row["brieflow_dominant_cluster"])
            else "—"
        )
        print(
            f"{row['funk_cluster']:>4d} {row['cell_class']:<11s} {row['n_genes']:>3d} "
            f"{row['preservation']:>6.0%} {row['tier']:<20s} "
            f"{row['funk_confidence']:<8s} {bl_cl_str:>5s} {row['brieflow_confidence']:<8s} "
            f"{row['description']}"
        )

    # -----------------------------------------------------------------------
    # Redistribution analysis (non-well-preserved text clusters)
    # -----------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("REDISTRIBUTION ANALYSIS (text-discussed clusters)")
    print("=" * 80)

    text_df = df[df["highlighted_in"] == "text"]
    redistributed = text_df[text_df["tier"] != "Well-preserved"]
    min_genes = 2

    redist_rows = []
    for _, row in redistributed.iterrows():
        cell_class = row["cell_class"]
        funk_cl_id = int(row["funk_cluster"])
        genes = get_funk_genes(funk_published, cell_class, funk_cl_id)

        matched = brieflow_data[cell_class][
            brieflow_data[cell_class]["gene_symbol_0"].isin(genes)
        ]
        dist = Counter(matched["cluster"].values)
        receivers = [
            (cl, cnt) for cl, cnt in dist.most_common() if cnt >= min_genes
        ]

        for bl_cl, n_genes_recv in receivers:
            redist_rows.append({
                "funk_cluster": funk_cl_id,
                "cell_class": cell_class,
                "funk_description": row["description"],
                "funk_confidence": row["funk_confidence"],
                "funk_n_genes": row["n_genes"],
                "brieflow_cluster": bl_cl,
                "n_genes_received": n_genes_recv,
                "brieflow_confidence": brieflow_conf[cell_class].get(bl_cl, "N/A"),
                "brieflow_process": brieflow_proc[cell_class].get(bl_cl, ""),
            })

    redist_df = pd.DataFrame(redist_rows)
    redist_df.to_csv(OUTPUT_DIR / "redistribution_detail.tsv", sep="\t", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'redistribution_detail.tsv'}")

    # Summary
    n_receivers = len(redist_df)
    n_recv_high = (redist_df["brieflow_confidence"] == "High").sum()
    n_recv_med = (redist_df["brieflow_confidence"] == "Medium").sum()
    n_recv_low = (redist_df["brieflow_confidence"] == "Low").sum()
    genes_in_high = redist_df.loc[
        redist_df["brieflow_confidence"] == "High", "n_genes_received"
    ].sum()
    genes_total = redist_df["n_genes_received"].sum()

    print(f"\n  Non-well-preserved text clusters: {len(redistributed)}")
    print(f"    Partially preserved: {(redistributed['tier'] == 'Partially preserved').sum()}")
    print(f"    Fragmented: {(redistributed['tier'] == 'Fragmented').sum()}")
    print(f"\n  Receiving Brieflow clusters (>={min_genes} genes): {n_receivers}")
    print(f"    High-conf:   {n_recv_high} ({100*n_recv_high/n_receivers:.0f}%)")
    print(f"    Medium-conf: {n_recv_med} ({100*n_recv_med/n_receivers:.0f}%)")
    print(f"    Low-conf:    {n_recv_low} ({100*n_recv_low/n_receivers:.0f}%)")
    print(f"\n  Genes in High-conf receivers: {genes_in_high}/{genes_total} "
          f"({100*genes_in_high/genes_total:.0f}%)")

    # -----------------------------------------------------------------------
    # Figure A: Preservation bars — all 57 clusters (supplemental)
    # -----------------------------------------------------------------------
    plot_preservation(df, OUTPUT_DIR / "cluster_preservation.png")

    # -----------------------------------------------------------------------
    # Figure B: Preservation bars — 15 text-discussed clusters
    # -----------------------------------------------------------------------
    plot_preservation(
        text_df, OUTPUT_DIR / "cluster_preservation_text.png",
        title="Brieflow Preservation of Funk Text-Discussed Clusters",
    )

    # -----------------------------------------------------------------------
    # Figure C: Redistribution Sankey — text clusters, non-well-preserved
    # -----------------------------------------------------------------------
    plot_fragmentation_sankey(
        redist_df, brieflow_conf, brieflow_proc,
        OUTPUT_DIR / "redistribution_sankey.png",
    )

    print("\n" + "=" * 80)
    print("VALIDATION COMPLETE")
    print("=" * 80)
    print(f"Outputs in: {OUTPUT_DIR}")


# ---------------------------------------------------------------------------
# Figure A: Preservation bar chart
# ---------------------------------------------------------------------------

TIER_COLORS = {
    "Well-preserved": "#2ca02c",
    "Partially preserved": "#ffbb78",
    "Fragmented": "#d62728",
}


def plot_preservation(df, output_path, title="Brieflow Preservation of Funk-Highlighted Clusters"):
    """Sorted horizontal bar chart of preservation for highlighted clusters."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    setup_plot_style()

    df_sorted = df.sort_values("preservation", ascending=True).reset_index(drop=True)
    n = len(df_sorted)

    fig, ax = plt.subplots(figsize=(14, max(6, n * 0.5)))

    y = np.arange(n)
    colors = [TIER_COLORS[t] for t in df_sorted["tier"]]

    ax.barh(y, df_sorted["preservation"], color=colors, alpha=0.85, height=0.85)

    # Threshold lines
    ax.axvline(0.50, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axvline(0.80, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)

    # Labels inside bars
    for i, (_, row) in enumerate(df_sorted.iterrows()):
        pres = row["preservation"]
        desc = row["description"]
        if row.get("cell_class") == "Mitotic" and "(mitotic)" not in desc.lower():
            desc = f"{desc} (mitotic)"
        label = f"{desc} ({row['n_genes']})"
        pct = f"{pres:.0%}"

        if pres >= 0.35:
            # Label inside bar, left-aligned
            ax.text(
                0.01, i, f"  {label}", va="center", ha="left",
                fontsize=11, color="white", fontweight="bold",
            )
            # Percentage at right edge of bar
            ax.text(
                pres - 0.01, i, pct, va="center", ha="right",
                fontsize=11, color="white", fontweight="bold",
            )
        else:
            # Short bars: label outside
            ax.text(
                pres + 0.01, i, f"{label}  {pct}", va="center", ha="left",
                fontsize=11, color="black",
            )

    ax.set_yticks([])
    ax.set_xlabel("Preservation (fraction in dominant Brieflow cluster)", fontsize=13)
    ax.set_xlim(0, 1.15)
    ax.set_title(title, fontweight="bold", fontsize=18, fontfamily="Arial")

    legend_elements = [
        Patch(facecolor=TIER_COLORS["Well-preserved"], label="Well-preserved (≥80%)"),
        Patch(facecolor=TIER_COLORS["Partially preserved"], label="Partially preserved (≥50%)"),
        Patch(facecolor=TIER_COLORS["Fragmented"], label="Fragmented (<50%)"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=12)

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"  Saved: {output_path}")


# ---------------------------------------------------------------------------
# Figure B: Fragmentation Sankey
# ---------------------------------------------------------------------------

CONF_NODE_COLORS = {
    "High": "rgb(30,140,30)",
    "Medium": "rgb(230,160,40)",
    "Low": "rgb(200,40,40)",
}
CONF_LINK_COLORS = {
    "High": "rgba(30,140,30,0.35)",
    "Medium": "rgba(230,160,40,0.35)",
    "Low": "rgba(200,40,40,0.35)",
}


def _shorten_process(proc, maxlen=60):
    """Lightly truncate process name, keeping as much as possible."""
    if len(proc) > maxlen:
        proc = proc[: maxlen - 3] + "..."
    return proc


def plot_fragmentation_sankey(frag_df, brieflow_conf, brieflow_proc, output_path,
                              min_genes=3):
    """Sankey diagram: fragmented Funk clusters → Brieflow receiving clusters.

    Args:
        frag_df: DataFrame from fragmentation analysis
        brieflow_conf: {cell_class: {cluster_id: confidence}} — cell_class-keyed
        brieflow_proc: {cell_class: {cluster_id: process}} — cell_class-keyed
        output_path: where to save the PNG
        min_genes: minimum genes per flow to include
    """
    df = frag_df[frag_df["n_genes_received"] >= min_genes].copy()

    # Resolve confidence using cell_class-aware lookup
    df["bl_conf"] = df.apply(
        lambda r: brieflow_conf[r["cell_class"]].get(
            int(r["brieflow_cluster"]), "N/A"
        ),
        axis=1,
    )
    df["bl_proc"] = df.apply(
        lambda r: brieflow_proc[r["cell_class"]].get(
            int(r["brieflow_cluster"]), "?"
        ),
        axis=1,
    )

    # Funk cluster info
    funk_info = df.groupby("funk_cluster").first()[
        ["cell_class", "funk_description", "funk_confidence", "funk_n_genes"]
    ].to_dict("index")

    funk_ids = sorted(df["funk_cluster"].unique())

    # Right-side nodes keyed by (cell_class, cluster) to avoid ID collisions
    bl_keys = sorted(
        df[["cell_class", "brieflow_cluster"]].drop_duplicates().itertuples(
            index=False
        )
    )

    # Build node lists
    nodes, node_colors = [], []

    funk_idx = {}
    for i, fk in enumerate(funk_ids):
        info = funk_info[fk]
        desc = info["funk_description"]
        if info["cell_class"] == "Mitotic" and "(mitotic)" not in desc.lower():
            desc = f"{desc} (mitotic)"
        nodes.append(desc)
        node_colors.append(CONF_NODE_COLORS[info["funk_confidence"]])
        funk_idx[fk] = i

    bl_idx = {}
    offset = len(funk_ids)
    for j, (cc, bl_cl) in enumerate(bl_keys):
        proc = _shorten_process(
            brieflow_proc[cc].get(int(bl_cl), "?")
        )
        conf = brieflow_conf[cc].get(int(bl_cl), "N/A")
        nodes.append(proc)
        node_colors.append(
            CONF_NODE_COLORS.get(conf, "rgb(150,150,150)")
        )
        bl_idx[(cc, bl_cl)] = offset + j

    # Build links
    sources, targets, values, link_colors = [], [], [], []
    for _, row in df.iterrows():
        sources.append(funk_idx[row["funk_cluster"]])
        targets.append(bl_idx[(row["cell_class"], row["brieflow_cluster"])])
        values.append(row["n_genes_received"])
        link_colors.append(
            CONF_LINK_COLORS.get(row["bl_conf"], "rgba(150,150,150,0.2)")
        )

    fig = go.Figure(
        data=[
            go.Sankey(
                arrangement="snap",
                node=dict(
                    pad=14,
                    thickness=20,
                    line=dict(color="black", width=0.5),
                    label=nodes,
                    color=node_colors,
                ),
                link=dict(
                    source=sources,
                    target=targets,
                    value=values,
                    color=link_colors,
                ),
            )
        ]
    )

    fig.update_layout(
        title=dict(
            text=(
                "<b>Fragmented Funk Clusters → Brieflow Receiving Clusters</b><br>"
                "<sup>Node color: "
                "<b style='color:rgb(30,140,30)'>■ High</b>  "
                "<b style='color:rgb(230,160,40)'>■ Medium</b>  "
                "<b style='color:rgb(200,40,40)'>■ Low</b> confidence</sup>"
            ),
            font_size=15,
        ),
        font=dict(size=11, family="Arial"),
        width=1800,
        height=950,
        margin=dict(l=5, r=5, t=80, b=10),
    )

    fig.write_image(str(output_path), scale=4)
    pdf_path = Path(output_path).with_suffix(".pdf")
    fig.write_image(str(pdf_path))
    print(f"  Saved: {output_path} + {pdf_path.name}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description="Cluster validation: Brieflow preservation of Funk clusters")
    parser.add_argument("--plots-only", action="store_true",
                        help="Regenerate figures from cached TSVs without recomputing")
    args = parser.parse_args()

    if args.plots_only:
        print("Loading cached TSVs...")
        df = pd.read_csv(OUTPUT_DIR / "cluster_tracking.tsv", sep="\t")
        redist_df = pd.read_csv(OUTPUT_DIR / "redistribution_detail.tsv", sep="\t")

        # Reconstruct brieflow_conf and brieflow_proc from redistribution_detail
        brieflow_conf = {}
        brieflow_proc = {}
        for _, row in redist_df.iterrows():
            cc = row["cell_class"]
            bl_cl = int(row["brieflow_cluster"])
            if cc not in brieflow_conf:
                brieflow_conf[cc] = {}
                brieflow_proc[cc] = {}
            brieflow_conf[cc][bl_cl] = row["brieflow_confidence"]
            brieflow_proc[cc][bl_cl] = row["brieflow_process"]

        text_df = df[df["highlighted_in"] == "text"]

        print("Regenerating figures...")
        plot_preservation(df, OUTPUT_DIR / "cluster_preservation.png")
        plot_preservation(
            text_df, OUTPUT_DIR / "cluster_preservation_text.png",
            title="Brieflow Preservation of Funk Text-Discussed Clusters",
        )
        plot_fragmentation_sankey(
            redist_df, brieflow_conf, brieflow_proc,
            OUTPUT_DIR / "redistribution_sankey.png",
        )
        print(f"\nDone. Outputs in: {OUTPUT_DIR}")
        return

    run_validation()


if __name__ == "__main__":
    main()
