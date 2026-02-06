"""
Validation Gene Cluster Tracking: Brieflow vs Funk et al. 2022

Two levels of analysis:

  Level 1 — Individual clusters: tracks co-clustering preservation for
            specific clusters highlighted in the paper text (Funk and Brieflow).

  Level 2 — Pathway groups: tracks broader functional categories that span
            multiple Funk clusters (e.g. "translation" = clusters 66,136,21,...).
            Measures how many pipeline clusters the group scatters into
            (normalized entropy) and mean within-group preservation.

Outputs:
    results/cluster/validation/validation_tracking.tsv        (level 1)
    results/cluster/validation/validation_gene_detail.tsv     (level 1)
    results/cluster/validation/validation_preservation.png    (level 1)
    results/cluster/validation/group_tracking.tsv             (level 2)
    results/cluster/validation/group_preservation.png         (level 2)

Usage:
    python cluster_validation.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
from pathlib import Path

from plot_style import setup_plot_style, COLORS, save_figure

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
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "cluster" / "validation"

# Funk supplementary CSV with published cluster assignments
FUNK_PUBLISHED_CSV = EXTERNAL_DIR / "1-s2.0-S0092867422013599-mmc4.csv"

# Ideal resolutions (matching cluster_benchmarks.py)
IDEAL_CONFIGS = {
    "Interphase": {"brieflow_k": 12, "funk_k": 10},
    "Mitotic": {"brieflow_k": 5, "funk_k": 9},
}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_funk_published():
    """Load Funk's published cluster assignments from supplementary CSV."""
    df = pd.read_csv(FUNK_PUBLISHED_CSV)
    # Standardize column access
    df = df.rename(
        columns={
            "\ufeffGene symbol": "Gene symbol",  # handle BOM
        }
    )
    # Use first column if 'Gene symbol' still not found
    if "Gene symbol" not in df.columns:
        df = df.rename(columns={df.columns[0]: "Gene symbol"})
    return df


def load_funk_reclustered(cell_class):
    """Load Funk re-clustered at ideal k."""
    k = IDEAL_CONFIGS[cell_class]["funk_k"]
    path = (
        EXTERNAL_DIR
        / "results"
        / "cluster"
        / f"Funk_{cell_class}_k{k}"
        / "phate_leiden_clustering.tsv"
    )
    df = pd.read_csv(path, sep="\t")
    return df


def load_brieflow_clustered(cell_class):
    """Load Brieflow clustering at ideal k."""
    k = IDEAL_CONFIGS[cell_class]["brieflow_k"]
    path = CLUSTER_DIR / cell_class / str(k) / "phate_leiden_clustering.tsv"
    df = pd.read_csv(path, sep="\t")
    return df


# ---------------------------------------------------------------------------
# Define validation gene groups
# ---------------------------------------------------------------------------


def build_validation_groups(funk_published_df, brieflow_interphase_df):
    """Build validation gene groups from published clusters and Brieflow clusters.

    Returns list of dicts with keys: name, source, cell_class, genes, description.
    """
    groups = []

    # --- Funk Interphase groups (from published supplementary Table S3) ---
    funk_interphase_groups = {
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
    }

    for cluster_id, description in funk_interphase_groups.items():
        genes = (
            funk_published_df[funk_published_df["Interphase cluster"] == cluster_id][
                "Gene symbol"
            ]
            .dropna()
            .tolist()
        )
        if genes:
            groups.append(
                {
                    "name": f"Funk_I_{cluster_id}",
                    "source": "Funk published",
                    "cell_class": "Interphase",
                    "original_cluster": cluster_id,
                    "genes": genes,
                    "description": description,
                }
            )

    # --- Funk Mitotic group: M109 / ZNF335 cluster ---
    m109_genes = (
        funk_published_df[
            funk_published_df["Mitotic cluster"].astype(str).str.strip() == "109"
        ]["Gene symbol"]
        .dropna()
        .tolist()
    )
    if m109_genes:
        groups.append(
            {
                "name": "Funk_M_109",
                "source": "Funk published",
                "cell_class": "Mitotic",
                "original_cluster": 109,
                "genes": m109_genes,
                "description": "ZNF335 / spindle / gamma-tubulin",
            }
        )

    # --- Brieflow Interphase groups (from Brieflow k=12 clustering) ---
    brieflow_interphase_groups = {
        5: "Ribosome biogenesis (C18orf21/NEPRO)",
        "10+61+87": "Mitochondrial clusters (respiratory chain)",
        60: "mRNA processing / vesicular trafficking (CCR4-NOT / exocyst)",
        95: "MYC-dependent transcription (ZMYND8/KAT2A/SETD2)",
    }

    for cluster_key, description in brieflow_interphase_groups.items():
        if isinstance(cluster_key, str) and "+" in cluster_key:
            # Combined clusters
            cluster_ids = [int(x) for x in cluster_key.split("+")]
            genes = (
                brieflow_interphase_df[
                    brieflow_interphase_df["cluster"].isin(cluster_ids)
                ]["gene_symbol_0"]
                .dropna()
                .tolist()
            )
            name = f"Brieflow_I_{cluster_key}"
        else:
            genes = (
                brieflow_interphase_df[
                    brieflow_interphase_df["cluster"] == cluster_key
                ]["gene_symbol_0"]
                .dropna()
                .tolist()
            )
            name = f"Brieflow_I_{cluster_key}"

        # Filter out control genes
        genes = [g for g in genes if not g.startswith("nontargeting")]

        if genes:
            groups.append(
                {
                    "name": name,
                    "source": "Brieflow",
                    "cell_class": "Interphase",
                    "original_cluster": cluster_key,
                    "genes": genes,
                    "description": description,
                }
            )

    return groups


# ---------------------------------------------------------------------------
# Level 2: Pathway group definitions
# ---------------------------------------------------------------------------

# Funk's broader functional categories (each maps to multiple published clusters)
FUNK_INTERPHASE_PATHWAY_GROUPS = {
    "I_translation": {
        "clusters": [66, 136, 21, 112, 216, 23, 15, 14, 203],
        "description": "Ribosome subunits, biogenesis, tRNA ligases, translation factors",
    },
    "I_transcription": {
        "clusters": [155, 8, 60, 199, 192, 45],
        "description": "Transcriptional regulators and co-factors",
    },
    "I_rna_processing": {
        "clusters": [52, 138, 110, 215, 157, 145, 217],
        "description": "Splicing, RNA modification, exosome",
    },
    "I_protein_degradation": {
        "clusters": [213, 106, 167, 200],
        "description": "Proteasome, ubiquitin-proteasome system",
    },
    "I_cell_cycle": {
        "clusters": [148, 95, 46, 218],
        "description": "Cytokinesis, cell cycle control",
    },
    "I_dna_replication_damage": {
        "clusters": [26, 13, 195, 179, 3],
        "description": "DNA replication, repair, damage response",
    },
    "I_chaperone_tubulin": {
        "clusters": [212],
        "description": "Chaperonin / tubulin folding",
    },
    "I_actin_adhesion": {
        "clusters": [29, 184],
        "description": "Actin cytoskeleton, integrins, adhesion",
    },
    "I_vesicle_trafficking": {
        "clusters": [201, 54, 140],
        "description": "Vesicle transport, ER-Golgi",
    },
    "I_nuclear_pore": {
        "clusters": [104],
        "description": "Nuclear pore complex, nuclear transport",
    },
}

FUNK_MITOTIC_PATHWAY_GROUPS = {
    "M_spindle_assembly": {
        "clusters": [205, 109, 214, 11],
        "description": "Tubulin, spindle bipolarity, augmin, chromosome alignment",
    },
    "M_translation": {
        "clusters": [17, 21, 0],
        "description": "Translation (mitotic)",
    },
    "M_mrna_splicing": {
        "clusters": [33, 6],
        "description": "mRNA splicing",
    },
    "M_protein_degradation": {
        "clusters": [88],
        "description": "Proteasome",
    },
    "M_dna_replication": {
        "clusters": [34],
        "description": "DNA replication",
    },
}


def build_pathway_groups(funk_published_df):
    """Build pathway-level gene groups from Funk's functional categories.

    Returns list of dicts with keys: name, cell_class, clusters, genes, description.
    """
    groups = []

    for group_name, info in FUNK_INTERPHASE_PATHWAY_GROUPS.items():
        genes = (
            funk_published_df[
                funk_published_df["Interphase cluster"].isin(info["clusters"])
            ]["Gene symbol"]
            .dropna()
            .tolist()
        )
        if genes:
            groups.append(
                {
                    "name": group_name,
                    "cell_class": "Interphase",
                    "funk_clusters": info["clusters"],
                    "n_funk_clusters": len(info["clusters"]),
                    "genes": genes,
                    "description": info["description"],
                }
            )

    for group_name, info in FUNK_MITOTIC_PATHWAY_GROUPS.items():
        # Mitotic cluster column may have mixed types; coerce to numeric
        mitotic_cl = pd.to_numeric(
            funk_published_df["Mitotic cluster"], errors="coerce"
        )
        genes = (
            funk_published_df[mitotic_cl.isin(info["clusters"])]["Gene symbol"]
            .dropna()
            .tolist()
        )
        if genes:
            groups.append(
                {
                    "name": group_name,
                    "cell_class": "Mitotic",
                    "funk_clusters": info["clusters"],
                    "n_funk_clusters": len(info["clusters"]),
                    "genes": genes,
                    "description": info["description"],
                }
            )

    return groups


def compute_group_metrics(genes, cluster_assignments, gene_col="gene_symbol"):
    """Compute group-level cluster spread metrics.

    Returns:
        n_found: genes found in this clustering
        n_clusters: number of distinct clusters the genes land in
        normalized_entropy: Shannon entropy / log(n_found), 0=all in one cluster, 1=uniform
        dominant_preservation: fraction in the single largest cluster
        distribution: Counter of cluster assignments
    """
    df = cluster_assignments[cluster_assignments[gene_col].isin(genes)]
    n_found = len(df)

    if n_found == 0:
        return 0, 0, np.nan, 0.0, Counter()

    distribution = Counter(df["cluster"].values)
    n_clusters = len(distribution)

    # Normalized Shannon entropy
    counts = np.array(list(distribution.values()), dtype=float)
    probs = counts / counts.sum()
    entropy = -np.sum(probs * np.log(probs))
    max_entropy = np.log(n_found) if n_found > 1 else 1.0
    normalized_entropy = entropy / max_entropy if max_entropy > 0 else 0.0

    dominant_count = distribution.most_common(1)[0][1]
    dominant_preservation = dominant_count / n_found

    return n_found, n_clusters, normalized_entropy, dominant_preservation, distribution


# ---------------------------------------------------------------------------
# Co-clustering preservation (Level 1)
# ---------------------------------------------------------------------------


def compute_preservation(genes, cluster_assignments, gene_col="gene_symbol"):
    """Compute co-clustering preservation for a gene group.

    Returns:
        preservation: fraction of genes in the dominant cluster
        dominant_cluster: the cluster containing the most genes
        n_found: number of genes found in the clustering
        n_total: total genes in the group
        distribution: Counter of cluster assignments
    """
    df = cluster_assignments[cluster_assignments[gene_col].isin(genes)]
    n_found = len(df)
    n_total = len(genes)

    if n_found == 0:
        return 0.0, None, 0, n_total, Counter()

    distribution = Counter(df["cluster"].values)
    dominant_cluster, dominant_count = distribution.most_common(1)[0]
    preservation = dominant_count / n_found

    return preservation, dominant_cluster, n_found, n_total, distribution


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------


def run_validation():
    """Run the full validation analysis."""
    print("=" * 80)
    print("VALIDATION GENE CLUSTER TRACKING")
    print("=" * 80)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    print("\nLoading data...")
    funk_published = load_funk_published()
    print(f"  Funk published: {len(funk_published)} genes")

    # Load re-clustered data
    cluster_data = {}
    for cell_class in ["Interphase", "Mitotic"]:
        funk_re = load_funk_reclustered(cell_class)
        brieflow = load_brieflow_clustered(cell_class)
        cluster_data[cell_class] = {
            "funk_reclustered": funk_re,
            "brieflow": brieflow,
        }
        fk = IDEAL_CONFIGS[cell_class]["funk_k"]
        bk = IDEAL_CONFIGS[cell_class]["brieflow_k"]
        print(
            f"  {cell_class}: Funk k={fk} ({len(funk_re)} genes), Brieflow k={bk} ({len(brieflow)} genes)"
        )

    # Build validation groups
    brieflow_interphase = cluster_data["Interphase"]["brieflow"]
    groups = build_validation_groups(funk_published, brieflow_interphase)
    print(f"\n  Validation groups: {len(groups)}")

    # Compute preservation for each group
    print("\nComputing co-clustering preservation...")
    summary_rows = []
    detail_rows = []

    for group in groups:
        cell_class = group["cell_class"]
        genes = group["genes"]
        cd = cluster_data[cell_class]

        # Funk re-clustered preservation
        funk_pres, funk_dom, funk_found, _, funk_dist = compute_preservation(
            genes, cd["funk_reclustered"], gene_col="gene_symbol"
        )

        # Brieflow preservation
        brieflow_pres, brieflow_dom, brieflow_found, _, brieflow_dist = (
            compute_preservation(genes, cd["brieflow"], gene_col="gene_symbol_0")
        )

        summary_rows.append(
            {
                "group": group["name"],
                "source": group["source"],
                "cell_class": cell_class,
                "description": group["description"],
                "n_genes": len(genes),
                "funk_k": IDEAL_CONFIGS[cell_class]["funk_k"],
                "funk_found": funk_found,
                "funk_preservation": funk_pres,
                "funk_dominant_cluster": funk_dom,
                "funk_n_clusters": len(funk_dist),
                "brieflow_k": IDEAL_CONFIGS[cell_class]["brieflow_k"],
                "brieflow_found": brieflow_found,
                "brieflow_preservation": brieflow_pres,
                "brieflow_dominant_cluster": brieflow_dom,
                "brieflow_n_clusters": len(brieflow_dist),
            }
        )

        # Per-gene detail
        funk_lookup = dict(
            zip(
                cd["funk_reclustered"]["gene_symbol"],
                cd["funk_reclustered"]["cluster"],
            )
        )
        brieflow_lookup = dict(
            zip(
                cd["brieflow"]["gene_symbol_0"],
                cd["brieflow"]["cluster"],
            )
        )

        # Funk published lookup
        if cell_class == "Interphase":
            pub_col = "Interphase cluster"
        else:
            pub_col = "Mitotic cluster"
        pub_lookup = dict(
            zip(
                funk_published["Gene symbol"],
                funk_published[pub_col],
            )
        )

        for gene in genes:
            detail_rows.append(
                {
                    "gene": gene,
                    "group": group["name"],
                    "description": group["description"],
                    "cell_class": cell_class,
                    "funk_published_cluster": pub_lookup.get(gene, None),
                    "funk_reclustered_cluster": funk_lookup.get(gene, None),
                    "brieflow_cluster": brieflow_lookup.get(gene, None),
                }
            )

    # Build DataFrames
    summary_df = pd.DataFrame(summary_rows)
    detail_df = pd.DataFrame(detail_rows)

    # Save tables
    summary_df.to_csv(OUTPUT_DIR / "validation_tracking.tsv", sep="\t", index=False)
    print(f"\n  Saved: {OUTPUT_DIR / 'validation_tracking.tsv'}")

    detail_df.to_csv(OUTPUT_DIR / "validation_gene_detail.tsv", sep="\t", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'validation_gene_detail.tsv'}")

    # Print summary
    print("\n" + "=" * 80)
    print("PRESERVATION SUMMARY")
    print("=" * 80)
    print(
        f"\n{'Group':<22s} {'Src':<8s} {'N':>4s}  "
        f"{'Funk Pres':>10s} {'Funk #Cl':>8s}  "
        f"{'BL Pres':>10s} {'BL #Cl':>8s}  {'Description'}"
    )
    print("-" * 110)
    for _, row in summary_df.iterrows():
        print(
            f"{row['group']:<22s} {row['source']:<8s} {row['n_genes']:>4d}  "
            f"{row['funk_preservation']:>10.1%} {row['funk_n_clusters']:>8d}  "
            f"{row['brieflow_preservation']:>10.1%} {row['brieflow_n_clusters']:>8d}  "
            f"{row['description']}"
        )

    # Plot level 1
    plot_preservation(summary_df, OUTPUT_DIR / "validation_preservation.png")

    # ===================================================================
    # Level 2: Pathway group analysis
    # ===================================================================
    print("\n" + "=" * 80)
    print("LEVEL 2: PATHWAY GROUP PRESERVATION")
    print("=" * 80)

    pathway_groups = build_pathway_groups(funk_published)
    print(f"\n  Pathway groups: {len(pathway_groups)}")

    group_rows = []
    for pg in pathway_groups:
        cell_class = pg["cell_class"]
        genes = pg["genes"]
        cd = cluster_data[cell_class]

        # Funk re-clustered
        funk_found, funk_ncl, funk_ent, funk_dom_pres, funk_dist = (
            compute_group_metrics(genes, cd["funk_reclustered"], gene_col="gene_symbol")
        )

        # Brieflow
        bl_found, bl_ncl, bl_ent, bl_dom_pres, bl_dist = compute_group_metrics(
            genes, cd["brieflow"], gene_col="gene_symbol_0"
        )

        group_rows.append(
            {
                "pathway_group": pg["name"],
                "cell_class": cell_class,
                "description": pg["description"],
                "n_genes": len(genes),
                "n_funk_published_clusters": pg["n_funk_clusters"],
                "funk_k": IDEAL_CONFIGS[cell_class]["funk_k"],
                "funk_found": funk_found,
                "funk_n_clusters": funk_ncl,
                "funk_normalized_entropy": funk_ent,
                "funk_dominant_preservation": funk_dom_pres,
                "brieflow_k": IDEAL_CONFIGS[cell_class]["brieflow_k"],
                "brieflow_found": bl_found,
                "brieflow_n_clusters": bl_ncl,
                "brieflow_normalized_entropy": bl_ent,
                "brieflow_dominant_preservation": bl_dom_pres,
            }
        )

    group_df = pd.DataFrame(group_rows)
    group_df.to_csv(OUTPUT_DIR / "group_tracking.tsv", sep="\t", index=False)
    print(f"  Saved: {OUTPUT_DIR / 'group_tracking.tsv'}")

    # Print summary
    print(
        f"\n{'Group':<25s} {'N':>5s} {'Funk#Pub':>8s}  "
        f"{'Funk#Cl':>7s} {'FunkEnt':>8s}  "
        f"{'BL#Cl':>7s} {'BL_Ent':>8s}  {'Description'}"
    )
    print("-" * 115)
    for _, row in group_df.iterrows():
        print(
            f"{row['pathway_group']:<25s} {row['n_genes']:>5d} {row['n_funk_published_clusters']:>8d}  "
            f"{row['funk_n_clusters']:>7d} {row['funk_normalized_entropy']:>8.3f}  "
            f"{row['brieflow_n_clusters']:>7d} {row['brieflow_normalized_entropy']:>8.3f}  "
            f"{row['description']}"
        )

    # Plot level 2
    plot_group_preservation(group_df, OUTPUT_DIR / "group_preservation.png")

    print("\n" + "=" * 80)
    print("VALIDATION COMPLETE")
    print("=" * 80)
    print(f"Outputs in: {OUTPUT_DIR}")


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot_preservation(summary_df, output_path):
    """Horizontal bar chart comparing preservation scores."""
    setup_plot_style()

    df = summary_df.copy()
    n = len(df)

    fig, ax = plt.subplots(figsize=(10, max(6, n * 0.45)))

    y = np.arange(n)
    bar_height = 0.35

    # Bars
    bars_funk = ax.barh(
        y + bar_height / 2,
        df["funk_preservation"],
        bar_height,
        label="Funk re-clustered",
        color=COLORS["funk"],
        alpha=0.85,
    )
    bars_brieflow = ax.barh(
        y - bar_height / 2,
        df["brieflow_preservation"],
        bar_height,
        label="Brieflow",
        color=COLORS["brieflow"],
        alpha=0.85,
    )

    # Value labels
    for bars in [bars_funk, bars_brieflow]:
        for bar in bars:
            width = bar.get_width()
            if width > 0:
                ax.text(
                    width + 0.01,
                    bar.get_y() + bar.get_height() / 2,
                    f"{width:.0%}",
                    va="center",
                    fontsize=9,
                )

    # Labels
    labels = [f"{row['group']}\n({row['n_genes']} genes)" for _, row in df.iterrows()]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Co-clustering Preservation (fraction in dominant cluster)")
    ax.set_xlim(0, 1.15)
    ax.invert_yaxis()
    ax.legend(loc="lower right", fontsize=10)
    ax.set_title(
        "Validation Gene Group Preservation Across Pipelines", fontweight="bold"
    )

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"  Saved: {output_path}")


def plot_group_preservation(group_df, output_path):
    """Two-panel plot for pathway group preservation: cluster count and entropy."""
    setup_plot_style()

    df = group_df.copy()
    n = len(df)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(5, n * 0.5)))

    y = np.arange(n)
    bar_height = 0.35
    labels = [
        f"{row['pathway_group']}\n({row['n_genes']} genes)" for _, row in df.iterrows()
    ]

    # --- Left panel: number of clusters ---
    ax1.barh(
        y + bar_height / 2,
        df["funk_n_clusters"],
        bar_height,
        label="Funk re-clustered",
        color=COLORS["funk"],
        alpha=0.85,
    )
    ax1.barh(
        y - bar_height / 2,
        df["brieflow_n_clusters"],
        bar_height,
        label="Brieflow",
        color=COLORS["brieflow"],
        alpha=0.85,
    )
    # Reference: number of original Funk published clusters
    for i, row in df.iterrows():
        ax1.plot(
            row["n_funk_published_clusters"],
            i,
            marker="d",
            color="black",
            markersize=7,
            zorder=5,
        )
    ax1.plot(
        [],
        [],
        marker="d",
        color="black",
        linestyle="none",
        label="Funk published #clusters",
    )

    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.set_xlabel("Number of Clusters")
    ax1.invert_yaxis()
    ax1.legend(loc="lower right", fontsize=9)
    ax1.set_title("Cluster Spread per Pathway Group", fontweight="bold")

    # --- Right panel: normalized entropy ---
    ax2.barh(
        y + bar_height / 2,
        df["funk_normalized_entropy"],
        bar_height,
        label="Funk re-clustered",
        color=COLORS["funk"],
        alpha=0.85,
    )
    ax2.barh(
        y - bar_height / 2,
        df["brieflow_normalized_entropy"],
        bar_height,
        label="Brieflow",
        color=COLORS["brieflow"],
        alpha=0.85,
    )

    ax2.set_yticks(y)
    ax2.set_yticklabels([], fontsize=9)
    ax2.set_xlabel("Normalized Entropy (0=tight, 1=uniform)")
    ax2.set_xlim(0, 1.0)
    ax2.invert_yaxis()
    ax2.legend(loc="lower right", fontsize=9)
    ax2.set_title("Cluster Assignment Entropy", fontweight="bold")

    fig.suptitle(
        "Pathway Group Preservation Across Pipelines",
        fontsize=14,
        fontweight="bold",
        y=1.01,
    )
    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"  Saved: {output_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    run_validation()
