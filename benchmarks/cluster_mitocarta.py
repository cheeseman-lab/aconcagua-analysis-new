"""
MitoCarta3.0 Pathway Boxes — Gene Names Colored by Cluster

For each pipeline (Brieflow, Funk), draws top-level MitoCarta3.0 pathway
boxes containing gene names colored by cluster assignment. Shows whether
unsupervised clusters recover known mitochondrial pathway structure.

Usage:
    python cluster_mitocarta.py
"""

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

import sys

sys.path.insert(0, str(Path(__file__).parent))
from plot_style import setup_plot_style, save_figure
from cluster_phate import (
    HIGHLIGHT_SETS,
    BRIEFLOW_CLUSTER_PALETTE,
    FUNK_CLUSTER_PALETTE,
    IDEAL_CONFIGS,
    load_features_for_gene_heatmap,
    load_phate_data,
)

BENCHMARKS_DIR = Path(__file__).parent
EXTERNAL_DIR = BENCHMARKS_DIR / "external"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "cluster" / "mitocarta"
MITOCARTA_PATH = EXTERNAL_DIR / "Human.MitoCarta3.0.xls"


def _get_top_level_pathway(gene_symbols):
    """Map genes to top-level MitoCarta3.0 pathway.

    Returns dict {gene: pathway_label}.
    """
    mc = pd.read_excel(MITOCARTA_PATH, sheet_name="A Human MitoCarta3.0")
    mc_map = dict(zip(mc["Symbol"], mc["MitoCarta3.0_MitoPathways"]))

    result = {}
    for gene in gene_symbols:
        pw = str(mc_map.get(gene, ""))
        if pw in ("nan", "0", ""):
            result[gene] = "Not in MitoCarta"
            continue
        first = pw.split("|")[0].strip()
        top = first.split(">")[0].strip()
        # Shorten long names
        top = top.replace("Mitochondrial central dogma", "Central Dogma")
        top = top.replace("Protein import, sorting and homeostasis",
                          "Protein Import & Homeostasis")
        top = top.replace("Mitochondrial dynamics and surveillance",
                          "Dynamics & Surveillance")
        result[gene] = top
    return result


def _build_color_map(gene_clusters, gene_list, palette):
    """Build {cluster_id: color} for clusters present in gene_list."""
    clusters_sorted = sorted(set(gene_clusters[g] for g in gene_list))
    return {cl: palette[i % len(palette)] for i, cl in enumerate(clusters_sorted)}


def plot_pathway_boxes(cell_class, highlight_key, pipeline, output_path):
    """Draw top-level MitoCarta pathway boxes with gene names colored by cluster.

    Args:
        pipeline: "brieflow" or "funk"
    """
    setup_plot_style()

    hl_set = HIGHLIGHT_SETS[cell_class][highlight_key]

    if pipeline == "brieflow":
        data = load_features_for_gene_heatmap(cell_class, hl_set, select_by="brieflow")
        gcl = data["bl_gene_clusters"]
        gene_list = data["gene_list"]
        palette = BRIEFLOW_CLUSTER_PALETTE
        hl_clusters = set(hl_set["brieflow_clusters"])
        pipeline_label = "Brieflow"
    else:
        # Load ALL genes from Funk's highlighted clusters (not pre-filtered)
        brieflow_df, funk_df, cfg = load_phate_data(cell_class)
        hl_clusters = set(hl_set["funk_clusters"])
        fk_hl_df = funk_df[funk_df["cluster"].isin(hl_clusters)]
        gene_list = fk_hl_df["gene_symbol"].tolist()
        gcl = dict(zip(funk_df["gene_symbol"], funk_df["cluster"]))
        palette = FUNK_CLUSTER_PALETTE
        pipeline_label = "Funk"

    # Only include genes belonging to highlighted clusters
    gene_list = [g for g in gene_list if gcl.get(g) in hl_clusters]
    color_map = {cl: palette[i % len(palette)]
                 for i, cl in enumerate(sorted(hl_clusters))}
    pathways = _get_top_level_pathway(gene_list)

    # Group genes by pathway
    groups = defaultdict(list)
    for gene in gene_list:
        groups[pathways[gene]].append(gene)

    # Sort groups by size
    sorted_groups = sorted(groups.items(), key=lambda x: -len(x[1]))

    # --- Layout ---
    # Each box: pathway label on top, gene names wrapped inside
    box_padding = 0.4
    gene_fontsize = 9
    header_fontsize = 12
    col_width = 3.2  # width per gene column inside box
    row_height = 0.45  # height per gene row
    max_cols = 4  # max gene columns per box
    box_gap = 0.8  # vertical gap between boxes

    # Compute box dimensions
    box_specs = []
    for pathway, genes in sorted_groups:
        n = len(genes)
        cols = min(n, max_cols)
        rows = int(np.ceil(n / cols))
        w = cols * col_width + 2 * box_padding
        h = rows * row_height + 2 * box_padding + 1.0  # extra for header
        box_specs.append({
            "pathway": pathway,
            "genes": sorted(genes),
            "cols": cols,
            "rows": rows,
            "width": w,
            "height": h,
        })

    # Stack boxes vertically
    total_width = max(b["width"] for b in box_specs) + 2
    y_cursor = 0.5
    for spec in box_specs:
        spec["y"] = y_cursor
        spec["x"] = (total_width - spec["width"]) / 2  # center
        y_cursor += spec["height"] + box_gap

    total_height = y_cursor + 0.5
    fig, ax = plt.subplots(figsize=(total_width, total_height * 0.65))

    for spec in box_specs:
        x0, y0 = spec["x"], spec["y"]
        w, h = spec["width"], spec["height"]

        # Draw box
        rect = mpatches.FancyBboxPatch(
            (x0, y0), w, h,
            boxstyle="round,pad=0.15",
            facecolor="#f8f8f8",
            edgecolor="black",
            linewidth=1.2,
        )
        ax.add_patch(rect)

        # Pathway header
        ax.text(
            x0 + w / 2, y0 + 0.5,
            f"{spec['pathway']} ({len(spec['genes'])})",
            ha="center", va="center",
            fontsize=header_fontsize, fontweight="bold",
        )

        # Gene names in grid
        for idx, gene in enumerate(spec["genes"]):
            col = idx % spec["cols"]
            row = idx // spec["cols"]
            gx = x0 + box_padding + col * col_width + col_width / 2
            gy = y0 + 1.0 + box_padding + row * row_height

            cl = gcl.get(gene)
            color = color_map.get(cl, "#999999")

            ax.text(
                gx, gy, gene,
                ha="center", va="center",
                fontsize=gene_fontsize, fontweight="bold",
                color=color,
            )

    # Legend
    legend_handles = []
    for cl in sorted(color_map.keys()):
        in_hl = cl in hl_clusters
        label = f"Cluster {cl}" + (" *" if in_hl else "")
        legend_handles.append(mpatches.Patch(color=color_map[cl], label=label))
    ax.legend(
        handles=legend_handles,
        loc="lower right",
        fontsize=10,
        framealpha=0.9,
        title=f"{pipeline_label} Clusters",
        title_fontsize=11,
    )

    ax.set_xlim(-0.5, total_width + 0.5)
    ax.set_ylim(-0.5, total_height + 0.5)
    ax.axis("off")
    fig.tight_layout()

    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def _layout_gene_names(ax, genes, gene_clusters, color_map, x0, y0, width, height,
                        rows_per_col=6, col_w_gene=1.10, row_h_gene=0.30,
                        fontsize=9):
    """Lay out gene names in column-major order within a rectangular zone.

    Genes fill top-to-bottom within each column, then wrap to the next column.
    """
    if not genes:
        return 0
    sorted_genes = sorted(genes)
    n = len(sorted_genes)
    rows = min(n, rows_per_col)

    # Center the gene block vertically
    block_h = rows * row_h_gene
    y_offset = y0 + (height - block_h) / 2

    for idx, gene in enumerate(sorted_genes):
        c = idx // rows_per_col   # column first
        r = idx % rows_per_col    # then row
        gx = x0 + c * col_w_gene + col_w_gene / 2
        gy = y_offset + r * row_h_gene + row_h_gene / 2
        color = color_map.get(gene_clusters.get(gene), "#999999")
        ax.text(gx, gy, gene, ha="center", va="center",
                fontsize=fontsize, fontweight="bold", color=color)
    return rows


def _get_overlap_data(cell_class, highlight_key):
    """Load data and split genes into Brieflow-only, overlap, and Funk-only sets.

    Returns (pw_data, bl_only, overlap, fk_only) where pw_data maps
    pathway name -> {"bl_only": [...], "overlap": [...], "fk_only": [...]}.
    """
    hl_set = HIGHLIGHT_SETS[cell_class][highlight_key]

    # Brieflow data
    data = load_features_for_gene_heatmap(cell_class, hl_set, select_by="brieflow")
    bl_gcl = data["bl_gene_clusters"]
    bl_hl = set(hl_set["brieflow_clusters"])
    bl_genes_set = set(g for g in data["gene_list"] if bl_gcl.get(g) in bl_hl)

    # Funk data
    _, funk_df, _ = load_phate_data(cell_class)
    fk_hl = set(hl_set["funk_clusters"])
    fk_hl_df = funk_df[funk_df["cluster"].isin(fk_hl)]
    fk_genes_set = set(fk_hl_df["gene_symbol"].tolist())

    # Split
    overlap = bl_genes_set & fk_genes_set
    bl_only = bl_genes_set - overlap
    fk_only = fk_genes_set - overlap

    # Map to MitoCarta pathways
    all_genes = sorted(bl_genes_set | fk_genes_set)
    pathways = _get_top_level_pathway(all_genes)

    # Group by pathway and zone
    pw_data = defaultdict(lambda: {"bl_only": [], "overlap": [], "fk_only": []})
    for g in bl_only:
        pw_data[pathways[g]]["bl_only"].append(g)
    for g in overlap:
        pw_data[pathways[g]]["overlap"].append(g)
    for g in fk_only:
        pw_data[pathways[g]]["fk_only"].append(g)

    return pw_data


def plot_combined(cell_class, highlight_key, output_path,
                  include_not_mitocarta=True):
    """GridSpec combined figure: genes colored by overlap status.

    Row 1: Central Dogma (full width).
    Row 2: OXPHOS, Metabolism, Protein Import, Dynamics side by side.
    Row 3 (optional): Not in MitoCarta.
    """
    import matplotlib.gridspec as gridspec

    setup_plot_style()

    pw_data = _get_overlap_data(cell_class, highlight_key)

    # --- Overlap-status colors ---
    BL_COLOR = "#0077BB"   # blue for Brieflow-only
    OV_COLOR = "#882255"   # purple for Both
    FK_COLOR = "#EE7733"   # orange for Funk-only

    # --- Layout parameters ---
    rows_per_col = 6
    col_w = 1.10
    row_h = 0.28
    gene_fs = 10
    title_fs = 14

    # --- Helpers ---
    def _cat_genes(name):
        d = pw_data.get(name, {"bl_only": [], "overlap": [], "fk_only": []})
        return sorted(d["bl_only"]), sorted(d["overlap"]), sorted(d["fk_only"])

    def _n_cols(genes):
        return int(np.ceil(len(genes) / rows_per_col)) if genes else 0

    def _total_cols(bl, ov, fk):
        return _n_cols(bl) + _n_cols(ov) + _n_cols(fk)

    def _draw_panel(ax, title, bl, ov, fk):
        """Draw gene names in a single grid, colored by overlap status."""
        # Build ordered list: BL genes, then Both, then FK — each sub-sorted
        gene_list = ([(g, BL_COLOR) for g in bl]
                     + [(g, OV_COLOR) for g in ov]
                     + [(g, FK_COLOR) for g in fk])
        n = len(gene_list)
        if n == 0:
            ax.set_title(title, fontsize=title_fs, fontweight="bold")
            ax.axis("off")
            return

        total_c = int(np.ceil(n / rows_per_col))
        rows = min(n, rows_per_col)
        w = total_c * col_w
        h = rows * row_h

        ax.set_xlim(-0.05, w + 0.05)
        ax.set_ylim(-0.05, h + 0.05)
        ax.set_title(title, fontsize=title_fs, fontweight="bold", pad=8)
        ax.axis("off")

        for idx, (gene, color) in enumerate(gene_list):
            c = idx // rows_per_col
            r = idx % rows_per_col
            gx = c * col_w + col_w / 2
            gy = h - (r * row_h + row_h / 2)
            ax.text(gx, gy, gene, ha="center", va="center",
                    fontsize=gene_fs, fontweight="bold", color=color)

    # --- Category data ---
    ROW2_NAMES = ["OXPHOS", "Metabolism", "Protein Import & Homeostasis",
                  "Dynamics & Surveillance"]
    cat_data = {}
    for name in ["Central Dogma"] + ROW2_NAMES + ["Not in MitoCarta"]:
        cat_data[name] = _cat_genes(name)

    # --- Width ratios for row 2 (proportional to total gene columns) ---
    row2_cols = [max(_total_cols(*cat_data[n]), 1) for n in ROW2_NAMES]

    # --- Figure dimensions ---
    # Scale: each gene column ≈ col_w inches, plus padding
    row2_total_cols = sum(row2_cols)
    row1_cols = max(_total_cols(*cat_data["Central Dogma"]), 1)
    nic_cols = max(_total_cols(*cat_data["Not in MitoCarta"]), 1)
    max_cols = max(row2_total_cols + 4, row1_cols, nic_cols)  # +4 for row2 gaps
    fig_w = max_cols * col_w + 1.5

    panel_h = rows_per_col * row_h + 0.7  # gene grid + title
    n_rows = 3 if include_not_mitocarta else 2
    height_ratios = [panel_h, panel_h]
    if include_not_mitocarta:
        height_ratios.append(panel_h)
    fig_h = sum(height_ratios) + 0.6  # legend

    fig = plt.figure(figsize=(fig_w, fig_h))
    outer = gridspec.GridSpec(n_rows, 1, figure=fig, height_ratios=height_ratios,
                              hspace=0.45)

    # Row 1: Central Dogma
    ax_cd = fig.add_subplot(outer[0])
    _draw_panel(ax_cd, "Central Dogma", *cat_data["Central Dogma"])

    # Row 2: four panels side by side
    inner = gridspec.GridSpecFromSubplotSpec(
        1, 4, subplot_spec=outer[1],
        width_ratios=row2_cols, wspace=0.4,
    )
    for i, name in enumerate(ROW2_NAMES):
        ax = fig.add_subplot(inner[i])
        _draw_panel(ax, name, *cat_data[name])

    # Row 3: Not in MitoCarta (optional)
    if include_not_mitocarta:
        ax_nic = fig.add_subplot(outer[2])
        _draw_panel(ax_nic, "Not in MitoCarta", *cat_data["Not in MitoCarta"])

    # --- Simple 3-item legend ---
    legend_handles = [
        mpatches.Patch(facecolor=BL_COLOR, label="Brieflow-only"),
        mpatches.Patch(facecolor=OV_COLOR, label="Both"),
        mpatches.Patch(facecolor=FK_COLOR, label="Funk-only"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=3,
               fontsize=11, framealpha=0.9, bbox_to_anchor=(0.5, 0.005))

    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_cluster_gene_list(cell_class, pipeline, cluster_id, highlight_genes,
                           output_path, n_cols=4):
    """Render a compact SVG with gene names from one cluster in a 4-column grid.

    MitoCarta genes are bolded. Genes in highlight_genes get a yellow background.
    Text is output as live text (not outlines) for easy Illustrator editing.
    """
    setup_plot_style()
    import matplotlib
    matplotlib.rcParams["svg.fonttype"] = "none"

    # Load data
    brieflow_df, funk_df, cfg = load_phate_data(cell_class)

    if pipeline == "brieflow":
        df = brieflow_df
        palette = BRIEFLOW_CLUSTER_PALETTE
        pipeline_label = "Brieflow"
    else:
        df = funk_df
        palette = FUNK_CLUSTER_PALETTE
        pipeline_label = "Funk"

    # Filter genes for this cluster
    cluster_df = df[df["cluster"] == cluster_id]
    genes = sorted(cluster_df["gene_symbol"].tolist())
    if not genes:
        print(f"No genes found for {pipeline_label} cluster {cluster_id}")
        return

    # MitoCarta membership
    pathways = _get_top_level_pathway(genes)
    mito_genes = {g for g in genes if pathways[g] != "Not in MitoCarta"}

    # Cluster color from palette
    all_clusters = sorted(df["cluster"].unique())
    cl_idx = list(all_clusters).index(cluster_id) if cluster_id in all_clusters else 0
    text_color = palette[cl_idx % len(palette)]

    # Layout parameters — minimal spacing
    n = len(genes)
    n_rows = int(np.ceil(n / n_cols))
    col_w = 0.68
    row_h = 0.18
    fontsize = 8
    pad_top = 0.35
    pad_bottom = 0.08
    pad_lr = 0.08

    fig_w = n_cols * col_w + 2 * pad_lr
    fig_h = n_rows * row_h + pad_top + pad_bottom

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)
    ax.axis("off")

    # Title
    ax.text(fig_w / 2, fig_h - pad_top / 2,
            f"{pipeline_label} Cluster {cluster_id}",
            ha="center", va="center", fontsize=10, fontweight="bold")

    # Gene names in column-major order
    highlight_set = set(highlight_genes)
    for idx, gene in enumerate(genes):
        col = idx // n_rows
        row = idx % n_rows
        gx = pad_lr + col * col_w + col_w / 2
        gy = fig_h - pad_top - row * row_h - row_h / 2

        is_mito = gene in mito_genes
        is_highlight = gene in highlight_set
        weight = "bold" if is_mito else "normal"

        if is_highlight:
            ax.text(gx, gy, gene, ha="center", va="center",
                    fontsize=fontsize, fontweight=weight, color=text_color,
                    bbox=dict(boxstyle="round,pad=0.08", facecolor="#FFFF00",
                              edgecolor="none", alpha=0.7))
        else:
            ax.text(gx, gy, gene, ha="center", va="center",
                    fontsize=fontsize, fontweight=weight, color=text_color)

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="MitoCarta3.0 pathway boxes with cluster-colored gene names"
    )
    parser.add_argument(
        "--cell-class", default="Interphase", choices=list(IDEAL_CONFIGS.keys()),
    )
    parser.add_argument("--highlight-key", default="mitochondrial")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cc = args.cell_class.lower()
    key = args.highlight_key

    for pipeline in ("brieflow", "funk"):
        plot_pathway_boxes(
            args.cell_class, key, pipeline,
            OUTPUT_DIR / f"mitocarta_{cc}_{key}_{pipeline}.pdf",
        )

    plot_combined(
        args.cell_class, key,
        OUTPUT_DIR / f"mitocarta_{cc}_{key}_combined.pdf",
        include_not_mitocarta=True,
    )
    plot_combined(
        args.cell_class, key,
        OUTPUT_DIR / f"mitocarta_{cc}_{key}_combined_mito_only.pdf",
        include_not_mitocarta=False,
    )

    # Single-cluster gene list SVGs
    highlight_map = {
        ("funk", 149): ["KRAS", "BRAF"],
        ("funk", 147): ["CAP1", "CFL1", "WDR1"],
    }
    cols_override = {
        ("brieflow", 22): 8,
    }
    hl_set = HIGHLIGHT_SETS["Interphase"]["mitochondrial"]
    cluster_specs = (
        [{"pipeline": "brieflow", "cluster": c} for c in hl_set["brieflow_clusters"]]
        + [{"pipeline": "funk", "cluster": c} for c in hl_set["funk_clusters"]]
    )
    for spec in cluster_specs:
        plot_cluster_gene_list(
            cell_class="Interphase",
            pipeline=spec["pipeline"],
            cluster_id=spec["cluster"],
            highlight_genes=highlight_map.get(
                (spec["pipeline"], spec["cluster"]), []
            ),
            output_path=OUTPUT_DIR / f"genelist_{spec['pipeline']}_cluster{spec['cluster']}.svg",
            n_cols=cols_override.get(
                (spec["pipeline"], spec["cluster"]), 4
            ),
        )


if __name__ == "__main__":
    main()
