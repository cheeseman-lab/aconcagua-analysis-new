"""Regenerate eval figures for aconcagua.

Re-renders evaluation PNGs and PDFs at 300 DPI with transparent backgrounds
and sans-serif fonts. Covers SBS, phenotype, merge, and aggregate eval.

Outputs to results/eval_figures/ mirroring the original brieflow_output/
directory structure.

Usage:
    sbatch replot_eval_figures.sh
"""

import gc
import glob
import random
import sys
from pathlib import Path

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.dataset as ds

sys.path.insert(0, "../brieflow/workflow")

from lib.shared.eval import plot_plate_heatmap
from lib.sbs.standardize_barcode_design import get_barcode_list
from lib.sbs.eval_mapping import (
    plot_mapping_vs_threshold,
    plot_read_mapping_heatmap,
    plot_cell_mapping_heatmap,
    plot_cell_metric_histogram,
    plot_gene_symbol_histogram,
    plot_barcode_prefix_matching,
)
from lib.merge.eval_merge import plot_sbs_ph_matching_heatmap, plot_cell_positions
from lib.aggregate.eval_aggregate import nas_summary, plot_feature_distributions
from lib.shared.file_utils import validate_dtypes

# --- Configuration ---
_arial_path = Path("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")
if _arial_path.exists():
    fm.fontManager.addfont(str(_arial_path))
    for variant in _arial_path.parent.glob("Arial*.ttf"):
        fm.fontManager.addfont(str(variant))

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Nimbus Sans", "Liberation Sans", "DejaVu Sans"],
    }
)

SAVE_KWARGS = dict(dpi=300, bbox_inches="tight", transparent=True)
SRC = Path("../analysis/brieflow_output")
OUT = Path("results/eval_figures")
PLATES = range(1, 9)
WELLS = ["A1", "A2", "A3", "B1", "B2", "B3"]
CHANNELS = ["DAPI", "TUBULIN", "GH2AX", "PHALLOIDIN"]
CELL_CLASSES = ["all", "Interphase", "Mitotic"]
CHANNEL_COMBO = "DAPI_TUBULIN_GH2AX_PHALLOIDIN"
SUBSET_SIZE = 100000

BARCODE_LIBRARY_FP = Path("../analysis/config/barcode_library.tsv")
BARCODE_TYPE = "simple"
SORT_BY = "count"
LIBRARY_BARCODE_COL = "prefix"
SBS_HEATMAP_SHAPE = "6W_sbs"
PH_HEATMAP_SHAPE = "6W_ph"
HEATMAP_PLATE = "6W"
MERGE_CMAP = "magma"
N_PCS = 10


def save_and_close(fig, path):
    """Save figure as PNG and PDF, then close to free memory."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"  Saving {path}")
    fig.savefig(path, **SAVE_KWARGS)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", transparent=True)
    plt.close(fig)


# =============================================================================
# SECTION 1: TSV-based heatmaps (fast, no parquets needed)
# =============================================================================
print("=" * 60)
print("SECTION 1: TSV-based heatmaps")
print("=" * 60)

for p in PLATES:
    tsv_path = SRC / f"sbs/eval/segmentation/P-{p}__cell_density_heatmap.tsv"
    if tsv_path.exists():
        df = pd.read_csv(tsv_path, sep="\t")
        fig, _ = plot_plate_heatmap(df, shape=SBS_HEATMAP_SHAPE, plate=HEATMAP_PLATE)
        save_and_close(fig, OUT / f"sbs/eval/segmentation/P-{p}__cell_density_heatmap.png")

for p in PLATES:
    tsv_path = SRC / f"phenotype/eval/segmentation/P-{p}__cell_density_heatmap.tsv"
    if tsv_path.exists():
        df = pd.read_csv(tsv_path, sep="\t")
        fig, _ = plot_plate_heatmap(df, shape=PH_HEATMAP_SHAPE, plate=HEATMAP_PLATE)
        save_and_close(fig, OUT / f"phenotype/eval/segmentation/P-{p}__cell_density_heatmap.png")

for p in PLATES:
    for suffix in ["cell_mapping_heatmap_one", "cell_mapping_heatmap_any"]:
        tsv_path = SRC / f"sbs/eval/mapping/P-{p}__{suffix}.tsv"
        if tsv_path.exists():
            df = pd.read_csv(tsv_path, sep="\t")
            fig, _ = plot_plate_heatmap(df, shape=SBS_HEATMAP_SHAPE, plate=HEATMAP_PLATE)
            save_and_close(fig, OUT / f"sbs/eval/mapping/P-{p}__{suffix}.png")

for p in PLATES:
    for suffix, shape in [
        ("sbs_to_ph_matching_rates", SBS_HEATMAP_SHAPE),
        ("ph_to_sbs_matching_rates", PH_HEATMAP_SHAPE),
    ]:
        tsv_path = SRC / f"merge/eval/P-{p}__{suffix}.tsv"
        if tsv_path.exists():
            df = pd.read_csv(tsv_path, sep="\t")
            fig, _ = plot_plate_heatmap(df, shape=shape, plate=HEATMAP_PLATE, cmap=MERGE_CMAP)
            save_and_close(fig, OUT / f"merge/eval/P-{p}__{suffix}.png")

for p in PLATES:
    for ch in CHANNELS:
        tsv_path = SRC / f"phenotype/eval/features/P-{p}__cell_{ch}_min_heatmap.tsv"
        if tsv_path.exists():
            df = pd.read_csv(tsv_path, sep="\t")
            fig, _ = plot_plate_heatmap(df, shape=PH_HEATMAP_SHAPE, plate=HEATMAP_PLATE)
            save_and_close(fig, OUT / f"phenotype/eval/features/P-{p}__cell_{ch}_min_heatmap.png")

print("Section 1 complete.\n")


# =============================================================================
# SECTION 2: SBS eval (parquet-based)
# =============================================================================
print("=" * 60)
print("SECTION 2: SBS eval (parquet-based)")
print("=" * 60)

df_barcode_library = pd.read_csv(BARCODE_LIBRARY_FP, sep="\t")
barcodes = get_barcode_list(df_barcode_library)

for p in PLATES:
    print(f"\nPlate {p}:")
    out_dir = OUT / "sbs/eval/mapping"

    reads_paths = sorted(glob.glob(str(SRC / f"sbs/parquets/P-{p}_W-*__reads.parquet")))
    cells_paths = sorted(glob.glob(str(SRC / f"sbs/parquets/P-{p}_W-*__cells.parquet")))
    sbs_info_paths = sorted(glob.glob(str(SRC / f"sbs/parquets/P-{p}_W-*__sbs_info.parquet")))

    if not reads_paths or not cells_paths:
        print(f"  Skipping plate {p} — missing parquets")
        continue

    reads = pd.concat([pd.read_parquet(f) for f in reads_paths], ignore_index=True)
    cells = pd.concat([pd.read_parquet(f) for f in cells_paths], ignore_index=True)
    sbs_info = pd.concat([pd.read_parquet(f) for f in sbs_info_paths], ignore_index=True)

    _, fig = plot_mapping_vs_threshold(reads, barcodes, "peak", num_thresholds=10)
    save_and_close(fig, out_dir / f"P-{p}__mapping_vs_threshold_peak.png")

    _, fig = plot_mapping_vs_threshold(reads, barcodes, "Q_min", num_thresholds=10)
    save_and_close(fig, out_dir / f"P-{p}__mapping_vs_threshold_qmin.png")

    fig = plot_read_mapping_heatmap(reads, barcodes, plate=HEATMAP_PLATE, shape=SBS_HEATMAP_SHAPE)
    save_and_close(fig, out_dir / f"P-{p}__read_mapping_heatmap.png")

    _, fig = plot_cell_mapping_heatmap(
        cells, sbs_info, barcodes, mapping_to="one",
        mapping_strategy="gene symbols", shape=SBS_HEATMAP_SHAPE,
        plate=HEATMAP_PLATE, return_summary=True,
    )
    save_and_close(fig, out_dir / f"P-{p}__cell_mapping_heatmap_one.png")

    _, fig = plot_cell_mapping_heatmap(
        cells, sbs_info, barcodes, mapping_to="any",
        mapping_strategy="gene symbols", shape=SBS_HEATMAP_SHAPE,
        plate=HEATMAP_PLATE, return_summary=True,
    )
    save_and_close(fig, out_dir / f"P-{p}__cell_mapping_heatmap_any.png")

    _, fig = plot_cell_metric_histogram(cells, sort_by=SORT_BY)
    save_and_close(fig, out_dir / f"P-{p}__cell_metric_histogram.png")

    _, fig = plot_gene_symbol_histogram(cells)
    save_and_close(fig, out_dir / f"P-{p}__gene_symbol_histogram.png")

    _, fig = plot_barcode_prefix_matching(reads, df_barcode_library, library_col=LIBRARY_BARCODE_COL)
    save_and_close(fig, out_dir / f"P-{p}__barcode_prefix_matching.png")

    del reads, cells, sbs_info
    gc.collect()

print("Section 2 complete.\n")


# =============================================================================
# SECTION 3: Merge eval (parquet-based)
# =============================================================================
print("=" * 60)
print("SECTION 3: Merge eval (parquet-based)")
print("=" * 60)

for p in PLATES:
    print(f"\nPlate {p}:")
    out_dir = OUT / "merge/eval"

    dedup_paths = sorted(
        glob.glob(str(SRC / f"merge/parquets/P-{p}_W-*__merge_deduplicated.parquet"))
    )
    if not dedup_paths:
        print(f"  Skipping plate {p} — no merge_deduplicated parquets")
        continue

    merge_deduplicated = validate_dtypes(
        pd.concat([pd.read_parquet(f) for f in dedup_paths], ignore_index=True)
    )
    if "global_i_0" in merge_deduplicated.columns:
        merge_deduplicated = merge_deduplicated.rename(
            columns={
                "global_i_0": "i_0", "global_j_0": "j_0",
                "global_i_1": "i_1", "global_j_1": "j_1",
            }
        )

    fig = plot_cell_positions(merge_deduplicated, title="All Cells by Channel Min")
    save_and_close(fig, out_dir / f"P-{p}__all_cells_by_channel_min.png")

    fig = plot_cell_positions(
        merge_deduplicated.query("channels_min==0"),
        title="Cells with Channel Min = 0", color="red",
    )
    save_and_close(fig, out_dir / f"P-{p}__cells_with_channel_min_0.png")

    formatted_paths = sorted(
        glob.glob(str(SRC / f"merge/parquets/P-{p}_W-*__merge_formatted.parquet"))
    )
    sbs_info_paths = sorted(
        glob.glob(str(SRC / f"sbs/parquets/P-{p}_W-*__sbs_info.parquet"))
    )
    ph_info_paths = sorted(
        glob.glob(str(SRC / f"phenotype/parquets/P-{p}_W-*__phenotype_info.parquet"))
    )

    if formatted_paths and sbs_info_paths and ph_info_paths:
        merge_formatted = validate_dtypes(
            pd.concat([pd.read_parquet(f) for f in formatted_paths], ignore_index=True)
        )
        sbs_info = validate_dtypes(
            pd.concat([pd.read_parquet(f) for f in sbs_info_paths], ignore_index=True)
        )
        ph_info = validate_dtypes(
            pd.concat([pd.read_parquet(f) for f in ph_info_paths], ignore_index=True)
        )
        merge_minimal = merge_formatted[
            ["plate", "well", "tile", "site", "cell_0", "cell_1", "distance"]
        ]
        sbs_summary, ax = plot_sbs_ph_matching_heatmap(
            merge_minimal, sbs_info.rename(columns={"cell": "cell_1"}),
            target="sbs", shape=SBS_HEATMAP_SHAPE, plate=HEATMAP_PLATE, return_summary=True,
        )
        save_and_close(ax.get_figure(), out_dir / f"P-{p}__sbs_to_ph_matching_rates.png")
        ph_summary, ax = plot_sbs_ph_matching_heatmap(
            merge_minimal, ph_info.rename(columns={"cell": "cell_0"}),
            target="phenotype", shape=PH_HEATMAP_SHAPE, plate=HEATMAP_PLATE, return_summary=True,
        )
        save_and_close(ax.get_figure(), out_dir / f"P-{p}__ph_to_sbs_matching_rates.png")
        del merge_formatted, sbs_info, ph_info, merge_minimal
    else:
        print(f"  Skipping matching rate heatmaps for plate {p} — missing parquets")

    del merge_deduplicated
    gc.collect()

print("Section 3 complete.\n")


# =============================================================================
# SECTION 4: Aggregate eval (parquet-based)
# =============================================================================
print("=" * 60)
print("SECTION 4: Aggregate eval (parquet-based)")
print("=" * 60)

for cell_class in CELL_CLASSES:
    prefix = f"CeCl-{cell_class}_ChCo-{CHANNEL_COMBO}"
    print(f"\n{prefix}:")
    out_dir = OUT / "aggregate/eval"

    split_paths = sorted(glob.glob(str(
        SRC / f"aggregate/parquets/P-*_W-*_CeCl-{cell_class}_ChCo-{CHANNEL_COMBO}__merge_data.parquet"
    )))
    aligned_path = SRC / f"aggregate/parquets/{prefix}__aligned.parquet"

    if not split_paths:
        print("  Skipping — no merge_data parquets")
        continue
    if not aligned_path.exists():
        print("  Skipping — no aligned parquet")
        continue

    merge_data = ds.dataset(split_paths, format="parquet")
    total_rows = merge_data.count_rows()
    n_sample = min(SUBSET_SIZE, total_rows)
    random_indices = np.random.choice(total_rows, size=n_sample, replace=False)
    random_indices.sort()
    merge_data = merge_data.scanner().take(random_indices).to_pandas(
        use_threads=True, memory_pool=None
    )

    nas_df, nas_fig = nas_summary(merge_data, vis_subsample=50000)
    if nas_fig is not None:
        save_and_close(nas_fig, out_dir / f"{prefix}__na_stats.png")
    else:
        print("  No NAs found — skipping na_stats figure")

    aligned_data = ds.dataset(str(aligned_path), format="parquet")
    total_rows = aligned_data.count_rows()
    n_sample = min(SUBSET_SIZE, total_rows)
    random_indices = np.random.choice(total_rows, size=n_sample, replace=False)
    random_indices.sort()
    aligned_data = aligned_data.scanner().take(random_indices).to_pandas(
        use_threads=True, memory_pool=None
    )

    random.seed(42)
    merge_feature_cols = [c for c in merge_data.columns if "cell_" in c and c.endswith("_mean")]
    if not merge_feature_cols:
        merge_feature_cols = [c for c in merge_data.columns if "nucleus_" in c and c.endswith("_mean")]
    pc_cols = [c for c in aligned_data.columns if c.startswith("PC_")]
    aligned_feature_cols = random.sample(pc_cols, k=min(len(merge_feature_cols), len(pc_cols)))

    feature_distributions_fig = plot_feature_distributions(
        merge_feature_cols, merge_data, aligned_feature_cols, aligned_data,
    )
    save_and_close(feature_distributions_fig, out_dir / f"{prefix}__feature_distributions.png")

    del merge_data, aligned_data
    gc.collect()

print("\nSection 4 complete.\n")


# =============================================================================
# SECTION 5: TVN summary TSV for source data (Fig3d)
# Loads only needed columns per well — low memory footprint.
# =============================================================================
print("=" * 60)
print("SECTION 5: TVN summary TSV")
print("=" * 60)

for cell_class in ["Interphase", "Mitotic"]:
    prefix = f"CeCl-{cell_class}_ChCo-{CHANNEL_COMBO}"
    print(f"\n{prefix}:")
    out_dir = OUT / "aggregate/eval"
    out_dir.mkdir(parents=True, exist_ok=True)

    split_paths = sorted(glob.glob(str(
        SRC / f"aggregate/parquets/P-*_W-*_CeCl-{cell_class}_ChCo-{CHANNEL_COMBO}__merge_data.parquet"
    )))
    aligned_path = SRC / f"aggregate/parquets/{prefix}__aligned.parquet"

    if not split_paths or not aligned_path.exists():
        print("  Skipping — missing parquets")
        continue

    summary_rows = []
    raw_feature_cols = None

    for p in split_paths:
        df = pd.read_parquet(p, columns=None)
        if raw_feature_cols is None:
            raw_feature_cols = [c for c in df.columns if "cell_" in c and c.endswith("_mean")]
            if not raw_feature_cols:
                raw_feature_cols = [c for c in df.columns if "nucleus_" in c and c.endswith("_mean")]
        plate = df["plate"].iloc[0]
        well = df["well"].iloc[0]
        for col in raw_feature_cols:
            vals = df[col].dropna()
            summary_rows.append({
                "cell_class": cell_class, "plate": plate, "well": well,
                "feature_name": col, "feature_type": "raw",
                "median": round(float(vals.median()), 4),
                "q25": round(float(vals.quantile(0.25)), 4),
                "q75": round(float(vals.quantile(0.75)), 4),
                "n_cells": len(vals),
            })
        del df
        gc.collect()

    pc_load_cols = ["plate", "well"] + [f"PC_{i}" for i in range(N_PCS)]
    aligned_df = pd.read_parquet(str(aligned_path), columns=pc_load_cols)
    for col in [c for c in pc_load_cols if c.startswith("PC_")]:
        for (plate, well), grp in aligned_df.groupby(["plate", "well"]):
            vals = grp[col].dropna()
            summary_rows.append({
                "cell_class": cell_class, "plate": plate, "well": well,
                "feature_name": col, "feature_type": "aligned_pc",
                "median": round(float(vals.median()), 4),
                "q25": round(float(vals.quantile(0.25)), 4),
                "q75": round(float(vals.quantile(0.75)), 4),
                "n_cells": len(vals),
            })
    del aligned_df
    gc.collect()

    tsv_out = out_dir / f"{prefix}__tvn_summary.tsv"
    pd.DataFrame(summary_rows).to_csv(tsv_out, sep="\t", index=False)
    print(f"  Saved: {tsv_out.name} ({len(summary_rows)} rows)")

print("\nSection 5 complete.\n")
print("=" * 60)
print("All done!")
print("=" * 60)
