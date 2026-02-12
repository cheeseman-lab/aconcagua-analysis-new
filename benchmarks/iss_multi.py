"""
ISS Multi-Construct Validation Benchmark

Validates brieflow's SBS pipeline using multi-construct ISS data with 3 known
barcodes. All wells contain the same pool of 3 constructs, providing a direct
test of barcode detection specificity with 3 sequencing cycles.

Data Source: benchmarks/results/iss_multi/
- 3 cycles (c1-c3), no skip
- 5 channels per ND2: [DAPI, CY3_T7, AF594_T7, CY5_T7, CY7_T7] → [DAPI, G, T, A, C]
- 3 wells (B03, B06, C03), 3 positions (FOVs) each
- Pool of 3 constructs with known 12-mer barcodes:
    ATACGGTTACTC → 3-base prefix: ATA
    TGATTAGAGGCA → 3-base prefix: TGA
    CATTACCTGCAA → 3-base prefix: CAT

Each ND2 file has shape (3, 5, 2960, 2960) = (Positions, Channels, Y, X).

Usage:
    python iss_multi.py                          # well B03, position 0
    python iss_multi.py --well B06 --position 1  # specific well/position
    python iss_multi.py --positions all           # all 3 positions
    python iss_multi.py --all-wells               # all wells, all positions
    python iss_multi.py --all-wells --positions all --verbose
    python iss_multi.py --cross-well              # mix reads across wells
    python iss_multi.py --plots-only              # regenerate plots from saved data
    python iss_multi.py --rerun --threshold 400   # re-run from cached arrays at new threshold
    python iss_multi.py --all-wells --positions all --rerun --threshold 400
"""

import sys
import argparse
from pathlib import Path

import nd2
import numpy as np
import pandas as pd
import tifffile

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------
BENCHMARKS_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BENCHMARKS_DIR.parent
BRIEFLOW_WORKFLOW = PROJECT_ROOT / "brieflow" / "workflow"
sys.path.insert(0, str(BRIEFLOW_WORKFLOW))

from lib.sbs.align_cycles import align_cycles
from lib.shared.log_filter import log_filter
from lib.sbs.compute_standard_deviation import compute_standard_deviation
from lib.sbs.find_peaks import find_peaks
from lib.sbs.max_filter import max_filter
from lib.shared.segment_cellpose import segment_cellpose
from lib.sbs.extract_bases import extract_bases
from lib.sbs.call_reads import call_reads
from lib.sbs.call_cells import call_cells

# ---------------------------------------------------------------------------
# PARAMS — pipeline configuration (matches aconcagua production)
# ---------------------------------------------------------------------------
PARAMS = {
    "channel_names": ["DAPI", "G", "T", "A", "C"],
    "bases": ["G", "T", "A", "C"],
    "alignment_method": "sbs_mean",
    "upsample_factor": 1,
    "max_filter_width": 3,
    "peak_width": 5,
    "threshold_peaks": 400,
    "call_reads_method": "median",
    "dapi_index": 0,
    "cyto_index": 4,
    "cyto_cycle_index": 0,
    "cellpose_model": "cyto3",
    "gpu": True,
    "segment_cells": False,
    "reconcile": "contained_in_cells",
    "q_min": 0,
    "sort_calls": "peak",
    "error_correct": False,
    "prefix_col": "prefix",
}

_EXTRA_CHANNEL_INDICES = [
    i for i, ch in enumerate(PARAMS["channel_names"]) if ch not in PARAMS["bases"]
]
_BASE_CHANNEL_INDICES = [
    i for i, ch in enumerate(PARAMS["channel_names"]) if ch in PARAMS["bases"]
]

# ---------------------------------------------------------------------------
# DATA_CONFIG
# ---------------------------------------------------------------------------
DATA_CONFIG = {
    "data_root": BENCHMARKS_DIR / "results" / "iss_multi",
    "output_dir": BENCHMARKS_DIR / "results" / "iss_multi",
    "cycle_dirs": ["c1", "c2", "c3"],
    "skip_cycles": [],
    "n_positions": 3,
    "well_seq": {
        "B03": "Seq0001",
        "B06": "Seq0004",
        "C03": "Seq0010",
    },
    "all_wells": ["B03", "B06", "C03"],
    "known_barcodes_12mer": {
        "ATACGGTTACTC": "construct_1",
        "TGATTAGAGGCA": "construct_2",
        "CATTACCTGCAA": "construct_3",
    },
    "known_prefixes": ["ATA", "TGA", "CAT"],
}

DF_BARCODE_LIBRARY = pd.DataFrame({
    "gene_symbol": ["construct_1", "construct_2", "construct_3"],
    "prefix": ["ATA", "TGA", "CAT"],
})

_ND2_CHANNEL_STR = "DAPI,CY3_T7,AF594_T7,CY5_T7,CY7_T7"
_SUMMARY_FILENAME = "iss_multi_summary.tsv"

_VERBOSE = False


def vprint(*args, **kwargs):
    if _VERBOSE:
        print(*args, **kwargs)


def _hamming(a: str, b: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(a, b))


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_position_across_cycles(well: str, position: int) -> list[np.ndarray]:
    """Load one position from one well across all SBS cycles."""
    root = DATA_CONFIG["data_root"]
    cycle_dirs = DATA_CONFIG["cycle_dirs"]
    skip = set(DATA_CONFIG["skip_cycles"])
    seq = DATA_CONFIG["well_seq"][well]

    images: list[np.ndarray] = []

    vprint(f"\n--- Data Loading ---")
    vprint(f"  Root:     {root}")
    vprint(f"  Well:     {well}  |  Position: {position}")
    vprint(f"  Seq:      {seq}")

    for idx, cdir in enumerate(cycle_dirs):
        cycle_num = idx + 1
        if idx in skip:
            vprint(f"  cycle {cycle_num}  — SKIPPED")
            continue

        nd2_name = f"Well{well}_Channel{_ND2_CHANNEL_STR}_{seq}.nd2"
        filepath = root / cdir / nd2_name
        if not filepath.exists():
            raise FileNotFoundError(f"Missing: {filepath}")

        vprint(f"  cycle {cycle_num} ({cdir})  loading: {filepath.name}")

        with nd2.ND2File(str(filepath)) as f:
            data = f.asarray()

        img = data[position]
        images.append(img)

        if _VERBOSE:
            means = [f"{img[c].mean():.1f}" for c in range(img.shape[0])]
            print(f"    shape={img.shape}  ch_means={means}")

    print(
        f"\nLoaded {len(images)} cycles for {well} pos {position}  "
        f"(shape: {images[0].shape})"
    )
    return images


# ---------------------------------------------------------------------------
# SBS pipeline
# ---------------------------------------------------------------------------
def run_sbs_pipeline(
    cycle_images: list[np.ndarray],
    well: str,
    position: int,
) -> dict:
    """Run the full SBS pipeline on one position's images."""
    p = PARAMS
    ch_names = p["channel_names"]
    extra_idx = _EXTRA_CHANNEL_INDICES
    wildcards = {"well": well, "tile": position}

    stack = np.stack(cycle_images, axis=0)
    print(f"\n1. Stack           shape={stack.shape}  dtype={stack.dtype}")
    vprint(f"   channels: {ch_names}")
    vprint(f"   extra_idx={extra_idx}  base_idx={_BASE_CHANNEL_INDICES}")

    aligned = align_cycles(
        stack,
        channel_order=ch_names,
        method=p["alignment_method"],
        upsample_factor=p["upsample_factor"],
        verbose=_VERBOSE,
    )
    print(f"2. Align           shape={aligned.shape}  method={p['alignment_method']}")

    loged = log_filter(aligned, sigma=1, skip_index=extra_idx)
    print(f"3. LoG filter      shape={loged.shape}  skipped={extra_idx}")

    std = compute_standard_deviation(loged, remove_index=extra_idx)
    print(f"4. Std dev         shape={std.shape}")

    peaks = find_peaks(std, width=p["peak_width"])
    n_above = int(np.sum(peaks > p["threshold_peaks"]))
    print(
        f"5. Find peaks      shape={peaks.shape}  "
        f"above threshold({p['threshold_peaks']})={n_above}"
    )

    maxed = max_filter(loged, width=p["max_filter_width"], remove_index=extra_idx)
    print(f"6. Max filter      shape={maxed.shape}")

    seg_cycle = aligned[p["cyto_cycle_index"]]
    print(
        f"7. Segmentation    cycle_idx={p['cyto_cycle_index']}  "
        f"dapi_idx={p['dapi_index']}  cyto_idx={p['cyto_index']}  "
        f"model={p['cellpose_model']}  gpu={p['gpu']}"
    )
    result = segment_cellpose(
        seg_cycle,
        dapi_index=p["dapi_index"],
        cyto_index=p["cyto_index"],
        nuclei_diameter=None,
        cell_diameter=None,
        cellpose_model=p["cellpose_model"],
        cellpose_kwargs=dict(
            nuclei_flow_threshold=0.4,
            nuclei_cellprob_threshold=0.0,
            cell_flow_threshold=1,
            cell_cellprob_threshold=0,
        ),
        gpu=p["gpu"],
        reconcile=p["reconcile"],
        cells=p["segment_cells"],
    )
    if p["segment_cells"]:
        nuclei, cells = result
    else:
        nuclei = result
        cells = nuclei  # Use nuclei as cells for downstream code compatibility

    n_nuclei = len(np.unique(nuclei)) - 1
    n_cells = len(np.unique(cells)) - 1
    print(f"   nuclei={n_nuclei}  cells={n_cells}")

    # Save arrays for cheap re-runs at different thresholds
    pipeline_dir = DATA_CONFIG["data_root"] / "pipeline"
    pos_dir = pipeline_dir / f"{well}_pos{position}"
    pos_dir.mkdir(parents=True, exist_ok=True)
    for name, arr in [("aligned", aligned), ("nuclei", nuclei), ("peaks", peaks), ("maxed", maxed)]:
        tifffile.imwrite(str(pos_dir / f"{name}.tiff"), arr)
    vprint(f"   Saved arrays to {pos_dir}")

    df_bases = extract_bases(
        peaks_data=peaks,
        max_filtered_data=maxed,
        cells_data=nuclei,
        threshold_peaks=p["threshold_peaks"],
        wildcards=wildcards,
        bases="".join(p["bases"]),
    )
    n_reads = df_bases["read"].nunique() if not df_bases.empty else 0
    print(f"8. Extract bases   reads={n_reads}")
    if _VERBOSE and not df_bases.empty:
        cycles_in = sorted(df_bases["cycle"].unique())
        print(f"   cycles in df: {cycles_in}  (n={len(cycles_in)})")

    df_reads = pd.DataFrame()
    df_cells = pd.DataFrame()
    if not df_bases.empty:
        df_reads = call_reads(
            df_bases,
            peaks_data=peaks,
            correction_only_in_cells=True,
            normalize_bases_first=True,
            method=p["call_reads_method"],
        )
        n_called = len(df_reads)
        print(f"9. Call reads      reads={n_called}")

        df_cells = call_cells(
            df_reads,
            df_barcode_library=DF_BARCODE_LIBRARY,
            prefix_col=p["prefix_col"],
            sort_calls=p["sort_calls"],
            error_correct=p["error_correct"],
        )
        n_mapped = len(df_cells)
        print(f"10. Call cells     cells={n_mapped}")

    # Save DataFrames under threshold dir
    t_dir = pipeline_dir / f"threshold_{p['threshold_peaks']}"
    t_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{well}_pos{position}"
    df_bases.to_csv(t_dir / f"{tag}_df_bases.tsv", sep="\t", index=False)
    df_reads.to_csv(t_dir / f"{tag}_df_reads.tsv", sep="\t", index=False)
    df_cells.to_csv(t_dir / f"{tag}_df_cells.tsv", sep="\t", index=False)
    vprint(f"   Saved DataFrames to {t_dir}")

    return {
        "aligned": aligned,
        "nuclei": nuclei,
        "cells": cells,
        "peaks": peaks,
        "maxed": maxed,
        "df_bases": df_bases,
        "df_reads": df_reads,
        "df_cells": df_cells,
        "n_cells": n_cells,
    }


# ---------------------------------------------------------------------------
# Rerun from cached arrays
# ---------------------------------------------------------------------------
def rerun_from_cache(well: str, position: int, threshold: int) -> dict:
    """Re-run extract_bases → call_reads → call_cells from cached arrays."""
    p = PARAMS
    pipeline_dir = DATA_CONFIG["data_root"] / "pipeline"
    pos_dir = pipeline_dir / f"{well}_pos{position}"

    aligned = tifffile.imread(str(pos_dir / "aligned.tiff"))
    nuclei = tifffile.imread(str(pos_dir / "nuclei.tiff"))
    peaks = tifffile.imread(str(pos_dir / "peaks.tiff"))
    maxed = tifffile.imread(str(pos_dir / "maxed.tiff"))

    n_cells = len(np.unique(nuclei)) - 1
    wildcards = {"well": well, "tile": position}

    print(f"  Loaded cached arrays from {pos_dir}")
    print(f"  nuclei={n_cells}  threshold={threshold}")

    df_bases = extract_bases(
        peaks_data=peaks,
        max_filtered_data=maxed,
        cells_data=nuclei,
        threshold_peaks=threshold,
        wildcards=wildcards,
        bases="".join(p["bases"]),
    )
    n_reads = df_bases["read"].nunique() if not df_bases.empty else 0
    print(f"  Extract bases: {n_reads} reads")

    df_reads = pd.DataFrame()
    df_cells = pd.DataFrame()
    if not df_bases.empty:
        df_reads = call_reads(
            df_bases,
            peaks_data=peaks,
            correction_only_in_cells=True,
            normalize_bases_first=True,
            method=p["call_reads_method"],
        )
        print(f"  Call reads: {len(df_reads)} reads")

        df_cells = call_cells(
            df_reads,
            df_barcode_library=DF_BARCODE_LIBRARY,
            prefix_col=p["prefix_col"],
            sort_calls=p["sort_calls"],
            error_correct=p["error_correct"],
        )
        print(f"  Call cells: {len(df_cells)} cells")

    # Save under new threshold dir
    t_dir = pipeline_dir / f"threshold_{threshold}"
    t_dir.mkdir(parents=True, exist_ok=True)
    tag = f"{well}_pos{position}"
    df_bases.to_csv(t_dir / f"{tag}_df_bases.tsv", sep="\t", index=False)
    df_reads.to_csv(t_dir / f"{tag}_df_reads.tsv", sep="\t", index=False)
    df_cells.to_csv(t_dir / f"{tag}_df_cells.tsv", sep="\t", index=False)

    return {
        "aligned": aligned,
        "nuclei": nuclei,
        "cells": nuclei,
        "peaks": peaks,
        "maxed": maxed,
        "df_bases": df_bases,
        "df_reads": df_reads,
        "df_cells": df_cells,
        "n_cells": n_cells,
    }


# ---------------------------------------------------------------------------
# Consensus base call
# ---------------------------------------------------------------------------
def consensus_base_call(
    peaks: np.ndarray,
    maxed: np.ndarray,
    threshold: float,
) -> dict:
    """Call a consensus barcode by aggregating signal across all peaks."""
    bases = PARAMS["bases"]
    yi, xi = np.where(peaks > threshold)
    n_peaks = len(yi)
    if n_peaks == 0:
        return {"barcode": "", "n_peaks": 0, "cycle_table": None}

    n_cycles, n_bases = maxed.shape[0], maxed.shape[1]
    spot_intensities = maxed[:, :, yi, xi]
    cycle_sums = spot_intensities.sum(axis=2).astype(float)
    per_peak_calls = np.argmax(spot_intensities, axis=1)

    rows = []
    barcode_chars = []
    for c in range(n_cycles):
        sums = cycle_sums[c]
        total = sums.sum()
        fracs = sums / total if total > 0 else sums
        winner = int(np.argmax(sums))
        base = bases[winner]
        barcode_chars.append(base)
        agreement = (per_peak_calls[c] == winner).sum() / n_peaks
        rows.append(
            {
                "cycle": c + 1,
                "base": base,
                **{f"{bases[b]}_frac": f"{fracs[b]:.3f}" for b in range(n_bases)},
                "agreement": f"{agreement:.3f}",
            }
        )

    return {
        "barcode": "".join(barcode_chars),
        "n_peaks": n_peaks,
        "cycle_table": pd.DataFrame(rows),
    }


# ---------------------------------------------------------------------------
# Raw barcode specificity
# ---------------------------------------------------------------------------
def raw_barcode_specificity(
    peaks: np.ndarray,
    maxed: np.ndarray,
    threshold: float,
    known_prefixes: list[str],
) -> dict:
    """Per-spot base calling specificity against known barcodes (no normalization)."""
    bases = PARAMS["bases"]
    yi, xi = np.where(peaks > threshold)
    n_peaks = len(yi)
    if n_peaks == 0:
        return {"n_reads": 0}

    n_cycles = maxed.shape[0]
    spot_intensities = maxed[:, :, yi, xi]
    raw_calls = np.argmax(spot_intensities, axis=1)

    known_set = set(known_prefixes)
    per_read_barcodes = []
    per_read_known = []
    for p_idx in range(n_peaks):
        bc = "".join(bases[raw_calls[c, p_idx]] for c in range(n_cycles))
        per_read_barcodes.append(bc)
        per_read_known.append(bc in known_set)

    bc_series = pd.Series(per_read_barcodes)
    bc_counts = bc_series.value_counts()

    per_barcode = {}
    for prefix in known_prefixes:
        count = bc_counts.get(prefix, 0)
        per_barcode[prefix] = {"count": int(count), "fraction": count / n_peaks}

    n_known = sum(per_read_known)
    unknown_bcs = {
        bc: int(cnt) for bc, cnt in bc_counts.items() if bc not in known_set
    }

    per_cycle_dist = []
    for c in range(n_cycles):
        counts = np.bincount(raw_calls[c], minlength=len(bases))
        per_cycle_dist.append({bases[b]: int(counts[b]) for b in range(len(bases))})

    return {
        "n_reads": n_peaks,
        "n_known": n_known,
        "n_unknown": n_peaks - n_known,
        "known_fraction": n_known / n_peaks,
        "per_barcode": per_barcode,
        "unknown_barcodes": unknown_bcs,
        "top_unknown": dict(sorted(unknown_bcs.items(), key=lambda x: -x[1])[:10]),
        "per_read_barcodes": per_read_barcodes,
        "per_cycle_distribution": per_cycle_dist,
    }


# ---------------------------------------------------------------------------
# Cross-well read calling
# ---------------------------------------------------------------------------
def call_reads_cross_well(
    per_well_results: dict[str, dict],
) -> dict[str, dict]:
    """Combine df_bases from all wells, call reads + call cells, split back."""
    p = PARAMS
    combined_parts = []
    read_offset = 0

    for well, res in per_well_results.items():
        df_b = res["df_bases"].copy()
        if df_b.empty:
            continue
        max_read = df_b["read"].max() + 1
        df_b["read"] = df_b["read"] + read_offset
        read_offset += max_read
        combined_parts.append(df_b)

    if not combined_parts:
        empty_reads = pd.DataFrame(
            columns=["read", "cell", "i", "j", "tile", "well",
                     "barcode", "Q_min", "peak"]
        )
        empty_cells = pd.DataFrame()
        return {w: {"df_reads": empty_reads.copy(), "df_cells": empty_cells.copy()}
                for w in per_well_results}

    df_combined = pd.concat(combined_parts, ignore_index=True)
    n_reads = df_combined["read"].nunique()
    print(f"\n   Cross-well mix: {len(per_well_results)} wells, {n_reads} total reads")

    max_i = int(df_combined["i"].max()) + 1
    max_j = int(df_combined["j"].max()) + 1
    dummy_peaks = np.zeros((max_i, max_j), dtype=np.uint16)

    df_reads_all = call_reads(
        df_combined,
        peaks_data=dummy_peaks,
        correction_only_in_cells=True,
        normalize_bases_first=True,
        method=p["call_reads_method"],
    )

    result = {}
    for well in per_well_results:
        mask = df_reads_all["well"] == well
        df_reads_well = df_reads_all[mask].copy()
        n = len(df_reads_well)

        df_cells_well = pd.DataFrame()
        if n > 0:
            top_bc = df_reads_well["barcode"].value_counts().index[0]
            top_ct = df_reads_well["barcode"].value_counts().iloc[0]
            print(f"   {well}: {n} reads, top={top_bc} ({top_ct}/{n}, {top_ct/n*100:.1f}%)")

            df_cells_well = call_cells(
                df_reads_well,
                df_barcode_library=DF_BARCODE_LIBRARY,
                prefix_col=p["prefix_col"],
                sort_calls=p["sort_calls"],
                error_correct=p["error_correct"],
            )
            print(f"   {well}: {len(df_cells_well)} cells called")

        result[well] = {"df_reads": df_reads_well, "df_cells": df_cells_well}

    return result


# ---------------------------------------------------------------------------
# Build summary table
# ---------------------------------------------------------------------------
def build_summary_table(
    well_data: dict[str, dict],
    known_prefixes: list[str],
) -> pd.DataFrame:
    """Build per-well summary table with all key metrics."""
    construct_names = list(DF_BARCODE_LIBRARY["gene_symbol"])
    rows = []

    for well in sorted(well_data.keys()):
        wd = well_data[well]
        det = wd.get("detection", {})
        n_seg = wd.get("n_cells_total", 0)
        n_mapped = det.get("n_mapped", 0)
        n_unique = det.get("n_uniquely_mapped", 0)
        per_con = det.get("per_construct", {})

        rows.append({
            "well": well,
            "n_positions": len(wd.get("raw_specificity", [])),
            "n_segmented": n_seg,
            "n_mapped": n_mapped,
            "n_uniquely_mapped": n_unique,
            "pct_mapped": n_mapped / n_seg * 100 if n_seg > 0 else 0,
            "pct_uniquely_mapped": n_unique / n_seg * 100 if n_seg > 0 else 0,
            **{name: per_con.get(name, 0) for name in construct_names},
        })

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def generate_plots(summary_df: pd.DataFrame, out_dir: Path) -> None:
    """Generate publication-quality plots from summary table."""
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(BENCHMARKS_DIR))
    from plot_style import setup_plot_style, FIGSIZE, save_figure

    setup_plot_style()
    out_dir.mkdir(parents=True, exist_ok=True)

    wells = summary_df["well"].values
    x = np.arange(len(wells))
    bw = 0.5

    # Three categories per well (stacking to n_segmented)
    n_unique = summary_df["n_uniquely_mapped"].values.astype(float)
    n_ambig = (summary_df["n_mapped"] - summary_df["n_uniquely_mapped"]).values.astype(float)
    n_no_bc = (summary_df["n_segmented"] - summary_df["n_mapped"]).values.astype(float)
    n_seg = summary_df["n_segmented"].values.astype(float)

    c_unique = "#27ae60"   # green
    c_ambig = "#f39c12"    # orange
    c_nobc = "#bdc3c7"     # gray

    labels = ["Uniquely mapped", "Ambiguous (multi-construct)", "No barcode"]

    # ── Plot 1: Counts (no percentages, just counts) ──
    fig, ax = plt.subplots(figsize=FIGSIZE["single"])
    ax.bar(x, n_unique, bw, label=labels[0], color=c_unique,
           edgecolor="white", linewidth=0.5)
    ax.bar(x, n_ambig, bw, bottom=n_unique, label=labels[1],
           color=c_ambig, edgecolor="white", linewidth=0.5)
    ax.bar(x, n_no_bc, bw, bottom=n_unique + n_ambig, label=labels[2],
           color=c_nobc, edgecolor="white", linewidth=0.5)

    for i in range(len(wells)):
        # Counts only — no percentages
        if n_unique[i] / n_seg[i] > 0.04:
            ax.text(x[i], n_unique[i] / 2, f"{int(n_unique[i]):,}",
                    ha="center", va="center", fontsize=9, fontweight="bold", color="white")
        if n_ambig[i] / n_seg[i] > 0.04:
            ax.text(x[i], n_unique[i] + n_ambig[i] / 2, f"{int(n_ambig[i]):,}",
                    ha="center", va="center", fontsize=9, fontweight="bold", color="white")
        if n_no_bc[i] / n_seg[i] > 0.04:
            ax.text(x[i], n_unique[i] + n_ambig[i] + n_no_bc[i] / 2, f"{int(n_no_bc[i]):,}",
                    ha="center", va="center", fontsize=9, fontweight="bold", color="#555")

    ax.set_ylabel("Cells")
    ax.set_title("Cell Mapping (counts)")
    ax.set_xticks(x)
    ax.set_xticklabels(wells)
    ax.set_xlabel("Well")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
    fig.tight_layout()
    save_figure(fig, out_dir / "cell_mapping_counts")
    plt.close(fig)
    print(f"Saved: {out_dir / 'cell_mapping_counts.png'}")

    # ── Plot 2: Percentages ──
    pct_unique = n_unique / n_seg * 100
    pct_ambig = n_ambig / n_seg * 100
    pct_nobc = n_no_bc / n_seg * 100

    fig, ax = plt.subplots(figsize=FIGSIZE["single"])
    ax.bar(x, pct_unique, bw, label=labels[0], color=c_unique,
           edgecolor="white", linewidth=0.5)
    ax.bar(x, pct_ambig, bw, bottom=pct_unique, label=labels[1],
           color=c_ambig, edgecolor="white", linewidth=0.5)
    ax.bar(x, pct_nobc, bw, bottom=pct_unique + pct_ambig, label=labels[2],
           color=c_nobc, edgecolor="white", linewidth=0.5)

    for i in range(len(wells)):
        if pct_unique[i] > 4:
            ax.text(x[i], pct_unique[i] / 2, f"{pct_unique[i]:.1f}%",
                    ha="center", va="center", fontsize=11, fontweight="bold", color="white")
        if pct_ambig[i] > 4:
            ax.text(x[i], pct_unique[i] + pct_ambig[i] / 2, f"{pct_ambig[i]:.1f}%",
                    ha="center", va="center", fontsize=10, fontweight="bold", color="white")
        if pct_nobc[i] > 4:
            ax.text(x[i], pct_unique[i] + pct_ambig[i] + pct_nobc[i] / 2,
                    f"{pct_nobc[i]:.1f}%",
                    ha="center", va="center", fontsize=10, fontweight="bold", color="#555")

    ax.set_ylabel("% of segmented cells")
    ax.set_title("Cell Mapping (percent)")
    ax.set_xticks(x)
    ax.set_xticklabels(wells)
    ax.set_xlabel("Well")
    ax.set_ylim(0, 105)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, frameon=False)
    fig.tight_layout()
    save_figure(fig, out_dir / "cell_mapping_pct")
    plt.close(fig)
    print(f"Saved: {out_dir / 'cell_mapping_pct.png'}")


# ---------------------------------------------------------------------------
# Print summary
# ---------------------------------------------------------------------------
def print_summary(
    well_data: dict[str, dict],
    mode: str = "standard",
) -> None:
    """Print comprehensive summary across all wells."""
    sep = "=" * 72
    known = DATA_CONFIG["known_prefixes"]

    print(f"\n{sep}")
    print(f"ISS MULTI-CONSTRUCT VALIDATION SUMMARY — {mode.upper()}")
    print(sep)

    print(f"  channels:        {PARAMS['channel_names']}")
    print(f"  bases:           {PARAMS['bases']}")
    print(f"  known barcodes:  {known}")
    print(f"  threshold:       {PARAMS['threshold_peaks']}")
    print(f"  random chance:   {len(known)}/64 = {len(known)/64*100:.1f}%")

    # Raw specificity
    print(f"\n{'─' * 72}")
    print("RAW BASE CALLING SPECIFICITY (argmax, no normalization)")
    print(f"{'─' * 72}")
    print(
        f"{'Well':<6} {'Pos':<5} {'Reads':>7} "
        + "".join(f"  {bc:>7}" for bc in known)
        + f"  {'Known%':>7}  {'Unknown%':>8}"
    )
    print("-" * (6 + 5 + 8 + 9 * len(known) + 9 + 9))

    total_raw_reads = 0
    total_raw_known = 0
    for well in sorted(well_data.keys()):
        for spec in well_data[well].get("raw_specificity", []):
            if spec.get("n_reads", 0) == 0:
                continue
            pos = spec.get("position", "?")
            n = spec["n_reads"]
            total_raw_reads += n
            total_raw_known += spec["n_known"]
            parts = [f"{well:<6} {pos:<5} {n:>7}"]
            for bc in known:
                info = spec["per_barcode"].get(bc, {"fraction": 0})
                parts.append(f"  {info['fraction']*100:>6.1f}%")
            parts.append(f"  {spec['known_fraction']*100:>6.1f}%")
            parts.append(f"  {(1-spec['known_fraction'])*100:>7.1f}%")
            print("".join(parts))

    if total_raw_reads > 0:
        print(
            f"\n  Raw total: {total_raw_reads} reads, "
            f"{total_raw_known} known ({total_raw_known/total_raw_reads*100:.1f}%), "
            f"{total_raw_reads-total_raw_known} unknown "
            f"({(total_raw_reads-total_raw_known)/total_raw_reads*100:.1f}%)"
        )

    # Cell-level mapping (from call_cells)
    print(f"\n{'─' * 72}")
    print("CELL-LEVEL MAPPING (after call_cells)")
    print(f"{'─' * 72}")

    for well in sorted(well_data.keys()):
        det = well_data[well].get("detection", {})
        n_seg = well_data[well].get("n_cells_total", 0)
        n_mapped = det.get("n_mapped", 0)
        n_unique = det.get("n_uniquely_mapped", 0)
        if n_seg == 0:
            print(f"  {well}: no cells")
            continue
        n_ambig = n_mapped - n_unique
        n_no_bc = n_seg - n_mapped
        print(f"\n  {well}: {n_seg} segmented")
        print(f"    Uniquely mapped: {n_unique:>6} ({n_unique/n_seg*100:.1f}%)")
        print(f"    Ambiguous:       {n_ambig:>6} ({n_ambig/n_seg*100:.1f}%)")
        print(f"    No barcode:      {n_no_bc:>6} ({n_no_bc/n_seg*100:.1f}%)")
        per_con = det.get("per_construct", {})
        if per_con:
            print(f"    Per-construct (uniquely mapped):")
            for name, count in per_con.items():
                print(f"      {name}: {count:>6} ({count/n_unique*100:.1f}%)" if n_unique > 0 else f"      {name}: {count:>6}")

    # Consensus
    print(f"\n{'─' * 72}")
    print("CONSENSUS BASE CALLING (aggregate signal — informational)")
    print(f"{'─' * 72}")
    for well in sorted(well_data.keys()):
        for cons in well_data[well].get("consensus", []):
            bc = cons["barcode"]
            pos = cons.get("position", "?")
            n_peaks = cons["n_peaks"]
            if not bc:
                continue
            match_str = "KNOWN" if bc in known else "not in known set"
            print(f"  {well} pos {pos}: {bc}  ({n_peaks} peaks)  [{match_str}]")
            if _VERBOSE and cons["cycle_table"] is not None:
                print(cons["cycle_table"].to_string(index=False))

    # Overall
    print(f"\n{sep}")
    print("OVERALL")
    total_seg = sum(well_data[w].get("n_cells_total", 0) for w in well_data)
    total_mapped = sum(well_data[w].get("detection", {}).get("n_mapped", 0) for w in well_data)
    total_unique = sum(well_data[w].get("detection", {}).get("n_uniquely_mapped", 0) for w in well_data)
    if total_seg > 0:
        print(f"  Segmented:         {total_seg}")
        print(f"  Uniquely mapped:   {total_unique} ({total_unique/total_seg*100:.1f}%)")
        print(f"  Ambiguous:         {total_mapped-total_unique} ({(total_mapped-total_unique)/total_seg*100:.1f}%)")
        print(f"  No barcode:        {total_seg-total_mapped} ({(total_seg-total_mapped)/total_seg*100:.1f}%)")
    if total_raw_reads > 0:
        print(f"  Raw known:         {total_raw_known} ({total_raw_known/total_raw_reads*100:.1f}%)")
    print(sep)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    global _VERBOSE

    parser = argparse.ArgumentParser(
        description="ISS Multi-Construct Validation Benchmark"
    )
    parser.add_argument("--well", type=str, default="B03")
    pos_group = parser.add_mutually_exclusive_group()
    pos_group.add_argument("--position", type=int, default=None)
    pos_group.add_argument("--positions", type=str, default=None)
    parser.add_argument("--all-wells", action="store_true")
    parser.add_argument("--cross-well", action="store_true")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--skip-save", action="store_true")
    parser.add_argument("--threshold", type=int, default=None)
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Regenerate plots from existing summary table (no pipeline rerun)",
    )
    parser.add_argument(
        "--rerun",
        action="store_true",
        help="Re-run from cached arrays (skip alignment/segmentation)",
    )
    args = parser.parse_args()

    _VERBOSE = args.verbose
    if args.threshold is not None:
        PARAMS["threshold_peaks"] = args.threshold

    out_dir = DATA_CONFIG["output_dir"]
    known = DATA_CONFIG["known_prefixes"]
    construct_names = list(DF_BARCODE_LIBRARY["gene_symbol"])

    # ── Plots-only mode: load saved summary and regenerate plots ──
    if args.plots_only:
        summary_path = out_dir / _SUMMARY_FILENAME
        if not summary_path.exists():
            print(f"ERROR: Summary table not found: {summary_path}")
            print("Run the full pipeline first to generate data.")
            sys.exit(1)
        summary_df = pd.read_csv(summary_path, sep="\t")
        print(f"Loaded summary table: {summary_path}")
        print(summary_df.to_string(index=False))
        generate_plots(summary_df, out_dir)
        return

    # ── Parse positions and wells ──
    if args.positions is not None:
        if args.positions == "all":
            positions = list(range(DATA_CONFIG["n_positions"]))
        else:
            positions = [int(p) for p in args.positions.split(",")]
    elif args.position is not None:
        positions = [args.position]
    else:
        positions = [0]

    if args.all_wells or args.cross_well:
        wells = DATA_CONFIG["all_wells"]
    else:
        wells = [args.well]

    mode_str = "rerun" if args.rerun else ("cross-well" if args.cross_well else "standard")
    print("=" * 72)
    print("ISS MULTI-CONSTRUCT VALIDATION")
    print("=" * 72)
    print(f"  wells:       {wells}")
    print(f"  positions:   {positions}")
    print(f"  known bc:    {known}  (from 12-mers, 3 cycles → 3 bases)")
    print(f"  mode:        {mode_str}")
    print(f"  data:        {DATA_CONFIG['data_root']}")
    print(f"  threshold:   {PARAMS['threshold_peaks']}")

    # ── Run pipeline (full or rerun from cache) ──
    all_results: dict[str, dict[int, dict]] = {w: {} for w in wells}
    all_consensus: dict[str, list[dict]] = {w: [] for w in wells}

    for pos in positions:
        for well in wells:
            print(f"\n{'─' * 72}")
            print(f"WELL {well}  POSITION {pos}")
            print(f"{'─' * 72}")

            if args.rerun:
                results = rerun_from_cache(well, pos, PARAMS["threshold_peaks"])
            else:
                images = load_position_across_cycles(well, pos)
                results = run_sbs_pipeline(images, well=well, position=pos)
            all_results[well][pos] = results

            cons = consensus_base_call(
                results["peaks"], results["maxed"], PARAMS["threshold_peaks"]
            )
            cons["position"] = pos
            cons["well"] = well
            all_consensus[well].append(cons)
            if cons["barcode"]:
                match_str = "KNOWN" if cons["barcode"] in known else "unknown"
                print(f"\n   CONSENSUS: {cons['barcode']}  ({cons['n_peaks']} peaks)  [{match_str}]")

    # ── Raw specificity ──
    well_data: dict[str, dict] = {
        w: {"consensus": all_consensus[w], "raw_specificity": []}
        for w in wells
    }

    for pos in positions:
        for well in wells:
            res = all_results[well][pos]
            spec = raw_barcode_specificity(
                res["peaks"], res["maxed"], PARAMS["threshold_peaks"], known
            )
            spec["position"] = pos
            well_data[well]["raw_specificity"].append(spec)

    # ── Aggregate cells per well ──
    for well in wells:
        well_data[well]["n_cells_total"] = sum(
            all_results[well][pos]["n_cells"] for pos in positions
        )

    # ── Collect call_cells results ──
    if args.cross_well:
        for pos in positions:
            tile_results = {w: all_results[w][pos] for w in wells}
            per_well = call_reads_cross_well(tile_results)
            for well in wells:
                cw = per_well[well]
                well_data[well].setdefault("all_df_cells", []).append(cw["df_cells"])
                well_data[well]["n_reads_total"] = (
                    well_data[well].get("n_reads_total", 0) + len(cw["df_reads"])
                )
    else:
        for pos in positions:
            for well in wells:
                res = all_results[well][pos]
                well_data[well].setdefault("all_df_cells", []).append(res["df_cells"])
                well_data[well]["n_reads_total"] = (
                    well_data[well].get("n_reads_total", 0) + len(res["df_reads"])
                )

    # ── Cell-level detection analysis (from call_cells) ──
    for well in wells:
        cell_dfs = [df for df in well_data[well].get("all_df_cells", []) if not df.empty]
        if not cell_dfs:
            well_data[well]["detection"] = {
                "n_cells_called": 0, "n_mapped": 0,
                "n_uniquely_mapped": 0, "per_construct": {},
            }
            continue

        df_all_cells = pd.concat(cell_dfs, ignore_index=True)
        n_called = len(df_all_cells)

        # Mapped = gene_symbol_0 is not NaN (top barcode hit a known construct)
        mapped = df_all_cells["gene_symbol_0"].notna()
        n_mapped = int(mapped.sum())

        # Uniquely mapped = mapped AND (no second barcode OR second agrees)
        gs0 = df_all_cells["gene_symbol_0"]
        gs1 = df_all_cells["gene_symbol_1"]
        uniquely_mapped = mapped & (gs1.isna() | (gs0 == gs1))
        n_uniquely_mapped = int(uniquely_mapped.sum())

        # Per-construct counts (uniquely mapped only)
        gene_counts = df_all_cells.loc[uniquely_mapped, "gene_symbol_0"].value_counts()
        per_construct = {}
        for name in construct_names:
            per_construct[name] = int(gene_counts.get(name, 0))

        well_data[well]["detection"] = {
            "n_cells_called": n_called,
            "n_mapped": n_mapped,
            "n_uniquely_mapped": n_uniquely_mapped,
            "per_construct": per_construct,
        }

    # ── Print summary ──
    mode = "cross-well" if args.cross_well else ("rerun" if args.rerun else "standard")
    print_summary(well_data, mode=mode)

    # ── Save summary table and plots ──
    if not args.skip_save:
        out_dir.mkdir(parents=True, exist_ok=True)

        summary_df = build_summary_table(well_data, known)
        summary_path = out_dir / _SUMMARY_FILENAME
        summary_df.to_csv(summary_path, sep="\t", index=False)
        print(f"\nSaved summary table: {summary_path}")
        print(summary_df.to_string(index=False))

        generate_plots(summary_df, out_dir)


if __name__ == "__main__":
    main()
