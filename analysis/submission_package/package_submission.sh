#!/bin/bash
# Packages submission figures and source data into a single zip.
# Run from anywhere — uses absolute paths derived from script location.
#
# Output: analysis/submission_package/aconcagua_submission.zip
#
# Figures that exist only on Mac (schematics, BioRender panels) are reported
# as "not found on HPC" — this is expected.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
FIGURES_DIR="$SCRIPT_DIR/figures"
EXCEL_DIR="$SCRIPT_DIR/excel"
ZIP_NAME="$SCRIPT_DIR/aconcagua_submission.zip"

declare -A FIGURE_LABEL
FIGURE_LABEL["P-1__cell_mapping_heatmap_one"]="2b"
FIGURE_LABEL["P-1__gene_symbol_histogram"]="2b"
FIGURE_LABEL["spotiflow_vs_standard_comparison_square"]="2c"
FIGURE_LABEL["P-1__cell_PHALLOIDIN_min_heatmap"]="2f"
FIGURE_LABEL["P-1__cell_density_heatmap"]="2f"
FIGURE_LABEL["P-3__sbs_to_ph_matching_rates"]="3b"
FIGURE_LABEL["P-3__cells_with_channel_min_0"]="3b"
FIGURE_LABEL["base_legend"]="3b"
FIGURE_LABEL["CeCl-all_ChCo-DAPI_TUBULIN_GH2AX_PHALLOIDIN__feature_distributions"]="3d"
FIGURE_LABEL["phate_scatter_interphase_mitochondrial"]="5a"
FIGURE_LABEL["pearson_correlation_interphase"]="5a"
FIGURE_LABEL["heatmap_interphase_mitochondrial_pca_brieflow"]="5c"
FIGURE_LABEL["heatmap_interphase_mitochondrial_legend"]="5c"
FIGURE_LABEL["heatmap_interphase_mitochondrial_raw_funk"]="5d"
FIGURE_LABEL["brieflow_rulegraph"]="s1"
FIGURE_LABEL["runtime_memory_comparison"]="s3ab"
FIGURE_LABEL["cell_mapping_pct"]="s3d"
FIGURE_LABEL["count_comparison"]="s4ad"
FIGURE_LABEL["runtime_comparison"]="s4e"
FIGURE_LABEL["cluster_scatter_hu_moments"]="s4fh"
FIGURE_LABEL["cluster_scatter_intensity"]="s4fh"
FIGURE_LABEL["cluster_scatter_radial_distribution"]="s4fh"
FIGURE_LABEL["tiles_phenotype_fast_only"]="s5bc"
FIGURE_LABEL["classifier"]="s6d"
FIGURE_LABEL["confusion_matrix"]="s6c"
FIGURE_LABEL["confidence"]="s6d"
FIGURE_LABEL["ideal_comparison"]="s7a"
FIGURE_LABEL["figure_a_high_confidence"]="s7b"
FIGURE_LABEL["figure_c_mean_genes"]="s7c"
FIGURE_LABEL["cluster_preservation_text"]="s8a"
FIGURE_LABEL["crosstalk"]="s8a"
FIGURE_LABEL["redistribution_sankey"]="s8b"
FIGURE_LABEL["novel"]="s9a"
FIGURE_LABEL["cluster_145"]="s9a"
FIGURE_LABEL["brieflow_schema"]="1"
FIGURE_LABEL["montage"]="2e"

TARGET_STEMS=(
    "P-1__cell_mapping_heatmap_one"
    "P-1__gene_symbol_histogram"
    "spotiflow_vs_standard_comparison_square"
    "P-1__cell_PHALLOIDIN_min_heatmap"
    "P-1__cell_density_heatmap"
    "P-3__sbs_to_ph_matching_rates"
    "P-3__cells_with_channel_min_0"
    "base_legend"
    "CeCl-all_ChCo-DAPI_TUBULIN_GH2AX_PHALLOIDIN__feature_distributions"
    "phate_scatter_interphase_mitochondrial"
    "pearson_correlation_interphase"
    "heatmap_interphase_mitochondrial_pca_brieflow"
    "heatmap_interphase_mitochondrial_legend"
    "heatmap_interphase_mitochondrial_raw_funk"
    "brieflow_rulegraph"
    "runtime_memory_comparison"
    "cell_mapping_pct"
    "count_comparison"
    "runtime_comparison"
    "cluster_scatter_hu_moments"
    "cluster_scatter_intensity"
    "cluster_scatter_radial_distribution"
    "tiles_phenotype_fast_only"
    "classifier"
    "confusion_matrix"
    "confidence"
    "ideal_comparison"
    "figure_a_high_confidence"
    "figure_c_mean_genes"
    "cluster_preservation_text"
    "crosstalk"
    "redistribution_sankey"
    "novel"
    "cluster_145"
    "brieflow_schema"
    "montage"
)

rm -rf "$FIGURES_DIR" "$EXCEL_DIR" "$ZIP_NAME"
mkdir -p "$FIGURES_DIR" "$EXCEL_DIR"

found=0
missing=()

for stem in "${TARGET_STEMS[@]}"; do
    match=$(find "$REPO_ROOT/benchmarks/results" \
        -name "${stem}.pdf" 2>/dev/null | head -1)
    if [ -n "$match" ]; then
        label="${FIGURE_LABEL[$stem]:-unknown}"
        dest="${label}__${stem}.pdf"
        cp "$match" "$FIGURES_DIR/${dest}"
        echo "  OK  ${dest}"
        ((found++)) || true
    else
        missing+=("$stem")
    fi
done

echo ""
cp "$REPO_ROOT/benchmarks/results/source_data.xlsx" "$EXCEL_DIR/" \
    2>/dev/null && echo "  OK  source_data.xlsx" || echo "WARN: source_data.xlsx not found"
for f in "$REPO_ROOT/benchmarks/results/cluster/tables/Supplementary_Data_"*.xlsx; do
    [ -f "$f" ] && cp "$f" "$EXCEL_DIR/" && echo "  OK  $(basename "$f")" || true
done

echo ""
echo "=== Summary ==="
echo "PDFs found:   $found / ${#TARGET_STEMS[@]}"
if [ ${#missing[@]} -gt 0 ]; then
    echo "Not on HPC (Mac-only or missing):"
    for m in "${missing[@]}"; do echo "    $m"; done
fi
echo "Excel files:  $(ls "$EXCEL_DIR/" | wc -l)"

zip -r "$ZIP_NAME" "$FIGURES_DIR/" "$EXCEL_DIR/"
echo ""
echo "Done: $ZIP_NAME ($(du -sh "$ZIP_NAME" | cut -f1))"
