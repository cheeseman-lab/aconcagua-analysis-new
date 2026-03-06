# Benchmarks

Compares methods at four pipeline stages, plus clustering evaluation and supplemental table generation.

## Segmentation (`segmentation.py`)

**Methods:** Cellpose (cyto3), Cellpose4 (CPSAM), StarDist
**Metrics:** Runtime, memory, cell/nuclei counts, retention rates

## Spot Calling (`spot_calling.py`)

**Methods:** Standard (peak detection), Spotiflow (deep learning)
**Metrics:** Runtime, memory, peaks, read mapping, cell assignment

## Feature Extraction (`feature_extraction.py`)

Two modes:
- **comparison** (default): Compare CP Measure vs CP Emulator pipelines
- **timing**: Time individual cp_measure functions to identify bottlenecks

**Metrics:** Runtime, memory, feature counts

## Merge (`merge.py`)

**Methods:** Fast (tile-by-tile alignment), Stitch (full well stitching)
**Metrics:** Runtime, alignment quality (distance), cell retention, barcode mapping

**Conclusion:** Fast approach is recommended. Stitch requires ~2TB RAM for triangle hash distance matrix (beyond cluster capacity), falls back to identity transform with degraded alignment.

## Clustering (`cluster_enrichment.py`)

**Comparison:** Brieflow vs Funk et al. 2022 original analysis
**Metrics:** STRING F1, CORUM/KEGG enrichment, Precision-Recall curves

### Supporting cluster analysis scripts

| Script | Purpose | Output Dir |
|--------|---------|------------|
| `cluster_overlap.py` | Jaccard similarity between Brieflow↔Funk clusters | `results/cluster/overlap/` |
| `cluster_validation.py` | Co-clustering preservation of known gene groups | `results/cluster/validation/` |
| `cluster_similarity.py` | Feature embedding Pearson correlation heatmaps | `results/cluster/similarity/` |

### MozzareLLM

| Script | Purpose | Output Dir |
|--------|---------|------------|
| `cluster_mozzarellm.py` | Publication figures A-C from LLM cluster summaries | `results/cluster/mozzarellm/` |
| `run_mozzarellm_funk.py` | Run MozzareLLM on Funk re-clustered data | `external/results/cluster/*/mozzarellm/` |

### Supplemental Tables (`cluster_tables.py`)

Generates all supplemental tables (S1-S9) and a merged Excel workbook:
- **S1:** All high-confidence clusters
- **S2:** Gene annotations for high-confidence clusters
- **S3:** Cross-pipeline Jaccard similarity (requires `cluster_overlap.py` output)
- **S4:** Funk cluster retention in Brieflow
- **S5:** Shuffled control comparison
- **S6-S9:** Comprehensive cluster and gene tables (Brieflow + Funk)

**Output:** `results/cluster/tables/` including `all_supplemental_tables.xlsx`

## Classifier (`classifier.py`)

XGBoost cell stage classifier for validating morphological feature quality.
- Trains on Interphase vs Mitotic labels using extracted features
- Generates confusion matrix and training composition plots

**Metrics:** Accuracy, precision, recall, F1, feature importance

## ISS Multi (`iss_multi.py`)

Multi-construct barcode validation for in-situ sequencing quality control.
- Validates barcode detection across wells and positions
- Analyzes cell-to-barcode mapping rates

**Metrics:** Barcode counts, cell mapping rates, per-well/position breakdowns

## Visualization (`visualization.py`)

Generate publication-quality visualizations of benchmark results. Requires Napari viewer — run locally, not via `run_benchmarks.sh`.

## Usage

```bash
# Run all benchmarks
bash run_benchmarks.sh

# Run specific benchmark type
bash run_benchmarks.sh --segmentation
bash run_benchmarks.sh --spot-calling
bash run_benchmarks.sh --feature-extraction
bash run_benchmarks.sh --feature-timing
bash run_benchmarks.sh --merge
bash run_benchmarks.sh --cluster
bash run_benchmarks.sh --classifier

# Run single method/mode directly
python segmentation.py --method cellpose
python feature_extraction.py --mode timing
python merge.py                # Analyze merge results from pipeline
python merge.py --plots-only   # Regenerate plots only
python cluster_enrichment.py              # Compare clustering with Funk et al. 2022
python cluster_enrichment.py --plots-only # Regenerate plots only
python classifier.py           # Train and evaluate cell stage classifier

# Cluster analysis scripts (run independently)
python cluster_overlap.py
python cluster_validation.py
python cluster_similarity.py
python cluster_mozzarellm.py --include-shuffled
python cluster_tables.py       # Generate S1-S9 tables + Excel

# Generate visualizations (requires Napari, run locally)
python visualization.py --type all
```

## Environments

| Environment | Used For |
|-------------|----------|
| `brieflow_aconcagua` | Main: Cellpose3, spot calling, feature extraction, merge, clustering |
| `brieflow_cellpose4` | Cellpose4 (CPSAM) segmentation |
| `brieflow_stardist` | StarDist segmentation |
| `brieflow_napari` | Visualization (optional) |
| `mozzarellm` | MozzareLLM cluster analysis (`run_mozzarellm_funk.py`) |

Results saved to `results/{benchmark}/`.
