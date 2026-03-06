# Benchmarks

Compares methods at four pipeline stages.

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

Analyzes merge results from Snakemake pipeline runs. Extracts timing from slurm logs and quality metrics from parquet outputs.

**Output directories:**
- `merge/` - Fast approach (default, `approach: fast` in config)
- `merge_stitch/` - Stitch approach (requires `approach: stitch` in config)

To generate stitch outputs for comparison, set `merge.approach: stitch` in `config/config.yml` and run the merge step.

**Conclusion:** Fast approach is recommended. Stitch requires ~2TB RAM for triangle hash distance matrix (beyond cluster capacity), falls back to identity transform with degraded alignment.

## Clustering (`cluster.py`)

**Comparison:** Brieflow vs Funk et al. 2022 original analysis
**Metrics:**
- STRING F1 score (protein-protein interaction enrichment)
- CORUM/KEGG enrichment (% clusters with significant pathway enrichment)
- Cluster size distributions
- Precision-Recall curves for pathway detection

Addresses reviewer feedback requesting rigorous quantitative benchmarking and head-to-head comparison with PR curves. Demonstrates improved clustering quality through:
- Enhanced biological coherence (CORUM/KEGG enrichment)
- Better protein interaction network capture (STRING F1)
- Systematic evaluation across both Interphase and Mitotic configurations

**Baseline data:** External comparison uses Funk et al. 2022 supplementary data (Table S2) from the original OPS manuscript.

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
- **segmentation**: Side-by-side segmentation method comparison
- **spots**: Spot detection overlay comparison
- **overlay**: Single image with segmentation overlay
- **panel**: Multi-panel grid of samples

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
python cluster.py              # Compare clustering with Funk et al. 2022
python cluster.py --plots-only # Regenerate plots only
python classifier.py           # Train and evaluate cell stage classifier

# Generate visualizations (requires Napari, run locally)
python visualization.py --type all
python visualization.py --type segmentation --sample P-1_W-A1_T-1
python visualization.py --type panel --rows 2 --cols 4
```

## Environments

| Environment | Used For |
|-------------|----------|
| `brieflow_aconcagua` | Main: Cellpose3, spot calling, feature extraction, merge, clustering |
| `brieflow_cellpose4` | Cellpose4 (CPSAM) segmentation |
| `brieflow_stardist` | StarDist segmentation |
| `brieflow_napari` | Visualization (optional) |

### Setting up brieflow_napari (for visualization)

```bash
# Create environment with all visualization dependencies (locally)
conda create -n brieflow_napari python=3.11 -y
conda activate brieflow_napari
pip install napari[all] tifffile numpy pandas matplotlib scikit-image seaborn

# Run visualizations (offscreen rendering enabled by default)
python visualization.py --type all
```

Results saved to `results/{benchmark}/`.
