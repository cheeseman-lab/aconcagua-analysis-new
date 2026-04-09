# Aconcagua OPS Screen Analysis (v2.0.0)

Reanalysis of the optical pooled CRISPR screen from [Funk et al. 2022 (Cell)](https://doi.org/10.1016/j.cell.2022.12.009) using [brieflow v1.4.10](https://github.com/cheeseman-lab/brieflow).

## Screen Overview

| Parameter | Value |
|-----------|-------|
| Cell line | HeLa-TetR-Cas9 (A7) |
| Genes targeted | 5,299 |
| Total guides | 20,445 (4 per gene) |
| Plates | 8 |
| Phenotype markers | DAPI, Tubulin (alpha-tubulin-FITC), γ-H2AX (AF594), Phalloidin (AF750) |

## Analysis

Raw data was processed end-to-end with [brieflow](https://github.com/cheeseman-lab/brieflow), covering preprocessing, SBS decoding, phenotype feature extraction, merge, aggregate, and clustering steps. Configuration notebooks and Slurm scripts are in `analysis/`.

## Data

Raw imaging data (ND2 format) and brieflow outputs are stored at `/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis/`. Source data used for figure generation is in `analysis/source_data/`.

## Citation

> *Citation will be added upon publication.*
