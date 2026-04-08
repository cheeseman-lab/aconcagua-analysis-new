#!/bin/bash
#SBATCH --job-name=mdiberna_replot_eval
#SBATCH --partition=20
#SBATCH --account=wibrusers
#SBATCH --mem=700G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=%x_%j.out

set -euo pipefail
export PYTHONUNBUFFERED=1

cd /lab/ops_analysis_ssd/cheeseman/aconcagua-analysis/benchmarks

eval "$(conda shell.bash hook)" && conda activate brieflow_aconcagua

echo "Starting replot_eval_figures.py at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Memory: 700G"
echo ""

python replot_eval_figures.py

echo ""
echo "Finished at $(date)"
