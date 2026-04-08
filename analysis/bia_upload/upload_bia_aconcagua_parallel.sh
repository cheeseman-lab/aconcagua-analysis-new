#!/bin/bash
#SBATCH --job-name=mdiberna_bia_aconcagua
#SBATCH --partition=20
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=21-00:00:00
#SBATCH --array=1-8
#SBATCH --output=/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis/analysis/bia_upload/logs/mdiberna_bia_aconcagua_plate%a_%j.out
#SBATCH --error=/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis/analysis/bia_upload/logs/mdiberna_bia_aconcagua_plate%a_%j.err

# =============================================================================
# Parallel upload of aconcagua to BioStudies — one job per plate
# =============================================================================

set -euo pipefail

PLATE="plate_${SLURM_ARRAY_TASK_ID}"
BIA_USER="bs-upload"
BIA_SERVER="fasp.ebi.ac.uk"
BIA_PORT=33001
BIA_BASE="/d9/610818-a855-441f-9bd4-5e5304a7ed1b-a30517/aconcagua"
export ASPERA_SCP_PASS="vsr5nW7Y"

ARCHIVE_BASE="/archive/cheeseman/ops_data/aconcagua"
LOG_DIR="/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis/analysis/bia_upload/logs"
RATE_LIMIT="2500M"

eval "$(conda shell.bash hook)" && conda activate aspera

echo "===================================================================="
echo "BioStudies Upload: aconcagua ${PLATE}"
echo "Started: $(date)"
echo "Rate:    ${RATE_LIMIT}"
echo "===================================================================="

which ascp
ascp --version 2>&1 | head -1

run_ascp() {
    local src="$1"
    local dest="$2"
    echo "Uploading: ${src} -> ${dest}"
    ascp \
        -k 1 \
        -Q \
        -T \
        -l "${RATE_LIMIT}" \
        -P "${BIA_PORT}" \
        -r \
        --file-manifest=text \
        --file-manifest-path="${LOG_DIR}" \
        "${src}" \
        "${BIA_USER}@${BIA_SERVER}:${dest}" 2>&1
}

# Upload input_ph for this plate
if [[ -d "${ARCHIVE_BASE}/input_ph/${PLATE}" ]]; then
    run_ascp "${ARCHIVE_BASE}/input_ph/${PLATE}" "${BIA_BASE}/input_ph/"
else
    echo "WARNING: ${ARCHIVE_BASE}/input_ph/${PLATE} not found, skipping"
fi

# Upload input_sbs for this plate
if [[ -d "${ARCHIVE_BASE}/input_sbs/${PLATE}" ]]; then
    run_ascp "${ARCHIVE_BASE}/input_sbs/${PLATE}" "${BIA_BASE}/input_sbs/"
else
    echo "WARNING: ${ARCHIVE_BASE}/input_sbs/${PLATE} not found, skipping"
fi

echo ""
echo "===================================================================="
echo "Plate ${PLATE} upload complete"
echo "Ended: $(date)"
echo "===================================================================="
