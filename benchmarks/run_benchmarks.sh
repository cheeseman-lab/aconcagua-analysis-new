#!/bin/bash
#
# Brieflow Benchmark Runner
#
# This script runs all benchmarks sequentially:
# - Segmentation: Cellpose, Cellpose4, StarDist
# - Spot Calling: Standard, Spotiflow
# - Feature Extraction: CP Measure vs CP Multichannel comparison
# - Feature Timing: Individual cp_measure function timing
# - Merge: Fast mode vs Stitch mode cell matching
# - Clustering: Brieflow vs Funk et al. 2022 comparison
#
# Usage: bash run_benchmarks.sh [OPTIONS]
#
# Options:
#   --segmentation        Only run segmentation benchmarks
#   --spot-calling        Only run spot calling benchmarks
#   --feature-extraction  Only run feature extraction comparison
#   --feature-timing      Only run feature function timing
#   --merge               Only run merge mode comparison (fast vs stitch)
#   --cluster             Only run clustering comparison (Brieflow vs Funk et al.)
#   --help                Show this help message

set -e  # Exit on error

# ============================================================================
# CONFIGURATION
# ============================================================================

MAIN_ENV="brieflow_aconcagua"
CELLPOSE_ENV="brieflow_cellpose4"
STARDIST_ENV="brieflow_stardist"

# ============================================================================
# SCRIPT SETUP
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
LOG_DIR="$SCRIPT_DIR/results/logs"
mkdir -p "$LOG_DIR"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Parse command line arguments
RUN_SEGMENTATION=true
RUN_SPOT_CALLING=true
RUN_FEATURE_EXTRACTION=true
RUN_FEATURE_TIMING=true
RUN_MERGE=true
RUN_CLUSTER=true

while [[ $# -gt 0 ]]; do
    case $1 in
        --segmentation)
            RUN_SPOT_CALLING=false
            RUN_FEATURE_EXTRACTION=false
            RUN_FEATURE_TIMING=false
            RUN_MERGE=false
            RUN_CLUSTER=false
            shift
            ;;
        --spot-calling)
            RUN_SEGMENTATION=false
            RUN_FEATURE_EXTRACTION=false
            RUN_FEATURE_TIMING=false
            RUN_MERGE=false
            RUN_CLUSTER=false
            shift
            ;;
        --feature-extraction)
            RUN_SEGMENTATION=false
            RUN_SPOT_CALLING=false
            RUN_FEATURE_TIMING=false
            RUN_MERGE=false
            RUN_CLUSTER=false
            shift
            ;;
        --feature-timing)
            RUN_SEGMENTATION=false
            RUN_SPOT_CALLING=false
            RUN_FEATURE_EXTRACTION=false
            RUN_MERGE=false
            RUN_CLUSTER=false
            shift
            ;;
        --merge)
            RUN_SEGMENTATION=false
            RUN_SPOT_CALLING=false
            RUN_FEATURE_EXTRACTION=false
            RUN_FEATURE_TIMING=false
            RUN_CLUSTER=false
            shift
            ;;
        --cluster)
            RUN_SEGMENTATION=false
            RUN_SPOT_CALLING=false
            RUN_FEATURE_EXTRACTION=false
            RUN_FEATURE_TIMING=false
            RUN_MERGE=false
            shift
            ;;
        --help)
            echo "Brieflow Benchmark Runner"
            echo ""
            echo "Usage: bash run_benchmarks.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --segmentation        Only run segmentation benchmarks"
            echo "  --spot-calling        Only run spot calling benchmarks"
            echo "  --feature-extraction  Only run feature extraction comparison"
            echo "  --feature-timing      Only run feature function timing"
            echo "  --merge               Only run merge mode comparison (fast vs stitch)"
            echo "  --cluster             Only run clustering comparison (Brieflow vs Funk et al.)"
            echo "  --help                Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
}

log_success() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] SUCCESS: $1"
}

check_conda_env() {
    local env_name=$1
    if conda env list | grep -q "^${env_name} "; then
        return 0
    else
        return 1
    fi
}

run_in_env() {
    local env_name=$1
    local script_name=$2
    local log_suffix=$3
    local extra_args=${4:-""}
    local log_file="$LOG_DIR/${script_name%.py}_${log_suffix}_${TIMESTAMP}.log"

    log_info "Running $script_name in environment: $env_name"
    if [ -n "$extra_args" ]; then
        log_info "Arguments: $extra_args"
    fi
    log_info "Log file: $log_file"

    (
        source ~/.bashrc
        conda activate "$env_name"
        cd "$SCRIPT_DIR"
        python "$script_name" $extra_args 2>&1 | tee "$log_file"
    )

    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        log_success "$script_name completed successfully in $env_name"
        return 0
    else
        log_error "$script_name failed in $env_name (see log for details)"
        return 1
    fi
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

echo "============================================================================"
echo "BRIEFLOW BENCHMARK RUNNER"
echo "============================================================================"
echo ""
log_info "Project root: $PROJECT_ROOT"
log_info "Benchmark directory: $SCRIPT_DIR"
log_info "Log directory: $LOG_DIR"
echo ""

# Verify main environment exists
if ! check_conda_env "$MAIN_ENV"; then
    log_error "Main environment '$MAIN_ENV' not found!"
    exit 1
fi

log_info "Main environment: $MAIN_ENV"
echo ""

OVERALL_SUCCESS=true

# ============================================================================
# SEGMENTATION BENCHMARKS
# ============================================================================

if [ "$RUN_SEGMENTATION" = true ]; then
    echo "============================================================================"
    echo "SEGMENTATION BENCHMARKS"
    echo "============================================================================"

    log_info "Running Cellpose 3 segmentation benchmark..."
    if ! run_in_env "$MAIN_ENV" "segmentation.py" "cellpose3" "--method cellpose"; then
        OVERALL_SUCCESS=false
    fi

    if check_conda_env "$STARDIST_ENV"; then
        log_info "Running StarDist segmentation benchmark..."
        if ! run_in_env "$STARDIST_ENV" "segmentation.py" "stardist" "--method stardist"; then
            OVERALL_SUCCESS=false
        fi
    else
        log_error "StarDist environment '$STARDIST_ENV' not found - skipping"
    fi

    if check_conda_env "$CELLPOSE_ENV"; then
        log_info "Running Cellpose 4 (CPSAM) segmentation benchmark..."
        if ! run_in_env "$CELLPOSE_ENV" "segmentation.py" "cellpose4" "--method cellpose4"; then
            OVERALL_SUCCESS=false
        fi
    else
        log_error "Cellpose 4 environment '$CELLPOSE_ENV' not found - skipping"
    fi

    log_success "Segmentation benchmarks completed"
fi

# ============================================================================
# SPOT CALLING BENCHMARKS
# ============================================================================

if [ "$RUN_SPOT_CALLING" = true ]; then
    echo ""
    echo "============================================================================"
    echo "SPOT CALLING BENCHMARKS"
    echo "============================================================================"

    log_info "Running spot calling benchmark..."
    if ! run_in_env "$MAIN_ENV" "spot_calling.py" "main"; then
        OVERALL_SUCCESS=false
    fi

    log_success "Spot calling benchmarks completed"
fi

# ============================================================================
# FEATURE EXTRACTION BENCHMARKS (COMPARISON)
# ============================================================================

if [ "$RUN_FEATURE_EXTRACTION" = true ]; then
    echo ""
    echo "============================================================================"
    echo "FEATURE EXTRACTION BENCHMARKS (COMPARISON)"
    echo "============================================================================"

    log_info "Running feature extraction comparison benchmark..."
    if ! run_in_env "$MAIN_ENV" "feature_extraction.py" "comparison" "--mode comparison"; then
        OVERALL_SUCCESS=false
    fi

    log_success "Feature extraction comparison completed"
fi

# ============================================================================
# FEATURE EXTRACTION BENCHMARKS (TIMING)
# ============================================================================

if [ "$RUN_FEATURE_TIMING" = true ]; then
    echo ""
    echo "============================================================================"
    echo "FEATURE EXTRACTION BENCHMARKS (TIMING)"
    echo "============================================================================"

    log_info "Running feature function timing benchmark..."
    if ! run_in_env "$MAIN_ENV" "feature_extraction.py" "timing" "--mode timing"; then
        OVERALL_SUCCESS=false
    fi

    log_success "Feature function timing completed"
fi

# ============================================================================
# MERGE BENCHMARKS (FAST vs STITCH)
# ============================================================================

if [ "$RUN_MERGE" = true ]; then
    echo ""
    echo "============================================================================"
    echo "MERGE BENCHMARKS (FAST vs STITCH)"
    echo "============================================================================"

    log_info "Running merge mode comparison benchmark..."
    if ! run_in_env "$MAIN_ENV" "merge.py" "main"; then
        OVERALL_SUCCESS=false
    fi

    log_success "Merge benchmarks completed"
fi

# ============================================================================
# CLUSTERING BENCHMARKS (BRIEFLOW vs FUNK ET AL. 2022)
# ============================================================================

if [ "$RUN_CLUSTER" = true ]; then
    echo ""
    echo "============================================================================"
    echo "CLUSTERING BENCHMARKS (BRIEFLOW vs FUNK ET AL. 2022)"
    echo "============================================================================"

    log_info "Running clustering comparison benchmark..."
    if ! run_in_env "$MAIN_ENV" "cluster.py" "main"; then
        OVERALL_SUCCESS=false
    fi

    log_success "Clustering benchmarks completed"
fi

# ============================================================================
# SUMMARY
# ============================================================================

echo ""
echo "============================================================================"
echo "BENCHMARK SUMMARY"
echo "============================================================================"
echo ""
log_info "Results directory: $SCRIPT_DIR/results/"
log_info "Log files: $LOG_DIR/"
echo ""

if [ "$OVERALL_SUCCESS" = true ]; then
    log_success "All benchmarks completed successfully!"
    exit 0
else
    log_error "Some benchmarks failed - check logs for details"
    exit 1
fi
