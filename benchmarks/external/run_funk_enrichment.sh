#!/bin/bash
# Run Funk enrichment analysis in parallel

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate brieflow_aconcagua

# Create log directory
mkdir -p ./logs

echo "================================================================================"
echo "Running Funk Enrichment Analysis in Parallel"
echo "================================================================================"
echo "Using conda environment: brieflow_aconcagua"
echo "Number of parallel jobs: 22 (all clustering results at once)"
echo ""

# Check if GNU parallel is available
if command -v parallel &> /dev/null; then
    echo "Using GNU parallel for parallelization"

    # Generate all combinations and run with parallel
    parallel -j 22 --progress \
        "python run_funk_enrichment.py {1} {2} > logs/funk_{1}_k{2}.log 2>&1" \
        ::: Interphase Mitotic \
        ::: {5..15}

else
    echo "GNU parallel not found, using background jobs (max 22 concurrent)"

    # Counter for background jobs
    job_count=0
    max_jobs=22

    # Run all combinations
    for cell_class in Interphase Mitotic; do
        for k in {5..15}; do
            # Wait if we've reached max concurrent jobs
            while [ $(jobs -r | wc -l) -ge $max_jobs ]; do
                sleep 1
            done

            # Start job in background
            echo "Starting: Funk_${cell_class}_k${k}"
            python run_funk_enrichment.py "$cell_class" "$k" \
                > "logs/funk_${cell_class}_k${k}.log" 2>&1 &

            job_count=$((job_count + 1))
        done
    done

    # Wait for all background jobs to complete
    echo ""
    echo "Waiting for all jobs to complete..."
    wait
fi

echo ""
echo "================================================================================"
echo "Collecting Results"
echo "================================================================================"

# Create summary table
echo ""
printf "%-30s %12s %10s %10s\n" "Run" "STRING F1" "CORUM %" "KEGG %"
echo "--------------------------------------------------------------------------------"

for cell_class in Interphase Mitotic; do
    for k in {5..15}; do
        run_name="Funk_${cell_class}_k${k}"
        log_file="logs/funk_${cell_class}_k${k}.log"

        if [ -f "$log_file" ]; then
            # Extract metrics from log file
            string_f1=$(grep "STRING F1:" "$log_file" | awk '{print $3}' || echo "N/A")
            corum_pct=$(grep "CORUM enriched:" "$log_file" | awk '{print $3}' || echo "N/A")
            kegg_pct=$(grep "KEGG enriched:" "$log_file" | awk '{print $3}' || echo "N/A")

            printf "%-30s %12s %10s %10s\n" "$run_name" "$string_f1" "$corum_pct" "$kegg_pct"
        else
            printf "%-30s %12s %10s %10s\n" "$run_name" "FAILED" "-" "-"
        fi
    done
done

echo "================================================================================"
echo ""

# Count successful runs
success_count=$(find results/cluster/Funk_*/CB-Real__global_metrics.json 2>/dev/null | wc -l)
echo "✓ Successfully processed $success_count/22 clustering results"
echo ""
echo "Logs saved to: $SCRIPT_DIR/logs/"
echo "================================================================================"
