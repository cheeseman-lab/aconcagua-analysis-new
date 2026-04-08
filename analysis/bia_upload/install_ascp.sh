#!/bin/bash
# Install IBM Aspera CLI (ascp) via conda in user space
# Run this on the head node before submitting the upload job

set -euo pipefail

ENV_NAME="aspera"

echo "Creating conda environment: ${ENV_NAME}"
conda create -n "${ENV_NAME}" -c hcc aspera-cli -y

# Verify installation
eval "$(conda shell.bash hook)" && conda activate "${ENV_NAME}"
echo "ascp version: $(ascp --version 2>&1 | head -1)"
echo "ascp path: $(which ascp)"
echo ""
echo "Installation complete. To activate: eval \"\$(conda shell.bash hook)\" && conda activate ${ENV_NAME}"
