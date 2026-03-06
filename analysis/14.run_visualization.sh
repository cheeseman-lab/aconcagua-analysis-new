#!/bin/bash

export BRIEFLOW_OUTPUT_PATH="brieflow_output/"
export CONFIG_PATH="config/config.yml"
export SCREEN_PATH="screen.yaml"
export PERTURBATION_LIBRARY_PATH="config/pool10_design.tsv"
export FEATURE_DOC_PATH="../brieflow/workflow/lib/external/CP_EMULATOR_FEATURES.md"

# Start Streamlit server, force bind to 0.0.0.0
exec streamlit run ../brieflow/visualization/Cluster_Analysis.py --server.address=0.0.0.0 "$@"