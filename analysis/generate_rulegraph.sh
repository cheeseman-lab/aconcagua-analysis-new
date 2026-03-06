#!/bin/bash

# Generate a rulegraph of the Snakefile
snakemake \
    --snakefile "../brieflow/workflow/Snakefile" \
    --configfile "config/config.yml" \
    --until all_cluster \
    --rulegraph | dot -Tsvg -o "../images/brieflow_rulegraph.svg"
