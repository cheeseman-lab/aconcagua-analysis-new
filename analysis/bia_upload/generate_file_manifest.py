"""
Generate BioImage Archive file list manifests from brieflow sample TSVs.
Converts local archive paths to BIA-relative paths.

Local:  /archive/cheeseman/ops_data/aconcagua/input_sbs/...
BIA:    aconcagua/input_sbs/...
"""

import pandas as pd
from pathlib import Path

CONFIG_DIR = Path("/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis/analysis/config")
OUT_DIR = Path("/lab/ops_analysis_ssd/cheeseman/aconcagua-analysis/analysis/bia_upload")

LOCAL_PREFIXES = [
    "/archive/cheeseman/ops_data/aconcagua/",
    "/lab/ops_data/aconcagua/",
]
BIA_PREFIX = "aconcagua/"

def convert_path(fp: str) -> str:
    for prefix in LOCAL_PREFIXES:
        if fp.startswith(prefix):
            return BIA_PREFIX + fp[len(prefix):]
    return fp

# --- SBS samples ---
sbs = pd.read_csv(CONFIG_DIR / "sbs_samples.tsv", sep="\t")
sbs["sample_fp"] = sbs["sample_fp"].apply(convert_path)
sbs["data_type"] = "sbs"

# --- Phenotype samples ---
ph = pd.read_csv(CONFIG_DIR / "phenotype_samples.tsv", sep="\t")
ph["sample_fp"] = ph["sample_fp"].apply(convert_path)
ph["data_type"] = "phenotype"

# --- Combine and save ---
combined = pd.concat([sbs, ph], ignore_index=True)
combined = combined.rename(columns={"sample_fp": "Files"})
combined.to_csv(OUT_DIR / "file_manifest_bia.tsv", sep="\t", index=False)
print(f"Combined manifest: {len(combined)} rows → {OUT_DIR / 'file_manifest_bia.tsv'}")
print(f"  SBS: {len(sbs)} rows")
print(f"  Phenotype: {len(ph)} rows")
print("\nSample rows:")
print(combined[["Files", "data_type"]].head(3).to_string())
print(combined[combined["data_type"] == "phenotype"][["Files", "data_type"]].head(3).to_string())
