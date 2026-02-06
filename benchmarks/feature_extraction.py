"""
Feature Extraction Benchmark

Two benchmark modes:
- comparison: Compare CP Measure vs CP Multichannel pipelines
- timing: Time individual cp_measure functions to identify bottlenecks

Usage:
    python feature_extraction.py                    # Run comparison (default)
    python feature_extraction.py --mode comparison  # Compare pipelines
    python feature_extraction.py --mode timing      # Time individual functions
    python feature_extraction.py --plots-only       # Regenerate plots including correlation analysis
"""

import os
import sys
import time
import random
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tifffile import imread
from pathlib import Path
from datetime import datetime
from scipy.stats import pearsonr
import psutil
import gc
from tqdm import tqdm

# Add brieflow to path (relative to benchmarks directory)
BENCHMARKS_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = BENCHMARKS_DIR.parent
sys.path.insert(0, str(PROJECT_ROOT / "brieflow" / "workflow"))

random.seed(42)

# Configure paths (relative to project root)
SOURCE_DIR = PROJECT_ROOT / "analysis" / "brieflow_output" / "phenotype" / "images"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "feature_extraction"

# Default channel configuration
CHANNEL_NAMES = ["dapi", "tubulin", "gh2ax", "phalloidin"]
MASKS = ["nucleus", "cell", "cytoplasm"]

# =============================================================================
# FEATURE MAPPING AND CORRELATION ANALYSIS
# =============================================================================


def get_feature_mappings():
    """
    Define explicit feature mappings between cp_emulator and cp_measure.
    Returns dict of {cluster_name: [(emulator_col, measure_col, expected_corr), ...]}
    """
    mappings = {}

    # For each mask-channel combination
    for mask in MASKS:
        for channel in CHANNEL_NAMES:
            prefix_emu = f"{mask}_{channel}"
            prefix_mea = f"{mask}_{channel}__"

            # INTENSITY FEATURES
            intensity_maps = [
                (
                    f"{prefix_emu}_int",
                    f"{prefix_mea}_Intensity_IntegratedIntensity",
                    "high",
                ),
                (f"{prefix_emu}_mean", f"{prefix_mea}_Intensity_MeanIntensity", "high"),
                (f"{prefix_emu}_std", f"{prefix_mea}_Intensity_StdIntensity", "high"),
                (f"{prefix_emu}_max", f"{prefix_mea}_Intensity_MaxIntensity", "high"),
                (f"{prefix_emu}_min", f"{prefix_mea}_Intensity_MinIntensity", "high"),
                (
                    f"{prefix_emu}_int_edge",
                    f"{prefix_mea}_Intensity_IntegratedIntensityEdge",
                    "low",
                ),
                (
                    f"{prefix_emu}_mean_edge",
                    f"{prefix_mea}_Intensity_MeanIntensityEdge",
                    "low",
                ),
                (
                    f"{prefix_emu}_std_edge",
                    f"{prefix_mea}_Intensity_StdIntensityEdge",
                    "low",
                ),
                (
                    f"{prefix_emu}_max_edge",
                    f"{prefix_mea}_Intensity_MaxIntensityEdge",
                    "low",
                ),
                (
                    f"{prefix_emu}_min_edge",
                    f"{prefix_mea}_Intensity_MinIntensityEdge",
                    "low",
                ),
                (
                    f"{prefix_emu}_mass_displacement",
                    f"{prefix_mea}_Intensity_MassDisplacement",
                    "medium",
                ),
                (
                    f"{prefix_emu}_lower_quartile",
                    f"{prefix_mea}_Intensity_LowerQuartileIntensity",
                    "high",
                ),
                (
                    f"{prefix_emu}_median",
                    f"{prefix_mea}_Intensity_MedianIntensity",
                    "high",
                ),
                (f"{prefix_emu}_mad", f"{prefix_mea}_Intensity_MADIntensity", "medium"),
                (
                    f"{prefix_emu}_upper_quartile",
                    f"{prefix_mea}_Intensity_UpperQuartileIntensity",
                    "high",
                ),
            ]
            mappings.setdefault("intensity", []).extend(intensity_maps)

            # LOCATION FEATURES
            location_maps = [
                (
                    f"{prefix_emu}_center_mass_r",
                    f"{prefix_mea}_Location_CenterMassIntensity_Y",
                    "low",
                ),
                (
                    f"{prefix_emu}_center_mass_c",
                    f"{prefix_mea}_Location_CenterMassIntensity_X",
                    "low",
                ),
                (
                    f"{prefix_emu}_max_location_r",
                    f"{prefix_mea}_Location_MaxIntensity_Y",
                    "low",
                ),
                (
                    f"{prefix_emu}_max_location_c",
                    f"{prefix_mea}_Location_MaxIntensity_X",
                    "low",
                ),
            ]
            mappings.setdefault("location", []).extend(location_maps)

            # RADIAL DISTRIBUTION FEATURES
            radial_maps = []
            for i in range(4):
                radial_maps.extend(
                    [
                        (
                            f"{prefix_emu}_frac_at_d_{i}",
                            f"{prefix_mea}_RadialDistribution_FracAtD_{i + 1}of4",
                            "low",
                        ),
                        (
                            f"{prefix_emu}_mean_frac_{i}",
                            f"{prefix_mea}_RadialDistribution_MeanFrac_{i + 1}of4",
                            "low",
                        ),
                        (
                            f"{prefix_emu}_radial_cv_{i}",
                            f"{prefix_mea}_RadialDistribution_RadialCV_{i + 1}of4",
                            "low",
                        ),
                    ]
                )
            mappings.setdefault("radial_distribution", []).extend(radial_maps)

            # HARALICK TEXTURE FEATURES (emulator averages directions, measure has per-direction)
            haralick_names = [
                "AngularSecondMoment",
                "Contrast",
                "Correlation",
                "Variance",
                "InverseDifferenceMoment",
                "SumAverage",
                "SumVariance",
                "SumEntropy",
                "Entropy",
                "DifferenceVariance",
                "DifferenceEntropy",
                "InfoMeas1",
                "InfoMeas2",
            ]
            haralick_maps = []
            for idx, name in enumerate(haralick_names):
                # Emulator has haralick_5_0..12, measure has per-direction (we use 00)
                haralick_maps.append(
                    (
                        f"{prefix_emu}_haralick_5_{idx}",
                        f"{prefix_mea}_{name}_3_00_256",
                        "medium",
                    )
                )
            mappings.setdefault("haralick", []).extend(haralick_maps)

            # WEIGHTED HU MOMENTS
            whu_maps = []
            for i in range(7):
                whu_maps.append(
                    (
                        f"{prefix_emu}_weighted_hu_moments_{i}",
                        f"{prefix_mea}_HuMoment_{i}",
                        "medium",
                    )
                )
            mappings.setdefault("weighted_hu_moments", []).extend(whu_maps)

    # SHAPE FEATURES (mask-only, not channel-dependent in emulator)
    for mask in MASKS:
        # Shape features use first channel in cp_measure naming
        channel = CHANNEL_NAMES[0]  # dapi
        prefix_mea = f"{mask}_{channel}__"

        shape_maps = [
            (f"{mask}_area", f"{prefix_mea}_Area", "high"),
            (f"{mask}_perimeter", f"{prefix_mea}_Perimeter", "high"),
            (f"{mask}_convex_area", f"{prefix_mea}_ConvexArea", "high"),
            (f"{mask}_form_factor", f"{prefix_mea}_FormFactor", "high"),
            (f"{mask}_solidity", f"{prefix_mea}_Solidity", "high"),
            (f"{mask}_extent", f"{prefix_mea}_Extent", "high"),
            (f"{mask}_euler_number", f"{prefix_mea}_EulerNumber", "high"),
            (f"{mask}_eccentricity", f"{prefix_mea}_Eccentricity", "high"),
            (f"{mask}_major_axis", f"{prefix_mea}_MajorAxisLength", "high"),
            (f"{mask}_minor_axis", f"{prefix_mea}_MinorAxisLength", "high"),
            (f"{mask}_orientation", f"{prefix_mea}_Orientation", "high"),
            (f"{mask}_compactness", f"{prefix_mea}_Compactness", "medium"),
            (f"{mask}_max_radius", f"{prefix_mea}_MaximumRadius", "high"),
            (f"{mask}_median_radius", f"{prefix_mea}_MedianRadius", "high"),
            (f"{mask}_mean_radius", f"{prefix_mea}_MeanRadius", "high"),
        ]
        mappings.setdefault("shape", []).extend(shape_maps)

        # FERET FEATURES
        feret_maps = [
            (f"{mask}_min_feret_diameter", f"{prefix_mea}_MinFeretDiameter", "high"),
            (f"{mask}_max_feret_diameter", f"{prefix_mea}_MaxFeretDiameter", "high"),
        ]
        mappings.setdefault("feret", []).extend(feret_maps)

        # BOUNDING BOX FEATURES
        bounding_maps = [
            (f"{mask}_bounds_0", f"{prefix_mea}_BoundingBoxMinimum_X", "high"),
            (f"{mask}_bounds_1", f"{prefix_mea}_BoundingBoxMaximum_X", "high"),
            (f"{mask}_bounds_2", f"{prefix_mea}_BoundingBoxMinimum_Y", "high"),
            (f"{mask}_bounds_3", f"{prefix_mea}_BoundingBoxMaximum_Y", "high"),
        ]
        mappings.setdefault("bounding", []).extend(bounding_maps)

        # HU MOMENTS
        hu_maps = []
        for i in range(7):
            hu_maps.append(
                (f"{mask}_hu_moments_{i}", f"{prefix_mea}_HuMoment_{i}", "high")
            )
        mappings.setdefault("hu_moments", []).extend(hu_maps)

        # ZERNIKE FEATURES
        zernike_pairs = [
            (0, 0),
            (1, 1),
            (2, 0),
            (2, 2),
            (3, 1),
            (3, 3),
            (4, 0),
            (4, 2),
            (4, 4),
            (5, 1),
            (5, 3),
            (5, 5),
            (6, 0),
            (6, 2),
            (6, 4),
            (6, 6),
            (7, 1),
            (7, 3),
            (7, 5),
            (7, 7),
            (8, 0),
            (8, 2),
            (8, 4),
            (8, 6),
            (8, 8),
            (9, 1),
            (9, 3),
            (9, 5),
            (9, 7),
            (9, 9),
        ]
        zernike_maps = []
        for n, m in zernike_pairs:
            zernike_maps.append(
                (f"{mask}_zernike_{n}_{m}", f"{prefix_mea}_Zernike_{n}_{m}", "high")
            )
        mappings.setdefault("zernike", []).extend(zernike_maps)

    return mappings


def load_feature_dataframes():
    """Load all feature CSVs from cp_measure and cp_multichannel directories."""
    cp_measure_dir = OUTPUT_DIR / "cp_measure"
    cp_multi_dir = OUTPUT_DIR / "cp_multichannel"

    measure_files = sorted(cp_measure_dir.glob("*__features.csv"))
    multi_files = sorted(cp_multi_dir.glob("*__features.csv"))

    if not measure_files or not multi_files:
        print("No feature files found!")
        return None, None

    # Load and concatenate all files
    measure_dfs = []
    multi_dfs = []

    for mf in measure_files:
        base = mf.name
        matching_multi = cp_multi_dir / base
        if matching_multi.exists():
            print(f"Loading {base}...")
            measure_dfs.append(pd.read_csv(mf))
            multi_dfs.append(pd.read_csv(matching_multi))

    if not measure_dfs:
        print("No matching file pairs found!")
        return None, None

    measure_df = pd.concat(measure_dfs, ignore_index=True)
    multi_df = pd.concat(multi_dfs, ignore_index=True)

    print(f"Loaded {len(measure_df)} rows from {len(measure_dfs)} files")
    print(f"CP Measure columns: {len(measure_df.columns)}")
    print(f"CP Multichannel columns: {len(multi_df.columns)}")

    return measure_df, multi_df


def compute_feature_correlations(measure_df, multi_df, mappings):
    """
    Compute Pearson correlations for each mapped feature pair.
    Returns DataFrame with correlation results.
    """
    results = []

    for cluster, feature_maps in mappings.items():
        for emu_col, mea_col, expected in feature_maps:
            # Check if columns exist
            emu_exists = emu_col in multi_df.columns
            mea_exists = mea_col in measure_df.columns

            if not emu_exists or not mea_exists:
                results.append(
                    {
                        "cluster": cluster,
                        "emulator_feature": emu_col,
                        "measure_feature": mea_col,
                        "expected_correlation": expected,
                        "pearson_r": np.nan,
                        "p_value": np.nan,
                        "n_valid": 0,
                        "emulator_exists": emu_exists,
                        "measure_exists": mea_exists,
                    }
                )
                continue

            # Get values
            x = multi_df[emu_col].values
            y = measure_df[mea_col].values

            # Filter valid (non-NaN, non-Inf) pairs
            valid_mask = ~(np.isnan(x) | np.isnan(y) | np.isinf(x) | np.isinf(y))
            x_valid = x[valid_mask]
            y_valid = y[valid_mask]

            if len(x_valid) < 10:
                results.append(
                    {
                        "cluster": cluster,
                        "emulator_feature": emu_col,
                        "measure_feature": mea_col,
                        "expected_correlation": expected,
                        "pearson_r": np.nan,
                        "p_value": np.nan,
                        "n_valid": len(x_valid),
                        "emulator_exists": True,
                        "measure_exists": True,
                    }
                )
                continue

            # Compute correlation
            try:
                r, p = pearsonr(x_valid, y_valid)
                results.append(
                    {
                        "cluster": cluster,
                        "emulator_feature": emu_col,
                        "measure_feature": mea_col,
                        "expected_correlation": expected,
                        "pearson_r": r,
                        "p_value": p,
                        "n_valid": len(x_valid),
                        "emulator_exists": True,
                        "measure_exists": True,
                    }
                )
            except Exception as e:
                results.append(
                    {
                        "cluster": cluster,
                        "emulator_feature": emu_col,
                        "measure_feature": mea_col,
                        "expected_correlation": expected,
                        "pearson_r": np.nan,
                        "p_value": np.nan,
                        "n_valid": len(x_valid),
                        "emulator_exists": True,
                        "measure_exists": True,
                        "error": str(e),
                    }
                )

    return pd.DataFrame(results)


def generate_correlation_plots(corr_df):
    """Generate correlation analysis plots."""
    from plot_style import setup_plot_style

    setup_plot_style()

    # Create output directory for correlation plots
    corr_output_dir = OUTPUT_DIR / "correlation_analysis"
    corr_output_dir.mkdir(exist_ok=True)

    # Filter to valid correlations
    valid_df = corr_df[corr_df["pearson_r"].notna()].copy()

    if len(valid_df) == 0:
        print("No valid correlations to plot!")
        return

    print(f"\nPlotting {len(valid_df)} valid feature correlations...")

    # 1. CLUSTER SUMMARY HEATMAP
    cluster_summary = (
        valid_df.groupby("cluster")
        .agg(
            {
                "pearson_r": ["mean", "std", "count"],
            }
        )
        .round(3)
    )
    cluster_summary.columns = ["mean_r", "std_r", "n_features"]
    cluster_summary = cluster_summary.sort_values("mean_r", ascending=False)

    fig, ax = plt.subplots(figsize=(10, 8))
    clusters = cluster_summary.index.tolist()
    means = cluster_summary["mean_r"].values
    stds = cluster_summary["std_r"].values
    counts = cluster_summary["n_features"].values

    colors = plt.cm.RdYlGn((means + 1) / 2)  # Map [-1, 1] to colormap
    bars = ax.barh(
        range(len(clusters)),
        means,
        xerr=stds,
        color=colors,
        edgecolor="black",
        linewidth=0.5,
        capsize=3,
    )

    ax.set_yticks(range(len(clusters)))
    ax.set_yticklabels([f"{c} (n={n})" for c, n in zip(clusters, counts)])
    ax.set_xlabel("Mean Pearson Correlation (r)")
    ax.set_title(
        "Feature Cluster Correlation: CP Emulator vs CP Measure", fontweight="bold"
    )
    ax.axvline(x=0, color="black", linestyle="-", linewidth=0.5)
    ax.axvline(
        x=0.7,
        color="green",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
        label="High (r=0.7)",
    )
    ax.axvline(
        x=0.4,
        color="orange",
        linestyle="--",
        linewidth=1,
        alpha=0.7,
        label="Medium (r=0.4)",
    )
    ax.set_xlim(-0.2, 1.1)
    ax.legend(loc="lower right")

    # Add value labels
    for i, (mean, std) in enumerate(zip(means, stds)):
        ax.text(mean + std + 0.02, i, f"{mean:.2f}", va="center", fontsize=9)

    plt.tight_layout()
    plt.savefig(
        corr_output_dir / "cluster_correlation_summary.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()
    print("  Saved cluster_correlation_summary.png")

    # 2. EXPECTED VS ACTUAL CORRELATION COMPARISON
    fig, ax = plt.subplots(figsize=(10, 6))

    expected_map = {"high": 0.8, "medium": 0.5, "low": 0.2}
    valid_df["expected_numeric"] = valid_df["expected_correlation"].map(expected_map)

    for expected, color in [("high", "green"), ("medium", "orange"), ("low", "red")]:
        subset = valid_df[valid_df["expected_correlation"] == expected]
        if len(subset) > 0:
            ax.scatter(
                subset["expected_numeric"] + np.random.normal(0, 0.02, len(subset)),
                subset["pearson_r"],
                alpha=0.5,
                s=30,
                c=color,
                label=f"Expected {expected}",
            )

    ax.plot([0, 1], [0, 1], "k--", alpha=0.5, label="Perfect match")
    ax.set_xlabel("Expected Correlation Level")
    ax.set_ylabel("Actual Pearson r")
    ax.set_title("Expected vs Actual Correlations", fontweight="bold")
    ax.set_xticks([0.2, 0.5, 0.8])
    ax.set_xticklabels(["Low", "Medium", "High"])
    ax.legend()
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.2, 1.1)

    plt.tight_layout()
    plt.savefig(
        corr_output_dir / "expected_vs_actual_correlation.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()
    print("  Saved expected_vs_actual_correlation.png")

    # 3. PER-CLUSTER DETAILED PLOTS
    for cluster in valid_df["cluster"].unique():
        cluster_data = valid_df[valid_df["cluster"] == cluster].copy()
        cluster_data = cluster_data.sort_values("pearson_r", ascending=True)

        if len(cluster_data) == 0:
            continue

        # Truncate feature names for display
        cluster_data["feature_short"] = cluster_data["emulator_feature"].apply(
            lambda x: x.split("_")[-1] if len(x) > 30 else x
        )

        fig_height = max(6, len(cluster_data) * 0.3)
        fig, ax = plt.subplots(figsize=(10, fig_height))

        colors = plt.cm.RdYlGn((cluster_data["pearson_r"].values + 1) / 2)
        bars = ax.barh(
            range(len(cluster_data)),
            cluster_data["pearson_r"].values,
            color=colors,
            edgecolor="black",
            linewidth=0.3,
        )

        ax.set_yticks(range(len(cluster_data)))
        ax.set_yticklabels(cluster_data["feature_short"].values, fontsize=8)
        ax.set_xlabel("Pearson Correlation (r)")
        ax.set_title(
            f"{cluster.upper()} Features: CP Emulator vs CP Measure", fontweight="bold"
        )
        ax.axvline(x=0, color="black", linestyle="-", linewidth=0.5)
        ax.set_xlim(-0.2, 1.1)

        plt.tight_layout()
        plt.savefig(
            corr_output_dir / f"cluster_{cluster}_correlations.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()

    print("  Saved per-cluster correlation plots")

    # 4. DISTRIBUTION HISTOGRAM
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(
        valid_df["pearson_r"], bins=30, edgecolor="black", alpha=0.7, color="steelblue"
    )
    ax.axvline(
        x=valid_df["pearson_r"].mean(),
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Mean: {valid_df['pearson_r'].mean():.2f}",
    )
    ax.axvline(
        x=valid_df["pearson_r"].median(),
        color="orange",
        linestyle="--",
        linewidth=2,
        label=f"Median: {valid_df['pearson_r'].median():.2f}",
    )

    ax.set_xlabel("Pearson Correlation (r)")
    ax.set_ylabel("Number of Features")
    ax.set_title("Distribution of Feature Correlations", fontweight="bold")
    ax.legend()

    plt.tight_layout()
    plt.savefig(
        corr_output_dir / "correlation_distribution.png", dpi=300, bbox_inches="tight"
    )
    plt.close()
    print("  Saved correlation_distribution.png")

    # 5. SUMMARY STATISTICS
    print("\n" + "=" * 60)
    print("CORRELATION ANALYSIS SUMMARY")
    print("=" * 60)
    print(f"\nTotal feature pairs analyzed: {len(corr_df)}")
    print(f"Valid correlations computed: {len(valid_df)}")
    print(f"Missing features: {len(corr_df) - len(valid_df)}")

    print("\nOverall correlation statistics:")
    print(f"  Mean r: {valid_df['pearson_r'].mean():.3f}")
    print(f"  Median r: {valid_df['pearson_r'].median():.3f}")
    print(f"  Std r: {valid_df['pearson_r'].std():.3f}")

    print("\nCorrelation by expected level:")
    for expected in ["high", "medium", "low"]:
        subset = valid_df[valid_df["expected_correlation"] == expected]
        if len(subset) > 0:
            print(
                f"  {expected.upper()}: mean={subset['pearson_r'].mean():.3f}, n={len(subset)}"
            )

    print("\nCorrelation by cluster:")
    print(cluster_summary.to_string())

    return cluster_summary


def generate_scatter_plots(corr_df, measure_df, multi_df, max_features_per_cluster=20):
    """
    Generate log-scale scatter plots for each feature cluster.

    Each plot shows all feature pairs within a cluster, colored by feature type,
    allowing visual inspection of correlation quality.
    """
    from plot_style import setup_plot_style, FIGSIZE, save_figure

    setup_plot_style()

    # Create output directory
    scatter_dir = OUTPUT_DIR / "scatter_plots"
    scatter_dir.mkdir(exist_ok=True)

    print("\nGenerating scatter plots...")

    def get_base_feature(feat):
        """Extract base feature name without mask/channel prefix."""
        parts = feat.split("_")
        if parts[0] in ["nucleus", "cell", "cytoplasm"]:
            if len(parts) > 2 and parts[1] in CHANNEL_NAMES:
                return "_".join(parts[2:])
            return "_".join(parts[1:])
        return feat

    clusters = corr_df["cluster"].unique()

    for cluster_name in clusters:
        # Get valid features for this cluster
        cluster_data = corr_df[
            (corr_df["cluster"] == cluster_name)
            & (corr_df["pearson_r"].notna())
            & (corr_df["emulator_exists"] == True)
            & (corr_df["measure_exists"] == True)
        ].copy()

        if len(cluster_data) == 0:
            continue

        # Sort by correlation and get base features
        cluster_data = cluster_data.sort_values("pearson_r", ascending=False)
        cluster_data["base_feature"] = cluster_data["emulator_feature"].apply(
            get_base_feature
        )
        unique_bases = cluster_data["base_feature"].unique()

        # Sample features - up to 3 examples per base feature type
        sampled_rows = []
        for base in unique_bases[:max_features_per_cluster]:
            base_rows = cluster_data[cluster_data["base_feature"] == base]
            sampled_rows.append(base_rows.head(3))

        if not sampled_rows:
            continue

        plot_data = pd.concat(sampled_rows)

        # Create figure with proper styling
        fig, ax = plt.subplots(figsize=FIGSIZE["square"])

        # Color palette for different features
        n_colors = min(len(unique_bases), max_features_per_cluster)
        colors = plt.cm.tab20(np.linspace(0, 1, n_colors))
        color_map = {
            base: colors[i]
            for i, base in enumerate(unique_bases[:max_features_per_cluster])
        }

        # Track which features have been added to legend
        legend_added = set()

        for _, row in plot_data.iterrows():
            emu_col = row["emulator_feature"]
            mea_col = row["measure_feature"]
            base = row["base_feature"]
            r = row["pearson_r"]

            if emu_col not in multi_df.columns or mea_col not in measure_df.columns:
                continue

            x = multi_df[emu_col].values
            y = measure_df[mea_col].values

            # Filter valid values (positive for log scale)
            valid = (
                (x > 0)
                & (y > 0)
                & ~np.isnan(x)
                & ~np.isnan(y)
                & ~np.isinf(x)
                & ~np.isinf(y)
            )
            x_valid = x[valid]
            y_valid = y[valid]

            if len(x_valid) < 10:
                continue

            # Subsample for plotting if too many points
            if len(x_valid) > 2000:
                idx = np.random.choice(len(x_valid), 2000, replace=False)
                x_valid = x_valid[idx]
                y_valid = y_valid[idx]

            # Add to legend only once per base feature
            label = f"{base} (r={r:.2f})" if base not in legend_added else None
            ax.scatter(
                x_valid,
                y_valid,
                alpha=0.4,
                s=15,
                c=[color_map[base]],
                label=label,
                edgecolors="none",
            )

            if label:
                legend_added.add(base)

        # Set log scale
        ax.set_xscale("log")
        ax.set_yscale("log")

        # Add 1:1 reference line
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        lims = [max(xlim[0], ylim[0]), min(xlim[1], ylim[1])]
        ax.plot(lims, lims, "k--", alpha=0.5, linewidth=1.5, label="1:1 line", zorder=0)

        ax.set_xlabel("CP Emulator Features")
        ax.set_ylabel("CP Measure Features")

        mean_r = cluster_data["pearson_r"].mean()
        n_pairs = len(cluster_data)
        ax.set_title(
            f"Log-Scale Correlation: {cluster_name}\n(Mean r = {mean_r:.3f}, n = {n_pairs} feature pairs)",
            fontweight="bold",
        )

        # Legend - position based on number of entries
        if len(legend_added) <= 8:
            ax.legend(loc="upper left", fontsize=9, markerscale=1.5)
        else:
            ax.legend(loc="upper left", fontsize=7, markerscale=1.2, ncol=1)

        save_figure(fig, scatter_dir / f"cluster_scatter_{cluster_name}.png")
        plt.close()

    print(f"  Saved scatter plots for {len(clusters)} clusters to {scatter_dir}")


def parse_args():
    parser = argparse.ArgumentParser(description="Feature extraction benchmarks")
    parser.add_argument(
        "--mode",
        type=str,
        choices=["comparison", "timing"],
        default="comparison",
        help="Benchmark mode: 'comparison' compares pipelines, 'timing' times individual functions",
    )
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Only regenerate plots from existing results",
    )
    return parser.parse_args()


def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024


# =============================================================================
# COMPARISON MODE - Compare CP Measure vs CP Multichannel pipelines
# =============================================================================


def setup_comparison_dirs():
    """Create directories for comparison mode."""
    dirs = [
        OUTPUT_DIR,
        OUTPUT_DIR / "images" / "aligned",
        OUTPUT_DIR / "images" / "nuclei",
        OUTPUT_DIR / "images" / "cells",
        OUTPUT_DIR / "images" / "cytoplasms",
        OUTPUT_DIR / "cp_measure",
        OUTPUT_DIR / "cp_multichannel",
    ]
    for d in dirs:
        d.mkdir(exist_ok=True, parents=True)


def prepare_test_images():
    """Find one complete set of images per plate and copy to benchmark directory."""
    aligned_dir = OUTPUT_DIR / "images" / "aligned"
    nuclei_dir = OUTPUT_DIR / "images" / "nuclei"
    cells_dir = OUTPUT_DIR / "images" / "cells"
    cytoplasm_dir = OUTPUT_DIR / "images" / "cytoplasms"

    pattern = str(SOURCE_DIR / "*__aligned.tiff")
    all_aligned_images = glob.glob(pattern)
    print(f"Found {len(all_aligned_images)} total aligned images")

    images_by_plate = {}
    for img_path in all_aligned_images:
        img_name = os.path.basename(img_path)
        if img_name.startswith("P-"):
            plate = img_name.split("_")[0].replace("P-", "")
            well = img_name.split("_")[1].replace("W-", "")
            tile = img_name.split("_")[2].replace("T-", "")

            if plate not in images_by_plate:
                images_by_plate[plate] = []

            base_name = f"P-{plate}_W-{well}_T-{tile}"
            nuclei_path = str(SOURCE_DIR / f"{base_name}__nuclei.tiff")
            cells_path = str(SOURCE_DIR / f"{base_name}__cells.tiff")
            cytoplasm_path = str(
                SOURCE_DIR / f"{base_name}__identified_cytoplasms.tiff"
            )

            if all(
                os.path.exists(p) for p in [nuclei_path, cells_path, cytoplasm_path]
            ):
                images_by_plate[plate].append(
                    {
                        "base_name": base_name,
                        "aligned_path": img_path,
                        "nuclei_path": nuclei_path,
                        "cells_path": cells_path,
                        "cytoplasm_path": cytoplasm_path,
                    }
                )

    selected_images = []
    for plate, images in images_by_plate.items():
        if images:
            image_info = random.choice(images)
            base_name = image_info["base_name"]

            dst_aligned = aligned_dir / f"{base_name}__aligned.tiff"
            dst_nuclei = nuclei_dir / f"{base_name}__nuclei.tiff"
            dst_cells = cells_dir / f"{base_name}__cells.tiff"
            dst_cytoplasm = cytoplasm_dir / f"{base_name}__identified_cytoplasms.tiff"

            print(f"Copying files for plate {plate}: {base_name}")

            for src, dst in [
                (image_info["aligned_path"], dst_aligned),
                (image_info["nuclei_path"], dst_nuclei),
                (image_info["cells_path"], dst_cells),
                (image_info["cytoplasm_path"], dst_cytoplasm),
            ]:
                if not dst.exists():
                    os.system(f"cp {src} {dst}")

            selected_images.append(
                {
                    "base_name": base_name,
                    "aligned_path": dst_aligned,
                    "nuclei_path": dst_nuclei,
                    "cells_path": dst_cells,
                    "cytoplasm_path": dst_cytoplasm,
                }
            )

    print(f"\nCopied {len(selected_images)} complete sets to benchmark directories")
    return selected_images


def run_comparison_benchmark(images):
    """Compare CP Measure vs CP Multichannel on test images."""
    from lib.phenotype.extract_phenotype_cp_measure import extract_phenotype_cp_measure
    from lib.phenotype.extract_phenotype_cp_multichannel import (
        extract_phenotype_cp_multichannel,
    )

    all_stats = []
    print(f"\nRunning comparison benchmark on {len(images)} images")

    for image_info in tqdm(images, desc="Processing images"):
        base_name = image_info["base_name"]
        plate = base_name.split("_")[0].replace("P-", "")
        well = base_name.split("_")[1].replace("W-", "")
        tile = base_name.split("_")[2].replace("T-", "")

        print(f"\nProcessing {base_name}")

        try:
            data = imread(image_info["aligned_path"])
            nuclei = imread(image_info["nuclei_path"])
            cells = imread(image_info["cells_path"])
            cytoplasms = imread(image_info["cytoplasm_path"])
        except Exception as e:
            print(f"  Error loading: {e}")
            continue

        # CP Measure
        print("  Running CP Measure...")
        try:
            gc.collect()
            start_mem = get_memory_usage()
            start_time = time.time()

            features = extract_phenotype_cp_measure(
                data_phenotype=data,
                nuclei=nuclei,
                cells=cells,
                cytoplasms=cytoplasms,
                channel_names=CHANNEL_NAMES,
            )

            elapsed = time.time() - start_time
            mem = max(0, get_memory_usage() - start_mem)

            features.to_csv(
                OUTPUT_DIR / "cp_measure" / f"{base_name}__features.csv", index=False
            )
            all_stats.append(
                {
                    "method": "cp_measure",
                    "plate": plate,
                    "well": well,
                    "tile": tile,
                    "feature_count": len(features.columns),
                    "row_count": len(features),
                    "runtime_seconds": elapsed,
                    "memory_mb": mem,
                    "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                }
            )
            del features
            gc.collect()
        except Exception as e:
            print(f"  CP Measure error: {e}")

        # CP Multichannel
        print("  Running CP Multichannel...")
        try:
            gc.collect()
            start_mem = get_memory_usage()
            start_time = time.time()

            features = extract_phenotype_cp_multichannel(
                data_phenotype=data,
                nuclei=nuclei,
                cells=cells,
                cytoplasms=cytoplasms,
                channel_names=CHANNEL_NAMES,
                wildcards={"plate": plate, "well": well, "tile": tile},
            )

            elapsed = time.time() - start_time
            mem = max(0, get_memory_usage() - start_mem)

            features.to_csv(
                OUTPUT_DIR / "cp_multichannel" / f"{base_name}__features.csv",
                index=False,
            )
            all_stats.append(
                {
                    "method": "cp_multichannel",
                    "plate": plate,
                    "well": well,
                    "tile": tile,
                    "feature_count": len(features.columns),
                    "row_count": len(features),
                    "runtime_seconds": elapsed,
                    "memory_mb": mem,
                    "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                }
            )
            del features
            gc.collect()
        except Exception as e:
            print(f"  CP Multichannel error: {e}")

        del data, nuclei, cells, cytoplasms
        gc.collect()

    if all_stats:
        results = pd.DataFrame(all_stats)
        results.to_csv(OUTPUT_DIR / "benchmark_results.csv", index=False)
        generate_comparison_summary(results)
        return results
    return None


def generate_comparison_summary(results_df, run_correlation=True):
    """Generate summary for comparison benchmark."""
    from plot_style import setup_plot_style, COLORS

    # Apply consistent plot styling
    setup_plot_style()

    summary = (
        results_df.groupby("method")
        .agg(
            {
                "feature_count": ["mean", "std"],
                "row_count": ["mean", "std"],
                "runtime_seconds": ["mean", "std"],
                "memory_mb": ["mean", "std"],
            }
        )
        .reset_index()
    )
    summary.to_csv(OUTPUT_DIR / "benchmark_summary.csv")

    # Calculate speedup
    cp_measure_time = results_df[results_df["method"] == "cp_measure"][
        "runtime_seconds"
    ].mean()
    cp_multi_time = results_df[results_df["method"] == "cp_multichannel"][
        "runtime_seconds"
    ].mean()
    speedup = cp_measure_time / cp_multi_time if cp_multi_time > 0 else 0

    # 1. IMPROVED RUNTIME COMPARISON PLOT
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Left: Bar chart with actual values
    ax = axes[0]
    method_colors = {
        "cp_measure": COLORS.get("cp_measure", "#8c564b"),
        "cp_multichannel": COLORS.get("cp_multichannel", "#e377c2"),
    }
    bars = ax.bar(
        ["CP Measure", "CP Multichannel"],
        [cp_measure_time, cp_multi_time],
        color=[method_colors["cp_measure"], method_colors["cp_multichannel"]],
        edgecolor="black",
        linewidth=1,
    )

    # Add value labels
    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{height:.1f}s",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=12,
        )

    ax.set_ylabel("Runtime (seconds)", fontsize=12)
    ax.set_title("Runtime Comparison", fontsize=14, fontweight="bold")
    ax.set_ylim(0, cp_measure_time * 1.15)

    # Right: Speedup visualization
    ax2 = axes[1]
    ax2.barh(
        ["Speedup Factor"],
        [speedup],
        color="forestgreen",
        edgecolor="black",
        height=0.5,
    )
    ax2.axvline(x=1, color="red", linestyle="--", linewidth=2, label="No speedup (1x)")
    ax2.set_xlabel("Speedup (x times faster)", fontsize=12)
    ax2.set_title(
        f"CP Multichannel is {speedup:.0f}x Faster", fontsize=14, fontweight="bold"
    )
    ax2.set_xlim(0, speedup * 1.15)
    ax2.annotate(
        f"{speedup:.0f}x",
        xy=(speedup, 0),
        xytext=(5, 0),
        textcoords="offset points",
        ha="left",
        va="center",
        fontweight="bold",
        fontsize=16,
        color="forestgreen",
    )

    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "runtime_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 2. MEMORY COMPARISON PLOT
    cp_measure_mem = results_df[results_df["method"] == "cp_measure"][
        "memory_mb"
    ].mean()
    cp_multi_mem = results_df[results_df["method"] == "cp_multichannel"][
        "memory_mb"
    ].mean()

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(
        ["CP Measure", "CP Multichannel"],
        [cp_measure_mem, cp_multi_mem],
        color=[method_colors["cp_measure"], method_colors["cp_multichannel"]],
        edgecolor="black",
        linewidth=1,
    )

    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{height:.1f} MB",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=11,
        )

    ax.set_ylabel("Memory Usage (MB)", fontsize=12)
    ax.set_title("Memory Usage Comparison", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "memory_comparison.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 3. FEATURE COUNT COMPARISON
    cp_measure_feat = results_df[results_df["method"] == "cp_measure"][
        "feature_count"
    ].mean()
    cp_multi_feat = results_df[results_df["method"] == "cp_multichannel"][
        "feature_count"
    ].mean()

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(
        ["CP Measure", "CP Multichannel"],
        [cp_measure_feat, cp_multi_feat],
        color=[method_colors["cp_measure"], method_colors["cp_multichannel"]],
        edgecolor="black",
        linewidth=1,
    )

    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{int(height):,}",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 5),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontweight="bold",
            fontsize=11,
        )

    ax.set_ylabel("Number of Features", fontsize=12)
    ax.set_title("Feature Count Comparison", fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / "feature_count_comparison.png", dpi=300, bbox_inches="tight"
    )
    plt.close()

    # Print summary
    print("\n" + "=" * 60)
    print("COMPARISON SUMMARY")
    print("=" * 60)
    for method in results_df["method"].unique():
        m = results_df[results_df["method"] == method]
        print(f"\n{method.upper()}")
        print(
            f"  Runtime: {m['runtime_seconds'].mean():.2f}s (+/- {m['runtime_seconds'].std():.2f})"
        )
        print(f"  Memory: {m['memory_mb'].mean():.1f} MB")
        print(f"  Features: {m['feature_count'].mean():.0f}")

    print(f"\nSPEEDUP: CP Multichannel is {speedup:.1f}x faster than CP Measure")

    # 4. RUN CORRELATION ANALYSIS
    if run_correlation:
        print("\n" + "=" * 60)
        print("RUNNING CORRELATION ANALYSIS")
        print("=" * 60)

        measure_df, multi_df = load_feature_dataframes()

        if measure_df is not None and multi_df is not None:
            # Get feature mappings
            mappings = get_feature_mappings()
            total_mappings = sum(len(v) for v in mappings.values())
            print(
                f"\nDefined {total_mappings} feature mappings across {len(mappings)} clusters"
            )

            # Compute correlations
            print("\nComputing feature correlations...")
            corr_df = compute_feature_correlations(measure_df, multi_df, mappings)

            # Save correlation results
            corr_df.to_csv(OUTPUT_DIR / "feature_correlations.csv", index=False)
            print("Saved correlation results to feature_correlations.csv")

            # Generate correlation summary plots
            generate_correlation_plots(corr_df)

            # Generate scatter plots for each cluster
            generate_scatter_plots(corr_df, measure_df, multi_df)
        else:
            print("Could not load feature dataframes for correlation analysis")


# =============================================================================
# TIMING MODE - Time individual cp_measure functions
# =============================================================================


def find_random_image():
    """Find one random complete image set."""
    aligned_files = list(SOURCE_DIR.glob("*__aligned.tiff"))
    random.shuffle(aligned_files)

    for aligned_path in aligned_files:
        base_name = aligned_path.name.replace("__aligned.tiff", "")
        nuclei_path = SOURCE_DIR / f"{base_name}__nuclei.tiff"
        cells_path = SOURCE_DIR / f"{base_name}__cells.tiff"

        if nuclei_path.exists() and cells_path.exists():
            return {
                "base_name": base_name,
                "aligned": aligned_path,
                "nuclei": nuclei_path,
                "cells": cells_path,
            }

    raise RuntimeError("No complete image sets found!")


def time_function(name, func, args):
    """Time a function and return results."""
    print(f"  {name}...", end=" ", flush=True)
    start = time.time()
    try:
        result = func(*args)
        elapsed = time.time() - start
        feature_count = len(result) if isinstance(result, dict) else 0
        print(f"{elapsed:.2f}s ({feature_count} features)")
        return {
            "name": name,
            "status": "SUCCESS",
            "time_seconds": elapsed,
            "feature_count": feature_count,
            "error": None,
        }
    except Exception as e:
        elapsed = time.time() - start
        print(f"ERROR: {str(e)[:50]}")
        return {
            "name": name,
            "status": "ERROR",
            "time_seconds": elapsed,
            "feature_count": 0,
            "error": str(e),
        }


def run_timing_benchmark():
    """Time individual cp_measure functions."""
    # Import cp_measure functions
    from cp_measure.core.measureobjectsizeshape import (
        get_sizeshape,
        get_zernike,
        get_ferret,
    )
    from cp_measure.core.measureobjectintensity import get_intensity
    from cp_measure.core.measureobjectintensitydistribution import get_radial_zernikes
    from cp_measure.core.measuregranularity import get_granularity
    from cp_measure.core.measuretexture import get_texture
    from cp_measure.core.measurecolocalization import (
        get_correlation_pearson,
        get_correlation_manders_fold,
        get_correlation_rwc,
        get_correlation_costes,
        get_correlation_overlap,
    )
    from cp_measure.multimask.measureobjectoverlap import measureobjectoverlap
    from cp_measure.multimask.measureobjectneighbors import measureobjectneighbors

    print("\nLoading test image...")
    image_set = find_random_image()
    print(f"Using: {image_set['base_name']}")

    aligned = imread(image_set["aligned"])
    nuclei = imread(image_set["nuclei"])
    cells = imread(image_set["cells"])

    n_nuclei = len(np.unique(nuclei)) - 1
    n_cells = len(np.unique(cells)) - 1
    print(f"  Shape: {aligned.shape}, Nuclei: {n_nuclei}, Cells: {n_cells}")

    ch0 = aligned[0]
    ch1 = aligned[1] if aligned.shape[0] > 1 else aligned[0]

    results = []

    # Core measurements
    print("\n--- CORE MEASUREMENTS ---")
    for name, func, args in [
        ("get_sizeshape", get_sizeshape, (nuclei, ch0)),
        ("get_zernike", get_zernike, (nuclei, ch0)),
        ("get_ferret", get_ferret, (nuclei, ch0)),
        ("get_intensity", get_intensity, (nuclei, ch0)),
        ("get_radial_zernikes", get_radial_zernikes, (nuclei, ch0)),
        ("get_granularity", get_granularity, (nuclei, ch0)),
        ("get_texture", get_texture, (nuclei, ch0)),
    ]:
        result = time_function(name, func, args)
        result["category"] = "core"
        results.append(result)

    # Correlation measurements
    print("\n--- CORRELATION MEASUREMENTS ---")
    for name, func, args in [
        ("get_correlation_pearson", get_correlation_pearson, (ch0, ch1, nuclei)),
        (
            "get_correlation_manders_fold",
            get_correlation_manders_fold,
            (ch0, ch1, nuclei),
        ),
        ("get_correlation_rwc", get_correlation_rwc, (ch0, ch1, nuclei)),
        ("get_correlation_costes", get_correlation_costes, (ch0, ch1, nuclei)),
        ("get_correlation_overlap", get_correlation_overlap, (ch0, ch1, nuclei)),
    ]:
        result = time_function(name, func, args)
        result["category"] = "correlation"
        results.append(result)

    # Multimask measurements
    print("\n--- MULTIMASK MEASUREMENTS ---")
    for name, func, args in [
        ("measureobjectoverlap", measureobjectoverlap, (nuclei, nuclei)),
        ("measureobjectneighbors", measureobjectneighbors, (nuclei, nuclei)),
    ]:
        result = time_function(name, func, args)
        result["category"] = "multimask"
        results.append(result)

    # Save results
    df = pd.DataFrame(results)
    df["timestamp"] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    df["image"] = image_set["base_name"]
    df["n_objects"] = n_nuclei

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_file = OUTPUT_DIR / "timing_results.csv"
    df.to_csv(output_file, index=False)

    # Summary
    print("\n" + "=" * 50)
    print("TIMING SUMMARY")
    print("=" * 50)

    successful = df[df["status"] == "SUCCESS"]
    print(f"\nTotal time: {successful['time_seconds'].sum():.2f}s")
    print(f"Total features: {successful['feature_count'].sum()}")
    print("\nBy function (sorted by time):")

    for _, row in df.sort_values("time_seconds", ascending=False).iterrows():
        status = (
            f"{row['time_seconds']:.2f}s" if row["status"] == "SUCCESS" else "ERROR"
        )
        print(f"  {row['name']}: {status}")

    print(f"\nResults saved to: {output_file}")
    return df


# =============================================================================
# MAIN
# =============================================================================


def main():
    args = parse_args()

    print("=" * 60)
    print(f"FEATURE EXTRACTION BENCHMARK - {args.mode.upper()} MODE")
    print("=" * 60)

    if args.plots_only:
        # Load existing results and regenerate plots
        results_path = OUTPUT_DIR / "benchmark_results.csv"
        if results_path.exists():
            print("Regenerating plots from existing results...")
            results_df = pd.read_csv(results_path)
            generate_comparison_summary(results_df)
        else:
            print(f"No results found at {results_path}. Run benchmark first.")
    elif args.mode == "comparison":
        setup_comparison_dirs()
        images = prepare_test_images()
        run_comparison_benchmark(images)
    else:  # timing
        run_timing_benchmark()

    print("\nBenchmark completed!")


if __name__ == "__main__":
    main()
