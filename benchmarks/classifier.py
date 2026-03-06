"""
Classifier Benchmark

Generates two panels:
1. Confusion matrix from the trained XGBoost classifier
2. Training data composition: cell counts by class across plates

Usage:
    python classifier.py
    python classifier.py --plots-only
"""

import sys
import json
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Import shared plot styling
BENCHMARKS_DIR = Path(__file__).parent.resolve()
sys.path.insert(0, str(BENCHMARKS_DIR))
from plot_style import setup_plot_style, save_figure

# Configure paths
PROJECT_ROOT = BENCHMARKS_DIR.parent
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
BRIEFLOW_OUTPUT = ANALYSIS_DIR / "brieflow_output"
OUTPUT_DIR = BENCHMARKS_DIR / "results" / "classifier"
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Classifier paths
RUN_DIR = BRIEFLOW_OUTPUT / "classifier" / "classifier" / "run_20260112_172222"
MODEL_DIR = RUN_DIR / "models" / "multiclass_xgb_k50"
REPORT_PATH = MODEL_DIR / "multiclass_xgb_k50_classification_report.json"
MODEL_PATH = MODEL_DIR / "multiclass_xgb_k50_model.dill"
TRAINING_DATA_PATH = (
    BRIEFLOW_OUTPUT
    / "classifier"
    / "training_dataset"
    / "cell_classifier_training_dataset_for_cell_stage_20260112_172158.parquet"
)

# Class mapping (from config.yml)
CLASS_MAPPING = {1: "Mitotic", 2: "Interphase"}


def parse_args():
    parser = argparse.ArgumentParser(description="Classifier benchmark plots")
    parser.add_argument(
        "--plots-only",
        action="store_true",
        help="Only regenerate plots (same behavior, kept for consistency)",
    )
    return parser.parse_args()


def plot_confusion_matrix(output_path):
    """Generate confusion matrix from classification report."""
    setup_plot_style()

    with open(REPORT_PATH) as f:
        report = json.load(f)

    class_names = ["Mitotic", "Interphase"]

    # Reconstruct confusion matrix from report
    # Class 0 = Mitotic (label 1), Class 1 = Interphase (label 2)
    mitotic_support = int(report["0"]["support"])
    interphase_support = int(report["1"]["support"])
    mitotic_recall = report["0"]["recall"]
    interphase_recall = report["1"]["recall"]

    # TP and FN for each class
    mitotic_tp = round(mitotic_recall * mitotic_support)
    mitotic_fn = mitotic_support - mitotic_tp
    interphase_tp = round(interphase_recall * interphase_support)
    interphase_fn = interphase_support - interphase_tp

    # Confusion matrix: rows = true, cols = predicted
    # [Mitotic predicted Mitotic, Mitotic predicted Interphase]
    # [Interphase predicted Mitotic, Interphase predicted Interphase]
    cm = np.array(
        [
            [mitotic_tp, mitotic_fn],
            [interphase_fn, interphase_tp],
        ]
    )

    # Normalize for display
    cm_normalized = cm.astype(float) / cm.sum(axis=1, keepdims=True)

    fig, ax = plt.subplots(figsize=(6, 6))

    # Plot heatmap
    im = ax.imshow(cm_normalized, interpolation="nearest", cmap="Blues", vmin=0, vmax=1)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Proportion", fontsize=12)

    # Add text annotations
    for i in range(2):
        for j in range(2):
            color = "white" if cm_normalized[i, j] > 0.5 else "black"
            ax.text(
                j,
                i,
                f"{cm[i, j]}\n({cm_normalized[i, j]:.1%})",
                ha="center",
                va="center",
                fontsize=14,
                fontweight="bold",
                color=color,
            )

    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(class_names, fontsize=13)
    ax.set_yticklabels(class_names, fontsize=13)
    ax.set_xlabel("Predicted", fontsize=14)
    ax.set_ylabel("True", fontsize=14)
    ax.set_title(
        f"Confusion Matrix (Accuracy: {report['accuracy']:.1%})",
        fontweight="bold",
        fontsize=14,
    )

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_training_composition(output_path):
    """Plot training data composition by class and plate."""
    setup_plot_style()

    df = pd.read_parquet(TRAINING_DATA_PATH, columns=["plate", "cell_stage"])
    df["class_name"] = df["cell_stage"].map(CLASS_MAPPING)

    # Count by plate and class
    counts = (
        df.groupby(["plate", "class_name"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=["Mitotic", "Interphase"], fill_value=0)
    )
    counts.index = [f"P-{p}" for p in counts.index]

    fig, ax = plt.subplots(figsize=(8, 6))

    x = np.arange(len(counts))
    width = 0.35

    mitotic_vals = counts["Mitotic"].values
    interphase_vals = counts["Interphase"].values

    bars1 = ax.bar(
        x - width / 2,
        mitotic_vals,
        width,
        label="Mitotic",
        color="#8856a7",
        edgecolor="white",
        linewidth=0.5,
    )
    bars2 = ax.bar(
        x + width / 2,
        interphase_vals,
        width,
        label="Interphase",
        color="#43a2ca",
        edgecolor="white",
        linewidth=0.5,
    )

    # Add value labels
    for bar in bars1:
        val = bar.get_height()
        if val > 0:
            ax.annotate(
                f"{int(val)}",
                xy=(bar.get_x() + bar.get_width() / 2, val),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
            )
    for bar in bars2:
        val = bar.get_height()
        if val > 0:
            ax.annotate(
                f"{int(val)}",
                xy=(bar.get_x() + bar.get_width() / 2, val),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=9,
                fontweight="bold",
            )

    ax.set_xlabel("Plate", fontsize=12)
    ax.set_ylabel("Number of Cells", fontsize=12)
    ax.set_title(
        f"Training Data Composition (n={len(df)})",
        fontweight="bold",
        fontsize=14,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(counts.index, fontsize=11)
    ax.legend(fontsize=11)

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def plot_training_totals(output_path):
    """Plot total training cells by class."""
    setup_plot_style()

    df = pd.read_parquet(TRAINING_DATA_PATH, columns=["cell_stage"])
    df["class_name"] = df["cell_stage"].map(CLASS_MAPPING)
    totals = df["class_name"].value_counts().reindex(["Mitotic", "Interphase"])

    fig, ax = plt.subplots(figsize=(5, 6))

    colors = ["#8856a7", "#43a2ca"]
    bars = ax.bar(
        totals.index,
        totals.values,
        color=colors,
        edgecolor="white",
        linewidth=0.5,
        width=0.5,
    )

    for bar, val in zip(bars, totals.values):
        ax.annotate(
            f"{int(val)}",
            xy=(bar.get_x() + bar.get_width() / 2, val),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=13,
            fontweight="bold",
        )

    ax.set_ylabel("Number of Cells", fontsize=12)
    ax.set_title(
        f"Total Training Cells (n={len(df)})",
        fontweight="bold",
        fontsize=14,
    )
    ax.tick_params(axis="x", labelsize=13)

    plt.tight_layout()
    save_figure(fig, output_path)
    plt.close()
    print(f"Saved: {output_path}")


def main():
    parse_args()

    print("=" * 60)
    print("CLASSIFIER BENCHMARK")
    print("=" * 60)

    print("\nGenerating confusion matrix...")
    plot_confusion_matrix(OUTPUT_DIR / "confusion_matrix.png")

    print("Generating training composition plot...")
    plot_training_composition(OUTPUT_DIR / "training_composition.png")

    print("Generating training totals plot...")
    plot_training_totals(OUTPUT_DIR / "training_totals.png")

    print("\n" + "=" * 60)
    print("BENCHMARK COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
