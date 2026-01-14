#!/usr/bin/env python3
"""
Validate alignment strategy across plates and wells.

Generates 48 alignment quality plots:
- 8 plates
- 6 wells per plate (A1, A2, A3, B1, B2, B3)
- 1 random tile per plate (same tile used for all wells in that plate)

Total: 8 × 6 = 48 plots
"""

import numpy as np
import yaml
from pathlib import Path
from tifffile import imread
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from lib.shared.file_utils import get_filename
from lib.shared.illumination_correction import apply_ic_field
from lib.phenotype.align_channels import align_phenotype_channels
from lib.shared.align import apply_custom_offsets


# Configuration
CONFIG_FILE = "config/config.yml"
OUTPUT_DIR = Path("alignment_validation")
OUTPUT_DIR.mkdir(exist_ok=True)

PLATES = [1, 2, 3, 4, 5, 6, 7, 8]
WELLS = ["A1", "A2", "A3", "B1", "B2", "B3"]
TILE_RANGE = (1, 1280)  # Inclusive

# Alignment configurations per plate
# Strategy: Custom offsets only (no automatic DAPI-to-GH2AX), with optional GFP alignment
ALIGNMENTS = {
    1: {
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-5, 12), 3: (-5, 12)}
    },
    2: {
        "steps": [
            {
                "target": 2,
                "source": 3,
                "riders": [4],
                "remove_channel": "source"
            }
        ],
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-5, 12), 3: (-5, 12)}
    },
    3: {
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-15, 15), 3: (-15, 15)}
    },
    4: {
        "steps": [
            {
                "target": 1,
                "source": 4,
                "riders": [],
                "remove_channel": "target"
            }
        ],
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (0, 7), 3: (0, 7)}
    },
    5: {
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-5, 12), 3: (-5, 12)}
    },
    6: {
        "steps": [
            {
                "target": 1,
                "source": 4,
                "riders": [],
                "remove_channel": "target"
            }
        ],
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (0, 12), 3: (0, 12)}
    },
    7: {
        "steps": [
            {
                "target": 1,
                "source": 4,
                "riders": [],
                "remove_channel": "target"
            }
        ],
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-5, 12), 3: (-5, 12)}
    },
    8: {
        "steps": [
            {
                "target": 1,
                "source": 4,
                "riders": [],
                "remove_channel": "target"
            }
        ],
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-12, 12), 3: (-12, 12)}
    }
}

# Alternative: Full two-round configuration
# Strategy: Custom offsets FIRST, then automatic alignment (DAPI-to-GH2AX + GFP)
ALIGNMENTS_FULL_TWO_ROUND = {
    1: {
        "target": 0,
        "source": 2,
        "riders": [3],
        "remove_channel": False,
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-5, 10), 3: (-5, 10)}
    },
    2: {
        "steps": [
            {
                "target": 0,
                "source": 2,
                "riders": [],
                "remove_channel": False
            },
            {
                "target": 2,
                "source": 3,
                "riders": [4],
                "remove_channel": "source"
            }
        ],
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-5, 12), 3: (-5, 12)}
    },
    3: {
        "target": 0,
        "source": 2,
        "riders": [3],
        "remove_channel": False,
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-12, 15), 3: (-12, 15)}
    },
    4: {
        "steps": [
            {
                "target": 0,
                "source": 2,
                "riders": [3],
                "remove_channel": False
            },
            {
                "target": 1,
                "source": 4,
                "riders": [],
                "remove_channel": "target"
            }
        ],
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (0, 7), 3: (0, 7)}
    },
    5: {
        "target": 0,
        "source": 2,
        "riders": [3],
        "remove_channel": False,
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-5, 12), 3: (-5, 12)}
    },
    6: {
        "steps": [
            {
                "target": 0,
                "source": 2,
                "riders": [3],
                "remove_channel": False
            },
            {
                "target": 1,
                "source": 4,
                "riders": [],
                "remove_channel": "target"
            }
        ],
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (0, 15), 3: (0, 15)}
    },
    7: {
        "steps": [
            {
                "target": 0,
                "source": 2,
                "riders": [3],
                "remove_channel": False
            },
            {
                "target": 1,
                "source": 4,
                "riders": [],
                "remove_channel": "target"
            }
        ],
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-5, 12), 3: (-5, 12)}
    },
    8: {
        "steps": [
            {
                "target": 0,
                "source": 2,
                "riders": [3],
                "remove_channel": False
            },
            {
                "target": 1,
                "source": 4,
                "riders": [],
                "remove_channel": "target"
            }
        ],
        "window": 2,
        "upsample_factor": 4,
        "custom_channel_offsets": {2: (-12, 12), 3: (-12, 12)}
    }
}

# Select which configuration to use
USE_FULL_TWO_ROUND = True  # Set to False to use ALIGNMENTS instead
ACTIVE_ALIGNMENTS = ALIGNMENTS_FULL_TWO_ROUND if USE_FULL_TWO_ROUND else ALIGNMENTS

# Visualization channels (DAPI, TUBULIN, GH2AX, PHALLOIDIN)
VIZ_CHANNELS = [0, 1, 2, 3]  # Channel indices
CROP_SIZE = 300


def apply_alignment(image, plate_config):
    """Apply alignment strategy: custom offsets first, then automatic alignment."""
    aligned = image.copy()

    # Step 1: Apply custom offsets (if defined)
    custom_offsets = plate_config.get("custom_channel_offsets", None)
    if custom_offsets:
        aligned = apply_custom_offsets(aligned, custom_offsets)
        print(f"    Applied custom offsets: {custom_offsets}")

    # Step 2: Apply automatic alignment (if steps defined)
    if "steps" in plate_config:
        window = plate_config.get("window", 2)
        upsample = plate_config.get("upsample_factor", 2)

        for step_num, step in enumerate(plate_config["steps"], 1):
            print(f"    Step {step_num}: target={step['target']}, source={step['source']}")
            aligned = align_phenotype_channels(
                aligned,
                target=step["target"],
                source=step["source"],
                riders=step["riders"],
                remove_channel=step["remove_channel"],
                upsample_factor=upsample,
                window=window,
                verbose=False,
            )

    return aligned


def visualize_alignment(image, channels, crop_size=300):
    """Create 16-location RGB composite visualization for alignment QC.

    First channel shown in grayscale, remaining 3 as RGB composite.
    Color fringing indicates misalignment.
    """
    _, height, width = image.shape

    # Define 16 crop locations (4x4 grid)
    margin = 50
    mid_y = (height - crop_size) // 2
    mid_x = (width - crop_size) // 2

    np.random.seed(42)  # For reproducible random locations
    locations = [
        # Row 1
        ("Top-Left", margin, margin),
        ("Top-Mid", margin, mid_x),
        ("Top-Right", margin, width - crop_size - margin),
        ("Rand-1", np.random.randint(margin, height - crop_size - margin),
         np.random.randint(margin, width - crop_size - margin)),
        # Row 2
        ("Left-Mid", mid_y, margin),
        ("TL-Quad", mid_y // 2, mid_x // 2),
        ("TR-Quad", mid_y // 2, mid_x + mid_x // 2),
        ("Rand-2", np.random.randint(margin, height - crop_size - margin),
         np.random.randint(margin, width - crop_size - margin)),
        # Row 3
        ("BL-Quad", mid_y + mid_y // 2, mid_x // 2),
        ("Center", mid_y, mid_x),
        ("BR-Quad", mid_y + mid_y // 2, mid_x + mid_x // 2),
        ("Right-Mid", mid_y, width - crop_size - margin),
        # Row 4
        ("Rand-3", np.random.randint(margin, height - crop_size - margin),
         np.random.randint(margin, width - crop_size - margin)),
        ("Bot-Left", height - crop_size - margin, margin),
        ("Bot-Mid", height - crop_size - margin, mid_x),
        ("Bot-Right", height - crop_size - margin, width - crop_size - margin),
    ]

    # Create figure
    fig = plt.figure(figsize=(20, 20))
    gs = GridSpec(4, 4, figure=fig, hspace=0.3, wspace=0.2)

    for idx, (loc_name, y_start, x_start) in enumerate(locations):
        y_end = y_start + crop_size
        x_end = x_start + crop_size

        # Grayscale base (first channel)
        gray_crop = image[channels[0], y_start:y_end, x_start:x_end]
        p2, p98 = np.percentile(gray_crop, [2, 98])
        gray_norm = np.clip((gray_crop - p2) / (p98 - p2 + 1e-8), 0, 1)

        # RGB overlay (channels 2-4)
        rgb = np.zeros((crop_size, crop_size, 3))
        for i, ch_idx in enumerate(channels[1:]):
            crop = image[ch_idx, y_start:y_end, x_start:x_end]
            p2, p98 = np.percentile(crop, [2, 98])
            crop_norm = np.clip((crop - p2) / (p98 - p2 + 1e-8), 0, 1)
            rgb[:, :, i] = crop_norm

        # Blend: 50% grayscale, 50% RGB
        rgba = np.zeros((crop_size, crop_size, 3))
        for i in range(3):
            rgba[:, :, i] = 0.5 * gray_norm + 0.5 * rgb[:, :, i]

        # Plot
        ax = fig.add_subplot(gs[idx // 4, idx % 4])
        ax.imshow(rgba)
        ax.set_title(f"{loc_name}", fontsize=9)
        ax.axis("off")

    return fig


def main():
    # Load config
    with open(CONFIG_FILE, "r") as f:
        config = yaml.safe_load(f)

    root_fp = Path(config["all"]["root_fp"])
    preprocess_fp = root_fp / "preprocess"

    # Print strategy
    strategy_name = "Full Two-Round (Custom + Automatic)" if USE_FULL_TWO_ROUND else "Custom Only"
    print(f"{'='*60}")
    print(f"Alignment Strategy: {strategy_name}")
    if USE_FULL_TWO_ROUND:
        print(f"  1. Apply custom offsets")
        print(f"  2. Run automatic DAPI-to-GH2AX alignment")
        print(f"  3. Run automatic GFP alignment (if applicable)")
    else:
        print(f"  1. Apply custom offsets")
        print(f"  2. Run automatic GFP alignment only (if applicable)")
    print(f"{'='*60}\n")

    # Generate one random tile per plate-well combination
    np.random.seed(123)  # Reproducible random tiles
    random_tiles = {}
    for plate in PLATES:
        for well in WELLS:
            random_tiles[(plate, well)] = np.random.randint(TILE_RANGE[0], TILE_RANGE[1] + 1)

    print(f"Processing {len(PLATES)} plates × {len(WELLS)} wells = {len(PLATES) * len(WELLS)} images")
    print(f"Each plate-well gets a unique random tile from range {TILE_RANGE}\n")

    results = []

    for plate in PLATES:
        plate_config = ACTIVE_ALIGNMENTS[plate]
        print(f"=== PLATE {plate} ===")

        for well in WELLS:
            tile = random_tiles[(plate, well)]
            print(f"  Processing {well}...")

            try:
                # Load raw image
                img_path = preprocess_fp / "images" / "phenotype" / get_filename(
                    {"plate": plate, "well": well, "tile": tile}, "image", "tiff"
                )

                if not img_path.exists():
                    print(f"    WARNING: Image not found: {img_path}")
                    continue

                raw_image = imread(str(img_path))

                # Load illumination correction
                ic_path = preprocess_fp / "ic_fields" / "phenotype" / get_filename(
                    {"plate": plate, "well": well}, "ic_field", "tiff"
                )

                if not ic_path.exists():
                    print(f"    WARNING: IC field not found: {ic_path}")
                    continue

                ic_field = imread(str(ic_path))

                # Apply illumination correction
                corrected = apply_ic_field(raw_image, correction=ic_field)

                # Apply alignment strategy
                aligned = apply_alignment(corrected, plate_config)

                # Generate visualization
                fig = visualize_alignment(aligned, VIZ_CHANNELS, crop_size=CROP_SIZE)

                # Save figure
                output_file = OUTPUT_DIR / f"plate{plate:02d}_{well}_tile{tile:04d}.png"
                fig.suptitle(
                    f"Plate {plate} | Well {well} | Tile {tile}\n"
                    f"Channels: DAPI (gray) + TUBULIN (R), GH2AX (G), PHALLOIDIN (B)",
                    fontsize=14, y=0.995
                )
                plt.savefig(output_file, dpi=150, bbox_inches='tight')
                plt.close(fig)

                print(f"    Saved: {output_file}")
                results.append((plate, well, tile, "success"))

            except Exception as e:
                print(f"    ERROR: {e}")
                results.append((plate, well, tile, f"error: {e}"))

    # Summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    success = sum(1 for r in results if r[3] == "success")
    print(f"Successfully processed: {success}/{len(results)} images")
    print(f"Output directory: {OUTPUT_DIR.absolute()}")

    # Save summary
    strategy_name = "Full Two-Round (Custom + Automatic)" if USE_FULL_TWO_ROUND else "Custom Only"
    summary_file = OUTPUT_DIR / "summary.txt"
    with open(summary_file, "w") as f:
        f.write("Alignment Validation Summary\n")
        f.write("="*60 + "\n\n")
        f.write(f"Strategy: {strategy_name}\n")
        f.write(f"  - Custom offsets applied first\n")
        if USE_FULL_TWO_ROUND:
            f.write(f"  - Then automatic DAPI-to-GH2AX + GFP alignment\n")
        else:
            f.write(f"  - Automatic alignment only for GFP (no DAPI-to-GH2AX)\n")
        f.write(f"\nPlates: {PLATES}\n")
        f.write(f"Wells: {WELLS}\n")
        f.write(f"Tile range: {TILE_RANGE}\n")
        f.write(f"One random tile per plate-well combination\n\n")
        f.write(f"Results:\n")
        for plate, well, tile, status in results:
            f.write(f"  Plate {plate}, {well}, Tile {tile}: {status}\n")
        f.write(f"\nTotal: {success}/{len(results)} successful\n")

    print(f"Summary saved: {summary_file}")


if __name__ == "__main__":
    main()
