"""
Shared plotting style configuration for all benchmarks.

Provides consistent, publication-quality plots with Helvetica font.
"""

import matplotlib.pyplot as plt
import matplotlib as mpl


def setup_plot_style():
    """Configure matplotlib for publication-quality plots with consistent styling.

    Call this at the start of any visualization function.
    """
    # Reset to default first
    plt.rcdefaults()

    # Font configuration - Helvetica-like fonts with fallbacks
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Nimbus Sans', 'Liberation Sans', 'Arial', 'DejaVu Sans'],
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11,
        'legend.fontsize': 11,
        'figure.titlesize': 16,
    })

    # Axes styling
    plt.rcParams.update({
        'axes.linewidth': 1.2,
        'axes.edgecolor': '#333333',
        'axes.labelcolor': '#333333',
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

    # Grid styling
    plt.rcParams.update({
        'axes.grid': False,
        'grid.alpha': 0.3,
        'grid.linewidth': 0.5,
    })

    # Tick styling
    plt.rcParams.update({
        'xtick.color': '#333333',
        'ytick.color': '#333333',
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'xtick.major.size': 4,
        'ytick.major.size': 4,
    })

    # Figure styling
    plt.rcParams.update({
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.facecolor': 'white',
        'savefig.edgecolor': 'white',
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.1,
    })

    # Legend styling
    plt.rcParams.update({
        'legend.frameon': True,
        'legend.framealpha': 0.9,
        'legend.edgecolor': '#cccccc',
        'legend.fancybox': False,
    })


# Standard color palettes for benchmarks
COLORS = {
    # Segmentation methods
    'cellpose': '#1f77b4',      # Blue
    'cellpose4': '#ff7f0e',     # Orange
    'stardist': '#2ca02c',      # Green

    # Spot calling methods
    'standard': '#d62728',      # Red
    'spotiflow': '#9467bd',     # Purple

    # Feature extraction methods
    'cp_measure': '#8c564b',    # Brown
    'cp_multichannel': '#e377c2',  # Pink

    # Merge methods
    'fast': '#1f77b4',          # Blue
    'stitch': '#ff7f0e',        # Orange
}

# Standard figure sizes
FIGSIZE = {
    'single': (6, 5),
    'wide': (10, 5),
    'square': (8, 8),
    'double': (12, 5),
    'quad': (10, 8),
}


def get_method_colors(methods, palette='default'):
    """Get consistent colors for a list of methods.

    Args:
        methods: List of method names
        palette: Color palette to use ('default' uses COLORS dict)

    Returns:
        List of colors matching the methods
    """
    if palette == 'default':
        return [COLORS.get(m.lower(), '#7f7f7f') for m in methods]
    else:
        import seaborn as sns
        return sns.color_palette(palette, len(methods))


def add_value_labels(ax, bars, fmt='{:.0f}', fontsize=10, offset=5):
    """Add value labels on top of bar chart bars.

    Args:
        ax: Matplotlib axes
        bars: Bar container from barplot
        fmt: Format string for values
        fontsize: Font size for labels
        offset: Vertical offset from bar top
    """
    for bar in bars:
        height = bar.get_height()
        ax.annotate(fmt.format(height),
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, offset),
                    textcoords="offset points",
                    ha='center', va='bottom',
                    fontsize=fontsize, fontweight='bold')


def save_figure(fig, path, dpi=300):
    """Save figure as high-DPI PNG.

    Args:
        fig: Matplotlib figure
        path: Path (with or without .png extension)
        dpi: Resolution (default 300)
    """
    from pathlib import Path
    path = Path(path)
    if path.suffix != '.png':
        path = path.with_suffix('.png')
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
