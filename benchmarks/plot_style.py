"""
Shared plotting style configuration for all benchmarks.

Provides consistent, publication-quality plots with Helvetica font.
"""

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


def _register_system_fonts():
    """Register system fonts so matplotlib can find Arial, Liberation Sans, etc."""
    font_dirs = [
        "/usr/share/fonts/truetype/msttcorefonts",
        "/usr/share/fonts/truetype/liberation",
        "/usr/share/fonts/truetype",
    ]
    for d in font_dirs:
        font_files = fm.findSystemFonts(fontpaths=[d])
        for f in font_files:
            try:
                fm.fontManager.addfont(f)
            except Exception:
                pass


def setup_plot_style():
    """Configure matplotlib for publication-quality plots with consistent styling.

    Call this at the start of any visualization function.
    """
    # Reset to default first
    plt.rcdefaults()

    # Register system fonts (conda env may not see them by default)
    _register_system_fonts()

    # Font configuration - Helvetica-like fonts with fallbacks
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial",
                "Nimbus Sans",
                "Liberation Sans",
                "DejaVu Sans",
            ],
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 11,
            "figure.titlesize": 16,
        }
    )

    # Axes styling
    plt.rcParams.update(
        {
            "axes.linewidth": 1.2,
            "axes.edgecolor": "#333333",
            "axes.labelcolor": "#333333",
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )

    # Grid styling
    plt.rcParams.update(
        {
            "axes.grid": False,
            "grid.alpha": 0.3,
            "grid.linewidth": 0.5,
        }
    )

    # Tick styling
    plt.rcParams.update(
        {
            "xtick.color": "#333333",
            "ytick.color": "#333333",
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 1.0,
            "ytick.major.width": 1.0,
            "xtick.major.size": 4,
            "ytick.major.size": 4,
        }
    )

    # Figure styling
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.edgecolor": "white",
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.1,
        }
    )

    # Legend styling
    plt.rcParams.update(
        {
            "legend.frameon": True,
            "legend.framealpha": 0.9,
            "legend.edgecolor": "#cccccc",
            "legend.fancybox": False,
        }
    )


# Standard color palettes for benchmarks
COLORS = {
    # Segmentation methods
    "cellpose": "#1f77b4",  # Blue
    "cellpose4": "#ff7f0e",  # Orange
    "stardist": "#2ca02c",  # Green
    # Spot calling methods
    "standard": "#E41A1C",  # Red
    "spotiflow": "#FFD700",  # Yellow
    # Feature extraction methods
    "cp_measure": "#8c564b",  # Brown
    "cp_multichannel": "#e377c2",  # Pink
    # Merge methods
    "fast": "#1f77b4",  # Blue
    "stitch": "#ff7f0e",  # Orange
    # Clustering comparison (colorblind-friendly)
    "brieflow": "#0077BB",  # Teal blue
    "funk": "#EE7733",  # Coral orange
}

# Standard figure sizes
FIGSIZE = {
    "single": (6, 5),
    "wide": (10, 5),
    "square": (8, 8),
    "double": (12, 5),
    "quad": (10, 8),
}


def get_method_colors(methods, palette="default"):
    """Get consistent colors for a list of methods.

    Args:
        methods: List of method names
        palette: Color palette to use ('default' uses COLORS dict)

    Returns:
        List of colors matching the methods
    """
    if palette == "default":
        return [COLORS.get(m.lower(), "#7f7f7f") for m in methods]
    else:
        import seaborn as sns

        return sns.color_palette(palette, len(methods))


def _significance_stars(p):
    """Convert p-value to star notation."""
    if p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return "ns"


def _add_significance_brackets(ax, data, x, y, groups, ylim=None):
    """Add pairwise Mann-Whitney U significance brackets between all group pairs."""
    from itertools import combinations
    from scipy.stats import mannwhitneyu

    pairs = list(combinations(range(len(groups)), 2))
    if not pairs:
        return

    # Get current y-axis range to place brackets above data
    y_max = data[y].max()
    y_min = data[y].min()
    y_range = y_max - y_min
    bracket_height = y_range * 0.02
    bracket_start = y_max + y_range * 0.05

    # If ylim is set, compress brackets to fit within the specified range
    if ylim is not None:
        ceiling = ylim[1]
        available = ceiling - y_max
        n_brackets = len(pairs)
        per_bracket = available / (n_brackets + 1.5)
        bracket_height = per_bracket * 0.3
        bracket_start = y_max + per_bracket * 0.5
        bracket_step = per_bracket
    else:
        bracket_step = y_range * 0.06

    for k, (i, j) in enumerate(pairs):
        vals_a = data[data[x] == groups[i]][y].dropna()
        vals_b = data[data[x] == groups[j]][y].dropna()
        if len(vals_a) < 2 or len(vals_b) < 2:
            continue

        _, p = mannwhitneyu(vals_a, vals_b, alternative="two-sided")
        stars = _significance_stars(p)

        y_pos = bracket_start + k * bracket_step
        # Draw bracket
        ax.plot([i, i, j, j], [y_pos, y_pos + bracket_height,
                                y_pos + bracket_height, y_pos],
                color="black", linewidth=1, clip_on=False)
        # Draw stars
        ax.text((i + j) / 2, y_pos + bracket_height, stars,
                ha="center", va="bottom", fontsize=8, color="black")

    # Adjust ylim to fit brackets (only when no explicit ylim)
    if ylim is None:
        new_top = bracket_start + len(pairs) * bracket_step + y_range * 0.05
        current_bottom = ax.get_ylim()[0]
        ax.set_ylim(current_bottom, new_top)


def box_strip(ax, data, x, y, palette, order=None, ylabel=None, title=None,
              fmt="auto", ylim=None):
    """Draw box-and-whisker + jittered strip plot on a single axes.

    Better than bar charts for small sample sizes (n~8) because it shows
    the full distribution: median, IQR (box), range (whiskers), and all
    individual data points. Pairwise Mann-Whitney U significance brackets
    are added automatically.

    Args:
        ax: Matplotlib axes
        data: DataFrame with the data
        x: Column name for categories (x-axis)
        y: Column name for values (y-axis)
        palette: Dict mapping category names to colors
        order: List of category names in display order
        ylabel: Y-axis label
        title: Axes title
        fmt: Format for median annotation — "auto", "pct", "int", or a format spec
        ylim: Tuple for y-axis limits
    """
    import seaborn as sns

    sns.boxplot(
        data=data, x=x, y=y, hue=x, order=order, palette=palette,
        linewidth=1.2, fliersize=0, ax=ax, legend=False,
        boxprops=dict(alpha=0.6),
        medianprops=dict(color="black", linewidth=1.5),
    )
    sns.stripplot(
        data=data, x=x, y=y, order=order, color="black",
        size=5, jitter=0.12, ax=ax, alpha=0.7, zorder=10, legend=False,
    )
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, fontweight="bold")
    ax.set_xlabel("")
    if ylim:
        ax.set_ylim(ylim)

    # Collect medians (returned so callers can print for publication)
    groups = order if order else sorted(data[x].unique())
    medians = {}
    for g in groups:
        vals = data[data[x] == g][y].dropna()
        if len(vals) > 0:
            medians[g] = vals.median()

    # Add significance brackets (pairwise Mann-Whitney U)
    _add_significance_brackets(ax, data, x, y, groups, ylim=ylim)
    # Re-apply explicit ylim if provided (overrides bracket adjustment)
    if ylim:
        ax.set_ylim(ylim)

    return medians


def print_summary_table(title, data, x, metrics, fmt_map=None):
    """Print a publication-ready summary table with median ± IQR for each metric.

    Args:
        title: Table title
        data: DataFrame with the data
        x: Column name for groups (rows of the table)
        metrics: List of (column_name, display_name) tuples
        fmt_map: Optional dict mapping column_name to format function.
                 Default formats: ".1f" for floats, ",.0f" for large numbers.
    """
    groups = sorted(data[x].unique())
    # Build header
    col_widths = [max(18, max(len(str(g)) for g in groups) + 2)]
    for _, display in metrics:
        col_widths.append(max(20, len(display) + 2))

    header = f"{'Method':<{col_widths[0]}}"
    for i, (_, display) in enumerate(metrics):
        header += f"  {display:>{col_widths[i+1]}}"

    print(f"\n{'=' * len(header)}")
    print(f"  {title}")
    print(f"{'=' * len(header)}")
    print(header)
    print("-" * len(header))

    for g in groups:
        vals_by_col = {}
        for col, _ in metrics:
            v = data[data[x] == g][col].dropna()
            vals_by_col[col] = v

        row = f"{str(g):<{col_widths[0]}}"
        for i, (col, _) in enumerate(metrics):
            v = vals_by_col[col]
            med = v.median()
            iqr_lo = v.quantile(0.25)
            iqr_hi = v.quantile(0.75)

            if fmt_map and col in fmt_map:
                cell = fmt_map[col](med, iqr_lo, iqr_hi)
            elif abs(med) >= 1000:
                cell = f"{med:,.0f} [{iqr_lo:,.0f}–{iqr_hi:,.0f}]"
            elif abs(med) >= 1:
                cell = f"{med:.1f} [{iqr_lo:.1f}–{iqr_hi:.1f}]"
            else:
                cell = f"{med:.2f} [{iqr_lo:.2f}–{iqr_hi:.2f}]"
            row += f"  {cell:>{col_widths[i+1]}}"
        print(row)

    print(f"{'=' * len(header)}")
    print("  Values: median [IQR]\n")


def add_value_labels(ax, bars, fmt="{:.0f}", fontsize=10, offset=5):
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
        ax.annotate(
            fmt.format(height),
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, offset),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=fontsize,
            fontweight="bold",
        )


def save_figure(fig, path, dpi=300, transparent=False):
    """Save figure, respecting the file extension (png, pdf, svg, etc.).

    Args:
        fig: Matplotlib figure
        path: Path with desired extension
        dpi: Resolution for raster formats (default 300)
        transparent: If True, save with transparent background
    """
    from pathlib import Path

    path = Path(path)
    if not path.suffix:
        path = path.with_suffix(".png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight", transparent=transparent)
