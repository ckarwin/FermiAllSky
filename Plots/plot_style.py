"""Shared figure style for the comparison plots.

Matches the counts and flux spectra (make_counts_plot.py, plot_flux.py):
matplotlib's default font, axis labels at 16 pt, tick labels at 14 pt and a
light dashed grid.
"""

LABEL_SIZE = 16
TICK_SIZE = 14
LEGEND_SIZE = 11
TITLE_SIZE = 14
GRID = dict(ls="--", color="grey", alpha=0.2)

# Colours for overlays of several fits, in a fixed order.
RUN_COLOURS = [
    "#1f77b4",
    "#d62728",
    "#2ca02c",
    "#ff7f0e",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#17becf",
]


def style_axes(ax, xlabel=None, ylabel=None, grid=True):
    """Apply the shared label, tick and grid style to one axes."""
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=LABEL_SIZE)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=LABEL_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_SIZE)
    if grid:
        ax.grid(which="major", **GRID)
