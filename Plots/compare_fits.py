"""Overlay the results of several inner-Galaxy fits.

Run: python compare_fits.py

Each fit is a plot folder made by make_all_plots.py (Plots/<plots_subdir>,
written at the end of run_3ML_master.py). Edit RUNS below. Makes, in
OUTPUT_DIR:

* gce_overlay.png/.pdf: the GCE spectrum of each fit (best fit and 68% band)
  over the Ajello et al. (2016) Figure 13 envelopes.
* bubbles_overlay.png/.pdf: the bubble spectra of each fit over Ackermann et
  al. (2017), converted to the fit ROI.
* loopI_overlay.png/.pdf: the Loop I spectrum of the all-sky fit used by each
  fit, over the published Fermi-LAT Loop I spectra.
* residuals_overlay.png/.pdf: the fractional counts residual (D - M)/M of
  each fit.

A plot is skipped (with a message) when no fit has the files it needs.
"""

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

PLOTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PLOTS_DIR))
sys.path.insert(0, str(PLOTS_DIR / "compare_bubbles"))
sys.path.insert(0, str(PLOTS_DIR / "compare_loopI"))
sys.path.insert(0, str(PLOTS_DIR.parent / "Main"))
from plot_style import GRID, LABEL_SIZE, LEGEND_SIZE, RUN_COLOURS, TICK_SIZE  # noqa: E402

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

# (plot folder, label): folders are taken from Plots/ unless absolute.
RUNS = [
    ("baseline_with_all_components_pl_bubbles", "Baseline"),
]

# Output folder (taken from Plots/ unless absolute).
OUTPUT_DIR = "compare_fits_output"

# Y range of the residual overlay.
RESIDUAL_YLIM = (-0.1, 0.1)


# ---------------------------------------------------------------------------


def _folder(path):
    path = Path(path)
    return path if path.is_absolute() else PLOTS_DIR / path


def _save(fig, outdir, name):
    fig.savefig(outdir / f"{name}.png", dpi=150, bbox_inches="tight")
    fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {outdir / name}.png/.pdf")


def gce_overlay(runs, outdir):
    """GCE spectra of all fits over the Figure 13 envelopes."""
    from compare_gce_ajello2016 import (
        LEGEND_SIZE as GCE_LEGEND_SIZE,
        draw_fig13_envelopes,
        style_gce_axes,
    )

    tables = []
    for (folder, label), colour in runs:
        path = folder / "gce_spectrum.csv"
        if not path.is_file():
            print(f"GCE overlay: no gce_spectrum.csv in {folder}, skipped")
            continue
        header = [line for line in path.read_text().splitlines() if line.startswith("#")]
        roi = next(line.split("=", 1)[1].strip() for line in header if "comparison_roi_deg" in line)
        scale = float(next(line.split("=", 1)[1] for line in header if "paper_scale" in line))
        tables.append((label, colour, np.loadtxt(path), roi, scale))
    if not tables:
        print("GCE overlay: no fit has gce_spectrum.csv, skipped")
        return
    if len({t[3] for t in tables}) > 1:
        raise ValueError("The fits use different GCE comparison regions; overlay is not meaningful.")
    roi = tuple(float(v) for v in tables[0][3].strip("()").split(","))

    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
    draw_fig13_envelopes(ax, tables[0][4])
    for label, colour, table, _, _ in tables:
        energy, best, low, high = table.T
        ax.fill_between(energy, low, high, color=colour, alpha=0.25, lw=0, zorder=9)
        ax.plot(energy, best, color=colour, lw=2.5, label=label, zorder=10)
    style_gce_axes(ax, roi)
    ax.legend(fontsize=GCE_LEGEND_SIZE, loc="best", frameon=False)
    _save(fig, outdir, "gce_overlay")


def bubbles_overlay(runs, outdir):
    """Bubble spectra of all fits over Ackermann et al. (2017)."""
    import compare_bubbles as cb

    present = [(f, l, c) for (f, l), c in runs
               if (f / "bubbles_low_flux.dat").is_file() and (f / "bubbles_high_flux.dat").is_file()]
    if not present:
        print("Bubbles overlay: no fit has bubbles_*_flux.dat, skipped")
        return
    paper = cb.paper_spectra()
    fig, axes = cb.comparison_figure()
    for ax, (tag, _, _, _) in zip(axes, cb.COMPONENTS):
        cb.plot_paper(ax, paper[tag])
        for folder, label, colour in present:
            energy, sed = cb.load_fit_spectrum(folder / f"bubbles_{tag}_flux.dat")
            ax.plot(energy, sed, color=colour, lw=2.5, label=label)
        ax.legend(fontsize=LEGEND_SIZE, frameon=False)
    _save(fig, outdir, "bubbles_overlay")


def loopI_overlay(runs, outdir):
    """Loop I spectra (of the all-sky fit used by each fit) over the published ones."""
    import compare_loopI as cl

    present = []
    for (folder, label), colour in runs:
        info_file = folder / "run_info.json"
        allsky = json.loads(info_file.read_text()).get("allsky_results_dir") if info_file.is_file() else None
        if allsky:
            present.append((allsky, label, colour))
        else:
            print(f"Loop I overlay: no all-sky results folder recorded in {folder}, skipped")
    if not present:
        print("Loop I overlay: no fit has an all-sky results folder, skipped")
        return
    references = cl.load_references()
    fig, ax = plt.subplots(figsize=(9, 7))
    cl.draw_references(ax, references)
    regions = []
    for allsky, label, colour in present:
        energy, curves, _ = cl.fit_curves(allsky, references)
        drawn = set()
        for c in curves:
            key = cl.region_text(c["ref"]["region"])
            if key in drawn:
                continue
            drawn.add(key)
            if key not in regions:
                regions.append(key)
            ls = cl.REGION_STYLES[regions.index(key) % len(cl.REGION_STYLES)]
            ax.plot(energy, c["total"], color=colour, ls=ls, lw=2)
        ax.plot([], [], color=colour, lw=2, label=label)
    for i, key in enumerate(regions):
        ax.plot([], [], color="0.3", ls=cl.REGION_STYLES[i % len(cl.REGION_STYLES)], lw=2,
                label=f"fit, {key}")
    cl.draw_leptonic_curve(ax)
    cl.style_spectrum_axes(ax)
    ax.set_title("Loop I: all-sky fits vs. published Fermi-LAT spectra", fontsize=13)
    ax.legend(fontsize=9, loc="lower left", frameon=False)
    _save(fig, outdir, "loopI_overlay")


def residuals_overlay(runs, outdir):
    """Fractional counts residuals (D - M)/M of all fits (make_counts_plot.py)."""
    present = [(f, l, c) for (f, l), c in runs if (f / "counts_spectrum.csv").is_file()]
    if not present:
        print("Residual overlay: no fit has counts_spectrum.csv, skipped")
        return
    fig, ax = plt.subplots(figsize=(9, 4.5))
    n = len(present)
    for i, (folder, label, colour) in enumerate(present):
        table = np.genfromtxt(folder / "counts_spectrum.csv", delimiter=",", names=True)
        energy = table["energy_MeV"]
        # Small horizontal offsets so the error bars of the fits do not overlap.
        shift = 10 ** (0.012 * (i - (n - 1) / 2))
        ax.errorbar(energy * shift, table["frac_resid"], yerr=table["frac_resid_err"], ls="",
                    marker="s", ms=5, color=colour, capsize=2, label=label)
    ax.axhline(y=0, ls="-", color="grey", lw=1)
    ax.set_xscale("log")
    ax.set_xlim(3e2, 0.1e6)
    ax.set_ylim(*RESIDUAL_YLIM)
    ax.set_xlabel("Energy [MeV]", fontsize=LABEL_SIZE)
    ax.set_ylabel("(D - M)/M", fontsize=LABEL_SIZE)
    ax.tick_params(axis="both", labelsize=TICK_SIZE)
    ax.grid(**GRID)
    ax.legend(fontsize=LEGEND_SIZE, frameon=False, ncol=2)
    _save(fig, outdir, "residuals_overlay")


def main():
    outdir = _folder(OUTPUT_DIR)
    outdir.mkdir(parents=True, exist_ok=True)
    runs = [((_folder(folder), label), RUN_COLOURS[i % len(RUN_COLOURS)])
            for i, (folder, label) in enumerate(RUNS)]
    for folder, _ in (r[0] for r in runs):
        if not folder.is_dir():
            raise FileNotFoundError(f"Plot folder not found: {folder}")
    for overlay in (gce_overlay, bubbles_overlay, loopI_overlay, residuals_overlay):
        overlay(runs, outdir)


if __name__ == "__main__":
    main()
