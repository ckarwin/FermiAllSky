"""Compare isotropic spectra from several all-sky fits with the official one.

Run: python compare_isotropic.py

Top: E^2 I(E) of each fit (isotropic_spectrum.txt) and of the official
Fermi-LAT isotropic template with its 1-sigma band. Bottom: the fractional
counts residual (data - model) / model of each all-sky fit, summed over its
retained sky, as in that fit's counts_spectrum.png (read from diagnostics.npz
in the same results folder).

Edit FITS below: one entry per all-sky fit, with the results folder
(Isotropic/<outdir>/bestfit_<hash>), a legend label and a colour.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plot_style import GRID, LABEL_SIZE, LEGEND_SIZE, TICK_SIZE, DPI

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

ISOTROPIC = "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/isotropic_fit"

# (results folder, label, colour)
FITS = [
    (f"{ISOTROPIC}/output_with_solar_and_loopI/bestfit_2eb3f75ed5", "Sun and Moon + Loop I", "magenta"),
    (f"{ISOTROPIC}/output_with_solar/bestfit_34e87a564c", "Sun and Moon", "cyan"),
    (f"{ISOTROPIC}/output_with_solar_and_loopI_separate/bestfit_33aa0a726a", "Sun and Moon + Loop I (separate)", "orange"),
]  # fmt: skip

# Official isotropic spectrum (energy [MeV], intensity, 1-sigma error).
OFFICIAL = "/project/ckarwin/astro/chris/my_envs/joint_analysis_main/share/fermitools/refdata/fermi/galdiffuse/iso_P8R3_SOURCE_V3_v1.txt"

# Output figure (also shown on screen); None to only show it.
OUTPUT = "compare_isotropic.png"

# Y range of the residual panel, e.g. (-0.05, 0.05); None = automatic.
RESIDUAL_YLIM = None


# ---------------------------------------------------------------------------


def fractional_residuals(results_dir):
    """Energy [MeV], (D - M)/M and its error from a fit's diagnostics.npz.

    Sums data and model over the retained (unmasked) sky, exactly as the
    all-sky fit's counts_spectrum.png does.
    """
    with np.load(Path(results_dir) / "diagnostics.npz", allow_pickle=False) as z:
        keep = z["keep"]
        edges = z["edges"]
        data = z["counts"][:, keep].sum(axis=1)
        model = z["model"][:, keep].sum(axis=1)
    energy = np.sqrt(edges[:-1] * edges[1:])
    frac = np.divide(data - model, model, out=np.full_like(model, np.nan, dtype=float),
                     where=model > 0)
    err = np.divide(np.sqrt(data), model, out=np.zeros_like(model, dtype=float), where=model > 0)
    return energy, frac, err


def main():
    fig, (ax, ar) = plt.subplots(2, 1, figsize=(8, 8), sharex=True,
                                 gridspec_kw={"height_ratios": [3, 1.3], "hspace": 0.05})

    emin, emax = np.inf, 0.0
    for results_dir, label, colour in FITS:
        energy, intensity = np.loadtxt(Path(results_dir) / "isotropic_spectrum.txt", unpack=True)
        ax.loglog(energy, energy**2 * intensity, "o-", color=colour, label=label)
        emin, emax = min(emin, energy.min()), max(emax, energy.max())
        e, frac, err = fractional_residuals(results_dir)
        ar.errorbar(e, frac, yerr=err, fmt="o", ms=4, color=colour, capsize=2)

    e_ref, i_ref, di_ref = np.loadtxt(OFFICIAL, unpack=True)
    ax.loglog(e_ref, e_ref**2 * i_ref, color="black", label="Official isotropic")
    ax.fill_between(e_ref, e_ref**2 * np.maximum(i_ref - di_ref, np.finfo(float).tiny),
                    e_ref**2 * (i_ref + di_ref), color="black", alpha=0.15,
                    label=r"Official $1\sigma$")

    ax.set_xlim(emin / 1.15, emax * 1.15)
    ax.set_ylabel(r"$E^2 I(E)$ [MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]", fontsize=LABEL_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_SIZE)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=LEGEND_SIZE, frameon=False)

    ar.axhline(0, color="grey", lw=1)
    ar.set_xscale("log")
    if RESIDUAL_YLIM is not None:
        ar.set_ylim(*RESIDUAL_YLIM)
    ar.set_xlabel("Energy [MeV]", fontsize=LABEL_SIZE)
    ar.set_ylabel("(D - M)/M", fontsize=LABEL_SIZE)
    ar.tick_params(axis="both", which="major", labelsize=TICK_SIZE)
    ar.grid(**GRID)

    if OUTPUT:
        fig.savefig(OUTPUT, dpi=DPI, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
