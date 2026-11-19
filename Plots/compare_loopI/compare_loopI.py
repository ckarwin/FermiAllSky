"""Compare the fitted Loop I spectrum with the Fermi-LAT points in Ackermann et al. (2025), Fig. 3.

Run: python compare_loopI.py
Requires: numpy, matplotlib, pyyaml (and LoopI/shell_geometry.py from this repository).
Edit the settings below to compare a different fit. See README.md.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
REFERENCE_DIR = SCRIPT_DIR / "reference"
LOOPI_DIR = SCRIPT_DIR.parents[1] / "LoopI"
sys.path.insert(0, str(LOOPI_DIR))
from shell_geometry import SHELLS, path_length  # noqa: E402

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

# Loop I S1 and S2 flux files written by Plots/plot_flux.py for your fit
# (E^2 dF/dE summed over the analysis region). Relative to this folder.
S1_FLUX_FILE = "loopI_SI_flux.dat"
S2_FLUX_FILE = "loopI_S2_flux.dat"
FIT_LABEL = "This fit (S1 + S2)"

# Analysis region of the fit (Galactic, deg): |l| <= ROI_HALF_WIDTH_DEG and
# |b| <= ROI_HALF_WIDTH_DEG around the Galactic centre.
ROI_HALF_WIDTH_DEG = 20.0

# The published intensity is the Loop I intensity at |b| > 10 deg, averaged
# over the area of Loop I, which the paper gives as about 3 sr. The fitted
# model is averaged over the brightest PAPER_SOLID_ANGLE_SR of its own sky at
# |b| > LATITUDE_CUT_DEG; SOLID_ANGLE_RANGE_SR gives the band in the plot.
LATITUDE_CUT_DEG = 10.0
PAPER_SOLID_ANGLE_SR = 3.0
SOLID_ANGLE_RANGE_SR = (2.5, 3.5)

# Energy range used for the mean fit/paper ratio [MeV].
RATIO_BAND_MEV = (300.0, 30000.0)

# Template grid (same as LoopI/build_templates.py): full-sky CAR, 0.5 deg pixels.
PIXEL_DEG = 0.5


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def load_paper_points():
    """Published Loop I points: energy [MeV], E^2 I and its 1-sigma range [MeV cm^-2 s^-1 sr^-1]."""
    data = yaml.safe_load((REFERENCE_DIR / "fig3_loopI_points.yaml").read_text())
    return (np.array(data["energy"]), np.array(data["e2_intensity"]),
            np.array(data["e2_intensity_low"]), np.array(data["e2_intensity_high"]))


def load_fit_spectrum(filename):
    """Read a plot_flux.py file: energy [MeV] and E^2 dF/dE [MeV cm^-2 s^-1] in the region."""
    data = np.loadtxt(SCRIPT_DIR / filename, comments="#", skiprows=3)
    return data[:, 0], data[:, 1]


def shell_maps():
    """Line-of-sight path length through S1 and S2 on a full-sky grid.

    Returns (maps[2, nb, nl], l, b, pixel solid angle). The absolute scale is
    irrelevant here: only ratios of integrals of each map are used.
    """
    nb, nl = int(round(180 / PIXEL_DEG)), int(round(360 / PIXEL_DEG))
    b = (np.arange(nb) - (nb - 1) / 2) * PIXEL_DEG
    l = -(np.arange(nl) - (nl - 1) / 2) * PIXEL_DEG
    b, l = np.meshgrid(b, l, indexing="ij")
    l = (l + 180) % 360 - 180
    half = np.deg2rad(PIXEL_DEG / 2)
    omega = np.deg2rad(PIXEL_DEG) * (np.sin(np.deg2rad(b) + half) - np.sin(np.deg2rad(b) - half))
    maps = np.array([path_length(l, b, shell) for shell in SHELLS])
    return maps, l, b, omega


# ---------------------------------------------------------------------------
# Conversion of the fit to the paper's quantity
# ---------------------------------------------------------------------------


def paper_footprint(weights, maps, b, omega, solid_angle):
    """Pixels forming the brightest `solid_angle` sr of the fitted Loop I map at |b| > cut.

    The fitted map is weights[0] * S1 + weights[1] * S2 (the shells share one
    spectrum, so its shape does not depend on energy).
    """
    model = np.tensordot(weights, maps, axes=1)
    candidates = np.flatnonzero(np.abs(b) > LATITUDE_CUT_DEG)
    order = candidates[np.argsort(model.ravel()[candidates])[::-1]]
    keep = order[np.cumsum(omega.ravel()[order]) <= solid_angle]
    footprint = np.zeros(model.size, bool)
    footprint[keep] = True
    return footprint.reshape(model.shape)


def fit_to_paper(e2f_roi, maps, roi, footprint, omega):
    """Convert region fluxes of S1 and S2 to the mean intensity over the footprint.

    For each shell i, the fitted intensity map is k_i(E) * map_i, with
    k_i = (region flux of shell i) / (integral of map_i over the region). The
    result is the mean of sum_i k_i * map_i over the footprint, per shell and in
    total [MeV cm^-2 s^-1 sr^-1].
    """
    area = omega[footprint].sum()
    parts = []
    for flux, shell in zip(e2f_roi, maps):
        k = flux / (shell * omega)[roi].sum()
        parts.append(k * (shell * omega)[footprint].sum() / area)
    return np.array(parts)


def interp_loglog(x, xp, fp):
    """Log-log interpolation; NaN outside xp."""
    out = np.exp(np.interp(np.log(x), np.log(xp), np.log(fp)))
    return np.where((x >= xp[0]) & (x <= xp[-1]), out, np.nan)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    paper_e, paper_y, paper_lo, paper_hi = load_paper_points()
    e1, f1 = load_fit_spectrum(S1_FLUX_FILE)
    e2, f2 = load_fit_spectrum(S2_FLUX_FILE)
    if not np.allclose(e1, e2):
        raise ValueError("S1 and S2 flux files must use the same energies")

    maps, l, b, omega = shell_maps()
    roi = (np.abs(l) <= ROI_HALF_WIDTH_DEG) & (np.abs(b) <= ROI_HALF_WIDTH_DEG)

    # Shell weights of the fitted map, from the region fluxes (any energy works,
    # since both shells share one spectrum; use the band-summed fluxes).
    weights = np.array([f.sum() / (m * omega)[roi].sum() for f, m in zip((f1, f2), maps)])

    curves = {}
    for omega_sr in (PAPER_SOLID_ANGLE_SR, *SOLID_ANGLE_RANGE_SR):
        footprint = paper_footprint(weights, maps, b, omega, omega_sr)
        curves[omega_sr] = fit_to_paper(np.array([f1, f2]), maps, roi, footprint, omega)
    s1, s2 = curves[PAPER_SOLID_ANGLE_SR]
    total = s1 + s2
    band = np.sort([curves[w].sum(axis=0) for w in SOLID_ANGLE_RANGE_SR], axis=0)

    # Fit / paper ratio at the paper energies inside the fit range.
    fit_at_paper = interp_loglog(paper_e, e1, total)
    ratio = fit_at_paper / paper_y
    in_band = (paper_e >= RATIO_BAND_MEV[0]) & (paper_e <= RATIO_BAND_MEV[1]) & np.isfinite(ratio)
    mean_ratio = np.exp(np.mean(np.log(ratio[in_band])))

    # Alternative definition, for reference: total model flux at |b| > cut / paper area.
    high = np.abs(b) > LATITUDE_CUT_DEG
    flux_ratio = sum(w * (m * omega)[high].sum() for w, m in zip(weights, maps)) / sum(
        w * (m * omega)[roi].sum() for w, m in zip(weights, maps))
    total_over_area = (f1 + f2) * flux_ratio / PAPER_SOLID_ANGLE_SR
    mean_ratio_total = np.exp(np.mean(np.log(
        (interp_loglog(paper_e, e1, total_over_area) / paper_y)[in_band])))

    # --- Tables -----------------------------------------------------------
    header = "energy_MeV\tS1\tS2\ttotal\ttotal_low\ttotal_high  [E^2 I, MeV cm^-2 s^-1 sr^-1]"
    np.savetxt(SCRIPT_DIR / "converted_fit.csv", np.c_[e1, s1, s2, total, band[0], band[1]],
               header=header, delimiter="\t")
    footprint = paper_footprint(weights, maps, b, omega, PAPER_SOLID_ANGLE_SR)
    summary = {
        "paper_solid_angle_sr": PAPER_SOLID_ANGLE_SR,
        "footprint_north_fraction": float(omega[footprint & (b > 0)].sum() / omega[footprint].sum()),
        "S1_share_of_footprint_intensity": float(np.mean(s1 / total)),
        "mean_fit_over_paper": float(mean_ratio),
        "mean_fit_over_paper_total_flux_definition": float(mean_ratio_total),
        "ratio_band_MeV": list(RATIO_BAND_MEV),
    }
    with open(SCRIPT_DIR / "summary.csv", "w") as handle:
        for key, value in summary.items():
            handle.write(f"{key},{value}\n")

    # --- Plot ---------------------------------------------------------------
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(7, 7), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1]})
    ax.errorbar(paper_e, paper_y, yerr=[paper_y - paper_lo, paper_hi - paper_y], fmt="o",
                color="red", ms=4, capsize=2,
                label="Fermi-LAT Loop I, |b| > 10° (Ackermann et al. 2025, Fig. 3)")
    ax.fill_between(e1, band[0], band[1], color="C0", alpha=0.25, lw=0,
                    label=f"{FIT_LABEL}, area {SOLID_ANGLE_RANGE_SR[0]}–{SOLID_ANGLE_RANGE_SR[1]} sr")
    ax.plot(e1, total, color="C0", lw=2, label=f"{FIT_LABEL}, area {PAPER_SOLID_ANGLE_SR} sr")
    ax.plot(e1, s1, color="C0", ls=":", lw=1.2, label="S1")
    ax.plot(e1, s2, color="C0", ls="--", lw=1.2, label="S2")
    curve = np.loadtxt(LOOPI_DIR / "leptonic_curve.txt")
    ax.plot(curve[:, 0], curve[:, 2], color="purple", lw=1, alpha=0.6,
            label="Leptonic model (paper; spectral shape of the fit)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(80, 1.2e5)
    ax.set_ylim(1e-6, 3e-3)
    ax.set_ylabel(r"$E^2\,I(E)$ [MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]")
    ax.set_title(f"Loop I at |b| > {LATITUDE_CUT_DEG:.0f}°: mean intensity over "
                 f"{PAPER_SOLID_ANGLE_SR} sr", fontsize=11)
    ax.legend(fontsize=7.5, loc="lower left")
    axr.errorbar(paper_e, ratio, yerr=[ratio - fit_at_paper / paper_hi, fit_at_paper / paper_lo - ratio],
                 fmt="o", color="k", ms=3, capsize=2)
    axr.axhline(1, color="gray", lw=0.8)
    axr.axhline(mean_ratio, color="C0", ls="--", lw=1,
                label=f"mean {mean_ratio:.2f} ({RATIO_BAND_MEV[0] / 1e3:g}–{RATIO_BAND_MEV[1] / 1e3:g} GeV)")
    axr.set_yscale("log")
    axr.set_ylabel("fit / paper")
    axr.set_xlabel("Energy [MeV]")
    axr.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(SCRIPT_DIR / "comparison.png", dpi=150)
    fig.savefig(SCRIPT_DIR / "comparison.pdf")
    print(f"Mean fit/paper ratio ({RATIO_BAND_MEV[0]:g}-{RATIO_BAND_MEV[1]:g} MeV): {mean_ratio:.3f}")
    print(f"  with the total-flux definition: {mean_ratio_total:.3f}")
    print("Wrote comparison.png, comparison.pdf, converted_fit.csv, summary.csv")


if __name__ == "__main__":
    main()
