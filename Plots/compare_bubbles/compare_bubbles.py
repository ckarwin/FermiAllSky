"""Compare the fitted Fermi-bubble spectra with Ackermann et al. (2017), Fig. 9a.

Run: python compare_bubbles.py
Requires: numpy, pandas, matplotlib, astropy, healpy, pyyaml.
Edit the settings below to compare a different fit. See README.md.
"""
from pathlib import Path

import healpy as hp
import matplotlib
import numpy as np
import pandas as pd
import yaml
from astropy.io import fits

#matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ----- Settings: change these filenames to use another fit -----
# Relative paths are resolved beside this script; absolute paths also work.
LOW_LAT_FLUX_FILE = "../baseline_with_all_components_pl_bubbles/bubbles_low_flux.dat"
HIGH_LAT_FLUX_FILE = "../baseline_with_all_components_pl_bubbles/bubbles_high_flux.dat"
FIT_LABEL = "my fit"

# Keep this geometry consistent with the ROI used to export your spectra.
ROI_LON_LIMIT_DEG = 20
ROI_LAT_LIMIT_DEG = 20
LATITUDE_SPLIT_DEG = 10
BAND_GEV = (1, 80)

# The published bubble template is zero inside the masked cores of the 200
# brightest 3FGL sources, and the published spectra exclude those holes. The
# fitted spectra cover the whole ROI, so fill the holes (with the mean of the
# unmasked template within this radius) when converting the paper to the ROI.
FILL_SOURCE_HOLES = True
HOLE_FILL_RADIUS_DEG = 2.0

SCRIPT_DIR = Path(__file__).resolve().parent
REFERENCE_DIR = SCRIPT_DIR / "reference"


def load_fit_spectrum(filename):
    """Read the exporter table; convert energy and E²dF/dE from MeV to GeV."""
    path = SCRIPT_DIR / filename
    table = pd.read_csv(path, skiprows=2, sep=r"\s+")
    energy = table["Energy[MeV]"].to_numpy() / 1000
    sed = table["energy_flux[MeV/cm^2/s]"].to_numpy() / 1000
    order = np.argsort(energy)
    energy, sed = energy[order], sed[order]
    if (not np.all(np.isfinite(energy)) or not np.all(np.isfinite(sed))
            or np.any(energy <= 0) or np.any(sed <= 0)
            or np.any(np.diff(energy) <= 0) or len(energy) < 2):
        raise ValueError(f"{path}: need positive, finite fluxes at distinct energies.")
    # calculate_flux() already sums over the ROI: do NOT apply template factors.
    return energy, sed


def load_reference():
    """Load the published spectra and maps used for the aperture conversion."""
    with (REFERENCE_DIR / "fig9a_SCA_bubbles_model_spectra.yaml").open() as stream:
        spectra = yaml.safe_load(stream)
    template = np.load(REFERENCE_DIR / "fig8_bubbles_template_2compSCA.npy")
    source_mask = np.load(REFERENCE_DIR / "fig2_mask_200_3FGL_sources.npy")
    nside = hp.npix2nside(len(template))
    lon, lat = hp.pix2ang(nside, np.arange(len(template)), lonlat=True)
    lon = (lon + 180) % 360 - 180
    roi = (abs(lon) < ROI_LON_LIMIT_DEG) & (abs(lat) < ROI_LAT_LIMIT_DEG)
    with fits.open(REFERENCE_DIR / "exposure.fits") as hdus:
        exposure = hdus[1].data.field(0).copy()
        exposure_energy = hdus[2].data.field(0).copy() / 1000
    filled = fill_source_holes(template, source_mask) if FILL_SOURCE_HOLES else template
    return spectra, template, filled, source_mask, lat, roi, exposure, exposure_energy


def fill_source_holes(template, source_mask):
    """Replace the template in masked source holes by the local unmasked mean.

    Each masked pixel gets the mean template value of the unmasked pixels
    within HOLE_FILL_RADIUS_DEG. Pixels without unmasked neighbours keep 0.
    """
    nside = hp.npix2nside(len(template))
    filled = template.copy()
    radius = np.radians(HOLE_FILL_RADIUS_DEG)
    for pixel in np.where(source_mask == 0)[0]:
        neighbours = hp.query_disc(nside, hp.pix2vec(nside, pixel), radius)
        neighbours = neighbours[source_mask[neighbours] > 0]
        if len(neighbours):
            filled[pixel] = template[neighbours].mean()
    return filled


def paper_roi_factors(energies, template, filled, source_mask, latitude_mask,
                      roi, exposure, exposure_energy):
    """Convert the paper's masked all-sky flux to the full rectangular ROI.

    Approximate intensity = released counts morphology / exposure.
    Factor = sum(hole-filled intensity in ROI) / sum(intensity over unmasked sky).
    The denominator matches the published spectra (whole sky, source cores
    excluded); the numerator is the whole ROI, including the filled holes.
    Both sums select the same latitude component. This is a model-based
    aperture conversion, not a PSF deconvolution.
    """
    factors = []
    for energy in energies:
        index = np.clip(np.searchsorted(exposure_energy, energy) - 1,
                        0, len(exposure_energy) - 2)
        weight = (np.log(energy) - np.log(exposure_energy[index])) / np.log(
            exposure_energy[index + 1] / exposure_energy[index])
        # Interpolate exposure linearly against log(energy).
        at_energy = (exposure[:, index] * (1 - weight)
                     + exposure[:, index + 1] * weight)
        intensity = np.maximum(template, 0) / at_energy * latitude_mask
        filled_intensity = np.maximum(filled, 0) / at_energy * latitude_mask
        factors.append(filled_intensity[roi].sum() / (intensity * source_mask).sum())
    return np.array(factors)


def interpolate_sed(energy, sample_energy, sample_sed):
    """Log-log interpolation; return NaN outside the sampled energy range."""
    return np.exp(np.interp(np.log(energy), np.log(sample_energy),
                            np.log(sample_sed), left=np.nan, right=np.nan))


def integrated_energy_flux(energy, sed):
    """Integrate E²dF/dE over dln(E), giving GeV cm^-2 s^-1.

    Return NaN if the input does not cover the whole requested energy band.
    The exporter's 'Total energy flux' header is not used.
    """
    grid = np.geomspace(*BAND_GEV, 2000)
    return np.trapezoid(interpolate_sed(grid, energy, sed), x=np.log(grid))


def main():
    spectra, template, filled, source_mask, lat, roi, exposure, exposure_energy = load_reference()
    paper_energy = np.array(spectra["Energies"])
    components = [
        ("low", "Low lat bubbles", abs(lat) < LATITUDE_SPLIT_DEG,
         LOW_LAT_FLUX_FILE, "#a32ab2"),
        ("high", "Bubbles", abs(lat) > LATITUDE_SPLIT_DEG,
         HIGH_LAT_FLUX_FILE, "#c23a3a"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True,
                             layout="constrained")
    reference_rows = []
    summary_rows = []

    for ax, (tag, key, latitude_mask, filename, color) in zip(axes, components):
        # 1. Convert only the published spectrum to our ROI.
        factors = paper_roi_factors(paper_energy, template, filled, source_mask,
                                    latitude_mask, roi, exposure, exposure_energy)
        paper_sed = np.array(spectra[key]) * factors
        paper_error = np.array(spectra[key + "_err"]) * factors

        # 2. Load the fit, which is already integrated over the ROI.
        fit_energy, fit_sed = load_fit_spectrum(filename)

        # 3. Plot both spectra in the same units and aperture.
        ax.errorbar(paper_energy, paper_sed, yerr=paper_error, fmt="o", ms=3,
                    color="black", capsize=2, label="Ackermann 2017: converted to ROI")
        ax.plot(fit_energy, fit_sed, color=color, lw=2,
                label=f"{FIT_LABEL}: ROI counts/exposure flux")
        ax.set(title=tag.capitalize() + "-latitude component", xscale="log",
               yscale="log", xlim=(.3, 100), ylim=(3e-9, 8e-7), xlabel="Energy [GeV]")
        ax.grid(alpha=.2)
        ax.legend(fontsize=8)

        # 4. Save the 10 GeV comparison and the energy-integrated band flux.
        band_name = f"{BAND_GEV[0]:g}_{BAND_GEV[1]:g}GeV"
        summary_rows.append({
            "component": tag,
            "paper_roi_factor_10GeV": float(np.interp(10, paper_energy, factors)),
            "paper_E2F_10GeV": interpolate_sed(10, paper_energy, paper_sed),
            "user_E2F_10GeV": interpolate_sed(10, fit_energy, fit_sed),
            f"paper_energy_flux_{band_name}": integrated_energy_flux(paper_energy, paper_sed),
            f"user_energy_flux_{band_name}": integrated_energy_flux(fit_energy, fit_sed),
        })
        for energy, sed, factor in zip(paper_energy, paper_sed, factors):
            reference_rows.append({"component": tag, "energy_GeV": energy,
                                   "paper_roi_E2F": sed, "aperture_factor": factor})

    axes[0].set_ylabel(r"$E^2 dF/dE$ [GeV cm$^{-2}$ s$^{-1}$]")
    fig.suptitle(
        rf"Integrated over $|l|<{ROI_LON_LIMIT_DEG}^\circ,\ |b|<{ROI_LAT_LIMIT_DEG}^\circ$"
        " (no point-source holes)")
    fig.savefig(SCRIPT_DIR / "comparison.png", dpi=180)
    fig.savefig(SCRIPT_DIR / "comparison.pdf")
    plt.show()
    plt.close(fig)
    pd.DataFrame(reference_rows).to_csv(SCRIPT_DIR / "converted_reference.csv", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(SCRIPT_DIR / "summary.csv", index=False)
    print(summary.to_string(index=False))
    if summary.isna().any().any():
        print("NaN/blank summary values mean the spectrum does not cover the requested energy range.")


if __name__ == "__main__":
    main()
