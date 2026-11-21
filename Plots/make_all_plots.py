"""Make every plot of one inner-Galaxy fit (Main/run_3ML_master.py).

Called at the end of run_3ML_master.py (``make_plots = True``):

    from make_all_plots import make_all_plots
    make_all_plots(plot_dir, savename, fermipy_outdir, ...)

or again later, from the files the fit left behind:

    python make_all_plots.py /path/to/Plots/<plot_dir>

(the folder holds run_info.json, written by the first call).

Everything goes to ``plot_dir`` (a folder inside Plots/):

=================================  ============================================
counts_spectra.png/.pdf            make_counts_plot.py: counts per component,
                                   (D - M)/M below
counts_spectrum.csv                the numbers of that plot
<name>_flux.dat                    plot_flux.py: ROI flux of each component
flux_spectra_<savename>.png/.pdf   plot_flux.py: E^2 dN/dE per component
count_map.png, model_map.png,      plot_spatial.py: summed counts and model
residuals.png                      maps, data - model
GC_<i>-<j>_residual.fits, ...      custom_client.py (3d residual maps):
residuals_<mode>_<i>-<j>.png       data - model, (data - model)/model and
residual_maps_<i>-<j>.png          (data - model)/sqrt(model) in energy-bin
                                   ranges i..j-1; one figure with all three
GCE_Ajello2016_comparison.pdf/png  Main/compare_gce_ajello2016.py
gce_spectrum.csv                   the fitted GCE curve of that plot
threeML_analysis_results.fits      the 3ML results (to redo the GCE plot)
bubbles_comparison.png/.pdf, ...   compare_bubbles/compare_bubbles.py
loopI_comparison.png/.pdf, ...     compare_loopI/compare_loopI.py (all-sky fit)
run_info.json                      inputs of this run (for compare_fits.py)
make_all_plots.log                 which steps worked
=================================  ============================================

Each step runs on its own: a failing step is logged and the others still
run, so a finished fit is never lost to a plotting error.

The flux tables are computed exactly as lat_binned_analysis's
``Plots.calculate_flux`` does (same pixel-by-pixel sum and exposure lookup,
same output format), but vectorized and with every cube read once.
"""

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np

PLOTS_DIR = Path(__file__).resolve().parent
REPO = PLOTS_DIR.parent
sys.path.insert(0, str(PLOTS_DIR))

# Energy-bin index ranges [start, stop) of the 3d residual maps (custom_client.py):
# 300 MeV - 1.28 GeV, 1.28 - 9.79 GeV, 9.79 - 100 GeV for 300 MeV - 100 GeV at 8 bins/decade.
RESIDUAL_BINS = [[0, 5], [5, 12], [12, 20]]

RESIDUAL_MODES = [
    ("diff", "_residual.fits", "data - model [counts]"),
    ("frac", "_fractional_residual.fits", "(data - model) / model"),
    ("sqrt", "_sqrt_residual.fits", r"(data - model) / $\sqrt{\rm model}$"),
]

# ---------------------------------------------------------------------------
# Components of the fit
# ---------------------------------------------------------------------------
# key: (file name pattern, folder, flux-table name, label in the flux plot)
#   folder "fermipy" = fermipy_outdir, "run" = output/ of the run directory.
# The file names are those used by make_counts_plot.py and plot_flux.py.
COMPONENTS = {
    "observed": ("ccube_00.fits", "fermipy", "observed", "total observed"),
    "model": ("mcube_fit_roi_{s}_00.fits", "fermipy", "model", "total model"),
    "point_sources": ("mcube_point_sources_{s}.fits", "run", "pcs", "point sources"),
    "isodiff": ("mcube_isodiff_{s}_00.fits", "fermipy", "isodiff", "isodiff"),
    "pion_decay": ("mcube_pion_decay_{s}_00.fits", "fermipy", "pion_decay", "pion_decay"),
    "ic": ("mcube_ic_{s}_00.fits", "fermipy", "ics", "ics"),
    "brem": ("mcube_brem_{s}_00.fits", "fermipy", "bremss", "bremss"),
    "GCE_NFW": ("mcube_GCE_NFW_{s}_00.fits", "fermipy", "gce", "GCE"),
    "bubbles_low": ("mcube_bubbles_low_{s}_00.fits", "fermipy", "bubbles_low", "Bubbles low"),
    "bubbles_high": ("mcube_bubbles_high_{s}_00.fits", "fermipy", "bubbles_high", "Bubbles high"),
    "Sun": ("mcube_Sun_{s}_00.fits", "fermipy", "sun", "Sun"),
    "Moon": ("mcube_Moon_{s}_00.fits", "fermipy", "moon", "Moon"),
    "Loop_I_S1": ("mcube_Loop_I_S1_{s}_00.fits", "fermipy", "loopI_SI", "LIS1"),
    "Loop_I_S2": ("mcube_Loop_I_S2_{s}_00.fits", "fermipy", "loopI_S2", "LIS2"),
}
OPTIONAL = {"GCE_NFW": "include_dm", "bubbles_low": "include_bubbles",
            "bubbles_high": "include_bubbles"}

# Line styles, as in make_counts_plot.py (key, kwargs).
COUNTS_STYLE = [
    ("brem", dict(ls="--", lw=2, color="purple", label="bremss")),
    ("ic", dict(ls="--", lw=2, color="cyan", label="ics")),
    ("pion_decay", dict(ls="--", lw=2, color="cornflowerblue", label="pion decay")),
    ("isodiff", dict(ls="-.", lw=2, color="orange", label="isodiff")),
    ("point_sources", dict(ls=":", lw=2, color="green", label="point sources")),
    ("GCE_NFW", dict(ls="-", marker="D", lw=2, color="purple", label="GCE")),
    ("bubbles_high", dict(ls="-", marker="o", lw=2, color="red", label="Bubble high")),
    ("bubbles_low", dict(ls="-", marker="o", lw=2, color="magenta", label="Bubble low")),
    ("Sun", dict(ls=":", lw=2, color="red", label="Sun")),
    ("Moon", dict(ls=":", lw=2, color="pink", label="Moon")),
    ("Loop_I_S1", dict(ls="-", marker="s", lw=2, color="tan", label="LIS1")),
    ("Loop_I_S2", dict(ls="-", marker="s", lw=2, color="brown", label="LIS2")),
]

# Line styles, as in plot_flux.py (key, kwargs).
FLUX_STYLE = [
    ("pion_decay", dict(color="cornflowerblue", ls="--", label="pion_decay")),
    ("ic", dict(color="cyan", ls="--", label="ics")),
    ("brem", dict(color="purple", ls="--", label="bremss")),
    ("isodiff", dict(color="darkorange", ls="--", label="isodiff")),
    ("point_sources", dict(color="green", ls="--", label="point sources")),
    ("GCE_NFW", dict(color="purple", ls="-", marker="D", label="GCE")),
    ("bubbles_low", dict(color="magenta", marker="o", ls="-", label="Bubbles low")),
    ("bubbles_high", dict(color="red", marker="o", ls="-", label="Bubbles high")),
    ("Sun", dict(color="red", ls=":", label="Sun")),
    ("Moon", dict(color="pink", ls=":", label="Moon")),
    ("Loop_I_S1", dict(color="tan", marker="s", ls="-", label="LIS1")),
    ("Loop_I_S2", dict(color="brown", marker="s", ls="-", label="LIS2")),
]


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


class _Log:
    """Run each step on its own and record the outcome."""

    def __init__(self, path):
        self.path = Path(path)
        self.lines = []

    def step(self, name, func, *args, **kwargs):
        start = time.time()
        print(f"[make_all_plots] {name} ...", flush=True)
        try:
            result = func(*args, **kwargs)
        except Exception:
            text = traceback.format_exc()
            print(f"[make_all_plots] {name} FAILED:\n{text}", flush=True)
            self.lines.append(f"FAILED  {name}\n{text}")
            result = None
        else:
            self.lines.append(f"ok      {name} ({time.time() - start:.1f} s)")
        self.path.write_text("\n".join(self.lines) + "\n")
        return result


def _close_all():
    import matplotlib.pyplot as plt

    plt.close("all")


def _component_files(info):
    """{key: Path} of the cubes present for this run (optional ones only if included)."""
    folders = {"fermipy": Path(info["fermipy_outdir"]), "run": Path(info["run_output"])}
    files = {}
    for key, (pattern, folder, _, _) in COMPONENTS.items():
        if key in OPTIONAL and not info.get(OPTIONAL[key], True):
            continue
        files[key] = folders[folder] / pattern.format(s=info["savename"])
    return files


def _read_cube(path):
    """(data [E, y, x], energy-bounds table, header) of a counts or model cube."""
    from astropy.io import fits

    with fits.open(path) as hdus:
        data = np.asarray(hdus[0].data, dtype=float)
        # The energy table is the second HDU, as calculate_flux assumes. Its
        # name depends on how the cube was written: EBOUNDS (gtbin, Fermipy
        # model maps) or BANDS (gammapy's default, e.g. the point-source cube).
        names = [hdu.name.upper() for hdu in hdus]
        ebounds = hdus[names.index("EBOUNDS") if "EBOUNDS" in names else 1].data
        ebounds = {"E_MIN": np.asarray(ebounds["E_MIN"], float),
                   "E_MAX": np.asarray(ebounds["E_MAX"], float)}
        header = hdus[0].header.copy()
    return data, ebounds, header


# ---------------------------------------------------------------------------
# Fluxes (lat_binned_analysis Plots.calculate_flux, vectorized)
# ---------------------------------------------------------------------------


class _ExposureLookup:
    """Exposure of each cube pixel, as calculate_flux looks it up.

    calculate_flux maps each cube pixel (x, y) to sky coordinates, then to
    exposure-map pixel coordinates, truncates them with int(), and uses the
    geometric mean of exposure planes E and E+1.
    """

    def __init__(self, exposure_file):
        from astropy.io import fits
        from astropy.wcs import WCS

        self.wcs = WCS(str(exposure_file))
        with fits.open(str(exposure_file)) as hdus:
            self.data = np.asarray(hdus[0].data, dtype=float)
        self._cache = {}

    def planes(self, cube_file, shape):
        from astropy.wcs import WCS

        key = (str(cube_file), shape)
        if key not in self._cache:
            ny, nx = shape
            y, x = np.mgrid[0:ny, 0:nx]
            lon, lat = WCS(str(cube_file)).celestial.wcs_pix2world(x, y, 0)
            px, py = self.wcs.celestial.wcs_world2pix(lon, lat, 0)
            ix, iy = np.trunc(px).astype(int), np.trunc(py).astype(int)
            nexp = self.data.shape[0]
            self._cache[key] = np.array(
                [np.sqrt(self.data[k][iy, ix] * self.data[k + 1][iy, ix]) for k in range(nexp - 1)]
            )
        return self._cache[key]


def roi_flux(cube_file, cube, ebounds, exposure, name, filetype, outdir):
    """Write <name>_flux.dat exactly as calculate_flux does; return its columns.

    filetype "ccube": energies in keV (divided by 1000); "fermipy_mcube":
    energies in MeV.
    """
    from astropy.io import ascii
    from astropy.table import Table

    scale = 1000.0 if filetype == "ccube" else 1.0
    e_min, e_max = ebounds["E_MIN"] / scale, ebounds["E_MAX"] / scale
    delta_e = e_max - e_min
    energy = np.sqrt(e_max * e_min)
    time_ = exposure.planes(cube_file, cube.shape[1:])[: cube.shape[0]]
    per_pixel = cube / time_
    counts_flux = per_pixel.sum(axis=(1, 2))
    dn_de = counts_flux / delta_e
    flux = dn_de * energy**2
    raw_counts = cube.sum(axis=(1, 2))

    data = Table(
        {
            "Energy[MeV]": energy,
            "energy_flux[MeV/cm^2/s]": flux,
            "counts_flux[ph/cm^2/s]": counts_flux,
            "differential_flux[ph/cm^2/s/MeV]": dn_de,
            "raw_counts[ph]": raw_counts,
        },
        names=[
            "Energy[MeV]",
            "energy_flux[MeV/cm^2/s]",
            "counts_flux[ph/cm^2/s]",
            "differential_flux[ph/cm^2/s/MeV]",
            "raw_counts[ph]",
        ],
    )
    data.meta["comments"] = [
        "Total counts flux [ph/cm^2/s]: " + str(counts_flux.sum()),
        "Total energy flux [MeV/cm^2/s]: " + str(flux.sum()),
    ]
    ascii.write(data, str(Path(outdir) / f"{name}_flux.dat"), delimiter="\t", overwrite=True)
    return energy, flux


# ---------------------------------------------------------------------------
# Plot steps
# ---------------------------------------------------------------------------


def load_cubes(files):
    """Read every cube once: {key: (data, ebounds, header)} for the files present."""
    cubes = {}
    for key, path in files.items():
        if Path(path).is_file():
            try:
                cubes[key] = _read_cube(path)
            except Exception as error:
                raise RuntimeError(f"Cannot read {path}: {error!r}") from error
        else:
            print(f"[make_all_plots] missing, skipped: {path}")
    for key in ("observed", "model"):
        if key not in cubes:
            raise FileNotFoundError(f"Required cube missing: {files[key]}")
    return cubes


def _energy_axis(cubes):
    """Bin centres and x error bars [MeV] from the counts cube (as the scripts do)."""
    eb = cubes["observed"][1]
    centre = np.sqrt(eb["E_MIN"] * eb["E_MAX"]) / 1000.0
    xerr = np.array([centre - eb["E_MIN"] / 1000.0, eb["E_MAX"] / 1000.0 - centre])
    return centre, xerr


def counts_plot(cubes, plot_dir):
    """counts_spectra.png/.pdf and counts_spectrum.csv (make_counts_plot.py)."""
    import matplotlib.pyplot as plt

    energy, xerr = _energy_axis(cubes)
    counts = {key: c[0].sum(axis=(1, 2)) for key, c in cubes.items()}
    data, model = counts["observed"], counts["model"]
    resid = (data - model) / model
    resid_err = np.sqrt(data) / model

    columns = {"energy_MeV": energy, "e_min_MeV": energy - xerr[0], "e_max_MeV": energy + xerr[1],
               "data": data, "model": model, "frac_resid": resid, "frac_resid_err": resid_err}
    columns.update({f"counts_{k}": v for k, v in counts.items() if k not in ("observed", "model")})
    np.savetxt(Path(plot_dir) / "counts_spectrum.csv", np.column_stack(list(columns.values())),
               delimiter=",", header=",".join(columns), comments="")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    ax1.loglog(energy, data, marker="s", color="black", ls="", label="total observed")
    ax1.errorbar(energy, data, xerr=xerr, yerr=np.sqrt(data), marker="s", color="black", ls="",
                 label="_nolabel_")
    ax1.loglog(energy, model, ls="-", lw=2, color="black", label="total predicted")
    for key, style in COUNTS_STYLE:
        if key in counts:
            ax1.loglog(energy, counts[key], **style)
    ax1.set_xlim(3e2, 0.1e6)
    ax1.set_ylim(1e2, 1e7)
    ax1.set_ylabel("Counts", fontsize=16)
    ax1.tick_params(axis="both", labelsize=14)
    ax1.legend(loc=1, ncol=2, frameon=False, fontsize=10)
    ax1.grid(ls="--", color="grey", alpha=0.2)

    ax2.semilogx(energy, resid, ls="", color="black", marker="s")
    ax2.errorbar(energy, resid, xerr=xerr, yerr=resid_err, ls="", color="black", marker="s")
    ax2.axhline(y=0, ls="-", color="grey", lw=1)
    ax2.set_ylim(-1, 1)
    ax2.set_xlabel("Energy [MeV]", fontsize=16)
    ax2.set_ylabel("(D - M)/M", fontsize=16)
    ax2.tick_params(axis="both", labelsize=14)
    ax2.grid(ls="--", color="grey", alpha=0.2)
    fig.savefig(Path(plot_dir) / "counts_spectra.png")
    fig.savefig(Path(plot_dir) / "counts_spectra.pdf")
    plt.close(fig)


def flux_tables(cubes, files, exposure_file, plot_dir):
    """<name>_flux.dat for every component (plot_flux.py / calculate_flux)."""
    exposure = _ExposureLookup(exposure_file)
    fluxes = {}
    for key, (data, ebounds, _) in cubes.items():
        name = COMPONENTS[key][2]
        filetype = "ccube" if key == "observed" else "fermipy_mcube"
        fluxes[key] = roi_flux(files[key], data, ebounds, exposure, name, filetype, plot_dir)
    return fluxes


def flux_plot(cubes, fluxes, savename, plot_dir):
    """flux_spectra_<savename>.png/.pdf (plot_flux.py)."""
    import matplotlib.pyplot as plt

    energy, xerr = _energy_axis(cubes)
    data = cubes["observed"][0].sum(axis=(1, 2))
    model = cubes["model"][0].sum(axis=(1, 2))
    resid = (data - model) / model
    resid_err = np.sqrt(data) / model

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 9), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})
    e, f = fluxes["observed"]
    ax1.loglog(e, f, color="black", ls="", marker="s", label="total observed")
    ax1.errorbar(e, f, xerr=xerr, color="black", ls="", marker="s", label="_NONE_")
    e, f = fluxes["model"]
    ax1.loglog(e, f, color="black", ls="-", label="total model")
    for key, style in FLUX_STYLE:
        if key in fluxes:
            ax1.loglog(*fluxes[key], **style)
    ax1.set_xlim(3e2, 0.1e6)
    ax1.set_ylim(1e-6, 1e0)
    ax1.set_ylabel(r"$\mathrm{E^2 dN/dE \  [MeV \ cm^{-2} \ s^{-1}]}$", fontsize=16)
    ax1.tick_params(axis="both", labelsize=14)
    ax1.legend(loc=1, ncol=3, frameon=False, fontsize=12)
    ax1.grid(ls="--", color="grey", alpha=0.2)

    ax2.semilogx(energy, resid, ls="", color="black", marker="s")
    ax2.errorbar(energy, resid, xerr=xerr, yerr=resid_err, ls="", color="black", marker="s")
    ax2.axhline(y=0, ls="-", color="grey", lw=1)
    ax2.set_ylim(-0.1, 0.1)
    ax2.set_xlabel("Energy [MeV]", fontsize=16)
    ax2.set_ylabel("(D - M)/M", fontsize=16)
    ax2.tick_params(axis="both", labelsize=14)
    ax2.grid(ls="--", color="grey", alpha=0.2)
    fig.savefig(Path(plot_dir) / f"flux_spectra_{savename}.png")
    fig.savefig(Path(plot_dir) / f"flux_spectra_{savename}.pdf")
    plt.close(fig)


def spatial_plots(files, cubes, inputs_yaml, plot_dir):
    """count_map.png, model_map.png and residuals.png (plot_spatial.py)."""
    import matplotlib.pyplot as plt
    from astropy.wcs import WCS
    from lat_binned_analysis import Plots

    instance_plots = Plots(str(inputs_yaml))
    instance_plots.plot_skymap(str(files["observed"]), str(Path(plot_dir) / "count_map.png"),
                               input_dim=3, make_sum=True, vmin=100, title="Data")
    _close_all()
    instance_plots.plot_skymap(str(files["model"]), str(Path(plot_dir) / "model_map.png"),
                               input_dim=3, make_sum=True, vmin=100, title="Model")
    _close_all()

    data, _, header = cubes["observed"]
    image = data.sum(axis=0) - cubes["model"][0].sum(axis=0)
    wcs_2d = WCS(header, naxis=2)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection=wcs_2d)
    this_max = np.amax(image)
    im = ax.imshow(image, origin="lower", cmap="seismic", vmin=-1 * this_max, vmax=this_max,
                   interpolation="none")
    cbar = plt.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("data - model [counts]", fontsize=12)
    ax.set_xlabel(r"Galactic Longitude [$\circ$]", fontsize=14)
    ax.set_ylabel(r"Galactic Latitude [$\circ$]", fontsize=14)
    ax.grid(color="grey", ls=":", alpha=0.5)
    ax.tick_params(axis="both", which="major", labelsize=14)
    plt.tight_layout()
    fig.savefig(Path(plot_dir) / "residuals.png")
    plt.close(fig)


def _band_title(ebounds_kev, start, stop):
    """'300 MeV - 1.28 GeV' for energy bins start..stop-1 (EBOUNDS in keV)."""

    def fmt(e_mev):
        return f"{e_mev:.3g} MeV" if e_mev < 1000 else f"{e_mev / 1000:.3g} GeV"

    stop = min(stop, len(ebounds_kev["E_MAX"]))
    return f"{fmt(ebounds_kev['E_MIN'][start] / 1000)} - {fmt(ebounds_kev['E_MAX'][stop - 1] / 1000)}"


def residual_maps(files, cubes, inputs_yaml, bins, plot_dir):
    """3d residual maps (custom_client.py) in three forms, per energy range.

    lat_binned_analysis writes the FITS maps; each is plotted on its own
    (plot_skymap) and the three forms are shown together in
    residual_maps_<i>-<j>.png.
    """
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.wcs import WCS
    from lat_binned_analysis import Analysis, Plots

    instance = Analysis(str(inputs_yaml))
    instance_plots = Plots(str(inputs_yaml))
    ebounds = cubes["observed"][1]
    written = {}
    for mode, _, label in RESIDUAL_MODES:
        outfiles = instance.residual_map2(None, None, str(files["observed"]), str(files["model"]),
                                          mode=mode, bins=bins, outdir=str(plot_dir))
        for (start, stop), path in zip(bins, outfiles):
            written[(mode, start, stop)] = path
            instance_plots.plot_skymap(
                path, str(Path(plot_dir) / f"residuals_{mode}_{start}-{stop}.png"),
                plot_type="residuals", colorbar_label=label,
                title=_band_title(ebounds, start, stop))
            _close_all()

    # One figure per energy range with the three forms side by side, in the
    # style of plot_skymap (seismic, symmetric colour scale). The scale is the
    # 99.5th percentile of |residual|, so a few extreme pixels do not hide the rest.
    for start, stop in bins:
        fig = plt.figure(figsize=(24, 6.6))
        for j, (mode, _, label) in enumerate(RESIDUAL_MODES):
            with fits.open(written[(mode, start, stop)]) as hdus:
                image = np.asarray(hdus[0].data, dtype=float)
                wcs_2d = WCS(hdus[0].header, naxis=2)
            ax = fig.add_subplot(1, 3, j + 1, projection=wcs_2d)
            finite = np.abs(image[np.isfinite(image)])
            limit = float(np.percentile(finite, 99.5)) if finite.size else 1.0
            limit = limit if limit > 0 else 1.0
            im = ax.imshow(image, origin="lower", cmap="seismic", vmin=-limit, vmax=limit,
                           interpolation="none")
            cbar = plt.colorbar(im, ax=ax, pad=0.02, fraction=0.046, extend="both")
            cbar.set_label(label, fontsize=12)
            cbar.ax.tick_params(labelsize=12)
            ax.set_xlabel(r"Galactic Longitude [$\circ$]", fontsize=14)
            ax.set_ylabel(r"Galactic Latitude [$\circ$]", fontsize=14)
            ax.grid(color="grey", ls=":", alpha=0.5)
            ax.tick_params(axis="both", which="major", labelsize=14)
        fig.subplots_adjust(wspace=0.55)
        fig.suptitle(_band_title(ebounds, start, stop), fontsize=16)
        fig.savefig(Path(plot_dir) / f"residual_maps_{start}-{stop}.png", bbox_inches="tight")
        plt.close(fig)


def gce_plot(fit_results, nfw_file, plot_dir):
    """GCE_Ajello2016_comparison.pdf/.png and gce_spectrum.csv."""
    sys.path.insert(0, str(REPO / "Main"))
    from compare_gce_ajello2016 import plot_gce_comparison

    fig = plot_gce_comparison(fit_results, nfw_file,
                              output=Path(plot_dir) / "GCE_Ajello2016_comparison.pdf",
                              data_output=Path(plot_dir) / "gce_spectrum.csv")
    fig.savefig(Path(plot_dir) / "GCE_Ajello2016_comparison.png", dpi=150, bbox_inches="tight")
    _close_all()


def bubbles_plot(plot_dir, label):
    """bubbles_comparison.png/.pdf etc. (compare_bubbles.py)."""
    sys.path.insert(0, str(PLOTS_DIR / "compare_bubbles"))
    import compare_bubbles

    compare_bubbles.main(low_file=Path(plot_dir) / "bubbles_low_flux.dat",
                         high_file=Path(plot_dir) / "bubbles_high_flux.dat",
                         label=label, outdir=plot_dir, prefix="bubbles_", show=False)


def loopI_plot(allsky_results_dir, plot_dir, label):
    """loopI_comparison.png/.pdf etc. (compare_loopI.py, from the all-sky fit)."""
    sys.path.insert(0, str(PLOTS_DIR / "compare_loopI"))
    import compare_loopI

    compare_loopI.main(allsky_dir=Path(allsky_results_dir), label=label, outdir=plot_dir,
                       prefix="loopI_", show=False)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------


def make_all_plots(
    plot_dir,
    savename,
    fermipy_outdir,
    run_output="output",
    allsky_results_dir=None,
    include_dm=True,
    include_bubbles=True,
    fit_results=None,
    nfw_file=None,
    residual_bins=RESIDUAL_BINS,
    inputs_yaml=PLOTS_DIR / "inputs.yaml",
    label="This fit",
):
    """Make all plots of one fit in `plot_dir` (created if needed).

    Parameters
    ----------
    plot_dir : str or Path
        Output folder, normally Plots/<name>.
    savename : str
        `savename` of run_3ML_master.py (part of the cube file names).
    fermipy_outdir : str or Path
        Fermipy working directory of the fit (counts, exposure, model cubes).
    run_output : str or Path
        output/ folder of the run directory (point-source cube).
    allsky_results_dir : str or Path, optional
        Results folder of the all-sky fit (Isotropic/<outdir>/bestfit_<hash>)
        for the Loop I comparison. Skipped if None.
    include_dm, include_bubbles : bool
        As in run_3ML_master.py; the GCE and bubble plots need them.
    fit_results : threeML results, optional
        ``jl.results``, for the GCE comparison; saved as
        threeML_analysis_results.fits so the plot can be made again later.
    nfw_file : str or Path, optional
        NFW template of the fit (GCE comparison).
    residual_bins : list of [start, stop]
        Energy-bin index ranges of the 3d residual maps.
    inputs_yaml : str or Path
        lat_binned_analysis configuration (Plots/inputs.yaml).
    label : str
        Name of this fit in the comparison legends.
    """
    import matplotlib.pyplot as plt

    plt.switch_backend("Agg")
    plot_dir = Path(plot_dir).resolve()
    plot_dir.mkdir(parents=True, exist_ok=True)
    info = dict(
        savename=savename,
        label=label,
        fermipy_outdir=str(Path(fermipy_outdir).resolve()),
        run_output=str(Path(run_output).resolve()),
        allsky_results_dir=None if allsky_results_dir is None else str(Path(allsky_results_dir).resolve()),
        include_dm=bool(include_dm),
        include_bubbles=bool(include_bubbles),
        nfw_file=None if nfw_file is None else str(Path(nfw_file).resolve()),
        residual_bins=[list(map(int, b)) for b in residual_bins],
        inputs_yaml=str(Path(inputs_yaml).resolve()),
        created=time.strftime("%Y-%m-%d %H:%M:%S"),
    )
    (plot_dir / "run_info.json").write_text(json.dumps(info, indent=2) + "\n")
    log = _Log(plot_dir / "make_all_plots.log")

    files = _component_files(info)
    cubes = log.step("read cubes", load_cubes, files)
    if cubes is not None:
        log.step("counts spectrum", counts_plot, cubes, plot_dir)
        fluxes = log.step("flux tables", flux_tables, cubes, files,
                          Path(info["fermipy_outdir"]) / "bexpmap_00.fits", plot_dir)
        if fluxes is not None:
            log.step("flux spectrum", flux_plot, cubes, fluxes, savename, plot_dir)
        log.step("spatial maps", spatial_plots, files, cubes, info["inputs_yaml"], plot_dir)
        log.step("3d residual maps", residual_maps, files, cubes, info["inputs_yaml"],
                 info["residual_bins"], plot_dir)
        _close_all()

    if include_dm:
        if fit_results is not None:
            results_file = plot_dir / "threeML_analysis_results.fits"
            log.step("save 3ML results", fit_results.write_to, str(results_file), overwrite=True)
        if fit_results is not None and nfw_file is not None:
            log.step("GCE comparison", gce_plot, fit_results, nfw_file, plot_dir)
        else:
            log.lines.append("skipped GCE comparison (needs fit_results and nfw_file)")
    if include_bubbles:
        log.step("bubbles comparison", bubbles_plot, plot_dir, label)
    if allsky_results_dir is not None:
        log.step("Loop I comparison", loopI_plot, info["allsky_results_dir"], plot_dir, label)
    else:
        log.lines.append("skipped Loop I comparison (no allsky_results_dir)")
    _close_all()
    log.path.write_text("\n".join(log.lines) + "\n")
    failed = [line for line in log.lines if line.startswith("FAILED")]
    print(f"[make_all_plots] done: {len(failed)} failed step(s); see {log.path}")
    return info


def _from_run_info(plot_dir):
    """Redo the plots of a folder made earlier (uses its run_info.json)."""
    plot_dir = Path(plot_dir).resolve()
    info = json.loads((plot_dir / "run_info.json").read_text())
    fit_results = None
    results_file = plot_dir / "threeML_analysis_results.fits"
    if info["include_dm"] and results_file.is_file():
        try:
            from threeML import load_analysis_results

            fit_results = load_analysis_results(str(results_file))
        except Exception as error:  # e.g. threeML not installed in this environment
            print(f"Cannot load {results_file} ({error}); the GCE comparison is skipped.")
    return make_all_plots(
        plot_dir,
        info["savename"],
        info["fermipy_outdir"],
        run_output=info["run_output"],
        allsky_results_dir=info["allsky_results_dir"],
        include_dm=info["include_dm"],
        include_bubbles=info["include_bubbles"],
        fit_results=fit_results,
        nfw_file=info["nfw_file"],
        residual_bins=info["residual_bins"],
        inputs_yaml=info["inputs_yaml"],
        label=info.get("label", "This fit"),
    )


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    _from_run_info(sys.argv[1])
