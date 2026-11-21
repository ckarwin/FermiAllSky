"""Compare the Loop I spectrum of the all-sky fit with published Fermi-LAT results.

Run: python compare_loopI.py
Requires: numpy, matplotlib, pyyaml (and LoopI/shell_geometry.py from this repository).
Edit the settings below to compare a different fit. See README.md.

The fitted Loop I sky is rebuilt from the per-shell spectra written by the
all-sky fit (Isotropic/), then averaged over the region that each published
result refers to. Each published result is one YAML file in reference/.

The functions are also used by make_all_plots.py (one fit, called after
run_3ML_master.py) and compare_fits.py (several fits overlaid).
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
sys.path.insert(0, str(SCRIPT_DIR.parent))
from plot_style import GRID, LABEL_SIZE, TICK_SIZE  # noqa: E402
from shell_geometry import SHELLS, path_length  # noqa: E402

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

# Results folder of the all-sky fit (Isotropic/<outdir>/bestfit_<hash> or
# fit_<hash>). It must contain loopI_s1_spectrum.txt and loopI_s2_spectrum.txt
# (Loop I mode "separate") or loopI_spectrum.txt (mode "combined"). Relative
# paths are taken from this folder; absolute paths also work.
ALLSKY_RESULTS_DIR = "example_allsky"
FIT_LABEL = "This work"

# Published results to plot: YAML files in reference/ (all of them by default).
REFERENCE_FILES = sorted(REFERENCE_DIR.glob("*.yaml"))

# Energy range used for the mean fit/published ratios [MeV]. A reference file
# can narrow it with its own `ratio_band_mev`.
RATIO_BAND_MEV = (300.0, 30000.0)

# Y range of the ratio panel. Points outside (e.g. the steep high-energy tail)
# are cut off.
RATIO_YLIM = (0.1, 2.0)

# Template grid (same as LoopI/build_templates.py): full-sky CAR, 0.5 deg pixels.
PIXEL_DEG = 0.5

# Reference region of the Loop I templates: the summed S1 + S2 map averages to
# 1 over all longitudes at |b| > TEMPLATE_REFERENCE_B_DEG (LoopI/build_templates.py).
TEMPLATE_REFERENCE_B_DEG = 10.0


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------


def load_reference(path):
    """One published result: metadata plus energy [MeV], E^2 I and its 1-sigma range."""
    data = yaml.safe_load(Path(path).read_text())
    for key in ("energy", "e2_intensity", "e2_intensity_low", "e2_intensity_high"):
        data[key] = np.asarray(data[key], float)
    return data


def load_fit_spectra(results_dir):
    """Per-shell spectral coefficients from the all-sky fit.

    Returns (energy [MeV], coefficients[n_shell_terms, n_energy], weights) where
    the fitted Loop I intensity is sum_k coefficient_k(E) * template_k(l, b)
    [ph cm^-2 s^-1 sr^-1 MeV^-1], and template_k is S1, S2 or S1 + S2 (rows of
    `weights` select the shells that each term multiplies).
    """
    results_dir = Path(results_dir)
    if not results_dir.is_absolute():
        results_dir = SCRIPT_DIR / results_dir
    separate = [results_dir / f"loopI_s{i}_spectrum.txt" for i in (1, 2)]
    combined = results_dir / "loopI_spectrum.txt"
    if all(p.is_file() for p in separate):
        tables = [np.loadtxt(p) for p in separate]
        weights = np.eye(2)
    elif combined.is_file():
        tables = [np.loadtxt(combined)]
        weights = np.ones((1, 2))
    else:
        raise FileNotFoundError(
            f"No Loop I spectrum files in {results_dir}: expected loopI_s1_spectrum.txt and "
            "loopI_s2_spectrum.txt, or loopI_spectrum.txt (written by Isotropic/fit_isotropic.py)"
        )
    energy = tables[0][:, 0]
    if not all(np.allclose(t[:, 0], energy) for t in tables):
        raise ValueError("Loop I spectrum files must use the same energies")
    return energy, np.array([t[:, 1] for t in tables]), weights


def shell_templates():
    """S1 and S2 templates exactly as LoopI/build_templates.py makes them.

    Returns (templates[2, nb, nl], l, b, pixel solid angle): line-of-sight
    path lengths, scaled so that S1 + S2 averages to 1 over |b| > 10 deg.
    """
    nb, nl = int(round(180 / PIXEL_DEG)), int(round(360 / PIXEL_DEG))
    b = (np.arange(nb) - (nb - 1) / 2) * PIXEL_DEG
    l = -(np.arange(nl) - (nl - 1) / 2) * PIXEL_DEG
    b, l = np.meshgrid(b, l, indexing="ij")
    l = (l + 180) % 360 - 180
    half = np.deg2rad(PIXEL_DEG / 2)
    omega = np.deg2rad(PIXEL_DEG) * (np.sin(np.deg2rad(b) + half) - np.sin(np.deg2rad(b) - half))
    raw = np.array([path_length(l, b, shell) for shell in SHELLS])
    ref = np.abs(b) > TEMPLATE_REFERENCE_B_DEG
    scale = np.sum(raw.sum(axis=0)[ref] * omega[ref]) / omega[ref].sum()
    return raw / scale, l, b, omega


# ---------------------------------------------------------------------------
# Fitted Loop I averaged over a published region
# ---------------------------------------------------------------------------


def region_pixels(region, solid_angle, fitted_map, b, omega):
    """Boolean pixel mask for a published region definition.

    * ``latitude``: all longitudes at |b| > b_min_deg.
    * ``flux_per_area``: the same pixels; the intensity is then rescaled to
      the total flux over them divided by ``solid_angle_sr`` (see
      `region_scale`).
    * ``footprint``: the brightest `solid_angle` sr of the fitted Loop I map at
      |b| > b_min_deg, standing in for a published Loop I template area.
    """
    high = np.abs(b) > region["b_min_deg"]
    if region["type"] in ("latitude", "flux_per_area"):
        return high
    if region["type"] != "footprint":
        raise ValueError(f"Unknown region type {region['type']!r}")
    candidates = np.flatnonzero(high)
    order = candidates[np.argsort(fitted_map.ravel()[candidates])[::-1]]
    keep = order[np.cumsum(omega.ravel()[order]) <= solid_angle]
    mask = np.zeros(fitted_map.size, bool)
    mask[keep] = True
    return mask.reshape(fitted_map.shape)


def mean_intensity(energy, coefficients, term_maps, mask, omega):
    """E^2 x mean fitted intensity over `mask`, per term and in total [MeV cm^-2 s^-1 sr^-1]."""
    area = omega[mask].sum()
    means = np.array([(m * omega)[mask].sum() / area for m in term_maps])
    parts = energy**2 * coefficients * means[:, None]
    return parts, parts.sum(axis=0)


def interp_loglog(x, xp, fp):
    """Log-log interpolation; NaN outside xp."""
    out = np.exp(np.interp(np.log(x), np.log(xp), np.log(fp)))
    return np.where((x >= xp[0]) & (x <= xp[-1]), out, np.nan)


def region_scale(region, mask, omega):
    """Factor applied to the mean intensity over `mask`.

    For ``flux_per_area`` the published value is the total Loop I flux at
    |b| > b_min_deg divided by ``solid_angle_sr``; otherwise it is the mean.
    """
    if region["type"] == "flux_per_area":
        return omega[mask].sum() / region["solid_angle_sr"]
    return 1.0


def region_text(region):
    """Short plot label for a region definition."""
    if region["type"] == "latitude":
        return f"|b| > {region['b_min_deg']:g}°"
    if region["type"] == "flux_per_area":
        return f"flux at |b| > {region['b_min_deg']:g}° / {region['solid_angle_sr']:g} sr"
    return f"brightest {region['solid_angle_sr']:g} sr at |b| > {region['b_min_deg']:g}° (assumed)"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


# Published results: one fixed colour and marker each, in file order.
REF_COLOURS = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e", "#8c564b", "#e377c2"]
REF_MARKERS = ["o", "s", "^", "D", "v", "P", "X"]
REGION_STYLES = ["-", "--", "-.", ":"]


def load_references(files=None):
    """All published results (REFERENCE_FILES by default)."""
    references = [load_reference(p) for p in (REFERENCE_FILES if files is None else files)]
    if not references:
        raise FileNotFoundError(f"No reference YAML files in {REFERENCE_DIR}")
    return references


def fit_curves(allsky_dir, references):
    """The fitted Loop I evaluated over each published region.

    Returns (energy [MeV], curves, rows): one curve dict per reference
    (total, parts, band, ratio, ...) and one summary row per reference.
    """
    energy, coefficients, weights = load_fit_spectra(allsky_dir)
    shells, l, b, omega = shell_templates()
    term_maps = np.tensordot(weights, shells, axes=1)
    # Map used to rank pixels for "footprint" regions. All terms share one
    # spectral shape, so the band-summed coefficients give the map's shape.
    fitted_map = np.tensordot((coefficients * energy**2).sum(axis=1), term_maps, axes=1)

    rows, curves = [], []
    for ref in references:
        region = ref["region"]
        mask = region_pixels(region, region.get("solid_angle_sr"), fitted_map, b, omega)
        parts, total = mean_intensity(energy, coefficients, term_maps, mask, omega)
        scale = region_scale(region, mask, omega)
        parts, total = parts * scale, total * scale
        band = None
        if region["type"] == "footprint" and "solid_angle_range_sr" in region:
            band = np.sort([
                mean_intensity(energy, coefficients, term_maps,
                               region_pixels(region, w, fitted_map, b, omega), omega)[1]
                for w in region["solid_angle_range_sr"]
            ], axis=0)
        fit_at_ref = interp_loglog(ref["energy"], energy, total)
        ratio = fit_at_ref / ref["e2_intensity"]
        # A reference may narrow the band (e.g. to the energies where it was fitted).
        lo_e, hi_e = ref.get("ratio_band_mev", RATIO_BAND_MEV)
        lo_e, hi_e = max(lo_e, RATIO_BAND_MEV[0]), min(hi_e, RATIO_BAND_MEV[1])
        in_band = (ref["energy"] >= lo_e) & (ref["energy"] <= hi_e) & np.isfinite(ratio)
        mean_ratio = float(np.exp(np.mean(np.log(ratio[in_band])))) if in_band.any() else np.nan
        curves.append(dict(ref=ref, total=total, parts=parts, band=band, ratio=ratio,
                           fit_at_ref=fit_at_ref, mean_ratio=mean_ratio, band_mev=(lo_e, hi_e)))
        rows.append({
            "reference": ref["label"],
            "region": region_text(region),
            "region_solid_angle_sr": float(omega[mask].sum()),
            "S1_share": float(np.mean(parts[0] / total)) if len(parts) == 2 else np.nan,
            "mean_fit_over_published": mean_ratio,
            "ratio_band_MeV": f"{lo_e:g}-{hi_e:g}",
            "n_points_in_ratio_band": int(in_band.sum()),
        })
    return energy, curves, rows


def draw_references(ax, references, axr=None, curves=None):
    """Published points on `ax` (and their fit/published ratios on `axr`).

    Points outside a result's own ratio band (e.g. not fitted bin by bin)
    are drawn hollow.
    """
    for i, ref in enumerate(references):
        col, mk = REF_COLOURS[i % len(REF_COLOURS)], REF_MARKERS[i % len(REF_MARKERS)]
        e, y, lo, hi = ref["energy"], ref["e2_intensity"], ref["e2_intensity_low"], ref["e2_intensity_high"]
        own = ref.get("ratio_band_mev", (0, np.inf))
        inside = (e >= own[0]) & (e <= own[1])
        for sel, face, label in ((inside, col, ref["label"]), (~inside, "white", None)):
            if not sel.any():
                continue
            ax.errorbar(e[sel], y[sel], yerr=[(y - lo)[sel], (hi - y)[sel]], fmt=mk, color=col,
                        mfc=face, ms=5.5, capsize=2, lw=1.2, label=label)
            if axr is not None:
                c = curves[i]
                rerr = np.array([c["ratio"] - c["fit_at_ref"] / hi, c["fit_at_ref"] / lo - c["ratio"]])
                axr.errorbar(e[sel], c["ratio"][sel], yerr=rerr[:, sel], fmt=mk, color=col, mfc=face,
                             ms=4.5, capsize=2, lw=1.2)
        if axr is not None:
            axr.hlines(curves[i]["mean_ratio"], *curves[i]["band_mev"], color=col, ls="--", lw=1.2)
    if any(np.any(np.asarray(r["energy"]) > r.get("ratio_band_mev", (0, np.inf))[1])
           for r in references):
        ax.plot([], [], "o", color="0.4", mfc="white", ms=5.5,
                label="hollow: not fitted bin by bin in that paper")


def draw_leptonic_curve(ax):
    """The leptonic model whose shape the fit uses."""
    curve = np.loadtxt(LOOPI_DIR / "leptonic_curve.txt")
    ax.plot(curve[:, 0], curve[:, 2], color="0.5", lw=1.2, ls=":",
            label="Leptonic model of Ackermann et al. 2025 (shape used in the fit)")


def style_spectrum_axes(ax, xlabel=True):
    """Axis ranges, labels and grid of the spectrum panel."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(80, 1.2e5)
    ax.set_ylim(1e-6, 3e-3)
    ax.set_ylabel(r"$E^2\,I(E)$ [MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]", fontsize=LABEL_SIZE)
    if xlabel:
        ax.set_xlabel("Energy [MeV]", fontsize=LABEL_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_SIZE)
    ax.grid(**GRID)


def main(allsky_dir=ALLSKY_RESULTS_DIR, label=FIT_LABEL, outdir=SCRIPT_DIR, prefix="", show=False):
    """Make the comparison for one all-sky fit.

    Outputs go to `outdir`, with `prefix` before each file name:
    comparison.png, converted_fit.csv and summary.csv.
    """
    outdir = Path(outdir)
    references = load_references()
    energy, curves, rows = fit_curves(allsky_dir, references)

    # --- Tables -----------------------------------------------------------
    header = "energy_MeV\t" + "\t".join(f"[{c['ref']['label']}]" for c in curves)
    np.savetxt(outdir / f"{prefix}converted_fit.csv",
               np.column_stack([energy] + [c["total"] for c in curves]),
               header=header + "  (E^2 I of the fit averaged over each published region, "
               "MeV cm^-2 s^-1 sr^-1)", delimiter="\t")
    ratios = np.array([r["mean_fit_over_published"] for r in rows])
    with open(outdir / f"{prefix}summary.csv", "w") as handle:
        keys = list(rows[0])
        handle.write(",".join(keys) + "\n")
        for r in rows:
            handle.write(",".join(f'"{r[k]}"' if isinstance(r[k], str) else f"{r[k]}" for k in keys)
                         + "\n")
        handle.write(f"# all-sky results: {allsky_dir}\n")
        handle.write(f"# ratio band [MeV]: {RATIO_BAND_MEV[0]:g}-{RATIO_BAND_MEV[1]:g}\n")
        finite = ratios[np.isfinite(ratios)]
        if finite.size:
            handle.write(f"# range of mean ratios over references: {finite.min():.3f}-{finite.max():.3f}\n")

    # --- Plot ---------------------------------------------------------------
    # The fit: one black curve per distinct region definition.
    fig, (ax, axr) = plt.subplots(2, 1, figsize=(9, 9), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1.2], "hspace": 0.05})
    draw_references(ax, references, axr=axr, curves=curves)
    drawn = {}
    for c in curves:
        key = region_text(c["ref"]["region"])
        if key not in drawn:
            drawn[key] = REGION_STYLES[len(drawn) % len(REGION_STYLES)]
            if c["band"] is not None:
                ax.fill_between(energy, *c["band"], color="0.6", alpha=0.3, lw=0)
            ax.plot(energy, c["total"], color="k", ls=drawn[key], lw=2, label=f"{label}, {key}")
    draw_leptonic_curve(ax)
    style_spectrum_axes(ax, xlabel=False)
    ax.legend(fontsize=9, loc="lower left", frameon=False)
    axr.axhline(1, color="0.3", lw=0.8)
    axr.set_yscale("log")
    axr.set_ylim(*RATIO_YLIM)
    axr.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
    axr.yaxis.set_minor_formatter(plt.FuncFormatter(
        lambda v, _: f"{v:g}" if round(v, 3) in (0.2, 0.3, 0.5, 2.0) else ""))
    axr.set_ylabel("fit / published", fontsize=LABEL_SIZE)
    axr.set_xlabel("Energy [MeV]", fontsize=LABEL_SIZE)
    axr.tick_params(axis="both", which="both", labelsize=TICK_SIZE)
    axr.grid(**GRID)
    axr.text(0.01, 0.05, "dashed: mean ratio over the energies it is computed for",
             transform=axr.transAxes, fontsize=10, color="0.3")
    fig.savefig(outdir / f"{prefix}comparison.png", dpi=150, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    for r in rows:
        print(f"{r['reference']}: mean fit/published = {r['mean_fit_over_published']:.3f} "
              f"({r['region']})")
    print(f"Wrote {prefix}comparison.png, {prefix}converted_fit.csv, "
          f"{prefix}summary.csv in {outdir}")
    return rows


if __name__ == "__main__":
    main()
