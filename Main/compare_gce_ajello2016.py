"""Compare the fitted GCE_NFW spectrum with Ajello et al. (2016), Figure 13.

Called by ``run_3ML_master.py`` after the fit when ``include_dm = True``:

    from compare_gce_ajello2016 import plot_gce_comparison
    plot_gce_comparison(jl.results, nfw_path)

It writes ``GCE_Ajello2016_comparison.pdf``: the best-fit E^2 dN/dE of the
``GCE_NFW`` source (exponentially cutoff power law) with a pointwise 68% band,
over the four Galactic-centre-excess envelopes of Ajello et al. 2016, ApJ 819,
44, Figure 13 (pulsar/OB-star cosmic-ray source models, intensity- or
index-scaled IEMs).

How the two are put on the same footing
---------------------------------------
Let T be the unit-integral NFW-squared template and f_R the fraction of T
inside a Galactic rectangle R (``comparison_roi``). The fitted spectral
amplitude refers to the whole template footprint, so

    F_fit(R)   = f_R * F_fit(footprint)
    F_paper(R) = f_R / f_15x15 * F_paper(15 x 15 deg)

This is a pure aperture conversion assuming the same gamma = 1 NFW-squared
morphology; it does not refit data or reproduce the paper's IEMs. The default
R is the paper's 15 x 15 deg region (plane included); use (-20, 20, -20, 20)
for the 40 x 40 deg analysis region.

Caveats
-------
* The paper envelopes are graphical estimates read from the hatching of the
  published PDF (page 19) at 41 energies (1.2-85 GeV), accurate to several
  percent (worse near steep or clipped edges). They are envelope edges, not
  published best-fit curves or author-supplied data. Lower edges that leave
  the paper's plotted range are NaN; no unplotted tails are invented.
* The fit band is a pointwise 16th-84th percentile band from threeML's aligned
  samples of (K, index, xc), so their correlations are kept. For a
  maximum-likelihood fit these samples come from the fitted covariance
  (Gaussian approximation), not a profile likelihood. Figure 13's envelopes
  come from normalization and index uncertainties, so the band definitions
  differ, and the spread among the four IEMs reflects foreground-model
  systematics that a statistical band does not include.
* Units follow the paper: MeV and MeV cm^-2 s^-1 for E^2 dN/dE
  (dN/dE per MeV = 1000 x dN/dE per keV).

Reference: Ajello et al. 2016, ApJ 819, 44, https://doi.org/10.3847/0004-637X/819/1/44
"""

import warnings

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel

# Embedded estimates of Fig.13 ENVELOPE EDGES, obtained from the published
# PDF's vector hatching. Columns: E[MeV], lower and upper E^2 dN/dE[MeV/cm2/s]
# for its 15x15 deg ROI. Approximate graphical data, NOT author-supplied tables.
# Hatching thickness/spacing limits accuracy to several percent (worse near
# steep edges). Black/green stop before the lower edge leaves the plotted range;
# no unseen tails are invented. No published best-fit center is implied.
FIG13_ENVELOPES = {
    "red": np.array(
        [
            [1200, 0.00017325, 0.00021024],
            [1334.86, 0.00017897, 0.00022058],
            [1484.89, 0.00018475, 0.00023128],
            [1651.77, 0.00019055, 0.00024227],
            [1837.41, 0.00019628, 0.00025347],
            [2043.91, 0.00020143, 0.00026419],
            [2273.61, 0.00020669, 0.00027532],
            [2529.14, 0.00021173, 0.00028646],
            [2813.38, 0.00021689, 0.00029803],
            [3129.57, 0.00022146, 0.00030908],
            [3481.29, 0.00022607, 0.00032044],
            [3872.55, 0.00022945, 0.00033034],
            [4307.77, 0.00023248, 0.00033989],
            [4791.91, 0.00023558, 0.00034984],
            [5330.46, 0.00023878, 0.00036018],
            [5929.53, 0.00024085, 0.000369],
            [6595.93, 0.00024143, 0.00037569],
            [7337.23, 0.00024111, 0.00038106],
            [8161.84, 0.00023845, 0.00038272],
            [9079.13, 0.0002355, 0.00038387],
            [10099.5, 0.00023019, 0.00038108],
            [11234.6, 0.00022495, 0.00037821],
            [12497.2, 0.00021664, 0.0003699],
            [13901.7, 0.00020842, 0.0003614],
            [15464.1, 0.00019691, 0.00034675],
            [17202, 0.00018514, 0.00033107],
            [19135.3, 0.00017163, 0.00031167],
            [21285.9, 0.00015776, 0.00029093],
            [23678.1, 0.00014272, 0.00026726],
            [26339.3, 0.00012675, 0.00024101],
            [29299.5, 0.00011117, 0.00021463],
            [32592.3, 9.4749e-05, 0.00018574],
            [36255.3, 8.018e-05, 0.00015959],
            [40329.9, 6.4972e-05, 0.00013129],
            [44862.5, 5.2568e-05, 0.00010784],
            [49904.4, 4.016e-05, 8.3637e-05],
            [55513.1, 3.0596e-05, 6.4686e-05],
            [61752, 2.1867e-05, 4.6928e-05],
            [68692.1, 1.5406e-05, 3.3561e-05],
            [76412.2, 1.12e-05, 2.2543e-05],
            [85000, 8.2766e-06, 1.4649e-05],
        ]
    ),
    "blue": np.array(
        [
            [1200, 0.00017758, 0.00020576],
            [1334.86, 0.00018117, 0.00020992],
            [1484.89, 0.00018483, 0.00022494],
            [1651.77, 0.00018856, 0.00022948],
            [1837.41, 0.00020204, 0.0002459],
            [2043.91, 0.00020612, 0.00025085],
            [2273.61, 0.00021028, 0.00026881],
            [2529.14, 0.00021455, 0.00027423],
            [2813.38, 0.00021886, 0.00027976],
            [3129.57, 0.00022328, 0.00029978],
            [3481.29, 0.00022778, 0.00030582],
            [3872.55, 0.00023239, 0.00031199],
            [4307.77, 0.00023706, 0.00033433],
            [4791.91, 0.00024187, 0.00034105],
            [5330.46, 0.00024675, 0.00034795],
            [5929.53, 0.00025172, 0.00035496],
            [6595.93, 0.00025681, 0.00036214],
            [7337.23, 0.00024945, 0.00036945],
            [8161.84, 0.00025448, 0.00037691],
            [9079.13, 0.00024717, 0.00038452],
            [10099.5, 0.00025217, 0.00037348],
            [11234.6, 0.00024493, 0.00038103],
            [12497.2, 0.00023789, 0.0003701],
            [13901.7, 0.00023108, 0.00035947],
            [15464.1, 0.00022445, 0.00034917],
            [17202, 0.00020757, 0.00033915],
            [19135.3, 0.00019195, 0.00032943],
            [21285.9, 0.00018644, 0.00030464],
            [23678.1, 0.00016415, 0.00028172],
            [26339.3, 0.0001518, 0.00026052],
            [29299.5, 0.00013366, 0.00024093],
            [32592.3, 0.00011769, 0.00021212],
            [36255.3, 0.00010361, 0.00017782],
            [40329.9, 8.2696e-05, 0.00014906],
            [44862.5, 6.9321e-05, 0.00013124],
            [49904.4, 5.5324e-05, 0.00010474],
            [55513.1, 4.4155e-05, 7.9594e-05],
            [61752, 3.1944e-05, 6.0477e-05],
            [68692.1, 2.311e-05, 4.5957e-05],
            [76412.2, 1.5919e-05, 3.1654e-05],
            [85000, 1.0964e-05, 2.1804e-05],
        ]
    ),
    "green": np.array(
        [
            [1200, 0.0001075, 0.00015922],
            [1334.86, 0.00011068, 0.00018994],
            [1484.89, 0.00011968, 0.00021573],
            [1651.77, 0.00012322, 0.00024502],
            [1837.41, 0.00012685, 0.00029229],
            [2043.91, 0.0001306, 0.00031608],
            [2273.61, 0.00013446, 0.00035899],
            [2529.14, 0.00013843, 0.00040772],
            [2813.38, 0.0001357, 0.0004409],
            [3129.57, 0.0001397, 0.00047676],
            [3481.29, 0.00013694, 0.00051554],
            [3872.55, 0.00012779, 0.00055749],
            [4307.77, 0.00012526, 0.00057395],
            [4791.91, 0.00011691, 0.00059093],
            [5330.46, 0.00010387, 0.00060837],
            [5929.53, 9.694e-05, 0.00059635],
            [6595.93, 8.2007e-05, 0.00058455],
            [7337.23, 7.2865e-05, 0.00054551],
            [8161.84, 5.8687e-05, 0.00048469],
            [9079.13, 4.7271e-05, 0.00043068],
            [10099.5, 3.6247e-05, 0.00036433],
            [11234.6, 2.7795e-05, 0.0003082],
            [12497.2, 1.9321e-05, 0.00023633],
            [13901.7, 1.343e-05, 0.00017255],
            [15464.1, 1.0299e-05, 0.00011995],
            [17202, 7.8973e-06, 8.3376e-05],
            [19135.3, np.nan, 5.0017e-05],
            [21285.9, np.nan, 3.0004e-05],
            [23678.1, np.nan, 1.8e-05],
            [26339.3, np.nan, 1.0281e-05],
        ]
    ),
    "black": np.array(
        [
            [1200, 0.00013244, 0.00017104],
            [1334.86, 0.00014698, 0.00020347],
            [1484.89, 0.00016311, 0.00023378],
            [1651.77, 0.00017484, 0.0002686],
            [1837.41, 0.00018741, 0.00029809],
            [2043.91, 0.00019403, 0.00033081],
            [2273.61, 0.00020089, 0.00036712],
            [2529.14, 0.00020089, 0.0003801],
            [2813.38, 0.00020089, 0.00040744],
            [3129.57, 0.00019403, 0.00040744],
            [3481.29, 0.00018741, 0.00042182],
            [3872.55, 0.00016887, 0.00040744],
            [4307.77, 0.00015217, 0.0003801],
            [4791.91, 0.00012792, 0.0003425],
            [5330.46, 0.00010753, 0.00030862],
            [5929.53, 8.7312e-05, 0.00025944],
            [6595.93, 6.6137e-05, 0.00020347],
            [7337.23, 4.8388e-05, 0.00015957],
            [8161.84, 3.1901e-05, 0.00011276],
            [9079.13, 2.1031e-05, 7.9686e-05],
            [10099.5, 1.2494e-05, 4.901e-05],
            [11234.6, np.nan, 3.0144e-05],
            [12497.2, np.nan, 1.5584e-05],
            [13901.7, np.nan, 8.0579e-06],
        ]
    ),
}

PAPER_ROI = (-7.5, 7.5, -7.5, 7.5)  # Galactic lmin, lmax, bmin, bmax [deg] of Figure 13

# Figure style (matches the counts and flux spectra in Plots/).
LABEL_SIZE = 16
TICK_SIZE = 14
LEGEND_SIZE = 11
TITLE_SIZE = 14

FIG13_LABELS = {
    "black": "Pulsars, index-scaled",
    "green": "OB stars, index-scaled",
    "red": "Pulsars, intensity-scaled",
    "blue": "OB stars, intensity-scaled",
}
FIG13_HATCHES = {"black": "---", "green": "\\\\", "red": "|||", "blue": "///"}
FIG13_FLOOR = 7.1e-6  # lower plotted limit of the original figure [MeV cm^-2 s^-1]


def draw_fig13_envelopes(ax, paper_scale=1.0):
    """Draw the four Figure 13 envelopes, scaled by `paper_scale`.

    `paper_scale` = (template fraction in the comparison ROI) / (fraction in
    the paper's 15 x 15 deg ROI); 1 for the paper's ROI.
    """
    floor = FIG13_FLOOR * paper_scale
    for color in ("red", "blue", "green", "black"):
        energy, lo, hi = FIG13_ENVELOPES[color].T
        lo, hi = lo * paper_scale, hi * paper_scale
        # NaN = lower boundary was outside the paper's plotted range.
        # Show only its visible portion, without inventing a faint-energy tail.
        visible_lo = np.where(np.isfinite(lo), lo, floor)
        ax.fill_between(energy, visible_lo, hi, color=color, alpha=0.10)
        ax.fill_between(
            energy,
            visible_lo,
            hi,
            facecolor="none",
            edgecolor=color,
            hatch=FIG13_HATCHES[color],
            linewidth=1,
            label=FIG13_LABELS[color],
        )
        ax.plot(energy, lo, color=color, lw=0.8, alpha=0.6)
        ax.plot(energy, hi, color=color, lw=0.8, alpha=0.6)
    return floor


def style_gce_axes(ax, comparison_roi=PAPER_ROI):
    """Axis scales, labels, title and grid of the GCE comparison."""
    ax.set(xscale="log", yscale="log", xlim=(1e3, 1e5), ylim=(1e-5, 2e-3))
    ax.set_xlabel("Energy [MeV]", fontsize=LABEL_SIZE)
    ax.set_ylabel(r"$E^2\,dN/dE$ [MeV cm$^{-2}$ s$^{-1}$]", fontsize=LABEL_SIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_SIZE)
    lmin, lmax, bmin, bmax = comparison_roi
    ax.set_title(
        f"GCE: $l$ = [{lmin:g}°, {lmax:g}°], $b$ = [{bmin:g}°, {bmax:g}°]", fontsize=TITLE_SIZE
    )
    ax.grid(which="major", ls="--", color="grey", alpha=0.2)


def plot_gce_comparison(
    fit_result,
    nfw_path,
    comparison_roi=PAPER_ROI,
    output="GCE_Ajello2016_comparison.pdf",
    boundary_subpixels=4,
    max_samples=5000,
    data_output=None,
):
    """Plot the fitted GCE_NFW spectrum against Ajello et al. (2016) Figure 13.

    Parameters
    ----------
    fit_result : threeML analysis results
        Results of the completed fit, e.g. ``jl.results``. The saved best-fit
        model is used, so later changes to the live model have no effect.
    nfw_path : str or Path
        The unit-integral NFW template used in the fit (``DM_Model/nfw_gamma1.fits``).
    comparison_roi : tuple
        Galactic (lmin, lmax, bmin, bmax) in degrees; no plane mask. Default:
        the paper's 15 x 15 deg region.
    output : str or Path
        Output figure file.
    boundary_subpixels : int
        Supersampling per axis for pixels on the ROI boundary.
    max_samples : int
        Maximum number of correlated parameter samples for the 68% band.
    data_output : str or Path, optional
        Also write the plotted fit (energy, best fit and 68% band in the
        comparison ROI) to this CSV file, e.g. for overlays of several fits.

    Returns
    -------
    matplotlib.figure.Figure
    """

    # Read the SAVED best fit, so later changes to the live model do not affect it.
    _fitted = fit_result.optimized_model
    _shape = _fitted.GCE_NFW.spectrum.main.shape
    _spatial = _fitted.GCE_NFW.spatial_shape
    if _spatial.K.free or not np.isclose(_spatial.K.value, 1.0):
        raise ValueError("This cell assumes the spatial K is fixed at 1.")
    if _shape.piv.free:
        raise ValueError("This cell assumes a fixed spectral pivot.")

    # Integrate the actual FITS map. The native WCS may be Galactic or equatorial;
    # pixel_to_skycoord handles the conversion before testing Galactic ROI cuts.
    with fits.open(nfw_path) as _hdul:
        _hdu = next(h for h in _hdul if h.data is not None and h.data.ndim == 2)
        _map = np.asarray(_hdu.data, dtype=float).copy()
        _wcs = WCS(_hdu.header).celestial
    if not np.all(np.isfinite(_map)) or np.any(_map < 0):
        raise ValueError("Expected a finite, nonnegative NFW intensity map.")
    _ny, _nx = _map.shape
    _y, _x = np.indices(_map.shape, dtype=float)

    # Exact spherical solid angles of parent pixels, split into two triangles.
    def _unit_vectors(x, y):
        return np.moveaxis(pixel_to_skycoord(x, y, _wcs).cartesian.xyz.value, 0, -1)

    def _triangle_area(a, b, c):
        numerator = np.abs(np.sum(a * np.cross(b, c), axis=-1))
        denominator = 1 + np.sum(a * b, axis=-1) + np.sum(b * c, axis=-1) + np.sum(c * a, axis=-1)
        return 2 * np.arctan2(numerator, denominator)

    _a = _unit_vectors(_x - 0.5, _y - 0.5)
    _b = _unit_vectors(_x + 0.5, _y - 0.5)
    _c = _unit_vectors(_x + 0.5, _y + 0.5)
    _d = _unit_vectors(_x - 0.5, _y + 0.5)
    _omega = _triangle_area(_a, _b, _c) + _triangle_area(_a, _c, _d)
    _total = np.sum(_map * _omega)
    if not np.isfinite(_total) or _total <= 0:
        raise ValueError("Invalid template integral.")
    if not np.isclose(_total, 1.0, rtol=0.02):
        raise ValueError(
            f"FITS integral is {_total:g}, expected ~1 sr-weighted. "
            "Check that this is the unit-integral map used in the fit."
        )

    def _roi_fraction(roi):
        lmin, lmax, bmin, bmax = roi
        if not (-180 <= lmin < lmax <= 180 and -90 < bmin < bmax < 90):
            raise ValueError("Use increasing Galactic bounds in [-180,180].")
        # Ensure requested comparison region is contained in the map footprint.
        ll, bb = np.meshgrid(np.linspace(lmin, lmax, 121), np.linspace(bmin, bmax, 121))
        xx, yy = skycoord_to_pixel(SkyCoord(ll * u.deg, bb * u.deg, frame="galactic"), _wcs)
        if np.any(~np.isfinite(xx + yy)) or np.any(
            (xx < -0.5) | (xx > _nx - 0.5) | (yy < -0.5) | (yy > _ny - 0.5)
        ):
            raise ValueError("Requested ROI extends beyond the NFW template footprint.")
        occupancy = np.zeros_like(_map)
        offsets = (np.arange(boundary_subpixels) + 0.5) / boundary_subpixels - 0.5
        for dy in offsets:
            for dx in offsets:
                sky = pixel_to_skycoord(_x + dx, _y + dy, _wcs).galactic
                lon = sky.l.wrap_at(180 * u.deg).deg
                lat = sky.b.deg
                occupancy += (lon >= lmin) & (lon <= lmax) & (lat >= bmin) & (lat <= bmax)
        occupancy /= boundary_subpixels**2
        return np.sum(_map * _omega * occupancy) / _total

    _paper_roi = (-7.5, 7.5, -7.5, 7.5)
    _fpaper = _roi_fraction(_paper_roi)
    _ftarget = _fpaper if tuple(comparison_roi) == _paper_roi else _roi_fraction(comparison_roi)
    if min(_fpaper, _ftarget) <= 0:
        raise ValueError("Empty ROI integral.")
    _paper_scale = _ftarget / _fpaper
    print(f"Template flux fraction in paper ROI: {_fpaper:.5f}")
    print(f"Template flux fraction in comparison ROI: {_ftarget:.5f}")
    print(f"Multiply your full-template spectrum by {_ftarget:.5f}.")
    print(f"Multiply Figure 13 fluxes by {_paper_scale:.5f}.")

    # Correlated samples supplied by threeML: SAME row selection for every parameter.
    # For MLE these approximate the likelihood using its fitted covariance; for
    # Bayesian results they are posterior samples. No independent parameter draws.
    _pars = [_shape.K, _shape.index, _shape.xc]
    _columns = [
        np.asarray(fit_result.get_variates(p.path), dtype=float) if p.free else None for p in _pars
    ]
    _lengths = [len(a) for a in _columns if a is not None]
    if not _lengths or len(set(_lengths)) != 1:
        raise ValueError("Need aligned threeML samples for free GCE parameters.")
    _n = _lengths[0]
    _rows = np.random.default_rng(2026).choice(_n, min(_n, max_samples), replace=False)
    _draws = np.column_stack(
        [
            a[_rows] if a is not None else np.full(len(_rows), p.value)
            for a, p in zip(_columns, _pars)
        ]
    )
    _valid = np.all(np.isfinite(_draws), axis=1) & (_draws[:, 0] > 0) & (_draws[:, 2] > 0)
    if not np.all(_valid):
        warnings.warn(
            f"Discarded {np.mean(~_valid):.1%} nonphysical draws. "
            "If appreciable, use a profile likelihood or posterior band."
        )
    _draws = _draws[_valid]
    if len(_draws) < 100:
        raise ValueError("Too few usable uncertainty samples.")

    _E_MeV = np.geomspace(1e3, 1e5, 250)
    _E_keV = 1e3 * _E_MeV
    _pivot = _shape.piv.value
    _K, _index, _xc = _draws.T
    # Convert dN/dE per keV to per MeV by multiplying by 1000.
    # Thus E_MeV^2 * [1000*dN/dE_keV] gives MeV cm^-2 s^-1.
    _convert = _ftarget * 1e3 * _E_MeV**2
    _best = _convert * np.asarray(_shape(_E_keV), dtype=float)
    _expected = (
        _shape.K.value * (_E_keV / _pivot) ** _shape.index.value * np.exp(-_E_keV / _shape.xc.value)
    )
    if not np.allclose(np.asarray(_shape(_E_keV)), _expected, rtol=1e-6, atol=0):
        raise ValueError("Spectrum is not the assumed Cutoff_powerlaw convention.")
    _samples = (
        _convert[None, :]
        * _K[:, None]
        * (_E_keV[None, :] / _pivot) ** _index[:, None]
        * np.exp(-_E_keV[None, :] / _xc[:, None])
    )
    _low, _high = np.percentile(_samples, [16, 84], axis=0)

    # Both spectra are now intrinsic fluxes in exactly the same Galactic rectangle.
    # Assuming the same standard gamma=1 NFW-squared morphology, this is a pure
    # aperture conversion. It does not refit the data or reproduce a different IEM.
    fig, ax = plt.subplots(figsize=(10, 7), constrained_layout=True)
    # Put a border around the plot
    for spine in ax.spines.values():
        spine.set_visible(True)
    draw_fig13_envelopes(ax, _paper_scale)
    ax.fill_between(
        _E_MeV,
        _low,
        _high,
        color="darkorange",
        alpha=0.35,
        label="Your fit: pointwise 68% band",
        zorder=9,
    )
    ax.plot(_E_MeV, _best, color="darkorange", lw=2.8, label="Your GCE_NFW best fit", zorder=10)
    style_gce_axes(ax, comparison_roi)
    ax.legend(fontsize=LEGEND_SIZE, loc="best", frameon=False)
    if data_output is not None:
        np.savetxt(
            data_output,
            np.column_stack([_E_MeV, _best, _low, _high]),
            header=(
                "GCE_NFW fit in the comparison ROI (E^2 dN/dE, MeV cm^-2 s^-1)\n"
                f"comparison_roi_deg = {tuple(comparison_roi)}\n"
                f"paper_scale = {_paper_scale:.8g}\n"
                "energy_MeV best low68 high68"
            ),
        )
    print(
        "Paper envelopes are graphical estimates; no published central curves "
        "or parameter covariance were reconstructed."
    )
    print(
        "Your band is pointwise and includes K/index/xc correlations. For MLE "
        "it is a covariance approximation, not a profile-likelihood interval."
    )
    fig.savefig(output, bbox_inches="tight")
    return fig
