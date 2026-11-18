"""Spectrum/count-spectrum plots and unsmoothed masked Mollweide residual maps.

Uses actual saved pixel coordinates; no RING/NESTED assumption or healpy import.
The display grid is a nearest-pixel rendering, not a rebinned counts dataset.
"""

from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree


def residuals(n, mu):
    """Return (data - model, (data - model)/model, (data - model)/sqrt(model)).

    Pixels with model = 0 get NaN for the two normalized residuals.
    """
    delta = n - mu
    fractional = np.divide(delta, mu, out=np.full_like(mu, np.nan, dtype=float), where=mu > 0)
    pearson = np.divide(delta, np.sqrt(mu), out=np.full_like(mu, np.nan, dtype=float), where=mu > 0)
    return delta, fractional, pearson


def vectors(l, b):
    """Unit vectors for Galactic (l, b) in degrees (used for nearest-pixel lookup)."""
    l, b = np.radians(l), np.radians(b)
    return np.column_stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])


def plot_results(directory, options=None):
    """Write all plots for one fit result folder.

    Outputs (in ``directory``):

    * ``isotropic_spectrum.png/.pdf`` : E^2 x isotropic intensity per bin.
      With profiles: points with 68% errors (TS >= 4) or 95% upper limits.
      Without profiles: best-fit points only.
    * ``counts_spectrum.png`` : summed counts over retained pixels, per
      component, with fractional residuals.
    * ``residuals_<band>.png`` : Mollweide maps of data - model, fractional and
      Pearson residuals for each band in ``options['bands_MeV']`` (whole energy
      bins whose centres fall in the band). Excluded pixels are gray.
    * ``plot_bands.txt`` : requested versus actual energy range of each band.

    ``options`` is the ``plots`` block of config.yaml: ``bands_MeV`` and
    ``map_nlon`` (number of longitude cells in the display grid).
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    options = options or {}
    directory = Path(directory)
    table = np.genfromtxt(directory / "isotropic_bins.csv", delimiter=",", names=True)
    table = np.atleast_1d(table)
    e = table["energy_MeV"]
    y = table["intensity"]
    factor = e**2 / 1000.0
    det = (table["iso_TS"] >= 4) & (y > 0)
    fig, ax = plt.subplots(figsize=(7, 4.8), layout="constrained")
    if np.any(det):
        ax.errorbar(
            e[det] / 1000.0,
            y[det] * factor[det],
            xerr=np.array([e[det] - table["emin_MeV"][det], table["emax_MeV"][det] - e[det]])
            / 1000.0,
            yerr=np.maximum(
                0,
                np.array(
                    [y[det] - table["intensity_lo"][det], table["intensity_hi"][det] - y[det]]
                ),
            )
            * factor[det],
            fmt="o",
            color="#155e75",
            label="Isotropic (68% profile)",
        )
    limits = (~det) & np.isfinite(table["intensity_ul95"]) & (table["intensity_ul95"] > 0)
    if np.any(limits):
        lim = table["intensity_ul95"][limits] * factor[limits]
        ax.errorbar(
            e[limits] / 1000.0,
            lim,
            yerr=0.25 * lim,
            uplims=True,
            fmt="v",
            color="#9a3412",
            label="95% upper limit",
        )
    unprofiled = ~np.isfinite(table["iso_TS"]) & (y > 0)
    if np.any(unprofiled):
        ax.plot(
            e[unprofiled] / 1000.0,
            y[unprofiled] * factor[unprofiled],
            "o-",
            color="#155e75",
            label="Best fit (uncertainties not computed)",
        )
    ax.set(
        xscale="log",
        yscale="log",
        xlabel="Energy [GeV]",
        ylabel=r"$E^2 I(E)$ [GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]",
    )
    ax.set_title(
        "Effective isotropic spectrum\nEnergy dispersion disabled for this component", fontsize=11
    )
    ax.grid(alpha=0.2, which="both")
    ax.legend(fontsize=8)
    fig.savefig(directory / "isotropic_spectrum.png", dpi=180)
    fig.savefig(directory / "isotropic_spectrum.pdf")
    plt.close(fig)
    with np.load(directory / "diagnostics.npz", allow_pickle=False) as z:
        keep = z["keep"]
        edges = z["edges"]
        n = z["counts"]
        mu = z["model"]
        e = np.sqrt(edges[:-1] * edges[1:]) / 1000.0
        data = n[:, keep].sum(axis=1)
        model = mu[:, keep].sum(axis=1)
        fig, (ax, ar) = plt.subplots(
            2,
            1,
            figsize=(7, 6),
            sharex=True,
            layout="constrained",
            gridspec_kw={"height_ratios": [3, 1]},
        )
        ax.errorbar(
            e, data, yerr=np.sqrt(data), fmt=".", color="black", label="Data (retained sky)"
        )
        ax.plot(e, model, label="Total model", color="#0f766e")
        ax.plot(e, z["fixed_catalog"][:, keep].sum(axis=1), "--", label="Fixed 4FGL")
        ax.plot(e, z["isotropic"][:, keep].sum(axis=1), "--", label="Effective isotropic")
        for name, m in zip(z["galactic_names"], z["galactic"]):
            ax.plot(e, m[:, keep].sum(axis=1), ":", label=str(name))
        ax.set(yscale="log", ylabel="Counts per energy bin")
        ax.legend(fontsize=8)
        ax.set_ylim(1e0, 1e7)
        delta, frac, _ = residuals(data, model)
        ar.errorbar(
            e,
            frac,
            yerr=np.divide(np.sqrt(data), model, out=np.zeros_like(model), where=model > 0),
            fmt=".",
        )
        ar.axhline(0, color="black", lw=0.8)
        ar.set(xscale="log", xlabel="Energy [GeV]", ylabel="(D−M)/M")
        for a in [ax, ar]:
            a.grid(alpha=0.2)
        fig.savefig(directory / "counts_spectrum.png", dpi=180)
        plt.close(fig)
        # Nearest original pixel on a regular lon/lat display grid. No smoothing.
        nx = int(options.get("map_nlon", 360))
        ny = max(30, nx // 2)
        xx = np.linspace(-np.pi, np.pi, nx + 1)
        yy = np.linspace(-np.pi / 2, np.pi / 2, ny + 1)
        x, y = np.meshgrid((xx[:-1] + xx[1:]) / 2, (yy[:-1] + yy[1:]) / 2)
        tree = cKDTree(vectors(z["l"], z["b"]))
        _, index = tree.query(vectors(-np.degrees(x.ravel()), np.degrees(y.ravel())))
        bands = options.get("bands_MeV", [[float(edges[0]), float(edges[-1])]])
        centers = np.sqrt(edges[:-1] * edges[1:])
        log = []
        for low, high in bands:
            bins = (centers >= low) & (centers < high)
            if not np.any(bins):
                continue
            chosen = np.flatnonzero(bins)
            actual_low, actual_high = edges[chosen[0]], edges[chosen[-1] + 1]
            dn = n[bins].sum(axis=0)
            dm = mu[bins].sum(axis=0)
            maps = residuals(dn, dm)
            title = f"{actual_low / 1000:.3g}–{actual_high / 1000:.3g} GeV (whole energy bins; unsmoothed)"
            fig = plt.figure(figsize=(10, 10), layout="constrained")
            for j, (values, label) in enumerate(
                zip(
                    maps,
                    [
                        "Data − model [counts]",
                        "(Data − model) / model",
                        "(Data − model) / sqrt(model)",
                    ],
                )
            ):
                values = values.copy()
                values[~keep] = np.nan
                grid = values[index].reshape(ny, nx)
                finite = np.abs(values[np.isfinite(values)])
                limit = max(float(np.percentile(finite, 99)) if len(finite) else 1.0, 1e-12)
                ax = fig.add_subplot(3, 1, j + 1, projection="mollweide")
                ax.set_facecolor("#d1d5db")
                im = ax.pcolormesh(
                    xx,
                    yy,
                    grid,
                    cmap="RdBu_r",
                    vmin=-limit,
                    vmax=limit,
                    shading="flat",
                    rasterized=True,
                )
                ax.grid(alpha=0.25)
                ax.set_title(label, fontsize=10)
                ticks = np.arange(-150, 151, 30)
                ax.set_xticks(np.radians(ticks))
                ax.set_xticklabels([f"{-v}°" for v in ticks], fontsize=8)
                fig.colorbar(im, ax=ax, shrink=0.75, pad=0.03, extend="both")
            fig.suptitle(
                title
                + "\nGalactic coordinates; gray = excluded; Pearson residual is not a TS significance",
                fontsize=11,
            )
            tag = f"{actual_low:.0f}_{actual_high:.0f}MeV"
            fig.savefig(directory / f"residuals_{tag}.png", dpi=160)
            plt.close(fig)
            log.append(
                f"{low:g}–{high:g} requested MeV -> {actual_low:g}–{actual_high:g} actual MeV"
            )
        (directory / "plot_bands.txt").write_text("\n".join(log) + "\n")
    print("Plots written to", directory)
