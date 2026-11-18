"""Extend GALPROP mapcubes in energy with per-pixel power laws.

With energy dispersion, Fermitools needs the model a little below emin and
above emax (true-energy margin). GALPROP cubes sometimes stop short of that.
``prepare`` therefore writes padded *copies* under ``<outdir>/padded_galprop/``
when needed. The original planes are copied unchanged; only new planes beyond
the original grid are added, each pixel following a power law fitted to its
last ``fit_points`` planes and anchored to the boundary value.

Settings come from ``galprop_extrapolation`` in config.yaml (see ``DEFAULTS``):

* ``high`` / ``low``  : allow padding above / below the cube (bool)
* ``fit_points``     : number of boundary planes used to fit the slope
* ``index_min/max``  : clip the fitted exponent to this range (dN/dE ~ E**index)
* ``index_shift``    : add this to every fitted exponent (sensitivity tests)
* ``max_factor``     : refuse to extrapolate more than this factor in energy

Padding decisions are recorded in ``<outdir>/galprop_extrapolation.json``.
"""

from pathlib import Path
import numpy as np

# Default padding settings; overridden by the galprop_extrapolation block.
DEFAULTS = dict(
    high=True,
    low=False,
    fit_points=3,
    index_min=-6.0,
    index_max=2.0,
    index_shift=0.0,
    max_factor=3.0,
)


def options(spec=None):
    """Merge user settings with ``DEFAULTS`` and validate them."""
    cfg = {**DEFAULTS, **(spec or {})}
    if set(cfg) - set(DEFAULTS):
        raise ValueError("Unknown galprop_extrapolation option")
    if not all(isinstance(cfg[k], bool) for k in ("high", "low")):
        raise ValueError("high and low must be YAML booleans")
    if not isinstance(cfg["fit_points"], int) or cfg["fit_points"] < 2:
        raise ValueError("fit_points must be an integer >= 2")
    if not all(
        np.isfinite(cfg[k]) for k in ("index_min", "index_max", "index_shift", "max_factor")
    ):
        raise ValueError("Extrapolation settings must be finite")
    if cfg["index_min"] >= cfg["index_max"] or cfg["max_factor"] < 1:
        raise ValueError("Invalid index bounds or max_factor")
    return cfg


def extend_arrays(energies, cube, lower, upper, spec=None):
    """Fit log I versus log E per pixel; anchor exactly at the boundary.

    Zero boundary pixels remain zero. A positive boundary with fewer than two
    positive fitting samples raises rather than inventing a spectral slope.
    Returned index is the exponent (negative for a falling spectrum).
    """
    cfg = options(spec)
    e = np.asarray(energies, dtype=float)
    a = np.asarray(cube)
    if (
        e.ndim != 1
        or len(e) < 2
        or not np.all(np.isfinite(e))
        or np.any(e <= 0)
        or np.any(np.diff(e) <= 0)
    ):
        raise ValueError("GALPROP energies must be finite, positive, strictly increasing")
    if a.ndim != 3 or a.shape[0] != len(e) or not np.all(np.isfinite(a)) or np.any(a < 0):
        raise ValueError("GALPROP cube must be finite, nonnegative, and match ENERGIES")
    if not (0 < lower < upper):
        raise ValueError("Invalid requested energy bounds")
    report = dict(original_range_MeV=[float(e[0]), float(e[-1])], settings=cfg, sides={})
    energies_out, arrays_out = [], []
    for side, needed, edge, target in [
        ("low", lower < e[0], 0, lower),
        ("high", upper > e[-1], -1, upper),
    ]:
        if not needed:
            continue
        factor = e[0] / lower if side == "low" else upper / e[-1]
        if not cfg[side] or factor > cfg["max_factor"]:
            raise ValueError(
                f"{side} padding required: actual {e[0]:g}–{e[-1]:g} MeV; "
                f"required {lower:g}–{upper:g} MeV; factor {factor:g}. "
                "Enable this side and set max_factor explicitly if appropriate."
            )
        n = min(cfg["fit_points"], len(e))
        sub = slice(0, n) if side == "low" else slice(-n, None)
        x = np.log(e[sub] / e[edge])[:, None, None]
        y = np.asarray(a[sub], float)
        valid = y > 0
        logy = np.log(np.where(valid, y, 1.0))
        count = valid.sum(axis=0)
        sx = (x * valid).sum(axis=0)
        sy = (logy * valid).sum(axis=0)
        den = count * (x * x * valid).sum(axis=0) - sx * sx
        active = a[edge] > 0
        good = (count >= 2) & (den > 0)
        if np.any(active & ~good):
            raise ValueError(f"{side}: positive boundary pixels lack two positive slope samples")
        slope = np.zeros(a.shape[1:], float)
        np.divide(count * (x * logy * valid).sum(axis=0) - sx * sy, den, out=slope, where=good)
        slope += cfg["index_shift"]
        clipped = active & ((slope < cfg["index_min"]) | (slope > cfg["index_max"]))
        slope = np.clip(slope, cfg["index_min"], cfg["index_max"])
        step = np.min(np.diff(np.log(e[sub])))
        nadd = max(1, int(np.ceil(np.log(factor) / step)))
        if side == "low":
            extra = np.geomspace(target, e[0], nadd + 1)[:-1]
        else:
            extra = np.geomspace(e[-1], target, nadd + 1)[1:]
        with np.errstate(over="raise", invalid="raise"):
            vals = np.asarray(a[edge], float)[None, :, :] * np.exp(
                np.log(extra / e[edge])[:, None, None] * slope[None, :, :]
            )
        info = dict(
            added_planes=nadd,
            zero_boundary_pixels=int((~active).sum()),
            clipped_pixels=int(clipped.sum()),
            index_percentiles=np.percentile(slope[active], [0, 16, 50, 84, 100]).tolist()
            if np.any(active)
            else [],
        )
        report["sides"][side] = info
        if side == "low":
            energies_out.append(extra)
            arrays_out.append(vals)
        else:
            high = (extra, vals)
    energies_out.append(e)
    arrays_out.append(a)
    if "high" in report["sides"]:
        energies_out.append(high[0])
        arrays_out.append(high[1])
    return np.concatenate(energies_out), np.concatenate(arrays_out), report


def pad_mapcube(source, destination, lower, upper, spec=None):
    """Pad a FITS mapcube so it covers [lower, upper] MeV.

    Returns ``(path_to_use, report)``. If no padding is needed the original
    path is returned unchanged; otherwise the padded copy is written to
    ``destination`` (the original file is never modified).
    """
    from astropy.io import fits

    if Path(source).resolve() == Path(destination).resolve():
        raise ValueError("Padding output must differ from the original file")
    with fits.open(source, memmap=False) as hdus:
        table = hdus["ENERGIES"]
        unit = table.columns[0].unit
        if unit and unit.strip().lower() != "mev":
            raise ValueError(f"ENERGIES must be in MeV, got {unit!r}: {source}")
        e = np.asarray(table.data.field(0), float).ravel()
        ee, aa, report = extend_arrays(e, hdus[0].data, lower, upper, spec)
        report["source"] = str(source)
        if not report["sides"]:
            report["output"] = str(source)
            return str(source), report
        # Preserve spatial WCS and other extensions. ENERGIES is authoritative;
        # remove stale energy-axis WCS because the extended grid may be irregular.
        header = hdus[0].header.copy()
        for key in ("CRVAL3", "CRPIX3", "CDELT3", "CTYPE3", "CUNIT3"):
            header.pop(key, None)
        header.add_history("Boundary PL padding; original energy planes preserved.")
        hdus[0] = fits.PrimaryHDU(aa, header=header)
        hdus[0].header["EXTEND"] = True
        idx = hdus.index_of("ENERGIES")
        hdus[idx] = fits.BinTableHDU.from_columns(
            [fits.Column(name=table.columns[0].name, format="D", unit="MeV", array=ee)],
            name="ENERGIES",
        )
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        hdus.writeto(destination, overwrite=True, checksum=True)
        report["output"] = str(destination)
    return str(destination), report
