#!/usr/bin/env python3
"""All-sky Fermi-LAT fit for the effective isotropic background spectrum.

Usage (run from this directory, inside a Fermipy/Fermitools environment):

    python fit_isotropic.py prepare config.yaml   # build LAT count templates (slow, once)
    python fit_isotropic.py fit     config.yaml   # masked Poisson fit of the templates (fast)
    python fit_isotropic.py plot    config.yaml   # remake plots from saved fit results
    python fit_isotropic.py run     config.yaml   # prepare + fit + plot in one go
    python fit_isotropic.py self-test             # synthetic checks; no LAT data needed

How it works
------------
1. ``prepare`` asks Fermipy for a full-sky HEALPix counts cube and, for every
   model component, its predicted-count cube ("template") at unit
   normalization. Templates include the LAT exposure, PSF and (except for the
   isotropic component) energy dispersion. Everything is saved to
   ``<outdir>/templates.npz``.
2. ``fit`` reads ``templates.npz``, applies the sky mask, and maximizes the
   binned Poisson likelihood for the linear model

       mu[k, p] = fixed[k, p] + a_k * T_iso[k, p] + sum_j A_j * T_j[k, p]

   where k is the energy bin and p the pixel. ``fixed`` is the sum of all
   catalog sources (held at 4FGL values), a_k is one isotropic amplitude per
   energy bin, and A_j is one energy-independent amplitude per free diffuse
   component (GALPROP pi0/IC/brem, Sun, Moon, Loop I, ...). Diffuse components
   marked ``free: false`` are added to ``fixed`` with their configured value.
3. The isotropic intensity in bin k is ``a_k * I_ref(E_k)``, where I_ref is the
   reference spectrum used to build the template and E_k the bin's geometric
   centre. Results are written to a ``bestfit_<hash>/`` (or ``fit_<hash>/``
   with profiles) folder inside ``outdir``.

Energy dispersion is ON for catalog, GALPROP and additional diffuse components
and OFF for the isotropic component, matching the standard LAT convention for
isotropic templates. Keep it OFF when you reuse ``isotropic_spectrum.txt``.

``self-test`` and fixed-shape ``fit``/``plot`` do not import Fermitools. An
analytic Loop I spectral-shape fit (see ``loopI.py``) does require Fermitools.
"""

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from scipy.optimize import minimize, brentq
from scipy.special import xlogy
import yaml
from joint_fit import JointPoisson
from galprop_padding import pad_mapcube, options as padding_options

# Version of the preparation/cache format. It is part of the preparation
# signature: changing it invalidates every existing templates.npz.
PIPELINE_VERSION = 6


# ---------------------------------------------------------------------------
# Small helpers: configuration lookups, isotropic reference spectrum, sky mask
# ---------------------------------------------------------------------------


def diffuse_components(c):
    """Return all non-isotropic diffuse components (GALPROP first, then extras).

    The order of this mapping defines the order of the diffuse templates in
    ``templates.npz`` (after the isotropic template, which is always first).
    """
    return {**c.get("galprop", {}), **c.get("additional_diffuse", {})}


def galprop_controls(c, names):
    """Return (free flags, normalizations) for the named diffuse components.

    Parameters
    ----------
    c : dict
        Configuration from :func:`read_config`.
    names : sequence of str
        Diffuse component names, in template order (``templates.npz`` order).

    Returns
    -------
    free : ndarray of bool
        True if the component's normalization is fitted.
    values : ndarray of float
        Fixed multiplier (``free: false``) or starting value (``free: true``).
    """
    free = []
    values = []
    for name in names:
        item = diffuse_components(c).get(str(name), {})
        item = item if isinstance(item, dict) else {}
        free.append(item.get("free", c.get("galprop_mode", "global") != "fixed"))
        values.append(item.get("normalization", 1.0))
    return np.array(free, bool), np.array(values, float)


def build_seed(c):
    """Build the isotropic reference spectrum ("seed") as an (E, I) table.

    The seed only sets the spectral shape *within* each energy bin; its overall
    normalization is arbitrary because every bin gets its own free amplitude.

    * ``isotropic.mode: powerlaw`` -> ``I(E) = intensity_at_pivot * (E/pivot_MeV)**index``
      tabulated on 256 log-spaced energies from emin/2 to 2*emax (wide enough
      for the energy-dispersion margins).
    * ``isotropic.mode: file`` -> two-column text file (MeV, ph cm^-2 s^-1 sr^-1 MeV^-1).

    Returns
    -------
    ndarray, shape (N, 2)
        Columns: energy [MeV], differential intensity [ph cm^-2 s^-1 sr^-1 MeV^-1].
    """
    spec = c.get("isotropic", {"mode": "file", "file": c.get("isotropic_seed")})
    if spec["mode"] == "file":
        table = np.loadtxt(spec["file"], ndmin=2)[:, :2]
        return extend_seed(table, c["emin"], c["emax"], c.get("allow_seed_extrapolation", False))
    require(spec["mode"] == "powerlaw", "isotropic.mode must be powerlaw or file.")
    index = float(spec.get("index", -2.3))
    pivot = float(spec.get("pivot_MeV", 1000.0))
    amplitude = float(spec.get("intensity_at_pivot", 1e-10))
    require(
        np.isfinite(index)
        and np.isfinite(pivot)
        and pivot > 0
        and np.isfinite(amplitude)
        and amplitude > 0,
        "Invalid analytic isotropic reference.",
    )
    energies = np.geomspace(c["emin"] / 2, c["emax"] * 2, 256)
    return np.column_stack([energies, amplitude * (energies / pivot) ** index])


def require(condition, message):
    """Raise ``ValueError(message)`` unless ``condition`` is true."""
    if not condition:
        raise ValueError(message)


def sky_mask(l, b, dec, spec):
    """Return True for pixels that are KEPT in the likelihood.

    A pixel (judged at its centre) is kept when all of these hold:

    * ``|b| >= plane_abs_b_min``            (removes the Galactic plane)
    * ``|dec| <= pole_abs_dec_max``         (optional celestial-pole cut; 90 = off)
    * NOT (``|l| <= bubbles_abs_l_max`` and ``|b| <= bubbles_abs_b_max``)
      (removes a rectangle around the Fermi bubbles / inner Galaxy)

    Longitudes are wrapped to [-180, 180) first, so ``|l|`` is measured from
    the Galactic centre.
    """
    l = (np.asarray(l) + 180.0) % 360.0 - 180.0
    b, dec = np.asarray(b), np.asarray(dec)
    return (
        (np.abs(b) >= spec["plane_abs_b_min"])
        & (np.abs(dec) <= spec["pole_abs_dec_max"])
        & ~((np.abs(l) <= spec["bubbles_abs_l_max"]) & (np.abs(b) <= spec["bubbles_abs_b_max"]))
    )


def seed_at(energy, table):
    """Evaluate a tabulated spectrum at ``energy`` by log-log interpolation.

    Raises if ``energy`` lies outside the table: no silent extrapolation.
    """
    require(
        np.all(np.isfinite(table[:, :2])) and np.all(table[:, :2] > 0),
        "Seed energies/intensities must be finite and positive.",
    )
    require(np.all(np.diff(table[:, 0]) > 0), "Seed energies must increase.")
    require(
        np.min(energy) >= table[0, 0] and np.max(energy) <= table[-1, 0],
        "Seed does not cover requested energies; no extrapolation allowed.",
    )
    return np.exp(np.interp(np.log(energy), np.log(table[:, 0]), np.log(table[:, 1])))


def extend_seed(table, emin, emax, allow=False):
    """Extend a file-based seed to cover [emin, emax] with endpoint power laws.

    Only used when ``allow_seed_extrapolation: true``. The added points use the
    slope of the two nearest table entries. They only shape the reference
    spectrum inside a bin; they are never treated as data.
    """
    table = np.asarray(table, float)
    seed_at(table[:, 0], table)
    if len(table) < 2:
        raise ValueError("Need at least two isotropic seed points.")
    result = table.copy()
    for energy, i, j in [(emin, 0, 1), (emax, -2, -1)]:
        if table[0, 0] <= energy <= table[-1, 0]:
            continue
        require(
            allow,
            "Seed does not cover requested range. Supply a wider seed or explicitly "
            "set allow_seed_extrapolation: true to use an endpoint power-law initial shape.",
        )
        slope = np.log(table[j, 1] / table[i, 1]) / np.log(table[j, 0] / table[i, 0])
        value = table[i, 1] * (energy / table[i, 0]) ** slope
        result = np.vstack([result, [energy, value]])
    return result[np.argsort(result[:, 0])]


def choose_edisp(emin, emax, binsperdec, setting):
    """Validate (or choose, for ``'auto'``) the Fermitools ``edisp_bins`` value.

    ``edisp_bins`` is the number of extra true-energy bins added on each side
    of the measured energy range for energy dispersion. Negative values use the
    standard Fermitools method in which dispersion is applied to the source
    spectrum (recommended for E >= 100 MeV, e.g. -2 as in 4FGL). Positive
    values use the spatial-map mixing method, required below 100 MeV.

    The extra bins must span at least 0.2 decade above 100 MeV (0.4 decade
    below), enough to cover LAT energy migration.

    ``'auto'`` returns -2 for the default 300 MeV-100 GeV, 8 bins/decade setup,
    and +4 when emin < 100 MeV.
    """
    ne = int(np.round(binsperdec * np.log10(emax / emin)))
    require(ne > 0, "No energy bins.")
    width = np.log10(emax / emin) / ne
    span = 0.4 if emin < 100 else 0.2
    needed = max(2, int(np.ceil(span / width - 1e-10)))
    if setting == "auto":
        setting = needed if emin < 100 or needed > 3 else -needed
    require(setting in [-3, -2, 2, 3, 4], "Use auto, -2/-3, or 2/3/4 for edisp_bins.")
    require(
        abs(setting) * width >= span - 1e-10,
        "Energy migration range too narrow. Increase edisp_bins up to 4 or use coarser energy bins.",
    )
    require(
        emin >= 100 or setting > 0, "Below 100 MeV use positive edisp_bins (spatial-map mixing)."
    )
    return setting


# ---------------------------------------------------------------------------
# Single-energy Poisson solver. The production fit uses JointPoisson
# (joint_fit.py); these simpler functions are kept for the self-tests.
# ---------------------------------------------------------------------------


def nll(n, mu):
    """Poisson negative log likelihood ``sum(mu - n ln mu)`` (drops ln n!).

    Returns ``inf`` for an impossible model (negative, or zero where n > 0).
    """
    if np.any(mu < 0) or np.any((mu == 0) & (n > 0)):
        return np.inf
    return float(np.sum(mu - xlogy(n, mu)))


def solve(n, fixed, templates, start=None):
    """Fit ``mu = fixed + a @ templates`` to counts ``n`` with ``a >= 0``.

    Parameters
    ----------
    n, fixed : ndarray, shape (npix,)
        Observed counts and fixed model counts.
    templates : ndarray, shape (ncomp, npix)
        Predicted counts of each fitted component at unit amplitude.
    start : array-like, optional
        Starting amplitudes (default: all ones).

    Returns
    -------
    a : ndarray
        Best-fit amplitudes.
    value : float
        Minimum negative log likelihood.
    cov : ndarray
        Inverse Hessian; all NaN if the Hessian is (near-)singular.
    """
    n, fixed, t = map(lambda x: np.asarray(x, dtype=float), (n, fixed, templates))
    require(
        n.ndim == 1 and fixed.shape == n.shape and t.ndim == 2 and t.shape[1] == n.size,
        "Inconsistent count/template shapes.",
    )
    require(
        all(np.all(np.isfinite(x)) and np.all(x >= 0) for x in (n, fixed, t)),
        "Counts and templates must be finite and nonnegative.",
    )
    support = fixed + t.sum(axis=0) > 0
    require(not np.any((~support) & (n > 0)), "Observed counts with no model support.")
    n, fixed, t = n[support], fixed[support], t[:, support]
    require(n.size > 0, "No supported pixels.")
    if t.shape[0] == 0:
        return np.empty(0), nll(n, fixed), np.empty((0, 0))
    require(np.all(t.sum(axis=1) > 0), "A fitted template has no support.")
    # Normalize the objective by total counts for stable stopping tolerances.
    scale = max(float(n.sum()), 1.0)

    def fg(a):
        mu = fixed + a @ t
        if np.any((mu <= 0) & (n > 0)):
            return np.inf, np.full(a.shape, -1e100)
        ratio = np.divide(n, mu, out=np.zeros_like(n), where=mu > 0)
        return nll(n, mu) / scale, (t @ (1.0 - ratio)) / scale

    starts = [np.ones(t.shape[0]) if start is None else np.asarray(start), np.full(t.shape[0], 0.3)]
    candidates = []
    for x in starts:
        r = minimize(
            fg,
            x,
            jac=True,
            bounds=[(0.0, None)] * len(x),
            method="L-BFGS-B",
            options={"ftol": 1e-14, "gtol": 1e-10, "maxiter": 3000, "maxls": 60},
        )
        g = fg(r.x)[1]
        projected = np.where((r.x < 1e-8) & (g > 0), 0.0, g)
        if r.success and np.isfinite(r.fun) and np.max(np.abs(projected)) < 2e-6:
            candidates.append(r)
    require(bool(candidates), "Optimizer failed convergence/KKT checks.")
    r = min(candidates, key=lambda r: r.fun)
    mu = fixed + r.x @ t
    w = np.divide(n, mu**2, out=np.zeros_like(n), where=mu > 0)
    h = (t * w) @ t.T
    # Never report a pseudo-inverse as a well-measured covariance.
    cov = np.full(h.shape, np.nan)
    if np.linalg.cond(h) < 1e12:
        cov = np.linalg.inv(h)
    return r.x, nll(n, mu), cov


def profile_iso(n, fixed, templates, best, minimum):
    """Profile-likelihood interval for the first amplitude (the isotropic one).

    The other amplitudes are re-fitted at every trial value.

    Returns
    -------
    lo, hi : float
        68% interval (delta NLL = 0.5); ``lo`` is 0 when the interval touches 0.
    ul95 : float
        One-sided 95% upper limit (delta NLL = 2.71/2 = 1.3528).
    ts : float
        Test statistic ``2 * [NLL(a=0) - NLL(best)]``.
    """

    def delta(a):
        support = fixed + a * templates[0] + templates[1:].sum(axis=0)
        if np.any((support == 0) & (n > 0)):
            return np.inf
        if not np.any(support):
            return -minimum
        _, f, _ = solve(n, fixed + a * templates[0], templates[1:], best[1:])
        return f - minimum

    a = best[0]
    d0 = delta(0.0)
    require(d0 >= -1e-3, "Profile found a lower likelihood minimum.")

    def upper(level):
        hi = max(1.0, 2.0 * a + 0.1)
        for _ in range(60):
            if delta(hi) >= level:
                return brentq(lambda x: delta(x) - level, a, hi, xtol=1e-7)
            hi *= 2.0
        raise ValueError("Could not bracket profile upper limit.")

    lo = 0.0 if d0 <= 0.5 or a < 1e-10 else brentq(lambda x: delta(x) - 0.5, 0.0, a)
    return lo, upper(0.5), upper(1.352771727), max(0.0, 2.0 * d0)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------


class UniqueKeyLoader(yaml.SafeLoader):
    """YAML loader that rejects duplicate keys.

    Plain ``yaml.safe_load`` silently keeps the *last* of two identical keys
    (e.g. two ``outdir:`` lines), which can send output somewhere unexpected.
    This loader raises an error with the line number instead. YAML merge keys
    (``<<: *anchor``) are still allowed.
    """

    def construct_mapping(self, node, deep=False):
        seen = set()
        for key_node, _ in node.value:
            if key_node.tag == "tag:yaml.org,2002:merge":
                continue
            key = self.construct_object(key_node, deep=deep)
            if key in seen:
                raise yaml.constructor.ConstructorError(
                    "while reading configuration",
                    node.start_mark,
                    f"duplicate key: {key}",
                    key_node.start_mark,
                )
            seen.add(key)
        return super().construct_mapping(node, deep=deep)


def read_config(path):
    """Read, validate and complete the wrapper configuration (``config.yaml``).

    What this does:

    * Resolves every path. Relative paths are relative to the YAML file's
      directory (not the shell's working directory); ``~`` and ``$VARS`` work.
    * Fills defaults: each GALPROP/additional component gets ``free`` and
      ``normalization`` if missing.
    * Expands the ``loopI`` block into ``additional_diffuse`` entries
      (``LoopI`` or ``LoopI_S1``/``LoopI_S2``), see :func:`loopI.configure`.
    * Resolves ``edisp_bins: auto`` to an integer (see :func:`choose_edisp`).
    * Checks values (energy range, time range, mask limits, names, ...).

    Returns
    -------
    dict
        The completed configuration used by ``prepare``, ``fit`` and ``plot``.
    """
    path = Path(path).resolve()
    c = yaml.load(path.read_text(), Loader=UniqueKeyLoader)
    # Old configs may still contain this switch; reuse is now always selective
    # (source by source) whenever reuse_from is set.
    c.pop("reuse_source_maps", None)
    # Retired switch: GALPROP units are no longer confirmed in the config.
    # Old configs that still contain it are accepted.
    c.pop("galprop_intensity_units_confirmed", None)
    # --- Input/output paths -------------------------------------------------
    for key in ["evfile", "scfile", "catalog", "extdir", "outdir"]:
        p = Path(os.path.expandvars(c[key])).expanduser()
        c[key] = str((path.parent / p).resolve() if not p.is_absolute() else p)
    if c.get("reuse_from"):
        c["reuse_from"] = str(
            (path.parent / Path(os.path.expandvars(c["reuse_from"])).expanduser()).resolve()
        )
    # --- GALPROP components (a bare string is accepted as just the file) ---
    for name, item in c["galprop"].items():
        if isinstance(item, str):
            item = {"file": item}
        item = dict(item)
        item["file"] = str(
            (path.parent / Path(os.path.expandvars(item["file"])).expanduser()).resolve()
        )
        item.setdefault("free", c.get("galprop_mode", "global") != "fixed")
        item.setdefault("normalization", 1.0)
        require(isinstance(item["free"], bool), "GALPROP free must be true/false.")
        require(
            np.isfinite(item["normalization"]) and item["normalization"] >= 0,
            "GALPROP normalization must be finite and nonnegative.",
        )
        c["galprop"][name] = item
    # --- Loop I: adds LoopI / LoopI_S1 / LoopI_S2 to additional_diffuse ---
    if c.get("loopI", {}).get("directory") is not None:
        p = Path(os.path.expandvars(c["loopI"]["directory"])).expanduser()
        c["loopI"]["directory"] = str((path.parent / p).resolve())
    from loopI import configure

    configure(c)
    # --- Additional diffuse components (Sun, Moon, Loop I, user cubes) ------
    extra = c.setdefault("additional_diffuse", {})
    require(isinstance(extra, dict), "additional_diffuse must be a mapping.")
    require(
        not (set(extra) & (set(c["galprop"]) | {"isodiff"})),
        "Additional source names must be unique and must not be isodiff.",
    )
    for name, item in extra.items():
        item = {"file": item} if isinstance(item, str) else dict(item)
        require(
            not (set(item) - {"file", "free", "normalization"}),
            "additional_diffuse accepts file, free, normalization only.",
        )
        item["file"] = str(
            (path.parent / Path(os.path.expandvars(item["file"])).expanduser()).resolve()
        )
        item.setdefault("free", True)
        item.setdefault("normalization", 1.0)
        require(isinstance(item["free"], bool), "Additional diffuse free must be true/false.")
        require(
            np.isfinite(item["normalization"]) and item["normalization"] >= 0,
            "Additional normalization must be finite and nonnegative.",
        )
        extra[name] = item
    # --- Remaining settings and sanity checks -------------------------------
    c["galprop_extrapolation"] = padding_options(c.get("galprop_extrapolation"))
    c.setdefault("isotropic", {"mode": "file", "file": c.get("isotropic_seed")})
    if c["isotropic"]["mode"] == "file":
        require(bool(c["isotropic"].get("file")), "Supply isotropic.file or select powerlaw mode.")
        c["isotropic"]["file"] = str(
            (path.parent / Path(os.path.expandvars(c["isotropic"]["file"])).expanduser()).resolve()
        )
    require(c["emin"] >= 20 and c["emax"] > c["emin"], "Require 20 <= emin < emax (MeV).")
    require(c["tmax"] > c["tmin"], "Invalid observation time range.")
    if "galprop_mode" in c:
        require(c["galprop_mode"] in ["global", "fixed"], "Invalid legacy galprop_mode.")
    c["edisp_bins"] = choose_edisp(
        c["emin"], c["emax"], c["binsperdec"], c.get("edisp_bins", "auto")
    )
    require(bool(c["galprop"]) and "isodiff" not in c["galprop"], "Invalid GALPROP names.")
    m = c["mask"]
    require(
        0 <= m["plane_abs_b_min"] < 90
        and 0 < m["pole_abs_dec_max"] <= 90
        and 0 <= m["bubbles_abs_l_max"] <= 180
        and 0 <= m["bubbles_abs_b_max"] <= 90,
        "Invalid mask boundaries.",
    )
    return c


def template_signature(c):
    """Hash of every setting that affects the prepared templates.

    ``prepare`` stores this hash in ``templates.npz``; ``fit`` recomputes it
    and refuses to run if it differs, so a fit can never silently use
    templates built from other data or model files.

    Included: data/selection/binning/response settings, input file paths, the
    size and modification time of each diffuse cube and Loop I asset, and
    ``PIPELINE_VERSION``.

    Excluded (they can change between fits of the same preparation): the
    mask, free/fixed flags and normalization values, plot settings,
    ``compute_profiles``, ``reuse_from`` and the ``loopI`` block itself (its
    effect is captured through the generated cube paths and stamps).
    """
    # Mask and fit mode can change without rerunning expensive response calculation.
    d = {
        k: v
        for k, v in c.items()
        if k
        not in [
            "mask",
            "galprop_mode",
            "plots",
            "compute_profiles",
            "reuse_from",
            "reuse_source_maps",
            "loopI",
        ]
    }
    if "galprop" in d:
        d["galprop"] = {
            name: item["file"] if isinstance(item, dict) else item
            for name, item in d["galprop"].items()
        }
    if "additional_diffuse" in d:
        d["additional_diffuse"] = {
            name: item["file"] for name, item in d["additional_diffuse"].items()
        }
    # Detect changed input cubes even when their filenames are unchanged.
    d["cube_stamps"] = {
        name: {
            "size": Path(item["file"]).stat().st_size,
            "mtime_ns": Path(item["file"]).stat().st_mtime_ns,
        }
        for name, item in diffuse_components(c).items()
        if Path(item["file"]).is_file()
    }
    from loopI import input_files

    d["loopI_input_stamps"] = {
        str(p): {"size": p.stat().st_size, "mtime_ns": p.stat().st_mtime_ns}
        if p.is_file()
        else None
        for p in input_files(c)
    }
    d["pipeline_version"] = PIPELINE_VERSION
    # The retired galprop_intensity_units_confirmed switch was always true in a
    # completed preparation and was part of this hash. Keep it fixed so
    # existing templates.npz files remain valid.
    d["galprop_intensity_units_confirmed"] = True
    return hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest()


def validate_extra_cube(path, lower, upper):
    """Check that an additional diffuse FITS file is a usable LAT mapcube.

    Requirements: 3D primary image (energy, lat, lon) with a celestial WCS, an
    ``ENERGIES`` table in MeV covering ``[lower, upper]`` (the analysis range
    plus energy-dispersion margins), finite nonnegative values, and (if BUNIT
    is present) differential photon intensity ph cm^-2 s^-1 sr^-1 MeV^-1.
    Unlike GALPROP cubes, these are never padded automatically.
    """
    from astropy.io import fits
    from astropy import units as u
    from astropy.wcs import WCS

    with fits.open(path, memmap=True) as h:
        e = np.asarray(h["ENERGIES"].data.field(0), float).ravel()
        unit = h["ENERGIES"].columns[0].unit
        require(not unit or u.Unit(unit) == u.MeV, "Additional cube energies must be MeV.")
        a = h[0].data
        require(
            a is not None and a.ndim == 3 and len(e) == a.shape[0], "Invalid additional mapcube."
        )
        require(
            len(e) > 1 and np.all(np.isfinite(e)) and np.all(e > 0) and np.all(np.diff(e) > 0),
            "Invalid additional cube energy grid.",
        )
        require(
            e[0] <= lower * (1 + 1e-5) and e[-1] >= upper * (1 - 1e-5),
            f"{path}: generate additional templates covering {lower:g} to {upper:g} MeV for edisp.",
        )
        require(WCS(h[0].header).has_celestial, "Additional cube requires celestial WCS.")
        for plane in a:
            require(
                np.all(np.isfinite(plane)) and np.all(plane >= 0), "Invalid additional intensity."
            )
        require(np.any(a > 0), "Additional intensity cube is entirely zero.")
        if h[0].header.get("BUNIT"):
            require(
                np.isclose(
                    u.Unit(h[0].header["BUNIT"]).to(u.ph / u.cm**2 / u.s / u.sr / u.MeV), 1.0
                ),
                "Additional cube must contain differential photon intensity, not counts or E2 intensity.",
            )


def prepare(c):
    """Build and save the full-sky count templates (``templates.npz``).

    Steps:

    1. Check inputs; build separate Loop I shell cubes if needed.
    2. Pad GALPROP cubes in energy if they do not cover the analysis range plus
       the energy-dispersion margin (copies go to ``padded_galprop/``;
       originals are never modified). Additional cubes must already cover it.
    3. Write the isotropic reference spectrum and a standard Fermipy config
       (``fermipy_config.yaml``) for a full-sky Galactic HEALPix analysis.
    4. Run Fermipy setup (selection, counts cube, livetime, exposure, source
       maps). With ``reuse_from``, compatible products are copied from an
       older output instead of being recomputed.
    5. Freeze all source parameters, set every diffuse normalization to 1,
       and extract the counts cube and one predicted-count template per
       diffuse component plus the summed catalog prediction (``fixed``).
    6. Self-check: templates must add up to Fermipy's total model, and changing
       each normalization by +/-20% must change Fermipy's own log-likelihood
       by the same amount as the explicit Poisson formula used by ``fit``.

    Refuses to run if ``<outdir>/fermipy`` already exists (use a new outdir).
    """
    from astropy.io import fits
    from fermipy.gtanalysis import GTAnalysis

    from loopI import ensure_cubes

    # --- 1. Inputs ------------------------------------------------------------
    ensure_cubes(c)
    galfiles = {name: item["file"] for name, item in c["galprop"].items()}
    for p in (
        [c[k] for k in ["evfile", "scfile", "catalog"]]
        + list(galfiles.values())
        + [v["file"] for v in c.get("additional_diffuse", {}).values()]
    ):
        require(Path(p).is_file(), "Missing input file: " + p)
    require(Path(c["extdir"]).is_dir(), "Missing extended-source template directory.")
    seed = build_seed(c)
    seed_at(np.array([c["emin"], c["emax"]]), seed)
    # True-energy margin needed by energy dispersion: |edisp_bins| bins plus one
    # extra bin of safety on each side of [emin, emax].
    nbins = int(np.round(c["binsperdec"] * np.log10(c["emax"] / c["emin"])))
    padding = 10 ** ((abs(c["edisp_bins"]) + 1) * np.log10(c["emax"] / c["emin"]) / nbins)
    out = Path(c["outdir"])
    out.mkdir(parents=True, exist_ok=True)
    work = out / "fermipy"
    require(
        not work.exists(),
        "Fermipy directory already exists. Use fit for cached templates, "
        "or a NEW outdir to prepare again (including after failed preparation).",
    )
    # --- 2. Energy padding of GALPROP cubes -----------------------------------
    padding_report = {}
    for index, (name, p) in enumerate(galfiles.items()):
        galfiles[name], padding_report[name] = pad_mapcube(
            p,
            out / "padded_galprop" / f"component_{index:02d}.fits",
            c["emin"] / padding,
            c["emax"] * padding,
            c.get("galprop_extrapolation"),
        )
        print(f"GALPROP {name}: {json.dumps(padding_report[name])}", flush=True)
    (out / "galprop_extrapolation.json").write_text(json.dumps(padding_report, indent=2))
    # Additional cubes already contain physical spectra. Require complete true-energy
    # coverage; never apply GALPROP extrapolation settings to solar/lunar profiles.
    for name, item in c.get("additional_diffuse", {}).items():
        validate_extra_cube(item["file"], c["emin"] / padding, c["emax"] * padding)
        galfiles[name] = item["file"]
    # --- 3. Isotropic reference spectrum and Fermipy configuration -----------
    # The isotropic source uses this file as a Fermitools FileFunction.
    seed_path = out / "isotropic_seed_used.txt"
    np.savetxt(seed_path, seed, fmt="%.16e")
    # Full sky (radius 180 deg, no ROI cut), Galactic HEALPix RING pixels.
    # Energy dispersion is on for every source except the isotropic one.
    cfg = {
        "data": {"evfile": c["evfile"], "scfile": c["scfile"]},
        "selection": {
            **{k: c[k] for k in ["emin", "emax", "tmin", "tmax", "zmax", "evclass", "evtype"]},
            "glon": 0.0,
            "glat": 0.0,
            "radius": 180.0,
            "filter": c.get("filter", "DATA_QUAL>0 && LAT_CONFIG==1"),
            "roicut": "no",
        },
        "binning": {
            "projtype": "HPX",
            "coordsys": "GAL",
            "roiwidth": None,
            "binsz": c.get("binsz", 1.0),
            "hpx_order": c["hpx_order"],
            "hpx_ordering_scheme": "RING",
            "binsperdec": c["binsperdec"],
        },
        "gtlike": {
            "irfs": c["irfs"],
            "edisp": True,
            "edisp_bins": c["edisp_bins"],
            "edisp_disable": ["isodiff"],
        },
        "model": {
            "catalogs": [c["catalog"]],
            "extdir": c["extdir"],
            "src_radius": 180.0,
            "galdiff": [{"name": name, "file": p} for name, p in galfiles.items()],
            "isodiff": [str(seed_path)],
        },
        "fileio": {"outdir": str(work), "usescratch": False},
        "logging": c.get("logging", {"chatter": 0, "verbosity": 1}),
    }
    (out / "fermipy_config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
    # --- 4. Fermipy setup (optionally reusing products from reuse_from) -----
    # Fermipy queries conda for the Fermitools version; tolerate its absence.
    import fermipy
    import warnings

    try:
        fermipy.get_st_version()
    except FileNotFoundError as exc:
        if exc.filename != "conda":
            raise
        warnings.warn("Skipping Fermitools version detection: conda unavailable.")
        fermipy.get_st_version = lambda: "unknown (conda unavailable)"
    reuse_state = None
    if c.get("reuse_from"):
        from reuse_products import stage_products, install_source_reuse

        reuse_state = stage_products(c["reuse_from"], work, cfg)
    gta = GTAnalysis(cfg)
    if reuse_state is not None:
        require(len(gta.components) == 1, "Reuse supports one component only.")
        install_source_reuse(gta.components[0], *reuse_state, out / "reuse_report.json")
    # Verify that no catalog sources disappeared through an ROI selection.
    with fits.open(c["catalog"]) as hdus:
        catalog_names = [str(x).strip() for x in hdus[1].data["Source_Name"]]
    for name in catalog_names:
        gta.roi.get_source_by_name(name)  # Raises if missing (also handles aliases).
    gta.setup(init_sources=False)
    # Record source definitions and referenced-file content for subsequent runs.
    from reuse_products import source_key
    import xml.etree.ElementTree as ET

    xmlpath = Path(gta.components[0].files["srcmdl"])
    hashes = {}
    fingerprints = {
        s.get("name"): source_key(s, xmlpath, hashes=hashes)
        for s in ET.parse(xmlpath).getroot().findall("source")
    }
    (out / "source_fingerprints.json").write_text(json.dumps(fingerprints, indent=2))
    require(len(gta.components) == 1, "This version supports one event selection only.")
    gta.free_sources(free=False)
    require(not gta.get_params(freeonly=True), "Some source parameters are still free.")
    # --- 5. Extract templates ------------------------------------------------
    # Template order: isodiff first, then GALPROP, then additional components.
    # Every template is computed at normalization 1 (scale 1).
    names = ["isodiff"] + list(diffuse_components(c))
    for name in names:
        gta.set_norm_scale(name, 1.0)
        gta.set_norm(name, 1.0, update_source=False)
    comp = gta.components[0]
    require(comp.like.logLike.use_edisp(), "Native likelihood reports energy dispersion disabled.")
    cmap = comp.counts_map()
    counts = np.asarray(cmap.data, dtype=float)
    npix = 12 * (2 ** c["hpx_order"]) ** 2
    require(counts.shape == (gta.enumbins, npix), "Count cube is not full-sky HEALPix.")
    require(np.allclose(counts, np.rint(counts)), "Expected integer observed counts.")
    # Derive coordinates from map geometry: no assumptions about RING/NESTED indexing.
    sky = cmap.geom.to_image().get_coord().skycoord
    l = np.asarray(sky.galactic.l.deg).ravel()
    b = np.asarray(sky.galactic.b.deg).ravel()
    dec = np.asarray(sky.icrs.dec.deg).ravel()
    require(l.size == npix, "Coordinate array does not match count pixels.")

    def prediction(name=None, exclude=None):
        """Full-sky predicted counts for one source (``name``), for all sources
        except ``exclude``, or for the whole model. Tiny negative values from
        numerical noise are set to zero; larger ones are an error."""
        m = comp.model_counts_map(name=name, exclude=exclude, use_mask=False)
        require(m.data.shape == counts.shape, "Model/count grid mismatch.")
        a = np.array(m.data, dtype=float, copy=True)
        require(
            np.all(np.isfinite(a)),
            f"Nonfinite predictions for {name!r}, exclude={exclude!r}: "
            f"NaN={np.isnan(a).sum()}, Inf={np.isinf(a).sum()}",
        )
        neg = a < 0
        if np.any(neg):
            minimum = float(a.min())
            mass = float(-a[neg].sum())
            positive = float(a[a > 0].sum())
            print(
                f"Negative predictions: min={minimum:g}, mass={mass:g}, "
                f"pixels={neg.sum()}, source={name!r}",
                flush=True,
            )
            require(
                minimum >= -1e-6 and mass <= 1e-5 * max(positive, 1.0),
                "Negative predictions exceed numerical-cleanup limits.",
            )
            a[neg] = 0.0
        return a

    templates = np.stack([prediction(name) for name in names])
    # `fixed` = everything that is not a diffuse template (4FGL point and
    # extended sources), held at catalog values.
    fixed = prediction(exclude=names)
    total = fixed + templates.sum(axis=0)
    np.testing.assert_allclose(total, prediction(), rtol=2e-4, atol=1e-5)
    # Installed-environment integration gate: compare likelihood CHANGES, so
    # parameter-independent Poisson constants cancel. Detect response/map mistakes.
    # --- 6. Native likelihood check ------------------------------------------
    checks = []
    baseline = float(gta.like())
    for j, name in enumerate(names):
        for factor in [0.8, 1.2]:
            gta.set_norm(name, factor, update_source=False)
            native = float(gta.like()) - baseline
            mu = total + (factor - 1.0) * templates[j]
            explicit = float(
                np.sum(
                    mu
                    - total
                    - xlogy(counts, np.divide(mu, total, out=np.ones_like(mu), where=total > 0))
                )
            )
            require(
                abs(native - explicit) < 0.05 + 0.002 * abs(native),
                f"Native likelihood validation failed for {name}: {native} vs {explicit}",
            )
            checks.append(
                {
                    "name": name,
                    "factor": factor,
                    "native_delta_nll": native,
                    "explicit_delta_nll": explicit,
                }
            )
            gta.set_norm(name, 1.0, update_source=False)
    # --- Save --------------------------------------------------------------
    meta = {
        "signature": template_signature(c),
        "config": c,
        "catalog_count": len(catalog_names),
        "native_likelihood_checks": checks,
        "versions": {},
    }
    for name in ["fermipy", "numpy", "scipy", "astropy", "gammapy"]:
        try:
            meta["versions"][name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            pass
    np.savez_compressed(
        out / "templates.npz",
        counts=counts,
        fixed=fixed,
        templates=templates,
        names=np.array(names),
        edges=gta.energies,
        l=l,
        b=b,
        dec=dec,
        seed=seed,
        metadata=json.dumps(meta),
    )
    (out / "preparation_checks.json").write_text(json.dumps(meta, indent=2))
    print("Prepared templates and passed native likelihood checks:", out, flush=True)


def result_dir(c):
    """Folder for the results of one fit configuration.

    Name: ``bestfit_<hash>`` (or ``fit_<hash>`` when ``compute_profiles`` is
    true). The hash covers the mask, free/fixed flags, normalizations and the
    ``loopI`` block, so a different fit setup never overwrites an older one,
    while plot-only changes map to the same folder.
    """
    controls = {
        name: {k: item.get(k, default) for k, default in [("free", True), ("normalization", 1.0)]}
        for name, item in diffuse_components(c).items()
        if isinstance(item, dict)
    }
    key_parts = [PIPELINE_VERSION, c["mask"], c.get("galprop_mode"), controls]
    if "loopI" in c:
        key_parts.append(c["loopI"])
    key = json.dumps(key_parts, sort_keys=True)
    prefix = "fit_" if c.get("compute_profiles", False) else "bestfit_"
    return Path(c["outdir"]) / (prefix + hashlib.sha256(key.encode()).hexdigest()[:10])


def fit(c):
    """Fit the prepared templates and write all result files.

    Uses ``<outdir>/templates.npz`` only; Fermipy is not needed unless an
    analytic Loop I spectral shape is fitted. Returns the result folder.

    Main outputs (see README): ``isotropic_bins.csv``,
    ``diffuse_normalizations.json``, ``covariance.npz``, ``fit_report.json``,
    ``diagnostics.npz`` and, when every bin is positive,
    ``isotropic_spectrum.txt``.
    """
    out = Path(c["outdir"])
    with np.load(out / "templates.npz", allow_pickle=False) as z:
        meta = json.loads(str(z["metadata"]))
        require(
            meta["signature"] == template_signature(c),
            "Preparation inputs or pipeline version changed. Use a NEW outdir and prepare again.",
        )
        # keep[p] is True for pixels used in the likelihood.
        keep = sky_mask(z["l"], z["b"], z["dec"], c["mask"])
        require(np.any(keep), "Mask excludes every pixel.")
        results = result_dir(c)
        require(
            not results.exists(),
            f"Results already exist: {results}. Move this directory aside to refit, or use plot.",
        )
        results.mkdir()
        # Energy bins: edges from Fermipy; reported energies are geometric centres.
        edges = z["edges"]
        centers = np.sqrt(edges[:-1] * edges[1:])
        ne = len(centers)
        seed = seed_at(centers, z["seed"])
        counts, fixed, templates = z["counts"], z["fixed"], z["templates"]
        from loopI import analytic_enabled, optimize_templates

        # Optional: fit an analytic Loop I spectral shape first. This replaces the
        # Loop I templates with ones refolded through the LAT response.
        spectral_report = None
        if analytic_enabled(c):
            templates, spectral_report = optimize_templates(
                c, z, keep, sys.modules[__name__], results
            )
        # Build the linear problem on retained pixels only:
        #   f : catalog + every FIXED diffuse component (times its normalization)
        #   t : [isotropic template, FREE diffuse templates...]
        f = fixed[:, keep].copy()
        all_t = templates[:, :, keep]
        free, gal = galprop_controls(c, z["names"][1:])
        f += np.einsum("j,jkp->kp", gal[~free], all_t[1:][~free])
        t = np.concatenate([all_t[:1], all_t[1:][free]], axis=0)
        # Parameters: one isotropic amplitude per energy bin, then one global
        # amplitude per free diffuse component (see joint_fit.JointPoisson).
        problem = JointPoisson(counts[:, keep], f, t)
        print(
            "Optimizing best-fit spectrum from cached templates (no Fermipy setup)...", flush=True
        )
        best, minimum = problem.solve(np.r_[np.ones(ne), gal[free]])
        np.savez_compressed(results / "best_fit_parameters.npz", best=best, minimum=minimum)
        print("Main optimization completed; parameters checkpointed.", flush=True)
        # Inverse Hessian of the Poisson likelihood (NaN if singular).
        covariance = problem.covariance(best)
        gal[free] = best[ne:]
        parameter_names = [f"isodiff_bin{k:03d}" for k in range(ne)]
        parameter_names += list(z["names"][1:][free])
        # Full-sky best-fit model (used for per-bin summaries and diagnostics).
        model = fixed + best[:ne, None] * templates[0] + np.einsum("j,jkp->kp", gal, templates[1:])
        # --- Per-energy-bin isotropic results -> isotropic_bins.csv -------------
        # intensity = fitted amplitude * reference intensity at the bin centre.
        rows = []
        print("Shared diffuse normalizations:", dict(zip(z["names"][1:], gal)), flush=True)
        for k, energy in enumerate(centers):
            if c.get("compute_profiles", False):
                print(f"Profiling isotropic bin {k + 1}/{ne}...", flush=True)
                lo, hi, ul, ts = problem.profile(k, best, minimum)
            else:
                lo = hi = ul = ts = np.nan
            n = counts[k, keep]
            mu = model[k, keep]
            row = dict(
                emin_MeV=edges[k],
                emax_MeV=edges[k + 1],
                energy_MeV=energy,
                intensity=seed[k] * best[k],
                intensity_lo=seed[k] * lo,
                intensity_hi=seed[k] * hi,
                intensity_ul95=seed[k] * ul,
                iso_norm=best[k],
                iso_TS=ts,
                counts=n.sum(),
                predicted_counts=mu.sum(),
                deviance=2 * (nll(n, mu) - nll(n, n)),
                boundary=bool(best[k] < 1e-7),
                covariance_valid=bool(np.all(np.isfinite(covariance))),
            )
            rows.append(row)
            print(f"{energy / 1000:.3g} GeV: I={row['intensity']:.5g}", flush=True)
        with (results / "isotropic_bins.csv").open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=rows[0].keys())
            w.writeheader()
            w.writerows(rows)
        np.savez_compressed(
            results / "covariance.npz",
            covariance=covariance,
            names=np.array(parameter_names),
            best=best,
        )
        # --- Diffuse normalizations -> *_normalizations.json --------------------
        gal_rows = []
        for j, name in enumerate(z["names"][1:]):
            item = {"name": str(name), "normalization": float(gal[j]), "fixed": not bool(free[j])}
            if not item["fixed"] and c.get("compute_profiles", False):
                lo, hi, ul, ts = problem.profile(ne + int(free[:j].sum()), best, minimum)
                item.update(lo68=lo, hi68=hi, ul95=ul, TS=ts)
            gal_rows.append(item)
        (results / "diffuse_normalizations.json").write_text(json.dumps(gal_rows, indent=2))
        (results / "galprop_normalizations.json").write_text(
            json.dumps([r for r in gal_rows if r["name"] in c.get("galprop", {})], indent=2)
        )
        (results / "additional_diffuse_normalizations.json").write_text(
            json.dumps(
                [r for r in gal_rows if r["name"] in c.get("additional_diffuse", {})], indent=2
            )
        )
        # Retain full-sky predictions; excluded pixels are diagnostic only and are
        # NEVER included in fit statistics. Plotting marks them as excluded.
        np.savez_compressed(
            results / "diagnostics.npz",
            counts=counts,
            model=model,
            fixed_catalog=fixed,
            isotropic=best[:ne, None] * templates[0],
            galactic=gal[:, None, None] * templates[1:],
            galactic_names=z["names"][1:],
            keep=keep,
            l=z["l"],
            b=z["b"],
            dec=z["dec"],
            edges=edges,
        )
        # --- Summary report ----------------------------------------------------
        report = {
            "loopI_spectral_fit": spectral_report,
            "config": c,
            "preparation": meta,
            "retained_pixels": int(keep.sum()),
            "solid_angle_sr": float(4 * np.pi * keep.mean()),
            "deviance": 2 * minimum,
            "diffuse_normalizations": gal_rows,
            "intensity_units": "ph cm^-2 s^-1 sr^-1 MeV^-1",
            "isotropic_convention": "effective background; apply_edisp=false on reuse",
            "error_method": (
                "Joint profile likelihood"
                if c.get("compute_profiles", False)
                else "Best-fit values only; no profile intervals, limits, or TS. Hessian covariance saved separately."
            ),
            "text_export": "withheld",
        }
        vals = np.array([r["intensity"] for r in rows])
        # Export the text spectrum only if every bin is positive (and, with
        # profiles, significant: TS >= 4 with a valid covariance).
        if np.all(vals > 0) and (
            not c.get("compute_profiles", False)
            or all(r["iso_TS"] >= 4 and r["covariance_valid"] for r in rows)
        ):
            np.savetxt(
                results / "isotropic_spectrum.txt",
                np.column_stack([centers, vals]),
                fmt="%.10e",
                header="Energy [MeV] intensity [ph cm^-2 s^-1 sr^-1 MeV^-1]\n"
                "Best-fit values; this file does not establish per-bin detections.\n"
                "Effective background: disable energy dispersion for this component.\n"
                "Log-log interpolation approximation; valid only between first and last energies.",
            )
            report["text_export"] = "isotropic_spectrum.txt"
        else:
            print("Text export withheld: weak/zero bin or singular covariance. Use CSV/limits.")
        if spectral_report and spectral_report["shape_free"]:
            report["error_method"] = (
                "Best fit only. Saved linear covariance is conditional on the fitted Loop I shape, not full joint covariance."
            )
        # Loop I coefficient spectra and amplitude correlations (if enabled).
        from loopI import export_spectrum, export_correlations

        export_correlations(c, parameter_names, covariance, results)
        export_spectrum(c, centers, gal_rows, spectral_report, results)
        (results / "fit_report.json").write_text(json.dumps(report, indent=2))
        print("Results:", results)
    return results


def plot(c):
    """Make the spectrum, count-spectrum and residual-map plots for the
    result folder selected by the config (the fit must already exist)."""
    from plots import plot_results

    plot_results(result_dir(c), c.get("plots", {}))


class Validation(unittest.TestCase):
    """Synthetic checks run by ``python fit_isotropic.py self-test``.

    They exercise the numerical fit (recovery of known amplitudes, profile
    intervals, boundaries, covariance), the sky mask, the reference spectrum,
    the cached fit/export path and plotting. They do not touch LAT data or
    Fermitools, so they validate the code, not the LAT response.
    """

    def test_reference_normalization_invariance(self):
        """Rescaling the isotropic template rescales its amplitude inversely."""
        c = {
            "emin": 20.0,
            "emax": 800000.0,
            "isotropic": {
                "mode": "powerlaw",
                "index": -2.3,
                "pivot_MeV": 1000.0,
                "intensity_at_pivot": 1e-10,
            },
        }
        reference = build_seed(c)
        np.testing.assert_allclose(
            seed_at([20.0, 1000.0, 800000.0], reference),
            1e-10 * (np.array([20.0, 1000.0, 800000.0]) / 1000) ** -2.3,
        )
        t = np.array([[2.0, 3.0, 4.0, 5.0]])
        n = np.array([4.0, 6.0, 8.0, 10.0])
        f = np.ones(4) * 0.5
        a, _, _ = solve(n, f, t)
        scaled, _, _ = solve(n, f, t * 37.0)
        np.testing.assert_allclose(scaled[0] * 37.0, a[0], rtol=1e-6)

    def test_independent_galprop_controls(self):
        """Free/fixed flags per GALPROP component are honoured in a cached fit."""
        rng = np.random.default_rng(12)
        with tempfile.TemporaryDirectory() as folder:
            c = {
                "outdir": folder,
                "mask": dict(
                    plane_abs_b_min=10,
                    pole_abs_dec_max=60,
                    bubbles_abs_l_max=30,
                    bubbles_abs_b_max=60,
                ),
                "galprop": {
                    "IC": {"file": "unused_ic", "free": True, "normalization": 0.6},
                    "pi0": {"file": "unused_pi0", "free": True, "normalization": 0.5},
                    "brem": {"file": "unused_brem", "free": False, "normalization": 1.35},
                },
            }
            t = np.concatenate([np.full((1, 2, 300), 5.0), rng.uniform(0.1, 12.0, (3, 2, 300))])
            f = np.ones((2, 300))
            counts = (
                f
                + np.array([1.7, 1.4])[:, None] * t[0]
                + np.einsum("j,jkp->kp", [0.8, 1.2, 1.35], t[1:])
            )
            np.savez_compressed(
                Path(folder) / "templates.npz",
                counts=counts,
                fixed=f,
                templates=t,
                names=np.array(["isodiff", "IC", "pi0", "brem"]),
                edges=[1000.0, 2000.0, 4000.0],
                l=np.full(300, 180.0),
                b=np.full(300, 30.0),
                dec=np.zeros(300),
                seed=[[1000.0, 1e-9], [4000.0, 1e-10]],
                metadata=json.dumps({"signature": template_signature(c)}),
            )
            original = template_signature(c)
            original_dir = result_dir(c)
            c["galprop"]["brem"]["free"] = True
            self.assertEqual(original, template_signature(c))
            self.assertNotEqual(original_dir, result_dir(c))
            c["galprop"]["brem"]["free"] = False
            directory = fit(c)
            result = json.loads((directory / "galprop_normalizations.json").read_text())
            np.testing.assert_allclose(
                [r["normalization"] for r in result], [0.8, 1.2, 1.35], rtol=1e-5
            )
            self.assertEqual([r["fixed"] for r in result], [False, False, True])
            with np.load(directory / "covariance.npz") as z:
                self.assertEqual(z["covariance"].shape, (4, 4))
                self.assertNotIn("brem", z["names"])

    def test_low_energy_settings(self):
        """edisp_bins auto-selection and seed extrapolation rules."""
        self.assertEqual(choose_edisp(300, 100000, 8, "auto"), -2)
        self.assertEqual(choose_edisp(20, 100000, 8, "auto"), 4)
        with self.assertRaises(ValueError):
            choose_edisp(20, 100000, 8, -2)
        with self.assertRaises(ValueError):
            choose_edisp(20, 100000, 16, "auto")
        table = np.array([[100.0, 1e-8], [1000.0, 1e-10]])
        extended = extend_seed(table, 20.0, 1000.0, True)
        self.assertAlmostEqual(extended[0, 1] / 2.5e-7, 1.0)
        with self.assertRaises(ValueError):
            extend_seed(table, 20.0, 1000.0, False)

    def test_joint_global_shape_and_derivatives(self):
        """Joint fit recovers truth; analytic gradient/Hessian match finite differences."""
        rng = np.random.default_rng(35)
        # Redistribute TRUE-energy Galactic templates before fitting their one
        # shared normalization. This is a toy response, not a LAT IRF test.
        response = np.array([[0.8, 0.15, 0], [0.2, 0.7, 0.2], [0, 0.15, 0.8]])
        true_templates = rng.uniform(1, 10, (2, 3, 100))
        folded = np.einsum("kt,jtp->jkp", response, true_templates)
        t = np.concatenate([np.full((1, 3, 100), 4.0), folded])
        f = np.full((3, 100), 2.0)
        truth = np.array([1.2, 1.7, 0.8, 0.75, 1.4])
        n = f + truth[:3, None] * t[0] + np.einsum("j,jkp->kp", truth[3:], t[1:])
        problem = JointPoisson(n, f, t)
        a, value = problem.solve()
        self.assertEqual(len(a), 5)  # three isotropic bins, TWO global Galactic norms.
        np.testing.assert_allclose(a, truth, rtol=1e-5)
        point = np.array([1.1, 1.2, 1.3, 0.9, 1.1])
        h = 1e-5
        numerical = []
        for i in range(len(point)):
            d = np.eye(len(point))[i] * h
            numerical.append(
                (problem.objective(point + d)[0] - problem.objective(point - d)[0]) / (2 * h)
            )
        np.testing.assert_allclose(problem.objective(point)[1], numerical, rtol=1e-5, atol=1e-5)
        cov = problem.covariance(a)
        nh = np.column_stack(
            [
                (
                    problem.objective(a + np.eye(len(a))[i] * h)[1]
                    - problem.objective(a - np.eye(len(a))[i] * h)[1]
                )
                / (2 * h)
                for i in range(len(a))
            ]
        )
        np.testing.assert_allclose(cov, np.linalg.inv(nh), rtol=1e-5, atol=1e-7)
        self.assertGreater(abs(cov[0, 1]), 1e-6)  # shared Galactic nuisance couples bins.
        lo, hi, ul, ts = problem.profile(0, a, value)
        self.assertLess(lo, truth[0])
        self.assertGreater(hi, truth[0])
        self.assertGreater(ul, hi)

    def test_joint_boundary(self):
        """Zero-count data give a zero amplitude and the expected upper interval."""
        problem = JointPoisson(np.zeros((2, 20)), np.ones((2, 20)), np.ones((1, 2, 20)))
        a, v = problem.solve()
        np.testing.assert_allclose(a, 0, atol=1e-8)
        lo, hi, ul, ts = problem.profile(0, a, v)
        self.assertAlmostEqual(hi, 0.5 / 20.0, places=6)

    def test_residual_definitions(self):
        """Count, fractional and Pearson residual formulas."""
        from plots import residuals

        d, f, p = residuals(np.array([12.0, 6.0, 0.0]), np.array([9.0, 4.0, 0.0]))
        np.testing.assert_allclose(d, [3, 2, 0])
        np.testing.assert_allclose(f[:2], [1 / 3, 1 / 2])
        np.testing.assert_allclose(p[:2], [1, 1])
        self.assertTrue(np.isnan(f[2]))

    def test_mask(self):
        """Plane, pole and bubble-rectangle cuts select the right pixels."""
        m = dict(
            plane_abs_b_min=10, pole_abs_dec_max=60, bubbles_abs_l_max=30, bubbles_abs_b_max=60
        )
        np.testing.assert_array_equal(
            sky_mask([0, 359, 180, 180, 180, 0], [30, -30, 0, 30, 30, 70], [0, 0, 0, 61, 0, 0], m),
            [False, False, False, False, True, True],
        )

    def test_seed_units(self):
        """Log-log interpolation of the reference spectrum; no extrapolation."""
        tab = np.array([[1000.0, 1e-9], [10000.0, 1e-11]])
        self.assertAlmostEqual(seed_at([np.sqrt(1e7)], tab)[0] / 1e-10, 1.0)
        with self.assertRaises(ValueError):
            seed_at([100.0], tab)

    def test_analytic_poisson(self):
        """Single-template fit matches the analytic Poisson solution."""
        t = np.array([[1.0, 2.0, 3.0, 4.0]])
        n = np.array([3.0, 7.0, 8.0, 14.0])
        a, v, cov = solve(n, np.zeros(4), t)
        self.assertAlmostEqual(a[0], n.sum() / t.sum(), places=5)
        self.assertAlmostEqual(cov[0, 0], a[0] / t.sum(), places=5)
        lo, hi, ul, ts = profile_iso(n, np.zeros(4), t, a, v)
        self.assertLess(lo, a[0])
        self.assertGreater(hi, a[0])
        self.assertGreater(ul, hi)

    def test_asimov_recovery_and_profile(self):
        """Noise-free data recover the true amplitudes; profile brackets truth."""
        rng = np.random.default_rng(41)
        t = np.vstack([np.full(800, 3.0), rng.uniform(0, 10, 800), rng.uniform(0, 4, 800)])
        truth = np.array([1.8, 0.75, 1.4])
        f = rng.uniform(0.1, 5, 800)
        n = f + truth @ t
        a, v, cov = solve(n, f, t)
        np.testing.assert_allclose(a, truth, rtol=2e-5)
        lo, hi, ul, _ = profile_iso(n, f, t, a, v)
        self.assertLess(lo, truth[0])
        self.assertGreater(hi, truth[0])
        self.assertTrue(np.all(np.linalg.eigvalsh(cov) > 0))
        self.assertGreater(ul, hi)

    def test_masked_contamination(self):
        """Masked pixels have no influence on the fit."""
        t = np.array([[2.0, 3.0, 4.0, 5.0]])
        n = np.array([4.0, 6.0, 8.0, 1e9])
        keep = np.array([1, 1, 1, 0], bool)
        a, _, _ = solve(n[keep], np.zeros(3), t[:, keep])
        self.assertAlmostEqual(a[0], 2.0, places=5)

    def test_boundary(self):
        """Boundary case: amplitude 0, 68% and 95% upper limits are analytic."""
        n = np.zeros(10)
        f = np.ones(10)
        t = np.ones((1, 10))
        a, v, _ = solve(n, f, t)
        self.assertAlmostEqual(a[0], 0.0, places=7)
        lo, hi, ul, ts = profile_iso(n, f, t, a, v)
        self.assertEqual(lo, 0.0)
        self.assertAlmostEqual(hi, 0.05, places=6)
        self.assertAlmostEqual(ul, 0.1352771727, places=6)

    def test_degenerate_templates(self):
        """Identical templates give a NaN (not a fake) covariance."""
        t = np.ones((2, 20))
        _, _, cov = solve(np.full(20, 4.0), np.ones(20), t)
        self.assertTrue(np.all(np.isnan(cov)))

    def test_invalid_model_support(self):
        """Counts where the model is zero are rejected."""
        with self.assertRaises(ValueError):
            solve(np.ones(3), np.zeros(3), np.zeros((1, 3)))

    def test_poisson_injections(self):
        """Pull distribution from Poisson simulations has mean ~0 and width ~1."""
        rng = np.random.default_rng(100)
        t = np.vstack([np.full(500, 4.0), rng.uniform(0, 15, 500)])
        f = np.full(500, 2.0)
        truth = np.array([1.6, 0.8])
        pulls = []
        for _ in range(40):
            n = rng.poisson(f + truth @ t)
            a, _, cov = solve(n, f, t)
            pulls.append((a[0] - truth[0]) / np.sqrt(cov[0, 0]))
        self.assertLess(abs(np.mean(pulls)), 0.5)
        self.assertTrue(0.6 < np.std(pulls) < 1.4)

    def test_cached_fit_and_intensity_export(self):
        """End-to-end cached fit, text export and plots on a toy templates.npz."""
        rng = np.random.default_rng(19)
        with tempfile.TemporaryDirectory() as folder:
            c = {
                "outdir": folder,
                "mask": dict(
                    plane_abs_b_min=10,
                    pole_abs_dec_max=60,
                    bubbles_abs_l_max=30,
                    bubbles_abs_b_max=60,
                ),
                "galprop_mode": "global",
            }
            t = np.stack(
                [np.full((2, 400), 8.0), rng.uniform(1, 12, (2, 400)), rng.uniform(0, 5, (2, 400))]
            )
            f = np.full((2, 400), 2.0)
            truth = np.array([1.7, 1.0, 1.0])
            counts = f + np.einsum("j,jkp->kp", truth, t)
            edges = np.array([1000.0, 2000.0, 4000.0])
            seed = np.array([[1000.0, 1e-9], [4000.0, 1e-10]])
            np.savez_compressed(
                Path(folder) / "templates.npz",
                counts=counts,
                fixed=f,
                templates=t,
                names=np.array(["isodiff", "GP1", "GP2"]),
                edges=edges,
                l=np.full(400, 180.0),
                b=np.full(400, 30.0),
                dec=np.zeros(400),
                seed=seed,
                metadata=json.dumps({"signature": template_signature(c)}),
            )
            for mode in ["global", "fixed"]:
                c["galprop_mode"] = mode
                fit(c)
            # Plot the last fit to exercise spectrum and map output without LAT dependencies.
            c["plots"] = {"map_nlon": 60}
            plot(c)
            self.assertTrue((result_dir(c) / "isotropic_spectrum.png").is_file())
            self.assertTrue(list(result_dir(c).glob("residuals_*.png")))
            for p in Path(folder).glob("bestfit_*/isotropic_spectrum.txt"):
                a = np.loadtxt(p)
                np.testing.assert_allclose(a[:, 1], 1.7 * seed_at(a[:, 0], seed), rtol=1e-5)
            self.assertEqual(len(list(Path(folder).glob("bestfit_*/isotropic_spectrum.txt"))), 2)


def main():
    """Command-line entry point; see the module docstring for the commands."""
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "command",
        choices=["self-test", "prepare", "fit", "plot", "run"],
        help="what to do (run = prepare, then fit, then plot)",
    )
    p.add_argument("config", nargs="?", help="wrapper config YAML (not needed for self-test)")
    a = p.parse_args()
    if a.command == "self-test":
        suite = unittest.defaultTestLoader.loadTestsFromTestCase(Validation)
        return 0 if unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful() else 1
    require(a.config is not None, "Supply config.yaml.")
    c = read_config(a.config)
    if a.command in ["prepare", "run"]:
        prepare(c)
    if a.command in ["fit", "run"]:
        fit(c)
    if a.command in ["plot", "run"]:
        plot(c)
    return 0


if __name__ == "__main__":
    sys.exit(main())
