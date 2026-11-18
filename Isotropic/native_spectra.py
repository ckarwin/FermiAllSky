"""Fermitools helpers for fitting an analytic Loop I spectral shape.

Used only when ``loopI.spectrum`` is ``powerlaw``, ``cutoff_powerlaw`` or
``logparabola``. Then ``prepare`` builds Loop I templates from a "seed" cube
with a fixed E^-2.4 spectrum, and ``fit`` changes the spectral shape by
attaching a multiplicative Fermitools spectrum (desired / seed) to that source
and asking Fermitools for new predicted counts. This keeps exposure, PSF and
energy dispersion exactly consistent with the rest of the model.
"""

import json
import numpy as np


def check_seed(path, s):
    """Check that a Loop I seed cube is what the spectral fit assumes.

    The cube must be a full-sky Galactic CAR map whose every energy plane is
    ``spatial_map * seed_flux * (E/pivot)**(-seed_index)``. Otherwise the
    multiplicative correction would not produce the requested spectrum.
    """
    from astropy.io import fits
    from astropy.wcs import WCS

    with fits.open(path, memmap=True) as h:
        e = np.asarray(h["ENERGIES"].data.field(0), float)
        data = h[0].data
        ny, nx = data.shape[1:]
        w = WCS(h[0].header).celestial
        if list(w.wcs.ctype) != ["GLON-CAR", "GLAT-CAR"] or not np.allclose(w.wcs.crval, [0, 0]):
            raise ValueError("Expected the supplied Galactic CAR seed cubes")
        y, x = np.indices((ny, nx))
        _, b = w.all_pix2world(x, y, 0)
        dx, dy = abs(w.wcs.cdelt[0]), abs(w.wcs.cdelt[1])
        om = np.deg2rad(dx) * (np.sin(np.deg2rad(b + dy / 2)) - np.sin(np.deg2rad(b - dy / 2)))
        np.testing.assert_allclose(om.sum(), 4 * np.pi, rtol=1e-6)
        flux = s["seed_flux"] * (e / s["pivot_MeV"]) ** (-s["seed_index"])
        spatial = np.asarray(data[len(e) // 2], float) / flux[len(e) // 2]
        for k, plane in enumerate(data):
            np.testing.assert_allclose(
                np.sum(plane * om),
                flux[k],
                rtol=2e-5,
                err_msg="Seed cube flux differs from spectral YAML seed parameters",
            )
            np.testing.assert_allclose(
                plane,
                flux[k] * spatial,
                rtol=2e-5,
                atol=0,
                err_msg="Seed cube is not spatially separable with the specified power law",
            )


class NativePredictor:
    """Return Fermitools-predicted counts of a Loop I source for a given shape.

    ``predictor(name, spec)`` sets the analytic spectrum ``spec`` on source
    ``name``, recomputes its model counts with the native LAT response and
    returns them on the retained pixels. Results are cached per source.
    """

    def __init__(self, gta, keep, shape):
        self.gta = gta
        self.comp = gta.components[0]
        self.keep = keep
        self.shape = shape
        self.cache = {}
        self.calls = 0

    def full(self, name, s):
        """Full-sky predicted counts for ``name`` with spectral shape ``s``."""
        from spectral_core import native_spectrum, ratio

        st, pars = native_spectrum(s)
        self.gta.set_source_spectrum(name, st, pars, update_source=False)
        self.gta.set_edisp_flag(name, True)
        # Confirm installed Fermitools evaluates the intended multiplicative ratio.
        import pyLikelihood as pyLike

        energies = np.geomspace(300, 100000, 12)
        actual = np.array(
            [self.comp.like[name].spectrum()(pyLike.dArg(float(e))) for e in energies]
        )
        np.testing.assert_allclose(actual, ratio(energies, s), rtol=1e-7, atol=1e-250)
        data = np.asarray(self.comp.model_counts_map(name=name, use_mask=False).data, float).copy()
        if data.shape != self.shape or np.any(~np.isfinite(data)) or data.min() < -1e-6:
            raise ValueError("Invalid native source prediction")
        data = np.maximum(data, 0.0)
        if not data[:, self.keep].sum() > 0:
            raise ValueError("Spectral source has no retained model support")
        self.calls += 1
        return data

    def __call__(self, name, s):
        """Predicted counts on retained pixels (cached until the shape changes)."""
        key = json.dumps(s, sort_keys=True)
        if name not in self.cache or self.cache[name][0] != key:
            self.cache[name] = (key, self.full(name, s)[:, self.keep])
            if self.calls % 10 == 0:
                print(f"Native shell prediction updates: {self.calls}", flush=True)
        return self.cache[name][1]
