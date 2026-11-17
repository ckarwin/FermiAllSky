#!/usr/bin/env python3
"""Make an intrinsic, unit-integral NFW annihilation template for LAT/3ML.

Dependencies: numpy, scipy, astropy
Example:
    python make_nfw_template.py --output nfw_gamma1.fits
    python make_nfw_template.py --gamma 1.2 --output nfw_gamma1p2.fits

The primary HDU is a 2D equatorial CEA (equal-area) WCS image in sr^-1.
Load with SpatialTemplate_2D(fits_file=..., frame='icrs'), fixing its K=1.
Attach a separate empirical spectral model to the ExtendedSource.

Density: rho(r) proportional to (r/rs)^(-gamma)*(1+r/rs)^(gamma-3).
Intensity: J(psi) = integral rho(r)^2 ds, from the observer to the halo edge.
No particle mass, cross section, absolute J-factor, exposure or PSF is included.
The arbitrary density normalization cancels when the map is normalized.
The map is zero outside its footprint. Spectral flux refers to that footprint,
not the whole halo or just the likelihood ROI. Choose a footprint that contains
the likelihood ROI plus an energy-dependent PSF margin; check edge stability.

The central cusp is pixel-averaged with enhanced Gauss-Legendre quadrature.
There is no physical core or arbitrary central-pixel replacement. The exact
central ray diverges, but the pixel integral is finite for gamma < 1.5.
For convergence checks, double --cusp-order and/or halve --binsz.
"""

import argparse
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.interpolate import PchipInterpolator


def make_nfw_template(output="nfw_gamma1.fits", *, gamma=1.0, rs=20.0,
                      r_sun=8.5, r_halo=250.0, width=60.0, binsz=0.1,
                      cusp_order=256, overwrite=False):
    """Write and return a FITS path; distances in kpc, angles in degrees.

    width is the square projection-plane width, not a Galactic l/b box.
    binsz is the pixel scale at the projection center. Both can be adjusted.
    The halo cutoff is spherical about the Galactic center.
    """
    output = Path(output).expanduser()
    if output.exists() and not overwrite:
        raise FileExistsError(f"{output} exists; use --overwrite to replace it")
    values = (gamma, rs, r_sun, r_halo, width, binsz)
    if not np.all(np.isfinite(values)):
        raise ValueError("All numerical inputs must be finite")
    if not 0 < gamma < 1.5:
        raise ValueError("Require 0 < gamma < 1.5 for an integrable central cusp")
    if not (rs > 0 and 0 < r_sun < r_halo):
        raise ValueError("Require rs > 0 and 0 < r_sun < r_halo")
    if not (0 < binsz <= 1 and 2*binsz <= width <= 80):
        raise ValueError("Require 0 < binsz <= 1 and 2*binsz <= width <= 80 deg")
    if int(cusp_order) != cusp_order or cusp_order < 32:
        raise ValueError("cusp_order must be an integer >= 32")
    cusp_order = int(cusp_order)

    center = SkyCoord(l=0, b=0, unit="deg", frame="galactic").icrs
    ra0, dec0 = center.ra.deg, center.dec.deg
    n = 2 * int(np.ceil(width / binsz / 2))
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---CEA", "DEC--CEA"]
    w.wcs.cunit = ["deg", "deg"]
    w.wcs.crval = [ra0, dec0]
    w.wcs.crpix = [(n+1)/2, (n+1)/2]
    w.wcs.cdelt = [-binsz, binsz]
    w.wcs.set_pv([(2, 1, 1.0)])  # CEA lambda=1: constant solid angle.
    w.wcs.radesys = "ICRS"

    def separation(x, y):
        ra, dec = w.all_pix2world(x, y, 0)
        dra = np.deg2rad(ra-ra0)
        dec = np.deg2rad(dec)
        d0 = np.deg2rad(dec0)
        hav = np.sin((dec-d0)/2)**2 + np.cos(dec)*np.cos(d0)*np.sin(dra/2)**2
        return np.rad2deg(2*np.arcsin(np.sqrt(np.clip(hav, 0, 1))))

    def los(psi_deg):
        # s - R_sun*cos(psi) = impact*sinh(t) resolves the cusp without
        # a coarse uniform LOS grid. Then r = impact*cosh(t), ds = r dt.
        psi = np.deg2rad(psi_deg)
        impact = r_sun*np.sin(psi)
        along = r_sun*np.cos(psi)
        t0 = np.arcsinh(-along/impact)
        t1 = np.arcsinh(np.sqrt(r_halo**2-impact**2)/impact)

        def integrand(t):
            radius = impact*np.cosh(t)
            q = radius/rs
            return radius*q**(-2*gamma)*(1+q)**(2*gamma-6)

        # Split at closest approach whenever it lies in front of the observer.
        kwargs = dict(epsabs=0, epsrel=1e-9, limit=200)
        if t0 < 0 < t1:
            return quad(integrand, t0, 0, **kwargs)[0] + quad(integrand, 0, t1, **kwargs)[0]
        return quad(integrand, t0, t1, **kwargs)[0]

    print(f"Computing NFW LOS profile: gamma={gamma:g}, rs={rs:g} kpc", flush=True)
    angles = np.geomspace(1e-10, 179.999, 2400)
    profile = np.array([los(angle) for angle in angles])
    interpolator = PchipInterpolator(np.log(angles), np.log(profile), extrapolate=False)

    def pixel_average(x, y, order):
        nodes, weights = leggauss(order)
        dx, dy = np.meshgrid(nodes/2, nodes/2)
        weights2 = np.outer(weights, weights).ravel()/4
        result = np.empty(len(x))
        chunk = max(1, 300000//(order*order))
        for start in range(0, len(x), chunk):
            stop = min(len(x), start+chunk)
            angle = separation(x[start:stop, None]+dx.ravel(),
                               y[start:stop, None]+dy.ravel())
            if np.any(angle <= angles[0]) or np.any(angle >= angles[-1]):
                raise ValueError("Pixel sample outside radial interpolation grid")
            result[start:stop] = np.exp(interpolator(np.log(angle))) @ weights2
        return result

    yy, xx = np.indices((n, n), dtype=float)
    radius = separation(xx.ravel(), yy.ravel())
    orders = np.full(n*n, 4, dtype=int)
    orders[radius < max(0.5, 4*binsz)] = 32
    orders[radius < 1.5*binsz] = cusp_order
    sky = np.empty(n*n)
    for order in np.unique(orders):
        chosen = orders == order
        print(f"Averaging {chosen.sum()} pixels with {order}x{order} quadrature", flush=True)
        sky[chosen] = pixel_average(xx.ravel()[chosen], yy.ravel()[chosen], int(order))
    sky = sky.reshape(n, n)
    solid_angle = np.deg2rad(binsz)**2  # Exact for CEA with lambda=1.
    if not np.all(np.isfinite(sky)) or np.any(sky <= 0):
        raise RuntimeError("Invalid template values")
    sky /= np.sum(sky, dtype=np.float64)*solid_angle

    header = w.to_header()
    header["BUNIT"] = ("sr-1", "Normalized intrinsic intensity")
    header["OBJECT"] = "GCE_NFW"
    header["PROFILE"] = "GNFW"
    header["LOSPOWER"] = (2, "Line-of-sight integral of density squared")
    header["GAMMA"] = float(gamma)
    header["RS_KPC"] = float(rs)
    header["RSUN_KPC"] = float(r_sun)
    header["RHALOKPC"] = float(r_halo)
    header["PIXSOLID"] = (solid_angle, "Pixel solid angle in sr")
    header["MAPWIDTH"] = (n*binsz, "Projection-plane width in degrees")
    header["CUSPORD"] = cusp_order
    header["MAPNORM"] = (float(sky.sum()*solid_angle), "Integral over map footprint")
    header.add_history("Pixel-averaged LOS integral of generalized NFW density squared.")
    header.add_history("Unit integral over footprint. No PSF, exposure, or spectrum.")
    header.add_history("Absolute density scale cancels; this is not a physical J-factor map.")
    output.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(sky.astype(np.float64), header).writeto(
        output, overwrite=overwrite, checksum=True)
    print(f"Wrote {output.resolve()} ({n}x{n}); integral = {sky.sum()*solid_angle:.12f}")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", default="nfw_gamma1.fits")
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--rs", type=float, default=20.0, help="Scale radius [kpc]")
    parser.add_argument("--r-sun", type=float, default=8.5, help="Solar radius [kpc]")
    parser.add_argument("--r-halo", type=float, default=250.0, help="Halo edge [kpc]")
    parser.add_argument("--width", type=float, default=60.0, help="Projection width [deg]")
    parser.add_argument("--binsz", type=float, default=0.1, help="Pixel scale [deg]")
    parser.add_argument("--cusp-order", type=int, default=256)
    parser.add_argument("--overwrite", action="store_true")
    make_nfw_template(**vars(parser.parse_args()))
