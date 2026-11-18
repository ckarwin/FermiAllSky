"""Geometry of the two Wolleben (2007) Loop I shells, S1 and S2.

Each shell is the region between two concentric spheres. ``SHELLS`` lists
(centre l [deg], centre b [deg], distance [pc], inner radius [pc],
outer radius [pc]) for S1 and S2.
"""

import numpy as np

# (l0, b0, distance, r_inner, r_outer) for S1 and S2
SHELLS = [(341.0, 3.0, 78.0, 62.0, 81.0), (332.0, 37.0, 95.0, 58.0, 82.0)]
EXAMPLE_MASK = dict(
    plane_abs_b_min=30.0, pole_abs_dec_max=90.0, bubbles_abs_l_max=45.0, bubbles_abs_b_max=90.0
)


def path_length(l, b, pars):
    """Forward-ray intersection with outer sphere minus inner sphere, in pc.
    Correctly includes observer inside S1's outer sphere (d=78 < Rout=81).
    Uniform volume emissivity, no magnetic-field weighting or extra masks.
    """
    l0, b0, d, ri, ro = pars
    l, b, l0, b0 = map(np.deg2rad, (l, b, l0, b0))
    q = d * (np.sin(b) * np.sin(b0) + np.cos(b) * np.cos(b0) * np.cos(l - l0))
    impact2 = np.maximum(d * d - q * q, 0.0)

    def inside(r):
        disc = r * r - impact2
        root = np.sqrt(np.maximum(disc, 0.0))
        return np.where(disc >= 0, np.maximum(q + root, 0) - np.maximum(q - root, 0), 0.0)

    return np.maximum(inside(ro) - inside(ri), 0.0)


def basic_mask(l, b, dec, spec):
    """Same sky mask as ``Isotropic/fit_isotropic.sky_mask`` (True = kept)."""
    expected = set(EXAMPLE_MASK)
    if set(spec) != expected:
        raise ValueError(
            "Mask has unsupported keys; supply --driver with its sky_mask, or --diagnostics."
        )
    l = (l + 180) % 360 - 180
    return (
        (abs(b) >= spec["plane_abs_b_min"])
        & (abs(dec) <= spec["pole_abs_dec_max"])
        & ~((abs(l) <= spec["bubbles_abs_l_max"]) & (abs(b) <= spec["bubbles_abs_b_max"]))
    )
