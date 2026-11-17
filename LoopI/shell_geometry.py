"""Uniform-volume Wolleben shell geometry and isotropic sky mask."""
import numpy as np
SHELLS = [(341.,3.,78.,62.,81.), (332.,37.,95.,58.,82.)]
EXAMPLE_MASK = dict(plane_abs_b_min=30., pole_abs_dec_max=90.,
                    bubbles_abs_l_max=45., bubbles_abs_b_max=90.)

def path_length(l,b,pars):
    """Forward-ray intersection with outer sphere minus inner sphere, in pc.
    Correctly includes observer inside S1's outer sphere (d=78 < Rout=81).
    Uniform volume emissivity, no magnetic-field weighting or extra masks.
    """
    l0,b0,d,ri,ro=pars
    l,b,l0,b0=map(np.deg2rad,(l,b,l0,b0))
    q=d*(np.sin(b)*np.sin(b0)+np.cos(b)*np.cos(b0)*np.cos(l-l0))
    impact2=np.maximum(d*d-q*q,0.)
    def inside(r):
        disc=r*r-impact2
        root=np.sqrt(np.maximum(disc,0.))
        return np.where(disc>=0,np.maximum(q+root,0)-np.maximum(q-root,0),0.)
    return np.maximum(inside(ro)-inside(ri),0.)

def basic_mask(l,b,dec,spec):
    expected=set(EXAMPLE_MASK)
    if set(spec)!=expected:
        raise ValueError('Mask has unsupported keys; supply --driver with its sky_mask, or --diagnostics.')
    l=(l+180)%360-180
    return ((abs(b)>=spec['plane_abs_b_min']) & (abs(dec)<=spec['pole_abs_dec_max']) &
            ~((abs(l)<=spec['bubbles_abs_l_max']) & (abs(b)<=spec['bubbles_abs_b_max'])))
