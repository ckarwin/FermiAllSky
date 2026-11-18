# Fermi-bubble spatial templates

Spatial templates of the Fermi bubbles for the inner-Galaxy fit in
[`Main/`](../Main/), split into a low-latitude and a high-latitude part so that
each can have its own spectrum.

| File | Part | Description |
|---|---|---|
| `bubbles_lowlat_ack2017_smooth.fits` | \|b\| ≲ 10° | Smoothed morphology (used by `Main/run_3ML_master.py`) |
| `bubbles_highlat_ack2017_smooth.fits` | \|b\| ≳ 10° | Smoothed morphology (used by `Main/run_3ML_master.py`) |
| `bubbles_lowlat_ack2017_approx.fits` | \|b\| ≲ 10° | Unsmoothed version |
| `bubbles_highlat_ack2017_approx.fits` | \|b\| ≳ 10° | Unsmoothed version |

## How they were made

The maps are derived from the bubble template in Figure 8 of
Ackermann et al. (2017), arXiv:1704.03910. The published counts map was divided
by the exposure at 1 GeV to approximate an intensity map. This is an
approximate reconstruction, not an official intensity release.

- **Smoothed maps**: the bubble outline was closed with a 3° morphological
  closing and its holes filled; zeros inside the outline are treated as
  missing; the map is then convolved with a 2° Gaussian (a morphology choice,
  not the LAT PSF). Low and high latitudes are separated by a smooth
  (cosine) transition between |b| = 9° and 11°.
- **Approximate maps**: published zeros are kept, no smoothing; the split is a
  sharp cut at |b| = 10°.

The provenance is also recorded in each file's FITS header (`HISTORY` cards).

## Conventions

- 400 × 400 pixels of 0.2° (80° × 80°), ICRS cylindrical equal-area projection
  centred on the Galactic centre.
- Units: sr⁻¹. Each map integrates to 1 over its footprint (header `SPATNORM`),
  so the attached spectral model gives the flux of the bubble part **inside
  the map footprint**. `ROIFRAC` in the header records the fraction of the
  template that falls inside the inner-Galaxy region (about 0.92–0.93 for
  the low-latitude and 0.33 for the high-latitude maps).
- No PSF deconvolution.

In 3ML/astromodels, load each map with
`SpatialTemplate_2D(fits_file=..., frame="icrs")`, fix `K` to 1 and give each
part its own spectrum (Main offers power law, log parabola or cutoff power law).

To use a bubble template in the isotropic fit instead, it must first be turned
into a 3D mapcube with a physical spectrum (see the
[isotropic guide](../Isotropic/README.md#additional-diffuse-components-sun-moon-others)).
The default isotropic fit masks the bubble region instead.
