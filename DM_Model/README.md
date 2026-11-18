# NFW (dark-matter-like) spatial templates

Spatial templates for the Galactic-centre excess component of the inner-Galaxy
fit in [`Main/`](../Main/). Each map is the line-of-sight integral of the
**square** of a generalized NFW density profile (the morphology expected from
dark-matter annihilation), normalized to a unit integral over the map.

| File | γ (inner slope) | Map |
|---|---|---|
| `nfw_gamma1.fits` | 1.0 (used by `Main/run_3ML_master.py`) | 60° × 60°, 0.1° pixels |
| `nfw_gamma1p2.fits` | 1.2 | 60° × 60°, 0.1° pixels |

Both use r_s = 20 kpc, R_sun = 8.5 kpc and a halo edge at 250 kpc.

## Conventions

- Density: ρ(r) ∝ (r/r_s)^−γ (1 + r/r_s)^(γ−3). Intensity: J(ψ) = ∫ ρ² ds along
  the line of sight from the Sun to the halo edge.
- Units: sr⁻¹. The map integrates to 1 over its footprint (header `MAPNORM`), so
  the spectral model attached to the source gives the flux **inside the map
  footprint**, not of the whole halo.
- Projection: equatorial (ICRS) cylindrical equal-area (CEA) centred on the
  Galactic centre, so every pixel has the same solid angle.
- No exposure, PSF, particle mass, cross section or absolute J-factor is
  included.

In 3ML/astromodels, load the map with
`SpatialTemplate_2D(fits_file=..., frame="icrs")`, fix its `K` to 1 and fit the
amplitude of the spectral model (an exponentially cutoff power law in Main).

## Making a new template

```bash
python make_nfw_template.py --output nfw_gamma1.fits
python make_nfw_template.py --gamma 1.2 --output nfw_gamma1p2.fits
```

| Option | Default | Meaning |
|---|---|---|
| `--gamma` | 1.0 | Inner slope; must be < 1.5 |
| `--rs` | 20 | Scale radius [kpc] |
| `--r-sun` | 8.5 | Sun–Galactic-centre distance [kpc] |
| `--r-halo` | 250 | Halo edge [kpc] |
| `--width` | 60 | Map width [deg]; keep it larger than the fit region plus a PSF margin |
| `--binsz` | 0.1 | Pixel size [deg] |
| `--cusp-order` | 256 | Quadrature order for the central pixels |
| `--overwrite` | off | Replace an existing file |

The central cusp is averaged over each pixel with high-order Gauss–Legendre
quadrature (no artificial core). To check convergence, double `--cusp-order`
or halve `--binsz` and compare. Needs NumPy, SciPy and Astropy.
