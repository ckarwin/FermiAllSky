# Isotropic extraction with global GALPROP normalizations — version 3

This version supersedes versions 1 and 2. It preserves each GALPROP component's
predicted spectrum and spatial distribution, fitting **one normalization per
component across the entire energy range**. All 4FGL parameters remain fixed.
The empirical isotropic component has one free amplitude per energy bin.
All amplitudes are fitted simultaneously in a masked Poisson likelihood.


## Matching your uploaded config and initializing isotropic emission

The shipped wrapper YAML preserves your FT1/spacecraft paths, tmin=239557417,
tmax=780439551, emin=300, emax=800000, binsperdec=8, evclass=128, evtype=3,
zmax=90, quality filter, IRFs and logging (chatter=0, verbosity=1).
It keeps the Galactic coordinate convention. The GC center no longer selects
an ROI because the preparation is full-sky.

The required spatial change is from a 40-degree WCS ROI to full-sky HEALPix.
The supplied hpx_order=8 corresponds to NSIDE=256, characteristic pixel scale
about 0.229 degrees, close to your 0.2-degree WCS bins. The preserved binsz=0.2
is only bookkeeping in HPX mode; hpx_order determines actual angular resolution.
All-sky source maps at this resolution can be large. Use order 6 for a first
smoke test in a separate outdir, then verify convergence with finer binning.
The output is a sibling Isotropic_Analysis_v3 directory, not your existing
ThreeML_Analysis directory.

Your original edisp_bins=-1 is changed to auto, resolving to -2 with your energy
range/binning. This covers the wider migration range discussed by the FSSC.
The isotropic component retains apply_edisp=false.

Default isotropic initialization:

```yaml
isotropic:
  mode: powerlaw
  index: -2.3
  pivot_MeV: 1000.0
  intensity_at_pivot: 1.0e-10
```

The reference intensity has units ph cm^-2 s^-1 sr^-1 MeV^-1. There is no external
measurement or prior in this mode. Every isotropic bin has an independent
nonnegative normalization. Rescaling the reference normalization just inversely
rescales those fitted coefficients, leaving physical fitted intensities unchanged
(up to numerical precision). The tests explicitly verify this invariance.
The assumed index controls shape WITHIN each finite energy bin; it does not tie
the amplitudes across bins. This within-bin assumption is needed to integrate
exposure/response over a bin. Repeat with other reference indices or finer bins
if the result is sensitive to it. A constant differential intensity is available
by setting index: 0.0, although a declining spectrum is a more practical basis
across the large energy range. The reference normalization should stay within a
reasonable numerical scale even though it has no physical constraint on the fit.

An optional external reference remains available:

```yaml
isotropic:
  mode: file
  file: /path/to/reference.txt
```

Independent Galactic controls (normalization need not equal 1):

```yaml
galprop:
  IC:   {file: /path/to/IC.fits,   free: true,  normalization: 1.0}
  pi0:  {file: /path/to/pi0.fits,  free: true,  normalization: 1.0}
  brem: {file: /path/to/brem.fits, free: false, normalization: 1.0}
```

Frozen contributions are added to the fixed counts prediction and omitted from
the optimizer/covariance. Their fixed values are recorded in the normalization
summary. Uncertainty in a frozen component is not propagated into the isotropic
errors. They still appear separately in the component plots.

## Energy dispersion and the meaning of the isotropic result

Energy dispersion is enabled in Fermipy/Fermitools for GALPROP and the fixed
4FGL catalog. It is deliberately disabled for the effective isotropic template,
following the FSSC convention for empirical isotropic background files:

    gtlike:
      edisp: true
      edisp_bins: <resolved setting>
      edisp_disable: [isodiff]

This is an effective, event-selection-dependent isotropic background spectrum,
including residual charged particles. It is NOT a deconvolved measurement of
pure astrophysical isotropic gamma rays. Reuse the extracted text file with
energy dispersion OFF for that component, while retaining it for GALPROP and
4FGL. Applying gamma-ray energy dispersion again to this effective template
would change the convention of the extracted model.

`edisp_bins: auto` selects enough bins for at least 0.2 dex of migration at
emin >=100 MeV, using the spectral correction (-2 or -3) when possible. For
emin <100 MeV it requires at least 0.4 dex and selects positive `edisp_bins`
(spatial-map mixing, which also handles the energy-dependent PSF). With the
provided 8 bins/decade, the settings are -2 at 300 MeV and +4 at 20 MeV.
FSSC documents that the positive mode is more computationally expensive and
its recommendation page describes it as still being validated. This script's
automatic low-energy choice is a modeling choice, not a certified 20 MeV
analysis prescription. You can also request a positive setting at higher energies.
Requests requiring more than the supported four migration bins fail explicitly;
reduce the number of analysis bins per decade if necessary.

GALPROP MapCubes must extend below and above the selected reconstructed-energy
range. The preflight conservatively requires abs(edisp_bins)+1 padded bins on
each side. No silent extrapolation of GALPROP is allowed. Positive edisp requires
fresh exposure/source maps with the extra energy coverage: old version-1 caches
cannot be used. The on-machine checks compare explicit Poisson likelihood
changes with the native energy-dispersion-enabled likelihood.

## Energy range

Default: **300 MeV to 800 GeV**, matching your uploaded configuration. Set `emin: 20.0` for a 20 MeV lower threshold.
The code accepts emin >=20 MeV. This is software support, not evidence that a
particular event selection/IRF, source catalog or foreground model gives an
accurate isotropic measurement at 20 MeV. The fit needs valid response and model
coverage across the padded true-energy range, and FT1 data that actually contain
the selected low-energy events. In particular, an event file already selected
above 100/300 MeV cannot be used for a 20 MeV fit.

Below 100 MeV, investigate broad-PSF spillover from the plane/bubbles, Earth-limb
contamination, catalog spectral extrapolation and diffuse-model uncertainty.
The code uses a uniform zmax=90 degrees by default; choose event/zenith cuts
appropriate to your low-energy analysis and repeat with alternative masks/cuts.
These systematic effects are not included in the statistical error bars.

If the isotropic seed file does not reach 20 MeV, provide a wider seed or set
`allow_seed_extrapolation: true`. That option extends only the seed with an
endpoint power law; the amplitude of every isotropic bin is still fitted.
It does not fabricate a measured low-energy spectrum. The exact seed used is
saved as `isotropic_seed_used.txt`. The within-bin seed shape remains an
assumption, including wherever it is extrapolated.

## Mask

Excluded at pixel centers:

- |Galactic b| <10 degrees;
- |ICRS Dec| >60 degrees;
- |wrapped Galactic l| <=30 degrees AND |b| <=60 degrees.

The plane and pole cuts follow the FSSC background-model page. The bubble
rectangle is a deliberately broad approximation, not an official bubble
contour. All boundaries are configurable. The mask acts on observed count
pixels only. Full-sky Galactic and catalog emission is convolved first, so
modeled flux from excluded directions can still contribute through the PSF.
This does not remove PSF leakage from any UNMODELED bubble/plane component.

## Inputs and environment

Run inside your working Fermipy/Fermitools environment with CALDB, NumPy,
SciPy, PyYAML, Astropy, Matplotlib, and Fermipy's compatible Gammapy/HEALPix
requirements. The plotting code itself does not require healpy. Avoid replacing
packages in an otherwise working LAT environment just for this script.

Edit `config.yaml`:

1. Supply genuinely all-sky FT1 data, or a text list of absolute FT1 paths,
   plus spacecraft data covering the observation. A GC-only download is invalid.
2. Set tmin/tmax in Fermi MET; the example dates are placeholders for your data.
3. Supply a local 4FGL FITS catalog and matching extended-source template folder.
4. Supply full-sky Fermitools WCS GALPROP FITS MapCubes with an ENERGIES extension
   in MeV. Intensities must be ph cm^-2 s^-1 sr^-1 MeV^-1, NOT E^2 dN/dE and
   NOT normalized to unit integral. Raw GALPROP HEALPix outputs need conversion.
   Set `galprop_intensity_units_confirmed: true` once checked.
5. The event class/type, IRFs, time range, data paths and logging now match your
   uploaded configuration. Catalog, extended-template and separate GALPROP paths
   were absent in that file and remain explicit /EDIT placeholders.
6. For each GALPROP entry, set `file`, `free` and `normalization`. Defaults:
   IC free, pi0 free, brem fixed at 1.0. For a free component, normalization is
   the optimizer's initial value; for a frozen component it is the fixed value.
   Each free component has just ONE fitted normalization over all energies.
   These per-component settings supersede the former galprop_mode switch.
7. No external isotropic file is needed in the default powerlaw mode. The
   generated reference is only a count-template basis; it is not a prior.

Do not add the official Galactic diffuse model alongside these GALPROP components.
That would double-count the foreground. YAML paths are resolved relative to the
YAML; paths inside FT1 list files should be absolute.

## Run and plot

After extracting the ZIP, work inside its `isotropic_fit` folder:

```bash
python fit_isotropic.py self-test
python fit_isotropic.py run config.yaml
```

`run` prepares the templates, fits, and plots. Alternatively:

```bash
python fit_isotropic.py prepare config.yaml
python fit_isotropic.py fit config.yaml
python fit_isotropic.py plot config.yaml
```

Replot from a notebook in that folder:

```python
from fit_isotropic import read_config, plot, result_dir
c = read_config('config.yaml')
plot(c)
print(result_dir(c))
```

Or plot a specific existing fit directory directly:

```python
from plots import plot_results
plot_results('/path/to/isotropic_run_v2/fit_...', {
    'bands_MeV': [[20, 100], [100, 300], [300, 1000], [1000, 10000]],
    'map_nlon': 360,
})
```

Plots use whole energy bins. Requested band boundaries are snapped by selecting
bins whose centers fall in the requested band; actual edge coverage is displayed
in plot titles and recorded in `plot_bands.txt`. There is no fractional splitting
of observed energy bins.

## Joint likelihood

For energy bin k and retained pixel p:

    mu_kp = fixed_4FGL_kp + a_iso,k*T_iso,kp + sum_j A_j*T_GALPROP,jkp
    -log L = sum_k,p [mu_kp - n_kp log(mu_kp)] + constant.

A_j is a single scalar shared by ALL energies. The Galactic templates already
include the requested energy-dispersion response; linearity permits scaling
them globally without recalculating the response each optimizer iteration.
The effective isotropic template has no energy redistribution, so its amplitudes
can be varied independently in reconstructed-energy bins. Catalog sources are
never optimized. `setup(init_sources=False)` skips expensive source diagnostics.

Isotropic intensity = fitted a_iso,k times the seed intensity at the geometric
bin center. Do not multiply/divide by 4*pi or the retained solid angle.
Profile intervals reoptimize every OTHER isotropic amplitude and every global
Galactic amplitude. Shared Galactic nuisance parameters create correlations
between isotropic bins. Their full covariance is saved.

## Outputs and plots

Preparation produces `templates.npz`, `preparation_checks.json`,
`fermipy_config.yaml`, and `isotropic_seed_used.txt`.
Each `fit_<hash>/` directory contains:

- `isotropic_bins.csv`: intensities, 68% profile interval endpoints, 95% upper
  limits, TS, counts and deviance per bin.
- `galprop_normalizations.json`: ONE normalization and profile interval per
  Galactic component (or the specified value for frozen components).
- `covariance.npz`: full joint amplitude covariance, best fit and parameter
  order (all isotropic bins followed by FREE global Galactic normalizations).
- `diagnostics.npz`: observed/model count cubes, separate fitted components,
  pixel coordinates, retained mask and energy edges. Predictions in excluded
  regions are extrapolations; they never contribute to reported fit statistics.
- `fit_report.json`: settings, selected area, assumptions and joint deviance.
- `isotropic_spectrum.png` and `.pdf`: E^2 I(E) in GeV cm^-2 s^-1 sr^-1,
  with profile error bars and upper limits for isotropic TS<4.
- `counts_spectrum.png`: masked data, fitted total, fixed catalog, isotropic and
  GALPROP counts; fractional residuals below. Data errors are sqrt(counts), not
  uncertainty in the fitted model prediction.
- `residuals_<band>.png`: Galactic Mollweide maps of data−model, (data−model)/model,
  and (data−model)/sqrt(model). Maps are unsmoothed; excluded pixels are gray.
  The Pearson residual is a diagnostic, NOT a likelihood significance/TS map.
  Colors clip at the 99th percentile of absolute retained-pixel residuals.
  The display grid samples nearest original pixels; it is not a new count binning.

If all bins have TS>=4, positive intensity and a nonsingular joint covariance,
`isotropic_spectrum.txt` is exported (MeV, differential intensity per MeV per sr).
Otherwise use the CSV/upper limits. The text file log-log interpolates bin-center
estimates, which approximates the piecewise seed-scaled extracted spectrum.
Its valid range is first to last bin CENTER. Extract beyond the intended reuse
interval, then reuse its interior. Reforward-fold the interpolated text spectrum
before precision use. Keep energy dispersion disabled for it on reuse.

## Reuse and validation limits

Changing the mask or per-component free/normalization controls reuses the prepared response products and
creates a different fit directory. Plot options can be changed without refitting.
A fit directory is never overwritten. Choose a NEW outdir if preparation fails,
if input file contents change, or if energy/event/response settings change.
The configuration signature includes the pipeline version, but does not hash
all large science input files. Preparation of thousands of all-sky source maps
can be costly. hpx_order=6 can be used as a smoke test in a separate outdir;
check finer angular binning before drawing scientific conclusions.

Local tests cover global-normalization recovery, a toy energy-redistribution
response, analytic gradients/Hessian, coupled profiles, Poisson recovery, mask
exclusion, covariance/boundaries, seed extrapolation, 20/300 MeV configuration,
end-to-end cached fitting/export, and generation of spectrum/residual plots.
See VALIDATION.txt for what actually ran. The tests are not LAT calibration tests.

On your machine, preparation checks catalog inclusion, frozen parameters,
full-sky dimensions, energy coverage, model-map sum closure, energy-dispersion
activation, and native-vs-explicit likelihood changes under normalization
perturbations. Real LAT response generation has NOT been run in the authoring
environment; Fermipy/Fermitools and your science data were unavailable.

## Sources

- https://fermipy.readthedocs.io/en/latest/_modules/fermipy/gtanalysis.html
- https://fermipy.readthedocs.io/en/latest/_modules/fermipy/roi_model.html
- https://fermi.gsfc.nasa.gov/ssc/data/access/lat/BackgroundModels.html
- https://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Pass8_edisp_usage.html
