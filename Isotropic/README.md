# All-sky isotropic fit

This folder measures the **effective isotropic gamma-ray background spectrum**
with Fermi-LAT data. It fits the whole high-latitude sky with a model made of
Galactic diffuse emission, catalog sources, the Sun, the Moon, Loop I and an
isotropic component, and reports the isotropic intensity in every energy bin.
The resulting spectrum (`isotropic_spectrum.txt`) can then be used as a fixed
background in other analyses, such as the inner-Galaxy fit in [`Main/`](../Main/).

- [How the fit works](#how-the-fit-works)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Commands](#commands)
- [Configuration reference](#configuration-reference)
- [What to rerun after a change](#what-to-rerun-after-a-change)
- [Output files](#output-files)
- [Using the isotropic spectrum elsewhere](#using-the-isotropic-spectrum-elsewhere)
- [Conventions and caveats](#conventions-and-caveats)
- [Testing](#testing)
- [Files in this folder](#files-in-this-folder)

## How the fit works

The analysis has two stages.

**1. Prepare (slow, run once).** Fermipy builds a full-sky HEALPix counts cube
and, for every model component, the counts that component predicts at unit
normalization (a *template*). The templates include the LAT exposure, PSF and
energy dispersion, so all instrument physics is handled by Fermitools.
Everything is saved to `templates.npz`.

**2. Fit (fast, repeat as often as you like).** The templates are combined in a
binned Poisson likelihood, using only pixels that pass the sky mask. For energy
bin *k* and pixel *p* the model is

```
counts[k, p] = catalog[k, p]                    (all 4FGL sources, fixed)
             + a_k  * isotropic[k, p]           (one free amplitude per energy bin)
             + sum_j A_j * diffuse_j[k, p]      (one amplitude per component, same at all energies)
```

- **Catalog sources** (4FGL point and extended sources) stay at catalog values.
- **Isotropic**: each energy bin has its own free amplitude `a_k`, so the
  isotropic spectrum is measured bin by bin. The reported intensity is
  `a_k × I_ref(E_k)`, where `I_ref` is a reference power law and `E_k` is the
  bin's geometric centre.
- **Diffuse components** (GALPROP π⁰, IC, bremsstrahlung, Sun, Moon, Loop I, and
  anything else you add) each have one normalization `A_j` shared by all
  energies. Each can be free or fixed.

All bins are fitted jointly because the diffuse normalizations couple them.
Amplitudes are constrained to be ≥ 0.

Templates are always computed for the **full sky** and the mask is applied
afterwards. Counts that leak into the retained region from bright masked areas
(through the PSF) are therefore modelled correctly.

## Requirements

| Task | Needs |
|---|---|
| `prepare`, `run` | Fermipy and Fermitools (pyLikelihood) with the LAT IRFs; event (FT1) and spacecraft (FT2) files; the 4FGL catalog and its extended-source templates; GALPROP mapcubes; Sun/Moon templates; Loop I templates |
| `fit`, `plot` | NumPy, SciPy, PyYAML, Matplotlib (Fermipy only if an analytic Loop I shape is fitted) |
| `self-test` | NumPy, SciPy, PyYAML, Matplotlib |
| Loop I template build | NumPy, SciPy, Astropy, Matplotlib, PyYAML |

The Sun and Moon templates are made with [`Sun_and_Moon/`](../Sun_and_Moon/).
The Loop I templates are made with [`LoopI/`](../LoopI/).

## Quick start

All commands are run from this folder in an activated Fermipy environment.

```bash
cd Isotropic

# 0. (Optional) check the installation with synthetic data: 16 tests should pass.
python fit_isotropic.py self-test

# 1. Edit config.yaml: data, catalog, GALPROP, Sun/Moon paths, outdir.
#    For a first preparation, delete the reuse_from line (or set it to null).

# 2. Build the Loop I templates once (they are not stored in git).
python ../LoopI/build_templates.py --config config.yaml

# 3. Prepare, fit and plot.
python fit_isotropic.py run config.yaml
```

On a cluster, submit [`slurm.sh`](slurm.sh) instead of step 3 (`sbatch slurm.sh`).
The full-sky preparation at `hpx_order: 8` is memory-hungry; the example job
requests 220 GB.

When the run finishes, the results are in `<outdir>/bestfit_<hash>/`. Start with
`isotropic_spectrum.png`, `diffuse_normalizations.json` and `isotropic_bins.csv`.

## Commands

```
python fit_isotropic.py <command> [config.yaml]
```

| Command | What it does | When to use it |
|---|---|---|
| `prepare config.yaml` | Checks inputs, builds the per-shell Loop I cubes if needed, runs Fermipy setup, extracts the templates, checks them against Fermipy's own likelihood, and writes `templates.npz`. | First run, or after changing data, response or model files (in a **new** `outdir`). |
| `fit config.yaml` | Loads `templates.npz`, applies the mask, runs the fit and writes the results. | After `prepare`, or after changing only the mask or free/fixed/normalization settings. |
| `plot config.yaml` | Remakes the plots from the saved fit selected by the config. | After a fit, or after changing plot settings. |
| `run config.yaml` | `prepare`, then `fit`, then `plot`. | A complete run into a fresh `outdir`. It does not resume a partly finished run. |
| `self-test` | 16 synthetic tests of the fit, profiles, mask, export and plotting. | Checking an installation. It does not test LAT data or the LAT response. |

Two safety rules prevent results from being silently overwritten:

- `prepare` stops if `<outdir>/fermipy/` already exists (also after a failed
  preparation). Use a new `outdir`, or delete the failed one.
- `fit` stops if its result folder already exists. Use `plot` to remake plots,
  or move the old result folder aside to refit.

`fit` also checks that `templates.npz` was prepared from the current inputs. If
any input file, path or data/response setting changed, it asks you to prepare
again in a new `outdir`.

## Configuration reference

[`config.yaml`](config.yaml) is the only file you normally edit. It is a wrapper
configuration; the code writes the actual Fermipy configuration for you
(`<outdir>/fermipy_config.yaml`).

Relative paths are resolved from the folder that contains the YAML file, not
from your shell's working directory. `~` and environment variables (`$HOME`)
are allowed. Each key may appear only once; a duplicated key is reported as an
error with its line number.

### Data, selection and response

| Key | Meaning | Example |
|---|---|---|
| `evfile` | FT1 event file, or a text file listing FT1 files | `filelist.txt` |
| `scfile` | FT2 spacecraft file | `lat_spacecraft_merged.fits` |
| `catalog` | Source catalog FITS file (all sources are included and fixed) | `gll_psc_v35.fit` (4FGL-DR4) |
| `extdir` | Folder with the catalog's extended-source templates | `Extended_14years/Templates/` |
| `evclass`, `evtype`, `irfs` | Event class, event type and matching IRFs | `128`, `3`, `P8R3_SOURCE_V3` |
| `tmin`, `tmax` | Time range in MET seconds | `239557417`, `780439551` |
| `emin`, `emax` | Energy range in MeV | `300`, `100000` |
| `binsperdec` | Energy bins per decade | `8` |
| `hpx_order` | HEALPix order of the sky grid; NSIDE = 2^order | `8` (NSIDE 256, ~0.23° pixels) |
| `binsz` | Kept for Fermipy bookkeeping; does **not** set the HEALPix resolution | `0.2` |
| `zmax` | Zenith-angle cut in degrees | `90` |
| `filter` | gtmktime filter | `DATA_QUAL>0 && LAT_CONFIG==1` |
| `edisp_bins` | Extra true-energy bins for energy dispersion: `auto`, `-2`, `-3`, `2`, `3` or `4` | `auto` (gives -2 here) |
| `logging` | Fermipy logging settings | `{chatter: 0, verbosity: 1}` |

About `edisp_bins`: negative values use the standard Fermitools method (as in
4FGL) and are used for `emin ≥ 100 MeV`. Below 100 MeV the code requires
positive values; `auto` then selects `+4`. The extra bins must span at least
0.2 decade (0.4 decade below 100 MeV).

For a quick, cheap test of the whole chain, use `hpx_order: 6` and a separate
`outdir`. Check that the results are stable at the full resolution before using
them for science.

### GALPROP components

```yaml
galprop:
  IC:   {file: /path/ics_mapcube.fits,        free: true,  normalization: 1.0}
  pi0:  {file: /path/pion_decay_mapcube.fits, free: true,  normalization: 1.0}
  brem: {file: /path/bremss_mapcube.fits,     free: false, normalization: 1.0}
```

- Any names can be used, except `isodiff`. At least one component is required.
- `free: true` fits one normalization for all energies, starting at
  `normalization`. `free: false` multiplies the cube by `normalization`.
- The files must be 3D mapcubes in differential intensity
  (ph cm⁻² s⁻¹ sr⁻¹ MeV⁻¹) with an `ENERGIES` table in MeV.

**Energy padding.** Energy dispersion needs the model slightly outside
`[emin, emax]`. If a GALPROP cube stops short, `prepare` writes a padded copy
(the original file is never changed) in which each pixel is extended with a
power law fitted to its last energy planes:

```yaml
galprop_extrapolation:
  high: true          # allow extension above the cube
  low: false          # allow extension below the cube
  fit_points: 3       # boundary planes used for each pixel's slope
  index_min: -6.0     # clip fitted exponents (dN/dE ∝ E^index) to this range
  index_max: 2.0
  index_shift: 0.0    # add to every exponent (for systematic tests)
  max_factor: 3.0     # never extrapolate by more than this factor in energy
```

If padding is needed on a side that is not enabled, or by more than
`max_factor`, `prepare` stops with a message. What was done is recorded in
`<outdir>/galprop_extrapolation.json`.

### Additional diffuse components (Sun, Moon, others)

```yaml
additional_diffuse:
  Sun:  {file: /path/sun_template.fits,  free: true, normalization: 1.0}
  Moon: {file: /path/moon_template.fits, free: true, normalization: 1.0}
```

These behave like GALPROP components (`free` defaults to `true` here). Each file
must be a 3D differential-intensity mapcube with a celestial WCS and an
`ENERGIES` table in MeV. Unlike GALPROP cubes these are **not** padded: they
must already cover the energy-dispersion margin (about 127 MeV to 237 GeV for
the default 300 MeV–100 GeV, `edisp_bins: -2` setup). The Sun/Moon generator in
[`Sun_and_Moon/`](../Sun_and_Moon/) produces files with this coverage.

You can add any other component in the same way, for example a Fermi-bubble
mapcube. The 2D maps in [`Bubbles_Model/`](../Bubbles_Model/) cannot be used
directly: they must first be multiplied by a spectrum and written as a 3D
mapcube. Do not add `LoopI`, `LoopI_S1` or `LoopI_S2` here; use the `loopI`
block instead.

### Loop I

```yaml
loopI:
  enabled: true
  mode: separate            # combined | separate
  spectrum: empirical       # empirical | powerlaw | cutoff_powerlaw | logparabola
  free_norm: true
  normalization: 0.14
  # directory: ../LoopI     # default: the repository-level LoopI/ folder
```

| Option | Meaning |
|---|---|
| `enabled: false` | Leave Loop I out of the model. |
| `mode: combined` | One source, `LoopI`: one amplitude multiplies the S1 + S2 sum. This is also used when `mode` is omitted. |
| `mode: separate` | Two sources, `LoopI_S1` and `LoopI_S2`, with independent amplitudes and the same spectral shape. With both amplitudes equal, the sum equals the combined model. |
| `free_norm` | `true`: fit the amplitude (one per shell in separate mode). `false`: hold it at `normalization`. |
| `normalization` | Starting value if free, fixed value otherwise. |
| `shells` | Separate mode only: per-shell overrides, e.g. `shells: {S1: {free_norm: true, normalization: 0.14}, S2: {free_norm: false, normalization: 0.2}}`. |
| `directory` | Folder with the Loop I templates, if not `../LoopI`. |
| `spectrum` | `empirical` uses the fixed reconstructed spectrum (no `parameters`). The analytic options let the spectral shape be fitted; see below. |

The templates are built once with `python ../LoopI/build_templates.py --config
config.yaml` (see [`LoopI/README.md`](../LoopI/README.md) for the model). In
separate mode, `prepare` then creates `loopI_S1_*_mapcube.fits` and
`loopI_S2_*_mapcube.fits` automatically.

**Fitting the Loop I spectral shape.** Choose an analytic spectrum and give its
parameters:

```yaml
loopI:
  enabled: true
  mode: separate
  spectrum: cutoff_powerlaw
  free_norm: true
  normalization: 0.14
  parameters:
    index:      {value: 2.4,   free: true, min: 0.5,  max: 5.0}
    cutoff_MeV: {value: 30000, free: true, min: 1000, max: 10000000}
compute_profiles: false
```

| Spectrum | Shape (E₀ = 1000 MeV) | Parameters |
|---|---|---|
| `powerlaw` | (E/E₀)^−index | `index` |
| `cutoff_powerlaw` | (E/E₀)^−index × exp(−(E−E₀)/cutoff) | `index`, `cutoff_MeV` |
| `logparabola` | (E/E₀)^−(alpha + beta ln(E/E₀)) | `alpha`, `beta` |

Indices are positive for falling spectra (`index: 2.4` means E^−2.4). Each
parameter accepts `value`, `free`, `min` and `max`; parameters you leave out get
defaults and are **free** unless you set `free: false`. S1 and S2 always share
the shape; their amplitudes can differ. During `fit`, every trial shape is
passed through the full LAT response with Fermitools, so this mode needs
Fermipy even for `fit`. Two restrictions apply:

- A free shape requires `compute_profiles: false`.
- The saved covariance of the amplitudes is conditional on the best-fit shape;
  there are no uncertainties on the shape parameters.

### Isotropic reference spectrum

```yaml
isotropic:
  mode: powerlaw
  index: -2.0
  pivot_MeV: 1000.0
  intensity_at_pivot: 1.0e-10   # ph cm^-2 s^-1 sr^-1 MeV^-1
allow_seed_extrapolation: false
```

The reference spectrum is built internally and saved as
`<outdir>/isotropic_seed_used.txt`. Its normalization does not matter (each bin
has a free amplitude); `index` only sets the spectral slope *within* each bin.
Alternatively, `isotropic: {mode: file, file: reference.txt}` reads a
two-column table (MeV, intensity). A table that does not cover
`[emin, emax]` is rejected unless `allow_seed_extrapolation: true`, which
extends it with endpoint power laws.

### Sky mask

```yaml
mask:
  plane_abs_b_min: 30.0     # keep |b| >= 30°
  pole_abs_dec_max: 90.0    # keep |dec| <= 90° (90 = no pole cut)
  bubbles_abs_l_max: 45.0   # remove |l| <= 45° ...
  bubbles_abs_b_max: 90.0   # ... where also |b| <= 90°
```

A pixel is kept if its centre passes all cuts. The example keeps |b| ≥ 30° and
removes the 90°-wide longitude band around the Galactic centre at all
latitudes, which excludes the Fermi bubbles and the 40° × 40° inner-Galaxy
region. The mask only affects `fit`; changing it does not require a new
preparation.

### Fit and plot options

| Key | Meaning |
|---|---|
| `compute_profiles` | `false` (default): best fit and Hessian covariance only. `true`: also profile-likelihood 68% intervals, 95% upper limits and TS for every isotropic bin and free normalization (slower; results go to `fit_<hash>/`). |
| `plots.bands_MeV` | Energy bands for residual maps, e.g. `[[300, 1000], [1000, 10000]]`. Each band uses the whole energy bins whose centres fall in it. |
| `plots.map_nlon` | Number of longitude cells of the residual-map display grid. |

### Output folder and reusing an earlier preparation

```yaml
outdir: output_with_solar_and_loopI_separate
reuse_from: output_with_solar_and_loopI     # omit or null for a first preparation
```

- `outdir` is where this run writes everything. Use a new `outdir` for every
  preparation.
- `reuse_from` points to the **root** of an older output (the folder that
  contains `fermipy_config.yaml` and `fermipy/`), never to its `fermipy/` or
  `bestfit_*` subfolder. It must differ from `outdir`.

With `reuse_from`, the data, selection, binning and IRF settings must be
identical to the old run; otherwise `prepare` stops and you need a fresh
preparation without `reuse_from`. When they match, the code copies the selected
events, counts cube, livetime cube and exposure map. It then compares every
source in the new model with the old one, including the contents of every file
the source refers to. Unchanged sources reuse their old source maps; new or
changed sources (for example a new GALPROP cube) are computed by Fermitools.
`<outdir>/reuse_report.json` lists which sources were reused and why the others
were regenerated. Old files are never modified.

## What to rerun after a change

| You changed... | Do this |
|---|---|
| The mask, or a free/fixed flag or normalization of a diffuse component | `fit`, then `plot` (same `outdir`) |
| Plot settings | `plot` |
| Analytic Loop I parameters, or one analytic spectrum for another | `fit`, then `plot` |
| Empirical ↔ analytic Loop I; combined ↔ separate; adding or removing a component | New `outdir` + `reuse_from` the old one, then `run` |
| A model file or its path (GALPROP, Sun, Moon, Loop I assets), the isotropic reference, or `galprop_extrapolation` | New `outdir` + `reuse_from` the old one, then `run` |
| Data, time or energy selection, binning, HEALPix order, event class/type, IRFs, `edisp_bins` | New `outdir` **without** `reuse_from`, then `run` |

The preparation records the size and modification time of every model cube, so
even overwriting a file under the same name is detected. Treat input files as
read-only and give changed models new file names.

## Output files

### Written by `prepare` (in `outdir`)

| File | Contents |
|---|---|
| `templates.npz` | Counts cube, one template per diffuse component, summed catalog prediction, energy edges, pixel coordinates (l, b, dec), reference spectrum and preparation metadata. This is everything `fit` needs. |
| `fermipy_config.yaml` | The Fermipy configuration that was generated and used. |
| `fermipy/` | Fermipy working folder: selected events, counts cube, livetime cube, exposure, `srcmdl_00.xml` (model definition) and `srcmap_00.fits` (+ numbered chunks: source maps for **all** catalog and diffuse sources). |
| `preparation_checks.json` | Template-versus-Fermipy likelihood checks, software versions, full configuration. |
| `source_fingerprints.json` | Fingerprints of every source, used by later `reuse_from` runs. |
| `isotropic_seed_used.txt` | The isotropic reference spectrum. |
| `galprop_extrapolation.json`, `padded_galprop/` | Energy-padding report and padded GALPROP copies (if needed). |
| `reuse_report.json` | Only with `reuse_from`: reused and regenerated sources. |

### Written by `fit` (in `outdir/bestfit_<hash>/` or `outdir/fit_<hash>/`)

The hash identifies the fit settings (mask, free/fixed flags, normalizations,
Loop I block), so different fit setups of one preparation sit side by side.

| File | Contents |
|---|---|
| `isotropic_spectrum.txt` | Energy [MeV] and isotropic intensity [ph cm⁻² s⁻¹ sr⁻¹ MeV⁻¹] at each bin centre. Written only if every bin is positive (with profiles also TS ≥ 4 in every bin and a valid covariance). |
| `isotropic_bins.csv` | Per bin: edges, centre, intensity, profile interval and 95% upper limit (NaN without profiles), amplitude, TS, counts, predicted counts, deviance, boundary flag, covariance flag. |
| `diffuse_normalizations.json` | Normalization of every diffuse component and whether it was fixed (plus profile results if requested). Also split into `galprop_normalizations.json` and `additional_diffuse_normalizations.json`. |
| `covariance.npz` | Covariance of all free parameters (`names` gives the order: `isodiff_bin000`, ..., then free diffuse components). |
| `best_fit_parameters.npz` | Best-fit parameter vector and minimum. |
| `fit_report.json` | Configuration, retained pixels and solid angle, total deviance, normalizations and the error method. |
| `diagnostics.npz` | Full-sky counts and best-fit model per component, plus the mask. |
| `loopI_s1_spectrum.txt`, `loopI_s2_spectrum.txt` (or `loopI_spectrum.txt`) | Fitted Loop I spectral coefficient: amplitude × reference spectrum. It multiplies that shell's spatial map; it is **not** the shell's sky-averaged intensity. |
| `loopI_correlations.json` | Correlation of each Loop I amplitude with every other parameter. |
| `loopI_spectral_fit.json`, `spectral_best_fit_parameters.npz` | Analytic Loop I only: fitted shape and Fermitools response checks. |

### Written by `plot` (same folder)

| File | Contents |
|---|---|
| `isotropic_spectrum.png/.pdf` | E² × isotropic intensity per bin. Without profiles, best-fit points only; with profiles, 68% errors (TS ≥ 4) or 95% upper limits. |
| `counts_spectrum.png` | Counts in the retained sky for data, total model and each component, with fractional residuals. |
| `residuals_<band>.png` | Maps of data − model, (data − model)/model and (data − model)/√model for each band. Gray = masked. The last one is a Pearson residual, not a significance. |
| `plot_bands.txt` | Requested versus actual (whole-bin) energy range of each band. |

## Using the isotropic spectrum elsewhere

`isotropic_spectrum.txt` is an **effective** background: it contains everything
isotropic at the fitted latitudes (extragalactic emission, unresolved sources
and residual cosmic rays in the event class), as seen with this event selection.

- Use it only with the **same event class, event type and IRFs**
  (here `P8R3_SOURCE_V3`, `evtype 3`).
- **Disable energy dispersion** for the isotropic component, as for the
  standard LAT isotropic templates.
- The values are intensities per steradian per MeV. Do not multiply by 4π or
  by the retained solid angle.
- The table runs from the first to the last bin **centre** (about 347 MeV to
  86 GeV for the default binning). A wider analysis range needs an explicit
  extrapolation; check how the receiving code handles it (`Main/run_3ML_master.py`
  extends the end-point slopes).
- The file gives best-fit values; it does not by itself establish a detection
  in every bin. Use `isotropic_bins.csv` (ideally with `compute_profiles: true`)
  for uncertainties.

To carry the whole all-sky model into the inner-Galaxy analysis, also copy the
fitted Sun, Moon and Loop I amplitudes from `diffuse_normalizations.json`. If
Loop I was fitted with an analytic shape, transfer that shape as well; the
amplitudes alone refer to that shape. The Main analysis names its Loop I
sources `Loop_I_S1`/`Loop_I_S2` and has its own `loopI_directory` setting.

## Conventions and caveats

- **Energy dispersion** is applied to catalog, GALPROP and additional diffuse
  components and disabled for the isotropic component.
- **Catalog sources** are fixed at their catalog parameters; their uncertainties
  are not propagated. The same holds for diffuse components with `free: false`.
- **Uncertainties.** Without profiles, `covariance.npz` holds the inverse
  Hessian. It is not meaningful for a parameter at the zero boundary (see the
  `boundary` column) and is set to NaN when the problem is degenerate. With
  profiles, intervals are profile-likelihood intervals in which all other
  amplitudes are re-fitted.
- **Mask edges** are applied at pixel centres. Model predictions in masked
  pixels (saved in `diagnostics.npz`) are extrapolations, not fit results.
- **Self-checks during `prepare`.** The templates must add up to Fermipy's full
  model, and changing each normalization by ±20% must change Fermipy's own
  log-likelihood by the same amount as the explicit formula used in `fit`.
  `prepare` stops if either check fails.

## Testing

```bash
python fit_isotropic.py self-test     # 16 synthetic tests; no LAT data or Fermitools
cd .. && python -m pytest -q tests    # repository tests (need Astropy and pytest)
```

These tests check the code (fit, profiles, mask, caching, reuse logic, Loop I
templates, plots). They do not validate the LAT response or a science result;
the `prepare` self-checks above do that on real data.

## Files in this folder

| File | Purpose |
|---|---|
| `fit_isotropic.py` | Command-line driver: configuration, `prepare`, `fit`, `plot`, `self-test`. |
| `config.yaml` | Example configuration (site-specific paths). |
| `joint_fit.py` | The joint Poisson likelihood, optimizer, covariance and profiles. |
| `loopI.py` | Turns the `loopI` block into sources; builds shell cubes; analytic shape fit. |
| `spectral_core.py`, `native_spectra.py` | Analytic Loop I spectra and their Fermitools evaluation. |
| `galprop_padding.py` | Energy padding of GALPROP cubes. |
| `reuse_products.py` | Reuse of products from an earlier preparation. |
| `plots.py` | Spectrum, count-spectrum and residual-map plots. |
| `slurm.sh` | Example Slurm job. |
