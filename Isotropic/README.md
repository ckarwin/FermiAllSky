# All-sky isotropic fit pipeline

The maintained entry point is [`fit_isotropic.py`](fit_isotropic.py). It prepares
full-sky LAT response products, then fits a **masked** Poisson likelihood with
one isotropic amplitude per energy bin and one amplitude across all energies
for each free diffuse component. Catalog point and extended sources stay fixed.
The current preparation/cache version is **6**.

## Quick start

Run these commands from the repository root, in your Fermipy/Fermitools environment:

```bash
python Isotropic/fit_isotropic.py self-test
# Edit Isotropic/config.yaml for your data and model paths first.
python LoopI/build_templates.py --config Isotropic/config.yaml
python Isotropic/fit_isotropic.py run Isotropic/config.yaml
```

For a **first preparation**, remove `reuse_from` or set it to `null`, and choose an
unused `outdir`. The supplied config is an example of a follow-up run: it points
to a previous output directory that is not included in the repository. Its data,
GALPROP, catalog, extended-source, and Sun/Moon paths are site-specific examples.
Change those paths to files available in your environment.

The Loop I build is needed once on a fresh checkout: the reconstruction inputs
are tracked, but the large generated FITS files are not. Skip it if Loop I is
disabled or the required cubes already exist. Rebuilding overwrites the generated
assets and invalidates affected preparation caches; do not rebuild during a fit.

`prepare` requires Fermipy, Fermitools/pyLikelihood, the LAT IRFs, event and
spacecraft data, a catalog and its extended templates. NumPy, SciPy, PyYAML,
Astropy, and Matplotlib support the local utilities. `self-test` and empirical
cached `fit`/`plot` do not require a live LAT response calculation. **Every analytic
Loop I fit uses Fermipy/Fermitools**, including an analytic shape with fixed
parameters. The regression suite additionally uses pytest:

```bash
python -m pytest -q tests
```

## What each command does

The following commands assume `cd Isotropic`; `config.yaml` is the wrapper
configuration, not raw Fermipy YAML.

| Command | Action | When to use it |
|---|---|---|
| `python fit_isotropic.py prepare config.yaml` | Validates inputs, generates separate Loop I cubes if requested, prepares LAT event/exposure/source-map products, checks the native likelihood, and saves `templates.npz`. | First run, or changed inputs/response in a **new** `outdir`. |
| `python fit_isotropic.py fit config.yaml` | Loads and checks the prepared cache, applies the mask, fits amplitudes and optional Loop I shape, and saves results. | Fit after preparation, or change fit-only controls. Does not automatically prepare stale inputs. |
| `python fit_isotropic.py plot config.yaml` | Reads the matching saved fit and writes spectrum/count/residual plots. | Regenerate plots or change plot settings without refitting. |
| `python fit_isotropic.py self-test` | Runs 16 synthetic numerical, cache/export, and plotting tests. No YAML is needed. | Check the local installation and numerical implementation; this is not a real LAT data test. |
| `python fit_isotropic.py run config.yaml` | Runs `prepare`, then `fit`, then `plot`. | A complete run into a fresh output directory. It is **not** an automatic resume command. |

Preparation refuses to reuse an existing `outdir/fermipy/` as its destination,
including one left by a failed preparation. A fit refuses to overwrite its result
directory. Use `plot` for existing results, move an incomplete result directory
aside before refitting, or select a new output for a new preparation.

## Output paths and reuse

```yaml
outdir: output_with_solar_and_loopI_separate
reuse_from: output_with_solar_and_loopI
```

`outdir` is the **new run's destination**. `reuse_from` is an **old run's root**,
containing `fermipy_config.yaml` and `fermipy/`, from which compatible products may
be copied. It is not the old `bestfit_<hash>/` directory or `fermipy/` itself.
They must differ. Old products are read-only; copies are used in the new run.

Relative paths resolve against the directory containing the YAML, independent of
the shell's working directory. Thus the example above writes under `Isotropic/`
when used with `Isotropic/config.yaml`. Absolute paths, `~`, and environment
variables are supported for configured input/output paths.

There must be **one `outdir` key**. The old example contained two; normal YAML
loading silently used the later `output_with_solar_and_loopI_separate` value.
The duplicate has been removed, preserving that effective destination. The
wrapper now rejects duplicate keys, including nested duplicates, with a line
number instead of silently changing your configuration.

With `reuse_from`, data, selection, binning, and response settings must match.
The code copies compatible selected events, count cube, livetime cube, and
exposure map. It checks source definitions and referenced-file content before
reusing individual source maps, including maps stored in numbered FITS chunks.
New, changed, missing, or invalid source maps are left for native generation.
`reuse_report.json` records those decisions. A changed event selection or IRF
requires a fresh preparation **without** `reuse_from`; reuse fails explicitly
rather than silently falling back. Data files should be immutable.

Version 6 also checks the underlying Loop I base cube, spatial maps, and spectral
metadata/table for changes. Existing version-5 `templates.npz` files need a new
preparation in a new `outdir`; `reuse_from` can still reuse compatible expensive
LAT products from that old run.

## Loop I directory and model options

The default asset directory is now the repository-level [`LoopI/`](../LoopI/),
next to `Isotropic/`. It is located relative to the driver, not the working
directory. **No filename is needed** in the normal integrated configuration:

```yaml
loopI:
  enabled: true
  mode: separate
  spectrum: empirical
  free_norm: true
  normalization: 0.14
  # directory: ../LoopI  # Optional override, relative to THIS YAML.
```

If you move the assets elsewhere, use `directory: /absolute/path/to/LoopI` (or an
appropriate relative path). The directory must contain the generated cubes,
spatial maps, `model_metadata.json`, and `empirical_intensity.txt`. Metadata,
normalization checks, and spectrum export all use this same directory.

The builder reads `loopI.directory` when present. An explicit `--outdir` overrides
that choice for the build only; make sure it agrees with the fit configuration:

```bash
python LoopI/build_templates.py --config Isotropic/config.yaml --outdir /path/to/LoopI
```

By default the builder reads `Isotropic/config.yaml` and writes beside itself.
It produces the combined empirical and power-law seed mapcubes plus the combined,
S1, and S2 spatial maps, metadata, spectrum table, and diagnostic plots. `prepare`
then creates the requested S1/S2 mapcubes from these inputs, preserving the common
spatial normalization. At unit amplitudes, S1+S2 equals the combined model.

| Option | Meaning |
|---|---|
| `enabled: false` | Omit the integrated Loop I sources. |
| `mode: combined` | One source, `LoopI`; one amplitude multiplies the fixed S1+S2 sum. This is also the default when `mode` is omitted. |
| `mode: separate` | Two sources, `LoopI_S1` and `LoopI_S2`; independent amplitudes, with a shared spectral shape. |
| `free_norm: true` | Fit the amplitude: **one per shell** in separate mode, one total in combined mode. |
| `free_norm: false` | Hold the amplitude at `normalization`. |
| `normalization: 0.14` | Starting amplitude when free; fixed multiplier otherwise. It does not rescale the input files. |
| `spectrum: empirical` | Use the reconstructed curve with fixed spectral shape; no `parameters` block. |

Optional shell-specific controls, nested inside `loopI`, override the shared
normalization defaults in separate mode:

```yaml
  shells:
    S1: {free_norm: true, normalization: 0.14}
    S2: {free_norm: false, normalization: 0.20}
```

The current driver also supports analytic shape fitting directly; no separate
spectral-fit script is required:

```yaml
loopI:
  enabled: true
  mode: separate
  spectrum: cutoff_powerlaw
  free_norm: true
  normalization: 0.14
  parameters:
    index: {value: 2.4, free: true, min: 0.5, max: 5.0}
    cutoff_MeV: {value: 30000, free: true, min: 1000, max: 10000000}
compute_profiles: false
```

`powerlaw` accepts `index`; `cutoff_powerlaw` accepts `index` and `cutoff_MeV`;
`logparabola` accepts `alpha` and `beta`. Analytic indices are positive for falling
spectra, e.g. `index: 2.4` means E^-2.4. Shape parameters are tied between S1/S2;
amplitudes can differ. The spectrum is normalized at 1000 MeV. Unspecified
analytic parameters default to free; set `free: false` explicitly to freeze them.

Analytic shapes are refolded through the native response during `fit`, before the
final amplitude fit. Free shape fits require `compute_profiles: false`; their
saved linear covariance is **conditional on the fitted shape**, not a full joint
shape/amplitude covariance. Fixed-shape normalization profiling is supported.

The empirical curve is a reconstruction, not an author-supplied flux table.
Its amplitude convention and spectral extrapolation limits are documented in
[`LoopI/MODEL_NOTES.md`](../LoopI/MODEL_NOTES.md).

Moving assets changes absolute input paths and may leave old XML/config files
pointing at the old location. Set the new directory and prepare a new output;
do not assume an existing analytic fit can still reopen the old Fermipy work
area. Reuse may safely regenerate any sources whose old references cannot be
verified. `Main/run_3ML_master.py` has its own `loopI_directory` setting; update it
too if you move assets beyond the default repository location.

## What is in the `src` file? What happens when a model changes?

There are several different source-related files:

| File | Contents |
|---|---|
| `fermipy/srcmdl_00.xml` | Source definitions for the LAT model: catalog point/extended sources and configured diffuse components. |
| `fermipy/srcmap_00.fits` and numbered chunks | Response-dependent source maps for catalog **and diffuse** sources; not a point-source-only file. |
| `templates.npz` | Full-sky counts, per-component predicted-count templates, fixed catalog contribution, energy grid, pixel coordinates, seed, and preparation metadata. |
| [`Main/fermi_src_model.yaml`](../Main/fermi_src_model.yaml) | Saved 3ML model. The checked-in model includes catalog sources **and** GALPROP pion-decay, IC, and bremsstrahlung components with file paths. |

The isotropic wrapper builds its model from **its own YAML**. It does not read
`Main/fermi_src_model.yaml` or automatically discover model files in repository
folders. `Main/run_3ML_master.py` loads that saved model by default
(`make_src_model = False`) and adds its other components separately. Changing an
isotropic YAML setting does not update this saved 3ML model. Rebuild or edit the
3ML model and set the matching fitted amplitudes/input paths when switching
GALPROP models in the main analysis.

To use a new GALPROP cube, update the corresponding `galprop.<name>.file`, choose
a new `outdir`, and run `prepare`/`run`. With compatible response settings,
`reuse_from` preserves unchanged catalog/component maps and recalculates changed
ones. `fit` alone never recalculates a changed GALPROP cube. Input cube size/mtime
changes invalidate the fit cache; selective source reuse additionally uses
content fingerprints. Do not edit science inputs while jobs are running or
preserve their timestamps to bypass these checks.

To include another diffuse component in the isotropic analysis, add it explicitly:

```yaml
additional_diffuse:
  Bubbles:
    file: /path/to/bubbles_intensity_mapcube.fits
    free: true
    normalization: 1.0
```

This must be a **3D differential-intensity mapcube** with celestial WCS and an
`ENERGIES` table in MeV, covering the analysis plus energy-dispersion margins.
The 2D morphology maps in [`Bubbles_Model/`](../Bubbles_Model/) need a physical
spectrum and conversion to such a cube before use here. The isotropic config
currently excludes a longitude band using the mask and does not include a
Bubbles source by default. A new bubble file is not picked up automatically just
because it exists. Add/update it and prepare a new output as above. Sun and Moon
likewise enter through explicit `additional_diffuse` paths; regenerate their
time-dependent templates when the observation/response settings change. Do not
also add `LoopI`, `LoopI_S1`, or `LoopI_S2` there when using integrated `loopI`.

| Change | Required action |
|---|---|
| Mask or diffuse fixed/free/amplitude controls | `fit`, then `plot`, using the same prepared output. |
| Plot bands or display resolution | `plot` only. |
| Analytic Loop I parameters or analytic spectral family | `fit`, then `plot` using the same power-law seed preparation. |
| Empirical ↔ analytic Loop I; combined ↔ separate; add/remove a component | New `outdir`, `prepare` then `fit`/`plot`; compatible LAT products may be reused. |
| Input cube/path, Loop I assets, isotropic reference, or extrapolation settings | New `outdir` and preparation. |
| Data, time/energy selection, pixelization, event class/type, IRF, or edisp settings | New `outdir`, fresh preparation without `reuse_from`. |

## Fit conventions, mask, and plotting

The example config uses 300–100000 MeV, 8 bins/decade, `hpx_order: 8`
(NSIDE=256), and `edisp_bins: auto` (resolved to -2 for these settings).
`binsz` is bookkeeping in this HEALPix analysis; `hpx_order` sets the sky grid.
For a cheaper native smoke test, use order 6 and a separate output, then check
angular convergence before interpreting the science result.

The mask retains |b| >= 30 degrees and |dec| <= 90 degrees and removes the region
|wrapped l| <= 45 degrees, |b| <= 90 degrees in the shipped config. Preparation
remains full-sky so response folding precedes masking. The likelihood and count
spectrum use **retained pixels**, not the full sky. The GC-centred 40°×40° region
is excluded by this mask. The Loop I overlap figure shows that rectangle only as
a diagnostic; it does not select the fit region. Predictions outside the mask
are extrapolations. `plots.bands_MeV` changes map energy bands, not the fit ROI.

GALPROP and additional diffuse normalizations are global across energies:
`free: false` freezes the supplied multiplier; `free: true` fits it. Uncertainty
in fixed components is not propagated. All catalog parameters remain fixed.
The isotropic reference is generated internally (default example index -2.0),
with a free nonnegative amplitude in each energy bin. Rescaling its arbitrary
reference intensity inversely rescales those amplitudes, preserving fitted
physical intensity. The reference index determines the within-bin shape.

The result is an **effective isotropic background**, including residual charged
particles. Energy dispersion is disabled for `isodiff` and enabled for other
components. Keep it disabled when reusing `isotropic_spectrum.txt`. Intensities
are per sr per MeV: do not multiply by 4π or the retained solid angle to use the
file as an isotropic intensity spectrum.

GALPROP padding can extend boundary spectra within the explicitly configured
`galprop_extrapolation` limits; originals are preserved and decisions are recorded
in `galprop_extrapolation.json`. Additional diffuse cubes are validated for energy
coverage but are not padded automatically. `edisp_bins: auto` uses positive
spatial-map mixing below 100 MeV (+4 at 20 MeV with 8 bins/decade). Low-energy
response/data validity still needs checking in the installed LAT environment.

## Results and diagnostics

Preparation writes `templates.npz`, `fermipy_config.yaml`,
`preparation_checks.json`, `source_fingerprints.json`, `isotropic_seed_used.txt`,
and native products under `fermipy/`. Reuse also writes `reuse_report.json`.

Fits write `bestfit_<hash>/` by default or `fit_<hash>/` with
`compute_profiles: true`. The hash identifies fit settings; changing only plot
options does not create a new fit. Useful outputs are:

| Output | Meaning |
|---|---|
| `isotropic_bins.csv` | Bin intensities, counts, deviance, and profile results when requested. |
| `diffuse_normalizations.json` | All fitted/fixed diffuse multipliers. GALPROP and additional components also have separate JSON summaries. |
| `best_fit_parameters.npz`, `covariance.npz` | Amplitude checkpoint, covariance, and parameter order. |
| `fit_report.json` | Configuration, fit convention, mask area, and validation/uncertainty information. |
| `diagnostics.npz` | Full-sky counts and fitted predictions with the retained mask. |
| `loopI_spectral_fit.json` | Analytic shape results and native response checks, when applicable. |
| `loopI_spectrum.txt` or `loopI_s1_spectrum.txt`/`loopI_s2_spectrum.txt`; `loopI_correlations.json` | Spectral coefficients and amplitude correlations. Each coefficient multiplies its corresponding spatial map; it is not that shell's sky-mean intensity. |
| `isotropic_spectrum.png/.pdf`, `counts_spectrum.png` | Intensity spectrum and retained-sky component count spectra. |
| `residuals_<band>.png`, `plot_bands.txt` | Unsmoothed residual maps and actual whole-bin energy coverage. Excluded pixels are gray; Pearson residuals are not TS significances. |

With `compute_profiles: false`, profile intervals, upper limits, and TS in the
CSV are NaN; the saved Hessian covariance is separate. With profiles enabled,
all other free amplitudes are reoptimized for each interval. Free Loop I shape
parameters currently preclude this profiling mode.

`isotropic_spectrum.txt` is exported when every best-fit intensity is positive;
with profiles enabled it additionally requires TS >= 4 in every bin and valid
covariance. Its absence is therefore not necessarily an I/O failure. Check the
CSV and report. The text spectrum interpolates bin-centre values; its stated
valid range is the first to last bin **centre**, not the full bin-edge range.
`Main/run_3ML_master.py` currently extrapolates endpoint slopes when its requested
range is wider, so inspect that choice when transferring a fit. Set its spectrum
path and fixed Sun/Moon/Loop I amplitudes to the outputs of the matching run.
The main script currently loads **empirical** Loop I cubes. After an analytic
Loop I fit, transfer the fitted spectral shape as well; applying those amplitudes
to the empirical cubes alone would use a different model.

## Validation and remaining environment requirements

See [`VALIDATION.txt`](VALIDATION.txt) for the commands actually run. Synthetic
checks cover the joint Poisson optimizer, profiles, mask, covariance/boundaries,
reference normalization, cached export/plots, full-resolution Loop I cube build,
shell sums and analytic seeds, path overrides, duplicate YAML rejection, cache
invalidation, selective source-map reuse, sparse-map handling, and solar GTIs.

Real LAT preparation and analytic shape fits additionally check source inclusion,
full-sky geometry, model sum closure, and native-versus-explicit likelihood
changes. They require the installed Fermipy/Fermitools and your data; the local
checks do not establish instrument calibration or validate a scientific fit.

The repository's `Main/` driver and Slurm scripts remain site-specific analysis
scripts: review external paths, output locations, environment activation, and
fitted amplitudes before execution. The main driver uses custom `GalpropMap` and
`FermipyLike` behavior supplied by your installed 3ML/astromodels environment.
`Sun_and_Moon/fit_isotropic.py` is a compatibility entry point to this maintained
driver; use `Isotropic/fit_isotropic.py` for new commands. The solar template
generator and its Slurm instructions live in [`Sun_and_Moon/`](../Sun_and_Moon/).
