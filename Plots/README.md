# Plots

Scripts for plotting the analysis results:

- **One inner-Galaxy fit** ([`Main/`](../Main/README.md)):
  [`make_all_plots.py`](#make_all_plotspy-all-plots-of-one-fit) makes every
  plot of a fit (count and flux spectra by component, sky maps, residual maps,
  and comparisons of the GCE, Fermi bubbles and Loop I with published results).
  `run_3ML_master.py` calls it at the end of each fit.
- **Several inner-Galaxy fits**:
  [`compare_fits.py`](#compare_fitspy-several-fits-overlaid) overlays their
  spectra and residuals, for example fits with different GALPROP models.
- **All-sky fits** ([`Isotropic/`](../Isotropic/README.md)):
  [`compare_isotropic.py`](#compare_isotropicpy-isotropic-spectra) compares
  their isotropic spectra with the official LAT template.

## Requirements

- NumPy, Matplotlib, Astropy, pandas, PyYAML and healpy.
- [`lat_binned_analysis`](https://github.com/ckarwin/Fermi-LAT_Binned_Analysis),
  for the sky maps and residual maps. Install it in the Fermitools environment:

  ```bash
  git clone https://github.com/ckarwin/Fermi-LAT_Binned_Analysis.git
  cd Fermi-LAT_Binned_Analysis
  pip install -e .
  ```

  The residual maps need a version whose `residual_map2` has the `mode`,
  `bins` and `outdir` options. Its settings are in [`inputs.yaml`](inputs.yaml)
  (map size, energy binning, data files, IRFs, time range).
- Fermipy and Fermitools, for the standalone `plot_flux.py` and
  `custom_client.py`.

## `make_all_plots.py`: all plots of one fit

### After a fit

`Main/run_3ML_master.py` calls it at the end of the fit when `make_plots = True`
(the default). The output folder is `Plots/<plots_subdir>`, set at the top of
`run_3ML_master.py` (default: the run's `savename`). The fit's own outputs are
found automatically.

Each plot is made on its own: if one fails (for example because a package is
missing), the error is written to `make_all_plots.log` and the other plots are
still made. The fit results are never affected.

### Again later

```bash
python make_all_plots.py /path/to/Plots/<plots_subdir>
```

This reads `run_info.json` in that folder (written by the first call) and makes
all plots again from the fit's files. The GCE comparison needs threeML (it
reloads the saved 3ML results).

### Output files

| File | Contents |
|---|---|
| `counts_spectra.png` | Counts per energy bin: data, total model and each component, with (D − M)/M below |
| `counts_spectrum.csv` | The numbers of that plot |
| `<name>_flux.dat` | Flux of each component in the region, per energy bin (`observed`, `model`, `pion_decay`, `ics`, `bremss`, `isodiff`, `pcs`, `gce`, `bubbles_low`, `bubbles_high`, `sun`, `moon`, `loopI_SI`, `loopI_S2`) |
| `flux_spectra_<savename>.png` | E² dN/dE of each component, with (D − M)/M below |
| `count_map.png`, `model_map.png` | Counts and model summed over energy |
| `residuals.png`, `residuals_smoothed.png` | Data − model summed over energy, unsmoothed and smoothed |
| `GC_<i>-<j>_residual.fits` | Data − model in energy bins i to j−1 |
| `GC_<i>-<j>_fractional_residual.fits` | (Data − model) / model |
| `GC_<i>-<j>_sqrt_residual.fits` | (Data − model) / √model |
| `GC_<i>-<j>_*_smoothed.fits` | The same three maps after smoothing |
| `residuals_<diff/frac/sqrt>_<i>-<j>.png`, `..._smoothed.png` | Each map, unsmoothed and smoothed |
| `residual_maps_<i>-<j>.png`, `..._smoothed.png` | The three maps of one energy range side by side, unsmoothed and smoothed |
| `GCE_Ajello2016_comparison.png`, `gce_spectrum.csv` | GCE spectrum vs. Ajello et al. (2016) ([Main](../Main/README.md#comparing-the-gce-with-ajello-et-al-2016)) |
| `correlation_matrix.png` | Correlation matrix of the free parameters after the fit |
| `correlation_matrix.csv`, `covariance_matrix.csv` | The numbers behind it, with full parameter names |
| `threeML_analysis_results.fits` | The 3ML results of the fit (used to redo the correlation matrix and GCE plots) |
| `bubbles_comparison.png`, `bubbles_*.csv` | Bubble spectra vs. Ackermann et al. (2017) ([`compare_bubbles/`](compare_bubbles/README.md)) |
| `loopI_comparison.png`, `loopI_*.csv` | Loop I of the all-sky fit vs. published spectra ([`compare_loopI/`](compare_loopI/README.md)) |
| `run_info.json` | Inputs of the fit, used to redo the plots and by `compare_fits.py` |
| `make_all_plots.log` | Which plots were made, and the errors of any that failed |

The GCE and bubble files are made only when the fit includes those components.

The correlation matrix is the fit's covariance matrix divided by the product
of the two parameters' errors, so every entry is between −1 and +1. The
parameters' values span many orders of magnitude, which would make a plot of
the covariance itself unreadable. Values near ±1 mean two parameters are
nearly degenerate (for example the GCE index and cutoff).

Its axis labels are short: the source name and the parameter name, e.g.
`GCE_NFW xc` for `GCE_NFW.spectrum.main.Cutoff_powerlaw.xc`. If two labels
would be the same, just enough of the path is added to tell them apart (e.g.
`pion_decay GalpropMap.K` and `pion_decay Powerlaw.K`). The CSV files keep
the full names.

### Settings

At the top of `make_all_plots.py`:

| Setting | Default | Meaning |
|---|---|---|
| `COUNTS_YLIM`, `COUNTS_RESID_YLIM` | `(1e2, 1e7)`, `(-1, 1)` | Y ranges of the counts spectrum and its (D − M)/M panel |
| `FLUX_YLIM`, `FLUX_RESID_YLIM` | `(1e-6, 1)`, `(-0.1, 0.1)` | Y ranges of the flux spectrum and its (D − M)/M panel |
| `RESIDUAL_SMOOTH_SIGMA_DEG` | `1.0` | Gaussian smoothing of the `*_smoothed` residual maps (σ in degrees; 0 for unsmoothed maps only) |
| `RESIDUAL_BINS` | `[[0, 5], [5, 12], [12, 20]]` | Energy-bin ranges of the residual maps |
| `PARAMETER_LABELS` | `{}` | Correlation-matrix axis labels to replace, `{full parameter name: label}` |
| `MAX_LABEL_LENGTH` | `28` | Longer correlation-matrix labels have their source name shortened with "…" |

To change a setting for a finished fit, edit it and run
`python make_all_plots.py Plots/<plots_subdir>`.

Every PNG is saved at 300 dpi: `DPI` in `plot_style.py`.

### Smoothing of the residual maps

Single 0.2° pixels hold few counts, so unsmoothed residual maps are dominated
by Poisson noise. The data and the model are smoothed with the same Gaussian
kernel before each residual is formed:

- data − model: the mean residual per pixel within the kernel, in counts per
  pixel (as unsmoothed);
- (data − model) / model: smoothed data over smoothed model, minus 1;
- (data − model) / √model: the smoothed excess divided by its Poisson
  uncertainty, so it reads as a significance (about ±1 where the model is
  right).

σ = 1° is about the LAT PSF (68% containment) at 1 GeV. The PSF is broader at
lower energies and narrower at higher energies, so features smaller than the
kernel are smoothed away in the high-energy maps, so each residual map is
also plotted unsmoothed.

The default residual-map energy ranges are 300 MeV–1.28 GeV, 1.28–9.79 GeV
and 9.79–100 GeV (energy bins 0–4, 5–11 and 12–19 for 300 MeV–100 GeV at 8
bins per decade). The side-by-side residual figure draws each map exactly as
its single-map plot.

The flux tables are computed exactly as `lat_binned_analysis`'s
`calculate_flux` computes them (same sum over pixels, exposure lookup and file
format), but much faster.

## `compare_fits.py`: several fits overlaid

Overlays the results of several fits made with `make_all_plots.py`. List them
at the top of the script, as plot folder and legend label:

```python
RUNS = [
    ("galprop_model_A", "GALPROP model A"),
    ("galprop_model_B", "GALPROP model B"),
]
OUTPUT_DIR = "compare_fits_output"
```

```bash
python compare_fits.py
```

It writes to `OUTPUT_DIR`:

| File | Contents |
|---|---|
| `gce_overlay.png` | GCE spectrum of each fit (best fit and 68% band) over the Ajello et al. (2016) results |
| `bubbles_overlay.png` | Bubble spectra of each fit over Ackermann et al. (2017) |
| `loopI_overlay.png` | Loop I spectrum of the all-sky fit used by each fit, over the published spectra |
| `residuals_overlay.png` | (D − M)/M of each fit |

A fit without a component (for example no GCE) is left out of that plot. The
inner-Galaxy fit keeps Loop I fixed at the all-sky result, so fits that use the
same all-sky fit have the same Loop I curve.

## `compare_isotropic.py`: isotropic spectra

Compares the isotropic spectra of several all-sky fits with the official LAT
isotropic template (`iso_P8R3_SOURCE_V3_v1.txt`) and its 1σ band. Below the
spectra, it shows each fit's (D − M)/M, summed over the sky that the fit
used, as in that fit's own counts spectrum. List the fits at the top of the
script, as results folder (`Isotropic/<outdir>/bestfit_<hash>`), label and
colour:

```python
FITS = [
    (".../output_with_solar/bestfit_34e87a564c", "Sun and Moon", "cyan"),
    ...
]
```

```bash
python compare_isotropic.py
```

The figure is shown on screen and saved as `OUTPUT` (default
`compare_isotropic.png`).

## Comparisons with published results

- [`compare_bubbles/`](compare_bubbles/README.md): fitted low- and
  high-latitude bubble spectra vs. Ackermann et al. (2017), Figure 9a,
  converted to the same 40° × 40° region.
- [`compare_loopI/`](compare_loopI/README.md): Loop I spectrum of the all-sky
  fit vs. published Fermi-LAT Loop I spectra, each averaged over its own sky
  region.

Both are run by `make_all_plots.py` and `compare_fits.py`, and can also be run
on their own.

## Standalone scripts

These are the scripts that `make_all_plots.py` brings together. They are still
useful for checking a single plot. Run them from this folder, after setting
`savename` and the folder of the fit's files at the top:

```python
OUTPUT_DIR = "../Main/output"
```

`OUTPUT_DIR` is the Fermipy working directory of the fit (`fermipy_outdir`),
which is assumed to also hold `mcube_point_sources_<savename>.fits` (written
to `output/` in the run directory).

| Script | Makes |
|---|---|
| `make_counts_plot.py` | `counts_spectra.png` |
| `plot_flux.py` | `<name>_flux.dat` (with `calculate_flux`), then `flux_spectra_<savename>.png`. It contains a `sys.exit()` after the flux calculations; remove it to make the plot. |
| `plot_spatial.py` | `count_map.png`, `model_map.png`, `residuals.png` (unsmoothed) |
| `custom_client.py` | Maps with `lat_binned_analysis`, including residual maps in energy ranges. The model XML written by Fermipy refers to `$(FERMIPY_DATA_DIR)`; replace it with the actual path before using it here. |

The component file names follow the source names of the inner-Galaxy fit (for
example `Loop_I_S1`). If you rename or add sources there, update
`COMPONENTS` in `make_all_plots.py` (and the standalone scripts) to match.

## Files in this folder

| File | Purpose |
|---|---|
| `make_all_plots.py` | All plots of one fit (called by `run_3ML_master.py`) |
| `compare_fits.py` | Overlays of several fits |
| `compare_isotropic.py` | Isotropic spectra of several all-sky fits vs. the official template |
| `compare_bubbles/` | Fermi-bubble spectra vs. Ackermann et al. (2017) |
| `compare_loopI/` | Loop I spectrum of the all-sky fit vs. published results |
| `plot_style.py` | Shared font sizes and colours of the comparison plots |
| `make_counts_plot.py`, `plot_flux.py`, `plot_spatial.py`, `custom_client.py` | Standalone versions of single plots |
| `inputs.yaml` | Settings for `lat_binned_analysis` |
