# Plots

Scripts for plotting the analysis results:

- **Inner-Galaxy fit** ([`Main/`](../Main/README.md)): count spectra, flux
  spectra and sky maps of the data, the model and the residuals, broken down by
  model component.
- **All-sky fit** ([`Isotropic/`](../Isotropic/README.md)): comparison of
  isotropic spectra from different runs with the official LAT isotropic
  template.

## Requirements

- Fermipy and Fermitools (for `plot_flux.py`).
- [`lat_binned_analysis`](https://github.com/ckarwin/Fermi-LAT_Binned_Analysis),
  used by `plot_flux.py`, `plot_spatial.py` and `custom_client.py`. Install it
  in the Fermitools environment:

  ```bash
  git clone https://github.com/ckarwin/Fermi-LAT_Binned_Analysis.git
  cd Fermi-LAT_Binned_Analysis
  pip install -e .
  ```
- NumPy, Matplotlib, Astropy and pandas.

## Inputs

The inner-Galaxy plots read the files written by that fit from `../output/`,
which should be the Fermipy working directory of that fit (`fermipy_outdir` in
`Main/run_3ML_master.py`). Edit the paths if your files are elsewhere.

| File | Written by the fit as |
|---|---|
| `ccube_00.fits` | Counts cube |
| `bexpmap_00.fits` | Exposure |
| `mcube_fit_roi_<savename>_00.fits` | Total model cube |
| `mcube_<component>_<savename>[_00].fits` | Model cube of each component (`pion_decay`, `ic`, `brem`, `isodiff`, `GCE_NFW`, `bubbles_low`, `bubbles_high`, `Sun`, `Moon`, Loop I) |
| `mcube_point_sources_<savename>.fits` | All sources except the diffuse components |

Set `savename` at the top of each script to the `savename` of the fit you want
to plot.

## Scripts

Run them from this folder.

### `compare_isotropic.py`: isotropic spectra

Plots E² I(E) of several `isotropic_spectrum.txt` files from all-sky runs
together with the official LAT isotropic template
(`iso_P8R3_SOURCE_V3_v1.txt`, from the Fermitools reference data) and its 1σ
band. Use it to compare runs with different foreground models, for example
with and without Loop I. Edit the file paths at the top to choose the runs. The
figure is shown on screen, not saved to a file. Needs only NumPy and
Matplotlib.

### `make_counts_plot.py`: count spectra

Sums the counts over the region in each energy bin and plots the data, the
total model and every component, with the fractional residuals
(data − model)/model below. Writes `counts_spectra.png`.

### `plot_flux.py`: flux spectra

The same breakdown in E² dN/dE (MeV cm⁻² s⁻¹). It works in two steps:

1. **Compute the fluxes.** The `calculate_flux` lines at the top divide each
   counts or model cube by the exposure and write `<label>_flux.dat` (for
   example `observed_flux.dat`, `pion_decay_flux.dat`). Uncomment the lines for
   the components you need.
2. **Plot.** The rest of the script reads the `.dat` files and writes
   `flux_spectra_<savename>.png`.

The script contains a `sys.exit()` after the flux calculations, so a run stops
after step 1. Remove it (and comment out the `calculate_flux` lines if the
`.dat` files already exist) to make the plot.

### `plot_spatial.py`: sky maps

Plots the counts map (`count_map.png`) and the total model map
(`model_map.png`), summed over energy, and the residual map data − model in
counts (`residuals.png`).

### `custom_client.py`: maps with `lat_binned_analysis`

Uses `lat_binned_analysis` to make counts, model and residual maps, including
fractional-residual maps in energy bands. Each step is a block of lines;
uncomment the ones you need. The settings for these tools are in
[`inputs.yaml`](inputs.yaml) (map size and projection, energy binning, data
files, IRFs and time range). The model XML written by Fermipy refers to
`$(FERMIPY_DATA_DIR)`; replace it with the actual path before using the XML
with these tools (see the comment at the top of the script).

## Notes

- Each script sets its own `savename`; `plot_spatial.py` currently uses a
  different one from the other two.
- The component file names follow the source names of the inner-Galaxy fit
  (for example `Loop_I_S1`). If you rename or add sources there, update the
  scripts to match.

## Files in this folder

| File | Purpose |
|---|---|
| `compare_isotropic.py` | Isotropic spectra of several all-sky runs vs. the official template |
| `make_counts_plot.py` | Count spectra by component, with residuals |
| `plot_flux.py` | Flux (E² dN/dE) spectra by component, with residuals |
| `plot_spatial.py` | Counts, model and residual sky maps |
| `custom_client.py` | Map making with `lat_binned_analysis` |
| `inputs.yaml` | Settings for `lat_binned_analysis` |
