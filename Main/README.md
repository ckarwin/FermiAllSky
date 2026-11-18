# Inner-Galaxy fit (3ML)

This folder fits the **40° × 40° region around the Galactic centre** with
Fermi-LAT data, using 3ML. The fit combines:

- a Galactic diffuse model from GALPROP (π⁰, inverse Compton, bremsstrahlung),
- the 4FGL-DR4 catalog sources,
- the foregrounds measured by the all-sky fit in [`Isotropic/`](../Isotropic/README.md):
  the isotropic background, the Sun, the Moon and Loop I,
- optionally, a Galactic-centre excess (GCE) with an NFW morphology and the
  Fermi bubbles.

Because it is built on 3ML, other instruments (INTEGRAL/SPI, COSI) can be added
to the same fit as additional plugins.

- [Requirements](#requirements)
- [Quick start](#quick-start)
- [The model](#the-model)
- [Settings](#settings)
- [Using the all-sky results](#using-the-all-sky-results)
- [Output files](#output-files)
- [Comparing the GCE with Ajello et al. (2016)](#comparing-the-gce-with-ajello-et-al-2016)
- [Changing the model](#changing-the-model)
- [Files in this folder](#files-in-this-folder)

## Requirements

- Fermipy and Fermitools with the LAT IRFs.
- The `diffuse-analysis` branches of
  [threeML](https://github.com/ckarwin/threeML/tree/diffuse-analysis) and
  [astromodels](https://github.com/ckarwin/astromodels/tree/diffuse-analysis).
  They provide the diffuse-emission features used here (GALPROP mapcubes as
  3ML sources, a custom isotropic spectrum, and the π⁰ spectral tilt).
- Astropy and Matplotlib (for the GCE comparison plot).
- Results from the other folders: the all-sky fit, the Sun/Moon and Loop I
  templates, and the NFW and bubble templates in `DM_Model/` and `Bubbles_Model/`.

## Quick start

```bash
cd Main
# 1. Edit the USER SETTINGS block at the top of run_3ML_master.py
#    (all-sky results, optional components, output name).
# 2. Check that config.yaml uses the same data selection as the all-sky fit.
# 3. Run (the fit takes many hours):
sbatch slurm.sh            # or: python run_3ML_master.py
```

Run the script from this folder. Results go to `output/` in the current
directory and to the Fermipy working directory (`fermipy_outdir`).

## The model

| Component | Name in the fit | Fitted parameters |
|---|---|---|
| 4FGL-DR4 sources within 30° of the GC | catalog names | none by default (see `free_point_source_radius`) |
| GALPROP π⁰ decay | `pion_decay` | normalization (+ spectral tilt with `fit_pi0_index`) |
| GALPROP inverse Compton | `ic` | normalization |
| GALPROP bremsstrahlung | `brem` | none (fixed at 1) |
| Galactic-centre excess (`include_dm`) | `GCE_NFW` | cutoff power law: `K`, `index`, `xc` |
| Fermi bubbles (`include_bubbles`) | `bubbles_low`, `bubbles_high` | one spectrum each (see `bubble_spectrum`) |
| Sun, Moon | `Sun`, `Moon` | none (all-sky values) |
| Loop I | `Loop_I_S1`, `Loop_I_S2` | none (all-sky values) |
| Isotropic | `isodiff` | none (all-sky spectrum, normalization 1) |

Spectral forms (energies in keV, pivot 1 GeV):

- Power law: dN/dE = K (E/piv)^index
- Log parabola: dN/dE = K (E/piv)^(alpha − beta ln(E/piv))
- Cutoff power law: dN/dE = K (E/piv)^index exp(−E/xc)
- Broken power law: dN/dE = K (E/piv)^alpha below the break energy `xb`, and
  K (xb/piv)^(alpha − beta) (E/piv)^beta above it (continuous at the break)

The GCE and bubble templates are normalized to 1, so their fitted spectrum is
the total flux of the component over the whole template (60° × 60° for the NFW
map, 80° × 80° for the bubbles), not only inside the 40° × 40° region.

The 4FGL and GALPROP part of the model is stored in `fermi_src_model.yaml` and
loaded by the script. To rebuild it, for example with a different GALPROP run,
set `make_src_model = True` and run once in the Fermipy environment.

## Settings

All settings are in the **USER SETTINGS** block at the top of
`run_3ML_master.py`.

| Setting | Meaning |
|---|---|
| `savename` | Label added to every output file name. |
| `make_src_model` | `True`: rebuild `fermi_src_model.yaml` from 4FGL-DR4 and the GALPROP files `pi0_file`, `ic_file`, `brem_file`. `False`: use the saved model. |
| `fit_pi0_index` | Also fit a power-law tilt of the π⁰ spectrum: I(E) = K (E / 1 GeV)^index × I_GALPROP(E), with `index` between −0.5 and 0.5 (positive = harder). |
| `include_dm`, `nfw_file` | Include the GCE and choose its NFW template. |
| `include_bubbles`, `bubble_spectrum` | Include the bubbles; spectrum `"powerlaw"`, `"log_parabola"`, `"cutoff_powerlaw"` or `"broken_powerlaw"` (break energy free between 1 and 100 GeV). |
| `fixed_templates`, `loopI_directory` | Files and fixed normalizations of the Sun, Moon and Loop I. |
| `isotropic_spectrum_file` | Isotropic spectrum from the all-sky fit. |
| `free_point_source_radius` | Free the normalization of every 4FGL point source within this distance of the GC (default 0°: none). |
| `fermipy_outdir` | Fermipy working directory. |

[`config.yaml`](config.yaml) is a standard Fermipy configuration: data files,
event selection, binning and instrument response. The source model comes from
the script, not from this file.

## Using the all-sky results

The fixed components must come from the all-sky fit you want to use:

1. **Isotropic spectrum.** Set `isotropic_spectrum_file` to
   `<outdir>/bestfit_<hash>/isotropic_spectrum.txt` from the all-sky fit. The
   script extends the table slightly beyond the analysis energy range when
   needed (`output/isotropic_spectrum_extrapolated.txt`).
2. **Sun, Moon and Loop I.** Use the same template files as the all-sky fit and
   copy their fitted normalizations from `diffuse_normalizations.json` into
   `fixed_templates`. Note the different names: `LoopI_S1`/`LoopI_S2` in the
   all-sky fit are `Loop_I_S1`/`Loop_I_S2` here. These entries are for the
   all-sky fit with separate shells and the empirical Loop I spectrum.
3. **Same data selection.** `config.yaml` must use the same event class, event
   type, IRFs, time range and zenith cut as the all-sky fit. The isotropic
   spectrum is only valid for that selection.
4. **Same GALPROP model** as in the all-sky fit (`fermi_src_model.yaml` stores
   the GALPROP file paths).

## Output files

In `output/`:

| File | Contents |
|---|---|
| `threeML_fit_parameters_<savename>.csv` | Best-fit values and errors of the free parameters |
| `threeML_likelihood_values_<savename>.csv` | Likelihood values |
| `threeML_best_fit_model_<savename>.yml` | Best-fit model (reload with `load_model`) |
| `mcube_point_sources_<savename>.fits` | Model counts of all sources except the diffuse components |
| `isotropic_spectrum_extrapolated.txt` | Isotropic spectrum used in the fit |

In the Fermipy working directory (`fermipy_outdir`):

| File | Contents |
|---|---|
| `fit_roi_<savename>.*` | The fitted region and total model map |
| residual maps and plots | Data minus model |
| `<component>_<savename>*` | Model-count map of each diffuse component |

In the current directory: `GCE_Ajello2016_comparison.pdf` (when `include_dm = True`).

To plot spectra and maps of these results by component, see [`Plots/`](../Plots/README.md).

## Comparing the GCE with Ajello et al. (2016)

With `include_dm = True`, the script finishes by plotting the fitted GCE
spectrum (best fit and 68% band) against the four GCE results of Ajello et al.
(2016), Figure 13, in the same region of the sky (by default the paper's
15° × 15° region). This is done by `plot_gce_comparison` in
[`compare_gce_ajello2016.py`](compare_gce_ajello2016.py), which can also be
used on its own:

```python
from compare_gce_ajello2016 import plot_gce_comparison

plot_gce_comparison(results, "../DM_Model/nfw_gamma1.fits",
                    comparison_roi=(-20, 20, -20, 20),   # 40 x 40 deg region
                    output="GCE_comparison_40x40.pdf")
```

Here `results` is a 3ML results object: `jl.results` after a fit, or one read
with `load_analysis_results` if you saved it with
`jl.results.write_to("output/results.fits")`. The paper curves are read from
the published figure and are accurate to a few percent; see the module
docstring for details.

## Changing the model

Keep these rules in mind when you add or change components:

- **Mapcube sources** (GALPROP, Sun, Moon, Loop I, or any other 3D intensity
  cube loaded with `GalpropMap`) are scaled by the spatial normalization `K`.
  Give them a fixed `Constant(k=1)` spectrum, or a `Powerlaw` tilt with `K`
  fixed at 1 and a fixed pivot (as for π⁰).
- **2D templates** (`SpatialTemplate_2D`, as for the GCE and bubbles) should
  integrate to 1, with `K` fixed at 1; the spectrum then carries the flux.
- **Energy grid of mapcubes.** The LAT likelihood reads mapcube energies from
  the file's `ENERGIES` table. When astromodels itself evaluates a mapcube (for
  example in a joint fit with SPI), it reads the energies from the cube's third
  axis as log10(E/MeV). GALPROP cubes follow that convention; the Sun, Moon and
  Loop I cubes do not and need converting before such use.
- **Energy dispersion** is applied to all components except the isotropic one.
  The number of extra energy bins is set by `gtlike.edisp_bins` in
  `config.yaml` (`-1` here; the all-sky fit uses `-2`).

## Files in this folder

| File | Purpose |
|---|---|
| `run_3ML_master.py` | Builds the model, runs the fit and writes the results |
| `config.yaml` | Fermipy configuration (data, selection, binning, IRFs) |
| `fermi_src_model.yaml` | Saved 4FGL + GALPROP model |
| `compare_gce_ajello2016.py` | GCE comparison with Ajello et al. (2016) |
| `GCE_Ajello2016_comparison.pdf` | Example comparison plot |
| `slurm.sh` | Example Slurm job |
