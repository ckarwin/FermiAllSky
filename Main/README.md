# Inner-Galaxy fit (3ML)

This folder fits the **40° × 40° region around the Galactic centre** with
Fermi-LAT data, using 3ML and its `FermipyLike` plugin. The Galactic diffuse
emission is described by GALPROP components with free normalizations. The
foregrounds measured by the all-sky fit ([`Isotropic/`](../Isotropic/README.md))
are held fixed: the isotropic background, the Sun, the Moon and Loop I. On top
of that, the fit can include a Galactic-centre excess (GCE) with an NFW
morphology and the Fermi bubbles.

The 3ML framework is what will later combine this LAT likelihood with
INTEGRAL/SPI (and COSI) in a joint fit: other instruments are added as further
plugins in the same `DataList`.

- [Requirements](#requirements)
- [Quick start](#quick-start)
- [The model](#the-model)
- [Settings](#settings)
- [Carrying over the all-sky results](#carrying-over-the-all-sky-results)
- [Output files](#output-files)
- [GCE comparison with Ajello et al. (2016)](#gce-comparison-with-ajello-et-al-2016)
- [Notes on the custom threeML/astromodels branches](#notes-on-the-custom-threemlastromodels-branches)
- [Files in this folder](#files-in-this-folder)

## Requirements

- Fermipy and Fermitools with the LAT IRFs.
- The `diffuse-analysis` branches of
  [threeML](https://github.com/ckarwin/threeML/tree/diffuse-analysis) and
  [astromodels](https://github.com/ckarwin/astromodels/tree/diffuse-analysis).
  They provide the `FermipyLike` options used here (`galactic_diffuse`,
  `isotropic_spectrum`, `skip_source_diagnostics`) and `GalpropMap`, which
  turns a LAT mapcube into an astromodels source. The script relies on
  features specific to these branches.
- Astropy and Matplotlib for the GCE comparison plot.
- Inputs from the other folders (see [Carrying over the all-sky results](#carrying-over-the-all-sky-results)).

## Quick start

```bash
cd Main
# 1. Edit the USER SETTINGS section at the top of run_3ML_master.py:
#    isotropic spectrum, Sun/Moon/Loop I files and normalizations,
#    optional components, output name.
# 2. Check that config.yaml uses the same data selection and IRFs as the all-sky fit.
# 3. Run (the full fit takes many hours):
sbatch slurm.sh            # or: python run_3ML_master.py
```

Run the script from this folder: the source model, `config.yaml` and the
comparison module are found next to the script, but the `output/` folder and
the comparison PDF are written to the current working directory.

## The model

| Component | Source name(s) | Spatial model | Spectrum | Fitted? |
|---|---|---|---|---|
| 4FGL-DR4 point and extended sources within 30° of the GC | catalog names | catalog | catalog | Fixed (point sources can be freed, see `free_point_source_radius`) |
| GALPROP π⁰ decay | `pion_decay` | GALPROP mapcube | from the mapcube (optionally × a power-law tilt) | Normalization free (+ tilt index with `fit_pi0_index`) |
| GALPROP inverse Compton | `ic` | GALPROP mapcube | from the mapcube | Normalization free |
| GALPROP bremsstrahlung | `brem` | GALPROP mapcube | from the mapcube | Fixed at 1 |
| Galactic-centre excess (`include_dm`) | `GCE_NFW` | NFW² template, `DM_Model/nfw_gamma1.fits` | cutoff power law | `K`, `index`, `xc` free |
| Fermi bubbles (`include_bubbles`) | `bubbles_low`, `bubbles_high` | `Bubbles_Model/` smoothed templates | power law, log parabola or cutoff power law, each independent | Spectral parameters free |
| Sun, Moon | `Sun`, `Moon` | mapcubes from `Sun_and_Moon/` | from the mapcube | Fixed at all-sky values |
| Loop I shells | `Loop_I_S1`, `Loop_I_S2` | mapcubes from `LoopI/` | empirical, from the mapcube | Fixed at all-sky values |
| Isotropic | `isodiff` | isotropic | all-sky `isotropic_spectrum.txt` | Fixed (normalization 1), energy dispersion off |

Spectral conventions (astromodels; energies in keV, pivot 1 GeV = 10⁶ keV):

- Cutoff power law: dN/dE = K (E/piv)^index exp(−E/xc)
- Log parabola: dN/dE = K (E/piv)^(alpha − beta ln(E/piv))
- Power law: dN/dE = K (E/piv)^index

For template sources (GCE, bubbles) the spatial template integrates to 1 and
its `K` is fixed at 1, so the spectrum is the total flux of the component over
the template's footprint (60° × 60° for the NFW map, 80° × 80° for the bubbles),
not only inside the 40° × 40° region.

The 4FGL + GALPROP part of the model is stored in `fermi_src_model.yaml`, which
the script loads. To rebuild it (for example with other GALPROP files), set
`make_src_model = True` and run once in the Fermipy environment.

## Settings

All settings are in the **USER SETTINGS** block at the top of
`run_3ML_master.py`:

| Setting | Meaning |
|---|---|
| `savename` | Label used in every output file name. |
| `make_src_model` | `True`: rebuild `fermi_src_model.yaml` from 4FGL-DR4 and the GALPROP files `pi0_file`, `ic_file`, `brem_file`. `False`: load the saved model. |
| `fit_pi0_index` | Fit a power-law tilt of the π⁰ spectrum: I(E) = K (E / 1 GeV)^index × I_GALPROP(E), with `index` free in [−0.5, 0.5] (positive = harder). Needs the `diffuse-analysis` `FermipyLike` from commit `ffad2d3` or later. |
| `include_dm`, `nfw_file` | Add the GCE component and choose its template. |
| `include_bubbles`, `bubble_spectrum` | Add the bubbles; `"powerlaw"`, `"log_parabola"` or `"cutoff_powerlaw"` (used for both components). |
| `loopI_directory`, `fixed_templates` | Files and fixed normalizations of the Sun, Moon and Loop I sources. |
| `isotropic_spectrum_file` | `isotropic_spectrum.txt` from the all-sky fit. |
| `free_point_source_radius` | Free the normalization of 4FGL point sources within this distance of the GC (default 0°: all fixed). Freed normalizations are bounded to 10⁻³–10³ times the catalog value. |
| `fermipy_outdir` | Fermipy working directory (overrides `fileio.outdir` in `config.yaml`). |

[`config.yaml`](config.yaml) is a standard Fermipy configuration for the data
selection, binning and response. The source model comes from the script.

## Carrying over the all-sky results

The fixed components must match the all-sky fit you are using:

1. **Isotropic spectrum.** Set `isotropic_spectrum_file` to
   `<outdir>/bestfit_<hash>/isotropic_spectrum.txt`. The script keeps every
   point and, because the table only spans the first to last energy-bin centre,
   adds one power-law point just below `emin` and just above `emax` (slopes from
   the two end points). The extended table is
   `output/isotropic_spectrum_extrapolated.txt`. The isotropic normalization is
   fixed at 1 and energy dispersion is disabled for it.
2. **Sun, Moon and Loop I.** Use the same files as the all-sky fit, and copy the
   fitted normalizations from `diffuse_normalizations.json` into
   `fixed_templates`. The names differ: the all-sky `LoopI_S1`/`LoopI_S2` are
   `Loop_I_S1`/`Loop_I_S2` here. The Loop I entries assume the separate-shell
   empirical model (`loopI_S1/S2_empirical_mapcube.fits`). If the all-sky fit
   used an analytic Loop I spectrum, its normalizations refer to that shape and
   cannot be applied to the empirical cubes.
3. **Data selection.** `config.yaml` must use the same event class, event type,
   IRFs, time range and zenith cut as the all-sky fit, because the isotropic
   spectrum is an effective background for that selection.
4. **GALPROP model.** Use the same GALPROP run as the all-sky fit
   (`fermi_src_model.yaml` stores the file paths).

## Output files

In `output/` (relative to the working directory):

| File | Contents |
|---|---|
| `threeML_fit_parameters_<savename>.csv` | Best-fit values and errors of all free parameters |
| `threeML_likelihood_values_<savename>.csv` | −log(likelihood) per plugin and total |
| `threeML_best_fit_model_<savename>.yml` | Best-fit 3ML model (reload with `load_model`) |
| `mcube_point_sources_<savename>.fits` | Model counts of everything except the diffuse components |
| `isotropic_spectrum_extrapolated.txt` | Isotropic table actually used (see above) |

In the Fermipy working directory (`fermipy_outdir`):

| File | Contents |
|---|---|
| `fit_roi_<savename>.*` | Fermipy ROI after the fit, with the total model map |
| residual-map files and plots | From `gta.residmap` and `make_residmap_plots` |
| `<component>_<savename>*` | Model-count map of each diffuse component (π⁰, IC, brem, isotropic, GCE, bubbles, Sun, Moon, Loop I) |
| `isodiff_<tag>.txt` | Copy of the isotropic table passed to Fermitools |

In the working directory: `GCE_Ajello2016_comparison.pdf` (when `include_dm = True`).

## GCE comparison with Ajello et al. (2016)

When `include_dm = True`, the script ends by calling `plot_gce_comparison` from
[`compare_gce_ajello2016.py`](compare_gce_ajello2016.py). It plots the
best-fit GCE spectrum and a pointwise 68% band over the four GCE envelopes of
Ajello et al. (2016), Figure 13, all converted to the same region of the sky
(default: the paper's 15° × 15° region). The paper envelopes are graphical
estimates from the published figure. The module docstring describes the
method and its caveats.

To use it elsewhere, for example on results loaded from disk:

```python
from compare_gce_ajello2016 import plot_gce_comparison
plot_gce_comparison(results, "../DM_Model/nfw_gamma1.fits",
                    comparison_roi=(-20, 20, -20, 20),     # 40 x 40 deg instead
                    output="GCE_comparison_40x40.pdf")
```

`results` is a threeML results object, such as `jl.results` or one read with
`load_analysis_results`. The script itself does not save that object; add
`jl.results.write_to("output/results.fits", overwrite=True)` after the fit if
you want to replot later without refitting.

## Notes on the custom threeML/astromodels branches

These points come from the `diffuse-analysis` branches and matter when the
model is changed or other instruments are added:

- **Mapcube sources** (`GalpropMap`: GALPROP, Sun, Moon, Loop I) take their
  amplitude from the spatial `K`. Their spectrum must be either a fixed
  `Constant(k=1)` or a `Powerlaw` tilt with `K` fixed at 1 and a fixed pivot.
  `FermipyLike` passes the FITS file to Fermitools (`MapCubeFunction`) and
  sends `K` as the normalization and the tilt as the Fermitools PowerLaw
  `Index` and `Scale` (pivot in MeV). Both use the same signed exponent.
- **Energies of mapcubes.** Fermitools reads the energies from the `ENERGIES`
  extension. The astromodels `GalpropMap.load_file` instead reads them from the
  third WCS axis, as log10(E/MeV) by default. GALPROP cubes follow that
  convention, but the Sun/Moon and Loop I cubes do not (the Loop I third axis
  is a placeholder). For this LAT-only fit this has no effect, because the LAT
  model is evaluated by Fermitools. It matters as soon as astromodels evaluates
  these sources itself, for example in the SPI/COSI joint fit or when plotting
  their spectra from the 3ML model.
- **Template sources** (`SpatialTemplate_2D`: GCE, bubbles) are passed to
  Fermitools as `SpatialMap` files with the spectrum computed by astromodels.
- **Energy dispersion.** `config.yaml` uses `edisp_bins: -1`; the all-sky fit
  uses `-2`. The value is kept to reproduce existing results. Energy
  dispersion is disabled for the isotropic component in both.

## Files in this folder

| File | Purpose |
|---|---|
| `run_3ML_master.py` | The fit: builds the model, runs the likelihood, writes results |
| `config.yaml` | Fermipy configuration (data, selection, binning, IRFs) |
| `fermi_src_model.yaml` | Saved 3ML model: 4FGL-DR4 sources within 30° of the GC and the GALPROP mapcubes |
| `compare_gce_ajello2016.py` | GCE spectrum comparison with Ajello et al. (2016), Figure 13 |
| `GCE_Ajello2016_comparison.pdf` | Example output of that comparison |
| `slurm.sh` | Example Slurm job |
