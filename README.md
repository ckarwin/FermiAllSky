# FermiAllSky

Fermi-LAT all-sky and inner-Galaxy analysis of Galactic diffuse gamma-ray
emission, built for a joint analysis with INTEGRAL/SPI (and later COSI) that
searches for positron annihilation in flight.

The repository has two kinds of folders. The **model folders** hold
templates for individual emission components, together with the scripts that
make them. They are built once and then reused, and new models can be added
to them at any time. The **analysis folders** run the fits that use those
templates, and [`Plots/`](Plots/README.md) makes figures of the results. Each
folder has its own guide:

| Folder | Type | What it contains / does | Main products |
|---|---|---|---|
| [`Sun_and_Moon/`](Sun_and_Moon/README.md) | Model | Full-sky Sun and Moon templates, made for the exact event selection of the analysis | `sun_template.fits`, `moon_template.fits` |
| [`LoopI/`](LoopI/README.md) | Model | Loop I foreground templates (Wolleben shells S1 and S2) | `loopI_*_mapcube.fits` |
| [`DM_Model/`](DM_Model/README.md) | Model | NFW (dark-matter-like) spatial templates for the Galactic-centre excess | `nfw_gamma*.fits` |
| [`Bubbles_Model/`](Bubbles_Model/README.md) | Model | Low- and high-latitude Fermi-bubble spatial templates | `bubbles_*.fits` |
| [`Isotropic/`](Isotropic/README.md) | Analysis | Masked all-sky fit with GALPROP diffuse emission, 4FGL catalog sources, Sun, Moon and Loop I: measures the isotropic spectrum and the diffuse and foreground normalizations | `isotropic_spectrum.txt`, `diffuse_normalizations.json` |
| [`Main/`](Main/README.md) | Analysis | Inner-Galaxy (40° × 40°) 3ML fit with GALPROP diffuse emission and 4FGL catalog sources, the all-sky results as fixed backgrounds, and optional Galactic-centre excess and Fermi bubbles | fitted parameters, model maps, GCE comparison plot |
| [`Plots/`](Plots/README.md) | Plotting | Count spectra, flux spectra and sky maps (data, model, residuals) of the inner-Galaxy fit, by component; comparisons with published Fermi-bubble spectra and with the official isotropic template | `counts_spectra.png`, `flux_spectra_*.png`, `residuals.png` |
| [`tests/`](tests/) | — | Regression tests (`python -m pytest -q tests`) | — |

How they connect:

1. **Build the models** you need. The NFW and bubble templates are already in
   the repository; the Loop I templates are built once with one command; the
   Sun and Moon templates are made from the event selection of the all-sky fit.
2. **Run the all-sky fit** (`Isotropic/`). Its model contains the GALPROP
   diffuse emission, the 4FGL catalog sources (held at their catalog values),
   and the Sun, Moon and Loop I models. It measures the isotropic background
   and the diffuse and foreground normalizations.
3. **Run the inner-Galaxy fit** (`Main/`). Its model contains the GALPROP
   diffuse emission and the 4FGL catalog sources, the all-sky results held
   fixed, and the GCE and bubble models.
4. **Plot the results** (`Plots/`): spectra and maps of the data, the model and
   each component, and comparisons with published bubble spectra and the
   official isotropic template.

## Typical workflow

```bash
# --- 1. Build the models ------------------------------------------------------
# Loop I (once; the large FITS files are not stored in git):
python LoopI/build_templates.py --config Isotropic/config.yaml

# Sun and Moon (Slurm array; see Sun_and_Moon/README.md). They need the event
# file selected by the all-sky preparation, so on a first pass prepare the
# all-sky fit without Sun/Moon, make the templates, then prepare again with
# them (reuse_from avoids repeating the expensive steps).
cd Sun_and_Moon
python -u generate_sun_moon_templates.py ../Isotropic/config.yaml \
    --events <outdir>/fermipy/ft1_00.fits --ltcube <outdir>/fermipy/ltcube_00.fits \
    --outdir solar_parallel --prepare-array
# ... submit the array, then run again with --finalize-array
cd ..

# NFW and bubble templates: already provided (DM_Model/, Bubbles_Model/).

# --- 2. All-sky fit ----------------------------------------------------------
cd Isotropic
python fit_isotropic.py run config.yaml

# --- 3. Inner-Galaxy fit -----------------------------------------------------
# In the USER SETTINGS of Main/run_3ML_master.py, set the isotropic spectrum
# and the fitted Sun/Moon/Loop I normalizations (diffuse_normalizations.json).
cd ../Main
sbatch slurm.sh

# --- 4. Plots ----------------------------------------------------------------
# Set savename in each script to that of the fit (see Plots/README.md).
cd ../Plots
python make_counts_plot.py
```

## Adding a new model

Put the new template (and the script that makes it) in the matching model
folder, or in a new folder, and document it in that folder's README. Then
point the analysis at it:

- **All-sky fit:** add a 3D intensity mapcube under `additional_diffuse` in
  `Isotropic/config.yaml` (see the [isotropic guide](Isotropic/README.md#additional-diffuse-components-sun-moon-others)).
- **Inner-Galaxy fit:** change the template path in the USER SETTINGS of
  `Main/run_3ML_master.py`, or add a source there (see
  [Changing the model](Main/README.md#changing-the-model)).

## Software

- Fermipy and Fermitools (with the LAT IRFs and CALDB) for every step that
  touches LAT data.
- For `Main/` and the joint fits: the `diffuse-analysis` branches of
  [threeML](https://github.com/ckarwin/threeML/tree/diffuse-analysis) and
  [astromodels](https://github.com/ckarwin/astromodels/tree/diffuse-analysis),
  which provide the `FermipyLike` options and the `GalpropMap` mapcube source
  used there. Everything else uses only Fermipy/Fermitools.
- For `Plots/`: [`lat_binned_analysis`](https://github.com/ckarwin/Fermi-LAT_Binned_Analysis).
- NumPy, SciPy, Astropy, Matplotlib and PyYAML for the utilities; pytest for the
  tests.

## Site-specific paths and files not in the repository

The configuration files and Slurm scripts contain absolute paths on the
Palmetto cluster (`/project/ckarwin/...`). Edit them for another system.

These inputs are not stored in git:

- LAT event and spacecraft files, the 4FGL catalog and its extended-source
  templates (standard Fermi data products).
- GALPROP mapcubes (π⁰, IC, bremsstrahlung).
- Sun and Moon templates (made in step 1).
- Loop I FITS cubes and maps (too large for git; rebuild with
  `LoopI/build_templates.py`, see [`LoopI/README.md`](LoopI/README.md)).
- The isotropic fit outputs referenced by `Main/run_3ML_master.py`.
