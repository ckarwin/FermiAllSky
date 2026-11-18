# FermiAllSky

Fermi-LAT all-sky and inner-Galaxy analysis of Galactic diffuse gamma-ray
emission, built for a joint analysis with INTEGRAL/SPI (and later COSI) that
searches for positron annihilation in flight.

The analysis is done in stages. Each folder handles one stage and has its own
guide:

| Step | Folder | What it does | Main output |
|---|---|---|---|
| 1 | [`Sun_and_Moon/`](Sun_and_Moon/README.md) | Full-sky Sun and Moon templates for the exact event selection | `sun_template.fits`, `moon_template.fits` |
| 2 | [`LoopI/`](LoopI/README.md) | Loop I foreground templates (Wolleben shells S1 and S2) | `loopI_*_mapcube.fits` |
| 3 | [`Isotropic/`](Isotropic/README.md) | Masked all-sky fit: measures the isotropic spectrum and fits the GALPROP, Sun, Moon and Loop I normalizations | `isotropic_spectrum.txt`, `diffuse_normalizations.json` |
| 4 | [`Main/`](Main/README.md) | Inner-Galaxy (40° × 40°) 3ML fit using the all-sky results as fixed backgrounds, with Fermi bubbles and a Galactic-centre excess component | fitted parameters, model maps, GCE comparison plot |

Supporting folders:

| Folder | Contents |
|---|---|
| [`DM_Model/`](DM_Model/README.md) | NFW spatial templates for the Galactic-centre excess, and the script that makes them |
| [`Bubbles_Model/`](Bubbles_Model/README.md) | Low- and high-latitude Fermi-bubble spatial templates |
| [`tests/`](tests/) | Regression tests (`python -m pytest -q tests`) |

## Typical workflow

```bash
# Step 3 needs Sun/Moon templates made from its own event selection, and step 1
# needs that selection. On a first pass you can prepare the isotropic analysis
# without Sun/Moon, make the templates from <outdir>/fermipy/ft1_00.fits, and
# then prepare again with them (reuse_from saves the expensive products).

# 1. Sun and Moon templates (Slurm array; see Sun_and_Moon/README.md)
cd Sun_and_Moon
python -u generate_sun_moon_templates.py ../Isotropic/config.yaml \
    --events <outdir>/fermipy/ft1_00.fits --ltcube <outdir>/fermipy/ltcube_00.fits \
    --outdir solar_parallel --prepare-array
# ... submit the array, then run again with --finalize-array
cd ..

# 2. Loop I templates (once per checkout)
python LoopI/build_templates.py --config Isotropic/config.yaml

# 3. All-sky isotropic fit
cd Isotropic
python fit_isotropic.py run config.yaml

# 4. Inner-Galaxy fit: in the USER SETTINGS of Main/run_3ML_master.py, set the
#    isotropic spectrum path and the fitted Sun/Moon/Loop I normalizations
#    (from diffuse_normalizations.json), then run it.
cd ../Main
sbatch slurm.sh
```

## Software

- Fermipy and Fermitools (with the LAT IRFs and CALDB) for every step that
  touches LAT data.
- For `Main/` and the joint fits: the `diffuse-analysis` branches of
  [threeML](https://github.com/ckarwin/threeML/tree/diffuse-analysis) and
  [astromodels](https://github.com/ckarwin/astromodels/tree/diffuse-analysis),
  which provide the `FermipyLike` options and the `GalpropMap` mapcube source
  used there. Everything else uses only Fermipy/Fermitools.
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
