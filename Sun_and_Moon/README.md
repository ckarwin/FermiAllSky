# Sun and Moon templates

The Sun and the Moon are bright, moving gamma-ray sources. Averaged over the
mission they form faint, extended emission along the ecliptic. This folder
makes **full-sky template mapcubes** of that emission for exactly your event
selection, using the Fermitools Solar System Tools:

- `sun_template.fits`: solar disk + inverse-Compton emission from the solar
  photon field
- `moon_template.fits`: lunar emission

Both are Galactic CAR differential-intensity cubes
(ph cm⁻² s⁻¹ sr⁻¹ MeV⁻¹, 0.5° pixels by default) with an `ENERGIES` table. They
keep their absolute intensity: in a fit, start each with one free
normalization equal to 1. They plug straight into the isotropic fit:

```yaml
# Isotropic/config.yaml
additional_diffuse:
  Sun:  {file: .../solar_parallel/sun_template.fits,  free: true, normalization: 1.0}
  Moon: {file: .../solar_parallel/moon_template.fits, free: true, normalization: 1.0}
```

## What you need

- An activated Fermitools environment (`gtltcube`, `gtltcubesun`, `gtltsumsun`,
  `gtexpcube2`, `gtexphpsun`, `gtsuntemp`) with NumPy, Astropy, Matplotlib and PyYAML.
- The **final selected event file** of your analysis (after gtselect and
  gtmktime). The isotropic preparation writes it as `<outdir>/fermipy/ft1_00.fits`.
- The spacecraft file (taken from the config unless `--spacecraft` is given).
- Optionally the ordinary livetime cube of the same selection
  (`<outdir>/fermipy/ltcube_00.fits`) to save time.
- The config: either `Isotropic/config.yaml` or a standard Fermipy YAML. The
  event type, IRFs, zenith cut, energy range, binning and `edisp_bins` are read
  from it.

The emission profiles are the FSSC baseline models `solar_profile_v2r0.fits`
and `lunar_profile_v2r0.fits`. They are downloaded into `<outdir>/profiles/` on
the first run. On a node without internet access, download them from the
[FSSC](https://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/data/latSolarTemplate/)
and pass them with `--sun-profile` and `--moon-profile`. These profiles describe
steady average emission; they do not model flares or solar-cycle changes.

## Running on a cluster (recommended)

The expensive step is `gtltcubesun`, which computes the livetime in the
moving Sun/Moon frame. The script splits it into time chunks (30 days for the
Sun, 180 days for the Moon by default) that can run as a Slurm array.

**1. Prepare the chunk list.** This writes `solar_parallel/array_plan.json` and
`solar_parallel/array_range.txt` without running any chunk:

```bash
python -u generate_sun_moon_templates.py /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/config.yaml \
    --events /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/output/fermipy/ft1_00.fits \
    --ltcube /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/output/fermipy/ltcube_00.fits \
    --outdir solar_parallel --prepare-array
```

**2. Submit one array task per chunk** (at most 1000 running at once). Activate
the Fermitools environment before submitting, and create the log folder
(Slurm does not create it):

```bash
mkdir -p palmetto_logs
sbatch --array="$(cat solar_parallel/array_range.txt)" submit_solar_chunks.sbatch solar_parallel/array_plan.json
```

Each task waits a random 0–60 s before starting and writes its own log to
`palmetto_logs/`. The script requests 4 hours and 64 GB per task; adjust the
`#SBATCH` lines (and add your account/partition) as needed. Finished chunks are
skipped if you resubmit, so a partly failed array can simply be submitted again.

**3. Combine the chunks and make the templates** once every task has finished.
This is what [`slurm.sh`](slurm.sh) runs (`sbatch slurm.sh`):

```bash
python -u generate_sun_moon_templates.py /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/config.yaml \
    --events /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/output/fermipy/ft1_00.fits \
    --ltcube /project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/output/fermipy/ltcube_00.fits \
    --outdir solar_parallel --finalize-array
```

It stops with an error if any chunk is missing or failed.

## Running in a single job

Leave out `--prepare-array`/`--finalize-array` and the script does every step
itself, one chunk after the other. This is simple but can take days for a
multi-year data set.

An interactive version of the same procedure is in
`generate_sun_moon_templates.ipynb`; the script is the maintained version
(it adds the Slurm-array mode).

## Options

| Option | Default | Meaning |
|---|---|---|
| `config` | — | Isotropic wrapper YAML or Fermipy YAML |
| `--events` | required | Final selected FT1 file |
| `--spacecraft` | config `scfile` | FT2 file or `@list` |
| `--ltcube` | computed | Ordinary livetime cube with the same GTIs and zenith cut |
| `--outdir` | `solar_lunar_templates` | Output folder; use a new one for different settings |
| `--component` | — | Component index/name for multi-component Fermipy YAMLs |
| `--sun-profile`, `--moon-profile` | FSSC v2r0 | Alternative emission-profile files |
| `--map-binsz` | 0.5 | Template pixel size in degrees (must tile the sky) |
| `--sun-chunk-days`, `--moon-chunk-days` | 30, 180 | Livetime chunk lengths |
| `--padding-bins` | \|edisp_bins\| + 2 | Extra analysis bins added below `emin` and above `emax` |

The templates cover the analysis energy range widened by `--padding-bins`
analysis bins on each side (about 94 MeV to 320 GeV for 300 MeV–100 GeV, 8 bins per
decade and `edisp_bins: -2`). This covers the energy-dispersion margin the
isotropic fit requires.

The output folder records all inputs and settings in `run_manifest.json`. If
you rerun with different settings in the same folder, the script stops; choose
a new `--outdir`. Rerunning with identical settings resumes and skips finished
steps.

## Output files (in `--outdir`)

| File | Contents |
|---|---|
| `sun_template.fits`, `moon_template.fits` | The templates |
| `sun_moon_preview.png` | Inner 60° × 60° at about 1 GeV |
| `sun_moon_spectra.png/.pdf`, `sun_moon_skymaps.png/.pdf` | All-sky spectra and Mollweide maps |
| `sun_template_spectrum.csv`, `moon_template_spectrum.csv` | All-sky integrated spectra |
| `validation.json`, `solar_diagnostics.json` | Shapes, energy ranges and units |
| `sun_moon_sources.xml` | Example Fermitools XML model with both templates |
| `ordinary_ltcube.fits`, `average_exposure.fits`, `*_ltcube.fits`, `*_moving_exposure.fits` | Intermediate products |
| `chunks/`, `logs/`, `profiles/`, `pfiles/` | Livetime chunks, tool logs, emission profiles, private parameter files |

## Files in this folder

| File | Purpose |
|---|---|
| `generate_sun_moon_templates.py` | Main script (all steps, plus array preparation/finalization) |
| `solar_array.py` | Worker run by each Slurm array task |
| `submit_solar_chunks.sbatch` | Slurm array job script |
| `slurm.sh` | Slurm job for the final `--finalize-array` step |
| `solar_diagnostics.py` | Diagnostic plots; also runnable on an existing output folder |
| `fit_isotropic.py` | Link to `Isotropic/fit_isotropic.py` (lets the generator read the wrapper config) |
| `generate_sun_moon_templates.ipynb` | Interactive notebook version (single job) |
