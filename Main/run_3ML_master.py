"""Inner-Galaxy Fermi-LAT fit with 3ML.

Fits a 40 x 40 deg region around the Galactic centre with the 3ML
``FermipyLike`` plugin. The model contains:

* 4FGL-DR4 sources within 30 deg of the GC (fixed at catalog values unless
  freed below),
* GALPROP pi0, IC and bremsstrahlung mapcubes (one free normalization each;
  brem fixed),
* optional Galactic-centre excess with an NFW-squared template and a cutoff
  power-law spectrum (``include_dm``),
* optional low- and high-latitude Fermi bubbles, each with its own spectrum
  (``include_bubbles``),
* Sun, Moon and Loop I (S1, S2) mapcubes FIXED at the all-sky fit values,
* the isotropic spectrum measured by the all-sky fit, FIXED, with energy
  dispersion disabled.

Requires the ``diffuse-analysis`` branches of threeML and astromodels
(custom ``FermipyLike`` and ``GalpropMap``). Run from the folder where the
outputs should go (this folder, or a folder with a copy of this script):

    python run_3ML_master.py        (or: sbatch slurm.sh)

Edit the USER SETTINGS section below. See README.md for details.
"""

# Imports. Keep the two star imports first, so the explicit imports below
# (np, u, Path, ...) are not shadowed by names exported by threeML/astromodels.
from threeML import *
from astromodels import *
import yaml
from pathlib import Path
import fermipy
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
import numpy as np

pd.set_option("display.float_format", "{:.6e}".format)

# Need to specify ST version (fermipy workaround):
fermipy.get_st_version = lambda: "2.4.0"

# =============================================================================
# USER SETTINGS
# =============================================================================

# --- Locations --------------------------------------------------------------------
# REPO: your clone of the FermiAllSky repository. The model folders
# (DM_Model/, LoopI/, Bubbles_Model/) and Main/compare_gce_ajello2016.py are
# read from here.
REPO = Path("/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis")  # EDIT

# RUN_DIR: folder with this run's config.yaml and fermi_src_model.yaml. The
# default is the folder this script is in, so a copy of the script uses the
# files next to it. Set an absolute path to use files from elsewhere.
RUN_DIR = Path(__file__).resolve().parent

# Name used in all output files of this run.
savename = "baseline_with_all_components_pl_bubbles"

# --- Source model --------------------------------------------------------------
# False: load the saved 3ML model fermi_src_model.yaml (4FGL + GALPROP).
# True: rebuild and save it first (only works in the fermipy env, not the joint env).
make_src_model = False

# GALPROP mapcubes; used only when make_src_model = True. The saved model
# already points to these files.
pi0_file = "/project/ckarwin/astro/chris/annihilation_in_flight/GALPROP_Sims/Run_9_test2/results/pion_decay_mapcube_57_SA100_F98_example"
ic_file = "/project/ckarwin/astro/chris/annihilation_in_flight/GALPROP_Sims/Run_9_test2/results/ics_anisotropic_mapcube_57_SA100_F98_example"
brem_file = "/project/ckarwin/astro/chris/annihilation_in_flight/GALPROP_Sims/Run_9_test2/results/bremss_mapcube_57_SA100_F98_example"

# Optional pi0 spectral tilt: I(E) = K * (E/1 GeV)**index * I_GALPROP(E),
# with the index free in [-0.5, 0.5] (positive = harder).
fit_pi0_index = False

# Galactic-centre excess: NFW-squared template (DM_Model/) with a cutoff power law.
include_dm = True
nfw_file = str(REPO / "DM_Model" / "nfw_gamma1.fits")

# Fermi bubbles (Bubbles_Model/): spectrum of EACH bubble component.
include_bubbles = True
bubble_spectrum = "powerlaw"  # "powerlaw", "log_parabola", "cutoff_powerlaw" or "broken_powerlaw"

# --- Results of the all-sky fit (Isotropic/) -----------------------------------
# Fixed Sun, Moon and Loop I mapcubes: (source name, file, normalization).
# Use the files and fitted normalizations from the all-sky fit
# (diffuse_normalizations.json). For moved Loop I assets, set loopI_directory
# to the same folder as loopI.directory in the isotropic config.
loopI_directory = REPO / "LoopI"
fixed_templates = [
    ("Sun",    "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Sun_and_Moon/solar_parallel/sun_template.fits",    1.052455442579017),  # EDIT path and norm
    ("Moon",   "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Sun_and_Moon/solar_parallel/moon_template.fits",   1.0927249143495301),  # EDIT path and norm
    ("Loop_I_S1", str(loopI_directory / "loopI_S1_empirical_mapcube.fits"), 0.028642697745440024),  # EDIT path and norm
    ("Loop_I_S2", str(loopI_directory / "loopI_S2_empirical_mapcube.fits"), 0.21487439356080634),
]  # fmt: skip

# Isotropic spectrum exported by the all-sky fit
# (MeV, ph cm^-2 s^-1 sr^-1 MeV^-1).
isotropic_spectrum_file = "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/isotropic_fit/output_with_solar_and_loopI_separate/bestfit_33aa0a726a/isotropic_spectrum.txt"

# --- Fit options ------------------------------------------------------------------
# Free the normalization of every 4FGL point source within this distance of
# the GC (0 deg: all point sources stay fixed).
free_point_source_radius = 0.0 * u.deg

# Fermipy working directory (overrides fileio.outdir in config.yaml).
fermipy_outdir = "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/output"

# =============================================================================
# 1. Source model: 4FGL + GALPROP
# =============================================================================

if make_src_model:
    print("Making fermi source model...")

    # Make 3ML source model:
    lat_catalog = FermiPySourceCatalog("4FGL-DR4")

    # Get sources within 30 degrees of GC:
    # Note: The lat_catalog object is always updated to the last query,
    # i.e, it will only contain sources within the search region.
    gc = SkyCoord(l=0 * u.deg, b=0 * u.deg, frame="galactic").icrs
    ra = gc.ra.deg
    dec = gc.dec.deg
    table = lat_catalog.cone_search(ra, dec, radius=30.0)

    # Print table to terminal:
    print("List of sources in model:")
    table

    # Get a 3ML model (a Model instance) from the table above,
    # where every source in the 4FGL becomes a Source instance.
    # Note: by default all parameters of all sources are fixed.
    model = lat_catalog.get_model()

    # Add Galactic diffuse to model:
    print("Adding Galactic diffuse to model...")

    # Mapcube sources carry their own spectrum: FermipyLike requires a fixed
    # Constant(k=1) spectrum and fits the spatial normalization K.
    dummy_spectrum = Constant()
    dummy_spectrum.k.value = 1.0
    dummy_spectrum.k.fix = True

    ########################
    # pi0:
    pi0_model = GalpropMap()
    pi0_model.load_file(pi0_file)
    pi0_model.display()

    # Free overall normalization:
    pi0_model.K.value = 1.0
    pi0_model.K.free = True
    pi0_model.K.min_value = 0.0
    pi0_model.K.max_value = 100.0

    # Add as extended source
    pi0_source = ExtendedSource(
        "pion_decay", spatial_shape=pi0_model, spectral_shape=dummy_spectrum
    )
    model.add_source(pi0_source)

    ########################
    # IC:
    ic_model = GalpropMap()
    ic_model.load_file(ic_file)
    ic_model.display()

    # Free overall normalization:
    ic_model.K.value = 1.0
    ic_model.K.free = True
    ic_model.K.min_value = 0.0
    ic_model.K.max_value = 100.0

    # Add as extended source
    ic_source = ExtendedSource("ic", spatial_shape=ic_model, spectral_shape=dummy_spectrum)
    model.add_source(ic_source)

    ########################
    # Brem:
    brem_model = GalpropMap()
    brem_model.load_file(brem_file)
    brem_model.display()

    # Free overall normalization (fixed again in section 4):
    brem_model.K.value = 1.0
    brem_model.K.free = True
    brem_model.K.min_value = 0.0
    brem_model.K.max_value = 100.0

    # Add as extended source
    brem_source = ExtendedSource("brem", spatial_shape=brem_model, spectral_shape=dummy_spectrum)
    model.add_source(brem_source)

    # Save model:
    model.save(str(RUN_DIR / "fermi_src_model.yaml"), overwrite=True)

    print("... finished making source model.")

# Load saved model:
if make_src_model == False:
    model = load_model(str(RUN_DIR / "fermi_src_model.yaml"))

# Optional pi0 spectral tilt (see USER SETTINGS). The mapcube amplitude stays
# the spatial K; the spectral Powerlaw is a dimensionless multiplier with K fixed
# at 1 and a fixed pivot, as FermipyLike requires. FermipyLike sends it to
# Fermitools as the PowerLaw Index and Scale of the MapCubeFunction source.
if fit_pi0_index:
    pi0_spatial = model.pion_decay.spatial_shape
    pi0_spatial.K.free = True
    pi0_scaling = Powerlaw()
    pi0_scaling.K = 1.0
    pi0_scaling.K.fix = True  # Keep only spatial K as the free normalization.
    pi0_scaling.piv = 1.0e6  # 1 GeV in keV
    pi0_scaling.piv.fix = True
    pi0_scaling.index = 0.0  # Positive hardens; negative softens the mapcube.
    pi0_scaling.index.bounds = (-0.5, 0.5)
    pi0_scaling.index.free = True
    model.remove_source("pion_decay")
    model.add_source(
        ExtendedSource(
            "pion_decay",
            spatial_shape=pi0_spatial,
            spectral_shape=pi0_scaling,
        )
    )

# =============================================================================
# 2. Optional components: GCE (NFW) and Fermi bubbles
# =============================================================================

if include_dm:
    # ------------------------------------------------------------
    # NFW annihilation morphology: fixed spatial template
    # ------------------------------------------------------------
    # Generated using DM_Model/make_nfw_template.py.
    # gamma=1.0 is the standard NFW case in Ajello et al. (2016).
    dm_spatial = SpatialTemplate_2D(
        fits_file=nfw_file,
        frame="icrs",
    )

    # The template integrates to 1 over its full footprint.
    # Fix its normalization; fit the spectral amplitude instead.
    dm_spatial.K = 1.0
    dm_spatial.K.fix = True

    # ------------------------------------------------------------
    # Empirical spectrum: exponentially cutoff power law
    #
    # dN/dE = K * (E / piv)**index * exp(-E / xc)
    #
    # Energies: keV
    # K: ph cm^-2 s^-1 keV^-1, integrated over the template
    # ------------------------------------------------------------
    dm_spectrum = Cutoff_powerlaw()

    dm_spectrum.piv = 1.0e6  # 1 GeV
    dm_spectrum.piv.fix = True

    # Initial guesses; all three parameters are fitted.
    dm_spectrum.K = 1.0e-12
    dm_spectrum.K.free = True

    dm_spectrum.index = -1.5
    dm_spectrum.index.bounds = (-4.0, 1.0)
    dm_spectrum.index.free = True

    dm_spectrum.xc = 3.0e6  # 3 GeV
    dm_spectrum.xc.bounds = (5.0e5, 1.0e8)  # 0.5–100 GeV
    dm_spectrum.xc.free = True

    dm = ExtendedSource(
        "GCE_NFW",
        spatial_shape=dm_spatial,
        spectral_shape=dm_spectrum,
    )

    # Replace a GCE_NFW source if the loaded model already contains one.
    if "GCE_NFW" in model.sources:
        model.remove_source("GCE_NFW")

    model.add_source(dm)

    print(model.GCE_NFW)

if include_bubbles:
    bubble_dir = REPO / "Bubbles_Model"

    def make_bubble_source(name, filename):
        """Bubble component: fixed unit-integral template + free spectrum."""
        template = bubble_dir / filename
        if not template.is_file():
            raise FileNotFoundError(template)

        # Fixed spatial template, normalized to unit solid-angle integral.
        spatial = SpatialTemplate_2D(
            fits_file=str(template),
            frame="icrs",
        )
        spatial.K = 1.0
        spatial.K.fix = True

        # Independent spectral parameters for each bubble component.
        if bubble_spectrum == "powerlaw":
            spectrum = Powerlaw()
            spectrum.index = -1.9
            spectrum.index.bounds = (-4.0, -0.5)
            spectrum.index.free = True
        elif bubble_spectrum == "log_parabola":
            # dN/dE = K * (E/piv)**(alpha - beta*ln(E/piv)).
            spectrum = Log_parabola()
            spectrum.alpha = -1.9
            spectrum.alpha.bounds = (-4.0, -0.5)
            spectrum.alpha.free = True
            spectrum.beta = 0.05
            spectrum.beta.bounds = (0.0, 1.0)  # beta=0 recovers a power law
            spectrum.beta.free = True
        elif bubble_spectrum == "cutoff_powerlaw":
            # dN/dE = K * (E/piv)**index * exp(-E/xc).
            spectrum = Cutoff_powerlaw()
            spectrum.index = -1.9
            spectrum.index.bounds = (-4.0, -0.5)
            spectrum.index.free = True
            spectrum.xc = 1.0e8  # 100 GeV in keV; initial value
            spectrum.xc.bounds = (1.0e6, 1.0e10)  # 1 GeV–10 TeV
            spectrum.xc.free = True
        elif bubble_spectrum == "broken_powerlaw":
            # dN/dE = K * (E/piv)**alpha                          for E < xb
            #       = K * (xb/piv)**(alpha - beta) * (E/piv)**beta for E >= xb
            # (continuous at the break; K is the flux at piv on the low side).
            spectrum = Broken_powerlaw()
            spectrum.alpha = -1.9
            spectrum.alpha.bounds = (-4.0, -0.5)
            spectrum.alpha.free = True
            spectrum.beta = -2.5
            spectrum.beta.bounds = (-5.0, -0.5)
            spectrum.beta.free = True
            spectrum.xb = 1.0e7  # 10 GeV in keV; initial value
            spectrum.xb.bounds = (1.0e6, 1.0e8)  # 1-100 GeV
            spectrum.xb.free = True
        else:
            raise ValueError(
                "bubble_spectrum must be 'powerlaw', 'log_parabola', "
                "'cutoff_powerlaw' or 'broken_powerlaw'"
            )

        spectrum.piv = 1e6  # 1 GeV in keV
        spectrum.piv.fix = True

        spectrum.K = 1e-12  # ph cm^-2 s^-1 keV^-1; initial value
        spectrum.K.bounds = (1e-20, 1e-8)
        spectrum.K.free = True

        return ExtendedSource(
            name,
            spatial_shape=spatial,
            spectral_shape=spectrum,
        )

    bubbles_low = make_bubble_source(
        "bubbles_low",
        "bubbles_lowlat_ack2017_smooth.fits",
    )

    bubbles_high = make_bubble_source(
        "bubbles_high",
        "bubbles_highlat_ack2017_smooth.fits",
    )

    model.add_source(bubbles_low)
    model.add_source(bubbles_high)

# =============================================================================
# 3. Fixed foregrounds from the all-sky fit: Sun, Moon, Loop I
# =============================================================================
# Mapcubes in ph cm^-2 s^-1 sr^-1 MeV^-1. FermipyLike passes the file to
# Fermitools (MapCubeFunction, energies from the ENERGIES extension) and only
# sends it the normalization K.

fixed_template_names = []
for name, template_file, normalization in fixed_templates:
    spatial = GalpropMap()
    spatial.load_file(template_file)
    spatial.K.value = normalization
    spatial.K.fix = True

    spectrum = Constant()
    spectrum.k.value = 1.0
    spectrum.k.fix = True
    model.add_source(ExtendedSource(name, spatial_shape=spatial, spectral_shape=spectrum))
    fixed_template_names.append(name)

# =============================================================================
# 4. Free and fixed parameters
# =============================================================================

print()
print("Setting fixed/free params...")
print()

center = SkyCoord(l=0 * u.deg, b=0 * u.deg, frame="galactic").icrs
radius = free_point_source_radius

for src_name, src in model.point_sources.items():
    src_coord = SkyCoord(
        ra=src.position.ra.value * u.deg,
        dec=src.position.dec.value * u.deg,
        frame="icrs",
    )

    if src_coord.separation(center) <= radius:
        spectral_shape = src.spectrum.main.shape

        for par_name, par in spectral_shape.parameters.items():
            if par.is_normalization:
                val = par.value

                if val > 0:
                    par.min_value = val / 1e3
                    par.max_value = val * 1e3
                    par.delta = 0.1 * val

                par.free = True
                print("Freed:", src_name, par_name)

# Freeze Brem component:
model.brem.GalpropMap.K.free = False

# Display model:
print("3ML model before gta setup:")
model.display()

# =============================================================================
# 5. LAT configuration and isotropic spectrum
# =============================================================================

with (RUN_DIR / "config.yaml").open() as stream:
    lat_config = yaml.safe_load(stream)

# The fitted isotropic is an effective background: do not apply edisp again.
# Preserve the existing edisp settings for all other sources.
configurations = [lat_config]
components = lat_config.get("components") or []
if isinstance(components, dict):
    components = list(components.values())
configurations.extend(components)
for configuration in configurations:
    gtlike = configuration.setdefault("gtlike", {})
    disabled = list(gtlike.get("edisp_disable") or [])
    if "isodiff" not in disabled:
        disabled.append("isodiff")
    gtlike["edisp_disable"] = disabled

# The all-sky table runs from the first to the last energy-bin CENTRE.
# Extend it to the LAT energy limits only if needed. Use the first two
# points for the low-energy power law, and the last two for the high-energy
# power law. Keep every original spectrum point unchanged.
table = np.loadtxt(isotropic_spectrum_file, usecols=(0, 1), ndmin=2)
if (
    len(table) < 2
    or not np.all(np.isfinite(table))
    or np.any(table <= 0)
    or np.any(np.diff(table[:, 0]) <= 0)
):
    raise ValueError("Isotropic table needs increasing positive energies and positive intensities")
selections = [dict(lat_config["selection"], **c.get("selection", {})) for c in components or [{}]]
emin = min(c["emin"] for c in selections) * 0.99
emax = max(c["emax"] for c in selections) * 1.01
extra = []
for energy in [emin, emax]:
    if energy < table[0, 0] or energy > table[-1, 0]:
        i, j = (0, 1) if energy < table[0, 0] else (-2, -1)
        slope = np.log(table[j, 1] / table[i, 1]) / np.log(table[j, 0] / table[i, 0])
        extra.append([energy, table[i, 1] * (energy / table[i, 0]) ** slope])
Path("output").mkdir(exist_ok=True)
if extra:
    table = np.vstack([table, extra])
    table = table[np.argsort(table[:, 0])]
    isotropic_spectrum_file = "output/isotropic_spectrum_extrapolated.txt"
    np.savetxt(isotropic_spectrum_file, table, fmt="%.18e")

# =============================================================================
# 6. Likelihood setup
# =============================================================================

# Define instance of fermipy plugin. GALPROP sources replace the standard LAT
# Galactic model; the table above replaces the standard isotropic template.
LAT = FermipyLike(
    "LAT",
    lat_config,
    galactic_diffuse=["ic", "pion_decay", "brem"],
    isotropic_spectrum=str(isotropic_spectrum_file),
    skip_source_diagnostics=True,
)

# Change the output directory from the default:
LAT.configuration["fileio"]["outdir"] = fermipy_outdir

# Look at configuration:
config = LAT.configuration
print()
print("Configuration:")
print(config)
print()

# Fermipy / NumPy 2.x fix:
np.set_printoptions(legacy="1.25")

# Setup analysis (i.e. the .setup method in fermipy)
plugins = DataList(LAT)  # add other plugins in the future, e.g. COSI, SPI, etc.
jl = JointLikelihood(model, plugins)
print("Setup works!")

# Access Fermipy GTAnalysis object
gta = LAT._gta

# For validating with simulated data (fit just the diffuse normalizations),
# simulate counts from the current model:
# gta.simulate_roi(randomize=False)

# Fix isodiff to the absolute intensity of the all-sky spectrum:
model.LAT_isodiff_Normalization.value = 1.0
model.LAT_isodiff_Normalization.free = False

# check model again:
print("3ML model after setup:")
model.display()

print("Free parameters passed to 3ML:")
for name, par in model.free_parameters.items():
    print(name, par.value, par.delta, par.min_value, par.max_value, par.fix)

# =============================================================================
# 7. Fit
# =============================================================================

jl.set_minimizer("minuit")
fit_results, likelihood_values = jl.fit()
print("Fit Works!")
print("Fit results:")
print(fit_results)

# =============================================================================
# 8. Save results, maps and plots
# =============================================================================

# Save 3ml results:
fit_results.to_csv(f"output/threeML_fit_parameters_{savename}.csv", sep="\t", index=True)
likelihood_values.to_csv(f"output/threeML_likelihood_values_{savename}.csv", sep="\t", index=True)

# Save the optimized/best-fit model
model.save(f"output/threeML_best_fit_model_{savename}.yml", overwrite=True)

# Write fermipy ROI to file (in the Fermipy working directory)
gta.write_roi(f"fit_roi_{savename}", save_model_map=True, make_plots=False)
maps = gta.residmap(f"fit_roi_{savename}", make_plots=False)

# Use map plotter for better control:
gta.plotter.make_residmap_plots(
    maps,
    roi=None,
    label_ts_threshold=1e99,
)

# Model-count maps of each diffuse component:
gta.write_model_map(f"pion_decay_{savename}", name="pion_decay")
gta.write_model_map(f"ic_{savename}", name="ic")
gta.write_model_map(f"brem_{savename}", name="brem")
gta.write_model_map(f"isodiff_{savename}", name="isodiff")
optional_diffuse_names = []
if include_dm:
    optional_diffuse_names.append("GCE_NFW")
if include_bubbles:
    optional_diffuse_names.extend(["bubbles_low", "bubbles_high"])
for name in optional_diffuse_names:
    gta.write_model_map(f"{name}_{savename}", name=name)

for name in fixed_template_names:
    gta.write_model_map(f"{name}_{savename}", name=name)

# Total model excluding diffuse components
ps_map = gta.model_counts_map(
    exclude=["pion_decay", "ic", "brem", "isodiff"] + optional_diffuse_names + fixed_template_names
)
ps_map.write(f"output/mcube_point_sources_{savename}.fits", overwrite=True)

# Compare best-fit GCE spectrum with Ajello et al. (2016), Figure 13
# (writes GCE_Ajello2016_comparison.pdf; see compare_gce_ajello2016.py):
if include_dm:
    import sys

    sys.path.insert(0, str(REPO / "Main"))  # works for copies of this script too
    from compare_gce_ajello2016 import plot_gce_comparison

    plot_gce_comparison(jl.results, nfw_file)
