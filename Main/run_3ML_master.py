# imports:
from threeML import *
from astromodels import *
import yaml
from pathlib import Path
import fermipy
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys
import pandas as pd
import numpy as np
HERE = Path(__file__).resolve().parent
REPO = HERE.parent

pd.set_option("display.float_format", "{:.6e}".format)

# Need to specify ST version (fermipy workaround): 
fermipy.get_st_version = lambda : "2.4.0"

# Note: This only works in the fermipy env. It does not work in the joint env.
make_src_model = False
if make_src_model:
    
    print("Making fermi source model...")

    # Make 3ML source model:
    lat_catalog = FermiPySourceCatalog("4FGL-DR4")

    # Get sources within 35 degrees of GC:
    # Note: The lat_catalog object is always updated to the last query,
    # i.e, it will only contain sources within the search region. 
    gc = SkyCoord(l=0*u.deg, b=0*u.deg, frame="galactic").icrs
    ra = gc.ra.deg
    dec = gc.dec.deg
    table = lat_catalog.cone_search(ra,dec,radius=30.0)

    # Print table to terminal:
    print("List of sources in model:")
    table

    # Get  a 3ML model (a Model instance) from the table above,
    # where every source in the 4FGL becomes a Source instance.
    # Note: by default all parameters of all sources are fixed.
    model = lat_catalog.get_model()

    # Add Galactic diffuse to model:
    print("Adding Galactic diffuse to model...")

    # Define dummy spectrum:
    dummy_spectrum = Constant()
    dummy_spectrum.k.value = 1.0
    dummy_spectrum.k.fix = True

    ########################
    # pi0:
    pi0_file = "/project/ckarwin/astro/chris/annihilation_in_flight/GALPROP_Sims/Run_9_test2/results/pion_decay_mapcube_57_SA100_F98_example"
    pi0_model = GalpropMap()
    pi0_model.load_file(pi0_file)
    pi0_model.display()

    # Free overall normalization:
    pi0_model.K.value = 1.0
    pi0_model.K.free = True
    pi0_model.K.min_value = 0.0
    pi0_model.K.max_value = 100.0

    # Add as extended source
    pi0_source = ExtendedSource("pion_decay", spatial_shape=pi0_model, spectral_shape=dummy_spectrum)
    model.add_source(pi0_source)

    ########################
    # IC+Brem:
    #ic_brem_file = "/project/ckarwin/astro/chris/annihilation_in_flight/GALPROP_Sims/Run_9_test2/results/IC_plus_Bremss_mapcube_57_SA100_F98_example.fits"
    #ic_brem_model = GalpropMap()
    #ic_brem_model.load_file(ic_brem_file)
    #ic_brem_model.display()

    # Free overall normalization:
    #ic_brem_model.K.value = 1.0
    #ic_brem_model.K.free = True
    #ic_brem_model.K.min_value = 0.0
    #ic_brem_model.K.max_value = 100.0

    # Add as extended source
    #ic_brem_source = ExtendedSource("ic_brem",spatial_shape=ic_brem_model, spectral_shape=dummy_spectrum)
    #model.add_source(ic_brem_source)

    ########################
    # IC:
    ic_file = "/project/ckarwin/astro/chris/annihilation_in_flight/GALPROP_Sims/Run_9_test2/results/ics_anisotropic_mapcube_57_SA100_F98_example"
    ic_model = GalpropMap()
    ic_model.load_file(ic_file)
    ic_model.display()

    # Free overall normalization:
    ic_model.K.value = 1.0
    ic_model.K.free = True
    ic_model.K.min_value = 0.0
    ic_model.K.max_value = 100.0

    # Add as extended source
    ic_source = ExtendedSource("ic",spatial_shape=ic_model, spectral_shape=dummy_spectrum)
    model.add_source(ic_source)


    ########################
    # Brem:
    brem_file = "/project/ckarwin/astro/chris/annihilation_in_flight/GALPROP_Sims/Run_9_test2/results/bremss_mapcube_57_SA100_F98_example"
    brem_model = GalpropMap()
    brem_model.load_file(brem_file)
    brem_model.display()

    # Free overall normalization:
    brem_model.K.value = 1.0
    brem_model.K.free = True
    brem_model.K.min_value = 0.0
    brem_model.K.max_value = 100.0

    # Add as extended source
    brem_source = ExtendedSource("brem",spatial_shape=brem_model, spectral_shape=dummy_spectrum)
    model.add_source(brem_source)

    # Save model:
    model.save(str(HERE / "fermi_src_model.yaml"), overwrite=True)

    print("... finished making source model.")

# Load saved model:
if make_src_model == False:
    model = load_model(str(HERE / "fermi_src_model.yaml"))

# Optional pi0 spectral tilt: I(E) = spatial K * (E/piv)**index * I_GALPROP(E).
# Requires the accompanying FermipyLike.py. False preserves the original model.
fit_pi0_index = False
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
    model.add_source(ExtendedSource(
        "pion_decay", spatial_shape=pi0_spatial, spectral_shape=pi0_scaling,
    ))

# Add DM component:
include_dm = True
if include_dm == True:
    from astromodels import (
        ExtendedSource,
        SpatialTemplate_2D,
        Cutoff_powerlaw,
    )   

    # ------------------------------------------------------------
    # NFW annihilation morphology: fixed spatial template
    # ------------------------------------------------------------
    # Generated using make_nfw_template.py.
    # gamma=1.0 is the standard NFW case in Ajello et al. (2016).
    nfw_file = str(REPO / "DM_Model" / "nfw_gamma1.fits")

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

    dm_spectrum.piv = 1.0e6           # 1 GeV
    dm_spectrum.piv.fix = True

    # Initial guesses; all three parameters are fitted.
    dm_spectrum.K = 1.0e-12
    dm_spectrum.K.free = True

    dm_spectrum.index = -1.5
    dm_spectrum.index.bounds = (-4.0, 1.0)
    dm_spectrum.index.free = True

    dm_spectrum.xc = 3.0e6            # 3 GeV
    dm_spectrum.xc.bounds = (5.0e5, 1.0e8)  # 0.5–100 GeV
    dm_spectrum.xc.free = True

    # ------------------------------------------------------------
    # Add to your existing model
    # ------------------------------------------------------------
    dm = ExtendedSource(
        "GCE_NFW",
        spatial_shape=dm_spatial,
        spectral_shape=dm_spectrum,
    )

    # Allows this cell to be rerun without duplicating the source.
    if "GCE_NFW" in model.sources:
        model.remove_source("GCE_NFW")

    model.add_source(dm)

    print(model.GCE_NFW)

# Add Fermi Bubbles:
include_bubbles = True
bubble_spectrum = "powerlaw"  # "powerlaw", "log_parabola", or "cutoff_powerlaw"
if include_bubbles == True:
    from pathlib import Path
    from astromodels import SpatialTemplate_2D, Powerlaw, Log_parabola, Cutoff_powerlaw, ExtendedSource

    bubble_dir = REPO / "Bubbles_Model"

    def make_bubble_source(name, filename):
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
            spectrum.xc = 1.0e8            # 100 GeV in keV; initial value
            spectrum.xc.bounds = (1.0e6, 1.0e10)  # 1 GeV–10 TeV
            spectrum.xc.free = True
        else:
            raise ValueError(
                "bubble_spectrum must be 'powerlaw', 'log_parabola', or 'cutoff_powerlaw'"
            )

        spectrum.piv = 1e6            # 1 GeV in keV
        spectrum.piv.fix = True

        spectrum.K = 1e-12           # ph cm^-2 s^-1 keV^-1; initial value
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

# Fixed Sun, Moon, and Loop I intensity mapcubes.
# Update paths and amplitudes for your all-sky fit. For moved Loop I assets,
# set this to the same directory as loopI.directory in the isotropic config.
loopI_directory = REPO / "LoopI"
fixed_templates = [
    ("Sun",    "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Sun_and_Moon/solar_parallel/sun_template.fits",    1.052455442579017),  # EDIT path and norm
    ("Moon",   "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Sun_and_Moon/solar_parallel/moon_template.fits",   1.0927249143495301),  # EDIT path and norm
    ("Loop_I_S1", str(loopI_directory / "loopI_S1_empirical_mapcube.fits"), 0.028642697745440024),  # EDIT path and norm
    ("Loop_I_S2", str(loopI_directory / "loopI_S2_empirical_mapcube.fits"), 0.21487439356080634)
    ]

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

# EDIT: spectrum exported by your all-sky fit (MeV, ph cm^-2 s^-1 sr^-1 MeV^-1).
isotropic_spectrum_file = "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/isotropic_fit/output_with_solar_and_loopI_separate/bestfit_33aa0a726a/isotropic_spectrum.txt"

# Fix/free parameters:
print()
print("Setting fixed/free params...")
print()

center = SkyCoord(l=0*u.deg, b=0*u.deg, frame="galactic").icrs
radius = 0.0 * u.deg

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
model.brem.GalpropMap.K.free=False

# Display model:
print("3ML model before gta setup:")
model.display()

with (HERE / "config.yaml").open() as stream:
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

# Extend the isotropic spectrum to the LAT energy limits only if needed.
# Use the first two points for the low-energy power law, and the last two
# for the high-energy power law. Keep every original spectrum point unchanged.
table = np.loadtxt(isotropic_spectrum_file, usecols=(0, 1), ndmin=2)
if (len(table) < 2 or not np.all(np.isfinite(table))
        or np.any(table <= 0) or np.any(np.diff(table[:, 0]) <= 0)):
    raise ValueError("Isotropic table needs increasing positive energies and positive intensities")
selections = [dict(lat_config["selection"], **c.get("selection", {}))
              for c in components or [{}]]
emin = min(c["emin"] for c in selections)* 0.99
emax = max(c["emax"] for c in selections) * 1.01
extra = []
for energy in [emin, emax]:
    if energy < table[0, 0] or energy > table[-1, 0]:
        i, j = (0, 1) if energy < table[0, 0] else (-2, -1)
        slope = np.log(table[j, 1] / table[i, 1]) / np.log(table[j, 0] / table[i, 0])
        extra.append([energy, table[i, 1] * (energy / table[i, 0]) ** slope])
if extra:
    table = np.vstack([table, extra])
    table = table[np.argsort(table[:, 0])]
    Path("output").mkdir(exist_ok=True)
    isotropic_spectrum_file = "output/isotropic_spectrum_extrapolated.txt"
    np.savetxt(isotropic_spectrum_file, table, fmt="%.18e")

# Define instance of fermipy plugin:
LAT = FermipyLike(
    "LAT", lat_config,
    galactic_diffuse=["ic", "pion_decay", "brem"],
    isotropic_spectrum=str(isotropic_spectrum_file),
    skip_source_diagnostics=True,
)

# Change the output directory from the default:
LAT.configuration["fileio"]["outdir"] = "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/output"

# Look at configuration:
config = LAT.configuration
print()
print("Configuration:")
print(config)
print()

# Fermipy / NumPy 2.x fix:
np.set_printoptions(legacy="1.25")

# Setup analysis (i.e. the .setup method in fermipy)
plugins = DataList(LAT) # add other plugins in the future, e.g. COSI, SPI, etc.
jl = JointLikelihood(model, plugins)
print("Setup works!")

# For validating with simulated data:
# In this case, fit just normalizations of diffuse components.

# Access Fermipy GTAnalysis object
gta = LAT._gta

# Simulate counts from current model
#gta.simulate_roi(randomize=False)

# fix isodiff to the absolute intensity of the all-sky spectrum:
model.LAT_isodiff_Normalization.value = 1.0
model.LAT_isodiff_Normalization.free = False

# check model again:
print("3ML model after setup:")
model.display()

print("Free parameters passed to 3ML:")
for name, par in model.free_parameters.items():
    print(name, par.value, par.delta, par.min_value, par.max_value, par.fix)

# Perform the fit:
jl.set_minimizer("minuit")
fit_results, likelihood_values = jl.fit()
print("Fit Works!")
print("Fit results:")
print(fit_results)

savename="baseline_with_all_components_pl_bubbles"

# Save 3ml results:
fit_results.to_csv(f"output/threeML_fit_parameters_{savename}.csv", sep="\t", index=True)
likelihood_values.to_csv(f"output/threeML_likelihood_values_{savename}.csv", sep="\t", index=True)

# Save the optimized/best-fit model
model.save(f"output/threeML_best_fit_model_{savename}.yml", overwrite=True)

# Write fermipy ROI to file
gta.write_roi(f"fit_roi_{savename}", save_model_map=True,make_plots=False)
maps = gta.residmap(f"fit_roi_{savename}", make_plots=False)

# Use map plotter for better control:
gta.plotter.make_residmap_plots(
    maps,
    roi=None,
    label_ts_threshold=1e99,
)

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
ps_map = gta.model_counts_map(exclude=["pion_decay", "ic", "brem", "isodiff"]
                             + optional_diffuse_names + fixed_template_names)
ps_map.write(f"output/mcube_point_sources_{savename}.fits",overwrite=True)

# Compare best-fit DM spectrum:
if include_dm == True:

    # Paste this entire cell AFTER the maximum-likelihood fit.
    # Requires numpy, matplotlib, astropy, and your existing threeML fit results.
    # Reference: Ajello et al. 2016, ApJ 819, 44, Figure 13 (NOT Figure 18).
    import warnings
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    # ---- Settings: change only these lines as needed --------------------------
    fit_result = jl.results       # use jl.results if your fit object is jl
    nfw_path = str(REPO / "DM_Model" / "nfw_gamma1.fits")
    # Common plotting region in GALACTIC coordinates, with NO plane mask.
    # Default: the paper's 15 x 15 degree region. For 40 x 40 use (-20,20,-20,20).
    comparison_roi = (-7.5, 7.5, -7.5, 7.5)  # lmin,lmax,bmin,bmax in degrees
    boundary_subpixels = 4            # supersampling for fractional boundary pixels
    max_samples = 5000
    # -------------------------------------------------------------------------

    # Read the SAVED best fit, so later changes to the live model do not affect it.
    _fitted = fit_result.optimized_model
    _shape = _fitted.GCE_NFW.spectrum.main.shape
    _spatial = _fitted.GCE_NFW.spatial_shape
    if _spatial.K.free or not np.isclose(_spatial.K.value, 1.):
        raise ValueError("This cell assumes the spatial K is fixed at 1.")
    if _shape.piv.free:
        raise ValueError("This cell assumes a fixed spectral pivot.")

    # Integrate the actual FITS map. The native WCS may be Galactic or equatorial;
    # pixel_to_skycoord handles the conversion before testing Galactic ROI cuts.
    with fits.open(nfw_path) as _hdul:
        _hdu = next(h for h in _hdul if h.data is not None and h.data.ndim == 2)
        _map = np.asarray(_hdu.data, dtype=float).copy()
        _wcs = WCS(_hdu.header).celestial
    if not np.all(np.isfinite(_map)) or np.any(_map < 0):
        raise ValueError("Expected a finite, nonnegative NFW intensity map.")
    _ny, _nx = _map.shape
    _y, _x = np.indices(_map.shape, dtype=float)

    # Exact spherical solid angles of parent pixels, split into two triangles.
    def _unit_vectors(x, y):
        return np.moveaxis(pixel_to_skycoord(x, y, _wcs).cartesian.xyz.value, 0, -1)
    def _triangle_area(a, b, c):
        numerator = np.abs(np.sum(a * np.cross(b, c), axis=-1))
        denominator = (1 + np.sum(a*b, axis=-1) + np.sum(b*c, axis=-1)
                       + np.sum(c*a, axis=-1))
        return 2 * np.arctan2(numerator, denominator)
    _a = _unit_vectors(_x-.5, _y-.5)
    _b = _unit_vectors(_x+.5, _y-.5)
    _c = _unit_vectors(_x+.5, _y+.5)
    _d = _unit_vectors(_x-.5, _y+.5)
    _omega = _triangle_area(_a, _b, _c) + _triangle_area(_a, _c, _d)
    _total = np.sum(_map * _omega)
    if not np.isfinite(_total) or _total <= 0:
        raise ValueError("Invalid template integral.")
    if not np.isclose(_total, 1., rtol=.02):
        raise ValueError(f"FITS integral is {_total:g}, expected ~1 sr-weighted. "
                         "Check that this is the unit-integral map used in the fit.")

    def _roi_fraction(roi):
        lmin,lmax,bmin,bmax = roi
        if not (-180 <= lmin < lmax <= 180 and -90 < bmin < bmax < 90):
            raise ValueError("Use increasing Galactic bounds in [-180,180].")
        # Ensure requested comparison region is contained in the map footprint.
        ll,bb = np.meshgrid(np.linspace(lmin,lmax,121), np.linspace(bmin,bmax,121))
        xx,yy = skycoord_to_pixel(SkyCoord(ll*u.deg,bb*u.deg,frame='galactic'),_wcs)
        if np.any(~np.isfinite(xx+yy)) or np.any((xx<-.5)|(xx>_nx-.5)|
                                                           (yy<-.5)|(yy>_ny-.5)):
            raise ValueError("Requested ROI extends beyond the NFW template footprint.")
        occupancy = np.zeros_like(_map)
        offsets = (np.arange(boundary_subpixels)+.5)/boundary_subpixels-.5
        for dy in offsets:
            for dx in offsets:
                sky = pixel_to_skycoord(_x+dx,_y+dy,_wcs).galactic
                lon = sky.l.wrap_at(180*u.deg).deg
                lat = sky.b.deg
                occupancy += ((lon>=lmin)&(lon<=lmax)&(lat>=bmin)&(lat<=bmax))
        occupancy /= boundary_subpixels**2
        return np.sum(_map*_omega*occupancy)/_total

    _paper_roi = (-7.5,7.5,-7.5,7.5)
    _fpaper = _roi_fraction(_paper_roi)
    _ftarget = _fpaper if tuple(comparison_roi)==_paper_roi else _roi_fraction(comparison_roi)
    if min(_fpaper,_ftarget) <= 0:
        raise ValueError("Empty ROI integral.")
    _paper_scale = _ftarget/_fpaper
    print(f"Template flux fraction in paper ROI: {_fpaper:.5f}")
    print(f"Template flux fraction in comparison ROI: {_ftarget:.5f}")
    print(f"Multiply your full-template spectrum by {_ftarget:.5f}.")
    print(f"Multiply Figure 13 fluxes by {_paper_scale:.5f}.")

    # Correlated samples supplied by threeML: SAME row selection for every parameter.
    # For MLE these approximate the likelihood using its fitted covariance; for
    # Bayesian results they are posterior samples. No independent parameter draws.
    _pars = [_shape.K, _shape.index, _shape.xc]
    _columns = [np.asarray(fit_result.get_variates(p.path),dtype=float) if p.free
                else None for p in _pars]
    _lengths = [len(a) for a in _columns if a is not None]
    if not _lengths or len(set(_lengths)) != 1:
        raise ValueError("Need aligned threeML samples for free GCE parameters.")
    _n = _lengths[0]
    _rows = np.random.default_rng(2026).choice(_n,min(_n,max_samples),replace=False)
    _draws = np.column_stack([a[_rows] if a is not None else
                             np.full(len(_rows),p.value) for a,p in zip(_columns,_pars)])
    _valid = np.all(np.isfinite(_draws),axis=1)&(_draws[:,0]>0)&(_draws[:,2]>0)
    if not np.all(_valid):
        warnings.warn(f"Discarded {np.mean(~_valid):.1%} nonphysical draws. "
                      "If appreciable, use a profile likelihood or posterior band.")
    _draws = _draws[_valid]
    if len(_draws)<100:
        raise ValueError("Too few usable uncertainty samples.")

    _E_MeV = np.geomspace(1e3,1e5,250)
    _E_keV = 1e3*_E_MeV
    _pivot = _shape.piv.value
    _K,_index,_xc = _draws.T
    # Convert dN/dE per keV to per MeV by multiplying by 1000.
    # Thus E_MeV^2 * [1000*dN/dE_keV] gives MeV cm^-2 s^-1.
    _convert = _ftarget * 1e3 * _E_MeV**2
    _best = _convert*np.asarray(_shape(_E_keV),dtype=float)
    _expected = (_shape.K.value*(_E_keV/_pivot)**_shape.index.value
                 *np.exp(-_E_keV/_shape.xc.value))
    if not np.allclose(np.asarray(_shape(_E_keV)),_expected,rtol=1e-6,atol=0):
        raise ValueError("Spectrum is not the assumed Cutoff_powerlaw convention.")
    _samples = (_convert[None,:]*_K[:,None]
                *(_E_keV[None,:]/_pivot)**_index[:,None]
                *np.exp(-_E_keV[None,:]/_xc[:,None]))
    _low,_high = np.percentile(_samples,[16,84],axis=0)

    # Embedded estimates of Fig.13 ENVELOPE EDGES, obtained from the published
    # PDF's vector hatching. Columns: E[MeV], lower and upper E^2 dN/dE[MeV/cm2/s]
    # for its 15x15 deg ROI. Approximate graphical data, NOT author-supplied tables.
    # Hatching thickness/spacing limits accuracy to several percent (worse near
    # steep edges). Black/green stop before the lower edge leaves the plotted range;
    # no unseen tails are invented. No published best-fit center is implied.
    _fig13 = {
        "red": np.array([
            [1200, 0.00017325, 0.00021024],
            [1334.86, 0.00017897, 0.00022058],
            [1484.89, 0.00018475, 0.00023128],
            [1651.77, 0.00019055, 0.00024227],
            [1837.41, 0.00019628, 0.00025347],
            [2043.91, 0.00020143, 0.00026419],
            [2273.61, 0.00020669, 0.00027532],
            [2529.14, 0.00021173, 0.00028646],
            [2813.38, 0.00021689, 0.00029803],
            [3129.57, 0.00022146, 0.00030908],
            [3481.29, 0.00022607, 0.00032044],
            [3872.55, 0.00022945, 0.00033034],
            [4307.77, 0.00023248, 0.00033989],
            [4791.91, 0.00023558, 0.00034984],
            [5330.46, 0.00023878, 0.00036018],
            [5929.53, 0.00024085, 0.000369],
            [6595.93, 0.00024143, 0.00037569],
            [7337.23, 0.00024111, 0.00038106],
            [8161.84, 0.00023845, 0.00038272],
            [9079.13, 0.0002355, 0.00038387],
            [10099.5, 0.00023019, 0.00038108],
            [11234.6, 0.00022495, 0.00037821],
            [12497.2, 0.00021664, 0.0003699],
            [13901.7, 0.00020842, 0.0003614],
            [15464.1, 0.00019691, 0.00034675],
            [17202, 0.00018514, 0.00033107],
            [19135.3, 0.00017163, 0.00031167],
            [21285.9, 0.00015776, 0.00029093],
            [23678.1, 0.00014272, 0.00026726],
            [26339.3, 0.00012675, 0.00024101],
            [29299.5, 0.00011117, 0.00021463],
            [32592.3, 9.4749e-05, 0.00018574],
            [36255.3, 8.018e-05, 0.00015959],
            [40329.9, 6.4972e-05, 0.00013129],
            [44862.5, 5.2568e-05, 0.00010784],
            [49904.4, 4.016e-05, 8.3637e-05],
            [55513.1, 3.0596e-05, 6.4686e-05],
            [61752, 2.1867e-05, 4.6928e-05],
            [68692.1, 1.5406e-05, 3.3561e-05],
            [76412.2, 1.12e-05, 2.2543e-05],
            [85000, 8.2766e-06, 1.4649e-05],
        ]),
        "blue": np.array([
            [1200, 0.00017758, 0.00020576],
            [1334.86, 0.00018117, 0.00020992],
            [1484.89, 0.00018483, 0.00022494],
            [1651.77, 0.00018856, 0.00022948],
            [1837.41, 0.00020204, 0.0002459],
            [2043.91, 0.00020612, 0.00025085],
            [2273.61, 0.00021028, 0.00026881],
            [2529.14, 0.00021455, 0.00027423],
            [2813.38, 0.00021886, 0.00027976],
            [3129.57, 0.00022328, 0.00029978],
            [3481.29, 0.00022778, 0.00030582],
            [3872.55, 0.00023239, 0.00031199],
            [4307.77, 0.00023706, 0.00033433],
            [4791.91, 0.00024187, 0.00034105],
            [5330.46, 0.00024675, 0.00034795],
            [5929.53, 0.00025172, 0.00035496],
            [6595.93, 0.00025681, 0.00036214],
            [7337.23, 0.00024945, 0.00036945],
            [8161.84, 0.00025448, 0.00037691],
            [9079.13, 0.00024717, 0.00038452],
            [10099.5, 0.00025217, 0.00037348],
            [11234.6, 0.00024493, 0.00038103],
            [12497.2, 0.00023789, 0.0003701],
            [13901.7, 0.00023108, 0.00035947],
            [15464.1, 0.00022445, 0.00034917],
            [17202, 0.00020757, 0.00033915],
            [19135.3, 0.00019195, 0.00032943],
            [21285.9, 0.00018644, 0.00030464],
            [23678.1, 0.00016415, 0.00028172],
            [26339.3, 0.0001518, 0.00026052],
            [29299.5, 0.00013366, 0.00024093],
            [32592.3, 0.00011769, 0.00021212],
            [36255.3, 0.00010361, 0.00017782],
            [40329.9, 8.2696e-05, 0.00014906],
            [44862.5, 6.9321e-05, 0.00013124],
            [49904.4, 5.5324e-05, 0.00010474],
            [55513.1, 4.4155e-05, 7.9594e-05],
            [61752, 3.1944e-05, 6.0477e-05],
            [68692.1, 2.311e-05, 4.5957e-05],
            [76412.2, 1.5919e-05, 3.1654e-05],
            [85000, 1.0964e-05, 2.1804e-05],
        ]),
        "green": np.array([
            [1200, 0.0001075, 0.00015922],
            [1334.86, 0.00011068, 0.00018994],
            [1484.89, 0.00011968, 0.00021573],
            [1651.77, 0.00012322, 0.00024502],
            [1837.41, 0.00012685, 0.00029229],
            [2043.91, 0.0001306, 0.00031608],
            [2273.61, 0.00013446, 0.00035899],
            [2529.14, 0.00013843, 0.00040772],
            [2813.38, 0.0001357, 0.0004409],
            [3129.57, 0.0001397, 0.00047676],
            [3481.29, 0.00013694, 0.00051554],
            [3872.55, 0.00012779, 0.00055749],
            [4307.77, 0.00012526, 0.00057395],
            [4791.91, 0.00011691, 0.00059093],
            [5330.46, 0.00010387, 0.00060837],
            [5929.53, 9.694e-05, 0.00059635],
            [6595.93, 8.2007e-05, 0.00058455],
            [7337.23, 7.2865e-05, 0.00054551],
            [8161.84, 5.8687e-05, 0.00048469],
            [9079.13, 4.7271e-05, 0.00043068],
            [10099.5, 3.6247e-05, 0.00036433],
            [11234.6, 2.7795e-05, 0.0003082],
            [12497.2, 1.9321e-05, 0.00023633],
            [13901.7, 1.343e-05, 0.00017255],
            [15464.1, 1.0299e-05, 0.00011995],
            [17202, 7.8973e-06, 8.3376e-05],
            [19135.3, np.nan, 5.0017e-05],
            [21285.9, np.nan, 3.0004e-05],
            [23678.1, np.nan, 1.8e-05],
            [26339.3, np.nan, 1.0281e-05],
        ]),
        "black": np.array([
            [1200, 0.00013244, 0.00017104],
            [1334.86, 0.00014698, 0.00020347],
            [1484.89, 0.00016311, 0.00023378],
            [1651.77, 0.00017484, 0.0002686],
            [1837.41, 0.00018741, 0.00029809],
            [2043.91, 0.00019403, 0.00033081],
            [2273.61, 0.00020089, 0.00036712],
            [2529.14, 0.00020089, 0.0003801],
            [2813.38, 0.00020089, 0.00040744],
            [3129.57, 0.00019403, 0.00040744],
            [3481.29, 0.00018741, 0.00042182],
            [3872.55, 0.00016887, 0.00040744],
            [4307.77, 0.00015217, 0.0003801],
            [4791.91, 0.00012792, 0.0003425],
            [5330.46, 0.00010753, 0.00030862],
            [5929.53, 8.7312e-05, 0.00025944],
            [6595.93, 6.6137e-05, 0.00020347],
            [7337.23, 4.8388e-05, 0.00015957],
            [8161.84, 3.1901e-05, 0.00011276],
            [9079.13, 2.1031e-05, 7.9686e-05],
            [10099.5, 1.2494e-05, 4.901e-05],
            [11234.6, np.nan, 3.0144e-05],
            [12497.2, np.nan, 1.5584e-05],
            [13901.7, np.nan, 8.0579e-06],
        ]),
    }

    # Both spectra are now intrinsic fluxes in exactly the same Galactic rectangle.
    # Assuming the same standard gamma=1 NFW-squared morphology, this is a pure
    # aperture conversion. It does not refit the data or reproduce a different IEM.
    _labels = {
        'black': 'Pulsars, index-scaled',
        'green': 'OB stars, index-scaled',
        'red': 'Pulsars, intensity-scaled',
        'blue': 'OB stars, intensity-scaled',
    }
    _hatches = {'black':'---','green':'\\\\','red':'|||','blue':'///'}
    _floor = 7.1e-6*_paper_scale  # lower plotted limit in the original figure
    fig, ax = plt.subplots(figsize=(10,7),constrained_layout=True)
    # Put a border around the plot
    for spine in ax.spines.values():
        spine.set_visible(True)
    for color in ('red','blue','green','black'):
        energy,lo,hi = _fig13[color].T
        lo,hi = lo*_paper_scale,hi*_paper_scale
        # NaN = lower boundary was outside the paper's plotted range.
        # Show only its visible portion, without inventing a faint-energy tail.
        visible_lo = np.where(np.isfinite(lo),lo,_floor)
        ax.fill_between(energy,visible_lo,hi,color=color,alpha=.10)
        ax.fill_between(energy,visible_lo,hi,facecolor='none',edgecolor=color,
                        hatch=_hatches[color],linewidth=1,label=_labels[color])
        ax.plot(energy,lo,color=color,lw=.8,alpha=.6)
        ax.plot(energy,hi,color=color,lw=.8,alpha=.6)
    ax.fill_between(_E_MeV,_low,_high,color='darkorange',alpha=.35,
                    label='Your fit: pointwise 68% band',zorder=9)
    ax.plot(_E_MeV,_best,color='darkorange',lw=2.8,label='Your GCE_NFW best fit',zorder=10)
    ax.set(xscale='log',yscale='log',xlim=(1e3,1e5),
           ylim=(_floor,max(8e-4*_paper_scale,1.1*np.nanmax(_high))),
           xlabel='Energy (MeV)',
           ylabel=r'$E^2\,dN/dE\quad[\mathrm{MeV\,cm^{-2}\,s^{-1}}]$')
    _lmin,_lmax,_bmin,_bmax = comparison_roi
    ax.set_title(f'GCE comparison: l = [{_lmin:g}, {_lmax:g}] deg, '
                 f'b = [{_bmin:g}, {_bmax:g}] deg')
    ax.grid(which='major',alpha=.18)
    ax.legend(fontsize=9,loc='best')
    plt.ylim(1e-5,2e-3)
    #plt.show()
    print("Paper envelopes are graphical estimates; no published central curves "
          "or parameter covariance were reconstructed.")
    print("Your band is pointwise and includes K/index/xc correlations. For MLE "
          "it is a covariance approximation, not a profile-likelihood interval.")
    fig.savefig('GCE_Ajello2016_comparison.pdf',bbox_inches='tight')
