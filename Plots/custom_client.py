#imports
from lat_binned_analysis import Analysis
from lat_binned_analysis import Plots
from lat_binned_analysis import pyLikelihood
import pandas as pd

#define instance of class:
iteration_instance = pyLikelihood('inputs.yaml')
instance = Analysis('inputs.yaml')
instance_plots = Plots('inputs.yaml')

# Note: Need to replace fermipy data directory in xml file:
# %s/$(FERMIPY_DATA_DIR)*/\/project\/ckarwin\/astro\/software\/fermi\/lib\/python3.11\/site-packages\/fermipy\/data/g

# counts map:
#instance.counts_map_small(evfile="../output/ft1_00.fits")
#instance_plots.plot_skymap("GC_cmap_small.fits", "cmap.png", vmin=100, title="Data")

# model map:
#instance.model_map("../output/srcmap_00.fits", "../output/fit_roi_00.xml",
#        "../output/ltcube_00.fits", "../output/bexpmap_00.fits", "model_map.fits",
#        "cmap")
#instance_plots.plot_skymap("model_map.fits", "model_map.png", vmin=100, title="Model")

# residual map:
#instance.residual_map("GC_cmap_small.fits", "model_map.fits", "residuals.fits")
#instance_plots.plot_skymap("residuals.fits", "residual_map.png", plot_type="residuals")

# model cube:
#instance.model_map("../output/srcmap_00.fits", "../output/fit_model_00.xml",
#        "../output/ltcube_00.fits", "../output/bexpmap_00.fits", "model_cube.fits",
#        "ccube")

# residual map 3d:
instance.residual_map2("GC_cmap_small.fits", "model_map.fits", "../output/ccube_00.fits", "../output/mcube_fit_roi_00.fits")
instance_plots.plot_skymap("GC_0-4_fractional_residual.fits", "residauls_1.png", \
        plot_type="residuals", title="300 MeV - 1.3 GeV")
instance_plots.plot_skymap("GC_4-19_fractional_residual.fits", "residuals_2.png", \
        plot_type="residuals", title="1.3 GeV - 104 GeV")
instance_plots.plot_skymap("GC_19-27_fractional_residual.fits", "residuals_3.png", \
        plot_type="residuals", title="104 GeV - 800 GeV")
