from lat_binned_analysis import Plots
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import matplotlib.pyplot as plt

savename = "baseline_with_all_components_pl_bubbles_gas_index_scaled"

# EDIT: where the inner-Galaxy fit (Main/run_3ML_master.py) wrote its files.
# FERMIPY_OUTDIR: Fermipy working directory (`fermipy_outdir` in run_3ML_master.py):
#   counts cube, exposure and model cubes.
FERMIPY_OUTDIR = "../output"

instance_plots = Plots('inputs.yaml')

ccube_file = f"{FERMIPY_OUTDIR}/ccube_00.fits"
instance_plots.plot_skymap(ccube_file, "count_map.png", input_dim=3, make_sum=True, vmin=100, title="Data")

mcube_file = f"{FERMIPY_OUTDIR}/mcube_fit_roi_{savename}_00.fits"
instance_plots.plot_skymap(mcube_file, "model_map.png", input_dim=3, make_sum=True, vmin=100, title="Model")

# plot residauls:

# Get total counts from ccube:
hdu = fits.open(ccube_file)
data = hdu[0].data
total_counts = data.sum(axis=0)
header = hdu[0].header

# Get total model counts:
hdu = fits.open(mcube_file)
data = hdu[0].data
total_model_counts = data.sum(axis=0)

# total spatial residuals:
resid = total_counts - total_model_counts

# plot:
image = resid
wcs_2d = WCS(header, naxis=2)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection=wcs_2d)

this_max = np.amax(image)
im = ax.imshow(image,origin="lower",cmap="seismic",\
            vmin=-1*this_max, vmax=this_max, interpolation="none")
colorbar_label="data - model [counts]"
cbar = plt.colorbar(im, ax=ax, pad=0.02)
cbar.set_label(colorbar_label, fontsize=12)

ax.set_xlabel("Galactic Longitude [$\circ$]", fontsize=14)
ax.set_ylabel("Galactic Latitude [$\circ$]", fontsize=14)
ax.grid(color="grey", ls=":", alpha=0.5)

ax.tick_params(axis='both', which='major', labelsize=14)

plt.tight_layout()
plt.savefig("residuals.png")
plt.show()

