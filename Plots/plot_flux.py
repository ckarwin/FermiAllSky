import matplotlib.pyplot as plt
from fermipy.gtanalysis import GTAnalysis
import pyLikelihood
from BinnedAnalysis import *
import pandas as pd
import math
import numpy as np
from lat_binned_analysis import Analysis, Plots
from astropy.io import fits
import sys

savename = "baseline_with_all_components_pl_bubbles"

# EDIT: output folder of the inner-Galaxy fit (Main/run_3ML_master.py), i.e. its
# fermipy_outdir, which also holds the point-source cube and 3ML results.
OUTPUT_DIR = "../Main/output"

instance = Analysis('inputs.yaml')
instance_plots = Plots('inputs.yaml')

# calc total observed flux:
#instance_plots.calculate_flux(f'{OUTPUT_DIR}/ccube_00.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"observed")

# calc total model flux:
#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_fit_roi_{savename}_00.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"model",filetype="fermipy_mcube")

# calc isodiff flux:
#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_isodiff_{savename}_00.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"isodiff",filetype="fermipy_mcube")

# calc galdiff flux:
#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_pion_decay_{savename}_00.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"pion_decay",filetype="fermipy_mcube")

#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_ic_{savename}_00.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"ics",filetype="fermipy_mcube")

#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_brem_{savename}_00.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"bremss",filetype="fermipy_mcube")

# point source flux:
#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_point_sources_{savename}.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"pcs",filetype="fermipy_mcube")

# GCE:
#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_GCE_NFW_{savename}.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"gce",filetype="fermipy_mcube")

# Bubbles low:
instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_bubbles_low_{savename}.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"bubbles_low",filetype="fermipy_mcube")

# Bubbles high:
instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_bubbles_high_{savename}.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"bubbles_high",filetype="fermipy_mcube")
sys.exit()

# Sun:
#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_Sun_{savename}.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"sun",filetype="fermipy_mcube")

# Moon:
#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_Moon_{savename}.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"moon",filetype="fermipy_mcube")

# Loop I S1:
#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_Loop_I_S1_{savename}.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"loopI_SI",filetype="fermipy_mcube")

# Loop I S2:
#instance_plots.calculate_flux(f'{OUTPUT_DIR}/mcube_Loop_I_S2_{savename}.fits',f'{OUTPUT_DIR}/bexpmap_00.fits',"loopI_S2",filetype="fermipy_mcube")

# Get total counts from ccube:
ccube_file = f"{OUTPUT_DIR}/ccube_00.fits"
hdu = fits.open(ccube_file)
data = hdu[0].data
total_counts = data.sum(axis=(1,2))

energy = hdu[1].data
energy_mean_list = np.sqrt(energy['E_MIN']*energy['E_MAX']) / 1000.0
err_low_list = energy_mean_list - energy['E_MIN']/1000.0
err_high_list = energy['E_MAX']/1000.0 - energy_mean_list
xerr = np.array([err_low_list,err_high_list])

# Get total model counts:
mcube_file = f"{OUTPUT_DIR}/mcube_fit_roi_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
total_model_counts = data.sum(axis=(1,2))

# Create subplots:
fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(9, 9),
    sharex=True,
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
)

df = pd.read_csv("observed_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="black",ls="",marker="s", label="total observed")
ax1.errorbar(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],xerr=xerr,color="black",ls="",marker="s", label="_NONE_")

# model:
df = pd.read_csv("model_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="black",ls="-",label="total model")

# galdiff:
df = pd.read_csv("pion_decay_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="cornflowerblue",ls="--",label="pion_decay")

df = pd.read_csv("ics_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="cyan",ls="--",label="ics")

df = pd.read_csv("bremss_flux.dat", skiprows=2, sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="purple",ls="--",label="bremss")

# isodiff:
df = pd.read_csv("isodiff_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="darkorange",ls="--",label="isodiff")

# point sources:
df = pd.read_csv("pcs_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="green",ls="--",label="point sources")

# GCE:
df = pd.read_csv("gce_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="purple",ls="-",marker="D",label="GCE")

# bubbles low:
df = pd.read_csv("bubbles_low_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="magenta",marker='o', ls="-",label="Bubbles low")

# bubbles high:
df = pd.read_csv("bubbles_high_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="red",marker='o',ls="-",label="Bubbles high")

# Sun:
df = pd.read_csv("sun_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="red",ls=":",label="Sun")

# Moon:
df = pd.read_csv("moon_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="pink",ls=":",label="Moon")

# LoopI SI:
df = pd.read_csv("loopI_SI_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="tan",marker="s",ls="-",label="LIS1")

# LoopI S2:
df = pd.read_csv("loopI_S2_flux.dat", skiprows=2,  sep=r'\s+')
ax1.loglog(df["Energy[MeV]"],df["energy_flux[MeV/cm^2/s]"],color="brown",marker="s",ls="-",label="LIS2")

# fractional residuals:
resid = (total_counts - total_model_counts) / total_model_counts
resid_err = (np.sqrt(total_counts)/total_model_counts)

# Plot kwargs:
ax1.set_xlim(3e2,0.1e6)
ax1.set_ylim(1e-6,1e0)
ax1.set_ylabel("$\mathrm{E^2 dN/dE \  [MeV \ cm^{-2} \ s^{-1}]}$", fontsize=16)
ax1.tick_params(axis='both', labelsize=14)
ax1.legend(loc=1, ncol=3, frameon=False, fontsize=12)
ax1.grid(ls="--",color="grey",alpha=0.2)

ax2.semilogx(energy_mean_list, resid, ls='', color='black', marker='s')
ax2.errorbar(energy_mean_list, resid, xerr=xerr, yerr=resid_err, ls='', color='black', marker='s')
ax2.axhline(y=0, ls='-',color='grey',lw=1)
ax2.set_ylim(-0.1,0.1)
ax2.set_xlabel("Energy [MeV]", fontsize=16)
ax2.set_ylabel("(D - M)/M", fontsize=16)
ax2.tick_params(axis='both', labelsize=14)
ax2.grid(ls="--",color="grey",alpha=0.2)
plt.savefig(f"flux_spectra_{savename}.png")
plt.show()

