import numpy as np
import matplotlib.pyplot as plt
from fermipy.gtanalysis import GTAnalysis
from astropy.io import fits

savename = "baseline_with_all_components_pl_bubbles"

# EDIT: where the inner-Galaxy fit (Main/run_3ML_master.py) wrote its files.
# FERMIPY_OUTDIR: Fermipy working directory (`fermipy_outdir` in run_3ML_master.py):
#   counts cube, exposure and model cubes.
# RUN_OUTPUT: the output/ folder in the directory run_3ML_master.py was run from:
#   point-source model cube and 3ML results.
FERMIPY_OUTDIR = "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Main/output"
RUN_OUTPUT = "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Main/output"

# Get total counts from ccube:
ccube_file = f"{FERMIPY_OUTDIR}/ccube_00.fits"
hdu = fits.open(ccube_file)
data = hdu[0].data
total_counts = data.sum(axis=(1,2))

energy = hdu[1].data
energy_mean_list = np.sqrt(energy['E_MIN']*energy['E_MAX']) / 1000.0
err_low_list = energy_mean_list - energy['E_MIN']/1000.0
err_high_list = energy['E_MAX']/1000.0 - energy_mean_list

# Get total model counts:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_fit_roi_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
total_model_counts = data.sum(axis=(1,2))

print("model counts:")
print(total_model_counts)
print()
print("data counts:")
print(total_counts)

# Get contribution from point sources:
mcube_file = f"{RUN_OUTPUT}/mcube_point_sources_{savename}.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
total_src_counts = data.sum(axis=(1,2))

# Get contribution from isodiff:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_isodiff_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
isodiff_counts = data.sum(axis=(1,2))

# Get contribution from pi0:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_pion_decay_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
pion_decay_counts = data.sum(axis=(1,2))

# Get contribution from ic:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_ic_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
ics_counts = data.sum(axis=(1,2))

# Get contribution from Brem:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_brem_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
bremss_counts = data.sum(axis=(1,2))

# Get contribution from GCE:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_GCE_NFW_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
gce_counts = data.sum(axis=(1,2))

# Get contribution from high-lat bubbles:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_bubbles_high_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
bubble_high_counts = data.sum(axis=(1,2))

# Get contribution from low-lat bubbles:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_bubbles_low_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
bubble_low_counts = data.sum(axis=(1,2))

# Get contribution from Sun:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_Sun_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
sun_counts = data.sum(axis=(1,2))

# Get contribution from Moon:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_Moon_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
moon_counts = data.sum(axis=(1,2))

# Get contribution from Loop I S1:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_Loop_I_S1_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
LIS1_counts = data.sum(axis=(1,2))

# Get contribution from Loop I S2:
mcube_file = f"{FERMIPY_OUTDIR}/mcube_Loop_I_S2_{savename}_00.fits"
hdu = fits.open(mcube_file)
data = hdu[0].data
LIS2_counts = data.sum(axis=(1,2))

# Create subplots:
fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(9, 9),
    sharex=True,
    gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.05}
)


# plot observed counts:
xerr = np.array([err_low_list,err_high_list])
yerr = np.sqrt(total_counts)
ax1.loglog(energy_mean_list, total_counts, marker='s', color='black', ls='', label='total observed')
ax1.errorbar(energy_mean_list, total_counts, xerr=xerr, yerr=yerr, marker='s', color='black', ls='', label='_nolabel_')

# plot model counts:
ax1.loglog(energy_mean_list,total_model_counts, ls='-', lw=2, color='black', label='total predicted')

# plot galdiff and iso:
ax1.loglog(energy_mean_list,bremss_counts, ls='--', lw=2, color='purple', label='bremss')
ax1.loglog(energy_mean_list,ics_counts, ls='--', lw=2, color='cyan', label='ics')
ax1.loglog(energy_mean_list,pion_decay_counts, ls='--', lw=2, color='cornflowerblue', label='pion decay')
ax1.loglog(energy_mean_list,isodiff_counts, ls='-.', lw=2, color='orange', label='isodiff')

# plot source counts:
ax1.loglog(energy_mean_list,total_src_counts, ls=':', lw=2, color='green', label='point sources')

# plot GCE:
ax1.loglog(energy_mean_list,gce_counts, ls='-', marker="D", lw=2, color='purple', label='GCE')

# bubble high:
ax1.loglog(energy_mean_list,bubble_high_counts, ls='-',marker="o", lw=2, color='red', label='Bubble high')

# bubble low:
ax1.loglog(energy_mean_list,bubble_low_counts, ls='-', marker="o", lw=2, color='magenta', label='Bubble low')

# sun:
ax1.loglog(energy_mean_list,sun_counts, ls=':', lw=2, color='red', label='Sun')

# moon:
ax1.loglog(energy_mean_list,moon_counts, ls=':', lw=2, color='pink', label='Moon')

# LoopI S1:
ax1.loglog(energy_mean_list,LIS1_counts, ls='-', marker='s', lw=2, color='tan', label='LIS1')

# LoopI S2:
ax1.loglog(energy_mean_list,LIS2_counts, ls='-', marker='s', lw=2, color='brown', label='LIS2')


# Sanity check:
# Plot sum:
#model_sum = galdiff_counts + isodiff_counts + total_src_counts
#ax1.loglog(energy_mean_list,model_sum, ls='--', lw=2, color='cyan', label='sum')

# fractional residuals:
resid = (total_counts - total_model_counts) / total_model_counts
resid_err = (np.sqrt(total_counts)/total_model_counts)

# Plot kwargs:
ax1.set_xlim(3e2,0.1e6)
ax1.set_ylim(1e2,1e7)
ax1.set_ylabel("Counts", fontsize=16)
ax1.tick_params(axis='both', labelsize=14)
ax1.legend(loc=1, ncol=2, frameon=False, fontsize=10)
ax1.grid(ls="--",color="grey",alpha=0.2)

ax2.semilogx(energy_mean_list, resid, ls='', color='black', marker='s')
ax2.errorbar(energy_mean_list, resid, xerr=xerr, yerr=resid_err, ls='', color='black', marker='s')
ax2.axhline(y=0, ls='-',color='grey',lw=1)
ax2.set_ylim(-1,1)
ax2.set_xlabel("Energy [MeV]", fontsize=16)
ax2.set_ylabel("(D - M)/M", fontsize=16)
ax2.tick_params(axis='both', labelsize=14)
ax2.grid(ls="--",color="grey",alpha=0.2)
plt.savefig("counts_spectra.png")
plt.show()
