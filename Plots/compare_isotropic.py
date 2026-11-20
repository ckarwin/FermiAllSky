import numpy as np
import matplotlib.pyplot as plt



# Your best-fit spectrum (Sun and Moon + Loop I).
E, I = np.loadtxt("/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/output_with_solar_and_loopI/bestfit_2eb3f75ed5/isotropic_spectrum.txt", unpack=True)


# Your best-fit spectrum (Sun and Moon).
fit_0_path = "/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/output_with_solar/bestfit_34e87a564c/"
E_0, I_0 = np.loadtxt(fit_0_path+"isotropic_spectrum.txt", unpack=True)

# Your best-fit spectrum (Sun and Moon + Loop I separate).
E_1, I_1 = np.loadtxt("/project/ckarwin/astro/chris/annihilation_in_flight/build_analysis_pipeline/Run_19/ThreeML_Analysis/Isotropic/output_with_solar_and_loopI_separate/bestfit_33aa0a726a/isotropic_spectrum.txt", unpack=True)

# Official isotropic spectrum.
E_ref, I_ref, dI_ref = np.loadtxt(
    "/project/ckarwin/astro/chris/my_envs/joint_analysis_main/share/fermitools/refdata/fermi/galdiffuse/iso_P8R3_SOURCE_V3_v1.txt", unpack=True
)

fig, ax = plt.subplots(figsize=(7, 5))

ax.loglog(
    E, E**2 * I,
    "o-", color="magenta", label="Sun and Moon + Loop I"
)

ax.loglog(
    E_0, E_0**2 * I_0,
    "o-", color="cyan", label="Sun and Moon"
)

ax.loglog(
    E_1, E_1**2 * I_1,
    "o-", color="orange", label="Sun and Moon + Loop I (separate)"
)

ax.loglog(
    E_ref, E_ref**2 * I_ref,
    color="black", label="Official isotropic"
)

ax.fill_between(
    E_ref,
    E_ref**2 * np.maximum(I_ref - dI_ref, np.finfo(float).tiny),
    E_ref**2 * (I_ref + dI_ref),
    color="black", alpha=0.15, label=r"Official $1\sigma$"
)

ax.set_xlim(E.min(), E.max())
ax.set_xlabel("Energy [MeV]")
ax.set_ylabel(
    r"$E^2 I(E)$ [MeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$]"
)
ax.grid(True, which="both", alpha=0.25)
ax.legend()
fig.tight_layout()
plt.show()
