# Fermi-bubble spectra vs. Ackermann et al. (2017)

`compare_bubbles.py` compares the low- and high-latitude Fermi-bubble spectra
from the inner-Galaxy fit ([`Main/`](../../Main/README.md)) with the published
spectra of Ackermann et al. (2017), ApJ 840, 43, Figure 9a
([arXiv:1704.03910](https://arxiv.org/abs/1704.03910)). Both are shown as
E² dF/dE integrated over the same region of the sky: the 40° × 40° analysis
region around the Galactic centre (|l| < 20°, |b| < 20°).

## Quick start

```bash
cd Plots/compare_bubbles
python compare_bubbles.py
```

Requirements: NumPy, pandas, Matplotlib, Astropy, healpy and PyYAML. The
published data are included in `reference/`.

After each inner-Galaxy fit, [`make_all_plots.py`](../README.md) runs this
comparison for that fit (`bubbles_comparison.png` in the fit's plot folder),
and [`compare_fits.py`](../README.md) overlays several fits.

## Inputs

**Your fit.** Two flux files made by [`../plot_flux.py`](../README.md) from the
fitted bubble model cubes, one per component (`bubbles_low_flux.dat`,
`bubbles_high_flux.dat`). Set their paths at the top of the script:

```python
LOW_LAT_FLUX_FILE = "../baseline_with_all_components_pl_bubbles/bubbles_low_flux.dat"
HIGH_LAT_FLUX_FILE = "../baseline_with_all_components_pl_bubbles/bubbles_high_flux.dat"
FIT_LABEL = "my fit"
```

Relative paths are taken from this folder; absolute paths also work. The two
`.dat` files bundled here are an example.

**Published data** (`reference/`, from the
[Fermi-LAT data release](https://www-glast.stanford.edu/pub_data/1220/GCexcess_Pass8_1704.03910_data/)
of the paper):

| File | Contents |
|---|---|
| `fig9a_SCA_bubbles_model_spectra.yaml` | Figure 9a spectra (GeV cm⁻² s⁻¹), integrated over the whole sky except the masked cores of the 200 brightest 3FGL sources |
| `fig8_bubbles_template_2compSCA.npy` | Bubble template (HEALPix, nside 128, counts), zero in the masked source cores |
| `fig2_mask_200_3FGL_sources.npy` | Point-source mask (1 = kept, 0 = masked) |
| `exposure.fits` | All-sky exposure (HEALPix, cm² s) at 27 energies |

## Outputs

Written to this folder (each run overwrites them):

| File | Contents |
|---|---|
| `comparison.png`, `comparison.pdf` | Both components: published spectra converted to the region, and your fit |
| `summary.csv` | Conversion factor, E² dF/dE at 10 GeV and energy flux over 1–80 GeV for each component |
| `converted_reference.csv` | Published spectra converted to the region, with the conversion factor at each energy |

## How the two are made comparable

**Your fit.** `plot_flux.py` divides the fitted model counts of each bubble
component by the exposure and sums over every pixel of the analysis region.
The `.dat` spectra are therefore already integrated over the region, so no
further factor is applied.

**The paper.** The published spectra cover a different area (the whole
template, without the source cores), so they are rescaled to the region:

```
factor(E) = Σ_region  T_filled(p) / exposure(p, E)
          / Σ_sky     T(p) × mask(p) / exposure(p, E)
```

- `T` is the published bubble template, a counts map, so `T / exposure`
  approximates its intensity.
- The denominator reproduces the area of the published spectra.
- The numerator covers the whole region, as your fit does. The published
  template is empty inside the masked source cores, so these holes are first
  filled with the mean of the surrounding template (within 2°,
  `HOLE_FILL_RADIUS_DEG`).
- Each component uses its own half of the template, split at |b| = 10° as in
  the paper.

At 10 GeV the factors are about 1.14 for the low-latitude and 0.28 for the
high-latitude component. Filling the holes matters for the low-latitude
component, whose brightest part lies near the Galactic centre where many
sources are masked. Without filling (`FILL_SOURCE_HOLES = False`) its factor
would be about 0.92, i.e. about 20% lower. The result changes by less than 1%
for fill radii between 1.5° and 3°.

**Energy flux.** The 1–80 GeV band flux is ∫ E² dF/dE d(ln E), with log-log
interpolation between the tabulated points. The "Total energy flux" in the
header of the `.dat` files is a sum over energy bins, not an integral, and is
not used.

## Limitations

- **Different bubble templates.** The fit uses the smoothed templates in
  [`Bubbles_Model/`](../../Bubbles_Model/README.md), with a smooth split
  between |b| = 9° and 11°; the paper uses its own template with a sharp split
  at 10°. Any difference in morphology moves flux between the components and in
  or out of the region.
- **Different models.** The two analyses use different Galactic diffuse models
  and GC excess components, which affects the bubble spectra.
- **Counts-based flux.** The fit spectra are model counts divided by exposure,
  so they include PSF spill-over across the region boundary and energy
  dispersion. Both are small above about 1 GeV.
- **Region shape.** The fit region is a 40° × 40° map in the AIT projection,
  which is close to, but not exactly, the rectangle |l|, |b| < 20° used for the
  paper.
- **Uncertainties.** The published error bars are statistical only. The
  comparison does not include the uncertainty of the hole filling or of the
  fit.

## Files in this folder

| File | Purpose |
|---|---|
| `compare_bubbles.py` | The comparison script |
| `reference/` | Published spectra, template, mask and exposure |
| `bubbles_low_flux.dat`, `bubbles_high_flux.dat` | Example flux files |
| `comparison.*`, `summary.csv`, `converted_reference.csv` | Outputs of the last run |
| `Bubbles_Model/` | Copy of the fit's bubble templates (not used by the script) |
