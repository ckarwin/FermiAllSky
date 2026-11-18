# Updated: ROI-integrated Fermi bubble comparison

Run `python compare_2017.py`. Dependencies: numpy pandas matplotlib astropy healpy pyyaml. Inputs are bundled. Outputs: comparison.png/pdf, converted_reference.csv, summary.csv.

## Switching fits

Edit just these settings near the top of `compare_2017.py`:

```python
LOW_LAT_FLUX_FILE = "bubbles_low_flux.dat"
HIGH_LAT_FLUX_FILE = "bubbles_high_flux.dat"
FIT_LABEL = "Your fit"
```

Use filenames relative to the script directory, or absolute paths. Each file should
have the same two-line header and columns as the bundled exporter `.dat` files.
Run `python compare_2017.py`; this overwrites the comparison plots and CSV files
beside the script. Change `FIT_LABEL` if you want a different legend label.
Plot limits are set in `main()` and can be adjusted for other fits.

The default integration band is 1–80 GeV. If another fit does not cover 10 GeV
or the whole integration band, the corresponding summary entries are blank
(printed as NaN), rather than extrapolated. The numerical results below describe
the bundled fit. Refactoring was checked against both original CSV outputs;
all numerical values are unchanged.

## Verified export convention

The supplied plot_flux.py calls Plots.calculate_flux(model_cube, exposure_cube, ..., filetype='fermipy_mcube'). Inspection of ckarwin/Fermi-LAT_Binned_Analysis, lat_binned_analysis/analysis.py (blob 9177020d573a72e3c76f4c42c075a873137638c1) shows that this method takes the spatial dimensions from the input counts cube and sums over every x,y pixel:

    E_center^2 / delta_E * sum_pixels(counts / sqrt(exposure_low*exposure_high))

Thus the .dat spectra are already integrated over the input cube footprint, i.e. the ROI. They are not full-template spectral fluxes. No factors 0.92956 or 0.32790 should be applied to them. The alternative full-template curves and ambiguous export options have been removed. The user's flux-generation and plotting script needs no spatial-normalization change; the comparison script has been corrected.

The exporter uses PSF-convolved model counts and bin-edge exposure estimates. Its flux is a counts/exposure estimate over the reconstructed ROI, not a PSF-deconvolved intrinsic template integral; energy dispersion can also affect it. The comparison to intrinsic model flux is therefore approximate, particularly at low energies and near ROI boundaries. The exact cube footprint was not supplied; this comparison retains the stated Galactic rectangle |l|<20, |b|<20 degrees. If the cube uses a projection with different boundaries, that geometrical approximation remains.

A separate exporter issue: its 'Total energy flux' header simply sums the E^2 dF/dE samples. This is NOT an energy-integrated flux. Do not use that header for a band-flux comparison. This package computes band flux via integral (E^2 dF/dE) dlnE instead. Spectral columns are unaffected by the misleading header.

## Reference conversion

Ackermann et al. 2017, https://arxiv.org/abs/1704.03910, Figure 9a. Official data: https://www-glast.stanford.edu/pub_data/1220/GCexcess_Pass8_1704.03910_data/ . The spectral table specifies all-sky flux except masked cores of 200 bright 3FGL sources.

Split the released HEALPix RING nside=128 bubble morphology at |b|=10 degrees. Approximate its intensity by dividing its counts morphology by the supplied exposure at each energy. Multiply published spectra and error bars by:

    sum_ROI(template/exposure) / sum_allsky(template/exposure * PS_mask)

The numerator covers the full rectangle without source holes, predicting bubble emission in masked locations from the model. At 10 GeV these conversion factors are 0.91854 low and 0.27804 high. Using exposure around 1 GeV instead changes them by less than 0.1%. This is a model-based aperture conversion, not a PSF deconvolution or reproduction of the original fit. Native pixel-center integration is used. User templates have modified morphology and a smooth 9--11 degree latitude partition, versus the paper's sharp split.

## Updated results

At 10 GeV, spatially integrated E^2 dF/dE [GeV cm^-2 s^-1]:
                  Paper in ROI        User ROI
Low latitude      1.697e-7             1.376e-7
High latitude     1.318e-7             5.132e-8

Energy-integrated flux over 1--80 GeV [GeV cm^-2 s^-1]:
Low latitude      7.415e-7             6.096e-7
High latitude     4.860e-7             2.512e-7

The user's low component is about 19% below the reference at 10 GeV, and high about 2.57 times lower. The high component has a harder shape. These are amplitude comparisons, not statistical significances. Published errors exclude aperture-conversion and foreground-model systematic uncertainty; no user fit uncertainties were supplied.

Band integrals use log-log interpolation within the supplied spectra; the 1--80 GeV range avoids extrapolation. No band-flux errors are inferred without covariance information.

Validation: script executed successfully with actual input files; regenerated plots inspected. Updated summary contains only the verified ROI convention.
