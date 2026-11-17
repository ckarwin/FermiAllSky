# Loop I input model

This is a nominal foreground model, not a new measurement or an exact reproduction of the published spatial analysis.

## Spectrum

The empirical option reconstructs the **purple dashed leptonic curve**, not the red LAT points, in Figure 3 (upper-right panel, PDF page 8) of Ackermann, Malyshev & Malyshev, *Diffuse Continuum Emission and Large Extended Sources at MeV Energies*, Space Science Reviews 221:92 (2025), DOI 10.1007/s11214-025-01215-0. The figure attributes the LAT points to Ackermann et al. (2014). Do not describe this reconstructed curve as an author-supplied numerical best-fit table.

`leptonic_curve.txt` contains the PDF vector path transformed using the logarithmic axes: E=0.1–1e6 MeV and E²I=1e-6–1 MeV cm^-2 s^-1 sr^-1. The path extends below the plot frame and ends at 113250 MeV. It is interpolated monotonically in log E–log I with PCHIP. The original PDF path has finite precision; the interpolation does not recover the authors' physical electron/radiation-field calculation or model uncertainties.

At 1 GeV, I is approximately 5.867e-10 ph cm^-2 s^-1 sr^-1 MeV^-1 (E²I=5.867e-4 MeV cm^-2 s^-1 sr^-1). The measured-fit energy range in the supplied configuration is 300–100000 MeV, within the extracted path. Above 113250 MeV an anchored cutoff tail fitted to the final six vector points supplies the response margin. This extrapolation is ours, not published data. The cube grid is 10–1e6 MeV with 32 intervals per decade. The very small high-energy tail can underflow to zero in float32; there is no positive artificial floor in the cubes. Do not extend the analysis above the extracted range without reconsidering this tail.

## Intensity and morphology

I(E,l,b) = A I_purple(E) [L_S1(l,b)+L_S2(l,b)] / C.

L is the forward line-of-sight length through a uniform-volume spherical shell. The geometries (centre l,b; distance; inner/outer radii in pc) are S1=(341,3;78;62,81), S2=(332,37;95;58,82). This assumes the same volume emissivity and spectrum for both shells; their relative strengths are fixed by geometry. S1 includes the observer inside its outer sphere and consequently has some full-sky support. No magnetic-field weighting is applied.

C is the solid-angle-weighted mean of L_S1+L_S2 over **all longitudes with |b|>10 degrees**. Thus the dimensionless combined spatial template averages to one in that explicit reference domain. A=1 makes its mean intensity equal to the purple curve. The model does not divide intensity by two or convert it to a total flux using the figure's approximate 3 sr sensitivity area. The two separate spatial maps are contributions to the sum; neither is independently normalized to mean one.

**Normalization limitation:** the figure's |b|>10 label does not specify the complete spatial weighting of the original extraction. Our all-longitude reference domain is a modelling convention, not a recovered paper footprint. Therefore A=1 is a fixed nominal amplitude, not a guaranteed reproduction of the paper's absolute flux for the Wolleben morphology. A matched published template/footprint would be needed for that equivalence. `free_norm: true` offers a single shared nuisance amplitude without changing the spectral shape or relative shell geometry.

## Mask overlap

The plot uses the mask in the uploaded config: |b|>=30 degrees and excludes |wrapped l|<=45 degrees at all latitudes. The 40x40 rectangle is explicitly **assumed centred on the Galactic centre**, since the uploaded all-sky code contains no nominal ROI centre. If it is M31-centred or elsewhere, this overlay is not that ROI.

Intrinsic flux fractions (not exposure/PSF-folded count fractions):

| Map | Retained isotropic sky | GC-centred 40x40 |
|---|---:|---:|
| S1 | 26.89% | 4.90% |
| S2 | 23.52% | 14.53% |
| Combined | 25.69% | 8.34% |

The assumed GC rectangle is entirely excluded from the isotropic fit, but Loop I contributes in both regions. Angular resolution is 0.5 degrees; these are pixel-centre estimates.


## Separate-shell extension

`loopI.mode: separate` replaces the combined source with two sources sharing the same spectral shape and common spatial normalization. Independent amplitudes multiply L_S1 and L_S2; unit amplitudes recover the combined model. Their ratio is fitted rather than inferred from a published gamma-ray shell ratio. This extends the combined geometric prescription; the main LAT studies checked here did not supply individual shell spectra. See `../PUBLICATION_NOTES.md` for the literature comparison and recommended robustness tests.
