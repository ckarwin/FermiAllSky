# Loop I spectrum vs. published Fermi-LAT results

`compare_loopI.py` compares the Loop I spectrum of the inner-Galaxy fit
([`Main/`](../../Main/README.md)) with the Fermi-LAT Loop I points shown in
Figure 3 (right) of Ackermann, Malyshev & Malyshev (2025), Space Sci. Rev.
221:92 ([doi:10.1007/s11214-025-01215-0](https://doi.org/10.1007/s11214-025-01215-0)).
Those points come from Ackermann et al. (2014), ApJ 793, 64.

The comparison is made in the paper's quantity: the Loop I **intensity**
(E² dN/dE per steradian) at |b| > 10°, averaged over the area of Loop I
(about 3 sr).

## Quick start

```bash
cd Plots/compare_loopI
python compare_loopI.py
```

Requirements: NumPy, Matplotlib and PyYAML. The script uses the Loop I shell
geometry from [`LoopI/shell_geometry.py`](../../LoopI/shell_geometry.py); it
does not need the large Loop I FITS files.

## Inputs

**Your fit.** The flux files of the two Loop I shells made by
[`../plot_flux.py`](../README.md) (`loopI_SI_flux.dat` for S1 and
`loopI_S2_flux.dat` for S2). Set their paths at the top of the script:

```python
S1_FLUX_FILE = "loopI_SI_flux.dat"
S2_FLUX_FILE = "loopI_S2_flux.dat"
FIT_LABEL = "This fit (S1 + S2)"
```

Relative paths are taken from this folder; absolute paths also work. The two
files bundled here are an example.

**Published data.** `reference/fig3_loopI_points.yaml` holds the 17 Loop I
points (112 MeV–27 GeV) with their error bars, in MeV cm⁻² s⁻¹ sr⁻¹. The paper
gives no table, so they were read from the figure's vector graphics with
`reference/extract_fig3_points.py`, which also checks its axis calibration
against the purple leptonic curve in `LoopI/leptonic_curve.txt`. The paper
itself is not included; to repeat the extraction, run
`python reference/extract_fig3_points.py /path/to/paper.pdf` (needs pdfplumber).

## Outputs

Written to this folder (each run overwrites them):

| File | Contents |
|---|---|
| `comparison.png`, `comparison.pdf` | Published points, your fit converted to the paper's quantity (S1, S2 and total, with a band for the area), the paper's leptonic model, and the ratio fit/paper |
| `converted_fit.csv` | Your fit as mean intensity over the Loop I area: S1, S2, total and the band |
| `summary.csv` | Mean fit/paper ratio and the details of the conversion |

## How the two are made comparable

The two spectra describe different things:

| | Published points | Your flux files |
|---|---|---|
| Quantity | Intensity, per sr | Flux, summed over the region |
| Area | Loop I at \|b\| > 10° (about 3 sr) | 40° × 40° around the Galactic centre, mostly \|b\| < 20° |
| Components | Loop I as one component | S1 and S2 separately |

They overlap only in the band 10° < |b| < 20°, so the conversion uses the
morphology of your Loop I model, the Wolleben (2007) shells:

1. For each shell, the region flux fixes the amplitude of that shell's map:
   amplitude = flux / ∫_region map dΩ.
2. The fitted Loop I map is the sum of the two scaled shell maps. Because the
   shells have different shapes, they are converted separately and then added;
   simply adding the two region fluxes and converting them with one factor
   would be wrong.
3. The Loop I area is taken as the brightest 3 sr of this map at |b| > 10°
   (`PAPER_SOLID_ANGLE_SR`), and the fitted intensity is the mean over that
   area. The band shows 2.5–3.5 sr (`SOLID_ANGLE_RANGE_SR`), because the paper
   only states "Ω ≈ 3 sr".

`summary.csv` also gives the ratio for a second definition: the total model
flux at |b| > 10° divided by 3 sr. The two definitions agree to about 10%.

The region is the rectangle |l|, |b| ≤ 20° on a 0.5° grid. For the bundled
example, the S1 and S2 flux files agree with the model computed this way to
within 1–3% between 0.3 and 4 GeV (5% up to 11 GeV).

## What the comparison tests

- **Normalization.** In the inner-Galaxy fit, Loop I is fixed at the
  normalizations of the all-sky fit, and its spectral shape is the paper's own
  leptonic model (the purple curve in the plot). So the shape agrees by
  construction; the comparison tests the overall normalization of the
  all-sky fit with the shell morphology.
- **Note on the Loop I normalization in `Isotropic/`.** There, the shell
  template is scaled to average 1 over all longitudes at |b| > 10° (10.4 sr),
  not over Loop I's 3 sr. A normalization of 1 therefore does not reproduce
  the published intensity; see [`LoopI/MODEL_NOTES.md`](../../LoopI/MODEL_NOTES.md).

## Limitations

- **Different Loop I templates.** The published spectrum used a template
  derived from the data (the spectral-component analysis of Ackermann et al.
  2014), not the Wolleben shells. The geometric shells spread emission over a
  larger area (S1 contains the observer), so the "brightest 3 sr" is only an
  approximation of the published area.
- **Approximate area.** The paper gives the area only as about 3 sr; the
  result scales roughly inversely with it (see the band).
- **Extrapolation.** Your fit constrains Loop I mainly in the region; the
  comparison relies on the shell geometry to extend it to |b| > 10°.
- **Different analyses.** Different data sets, event selections, IRFs and
  models of the Galactic diffuse emission, the bubbles and the isotropic
  background.
- **Counts-based flux.** The flux files are model counts divided by
  exposure, which includes energy dispersion and PSF effects (small, except
  for the steep high-energy tail).
- **Digitized data.** The published points are read from the figure's vector
  graphics, which is precise to a fraction of a percent; the error bars are as
  drawn in the figure.

## Files in this folder

| File | Purpose |
|---|---|
| `compare_loopI.py` | The comparison script |
| `reference/fig3_loopI_points.yaml` | Published Loop I points |
| `reference/extract_fig3_points.py` | Re-reads the points from the paper's PDF |
| `loopI_SI_flux.dat`, `loopI_S2_flux.dat` | Example flux files (S1 and S2) |
| `comparison.*`, `converted_fit.csv`, `summary.csv` | Outputs of the last run |
