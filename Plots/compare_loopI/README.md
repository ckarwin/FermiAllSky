# Loop I spectrum vs. published results

`compare_loopI.py` compares the Loop I spectrum of the all-sky fit
([`Isotropic/`](../../Isotropic/README.md)) with published Fermi-LAT Loop I
spectra. Each published spectrum refers to its own sky region, so the fitted
Loop I sky is averaged over that region before comparing. Showing several
published results together gives an idea of the systematic uncertainty.

## Quick start

```bash
cd Plots/compare_loopI
python compare_loopI.py
```

Requirements: NumPy, Matplotlib and PyYAML. The script rebuilds the Loop I
shell templates from [`LoopI/shell_geometry.py`](../../LoopI/shell_geometry.py);
it does not need the large Loop I FITS files.

After each inner-Galaxy fit, [`make_all_plots.py`](../README.md) runs this
comparison for the all-sky fit that the inner-Galaxy fit used
(`loopI_comparison.png` in the fit's plot folder), and
[`compare_fits.py`](../README.md) overlays several fits.

## Inputs

**Your fit.** The results folder of an all-sky fit
(`<outdir>/bestfit_<hash>` or `<outdir>/fit_<hash>`). The fit writes the Loop I
spectra there: `loopI_s1_spectrum.txt` and `loopI_s2_spectrum.txt` (Loop I
`mode: separate`) or `loopI_spectrum.txt` (`mode: combined`). Set the folder at
the top of the script:

```python
ALLSKY_RESULTS_DIR = "/path/to/Isotropic/output_.../bestfit_xxxxxxxxxx"
FIT_LABEL = "This work"
```

Relative paths are taken from this folder. `example_allsky/` is an example: the
empirical spectrum scaled by the S1 and S2 normalizations that `Main/` uses
(0.0286 and 0.2149).

The all-sky fit, not the inner-Galaxy fit, is the right input: `Main/` keeps
Loop I fixed at the all-sky values, so its Loop I spectra carry no new
information.

**Published results.** Every YAML file in `reference/` is plotted (set
`REFERENCE_FILES` to choose). See [Published results](#published-results) for
the current set and [Adding a published result](#adding-a-published-result) for
the format.

## Outputs

Written to this folder (each run overwrites them):

| File | Contents |
|---|---|
| `comparison.png`, `comparison.pdf` | Published points (one colour each); the fit averaged over each published region (black, one line style per region); the ratio fit/published |
| `converted_fit.csv` | E² I of the fit averaged over each published region |
| `summary.csv` | For each published result: region, its solid angle, S1 share, mean fit/published ratio and the energies used for it |

## How the comparison works

The all-sky fit writes, for each Loop I term, a spectral coefficient c(E) that
multiplies the dimensionless shell template. The fitted Loop I intensity is

```
I(E, l, b) = c_S1(E) × S1(l, b) + c_S2(E) × S2(l, b)
```

with S1 and S2 the Wolleben (2007) shell templates exactly as
`LoopI/build_templates.py` makes them (S1 + S2 averages to 1 over all
longitudes at |b| > 10°). In combined mode there is one term, c(E) × (S1 + S2).

The script averages this map over the region of each published result:

| Region type | Region | Used for |
|---|---|---|
| `latitude` | All longitudes at \|b\| > `b_min_deg` | Results quoted as an average over a latitude cut |
| `flux_per_area` | Total flux at \|b\| > `b_min_deg`, divided by `solid_angle_sr` | Results quoted as a total flux spread over a nominal Loop I area |
| `footprint` | The brightest `solid_angle_sr` of the fitted map at \|b\| > `b_min_deg` | Results quoted as an average over a published Loop I template area |

S1 and S2 have different shapes, so they are averaged separately and then
added.

## Published results

| File | Result | Template | Region |
|---|---|---|---|
| `reference/ackermann2014.yaml` | Ackermann et al. 2014, ApJ 793, 64, Fig. 15 (as replotted in Ackermann, Malyshev & Malyshev 2025, Space Sci. Rev. 221:92, Fig. 3) | Derived from the data (spectral component analysis) | Total flux at \|b\| > 10° divided by 3 sr |
| `reference/igrb2015_model{A,B,C}.yaml` | Ackermann et al. 2015, ApJ 799, 86 (arXiv:1410.3696), Figs. 4–6, "Loop I / Local Loop", foreground models A, B and C | Wolleben (2007) shells, as here | Mean over all longitudes at \|b\| > 20° |

None of these papers tabulates its Loop I spectrum, so the points were read
from the vector graphics of the figures:

- `reference/extract_fig3_points.py /path/to/2025_paper.pdf` re-creates
  `ackermann2014.yaml` and checks its axis calibration against
  `LoopI/leptonic_curve.txt`.
- `reference/extract_igrb2015_points.py /path/to/1410.3696v1.pdf` re-creates
  the three `igrb2015_model*.yaml` files. Each energy bin is a band in the
  figure; the value is read at the bin's logarithmic centre, with the band as
  the error bar. The axes are calibrated from the labelled ticks of each panel.

Both need pdfplumber and the exact PDF.

**The 2015 results are the most direct comparison**: same template and a
well-defined region. Their three foreground models show the systematic spread
from the Galactic diffuse model alone. In that paper Loop I is fitted bin by
bin only below 12.8 GeV; above, its normalization is fixed (the bands are the
same for all three models). Those points are drawn hollow and left out of the
mean ratio (`ratio_band_mev` in the YAML files).

**What the 2014 points mean.** Fig. 15 of the 2014 paper shows the spectra of
all components averaged over the whole analysis region, |b| > 10° (10.4 sr).
The 2025 paper replots the Loop I spectrum as an intensity over Loop I's own
area, Ω ≈ 3 sr: its points are the 2014 values × 10.4 / 3. Reading Fig. 15
directly confirms this (the 2025 points are 3.5 × higher, which corresponds to
2.96 sr). So the 2014 points are the total Loop I flux at |b| > 10°
divided by 3 sr, and the fit is evaluated in the same way. This uses no
assumption about the shape of the 2014 template, and does not depend on how
the fit splits Loop I between S1 and S2.

The 2025 error bars (about ±18% at low energy, growing above a few GeV) are
larger than the statistical errors of the 2014 baseline fit in Fig. 15. They
most likely include systematic uncertainties, which the 2014 paper estimated
by varying the foreground models and templates but did not publish for
Loop I. This is why the 2025 points are used here.

**The 2014 and 2015 results differ by about 30–50%** (the fit-to-published
ratios differ by that much), a measure of the
systematic uncertainty from how Loop I is modelled (template derived from
the data vs. Wolleben shells, different foreground models).

Other Loop I results checked:

| Result | Why not included |
|---|---|
| Jóhannesson & Porter 2021 (arXiv:2104.13708), Table 4 | Mean of the whole smooth emission minus isotropic (includes large-scale inverse Compton), over a region given only as a figure outline; not Loop I alone |
| Ajello et al. 2016; Ackermann et al. 2017; DAMPE 2025; Totani 2025; Cholis et al. 2026 | Use the same Wolleben shells, but publish no Loop I spectrum or normalization |
| Casandjian & Grenier 2009; Acero et al. 2016 | Radio (408 MHz) templates and patches; no Loop I spectrum given |

## Adding a published result

Add one YAML file to `reference/`:

```yaml
label: Author et al. 20XX (template)        # legend text
citation: Journal reference, figure or table
template: Which Loop I template the paper used
region:
  type: latitude                             # or: flux_per_area, footprint
  b_min_deg: 20.0
  # flux_per_area and footprint:
  # solid_angle_sr: 3.0
  # footprint only:
  # solid_angle_range_sr: [2.5, 3.5]         # optional band
energy: [...]                                # MeV
e2_intensity: [...]                          # E^2 dN/dE per sr, MeV cm^-2 s^-1 sr^-1
e2_intensity_low: [...]                      # 1-sigma lower end
e2_intensity_high: [...]                     # 1-sigma upper end
```

Use one file per data set (for example one per foreground model).

## What the comparison tests

- **Mainly the normalization.** The fit uses the leptonic spectral shape of
  Ackermann et al. (2025), which was fitted to the 2014 points, so the shape
  agrees with those by construction. The 2015 points (a different template and
  analysis) are a more independent check of the shape.
- **Shell weights matter for the 2015 comparison.** The \|b\| > 20° average
  depends on the relative S1 and S2 amplitudes, which are free in this fit; the
  2015 paper used one template.
- **Consistency, not a measurement.** Loop I is a foreground that trades off
  against the inverse Compton, gas and isotropic components, and the all-sky
  mask removes |l| < 45° at all latitudes (including much of the bright North
  Polar Spur). Different templates with different background models can give
  Loop I amplitudes a factor of a few apart.
- **Different analyses.** Data sets, event selections, IRFs and models of the
  Galactic diffuse emission, the bubbles and the isotropic background all
  differ.
- **Statistical errors of the fit are not shown.** They are small compared with
  the differences between analyses.

## Files in this folder

| File | Purpose |
|---|---|
| `compare_loopI.py` | The comparison script |
| `reference/*.yaml` | Published Loop I spectra |
| `reference/extract_fig3_points.py` | Re-reads the 2014 points from the 2025 paper's PDF |
| `reference/extract_igrb2015_points.py` | Re-reads the 2015 points from that paper's PDF |
| `example_allsky/` | Example all-sky Loop I spectra (S1 and S2) |
| `comparison.*`, `converted_fit.csv`, `summary.csv` | Outputs of the last run |
