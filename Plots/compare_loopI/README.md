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
| `comparison.png`, `comparison.pdf` | Published points; the fit averaged over each published region (same colour); the ratio fit/published |
| `converted_fit.csv` | E² I of the fit averaged over each published region |
| `summary.csv` | For each published result: region, its solid angle, S1 share, mean fit/published ratio (0.3–30 GeV) |

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
| `footprint` | The brightest `solid_angle_sr` of the fitted map at \|b\| > `b_min_deg` | Results quoted as an average over a published Loop I template area |

S1 and S2 have different shapes, so they are averaged separately and then
added.

## Published results

| File | Result | Template | Region | Status |
|---|---|---|---|---|
| `reference/ackermann2014.yaml` | Ackermann et al. 2014, ApJ 793, 64, Fig. 15 (as replotted in Ackermann, Malyshev & Malyshev 2025, Space Sci. Rev. 221:92, Fig. 3) | Data-derived (spectral component analysis) | Mean over the Loop I template at \|b\| > 10°; the area, about 3 sr, is **assumed** | Included |

The 2014 points were read from the vector graphics of the 2025 figure with
`reference/extract_fig3_points.py` (which checks its axis calibration against
`LoopI/leptonic_curve.txt`). To repeat it: `python reference/extract_fig3_points.py /path/to/paper.pdf`
(needs pdfplumber).

**The 3 sr area is an assumption.** Neither paper states over which area the
2014 Loop I intensity is averaged. The 2014 paper gives the bubble spectrum as
the mean over the bubble template, and the 2025 paper gives Loop I an area of
about 3 sr, so the Loop I points are taken as the mean over a ~3 sr template.
The 2014 template is derived from the data, not from the shells, so "the
brightest 3 sr of the shell model" only stands in for it. The shaded band shows
2.5–3.5 sr. Averaging over all of |b| > 10° (10.4 sr) instead would lower the
fit by about a factor of 3.

Other Loop I results checked:

| Result | Why not included |
|---|---|
| Ackermann et al. 2015, ApJ 799, 86 (IGRB), Figs. 4–6 | Same Wolleben template, averaged over \|b\| > 20°, for foreground models A, B, C: a direct comparison. Only shown in figures; to be added once digitized from the PDF |
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
  type: latitude                             # or: footprint
  b_min_deg: 20.0
  # footprint only:
  # solid_angle_sr: 3.0
  # solid_angle_range_sr: [2.5, 3.5]         # optional band
energy: [...]                                # MeV
e2_intensity: [...]                          # E^2 dN/dE per sr, MeV cm^-2 s^-1 sr^-1
e2_intensity_low: [...]                      # 1-sigma lower end
e2_intensity_high: [...]                     # 1-sigma upper end
```

Use one file per data set (for example one per foreground model).

## What the comparison tests

- **Normalization only.** The fit uses the leptonic spectral shape of
  Ackermann et al. (2025), which was fitted to the 2014 points, so the shapes
  agree by construction.
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
| `example_allsky/` | Example all-sky Loop I spectra (S1 and S2) |
| `comparison.*`, `converted_fit.csv`, `summary.csv` | Outputs of the last run |
