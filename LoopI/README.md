# Loop I templates

Loop I is a large, nearby structure seen at high Galactic latitude in radio and
gamma rays. This folder builds the Loop I foreground templates used by the
isotropic fit ([`Isotropic/`](../Isotropic/)) and by the inner-Galaxy analysis
([`Main/`](../Main/)).

The model uses the two shells S1 and S2 of Wolleben (2007):

```
I(E, l, b) = A × I_curve(E) × [L_S1(l, b) + L_S2(l, b)] / C
```

- `L_S1`, `L_S2`: line-of-sight length through each shell, assuming uniform
  emissivity (geometry in [`shell_geometry.py`](shell_geometry.py)).
- `C`: makes the summed map average to 1 over all longitudes at |b| > 10°.
- `I_curve(E)`: an empirical leptonic spectrum reconstructed from a published
  figure ([`leptonic_curve.txt`](leptonic_curve.txt)); 5.87 × 10⁻¹⁰
  ph cm⁻² s⁻¹ sr⁻¹ MeV⁻¹ at 1 GeV.
- `A`: the amplitude fitted in the isotropic analysis (`loopI.normalization`).

The isotropic fit can use the shells together (`mode: combined`) or with
independent amplitudes (`mode: separate`), and can replace the empirical
spectrum by a fitted power law, cutoff power law or log parabola. See the
[Loop I section of the isotropic guide](../Isotropic/README.md#loop-i). The
scientific assumptions and limitations of the model are described in
[`MODEL_NOTES.md`](MODEL_NOTES.md).

## Building the templates

The large FITS files are **not stored in git**. Build them once on a fresh
checkout (needs NumPy, SciPy, Astropy, Matplotlib and PyYAML):

```bash
python LoopI/build_templates.py --config Isotropic/config.yaml
```

The files are written to `loopI.directory` from the config if it is set, or to
this folder otherwise. Use `--outdir DIR` to choose another folder, and point
`loopI.directory` in the isotropic config (and `loopI_directory` in
`Main/run_3ML_master.py`) to the same place.

The build is deterministic, so rebuilding gives the same templates. It does,
however, update the files' modification times, and the isotropic fit then asks
for a new preparation (in a new `outdir`, with `reuse_from` the old one). Do
not rebuild while a fit that uses these files is running.

## Files

Tracked in git:

| File | Contents |
|---|---|
| `build_templates.py` | Builds all template files below |
| `shell_geometry.py` | Shell parameters and line-of-sight path lengths |
| `leptonic_curve.txt` | Reconstructed spectrum: energy [MeV], intensity, E²×intensity |
| `extract_curve.py` | Re-creates `leptonic_curve.txt` from the paper's PDF (needs PyMuPDF and the PDF; not normally needed) |
| `empirical_intensity.txt`, `model_metadata.json` | Spectrum on the cube energy grid and model metadata (rewritten by the build) |
| `loopI_spectrum.png`, `loopI_overlap.png` | Diagnostic plots: spectrum, and shells with the isotropic mask and the 40° × 40° inner-Galaxy region (rewritten by the build) |
| `MODEL_NOTES.md` | Model assumptions and limitations |

Generated (not in git):

| File | Made by | Contents |
|---|---|---|
| `loopI_combined_spatial.fits`, `loopI_S1_spatial.fits`, `loopI_S2_spatial.fits` | `build_templates.py` | Dimensionless 0.5° Galactic CAR maps; S1 + S2 = combined |
| `loopI_empirical_mapcube.fits` | `build_templates.py` | Combined map × empirical spectrum (10 MeV–1 TeV) |
| `loopI_powerlaw_seed_mapcube.fits` | `build_templates.py` | Combined map × E^−2.4 seed spectrum (for analytic spectral fits) |
| `loopI_S1_empirical_mapcube.fits`, `loopI_S2_empirical_mapcube.fits` (and `*_powerlaw_seed_*`) | isotropic `prepare` in separate mode | Per-shell cubes; also used by `Main/run_3ML_master.py` |
