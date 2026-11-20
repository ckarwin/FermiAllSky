"""Re-create ackermann2014.yaml from the published figure (not normally needed).

    python extract_fig3_points.py /path/to/paper.pdf

Reads the red Loop I data points and their error bars from Figure 3 (upper
right panel, PDF page 8) of Ackermann, Malyshev & Malyshev, Space Sci. Rev.
221:92 (2025) as vector graphics, and converts them to physical units with the
figure's logarithmic axes (E = 0.1-1e6 MeV, E^2 I = 1e-6-1 MeV cm^-2 s^-1 sr^-1).
The axis calibration is the one used by LoopI/extract_curve.py; it is checked
here by recovering the purple leptonic curve and comparing it with
LoopI/leptonic_curve.txt. Requires pdfplumber and that exact PDF.
"""

import sys
from pathlib import Path

import numpy as np
import pdfplumber
import yaml

HERE = Path(__file__).resolve().parent
# Axis calibration of the upper-right panel (page coordinates in points, y down).
X0, X1 = 246.95269775, 389.09613037  # E = 1e-1 and 1e6 MeV
Y0, Y1 = 60.3150711, 158.09176636  # E^2 I = 1e0 and 1e-6 MeV cm^-2 s^-1 sr^-1


def energy(x):
    """Page x coordinate -> energy [MeV]."""
    return 10 ** (-1 + 7 * (x - X0) / (X1 - X0))


def e2i(y):
    """Page y coordinate -> E^2 I [MeV cm^-2 s^-1 sr^-1]."""
    return 10 ** (-6 * (y - Y0) / (Y1 - Y0))


def is_color(c, ref):
    return isinstance(c, (tuple, list)) and len(c) == 3 and all(abs(a - b) < 0.003 for a, b in zip(c, ref))


def main(pdf_path):
    page = pdfplumber.open(pdf_path).pages[7]

    # Calibration check: the purple leptonic curve must match LoopI/leptonic_curve.txt.
    curve = [o for o in page.curves if len(o["pts"]) == 74 and o["top"] < 150
             and is_color(o.get("stroking_color"), (0.5804, 0.40392, 0.74118))]
    reference = np.loadtxt(HERE.parents[2] / "LoopI" / "leptonic_curve.txt")
    if len(curve) != 1:
        raise ValueError("PDF does not match the expected figure; recalibrate the axes")
    pts = np.array(curve[0]["pts"], float)
    if np.max(np.abs(np.log10(e2i(pts[:, 1]) / reference[:, 2]))) > 1e-4:
        raise ValueError("Axis calibration does not reproduce LoopI/leptonic_curve.txt")

    # Red filled markers and red error bars inside the panel; the first marker
    # (top of the panel) is the legend symbol and is skipped.
    def in_panel(o):
        return 240 < o["x0"] and o["x1"] < 395 and 55 < o["top"] and o["bottom"] < 165

    red = (1.0, 0.0, 0.0)
    marks = sorted((o for o in page.curves if in_panel(o) and o.get("fill")
                    and is_color(o.get("non_stroking_color"), red)), key=lambda o: o["x0"])[1:]
    bars = sorted((o for o in page.lines if in_panel(o)
                   and is_color(o.get("stroking_color"), red)), key=lambda o: o["x0"])[1:]
    if len(marks) != 17 or len(bars) != 17:
        raise ValueError(f"Expected 17 points, found {len(marks)} markers and {len(bars)} bars")

    rows = []
    for m, b in zip(marks, bars):
        x = (m["x0"] + m["x1"]) / 2
        y = (m["top"] + m["bottom"]) / 2
        if abs(x - b["x0"]) > 0.01:
            raise ValueError("Marker and error bar do not line up")
        rows.append([float(energy(x)), float(e2i(y)), float(e2i(b["bottom"])), float(e2i(b["top"]))])
    rows = np.array(rows)

    out = {
        "label": "Ackermann et al. 2014 (SCA template)",
        "citation": "Ackermann et al. 2014, ApJ 793, 64 (Fig. 15); as replotted in Ackermann, "
        "Malyshev & Malyshev 2025, Space Sci. Rev. 221:92, Fig. 3 (right)",
        "template": "Data-derived Loop I template (spectral component analysis, 1 sigma threshold)",
        "region": {
            "type": "flux_per_area",
            "b_min_deg": 10.0,
            "solid_angle_sr": 3.0,
            "note": "Ackermann et al. 2014, Fig. 15 gives the Loop I (SCA) intensity averaged over all "
            "of |b| > 10 deg (10.4 sr). Ackermann et al. 2025 rescaled it to Loop I's area of about "
            "3 sr: the points are the total |b| > 10 deg flux divided by 3 sr (checked against Fig. 15: "
            "ratio about 3.5, i.e. 2.96 sr).",
        },
        "source": "Figure 3 (right), red circles, of Ackermann, Malyshev & Malyshev, "
        "Space Sci. Rev. 221:92 (2025), doi:10.1007/s11214-025-01215-0; data from "
        "Ackermann et al. 2014, ApJ 793, 64",
        "comment": "Loop I intensity at |b| > 10 deg (Fermi-LAT). Read from the PDF vector "
        "graphics (marker centres and error-bar ends); not an author-supplied table.",
        "energy_unit": "MeV",
        "data_unit": "MeV cm^-2 s^-1 sr^-1 (E^2 dN/dE per steradian)",
        "energy": rows[:, 0].round(4).tolist(),
        "e2_intensity": rows[:, 1].tolist(),
        "e2_intensity_low": rows[:, 2].tolist(),
        "e2_intensity_high": rows[:, 3].tolist(),
    }
    (HERE / "ackermann2014.yaml").write_text(yaml.safe_dump(out, sort_keys=False))
    print(f"Wrote {len(rows)} points to {HERE / 'ackermann2014.yaml'}")


if __name__ == "__main__":
    main(sys.argv[1])
