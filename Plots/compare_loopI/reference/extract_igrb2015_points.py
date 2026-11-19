"""Re-create the igrb2015_model{A,B,C}.yaml files from the published figures (not normally needed).

    python extract_igrb2015_points.py /path/to/1410.3696v1.pdf

Reads the "Loop I / Local Loop" component (cyan bands) from the left panels of
Figures 4, 5 and 6 (foreground models A, B and C) of Ackermann et al. 2015,
ApJ 799, 86 (arXiv:1410.3696v1; PDF pages 28 and 29) as vector graphics.
Each energy bin is drawn as a band whose upper and lower edges give the
1-sigma range; the value at the bin's logarithmic centre is read from the
edges, and the central value is the logarithmic midpoint of the band.

The axes are calibrated from the major tick marks of each panel, whose
decades are read from the tick labels. Requires pdfplumber and that PDF.
"""

import sys
from pathlib import Path

import numpy as np
import pdfplumber
import yaml

HERE = Path(__file__).resolve().parent
LOOP_I_COLOUR = (0.4, 0.8, 0.8)  # cyan "Loop I / Local Loop" bands
FIGURES = [  # (PDF page, panel index from the top of the page, model, figure)
    (28, 0, "A", 4),
    (29, 0, "B", 5),
    (29, 1, "C", 6),
]


def left_panels(page):
    """Plot frames of the left panels, from top to bottom: (x0, x1, top, bottom)."""
    frames = {
        (round(r["x0"], 2), round(r["x1"], 2), round(r["top"], 2), round(r["bottom"], 2))
        for r in page.rects
        if r["x0"] < page.width / 2 and 150 < r["width"] < 220 and r["height"] > 150
    }
    return sorted(frames, key=lambda f: f[2])


def major_ticks(page, frame):
    """Major tick positions on the x axis (page x) and y axis (page top) of a panel."""
    x0, x1, top, bottom = frame
    xs, ys = set(), set()
    for line in page.lines:
        if line.get("stroking_color") != (0.0, 0.0, 0.0):
            continue
        (ax, ay), (bx, by) = line["pts"][0], line["pts"][-1]
        if abs(ax - bx) < 0.01 and x0 - 0.1 <= ax <= x1 + 0.1 and abs(max(ay, by) - bottom) < 0.1:
            if abs(ay - by) > 4:
                xs.add(round(ax, 3))
        if abs(ay - by) < 0.01 and top - 0.1 <= ay <= bottom + 0.1 and abs(min(ax, bx) - x0) < 0.1:
            if abs(ax - bx) > 5:
                ys.add(round(ay, 3))
    return sorted(xs), sorted(ys)


def tick_decades(page, frame, xs, ys):
    """Labelled major ticks and their decades (log10), read from the exponent labels."""
    x0, x1, top, bottom = frame
    words = page.extract_words(extra_attrs=["size"])  # keeps superscripts separate
    exps = [w for w in words if w["text"].lstrip("-").isdigit() and len(w["text"]) <= 2
            and w["text"] != "10"]
    xt, xdec = [], []
    for x in xs:  # exponent label just below the axis, right of "10"
        cand = [w for w in exps if bottom < w["top"] < bottom + 12 and 0 < w["x0"] - x < 6]
        if cand:  # frame edges also look like ticks but carry no label
            xt.append(x)
            xdec.append(int(cand[0]["text"]))
    yt, ydec = [], []
    for y in ys:  # exponent label just left of the axis (not the rotated axis title)
        cand = [w for w in exps if x0 - 8 < w["x1"] < x0 and abs((w["top"] + w["bottom"]) / 2 - y) < 6]
        if cand:
            yt.append(y)
            ydec.append(int(cand[0]["text"]))
    return np.array(xt), np.array(xdec), np.array(yt), np.array(ydec)


def loop_i_bins(page, frame):
    """Loop I band polygons in a panel: list of (x_left, x_right, edge, other_edge)."""
    x0, x1, top, bottom = frame
    out = []
    for c in page.curves:
        if c.get("stroking_color") != LOOP_I_COLOUR or c.get("non_stroking_color") != "P02":
            continue
        if not (x0 - 0.1 <= c["x0"] and c["x1"] <= x1 + 0.1 and top <= c["top"] and c["bottom"] <= bottom):
            continue
        pts = np.array(c["pts"][:6])  # one edge left->right, the other right->left
        if len(c["pts"]) < 6 or not np.allclose(pts[[0, 5], 0], pts[0, 0]):
            continue  # legend swatch (a rectangle), not a bin
        edge_a = pts[:3]
        edge_b = pts[3:6][::-1]
        out.append((pts[0, 0], pts[2, 0], edge_a, edge_b))
    return sorted(out, key=lambda b: b[0])


def extract(pdf_path):
    pdf = pdfplumber.open(pdf_path)
    for page_no, panel, model, figure in FIGURES:
        page = pdf.pages[page_no - 1]
        frame = left_panels(page)[panel]
        xs, ys = major_ticks(page, frame)
        xs, xdec, ys, ydec = tick_decades(page, frame, xs, ys)
        assert len(xs) >= 3 and len(ys) >= 3, "too few labelled ticks"
        # Linear fits of decade vs. page coordinate; residuals check the calibration.
        px = np.polyfit(xs, xdec, 1)
        py = np.polyfit(ys, ydec, 1)
        assert np.max(np.abs(np.polyval(px, xs) - xdec)) < 1e-3, "x calibration"
        assert np.max(np.abs(np.polyval(py, ys) - ydec)) < 1e-3, "y calibration"
        rows = []
        for xl, xr, edge_a, edge_b in loop_i_bins(page, frame):
            xc = (xl + xr) / 2  # logarithmic bin centre (log x axis)
            up = np.interp(xc, edge_a[:, 0], edge_a[:, 1])
            lo = np.interp(xc, edge_b[:, 0], edge_b[:, 1])
            log_lo, log_hi = sorted((np.polyval(py, up), np.polyval(py, lo)))
            rows.append([10 ** np.polyval(px, xc), 10 ** ((log_hi + log_lo) / 2),
                         10 ** log_lo, 10 ** log_hi,
                         10 ** np.polyval(px, xl), 10 ** np.polyval(px, xr)])
        rows = np.array(rows)
        print(f"Model {model}: {len(rows)} bins, {rows[0, 4]:.4g}-{rows[-1, 5]:.4g} MeV")
        out = {
            "label": f"Ackermann et al. 2015, model {model} (Wolleben template)",
            "citation": f"Ackermann et al. 2015, ApJ 799, 86 (arXiv:1410.3696), Fig. {figure} (left), "
            "'Loop I / Local Loop'",
            "template": "Geometric Loop I / Local Loop template (Wolleben 2007)",
            "region": {
                "type": "latitude",
                "b_min_deg": 20.0,
                "note": "Average intensity over all longitudes at |b| > 20 deg.",
            },
            # Loop I is fitted bin by bin only below 12.8 GeV; above, the paper fixes the
            # normalizations (the drawn bands are identical for models A, B and C there).
            "ratio_band_mev": [100.0, 12800.0],
            "source": "Read from the PDF vector graphics of the cyan bands (value at the logarithmic "
            "bin centre; central value = logarithmic midpoint of the band); not an author-supplied table.",
            "energy_unit": "MeV",
            "data_unit": "MeV cm^-2 s^-1 sr^-1 (E^2 dN/dE per steradian)",
            "energy": rows[:, 0].round(3).tolist(),
            "e2_intensity": rows[:, 1].tolist(),
            "e2_intensity_low": rows[:, 2].tolist(),
            "e2_intensity_high": rows[:, 3].tolist(),
            "energy_bin_low": rows[:, 4].round(3).tolist(),
            "energy_bin_high": rows[:, 5].round(3).tolist(),
        }
        path = HERE / f"igrb2015_model{model}.yaml"
        path.write_text(yaml.safe_dump(out, sort_keys=False, width=100))
        print(f"Wrote {path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    extract(sys.argv[1])
