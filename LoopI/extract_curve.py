"""Re-create leptonic_curve.txt from the published figure (not normally needed).

    python extract_curve.py /path/to/paper.pdf

Reads the purple dashed curve of Figure 3 (upper right, PDF page 8) of
Ackermann, Malyshev & Malyshev, Space Sci. Rev. 221:92 (2025) as a vector path
and converts it to physical units with the figure's logarithmic axes. Requires
PyMuPDF (``fitz``) and that exact PDF; the axis calibration is hard-coded.
"""

import sys
from pathlib import Path
import fitz
import numpy as np

p = fitz.open(sys.argv[1])[7]
paths = [
    d
    for d in p.get_drawings()
    if d["color"]
    and abs(d["color"][0] - 0.5804) < 0.002
    and len(d["items"]) == 73
    and d["rect"].y0 < 150
]
if len(paths) != 1:
    raise ValueError("PDF does not match the supplied review; recalibrate axes")
d = paths[0]
pts = np.array([tuple(d["items"][0][1])] + [tuple(i[-1]) for i in d["items"]])
e = 10 ** (-1 + 7 * (pts[:, 0] - 246.95269775) / (389.09613037 - 246.95269775))
sed = 10 ** (-6 * (pts[:, 1] - 60.3150711) / (158.09176636 - 60.3150711))
np.savetxt(
    Path(__file__).parent / "leptonic_curve.txt",
    np.c_[e, sed / e**2, sed],
    header="Energy_MeV differential_intensity_ph_cm-2_s-1_sr-1_MeV-1 E2I_MeV_cm-2_s-1_sr-1; extracted PDF vector path, not author numerical table",
)
