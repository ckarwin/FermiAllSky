"""Make the analysis folders importable by the tests (Isotropic/ takes priority)."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
for directory in ("Sun_and_Moon", "LoopI", "DM_Model", "Isotropic"):
    sys.path.insert(0, str(ROOT / directory))
