#!/usr/bin/env python3
"""Compatibility entry point; the maintained pipeline lives in Isotropic/.

The Sun/Moon generator imports read_config here. Delegate both that import and
legacy CLI commands to the same driver so fixes and cache rules cannot diverge.
"""
import importlib.util
from pathlib import Path
import sys

_DRIVER = Path(__file__).resolve().parents[1] / 'Isotropic' / 'fit_isotropic.py'
sys.path.insert(0, str(_DRIVER.parent))
_spec = importlib.util.spec_from_file_location('_fermi_isotropic_driver', _DRIVER)
_driver = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _driver
_spec.loader.exec_module(_driver)
read_config = _driver.read_config

if __name__ == '__main__':
    sys.exit(_driver.main())
