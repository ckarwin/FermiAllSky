#!/usr/bin/env python3
"""Thin link to ``Isotropic/fit_isotropic.py``.

``generate_sun_moon_templates.py`` uses ``read_config`` from here to read the
isotropic wrapper config. This file simply loads the real driver, so there is
only one implementation. Running it directly behaves exactly like
``Isotropic/fit_isotropic.py``; prefer that script for fits.
"""

import importlib.util
from pathlib import Path
import sys

_DRIVER = Path(__file__).resolve().parents[1] / "Isotropic" / "fit_isotropic.py"
sys.path.insert(0, str(_DRIVER.parent))
_spec = importlib.util.spec_from_file_location("_fermi_isotropic_driver", _DRIVER)
_driver = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _driver
_spec.loader.exec_module(_driver)
read_config = _driver.read_config

if __name__ == "__main__":
    sys.exit(_driver.main())
