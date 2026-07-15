# -*- coding: utf-8 -*-

# Version information
try:
    from ._version import __version__
except ImportError:
    # Package not installed or built, use placeholder
    __version__ = "unknown"

import sys as _sys

from .legacy import pfsetc
from .sim import dm_utils, pfsspec

# Keep the pre-v2.0 public import paths working (deprecated, may be removed
# in a future release): `import pfsspecsim.pfsspec` /
# `from pfsspecsim.pfsspec import Pfsspec` etc. resolve to the real modules
# below (now living in `pfsspecsim.legacy`/`pfsspecsim.sim`).
_sys.modules[__name__ + ".pfsetc"] = pfsetc
_sys.modules[__name__ + ".pfsspec"] = pfsspec
_sys.modules[__name__ + ".dm_utils"] = dm_utils
