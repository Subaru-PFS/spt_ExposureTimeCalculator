"""Pure-Python PFS exposure time calculator engine.

Ports the C engine (vendored, frozen, at `legacy/c_src/gsetc.c`) to
numpy/scipy/astropy.
"""

from __future__ import annotations

from .engine import EtcResults, run_etc, run_etc_files
from .params import EtcParams, MagSpec, calc_obscuration, load_params, resolve_degrade

__all__ = [
    "EtcParams",
    "EtcResults",
    "MagSpec",
    "calc_obscuration",
    "load_params",
    "resolve_degrade",
    "run_etc",
    "run_etc_files",
]
