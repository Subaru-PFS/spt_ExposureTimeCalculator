"""Pure-Python PFS exposure time calculator engine.

Ports the C engine `src/gsetc.c` to numpy/scipy/astropy. See the project
plan (`PLAN-pure-python-etc.md`) for the module breakdown and task order.
`run_etc` / `EtcResults` land with `engine.py` in a later task; for now
this re-exports what already exists.
"""

from __future__ import annotations

from .params import EtcParams, MagSpec, calc_obscuration, load_params, resolve_degrade

__all__ = [
    "EtcParams",
    "MagSpec",
    "calc_obscuration",
    "load_params",
    "resolve_degrade",
]
