"""PFS spectral simulator: the modern snake_case API plus its underlying
legacy engine.

`SimSpecParams`/`load_params`/`run_sim_spec` (re-exported here) are the
modern entry points (see `pfsspecsim.sim.params`/`pfsspecsim.sim.engine`);
`pfsspecsim.sim.pfsspec.Pfsspec` is the underlying simulator engine
`run_sim_spec` drives, and `pfsspecsim.sim.dm_utils` builds the
`pfs.datamodel` objects it writes.
"""

from __future__ import annotations

from .engine import run_sim_spec
from .params import SimSpecParams, load_params

__all__ = [
    "SimSpecParams",
    "load_params",
    "run_sim_spec",
]
