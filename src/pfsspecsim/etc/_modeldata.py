"""Loader for the extracted C lookup tables (``etc/data/modeldata.npz``).

The archive was produced by `tools/extract_modeldata.py` (T1) from
`src/gsetc.c` / `src/modeldata.h`; see `etc/data/README.md` for provenance,
shapes, units, and upstream references for each array.
"""

from __future__ import annotations

import functools
from importlib import resources
from typing import NamedTuple

import numpy as np


class ModelData(NamedTuple):
    """Container for the extracted lookup tables.

    Attributes
    ----------
    uves_lambda, uves_int : ndarray, shape (2816,)
        UVES sky emission-line atlas: air wavelength [nm], intensity
        [1e-12 erg/m2/s/arcsec2].
    oh_data : ndarray, shape (698, 2)
        OH airglow line list: [vacuum wavelength nm, intensity].
    atm_trans_kp : ndarray, shape (40001,)
        Kitt Peak atmospheric transmission table, 0.025 nm grid from 500nm.
    mk_trans_3mm : ndarray, shape (30001,)
        Mauna Kea (3mm precipitable water) transmission table, 0.02 nm grid
        from 900nm.
    dust_norm : ndarray, shape (201,)
        Galactic dust extinction normalization curve (Draine).
    si_index_table : ndarray, shape (181,)
        Silicon refractive index table, 5nm grid from 200nm.
    """

    uves_lambda: np.ndarray
    uves_int: np.ndarray
    oh_data: np.ndarray
    atm_trans_kp: np.ndarray
    mk_trans_3mm: np.ndarray
    dust_norm: np.ndarray
    si_index_table: np.ndarray


@functools.lru_cache(maxsize=1)
def load_modeldata() -> ModelData:
    """Load and cache the extracted modeldata tables for this process.

    The archive is packaged as ``pfsspecsim.etc.data/modeldata.npz``. The
    returned arrays are marked read-only so a single cached load can safely
    be shared across callers (and across forked worker processes) without
    risking in-place mutation corrupting the shared cache.
    """
    data_traversable = resources.files("pfsspecsim.etc.data") / "modeldata.npz"
    with resources.as_file(data_traversable) as npz_path:
        with np.load(npz_path, allow_pickle=False) as data:
            arrays = {name: data[name] for name in ModelData._fields}
    for array in arrays.values():
        array.flags.writeable = False
    return ModelData(**arrays)
