"""Galactic dust extinction ratio A_lambda / E(B-V).

Pure-Python port of ``gsGalactic_Alambda__EBV`` (``src/gsetc.c:343-385``):
linear interpolation of the Milky Way R_V=3.1 dust extinction curve of
Weingartner & Draine (kext_albedo_WD_MW_3.1_60_D03.all), extracted to
``modeldata.npz`` as ``dust_norm`` (201 points, log-spaced in wavelength
from 10 um down to 0.1 um: ``lambda_um = 10 * 10**(-i/100)``).

The function accepts a scalar or ndarray wavelength in nm and is vectorized
with numpy; scalar input yields a Python float. Wavelengths outside
0.1-10 um raise :class:`ValueError` (the C code prints to stderr and calls
``exit(1)``).
"""

from __future__ import annotations

import numpy as np

from ._modeldata import load_modeldata

# Valid range of the dust extinction table, in nm (0.1-10 um; gsetc.c:372).
_DUST_LAMBDA_MIN_NM = 100.0
_DUST_LAMBDA_MAX_NM = 10000.0


def alambda_over_ebv(lam_nm):
    """A_lambda / E(B-V) for Milky Way dust (gsGalactic_Alambda__EBV, gsetc.c:343-385).

    Parameters
    ----------
    lam_nm : float or array_like
        Wavelength in nm; legal range 100-10000 nm (0.1-10 um).

    Returns
    -------
    float or ndarray
        Extinction-to-reddening ratio A_lambda / E(B-V), in magnitudes.
    """
    lam = np.asarray(lam_nm, dtype=np.float64)
    if np.any(lam < _DUST_LAMBDA_MIN_NM) or np.any(lam > _DUST_LAMBDA_MAX_NM):
        raise ValueError(
            "Wavelength out of range for dust extinction law: must be in "
            f"[{_DUST_LAMBDA_MIN_NM}, {_DUST_LAMBDA_MAX_NM}] nm"
        )
    norm = load_modeldata().dust_norm

    lam_um = lam / 1e3  # nm -> um (gsetc.c:370)
    # Table index is log-spaced, 100 points per dex from 10 um downward
    # (gsetc.c:378-382).
    lset = 100.0 / np.log(10.0) * np.log(10.0 / lam_um)
    il = np.clip(np.floor(lset).astype(np.int64), 0, 199)
    lf = lset - il

    result = norm[il] + lf * (norm[il + 1] - norm[il])
    return result if result.ndim else float(result)
