"""Silicon material properties (absorption length, real refractive index).

Pure-Python port of the material-property routines in ``src/gsetc.c``:

* :func:`si_abslength` <- ``gsOP_Si_abslength`` (gsetc.c:239-275),
  after Rajkanan, Singh & Shewchun (1979), Solid-State Electronics 22, 793.
* :func:`si_index_real` <- ``gsOP_Si_indexreal`` (gsetc.c:293-330),
  4-point Lagrange cubic interpolation of the Don Groom / LBNL silicon
  optical-constants table (extracted to ``modeldata.npz`` as
  ``si_index_table``, 5 nm grid starting at 200 nm).

Both functions accept a scalar or ndarray wavelength in nm and are
vectorized with numpy; scalar input yields a Python float. Wavelengths
outside the tabulated range raise :class:`ValueError` (the C code prints
to stderr and calls ``exit(1)``).
"""

from __future__ import annotations

import numpy as np

from ._modeldata import load_modeldata
from .constants import HC_EV_NM, K_B_EV_PER_K

# Rajkanan, Singh & Shewchun (1979) fit coefficients for silicon absorption
# (algorithm-specific literals, verbatim from gsetc.c:244-252).
_BETA = 7.021e-4  # eV/K
_GAMMA = 1108.0  # K
_E_G0 = (1.1557, 2.5)  # eV (indirect band gaps at T=0)
_E_GD0 = 3.2  # eV (direct band gap at T=0)
_E_P = (1.827e-2, 5.773e-2)  # eV (phonon energies)
_C = (5.5, 4.0)  # dimensionless
_A = (323.1, 7237.0)  # cm^-1 eV^-2
_A_D = 1.052e6  # cm^-1 eV^-2

# Legal wavelength range of the silicon data, in nm (the C code converts to
# meters and checks 2e-7 .. 1.1e-6 m; gsetc.c:256-260, 311-315).
_SI_LAMBDA_MIN_NM = 200.0
_SI_LAMBDA_MAX_NM = 1100.0

# si_index_table grid: lambda = 200 nm + 5 nm * i, i = 0..180 (gsetc.c:318).
_SI_TABLE_LAMBDA0_NM = 200.0
_SI_TABLE_DLAMBDA_NM = 5.0


def _check_si_range(lam_nm: np.ndarray) -> None:
    """Raise ValueError if any wavelength is outside 200-1100 nm."""
    if np.any(lam_nm < _SI_LAMBDA_MIN_NM) or np.any(lam_nm > _SI_LAMBDA_MAX_NM):
        raise ValueError(
            "Si data out of range: wavelength must be in "
            f"[{_SI_LAMBDA_MIN_NM}, {_SI_LAMBDA_MAX_NM}] nm"
        )


def si_abslength(lam_nm, temperature):
    """Silicon absorption length in microns (gsOP_Si_abslength, gsetc.c:239-275).

    Parameters
    ----------
    lam_nm : float or array_like
        Wavelength in nm; legal range 200-1100 nm.
    temperature : float
        Detector temperature in K.

    Returns
    -------
    float or ndarray
        Photon absorption length in silicon, in microns.
    """
    lam = np.asarray(lam_nm, dtype=np.float64)
    _check_si_range(lam)
    temp = float(temperature)

    hnu = HC_EV_NM / lam  # photon energy in eV
    # Temperature-dependent band-gap shift (same for direct and indirect).
    band_shift = _BETA * temp * temp / (temp + _GAMMA)
    e_gd = _E_GD0 - band_shift

    # Direct-gap contribution.
    de_d = hnu - e_gd
    alpha = np.where(de_d > 0.0, _A_D * np.sqrt(np.clip(de_d, 0.0, None)), 0.0)

    # Indirect (phonon-assisted) contributions: 2 phonons x 2 gaps, each with
    # phonon-absorption and phonon-emission terms (gsetc.c:267-272).
    kt = K_B_EV_PER_K * temp
    for i in range(2):
        occ_abs = 1.0 / (np.exp(_E_P[i] / kt) - 1.0)  # phonon absorption
        occ_emi = 1.0 / (1.0 - np.exp(-_E_P[i] / kt))  # phonon emission
        for j in range(2):
            e_g = _E_G0[j] - band_shift
            de = hnu - e_g + _E_P[i]
            alpha = alpha + np.where(de > 0.0, _C[i] * _A[j] * de * de * occ_abs, 0.0)
            de = hnu - e_g - _E_P[i]
            alpha = alpha + np.where(de > 0.0, _C[i] * _A[j] * de * de * occ_emi, 0.0)

    result = 1e4 / alpha  # alpha in cm^-1 -> length in microns
    return result if result.ndim else float(result)


def si_index_real(lam_nm):
    """Real part of the silicon refractive index (gsOP_Si_indexreal, gsetc.c:293-330).

    4-point Lagrange cubic interpolation on the tabulated index (5 nm grid
    from 200 nm), exact at the interior table nodes. Temperature dependence
    is vestigial in the C code (room-temperature table) and is dropped here.

    Parameters
    ----------
    lam_nm : float or array_like
        Wavelength in nm; legal range 200-1100 nm.

    Returns
    -------
    float or ndarray
        Real refractive index n of silicon.
    """
    lam = np.asarray(lam_nm, dtype=np.float64)
    _check_si_range(lam)
    table = load_modeldata().si_index_table

    x = (lam - _SI_TABLE_LAMBDA0_NM) / _SI_TABLE_DLAMBDA_NM
    xint = np.clip(np.floor(x).astype(np.int64), 1, 178)  # gsetc.c:319-321
    xfrac = x - xint

    # 4-point Lagrange cubic through table[xint-1 .. xint+2] (gsetc.c:324-329).
    f_m1 = table[xint - 1]
    f_0 = table[xint]
    f_p1 = table[xint + 1]
    f_p2 = table[xint + 2]
    result = (
        -xfrac * (xfrac - 1.0) * (xfrac - 2.0) / 6.0 * f_m1
        + (xfrac * xfrac - 1.0) * (xfrac - 2.0) / 2.0 * f_0
        - xfrac * (xfrac + 1.0) * (xfrac - 2.0) / 2.0 * f_p1
        + xfrac * (xfrac * xfrac - 1.0) / 6.0 * f_p2
    )
    return result if result.ndim else float(result)
