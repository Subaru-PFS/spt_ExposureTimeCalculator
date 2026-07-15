"""Air<->vacuum wavelength conversion, atmospheric continuum absorption,
transmission, and airmass.

Pure-Python port of the atmosphere-related routines in ``src/gsetc.c``:

* :func:`n_air` <- ``gs_n_air`` (gsetc.c:209-211), the Edlen (1953) index of
  refraction of air.
* :func:`vac2air` / :func:`air2vac` <- ``gs_vac2air`` / ``gs_air2vac``
  (gsetc.c:214-228); ``air2vac`` is a fixed 10-iteration fixed-point solve,
  verbatim (not a closed-form inverse).
* :func:`airmass` <- the ``1/sqrt(1-0.96*sin^2(za))`` approximation used at
  gsetc.c:759 (noise/throughput calculations elsewhere in the engine).
* :func:`cont_opacity` <- ``gsAtmContOp`` (gsetc.c:390-425); only the
  Mauna Kea median-extinction model (``(skytype>>8)&0xf == 0x0``) is
  implemented, per the task brief.
* :func:`transmission` <- ``gsAtmTrans`` (gsetc.c:431-491): a continuum term
  ``10**(-0.4*k/cos(za))`` -- note this is ``1/cos(za)``, *not* the 0.96
  airmass approximation used elsewhere in the engine -- times a line
  absorption factor selected by bits 12-15 of ``sky_type``: ``0x0`` uses the
  Kitt Peak transmission table only (with wavelengths below 500 nm passed
  through unattenuated, a verbatim quirk of the C index-clamp logic); ``0x1``
  uses the Kitt Peak table at and below 900 nm and switches to the Mauna Kea
  3mm precipitable-water table strictly above 900 nm (gsetc.c:462:
  ``if (lambda>900)``, so 900 nm itself still uses the Kitt Peak table).

All functions accept a scalar or array_like wavelength/angle and are
vectorized with numpy; scalar input yields a Python float. ``sky_type`` may
be passed as an ``int`` or as a hex string (as stored in
``EtcParams.sky_type``, e.g. ``"11006"``).
"""

from __future__ import annotations

import numpy as np

from ._modeldata import load_modeldata
from .constants import DEG_TO_RAD

# --- sky_type bit-field helpers --------------------------------------------
#
# `EtcParams.sky_type` is a hex string (gsetc.c's stdin-driven `skytype`
# bitmask); T2's params.py validates it is parseable as hex but does not
# split it into fields, so we parse the bits we need locally, matching the
# nibble layout used by gsAtmContOp/gsAtmTrans (gsetc.c:396, 440).

_OPACITY_MODEL_SHIFT = 8
_LINE_ABSORPTION_SHIFT = 12
_NIBBLE_MASK = 0xF


def _sky_type_int(sky_type: int | str) -> int:
    """Coerce `sky_type` (hex string or int) to a plain int bitmask."""
    if isinstance(sky_type, str):
        return int(sky_type, 16)
    return int(sky_type)


def _opacity_model(sky_type: int | str) -> int:
    return (_sky_type_int(sky_type) >> _OPACITY_MODEL_SHIFT) & _NIBBLE_MASK


def _line_absorption_model(sky_type: int | str) -> int:
    return (_sky_type_int(sky_type) >> _LINE_ABSORPTION_SHIFT) & _NIBBLE_MASK


# --- Conversion functions (gsetc.c:206-228) --------------------------------


def n_air(lam_vac_nm):
    """Index of refraction of air (Edlen 1953 formula; gs_n_air, gsetc.c:209-211).

    Parameters
    ----------
    lam_vac_nm : float or array_like
        Vacuum wavelength in nm.

    Returns
    -------
    float or ndarray
        Index of refraction of air at `lam_vac_nm`.
    """
    lam = np.asarray(lam_vac_nm, dtype=np.float64)
    inv_lam2 = 1e6 / lam / lam
    result = (
        1.000064328
        + 0.0294981 / (146.0 - inv_lam2)
        + 0.0002554 / (41.0 - inv_lam2)
    )
    return result if result.ndim else float(result)


def vac2air(lam_vac_nm):
    """Vacuum -> air wavelength conversion (gs_vac2air, gsetc.c:214-216)."""
    lam = np.asarray(lam_vac_nm, dtype=np.float64)
    result = lam / n_air(lam)
    return result if result.ndim else float(result)


def air2vac(lam_air_nm):
    """Air -> vacuum wavelength conversion (gs_air2vac, gsetc.c:219-228).

    Verbatim port of the C routine's fixed 10-iteration fixed-point solve
    (not a closed-form inverse of `vac2air`): each iteration nudges the
    current vacuum-wavelength estimate by the residual
    ``lambda_air - vac2air(lambda_vac)``.
    """
    lam_air = np.asarray(lam_air_nm, dtype=np.float64)
    lam_vac = lam_air.copy()
    for _ in range(10):
        lam_vac = lam_vac + (lam_air - vac2air(lam_vac))
    return lam_vac if lam_vac.ndim else float(lam_vac)


def airmass(zenith_ang_deg):
    """Airmass from zenith angle, ``1/sqrt(1-0.96*sin^2(za))`` (gsetc.c:759).

    Parameters
    ----------
    zenith_ang_deg : float or array_like
        Zenith angle in degrees.

    Returns
    -------
    float or ndarray
        Airmass.
    """
    za = np.asarray(zenith_ang_deg, dtype=np.float64) * DEG_TO_RAD
    result = 1.0 / np.sqrt(1.0 - 0.96 * np.sin(za) ** 2)
    return result if result.ndim else float(result)


# --- Continuum atmospheric opacity (gsAtmContOp, gsetc.c:390-425) ---------
#
# Mauna Kea median-extinction model (opacity model 0x0): Boulade 1988 (CFHT
# Bulletin 19, 16) below 400nm, CFHT Observer's Manual above. The C code
# expresses this as a cascade of unconditional `if (lambda<T) k = ...;`
# statements with strictly decreasing thresholds T, so for any given
# lambda the *last* (most restrictive / smallest-T) true condition wins;
# this is equivalent to a standard piecewise assignment keyed on the
# smallest threshold each wavelength falls under, which is how it's
# expressed below via `np.select` (most-restrictive condition first).
def _kp_opacity(lam_nm: np.ndarray) -> np.ndarray:
    conditions = [
        lam_nm < 310.0,
        lam_nm < 320.0,
        lam_nm < 340.0,
        lam_nm < 360.0,
        lam_nm < 380.0,
        lam_nm < 400.0,
        lam_nm < 450.0,
        lam_nm < 500.0,
        lam_nm < 550.0,
        lam_nm < 600.0,
        lam_nm < 650.0,
        lam_nm < 700.0,
        lam_nm < 800.0,
        lam_nm < 900.0,
    ]
    choices = [
        1.37,
        1.37 - 0.55 * (lam_nm - 310.0) / 10.0,
        0.82 - 0.31 * (lam_nm - 320.0) / 20.0,
        0.51 - 0.14 * (lam_nm - 340.0) / 20.0,
        0.37 - 0.07 * (lam_nm - 360.0) / 20.0,
        0.30 - 0.05 * (lam_nm - 380.0) / 20.0,
        0.25 - 0.08 * (lam_nm - 400.0) / 50.0,
        0.17 - 0.04 * (lam_nm - 450.0) / 50.0,
        0.13 - 0.01 * (lam_nm - 500.0) / 50.0,
        0.12 - 0.01 * (lam_nm - 550.0) / 50.0,
        0.11 - 0.00 * (lam_nm - 600.0) / 50.0,
        0.11 - 0.01 * (lam_nm - 650.0) / 50.0,
        0.10 - 0.03 * (lam_nm - 700.0) / 100.0,
        0.07 - 0.02 * (lam_nm - 800.0) / 100.0,
    ]
    return np.select(conditions, choices, default=0.05)


def cont_opacity(lam_nm, sky_type: int | str):
    """Continuum atmospheric opacity in mag/airmass (gsAtmContOp, gsetc.c:390-425).

    Only opacity model 0x0 (bits 8-11 of `sky_type`), the Mauna Kea median
    extinction curve, is implemented -- the only model gsetc.c itself
    implements.

    Parameters
    ----------
    lam_nm : float or array_like
        Wavelength in nm.
    sky_type : int or str
        Sky-model bitmask (as an int, or a hex string as stored in
        `EtcParams.sky_type`).

    Returns
    -------
    float or ndarray
        Opacity `k` in magnitudes per airmass.
    """
    model = _opacity_model(sky_type)
    if model != 0x0:
        raise NotImplementedError(
            f"cont_opacity: opacity model {model:#x} is not implemented "
            "(only model 0x0, Mauna Kea median extinction, is)"
        )
    lam = np.asarray(lam_nm, dtype=np.float64)
    result = _kp_opacity(lam)
    return result if result.ndim else float(result)


# --- Line absorption tables (Kitt Peak / Mauna Kea 3mm) --------------------
#
# Grids per etc/data/README.md (T1): atm_trans_kp is a 40001-point table on
# a 0.025nm grid starting at 500nm (500-1500nm); mk_trans_3mm is a
# 30001-point table on a 0.02nm grid starting at 900nm (900-1500nm).
_KP_LAMBDA0_NM = 500.0
_KP_DLAMBDA_NM = 0.025
_MK_LAMBDA0_NM = 900.0
_MK_DLAMBDA_NM = 0.02


def _kp_line_trans(lam_nm: np.ndarray) -> np.ndarray:
    """Kitt Peak line-absorption transmission factor.

    Below 500nm the C code's index clamp (`xint>=0` guard, gsetc.c:447,
    479) leaves the transmission unmultiplied, i.e. factor 1 -- a verbatim
    quirk, not a physical statement that there's no absorption there.
    Above the table's upper edge the C code exits with an error; we raise
    `ValueError` instead, using the table's actual last tabulated
    wavelength (`lam_max`, derived from `kp.size`) as the cutoff. This is
    marginally more permissive than the C guard (`xint>39998`, gsetc.c:451),
    which stops one grid step early and so never actually reads the
    table's last element -- a ~0.025nm sliver just below 1500nm that no
    physically reachable input hits, so the difference is not observable
    in practice.
    """
    kp = load_modeldata().atm_trans_kp
    lam_max = _KP_LAMBDA0_NM + _KP_DLAMBDA_NM * (kp.size - 1)
    if np.any(lam_nm > lam_max):
        raise ValueError(
            "transmission: line absorption model not valid in the IR "
            f"(lambda > {lam_max} nm)"
        )
    lam_grid = _KP_LAMBDA0_NM + _KP_DLAMBDA_NM * np.arange(kp.size)
    factor = np.ones_like(lam_nm)
    below = lam_nm < _KP_LAMBDA0_NM
    factor[~below] = np.interp(lam_nm[~below], lam_grid, kp)
    return factor


def _mk_line_trans(lam_nm: np.ndarray) -> np.ndarray:
    """Mauna Kea 3mm-precipitable-water line-absorption transmission factor."""
    mk = load_modeldata().mk_trans_3mm
    lam_max = _MK_LAMBDA0_NM + _MK_DLAMBDA_NM * (mk.size - 1)
    if np.any(lam_nm > lam_max):
        raise ValueError(
            "transmission: line absorption model not valid in the IR "
            f"(lambda > {lam_max} nm)"
        )
    lam_grid = _MK_LAMBDA0_NM + _MK_DLAMBDA_NM * np.arange(mk.size)
    return np.interp(lam_nm, lam_grid, mk)


def _line_absorption(lam_nm: np.ndarray, sky_type: int | str) -> np.ndarray:
    model = _line_absorption_model(sky_type)
    if model == 0x0:
        return _kp_line_trans(lam_nm)
    if model == 0x1:
        factor = np.empty_like(lam_nm)
        hi = lam_nm > _MK_LAMBDA0_NM
        factor[hi] = _mk_line_trans(lam_nm[hi])
        factor[~hi] = _kp_line_trans(lam_nm[~hi])
        return factor
    raise ValueError(
        f"transmission: unrecognized atmospheric line absorption model {model:#x}"
    )


def transmission(lam_nm, zenith_ang_deg, sky_type: int | str):
    """Atmospheric transmission (gsAtmTrans, gsetc.c:431-491).

    Combines the continuum extinction ``10**(-0.4*k/cos(za))`` (note:
    ``1/cos(za)``, *not* the 0.96 airmass approximation of :func:`airmass`
    -- gsetc.c:437) with a line-absorption factor selected by bits 12-15 of
    `sky_type` (see :func:`_line_absorption`).

    Parameters
    ----------
    lam_nm : float or array_like
        Wavelength in nm.
    zenith_ang_deg : float
        Zenith angle in degrees (scalar, as in gsetc.c's `OBS_ATTRIB`).
    sky_type : int or str
        Sky-model bitmask (as an int, or a hex string as stored in
        `EtcParams.sky_type`).

    Returns
    -------
    float or ndarray
        Atmospheric transmission (dimensionless, 0-1).
    """
    lam = np.asarray(lam_nm, dtype=np.float64)
    za_rad = float(zenith_ang_deg) * DEG_TO_RAD
    k = cont_opacity(lam, sky_type)
    trans = 10.0 ** (-0.4 * k / np.cos(za_rad))
    trans = trans * _line_absorption(lam, sky_type)
    return trans if trans.ndim else float(trans)
