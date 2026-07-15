"""Night-sky emission-line tables, sky continuum model(s), and moonlight.

Pure-Python port of the sky-model portions of `gsGetNoise` in
`src/gsetc.c`:

* :func:`load_sky_lines` <- the UVES air->vacuum conversion done inline at
  gsetc.c:795 (``lambda = gs_air2vac(gsSKY_UVES_LAMBDA[iline]);``), packaging
  the extracted UVES + OH line tables (`_modeldata`, T1) for T8/noise.py.
* :func:`sky_continuum` <- the ``obs->skytype & 0xf`` sky-continuum-model
  switch (gsetc.c:889-953): models 0x0-0x6, including model 0x6's
  lambda<=800nm linear correction (gsetc.c:944).
* :func:`moon_continuum` <- the Krisciunas & Schaefer (1991) moonlight
  model (gsetc.c:956-1004), selected by ``(obs->skytype>>4)&0xf == 0x0``.

Both `sky_continuum` and `moon_continuum` are vectorized over wavelength
`lam_nm`; the observation-level angles/phase they take (airmass inputs,
moon geometry) are scalars, matching gsetc.c's `OBS_ATTRIB` struct (one set
of angles per run, many wavelengths per spectrum).

Line-list vs. continuum/moonlight split: the UVES and OH line contributions
are rescaled by *different* reference airmasses in gsGetNoise (UVES to 1.1,
gsetc.c:805-813; OH to 1.0, gsetc.c:842-852) and go through
`gsFracTrace`/`gsSpectroDist` (T6/psf.py, spectrograph-arm-specific), so
`load_sky_lines` only loads and unit-converts the raw line lists -- it does
not apply airmass scaling, effective area, or PSF distribution; that is
T8/noise.py's job, per-arm.
"""

from __future__ import annotations

import functools
from typing import NamedTuple

import numpy as np

from ._modeldata import load_modeldata
from .atmosphere import _sky_type_int, air2vac, cont_opacity
from .constants import DEG_TO_RAD, PHOTONS_PER_ERG_1NM

# --- sky_type bit-field helpers ---------------------------------------
#
# Same nibble layout as atmosphere.py's _opacity_model/_line_absorption_model
# (gsetc.c's `skytype` bitmask), for the two nibbles this module reads:
# continuum model (bits 0-3, gsetc.c:892) and moonlight model (bits 4-7,
# gsetc.c:962). The sky-line-model nibble (bits 16-19, gsetc.c:786) is not
# parsed here -- `load_sky_lines` deliberately leaves that switch to its
# caller (see its docstring); noise.py (T8), the only caller that needs it,
# defines its own `_line_model` helper.

_CONTINUUM_MODEL_MASK = 0xF
_MOON_MODEL_SHIFT = 4
_NIBBLE_MASK = 0xF


def _continuum_model(sky_type: int | str) -> int:
    return _sky_type_int(sky_type) & _CONTINUUM_MODEL_MASK


def _moon_model(sky_type: int | str) -> int:
    return (_sky_type_int(sky_type) >> _MOON_MODEL_SHIFT) & _NIBBLE_MASK


# --- Sky emission-line lists (gsetc.c:786-866) -----------------------------


class SkyLines(NamedTuple):
    """UVES + OH sky emission-line lists, ready for per-arm processing.

    Kept as two separate line lists (rather than concatenated) because
    gsGetNoise rescales them by different reference airmasses -- UVES to
    1.1 (gsetc.c:805-813: ``count *= airmass/1.1 * exp(...)``), OH to 1.0
    (gsetc.c:842-852: ``count *= airmass * exp(...) * exp((14.8-15.8)/1.086)``)
    -- a distinction T8 (noise.py) needs to preserve; this module does not
    apply either rescaling.

    Attributes
    ----------
    uves_lambda_vac_nm : ndarray, shape (2816,)
        UVES sky emission-line wavelengths, air->vacuum converted (nm)
        via `atmosphere.air2vac` (gsetc.c:795).
    uves_intensity : ndarray, shape (2816,)
        UVES line intensities, 1e-12 erg/m2/s/arcsec2 (gsetc.c:805-806),
        unmodified from the extracted table. A couple of entries are
        slightly negative (noise-floor artifacts in the UVES atlas);
        gsetc.c does not clip these here, only the derived per-pixel
        `count` further downstream (``if (count<0) count=0;``,
        gsetc.c:807) -- a T8/noise.py responsibility, not replicated in
        this loader.
    oh_lambda_vac_nm : ndarray, shape (698,)
        OH airglow line wavelengths (nm); already vacuum in the extracted
        table (T1), so no conversion is applied.
    oh_intensity : ndarray, shape (698,)
        OH line intensities, 1e-12 erg/m2/s/arcsec2 (gsetc.c:836-837).
    """

    uves_lambda_vac_nm: np.ndarray
    uves_intensity: np.ndarray
    oh_lambda_vac_nm: np.ndarray
    oh_intensity: np.ndarray


@functools.lru_cache(maxsize=1)
def load_sky_lines() -> SkyLines:
    """Load the UVES + OH sky emission-line lists (gsetc.c:786-866).

    Applies the UVES air->vacuum conversion gsGetNoise does inline per-line
    (gsetc.c:795); the OH table is already tabulated in vacuum wavelengths
    (T1) and is passed through unmodified.

    This function only loads and unit-converts the two line lists -- it
    does not implement the ``(obs->skytype>>16)&0xf`` sky-line-model switch
    (gsetc.c:786): model 0x0 omits sky lines entirely ("for testing only",
    gsetc.c:788-790) while model 0x1 (the only other case gsetc.c
    implements) always includes *both* UVES and OH together (gsetc.c:793,
    828 -- there is no mix-and-match). Callers that need to honor
    `sky_type`'s line-model bits (T8/noise.py) should decide whether to
    call this function / use its result themselves (e.g. skip it, or
    discard the result, when that nibble is 0x0) rather than have this
    loader silently return empty arrays for a model it was never asked
    about.

    Returns
    -------
    SkyLines
        Line-list namedtuple (read-only arrays); see `SkyLines` for field
        documentation. Cached process-wide (`functools.lru_cache`), like
        `_modeldata.load_modeldata`.
    """
    data = load_modeldata()
    uves_lambda_vac_nm = np.asarray(air2vac(data.uves_lambda), dtype=np.float64)
    uves_intensity = np.array(data.uves_int, dtype=np.float64)
    oh_lambda_vac_nm = np.array(data.oh_data[:, 0], dtype=np.float64)
    oh_intensity = np.array(data.oh_data[:, 1], dtype=np.float64)
    for array in (uves_lambda_vac_nm, uves_intensity, oh_lambda_vac_nm, oh_intensity):
        array.flags.writeable = False
    return SkyLines(
        uves_lambda_vac_nm=uves_lambda_vac_nm,
        uves_intensity=uves_intensity,
        oh_lambda_vac_nm=oh_lambda_vac_nm,
        oh_intensity=oh_intensity,
    )


# --- Sky continuum (part of the switch inside gsGetNoise, gsetc.c:889-953) -


def sky_continuum(lam_nm, airmass, trans, sky_type: int | str):
    """Sky continuum brightness, models 0x0-0x6 (gsetc.c:889-953).

    Bits 0-3 of `sky_type` select the model (``obs->skytype & 0xf``,
    gsetc.c:892). Returns the continuum in photons/s/m2/arcsec2/nm,
    *already* rescaled by `airmass` and `trans` -- exactly the value
    gsetc.c's `continuum` variable holds just before Moonlight
    (:func:`moon_continuum`) is added to it at gsetc.c:1003. (In the C
    source this rescaling -- ``continuum *= airmass * trans;`` -- is
    written inline within each non-0x0 case; hoisting it out to apply
    uniformly, including to the always-zero 0x0 case, is numerically
    identical for finite inputs.)

    Parameters
    ----------
    lam_nm : float or array_like
        Wavelength in nm.
    airmass : float
        Airmass, e.g. from `atmosphere.airmass` (gsetc.c:759).
    trans : float or array_like
        Atmospheric transmission at `lam_nm`, broadcastable with `lam_nm`.
        In gsetc.c this is a PSF-weighted average of `gsAtmTrans`
        (`atmosphere.transmission`) over the spectrograph PSF window
        (gsetc.c:876-886) -- sky.py has no access to the spectrograph PSF
        (T6/psf.py), so it is a required input here rather than computed
        internally; a caller not modeling the PSF smearing may simply pass
        ``atmosphere.transmission(lam_nm, zenith_ang, sky_type)``.
    sky_type : int or str
        Sky-model bitmask (int, or hex string as stored in
        `EtcParams.sky_type`).

    Returns
    -------
    float or ndarray
        Sky continuum, photons/s/m2/arcsec2/nm.
    """
    lam = np.asarray(lam_nm, dtype=np.float64)
    model = _continuum_model(sky_type)

    if model == 0x0:
        # No sky continuum at all, for testing only (gsetc.c:894-896).
        continuum = np.zeros_like(lam)

    elif model in (0x1, 0x2):
        # UVES continuum, typical airmass ~1.1, supplemented with IR
        # (gsetc.c:899-907).
        continuum = np.select(
            [
                lam < 375.0,
                lam < 483.0,
                lam < 580.0,
                lam < 674.5,
                lam < 858.0,
            ],
            [0.17, 0.14, 0.09, 0.10, 0.08],
            default=0.07,
        )
        continuum = continuum / 1.1
        if model == 0x1:
            ir_correction = 0.4669 * (1000.0 / lam) * (2.0 * lam / 1000.0 - 0.5)
            continuum = np.where(lam > 1040.0, ir_correction, continuum)
        continuum = continuum * 1e-11 * PHOTONS_PER_ERG_1NM * lam

    elif model == 0x3:
        # Fit to Jim Gunn spectrum (gsetc.c:909-914).
        mag = (
            21.55
            + (lam - 600.0) * np.where(lam > 600.0, 5e-5, -6e-3)
            - 0.55 * np.exp(-0.005 * (lam - 594.0) ** 2)
            - 0.175 * (1.0 + np.tanh(375.0 - lam))
            - 6.14656e9 / lam**4
        )
        continuum = 0.01089 * 10.0 ** (0.4 * (22.5 - mag)) * 1e6 / lam**2
        continuum = continuum * 1e-11 * PHOTONS_PER_ERG_1NM * lam

    elif model == 0x4:
        # Fit to Big Boss proposal spectrum (gsetc.c:917-921).
        continuum = 0.035 * np.sqrt(1000.0 / lam) + 0.045 * np.exp(
            -0.005 * (lam - 594.0) ** 2
        )
        continuum = continuum * 10.0 ** (0.4 * (0.05 + 6.14656e9 / lam**4))
        continuum = continuum * 1e-11 * PHOTONS_PER_ERG_1NM * lam

    elif model in (0x5, 0x6):
        # Fit to Jim Gunn spectrum, modified (gsetc.c:927-947); model 0x6
        # (the default) additionally applies a lambda<=800nm linear
        # correction (gsetc.c:944).
        mag = (
            np.where(lam > 600.0, 24.316, 27.166)
            + np.where(lam > 600.0, -5.199e-03, -1.419e-02) * lam
            + np.where(lam > 600.0, 1.465e-06, 8.541e-06) * lam**2
            - 0.55 * np.exp(-0.005 * (lam - 594.0) ** 2)
            - 6.14656e9 / lam**4
        )
        continuum = 0.01089 * 10.0 ** (0.4 * (22.5 - mag)) * 1e6 / lam**2
        if model == 0x6:
            low_lam_correction = -1.02278215e-03 * lam + 1.77400498e00
            continuum = continuum * np.where(lam > 800.0, 1.0, low_lam_correction)
        continuum = continuum * 1e-11 * PHOTONS_PER_ERG_1NM * lam

    else:
        raise ValueError(f"sky_continuum: illegal sky continuum model {model:#x}")

    continuum = continuum * airmass * trans
    return continuum if continuum.ndim else float(continuum)


# --- Moonlight (Krisciunas & Schaefer 1991), gsetc.c:956-1004 --------------


def moon_continuum(
    lam_nm,
    moon_zenith_ang,
    moon_target_ang,
    moon_phase,
    zenith_ang,
    sky_type: int | str,
):
    """Moonlight continuum brightness (Krisciunas & Schaefer 1991), gsetc.c:956-1004.

    Returns 0 when the Moon is below the horizon: `moon_zenith_ang >= 90`
    (the negation of gsetc.c:956's guard, ``if (obs->lunarZA<90)``). Bits
    4-7 of `sky_type` select the moonlight model
    (``(obs->skytype>>4)&0xf``, gsetc.c:962); only model 0x0 (Krisciunas &
    Schaefer) is implemented, the only one gsetc.c itself implements.

    The returned value is *not* rescaled by airmass or atmospheric
    transmission: gsetc.c adds it directly to the already-rescaled
    :func:`sky_continuum` value (``continuum += lunar_cont;``, gsetc.c:1003)
    -- unlike the sky continuum models, moonlight's own airmass-like
    dependence enters only through the two ``10**(-0.4*kV/sqrt(...))``
    extinction factors below (gsetc.c:987-988), not a separate multiply.

    Parameters
    ----------
    lam_nm : float or array_like
        Wavelength in nm.
    moon_zenith_ang : float
        Angle from the Moon to the zenith, in degrees (`obs->lunarZA`);
        >=90 means the Moon is below the horizon.
    moon_target_ang : float
        Angle from the Moon to the line of sight, in degrees
        (`obs->lunarangle`).
    moon_phase : float
        Lunar phase: 0.0 (new) -> 0.5 (full) -> 1.0 (new) (`obs->lunarphase`);
        only the fractional part is used (gsetc.c:959).
    zenith_ang : float
        Target's angle from the zenith, in degrees (`obs->zenithangle`).
    sky_type : int or str
        Sky-model bitmask (int, or hex string as stored in
        `EtcParams.sky_type`).

    Returns
    -------
    float or ndarray
        Moonlight continuum contribution, photons/s/m2/arcsec2/nm.

    Note
    ----
    gsetc.c:960 also computes ``k = gsAtmContOp(obs,lambda,flags)`` inside
    this block, but the result is never used anywhere in the function (only
    `kV`, evaluated at 550nm below, feeds the Krisciunas & Schaefer
    formula) -- dead code in the C source, intentionally not reproduced
    here.
    """
    lam = np.asarray(lam_nm, dtype=np.float64)

    if float(moon_zenith_ang) >= 90.0:
        result = np.zeros_like(lam)
        return result if result.ndim else float(result)

    model = _moon_model(sky_type)
    if model != 0x0:
        raise NotImplementedError(
            f"moon_continuum: moonlight model {model:#x} is not implemented "
            "(only model 0x0, Krisciunas & Schaefer 1991, is)"
        )

    lunarphase = moon_phase - np.floor(moon_phase)
    kV = cont_opacity(550.0, sky_type)

    # Solar spectrum rescaling (CFHT Redeye manual color model), normalized
    # to 550nm (gsetc.c:975-978).
    scale_rs = (np.exp(2480.0 / 550.0) - 1.0) / (np.exp(2480.0 / lam) - 1.0) * (
        lam / 550.0
    ) ** -7.0
    scale_ms = (np.exp(2480.0 / 550.0) - 1.0) / (np.exp(2480.0 / lam) - 1.0) * (
        lam / 550.0
    ) ** -4.3

    # V-band moonlight brightness model (gsetc.c:980-988).
    alpha = 360.0 * abs(lunarphase - 0.5)
    i_star = 10.0 ** (-0.4 * (3.84 + 0.026 * alpha + 4e-9 * alpha**4))
    if alpha < 7.0:
        i_star *= 1.35 - 0.05 * alpha

    moon_target_rad = float(moon_target_ang) * DEG_TO_RAD
    moon_za_rad = float(moon_zenith_ang) * DEG_TO_RAD
    target_za_rad = float(zenith_ang) * DEG_TO_RAD

    f1 = 2.29e5 * (1.06 + np.cos(moon_target_rad) ** 2)
    f2 = 10.0 ** (6.15 - float(moon_target_ang) / 40.0)

    b_moon = (
        (f1 * scale_rs + f2 * scale_ms)
        * i_star
        * 10.0 ** (-0.4 * kV / np.sqrt(1.0 - 0.96 * np.sin(moon_za_rad) ** 2))
        * (
            1.0
            - 10.0 ** (-0.4 * kV / np.sqrt(1.0 - 0.96 * np.sin(target_za_rad) ** 2))
        )
    )

    # Continuum conversion at V band (gsetc.c:990-994): 3.408e10 = nanoLambert
    # brightness of a 0th-magnitude object per arcsec^2; 5.48e10 = photons
    # per second per m^2 per ln(lambda) from a 0th-magnitude object. Both are
    # algorithm-defining literals from the C source, not derivable physical
    # constants (Vega-AB correction neglected, per the C comment).
    lunar_cont = b_moon / 3.408e10 * 5.48e10 / 550.0

    return lunar_cont if lunar_cont.ndim else float(lunar_cont)
