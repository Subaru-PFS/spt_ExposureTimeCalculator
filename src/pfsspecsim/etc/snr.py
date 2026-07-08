"""Signal and various SNR computations (line, [OII] doublet, continuum).

Pure-Python port of `gsGetSignal`, `gsGetSNR`, `gsGetSNR_Single`,
`gsGetSNR_OII`, and `gsGetSNR_Continuum` (`src/gsetc.c:1063-1446`).

All functions here compute a *single-exposure* SNR (or signal), exactly like
their C counterparts -- the ``sqrt(n_exp)`` "coadd `n_exp` exposures"
rescaling seen throughout `gsetc.c`'s `main` (e.g. gsetc.c:2033-2034,
2069-2070, 2107) is applied by the caller (the T10 engine), not here.

Public API
----------
* :func:`get_signal` <- `gsGetSignal` (gsetc.c:1071-1118): the per-feature
  signal vector for one spectral line (flux `F`, velocity dispersion
  `sigma_v`), vectorized over an array of lines/redshifts. Rather than the
  C function's full `Npix`-length `Signal` array (mostly zero), this
  returns only the `NP_WIN`-pixel window that can ever be nonzero, plus the
  window's reference (leftmost) pixel index -- see the "Vectorization"
  section below.
* :func:`snr_line` <- `gsGetSNR` (gsetc.c:1130-1163): reduces a `get_signal`
  window against an externally supplied per-pixel noise vector to an SNR
  (`snr_type` 0 = 1D optimal, 1 = uniform matched filter).
* :func:`snr_single` <- `gsGetSNR_Single` (gsetc.c:1174-1257): as
  `snr_line`, but adds the emission line's host-continuum shot noise
  (magnitude `mag`) to the noise vector first.
* :func:`snr_oii` <- `gsGetSNR_OII` (gsetc.c:1273-1360): SNR for the [OII]
  doublet (two lines, flux-ratio `ROII`, both with continuum shot noise
  added), `snr_type` 2 additionally combining the two lines into a single
  "optimal" reduction over their *union* pixel window (see below).
* :func:`snr_continuum` <- `gsGetSNR_Continuum` (gsetc.c:1367-1446): the
  full-pixel-vectorized continuum SNR curve for a magnitude (or
  `MagSpec`-evaluated per-pixel-wavelength magnitude) source; returns all 6
  of the C function's output curves as a :class:`ContinuumSNR`.

`params: EtcParams` supplies the observing-condition/instrument-setup
scalars that are constant across an entire engine run in gsetc.c (`decent`
== `params.fiber_offset`, `fieldang` == `params.field_ang`, `t_exp` ==
`params.exp_time`, `zenith_ang`, `sky_type`, `galactic_ext` == EBV,
`seeing_fwhm_800` == `params.seeing`, `hgcdte_sutr`) -- mirroring
`noise.compute_noise_arm(params, spectro, i_arm)`'s calling convention.
`r_eff` (and, for `snr_oii`, `F`/`ROII`/`src_cont`) remain explicit
arguments rather than being read from `params`: the [OII] catalog use case
(gsetc.c:2147-2166, a later task) calls `gsGetSNR_OII` with a *per-catalog-row*
`r_eff`/`F`/`ROII`/`src_cont`, distinct from `params.reff`/`params.line_flux`.

Vectorization
-------------
`get_signal`'s C original allocates a full `Npix`-length `Signal` array per
call but only ever populates an `NP_WIN=32`-pixel window of it (gsetc.c:1076
comment); summing `Signal**2/Noise` (or any of the other per-pixel
reductions used downstream) over the full array is therefore identical to
summing over just that window, since `Signal` is exactly zero everywhere
else. `get_signal` here returns only the `(Nz, NP_WIN)` window (`iref` is
the window's reference/leftmost pixel, one per line) and `snr_line`/
`snr_single` index the (shared, full-length) noise vector with
`iref[:, None] + arange(NP_WIN)` to build the matching `(Nz, NP_WIN)` noise
window -- the `(Nz, Npix)` tensor implied by a naive vectorization of
`Signal` is never materialized.

`snr_oii`'s `snr_type == 2` branch needs `(Signal0 + Signal1)**2 / Noise`
summed over *both* lines' pixel ranges, but the two lines' `NP_WIN` windows
are independently positioned (`iref0 != iref1` in general -- the doublet
separation, projected to pixels, is nonzero) and so cannot simply be added
window-against-window. Both windows are deposited into a shared,
`_OII_WIN = 48`-pixel-wide *union* window instead (`iref_u = min(iref0,
iref1)`, each line's contribution placed at its own offset from `iref_u`);
see the `_OII_MAX_OFFSET` derivation in code for why 48 pixels always
suffices given each line's own `iref` is already clipped to
`[0, Npix - NP_WIN]` by `get_signal`.

Like `psf.py` (see its module docstring re: `constants.L_CHUNK`), this
module does not chunk internally over large line/redshift arrays; callers
sweeping e.g. the full [OII] curve or catalog redshift grids (potentially
several thousand entries -- gsetc.c:2031, 2067, 2129) should slice their
input arrays into `constants.Z_CHUNK`-sized pieces themselves if peak memory
(dominated by `psf.spectro_mtf`'s `(Nz, 1000)` MTF array, transitively
computed inside `get_signal`/`frac_trace`/`effective_area`) matters.

Preserved quirks (not "fixed" -- see task brief)
-------------------------------------------------
* `snr_single`'s object-continuum transmission is a 41-point Gaussian
  average evaluated at a *fixed* wavelength (`lambda`) for every
  quadrature point (gsetc.c:1213-1218), unlike `get_signal`'s average
  (evaluated at `lambda*(1 + x*sigma_v/c)`, i.e. genuinely smeared by the
  line's velocity dispersion). Averaging a constant over any weights
  reproduces that constant, so this is ported as the mathematically
  equivalent direct `atmosphere.transmission(lambda, ...)` call rather than
  a literal (but pointless) 41-point loop -- see :func:`snr_single`.
* `snr_continuum` evaluates each pixel's wavelength at the pixel's *left
  edge* (`lmin + dl*ipix`, gsetc.c:1401), whereas `noise.compute_noise_arm`'s
  sky continuum grid uses the pixel *center* (`lmin + dl*(ipix+0.5)`,
  gsetc.c:882) -- a genuine cross-function inconsistency in gsetc.c, kept
  verbatim.
* The [OII]/single-line curves' diagnostic `gsGeometricThroughput` printf
  columns in gsetc.c's `main` use `fieldang=0` for the OII curve
  (gsetc.c:2049) but the real `fieldang` for the single-line curve
  (gsetc.c:2084); the [OII] catalog loop reads (but never uses) an input
  `sigma` column, always passing a hardcoded `sigma_v=70` instead
  (gsetc.c:2154). Both are properties of the `main`-loop driver, not of any
  function in this module, and are noted here only for traceability -- they
  belong to the (later) engine task that ports that driver loop.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from . import psf
from .atmosphere import transmission as atm_transmission
from .config import Spectrograph
from .constants import (
    AB_ZEROPOINT_CGS,
    C_KM_PER_S,
    C_NM_PER_S,
    NP_WIN,
    OII_LAMBDA,
    PHOTONS_PER_ERG_1NM,
)
from .extinction import alambda_over_ebv
from .noise import sample_factor_for_arm, smoothed_transmission
from .params import EtcParams

#: 41-point Gaussian quadrature grid shared by `get_signal`'s and
#: `snr_oii`'s object-transmission averages (gsetc.c:1100-1105/1308-1312:
#: ``for(x=-4;x<4.01;x+=.2)``, 41 points). Built via `arange`, not repeated
#: `+=0.2` accumulation, to avoid float drift relative to the C loop.
_TRANS_AVG_X = np.arange(41, dtype=np.float64) * 0.2 - 4.0
_TRANS_AVG_W = np.exp(-0.5 * _TRANS_AVG_X**2)
_TRANS_AVG_WSUM = _TRANS_AVG_W.sum()

#: Absolute-count threshold gating the uniform-matched-filter SNR's
#: division (gsetc.c:1147/1155/1249/1257/1345/1353: ``numer>=0.001? ... :
#: 0``) -- an algorithm-defining literal from gsetc.c, not a derived
#: constant.
_SNR_UNIFORM_FLOOR = 0.001

#: [OII] doublet union-window width for `snr_oii`'s `snr_type == 2` (see
#: module docstring); not a gsetc.c literal (the C function materializes
#: the full `Npix`-length `Signal0`/`Signal1`/`myNoise` arrays instead) but
#: a vectorization parameter of this port.
_OII_WIN = 48
#: Max allowed |iref0 - iref1|; `snr_oii` raises if this is ever exceeded
#: rather than silently truncating a line's window (see module docstring
#: derivation). Verified empirically <=12 pixels across all packaged
#: PFS.*.dat configs over the full z=0.1..2.7627 sweep.
_OII_MAX_OFFSET = _OII_WIN - NP_WIN


def _reduce_window(
    signal_window: np.ndarray, noise_window: np.ndarray, snr_type: int
) -> np.ndarray:
    """Reduce a `(Nz, W)` signal/noise window pair to an SNR.

    `gsGetSNR`'s two `snrType` branches (gsetc.c:1139-1155), factored out
    so `snr_line`, `snr_single`, and `snr_oii` (whose "noise" differs --
    bare `Noise`, `Noise + continuum`, or an [OII] union-window slice) can
    all share the reduction.
    """
    if snr_type == 0:
        return np.sqrt(np.sum(signal_window**2 / noise_window, axis=-1))
    if snr_type == 1:
        numer = np.sum(signal_window**2, axis=-1)
        denom = np.sum(signal_window**2 * noise_window, axis=-1)
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = numer / np.sqrt(denom)
        return np.where(numer >= _SNR_UNIFORM_FLOOR, ratio, 0.0)
    raise ValueError(f"snr_type must be 0 (optimal) or 1 (uniform), got {snr_type!r}")


def _velocity_smeared_transmission(
    lam: np.ndarray, sigma_v: np.ndarray, zenith_ang: float, sky_type: int | str
) -> np.ndarray:
    """`get_signal`'s object-transmission average (gsetc.c:1100-1105).

    41-point Gaussian average of the atmospheric transmission over
    ``lambda*(1 + x*sigma_v/c)`` for ``x`` in `_TRANS_AVG_X`.
    """
    offset = sigma_v / C_KM_PER_S  # (Nz,), fractional velocity shift per unit x
    lam_grid = lam[:, None] * (1.0 + _TRANS_AVG_X[None, :] * offset[:, None])  # (Nz,41)
    trans_grid = atm_transmission(lam_grid, zenith_ang, sky_type)
    return np.sum(trans_grid * _TRANS_AVG_W[None, :], axis=1) / _TRANS_AVG_WSUM


def _oii_doublet_transmission(
    lam0: np.ndarray, lam1: np.ndarray, zenith_ang: float, sky_type: int | str
) -> np.ndarray:
    """`snr_oii`'s object-continuum transmission average (gsetc.c:1308-1312).

    41-point Gaussian average of the atmospheric transmission over
    ``lambda0 + (0.5 + 0.5*x)*(lambda1-lambda0)`` for ``x`` in
    `_TRANS_AVG_X` -- centered between, and spanning somewhat beyond, the
    two doublet components (*not* the same "smear by velocity dispersion"
    scheme as `_velocity_smeared_transmission`).
    """
    diff = lam1 - lam0  # (Nz,)
    lam_grid = lam0[:, None] + (0.5 + 0.5 * _TRANS_AVG_X[None, :]) * diff[:, None]
    trans_grid = atm_transmission(lam_grid, zenith_ang, sky_type)
    return np.sum(trans_grid * _TRANS_AVG_W[None, :], axis=1) / _TRANS_AVG_WSUM


def _continuum_counts_at(
    params: EtcParams,
    spectro: Spectrograph,
    i_arm: int,
    lam_eval: np.ndarray,
    r_eff: float,
    src_cont: np.ndarray,
    trans: np.ndarray,
    apply_sample_factor: bool = True,
) -> np.ndarray:
    """Detected per-pixel continuum counts from a source of flux density
    `src_cont` (erg/cm2/s/Hz) at `lam_eval`, shared by `snr_single`,
    `snr_oii`, and `snr_continuum` (gsetc.c:1219-1231/1310-1320/1428-1435
    -- the same formula, evaluated at each call site's own `lambda`/`ll`
    and `trans`).

    Includes the per-Hz -> per-pixel conversion (`C_NM_PER_S*dl/lambda**2`)
    and, if `apply_sample_factor` (the default), the HgCdTe-SUTR x1.2 gate
    (`sample_factor_for_arm`) -- `snr_single`/`snr_oii` fold this directly
    into `counts` (gsetc.c:1227-1229/1317-1319: ``counts *= 1.2``), but
    `snr_continuum` keeps its `sample_factor` separate from the returned
    `.counts` curve (gsetc.c:1440-1441), so it passes
    `apply_sample_factor=False` and applies the factor itself afterwards.
    """
    dl = spectro.dl[i_arm]
    counts = (
        src_cont
        * trans
        * 10.0 ** (-0.4 * alambda_over_ebv(lam_eval) * params.galactic_ext)
        * psf.geometric_throughput(
            spectro,
            lam_eval,
            r_eff,
            params.fiber_offset,
            params.field_ang,
            params.seeing,
        )
        * psf.frac_trace(spectro, i_arm, lam_eval, tr=0)
        * PHOTONS_PER_ERG_1NM
        * lam_eval
        * params.exp_time
        * psf.effective_area(spectro, i_arm, lam_eval, params.field_ang)
        * 1e4
    )
    if apply_sample_factor:
        counts = counts * sample_factor_for_arm(spectro, i_arm, params.hgcdte_sutr)
    counts = counts * C_NM_PER_S * dl / (lam_eval**2)
    return counts


def get_signal(
    params: EtcParams,
    spectro: Spectrograph,
    i_arm: int,
    lam,
    F,
    sigma_v,
    r_eff: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-feature signal window, `gsGetSignal` (gsetc.c:1071-1118).

    Parameters
    ----------
    params : EtcParams
    spectro : Spectrograph
    i_arm : int
        Internal 0-based arm index.
    lam : array_like, shape (Nz,)
        Feature wavelength(s), nm.
    F : array_like or float, shape (Nz,) or scalar
        Flux, erg/cm2/s.
    sigma_v : array_like or float, shape (Nz,) or scalar
        Velocity dispersion, km/s.
    r_eff : float
        Source effective radius, arcsec (0 for a point source).

    Returns
    -------
    signal_window : ndarray, shape (Nz, NP_WIN)
        Detected counts per pixel within each feature's `NP_WIN`-pixel
        window (all-zero for a line whose window falls entirely outside
        `[-NP_WIN, Npix)`, gsetc.c:1092-1093).
    iref : ndarray, shape (Nz,), int64
        Reference (leftmost) pixel index of each window, clipped to
        `[0, Npix - NP_WIN]` (gsetc.c:1094-1095) regardless of validity.
    """
    lam = np.atleast_1d(np.asarray(lam, dtype=np.float64))
    nz = lam.size
    F = np.broadcast_to(np.asarray(F, dtype=np.float64), (nz,))
    sigma_v = np.broadcast_to(np.asarray(sigma_v, dtype=np.float64), (nz,))

    npix = int(spectro.npix[i_arm])
    lmin = spectro.lmin[i_arm]
    dl = spectro.dl[i_arm]

    pos = (lam - lmin) / dl
    iref_raw = np.floor(pos - NP_WIN / 2.0).astype(np.int64)
    valid = (iref_raw >= -NP_WIN) & (iref_raw < npix)
    iref = np.clip(iref_raw, 0, npix - NP_WIN)

    trans = _velocity_smeared_transmission(
        lam, sigma_v, params.zenith_ang, params.sky_type
    )

    counts = (
        F
        * trans
        * 10.0 ** (-0.4 * alambda_over_ebv(lam) * params.galactic_ext)
        * psf.geometric_throughput(
            spectro, lam, r_eff, params.fiber_offset, params.field_ang, params.seeing
        )
        * psf.frac_trace(spectro, i_arm, lam, tr=0)
        * PHOTONS_PER_ERG_1NM
        * lam
        * params.exp_time
        * psf.effective_area(spectro, i_arm, lam, params.field_ang)
        * 1e4
    )

    sigma_pix = sigma_v / C_KM_PER_S * lam / dl
    fr = psf.spectro_dist(
        spectro, i_arm, lam, pos - iref, sigma_pix, NP_WIN
    )  # (Nz, NP_WIN)

    signal_window = fr * counts[:, None]
    signal_window = np.where(valid[:, None], signal_window, 0.0)
    return signal_window, iref


def snr_line(
    params: EtcParams,
    spectro: Spectrograph,
    i_arm: int,
    lam,
    F,
    sigma_v,
    r_eff: float,
    noise: np.ndarray,
    snr_type: int = 0,
):
    """SNR for a spectral feature given a precomputed noise vector,
    `gsGetSNR` (gsetc.c:1130-1163).

    Parameters
    ----------
    params, spectro, i_arm, lam, F, sigma_v, r_eff
        As `get_signal`.
    noise : ndarray, shape (Npix,)
        Per-pixel noise variance for this arm (`ArmNoise.noise`).
    snr_type : int
        0 = 1D optimal, 1 = uniform matched filter.

    Returns
    -------
    float or ndarray, shape (Nz,)
        SNR (scalar if `lam` was scalar).
    """
    lam_arr = np.asarray(lam, dtype=np.float64)
    scalar_input = lam_arr.ndim == 0

    signal_window, iref = get_signal(params, spectro, i_arm, lam_arr, F, sigma_v, r_eff)
    noise = np.asarray(noise, dtype=np.float64)
    noise_window = noise[iref[:, None] + np.arange(NP_WIN)[None, :]]
    snr = _reduce_window(signal_window, noise_window, snr_type)
    return float(snr[0]) if scalar_input else snr


def snr_single(
    params: EtcParams,
    spectro: Spectrograph,
    i_arm: int,
    mag,
    lam,
    F,
    sigma_v,
    r_eff: float,
    noise: np.ndarray,
    snr_type: int = 0,
):
    """SNR for a single emission line atop its host continuum,
    `gsGetSNR_Single` (gsetc.c:1174-1257).

    The legacy C function's ``mag==-99.9``-flagged, hand-rolled per-pixel
    magnitude-file interpolation (gsetc.c:1187-1210) is superseded by
    `params.MagSpec`: callers pass the already-resolved AB magnitude(s) at
    `lam` (e.g. ``MagSpec(params.mag, params.mag_file)(lam)``) via `mag`
    directly, scalar or per-line (`mag` may vary with `lam`/z when driven
    by a wavelength-dependent `mag_file`).

    Parameters
    ----------
    params, spectro, i_arm, lam, F, sigma_v, r_eff, noise, snr_type
        As `snr_line`.
    mag : array_like or float, shape (Nz,) or scalar
        AB magnitude of the host continuum at `lam`.

    Returns
    -------
    float or ndarray, shape (Nz,)
        SNR (scalar if `lam` was scalar).
    """
    lam_arr = np.asarray(lam, dtype=np.float64)
    scalar_input = lam_arr.ndim == 0
    lam_arr = np.atleast_1d(lam_arr)
    nz = lam_arr.size
    mag_arr = np.broadcast_to(np.asarray(mag, dtype=np.float64), (nz,))

    # QUIRK (gsetc.c:1213-1218), preserved verbatim -- not "fixed": the
    # 41-point Gaussian transmission average here evaluates
    # `gsAtmTrans(obs,lambda,flags)` -- the *same*, unshifted `lambda` --
    # at every quadrature point, unlike get_signal's genuinely
    # velocity-smeared average. The weighted average of a constant is that
    # constant, so this collapses to the plain pointwise transmission;
    # ported as the equivalent direct call rather than a literal (and
    # pointless) 41-point loop over an x-independent integrand. See the
    # module docstring.
    trans = atm_transmission(lam_arr, params.zenith_ang, params.sky_type)

    src_cont = AB_ZEROPOINT_CGS * 10.0 ** (-0.4 * mag_arr)
    counts = _continuum_counts_at(
        params, spectro, i_arm, lam_arr, r_eff, src_cont, trans
    )

    signal_window, iref = get_signal(params, spectro, i_arm, lam_arr, F, sigma_v, r_eff)
    noise = np.asarray(noise, dtype=np.float64)
    noise_window = noise[iref[:, None] + np.arange(NP_WIN)[None, :]] + counts[:, None]
    snr = _reduce_window(signal_window, noise_window, snr_type)
    return float(snr[0]) if scalar_input else snr


def snr_oii(
    params: EtcParams,
    spectro: Spectrograph,
    i_arm: int,
    z,
    F,
    sigma_v,
    r_eff: float,
    src_cont,
    ROII,
    noise: np.ndarray,
    snr_type: int = 2,
):
    """SNR for the [OII] 3727 doublet, `gsGetSNR_OII` (gsetc.c:1273-1360).

    Parameters
    ----------
    params, spectro, i_arm, sigma_v, r_eff, noise
        As `snr_line` (`sigma_v`/`r_eff` are shared by both doublet lines).
    z : array_like or float, shape (Nz,) or scalar
        Redshift.
    F : array_like or float, shape (Nz,) or scalar
        Total (both lines) doublet flux, erg/cm2/s.
    src_cont : array_like or float, shape (Nz,) or scalar
        Host-continuum flux density at the doublet, erg/cm2/s/Hz.
    ROII : array_like or float, shape (Nz,) or scalar
        Flux ratio of the two lines; clamped to [0.667, 3.87]
        (gsetc.c:1290-1291) -- the ETC always runs with `ROII=1`, per the
        comment there, but the clamp is preserved for API fidelity.
    snr_type : int
        0 = 1D optimal (brighter line), 1 = uniform matched filter
        (brighter line), 2 = combined 1D optimal of both lines (the only
        `snrType` gsetc.c's `main` actually drives for [OII], gsetc.c:2022).

    Returns
    -------
    float or ndarray, shape (Nz,)
        SNR (scalar if `z` was scalar).
    """
    z_arr = np.asarray(z, dtype=np.float64)
    scalar_input = z_arr.ndim == 0
    z_arr = np.atleast_1d(z_arr)
    nz = z_arr.size

    lam0 = OII_LAMBDA[0] * (1.0 + z_arr)
    lam1 = OII_LAMBDA[1] * (1.0 + z_arr)

    roii_arr = np.clip(
        np.broadcast_to(np.asarray(ROII, dtype=np.float64), (nz,)), 0.667, 3.87
    )
    frac0 = roii_arr / (1.0 + roii_arr)
    frac1 = 1.0 / (1.0 + roii_arr)

    f_arr = np.broadcast_to(np.asarray(F, dtype=np.float64), (nz,))
    src_cont_arr = np.broadcast_to(np.asarray(src_cont, dtype=np.float64), (nz,))

    ll = 0.5 * (lam0 + lam1)
    trans = _oii_doublet_transmission(lam0, lam1, params.zenith_ang, params.sky_type)
    counts = _continuum_counts_at(
        params, spectro, i_arm, ll, r_eff, src_cont_arr, trans
    )

    noise = np.asarray(noise, dtype=np.float64)
    npix = int(spectro.npix[i_arm])

    signal0, iref0 = get_signal(
        params, spectro, i_arm, lam0, frac0 * f_arr, sigma_v, r_eff
    )
    signal1, iref1 = get_signal(
        params, spectro, i_arm, lam1, frac1 * f_arr, sigma_v, r_eff
    )

    if snr_type in (0, 1):
        win = np.arange(NP_WIN)[None, :]
        noise_win0 = noise[iref0[:, None] + win] + counts[:, None]
        noise_win1 = noise[iref1[:, None] + win] + counts[:, None]
        snr0 = _reduce_window(signal0, noise_win0, snr_type)
        snr1 = _reduce_window(signal1, noise_win1, snr_type)
        snr = np.maximum(snr0, snr1)
    elif snr_type == 2:
        iref_u = np.minimum(np.minimum(iref0, iref1), npix - _OII_WIN)
        offset0 = iref0 - iref_u
        offset1 = iref1 - iref_u
        if np.any(offset0 > _OII_MAX_OFFSET) or np.any(offset1 > _OII_MAX_OFFSET):
            raise ValueError(
                "snr_oii: [OII] doublet pixel separation exceeds the "
                f"{_OII_WIN}-pixel union window (max offset {_OII_MAX_OFFSET}); "
                "increase _OII_WIN for this spectrograph configuration"
            )

        rows = np.arange(nz)[:, None]
        combined = np.zeros((nz, _OII_WIN), dtype=np.float64)
        win = np.arange(NP_WIN)[None, :]
        combined[rows, offset0[:, None] + win] += signal0
        combined[rows, offset1[:, None] + win] += signal1

        noise_window = (
            noise[iref_u[:, None] + np.arange(_OII_WIN)[None, :]] + counts[:, None]
        )
        snr = np.sqrt(np.sum(combined**2 / noise_window, axis=1))
    else:
        raise ValueError(f"snr_type must be 0, 1, or 2, got {snr_type!r}")

    return float(snr[0]) if scalar_input else snr


@dataclasses.dataclass
class ContinuumSNR:
    """Per-pixel continuum SNR curve, `gsGetSNR_Continuum`'s 6 output
    arrays (gsetc.c:1437-1443).

    Attributes
    ----------
    snr : ndarray, shape (Npix,)
        ``counts / sqrt(sample_factor*counts + noise)``.
    counts : ndarray, shape (Npix,)
        Detected continuum counts per pixel.
    noise : ndarray, shape (Npix,)
        ``sample_factor*counts + noise`` (total per-pixel variance).
    mag : ndarray, shape (Npix,)
        AB magnitude used at each pixel (echoes the input, broadcast).
    trans : ndarray, shape (Npix,)
        ``counts / src_cont`` (net system throughput, dimensionless).
    sample_factor : ndarray, shape (Npix,)
        The (pixel-independent, but returned per-pixel to match the C
        output layout) HgCdTe-SUTR sampling factor.
    """

    snr: np.ndarray
    counts: np.ndarray
    noise: np.ndarray
    mag: np.ndarray
    trans: np.ndarray
    sample_factor: np.ndarray


def snr_continuum(
    params: EtcParams,
    spectro: Spectrograph,
    i_arm: int,
    mag,
    r_eff: float,
    noise: np.ndarray,
) -> ContinuumSNR:
    """Continuum SNR curve over every pixel of one arm, `gsGetSNR_Continuum`
    (gsetc.c:1367-1446).

    Parameters
    ----------
    params, spectro, i_arm, r_eff
        As `snr_line`.
    mag : array_like or float, shape (Npix,) or scalar
        AB magnitude of the source continuum at each pixel's wavelength
        (see `snr_single`'s docstring re: `MagSpec` superseding the legacy
        `mag==-99.9` file-interpolation flag).
    noise : ndarray, shape (Npix,)
        Per-pixel noise variance for this arm (`ArmNoise.noise`).

    Returns
    -------
    ContinuumSNR
    """
    npix = int(spectro.npix[i_arm])
    dl = spectro.dl[i_arm]
    lmin = spectro.lmin[i_arm]

    # QUIRK (gsetc.c:1401), preserved verbatim -- not "fixed": each pixel's
    # wavelength here is its *left edge* (`lmin + dl*ipix`), unlike
    # `noise.compute_noise_arm`'s sky-continuum grid, which uses the pixel
    # *center* (`lmin + dl*(ipix+0.5)`, gsetc.c:882) -- a genuine
    # cross-function inconsistency in gsetc.c itself. See module docstring.
    lam_pix = lmin + dl * np.arange(npix, dtype=np.float64)

    mag_arr = np.broadcast_to(np.asarray(mag, dtype=np.float64), (npix,)).astype(
        np.float64
    )

    sample_factor = sample_factor_for_arm(spectro, i_arm, params.hgcdte_sutr)

    trans = smoothed_transmission(
        spectro, i_arm, lam_pix, params.zenith_ang, params.sky_type
    )

    src_cont = AB_ZEROPOINT_CGS * 10.0 ** (-0.4 * mag_arr)
    counts = _continuum_counts_at(
        params,
        spectro,
        i_arm,
        lam_pix,
        r_eff,
        src_cont,
        trans,
        apply_sample_factor=False,
    )

    noise = np.asarray(noise, dtype=np.float64)
    out_noise = sample_factor * counts + noise
    with np.errstate(invalid="ignore", divide="ignore"):
        out_snr = counts / np.sqrt(out_noise)

    return ContinuumSNR(
        snr=out_snr,
        counts=counts,
        noise=out_noise,
        mag=mag_arr,
        trans=counts / src_cont,
        sample_factor=np.full(npix, sample_factor),
    )
