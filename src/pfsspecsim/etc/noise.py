"""Per-arm noise vector construction.

Pure-Python port of ``gsGetNoise`` (`src/gsetc.c:730-1059`), the single most
important regression gate in the whole port: for a given spectrograph arm it
builds the per-pixel noise variance (``Noise``, e^2/pix/exposure) and the
sky-only counts vector (``SkyMod``, here called ``sky``, e/pix/exposure)
that ``tests/master_results/noise_omp.dat`` records for comparison.

Pipeline (mirrors the C function's structure exactly):

1. **Sky emission lines** (gsetc.c:786-871): UVES (optical, air->vacuum
   converted, referenced to airmass 1.1) and OH (near-IR, referenced to
   airmass 1.0, with an extra ``exp((14.8-15.8)/1.086)`` brightness-level
   rescaling) line lists, each deposited into a ``SP_PSF_LEN`` (32-pixel)
   window per line via `psf.spectro_dist`, weighted by the line's total
   captured flux (`psf.frac_trace` with ``tr=1``, i.e. including the two
   neighboring traces) and effective area (`psf.effective_area`).
2. **Sky continuum** (gsetc.c:876-1009): a per-pixel-center-wavelength
   continuum model (`sky.sky_continuum`) plus moonlight
   (`sky.moon_continuum`), remapped through a PSF-smoothed atmospheric
   transmission (:func:`smoothed_transmission`), then converted to counts
   the same way as the line contributions.
3. **Sky-subtraction systematic** (gsetc.c:1019-1026): a floor proportional
   to the square of the local (3-pixel-max) sky level.
4. **Diffuse stray light** (gsetc.c:1029-1034): a small fraction of the
   mean sky level, spread uniformly across the arm.
5. **Dark current + read noise** (gsetc.c:1037-1040).
6. **Empirical adjustment** (gsetc.c:1042-1055): a fixed per-arm scale
   factor (`constants.ADJUST_NOISE_LR`/`ADJUST_NOISE_MR`), fit to past
   on-sky observations.

All per-pixel accumulation happens in `Noise`; `sky` is captured *before*
steps 3-6 are added (gsetc.c:1012-1015: ``sky[ipix]=Noise[ipix]/sample_factor``
right after the line+continuum loops, and never touched again) and is
*never* multiplied by the empirical adjustment factor -- both are verbatim
quirks of the C source, not simplifications made here.

``sample_factor`` (1.0, or 1.2 for HgCdTe SUTR arms when
`EtcParams.hgcdte_sutr` is set, gsetc.c:753-756) is a pre-factor on the
*expected Poisson mean* fed to the noise variance -- it multiplies every
photon-counting contribution (lines, continuum, diffuse stray, dark
current) but not read noise or the sky-subtraction systematic, and it also
divides back out of `sky` (gsetc.c:1013) so that `sky` is a photon-count
estimate, not a variance.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from . import psf
from ._parallel import map_arms
from .atmosphere import _sky_type_int, airmass, cont_opacity
from .atmosphere import transmission as atm_transmission
from .config import Spectrograph, field_interp
from .constants import (
    ADJUST_NOISE_LR,
    ADJUST_NOISE_MR,
    ARCSEC_PER_URAD,
    MAG_TO_NEP_APPROX,
    PHOTONS_PER_ERG_1NM,
    SP_PSF_LEN,
)
from .params import EtcParams
from .sky import load_sky_lines, moon_continuum, sky_continuum

#: Reference airmass the UVES sky-line atlas is normalized to (gsetc.c:810-813).
_UVES_REF_AIRMASS = 1.1
#: Reference airmass the OH airglow line list is normalized to (gsetc.c:847-852).
_OH_REF_AIRMASS = 1.0
#: OH brightness-level rescaling: the OH table is tabulated at 14.8 mag/as2
#: (Vega); the desired sky level is 15.8 mag/as2 (UKIRT User Guide),
#: gsetc.c:848-852.
_OH_BRIGHTNESS_SCALE = np.exp((14.8 - 15.8) / MAG_TO_NEP_APPROX)

#: Bit shift/mask for the sky-line-model nibble of `sky_type` (bits 16-19,
#: gsetc.c:786): ``0x0`` = no sky lines at all (testing only), ``0x1`` =
#: UVES + OH together (no mix-and-match; the only non-trivial model
#: gsetc.c implements).
_LINE_MODEL_SHIFT = 16
_NIBBLE_MASK = 0xF


def _line_model(sky_type: int | str) -> int:
    return (_sky_type_int(sky_type) >> _LINE_MODEL_SHIFT) & _NIBBLE_MASK


@dataclasses.dataclass
class ArmNoise:
    """Per-arm noise/sky vectors, `gsGetNoise`'s per-call output.

    Attributes
    ----------
    lam : ndarray, shape (Npix,)
        Pixel-center wavelength, nm (``lmin + (ipix+0.5)*dl``).
    noise : ndarray, shape (Npix,)
        Noise variance, e^2/pix/exposure (``Noise`` in gsetc.c).
    sky : ndarray, shape (Npix,)
        Sky-only counts, e/pix/exposure (``SkyMod`` in gsetc.c) -- captured
        before the sky-subtraction-systematic/diffuse-stray/dark-read/
        empirical-adjustment terms are folded into `noise`.
    sample_factor : float
        Poisson pre-sampling factor for this arm (1.0, or 1.2 for HgCdTe
        SUTR arms; gsetc.c:753-756). Engine (T10) needs this to combine
        the source-signal Poisson variance with `noise` consistently.
    """

    lam: np.ndarray
    noise: np.ndarray
    sky: np.ndarray
    sample_factor: float


@dataclasses.dataclass
class NoiseResult:
    """Noise/sky vectors for every spectrograph arm, indexed by internal `ia`."""

    arms: list[ArmNoise]

    def __getitem__(self, ia: int) -> ArmNoise:
        return self.arms[ia]

    def __len__(self) -> int:
        return len(self.arms)


def sample_factor_for_arm(
    spectro: Spectrograph, i_arm: int, hgcdte_sutr: bool
) -> float:
    """Poisson pre-sampling factor for one arm (gsetc.c:753-756).

    1.2 for HgCdTe (`Dtype[i_arm] == 1`) arms when `hgcdte_sutr` is set
    (the C engine's ``#ifdef HGCDTE_SUTR`` compile flag, now a runtime
    `EtcParams` field); 1.0 otherwise.
    """
    if hgcdte_sutr and int(spectro.Dtype[i_arm]) == 1:
        return 1.2
    return 1.0


def _fiber_radius_arcsec(spectro: Spectrograph, fieldang: float) -> float:
    """Fiber radius on sky, in arcsec (gsetc.c:767-780).

    ``EFL`` (effective focal length) is the same 5-node field-angle
    interpolation as `psf.geometric_throughput`'s `efl`/`sigma`
    (`config.field_interp`), reused verbatim here (gsetc.c:768-779 is a
    third, textually distinct, open-coding of the identical formula).
    """
    efl = field_interp(spectro.EFL, fieldang, spectro.rfov)
    return spectro.fiber_ent_rad / efl * ARCSEC_PER_URAD


def smoothed_transmission(
    spectro: Spectrograph,
    i_arm: int,
    lam_pix: np.ndarray,
    zenith_ang: float,
    sky_type: int | str,
    mtf: np.ndarray | None = None,
) -> np.ndarray:
    """PSF-smoothed atmospheric transmission at each pixel (gsetc.c:880-887).

    The raw `atmosphere.transmission` line-absorption features are much
    narrower than the spectrograph PSF; gsetc.c remaps the transmission
    through the (sigma=0) spectrograph PSF response `FR` before using it to
    rescale the sky continuum model, by resampling `FR` (computed on the
    ``SP_PSF_LEN``-pixel grid) onto a 5x finer sub-pixel grid via integer
    division (``FR[j/5]``, C truncating division) and taking a weighted
    average of the transmission over that grid.

    Parameters
    ----------
    spectro : Spectrograph
    i_arm : int
        Internal 0-based arm index.
    lam_pix : ndarray, shape (Npix,)
        Pixel-center wavelengths, nm.
    zenith_ang : float
        Zenith angle, degrees.
    sky_type : int or str
        Sky-model bitmask.
    mtf : ndarray, shape (Npix, Nu), optional
        Precomputed `psf.spectro_mtf(spectro, i_arm, lam_pix, psf.U_GRID)`;
        computed internally if omitted. Callers that also need
        `psf.frac_trace` at the same `(spectro, i_arm, lam_pix)` (as
        `compute_noise_arm` does, for the continuum-count conversion)
        should pass it in to share the one expensive MTF evaluation.

    Returns
    -------
    ndarray, shape (Npix,)
        Transmission, PSF-smoothed and averaged, at each pixel.
    """
    lam_pix = np.atleast_1d(np.asarray(lam_pix, dtype=np.float64))
    dl = spectro.dl[i_arm]
    n = SP_PSF_LEN
    pos = n / 2 - 0.5  # 15.5: centers the SP_PSF_LEN-pixel window.

    fr = psf.spectro_dist(spectro, i_arm, lam_pix, pos, 0.0, n, mtf=mtf)  # (Npix, 32)

    j = np.arange(5 * n)  # 0..159
    weight_idx = j // 5  # C's `FR[j/5]` truncating integer division.
    offsets_nm = (0.2 * j - n / 2 + 0.5) * dl  # (160,)
    lam_grid = lam_pix[:, None] + offsets_nm[None, :]  # (Npix, 160)

    trans_grid = atm_transmission(lam_grid, zenith_ang, sky_type)  # (Npix, 160)
    fr_weighted = fr[:, weight_idx]  # (Npix, 160)

    num = np.sum(fr_weighted * trans_grid, axis=1)
    den = np.sum(fr_weighted, axis=1)
    return num / den


def _line_deposit_positions(pos: np.ndarray, npix: int) -> np.ndarray:
    """Reference (leftmost) pixel index for a `SP_PSF_LEN`-window deposit.

    ``iref = floor(pos - (SP_PSF_LEN/2 - 0.5))``, clipped to
    ``[0, Npix - SP_PSF_LEN]`` (gsetc.c:815-817 for UVES; gsetc.c:857-859
    for OH -- the same formula, ``floor(pos-16+0.5) == floor(pos-15.5)``,
    written slightly differently in the two C call sites).
    """
    iref = np.floor(pos - (SP_PSF_LEN / 2 - 0.5)).astype(np.int64)
    return np.clip(iref, 0, npix - SP_PSF_LEN)


def _line_counts(
    spectro: Spectrograph,
    i_arm: int,
    lam_vac: np.ndarray,
    intensity: np.ndarray,
    fieldang: float,
    rad: float,
    t_exp: float,
    am: float,
    sky_type: int | str,
    ref_airmass: float,
    brightness_scale: float,
    mtf: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-line detected photon count, before spatial deposit (gsetc.c:805-813, 842-852).

    Common to both the UVES and OH line lists, which differ only in their
    reference airmass and (for OH only) an extra brightness-level rescaling.

    Returns
    -------
    count : ndarray, shape (L,)
        Detected counts per line, post ``count<0 -> 0`` clamp and airmass/
        extinction/brightness rescaling.
    frac : ndarray, shape (L,)
        `psf.frac_trace(..., tr=1)` for these lines (also needed by the
        caller to build the deposit's MTF-sharing `spectro_dist` call).
    """
    lam_vac = np.asarray(lam_vac, dtype=np.float64)
    intensity = np.asarray(intensity, dtype=np.float64)
    if mtf is None:
        mtf = psf.spectro_mtf(spectro, i_arm, lam_vac, psf.U_GRID)
    frac = psf.frac_trace(spectro, i_arm, lam_vac, tr=1, mtf=mtf)
    aeff = psf.effective_area(spectro, i_arm, lam_vac, fieldang)

    count = (
        intensity
        * lam_vac
        * 1e-12
        * PHOTONS_PER_ERG_1NM
        * frac
        * aeff
        * t_exp
        * np.pi
        * rad**2
    )
    count = np.where(count < 0.0, 0.0, count)

    count = (
        count
        * (am / ref_airmass)
        * np.exp(-cont_opacity(lam_vac, sky_type) * am / MAG_TO_NEP_APPROX)
        * brightness_scale
    )
    return count, frac


def _deposit_lines(
    noise: np.ndarray,
    spectro: Spectrograph,
    i_arm: int,
    lam_vac: np.ndarray,
    intensity: np.ndarray,
    fieldang: float,
    rad: float,
    t_exp: float,
    am: float,
    sky_type: int | str,
    ref_airmass: float,
    brightness_scale: float,
    sample_factor: float,
) -> None:
    """Deposit one sky-line list's contribution into `noise`, in place.

    Ports the per-line-list loop body shared by UVES (gsetc.c:794-830) and
    OH (gsetc.c:833-864): selects lines within range of this arm
    (gsetc.c:797/838), computes their counts (:func:`_line_counts`), and
    spreads each line's count over a `SP_PSF_LEN`-pixel window
    (:func:`_line_deposit_positions` + `psf.spectro_dist`) via
    `numpy.add.at` (safe under repeated/overlapping indices, unlike plain
    fancy-index assignment).
    """
    npix = int(spectro.npix[i_arm])
    lmin = spectro.lmin[i_arm]
    dl = spectro.dl[i_arm]

    lam_vac = np.asarray(lam_vac, dtype=np.float64)
    pos = (lam_vac - lmin) / dl
    in_range = (pos > -(SP_PSF_LEN / 2)) & (pos < npix + SP_PSF_LEN / 2 - 1)
    if not np.any(in_range):
        return

    lam_m = lam_vac[in_range]
    pos_m = pos[in_range]
    intensity_m = np.asarray(intensity, dtype=np.float64)[in_range]

    # One spectro_mtf evaluation, reused for both frac_trace(tr=1) (the
    # total-captured-flux normalization) and the SP_PSF_LEN-pixel deposit
    # window below (gsetc.c reuses the same GsSpectroDist-underlying MTF
    # implicitly by calling gsFracTrace then gsSpectroDist separately with
    # the same lambda; here we make the shared work explicit).
    mtf = psf.spectro_mtf(spectro, i_arm, lam_m, psf.U_GRID)

    count, _frac = _line_counts(
        spectro,
        i_arm,
        lam_m,
        intensity_m,
        fieldang,
        rad,
        t_exp,
        am,
        sky_type,
        ref_airmass,
        brightness_scale,
        mtf=mtf,
    )

    iref = _line_deposit_positions(pos_m, npix)
    fr = psf.spectro_dist(spectro, i_arm, lam_m, pos_m - iref, 0.0, SP_PSF_LEN, mtf=mtf)

    j = np.arange(SP_PSF_LEN)
    np.add.at(noise, iref[:, None] + j[None, :], (count * sample_factor)[:, None] * fr)


def _continuum_counts(
    spectro: Spectrograph,
    i_arm: int,
    params: EtcParams,
    rad: float,
    am: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-pixel sky-continuum(+moonlight) detected counts (gsetc.c:876-1009).

    Returns
    -------
    lam_pix : ndarray, shape (Npix,)
        Pixel-center wavelengths, nm.
    count : ndarray, shape (Npix,)
        Detected counts per pixel (not yet scaled by `sample_factor`).
    """
    npix = int(spectro.npix[i_arm])
    lmin = spectro.lmin[i_arm]
    dl = spectro.dl[i_arm]
    lam_pix = lmin + (np.arange(npix, dtype=np.float64) + 0.5) * dl

    # Shared MTF for the whole pixel grid: reused by smoothed_transmission
    # (the SP_PSF_LEN-window remap) and frac_trace (the total captured-flux
    # fraction) below -- (Npix, 1000) f8, ~32MB peak for Npix=4096
    # (constants.L_CHUNK), per the task brief's memory-discipline note.
    mtf = psf.spectro_mtf(spectro, i_arm, lam_pix, psf.U_GRID)

    trans = smoothed_transmission(
        spectro, i_arm, lam_pix, params.zenith_ang, params.sky_type, mtf=mtf
    )
    continuum = sky_continuum(lam_pix, am, trans, params.sky_type)
    if params.moon_zenith_ang < 90.0:
        continuum = continuum + moon_continuum(
            lam_pix,
            params.moon_zenith_ang,
            params.moon_target_ang,
            params.moon_phase,
            params.zenith_ang,
            params.sky_type,
        )

    aeff = psf.effective_area(spectro, i_arm, lam_pix, params.field_ang)
    frac = psf.frac_trace(spectro, i_arm, lam_pix, tr=1, mtf=mtf)

    count = continuum * dl * aeff * params.exp_time * np.pi * rad**2 * frac
    return lam_pix, count


def compute_noise_arm(params: EtcParams, spectro: Spectrograph, i_arm: int) -> ArmNoise:
    """Noise/sky vectors for one spectrograph arm, `gsGetNoise` (gsetc.c:730-1059).

    Parameters
    ----------
    params : EtcParams
        Resolved ETC input parameters.
    spectro : Spectrograph
        Parsed spectrograph configuration (already built by the caller --
        `resolve_degrade`/`load_spectrograph_config` are the caller's
        responsibility, not this function's).
    i_arm : int
        Internal 0-based arm index.

    Returns
    -------
    ArmNoise
    """
    npix = int(spectro.npix[i_arm])
    samp = sample_factor_for_arm(spectro, i_arm, params.hgcdte_sutr)
    am = airmass(params.zenith_ang)
    rad = _fiber_radius_arcsec(spectro, params.field_ang)

    noise = np.zeros(npix, dtype=np.float64)

    line_model = _line_model(params.sky_type)
    if line_model == 0x1:
        lines = load_sky_lines()
        _deposit_lines(
            noise,
            spectro,
            i_arm,
            lines.uves_lambda_vac_nm,
            lines.uves_intensity,
            params.field_ang,
            rad,
            params.exp_time,
            am,
            params.sky_type,
            ref_airmass=_UVES_REF_AIRMASS,
            brightness_scale=1.0,
            sample_factor=samp,
        )
        _deposit_lines(
            noise,
            spectro,
            i_arm,
            lines.oh_lambda_vac_nm,
            lines.oh_intensity,
            params.field_ang,
            rad,
            params.exp_time,
            am,
            params.sky_type,
            ref_airmass=_OH_REF_AIRMASS,
            brightness_scale=_OH_BRIGHTNESS_SCALE,
            sample_factor=samp,
        )
    elif line_model != 0x0:
        raise ValueError(f"compute_noise: illegal sky line model {line_model:#x}")

    lam_pix, cont_count = _continuum_counts(spectro, i_arm, params, rad, am)
    noise += cont_count * samp

    # sky (SkyMod) is captured here -- before the systematic/stray/dark/read
    # terms below, and never touched by the empirical adjustment at the end
    # (gsetc.c:1012-1015).
    sky = noise / samp

    # Sky-subtraction systematic: floor proportional to the local (3-pixel,
    # edge-clipped) max sky level (gsetc.c:1019-1026). `sysfrac` itself
    # already carries the sqrt(n_exp) "don't average down" rescaling
    # (gsetc.c:1834), applied once at the parameter level, not per-arm here.
    sky_left = np.concatenate(([-np.inf], sky[:-1]))
    sky_right = np.concatenate((sky[1:], [-np.inf]))
    sky_sysref = np.maximum(np.maximum(sky, sky_left), sky_right)
    sysfrac = params.sky_sub_floor * np.sqrt(params.exp_num)
    noise = noise + (sysfrac**2) * sky_sysref**2

    # Diffuse stray light background (gsetc.c:1029-1034).
    stray_ref = (
        sky.sum()
        * spectro.width[i_arm]
        * spectro.pix[i_arm]
        / spectro.sep[i_arm]
        / npix
    )
    noise = noise + params.diffuse_stray * stray_ref * samp

    # Dark current & read noise (gsetc.c:1037-1040).
    var = (
        spectro.dark[i_arm] * params.exp_time * samp + spectro.read[i_arm] ** 2
    ) * spectro.width[i_arm]
    noise = noise + var

    # Empirical adjustment from past observations, indexed by the internal
    # arm index `ia` (gsetc.c:1042-1055) -- *not* the output arm id.
    adjust = ADJUST_NOISE_MR if spectro.MR else ADJUST_NOISE_LR
    noise = noise * adjust[i_arm]

    return ArmNoise(lam=lam_pix, noise=noise, sky=sky, sample_factor=samp)


def compute_noise(params: EtcParams, spectro: Spectrograph) -> NoiseResult:
    """Noise/sky vectors for every spectrograph arm (`gsGetNoise`, called
    once per arm by gsetc.c's `main`, gsetc.c:2004-2006).

    Arms are independent, so this maps `compute_noise_arm` across
    `params.n_workers` threads (`_parallel.map_arms`), collecting results in
    `ia` order for bit-identical output regardless of `n_workers`. When the
    sky-line model is enabled (`_line_model(params.sky_type) == 0x1`, the
    same gate `compute_noise_arm` itself checks), `load_sky_lines()` is
    warmed up once here, serially, before the parallel region: it is
    `lru_cache`d, and doing the first (file-loading) call from a single
    thread avoids relying on `lru_cache`'s own thread-safety for the
    concurrent first calls that would otherwise all race on a cache miss.
    """
    if _line_model(params.sky_type) == 0x1:
        load_sky_lines()
    return NoiseResult(
        arms=map_arms(
            lambda ia: compute_noise_arm(params, spectro, ia),
            spectro.N_arms,
            params.n_workers,
        )
    )
