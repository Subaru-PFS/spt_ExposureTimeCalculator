"""Geometric throughput, effective area, spectrograph MTF, spectral PSF
distribution, and fiber trace fraction.

Pure-Python port of the throughput / PSF routines in ``src/gsetc.c``:

* :func:`geometric_throughput` <- ``gsGeometricThroughput`` (gsetc.c:504-570):
  encircled energy of a point/extended source in a fiber, combining the
  Kolmogorov seeing MTF, a leading-order diffraction term, and an
  exponential-profile galaxy convolved with a circular fiber aperture
  (Hankel-transform overlap integral, 50-point midpoint quadrature).
* :func:`effective_area` <- ``gsAeff`` (gsetc.c:576-615): geometric area x
  vignetting x wavelength-dependent throughput (piecewise-linear
  interpolation of the per-arm throughput grid, clamped at the endpoints).
* :func:`spectro_mtf` <- ``gsSpectroMTF`` (gsetc.c:621-667): the 1D Fourier
  transform (real part only) of the spectrograph PSF in the dispersion
  direction, as a function of wavenumber `u` (cycles/pixel): fiber-aperture
  Airy factor, pixel sinc, Gaussian camera spot, Si detector defocus
  (6-point depth average; Si arms only), and grating scattering.
* :func:`spectro_dist` <- ``gsSpectroDist`` (gsetc.c:673-690): the fraction
  of a spectral feature's flux landing in each pixel of an `N`-pixel
  window, i.e. the inverse Fourier transform of `spectro_mtf` (times an
  extra Gaussian smearing term) evaluated on the C code's fixed 1000-point
  `u`-grid.
* :func:`frac_trace` <- ``gsFracTrace`` (gsetc.c:696-710): total fraction of
  flux captured in a fixed-width extraction window, optionally including
  the wings of the two adjacent traces (`tr=1`).

Bessel functions ``getJ0``/``getJ1`` (gsetc.c:115-152, Abramowitz & Stegun
approximations) are replaced by ``scipy.special.j0``/``j1``, per the task
brief; ``geterf`` (gsetc.c) is never called from any of the above and is not
ported.

Vectorization
-------------
Every function accepts (and, apart from `spectro_mtf`'s `u` argument, is
*meant* to be called with) a wavelength array of length `L` -- e.g. one
entry per emission line in a sky-line list. `spectro_mtf` returns a dense
`(L, Nu)` MTF array (`Nu=1000`, the C code's fixed quadrature grid,
`u = (arange(1000)+0.5)*0.005`). `spectro_dist` turns each row of that
array into an `(L, N)` array of per-pixel fractions using the
angle-addition identity ``cos(2*pi*u*(pos-ip)) = cos(2*pi*u*pos)*cos(2*pi*u*ip)
+ sin(2*pi*u*pos)*sin(2*pi*u*ip)``, which factors the C code's `Nu`-point
sum into two `(L, Nu) @ (Nu, N)` matrix products -- the `(Nlines, Nu, Npix)`
tensor implied by a naive vectorization of the C triple loop is never
materialized (see the memory-discipline note in the task brief). Likewise
the 6-depth-point Si-defocus loop in `spectro_mtf` accumulates into an
`(L, Nu)` pair of running sums rather than materializing `(L, 6, Nu)`.

When `pos` and `sigma` are both scalars (one shared feature position and
smearing width for every wavelength -- e.g. `frac_trace`'s window centers,
or `noise.smoothed_transmission`'s fixed half-window `pos`), the identity
factors further: every lambda-independent term folds into a single
`(Nu, N)` weight matrix ``W[iu, ip] = 2*du * exp(-2*pi^2*sigma^2*u^2) *
cos(2*pi*u*(pos-ip))`` and the whole distribution collapses to one
``mtf @ W`` product -- no `(L, Nu)` trig/exp evaluation at all.
`frac_trace` goes one step more: its sum over the `N`-pixel window and the
`2*tr+1` trace positions is likewise lambda-independent, so the total
captured fraction is a single matrix-vector product ``mtf @ w`` against a
precomputable `(Nu,)` weight vector. Both closed forms are exact
reassociations of the same sums (identical up to float rounding, ~1e-15).

For very large `L` (e.g. an all-sky-lines call with `L` ~ a few thousand),
callers can precompute `spectro_mtf` in `L`-chunks (see `constants.L_CHUNK`)
and pass each chunk's MTF array in via `spectro_dist`'s/`frac_trace`'s `mtf`
argument, rather than materializing `(L, Nu)` for the whole line list at
once; this module does not chunk internally.
"""

from __future__ import annotations

import numpy as np
from scipy.special import j0, j1

from .config import Spectrograph, field_interp
from .constants import ARCSEC_PER_URAD, RAT_HL_SL_EXP
from .materials import si_abslength, si_index_real

# --- gsSpectroMTF / gsSpectroDist quadrature grid (gsetc.c:681) ------------
_MTF_NU = 1000
_MTF_DU = 0.005
#: Fixed `u`-grid (cycles/pixel) used by `spectro_dist`/`frac_trace` when the
#: caller does not supply a precomputed `mtf`; midpoint rule, `Nu=1000`,
#: `du=0.005` (gsetc.c:681).
U_GRID = (np.arange(_MTF_NU) + 0.5) * _MTF_DU

# --- gsGeometricThroughput quadrature grid (gsetc.c:513-514) ---------------
_EE_NU = 50
_EE_DU = 0.024

# --- Si defocus model (gsetc.c:640-661): fixed 3-slab / 6-depth-point model
_SI_DEFOCUS_N = 3


def _si_defocus_d0(mfp: np.ndarray, ddepth: float) -> np.ndarray:
    """First depth sample `d0` for the Si defocus model (gsetc.c:654).

    ``d0 = mfp<1e4*ddepth ? mfp*(1-(1+ddepth/mfp)*exp(-ddepth/mfp))/(1-exp(-ddepth/mfp))
    : 0.5*ddepth``. Computed with a safe denominator substitute so that the
    (unused, per the C ternary) small-`ddepth/mfp` branch never triggers a
    0/0 warning when `mfp` is very large.
    """
    ratio = ddepth / mfp
    denom = 1.0 - np.exp(-ratio)
    safe_denom = np.where(denom > 0.0, denom, 1.0)
    small_ratio_branch = mfp * (1.0 - (1.0 + ratio) * np.exp(-ratio)) / safe_denom
    return np.where(mfp < 1e4 * ddepth, small_ratio_branch, 0.5 * ddepth)


def spectro_mtf(spectro: Spectrograph, i_arm: int, lam, u) -> np.ndarray:
    """Spectrograph PSF MTF (real part), `gsSpectroMTF` (gsetc.c:621-667).

    Parameters
    ----------
    spectro : Spectrograph
    i_arm : int
        Internal 0-based arm index.
    lam : float or array_like, shape (L,)
        Wavelength in nm.
    u : float or array_like, shape (Nu,)
        Wavenumber in cycles/pixel.

    Returns
    -------
    ndarray, shape (L, Nu)
        MTF at each (wavelength, wavenumber) pair.
    """
    lam = np.atleast_1d(np.asarray(lam, dtype=np.float64))
    u = np.atleast_1d(np.asarray(u, dtype=np.float64))

    # Fiber size (u-only). `np.where` evaluates both branches elementwise,
    # so guard the denominator at u=0 to avoid a spurious 0/0 warning (the
    # divide-branch result there is discarded by `np.where` regardless).
    d_spot = spectro.diam[i_arm] / spectro.pix[i_arm]
    x = np.pi * d_spot * u
    x_safe = np.where(x == 0.0, 1.0, x)
    fiber = np.where(np.abs(d_spot * u) > 1e-6, 2.0 * j1(x_safe) / x_safe, 1.0)

    # Pixelization (u-only).
    pu = np.pi * u
    pu_safe = np.where(pu == 0.0, 1.0, pu)
    pixel = np.where(np.abs(u) > 1e-9, np.sin(pu_safe) / pu_safe, 1.0)

    # Camera spot size, Gaussian (u-only).
    sigma_cam = spectro.rms_cam[i_arm] / spectro.pix[i_arm]
    gauss = np.exp(-2.0 * np.pi**2 * sigma_cam**2 * u**2)

    mtf = (fiber * pixel * gauss)[None, :]  # (1, Nu) -> broadcasts to (L, Nu)

    # Defocus in the detector -- Si only (depends on both lambda and u).
    if spectro.Dtype[i_arm] == 0:
        n_depth = _SI_DEFOCUS_N
        ddepth = spectro.thick[i_arm] / n_depth
        thick = spectro.thick[i_arm]
        fratio = spectro.fratio[i_arm]
        pix = spectro.pix[i_arm]

        mfp = np.atleast_1d(si_abslength(lam, spectro.temperature[i_arm]))  # (L,)
        nsi = np.atleast_1d(si_index_real(lam))  # (L,)
        d0 = _si_defocus_d0(mfp, ddepth)  # (L,)

        numer = np.zeros((lam.size, u.size))
        denom = np.zeros(lam.size)
        for i in range(2 * n_depth):
            depth = d0 + ddepth * i  # (L,)
            contrib = np.exp(-ddepth * i / mfp) * np.where(depth > thick, 0.3, 1.0)
            rspot = depth / nsi / 2.0 / fratio / pix  # (L,)
            arg = 2.0 * np.pi * rspot[:, None] * u[None, :]  # (L, Nu)
            denom += contrib
            numer += contrib[:, None] * (
                0.1666666667 * np.cos(0.866025404 * arg)
                + 0.5 * np.cos(0.5 * arg)
                + 0.3333333333
            )
        mtf = mtf * (numer / denom[:, None])

    # Scattering from grating (depends on both lambda and u).
    grating = np.exp(
        -lam[:, None] / spectro.dl[i_arm] / spectro.nline[i_arm] * np.abs(u)[None, :]
    )
    mtf = mtf * grating

    return mtf


def spectro_dist(
    spectro: Spectrograph,
    i_arm: int,
    lam,
    pos,
    sigma,
    N: int,
    mtf: np.ndarray | None = None,
) -> np.ndarray:
    """Per-pixel flux fraction of a spectral feature, `gsSpectroDist`
    (gsetc.c:673-690).

    Parameters
    ----------
    spectro : Spectrograph
    i_arm : int
        Internal 0-based arm index.
    lam : float or array_like, shape (L,)
        Wavelength in nm (one entry per independent feature/line).
    pos : float or array_like, shape (L,) or scalar
        Feature center, in pixels.
    sigma : float or array_like, shape (L,) or scalar
        Extra Gaussian smearing width, in pixels.
    N : int
        Window width, in pixels; output covers pixels `0 .. N-1`.
    mtf : ndarray, shape (L, Nu), optional
        Precomputed `spectro_mtf(spectro, i_arm, lam, U_GRID)`; computed
        internally if omitted. Callers looping over the same `(spectro,
        i_arm, lam)` with varying `pos`/`sigma` (e.g. `frac_trace`) should
        pass this in to avoid recomputing the MTF each time.

    Returns
    -------
    ndarray, shape (L, N)
        Fraction of flux landing in each of the `N` pixels.

    Notes
    -----
    The C code's `fr[ip] = sum_u 2*du*cos(2*pi*u*(pos-ip))*mtf1d` is
    rewritten via the cosine angle-addition identity as two `(L, Nu) @
    (Nu, N)` matrix products, so the underlying `(L, Nu, N)` triple-loop
    tensor is never materialized. When `pos` and `sigma` are both scalars,
    every lambda-independent factor folds into one `(Nu, N)` weight matrix
    and the result is the single product `mtf @ W` -- an exact
    reassociation of the same sum (see module docstring).
    """
    lam = np.atleast_1d(np.asarray(lam, dtype=np.float64))
    L = lam.size
    scalar_pos = np.ndim(pos) == 0
    scalar_sigma = np.ndim(sigma) == 0

    u = U_GRID
    if mtf is None:
        mtf = spectro_mtf(spectro, i_arm, lam, u)
    else:
        assert mtf.shape == (L, u.size), (
            f"spectro_dist: precomputed mtf has shape {mtf.shape}, expected "
            f"{(L, u.size)} (L=lam.size, Nu=U_GRID.size)"
        )

    ip = np.arange(N, dtype=np.float64)
    cm = np.cos(2.0 * np.pi * np.outer(u, ip))  # (Nu, N), shared across lines
    sm = np.sin(2.0 * np.pi * np.outer(u, ip))  # (Nu, N)

    if scalar_pos and scalar_sigma:
        # Fast path (see module docstring): one shared (pos, sigma) for all
        # lambda -> fold the trig/smearing terms into a lambda-independent
        # (Nu, N) weight matrix and collapse to a single matmul.
        ang_u = 2.0 * np.pi * u * float(pos)  # (Nu,)
        w = np.cos(ang_u)[:, None] * cm + np.sin(ang_u)[:, None] * sm  # (Nu, N)
        if float(sigma) != 0.0:
            w *= np.exp(-2.0 * np.pi**2 * float(sigma) ** 2 * u**2)[:, None]
        return 2.0 * _MTF_DU * (mtf @ w)  # (L, N)

    pos = np.broadcast_to(np.asarray(pos, dtype=np.float64), (L,))
    if scalar_sigma:
        # Per-u smearing factor is shared by every lambda: apply it as a
        # (Nu,)-row scale instead of a full (L, Nu) exp evaluation.
        if float(sigma) != 0.0:
            M = mtf * np.exp(-2.0 * np.pi**2 * float(sigma) ** 2 * u**2)[None, :]
        else:
            M = mtf
    else:
        sigma = np.broadcast_to(np.asarray(sigma, dtype=np.float64), (L,))
        M = mtf * np.exp(
            -2.0 * np.pi**2 * sigma[:, None] ** 2 * u[None, :] ** 2
        )  # (L,Nu)

    ang = 2.0 * np.pi * u[None, :] * pos[:, None]  # (L, Nu)
    fr = 2.0 * _MTF_DU * ((M * np.cos(ang)) @ cm + (M * np.sin(ang)) @ sm)  # (L, N)
    return fr


def frac_trace(
    spectro: Spectrograph,
    i_arm: int,
    lam,
    tr: int,
    mtf: np.ndarray | None = None,
) -> np.ndarray:
    """Fraction of flux captured in the extraction window, `gsFracTrace`
    (gsetc.c:696-710).

    Sums `spectro_dist` over the `spectro.width[i_arm]`-pixel window,
    centered at `0.5*(N-1)`, plus (if `tr>=1`) the `tr` adjacent traces on
    each side, offset by `j*spectro.sep[i_arm]/spectro.pix[i_arm]` pixels
    for `j` in `-tr .. tr`.

    Parameters
    ----------
    spectro : Spectrograph
    i_arm : int
        Internal 0-based arm index.
    lam : float or array_like, shape (L,)
        Wavelength in nm.
    tr : int
        Number of adjacent traces (each side) to include; `tr=0` is the
        correct trace only.
    mtf : ndarray, shape (L, Nu), optional
        Precomputed `spectro_mtf(spectro, i_arm, lam, U_GRID)`; computed
        once internally (and reused across the `2*tr+1` `spectro_dist`
        calls) if omitted.

    Returns
    -------
    ndarray, shape (L,)
        Fraction of flux captured, summed over the window and all traces.

    Notes
    -----
    The window/trace double sum is lambda-independent (each trace's `pos`
    is one scalar for every wavelength), so it folds into a single `(Nu,)`
    weight vector ``w(u) = 2*du * sum_j [cos(2*pi*u*pos_j)*sum_ip
    cos(2*pi*u*ip) + sin(2*pi*u*pos_j)*sum_ip sin(2*pi*u*ip)]`` and the
    total captured fraction is the single matrix-vector product `mtf @ w`
    -- an exact reassociation of the per-trace `spectro_dist` window sums
    (see module docstring).
    """
    lam = np.atleast_1d(np.asarray(lam, dtype=np.float64))
    N = int(spectro.width[i_arm])
    if mtf is None:
        mtf = spectro_mtf(spectro, i_arm, lam, U_GRID)
    else:
        assert mtf.shape == (lam.size, U_GRID.size), (
            f"frac_trace: precomputed mtf has shape {mtf.shape}, expected "
            f"{(lam.size, U_GRID.size)} (L=lam.size, Nu=U_GRID.size)"
        )

    u = U_GRID
    ip = np.arange(N, dtype=np.float64)
    cm_sum = np.cos(2.0 * np.pi * np.outer(u, ip)).sum(axis=1)  # (Nu,)
    sm_sum = np.sin(2.0 * np.pi * np.outer(u, ip)).sum(axis=1)  # (Nu,)

    sep_pix = spectro.sep[i_arm] / spectro.pix[i_arm]
    w = np.zeros(u.size)
    for j in range(-tr, tr + 1):
        pos = 0.5 * (N - 1) + j * sep_pix
        ang_u = 2.0 * np.pi * u * pos
        w += np.cos(ang_u) * cm_sum + np.sin(ang_u) * sm_sum
    return 2.0 * _MTF_DU * (mtf @ w)


def geometric_throughput(
    spectro: Spectrograph,
    lam,
    r_eff: float,
    decent: float,
    fieldang: float,
    seeing_fwhm_800: float,
):
    """Encircled energy in the fiber, `gsGeometricThroughput` (gsetc.c:504-570).

    Combines the Kolmogorov seeing MTF (1/e at `u=uscale`), a leading-order
    telescope-diffraction term, decenter (via `j0`), and an exponential
    surface-brightness profile (scale radius `r_eff / RAT_HL_SL_EXP`)
    convolved with the circular fiber aperture, via a 50-point midpoint
    Hankel-transform quadrature (`u = (arange(50)+0.5)*0.024`,
    gsetc.c:513-514).

    The C signature also takes an unused `flags` argument (gsetc.c:505),
    which is dropped here.

    Parameters
    ----------
    spectro : Spectrograph
    lam : float or array_like, shape (L,)
        Wavelength in nm.
    r_eff : float
        Target effective radius, arcsec (0 for a point source).
    decent : float
        Fiber decenter from the target, arcsec.
    fieldang : float
        Field angle, degrees (0 on axis).
    seeing_fwhm_800 : float
        Seeing FWHM at 800nm, arcsec.

    Returns
    -------
    float or ndarray
        Encircled energy (dimensionless fraction), same shape as `lam`.
    """
    lam_arr = np.asarray(lam, dtype=np.float64)
    scalar_input = lam_arr.ndim == 0
    lam_arr = np.atleast_1d(lam_arr)

    sigma = field_interp(spectro.rms_spot, fieldang, spectro.rfov)
    efl = field_interp(spectro.EFL, fieldang, spectro.rfov)
    sigma *= ARCSEC_PER_URAD / efl  # microns -> arcsec

    R = spectro.fiber_ent_rad / efl * ARCSEC_PER_URAD  # fiber radius, arcsec

    u = (np.arange(_EE_NU) + 0.5) * _EE_DU  # (Nu,)
    k = 2.0 * np.pi * u  # (Nu,)

    theta_d = (
        0.001 * lam_arr / spectro.D_outer / (1.0 - spectro.centobs) * ARCSEC_PER_URAD
    )  # (L,)
    uscale = 0.465 / seeing_fwhm_800 * (lam_arr / 800.0) ** 0.2  # (L,)
    rs = r_eff / RAT_HL_SL_EXP  # scalar, galaxy scale length

    # Telescope PSF (leading-order diffraction).
    seeing_term = np.exp(-(k**2) * sigma**2 / 2.0)  # (Nu,)
    decent_term = j0(k * decent)  # (Nu,)
    kolmogorov_term = np.exp(-((u[None, :] / uscale[:, None]) ** (5.0 / 3.0)))  # (L,Nu)
    diffraction_term = np.exp(-4.0 / np.pi * u[None, :] * theta_d[:, None])  # (L,Nu)
    g_tilde = (
        seeing_term[None, :] * decent_term[None, :] * kolmogorov_term * diffraction_term
    )

    # Galaxy profile.
    f_tilde = (1.0 + (k * rs) ** 2) ** -1.5  # (Nu,)

    ee_terms = (
        R
        * 2.0
        * np.pi
        * _EE_DU
        * j1(2.0 * np.pi * u * R)[None, :]
        * f_tilde[None, :]
        * g_tilde
    )
    ee = ee_terms.sum(axis=1)
    return float(ee[0]) if scalar_input else ee


def effective_area(spectro: Spectrograph, i_arm: int, lam, fieldang: float):
    """Effective collecting area, `gsAeff` (gsetc.c:576-615).

    Geometric area (with central obscuration) x field-angle-dependent
    vignetting x wavelength-dependent system throughput (piecewise-linear
    interpolation of the arm's throughput grid, clamped to the endpoint
    values outside its range).

    Parameters
    ----------
    spectro : Spectrograph
    i_arm : int
        Internal 0-based arm index.
    lam : float or array_like
        Wavelength in nm.
    fieldang : float
        Field angle, degrees (0 on axis).

    Returns
    -------
    float or ndarray
        Effective area, m^2, same shape as `lam`.
    """
    lam_arr = np.asarray(lam, dtype=np.float64)
    scalar_input = lam_arr.ndim == 0
    lam_arr = np.atleast_1d(lam_arr)

    aeff = np.pi / 4.0 * spectro.D_outer**2 * (1.0 - spectro.centobs**2)

    vig = field_interp(spectro.vignette, fieldang, spectro.rfov)
    aeff = aeff * vig

    imin = int(spectro.istart[i_arm])
    imax = int(spectro.istart[i_arm + 1])
    # np.interp clamps to fp[0]/fp[-1] outside [xp[0], xp[-1]], matching the
    # C code's explicit <=/>= endpoint branches (gsetc.c:598-611).
    thr = np.interp(lam_arr, spectro.l[imin:imax], spectro.T[imin:imax])
    aeff = aeff * thr

    return float(aeff[0]) if scalar_input else aeff
