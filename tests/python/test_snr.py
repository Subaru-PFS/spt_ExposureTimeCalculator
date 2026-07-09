"""Tests for pfsspecsim.etc.snr (task T9).

Ports of gsGetSignal (gsetc.c:1071-1118), gsGetSNR (gsetc.c:1130-1163),
gsGetSNR_Single (gsetc.c:1174-1257), gsGetSNR_OII (gsetc.c:1273-1360), and
gsGetSNR_Continuum (gsetc.c:1367-1446).

In addition to physical-sanity checks (optimal >= uniform SNR, SNR
proportional to line flux, an essentially-invisible-magnitude source having
~zero continuum SNR), this file includes literal, scalar, loop-based
Python transcriptions of gsGetSignal and gsGetSNR_OII's `snrType == 2`
branch, and cross-checks the vectorized implementation against them to
~1e-10 -- the same "independent oracle" pattern used in test_psf.py and
test_noise.py, which caught real transcription errors in earlier tasks.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from pfsspecsim.etc import psf, snr
from pfsspecsim.etc.atmosphere import transmission as atm_transmission
from pfsspecsim.etc.config import find_config_file, load_spectrograph_config
from pfsspecsim.etc.constants import (
    C_KM_PER_S,
    C_NM_PER_S,
    NP_WIN,
    OII_LAMBDA,
    PHOTONS_PER_ERG_1NM,
)
from pfsspecsim.etc.extinction import alambda_over_ebv
from pfsspecsim.etc.noise import compute_noise_arm
from pfsspecsim.etc.params import EtcParams

# Protected fixture; never modify (see task brief).
LEGACY_CONFIG_FIXTURE = Path(__file__).resolve().parents[1] / "PFS.20211220.dat"

ARM_BLUE = 0  # Si (Dtype=0)
ARM_NIR = 2  # HgCdTe (Dtype=1)


@pytest.fixture(scope="module")
def spectro():
    # The "real" packaged config named in the task brief.
    return load_spectrograph_config(find_config_file("20240714"), degrade=1.0)


@pytest.fixture(scope="module")
def params():
    return EtcParams(
        seeing=0.8,
        zenith_ang=45.0,
        galactic_ext=0.0,
        field_ang=0.675,
        fiber_offset=0.0,
        moon_zenith_ang=135.0,  # Moon below horizon: disabled.
        exp_time=900.0,
        exp_num=8,
        mag=22.5,
        sky_type="11005",
        degrade=1.0,
        obsc_fov_dep=False,
        hgcdte_sutr=True,
    )


@pytest.fixture(scope="module")
def arm_noise(params, spectro):
    return {ia: compute_noise_arm(params, spectro, ia) for ia in range(spectro.N_arms)}


# --- Scalar (loop-based) transcription of gsGetSignal / gsGetSNR_OII ------
#
# Deliberately re-derives the C assembly logic (the iref window arithmetic,
# the 41-point trans averages, the counts formula, the full-Npix myNoise
# array) from scratch in scalar Python, while still calling
# psf.geometric_throughput/frac_trace/effective_area/spectro_dist and
# atmosphere.transmission/extinction.alambda_over_ebv directly -- those are
# already independently unit-tested (T3/T5/T6) building blocks, not the
# thing under test here.


def _scalar_get_signal(spectro_, params_, i_arm, lam, F, sigma_v, r_eff):
    """Scalar transcription of gsGetSignal (gsetc.c:1071-1118).

    Returns the *full* `(Npix,)` Signal list (mirroring the C array) and
    the window's `iref`, or `iref=None` if the C function's early-return
    (gsetc.c:1092-1093) fires (Signal all zero).
    """
    npix = int(spectro_.npix[i_arm])
    lmin = float(spectro_.lmin[i_arm])
    dl = float(spectro_.dl[i_arm])
    signal = [0.0] * npix

    pos = (lam - lmin) / dl
    iref = math.floor(pos - NP_WIN / 2.0)
    if iref < -NP_WIN or iref >= npix:
        return signal, None
    if iref < 0:
        iref = 0
    if iref > npix - NP_WIN:
        iref = npix - NP_WIN

    trans = 0.0
    den = 0.0
    x = -4.0
    while x < 4.01:
        trans += atm_transmission(
            lam * (1.0 + x * sigma_v / C_KM_PER_S), params_.zenith_ang, params_.sky_type
        ) * math.exp(-0.5 * x * x)
        den += math.exp(-0.5 * x * x)
        x += 0.2
    trans /= den

    counts = (
        F
        * trans
        * 10.0 ** (-0.4 * float(alambda_over_ebv(lam)) * params_.galactic_ext)
        * psf.geometric_throughput(
            spectro_,
            lam,
            r_eff,
            params_.fiber_offset,
            params_.field_ang,
            params_.seeing,
        )
        * float(psf.frac_trace(spectro_, i_arm, np.array([lam]), tr=0)[0])
        * PHOTONS_PER_ERG_1NM
        * lam
        * params_.exp_time
        * float(psf.effective_area(spectro_, i_arm, lam, params_.field_ang))
        * 1e4
    )

    sigma_pix = sigma_v / C_KM_PER_S * lam / dl
    fr = psf.spectro_dist(
        spectro_, i_arm, np.array([lam]), pos - iref, sigma_pix, NP_WIN
    )[0]
    for ip in range(NP_WIN):
        signal[iref + ip] = fr[ip] * counts
    return signal, iref


def _scalar_snr_oii_type2(
    spectro_, params_, i_arm, z, F, sigma_v, r_eff, src_cont, roii, noise
):
    """Scalar transcription of gsGetSNR_OII's `snrType == 2` branch
    (gsetc.c:1273-1360, combining formula at 1345-1356)."""
    lam0 = OII_LAMBDA[0] * (1.0 + z)
    lam1 = OII_LAMBDA[1] * (1.0 + z)

    roii = min(max(roii, 0.667), 3.87)
    frac0 = roii / (1.0 + roii)
    frac1 = 1.0 / (1.0 + roii)

    ll = (lam0 + lam1) / 2.0
    trans = 0.0
    den = 0.0
    x = -4.0
    while x < 4.01:
        trans += atm_transmission(
            lam0 + (0.5 + 0.5 * x) * (lam1 - lam0), params_.zenith_ang, params_.sky_type
        ) * math.exp(-0.5 * x * x)
        den += math.exp(-0.5 * x * x)
        x += 0.2
    trans /= den

    dl = float(spectro_.dl[i_arm])
    counts = (
        src_cont
        * trans
        * 10.0 ** (-0.4 * float(alambda_over_ebv(ll)) * params_.galactic_ext)
        * psf.geometric_throughput(
            spectro_, ll, r_eff, params_.fiber_offset, params_.field_ang, params_.seeing
        )
        * float(psf.frac_trace(spectro_, i_arm, np.array([ll]), tr=0)[0])
        * PHOTONS_PER_ERG_1NM
        * ll
        * params_.exp_time
        * float(psf.effective_area(spectro_, i_arm, ll, params_.field_ang))
        * 1e4
    )
    if params_.hgcdte_sutr and int(spectro_.Dtype[i_arm]) == 1:
        counts *= 1.2
    counts *= C_NM_PER_S * dl / (ll * ll)

    my_noise = [float(n) + counts for n in noise]

    signal0, _iref0 = _scalar_get_signal(
        spectro_, params_, i_arm, lam0, frac0 * F, sigma_v, r_eff
    )
    signal1, _iref1 = _scalar_get_signal(
        spectro_, params_, i_arm, lam1, frac1 * F, sigma_v, r_eff
    )

    snr_sq = 0.0
    for ipix in range(len(my_noise)):
        s = signal0[ipix] + signal1[ipix]
        snr_sq += s * s / my_noise[ipix]
    return math.sqrt(snr_sq)


class TestGetSignalScalarTranscription:
    @pytest.mark.parametrize(
        "i_arm,lam,sigma_v,r_eff",
        [
            (ARM_BLUE, 500.0, 70.0, 0.3),
            (ARM_BLUE, 559.0, 150.0, 0.0),
            (ARM_NIR, 1000.0, 70.0, 0.5),
        ],
    )
    def test_matches_scalar_reference(
        self, spectro, params, i_arm, lam, sigma_v, r_eff
    ):
        F = 1e-16
        signal_window, iref = snr.get_signal(
            params, spectro, i_arm, lam, F, sigma_v, r_eff
        )
        assert signal_window.shape == (1, NP_WIN)

        ref_signal, ref_iref = _scalar_get_signal(
            spectro, params, i_arm, lam, F, sigma_v, r_eff
        )
        assert ref_iref is not None
        assert int(iref[0]) == ref_iref
        np.testing.assert_allclose(
            signal_window[0],
            ref_signal[ref_iref : ref_iref + NP_WIN],
            rtol=1e-10,
            atol=1e-30,
        )

    def test_out_of_range_gives_all_zero_signal(self, spectro, params):
        # Far outside [-NP_WIN, Npix) in pixel units for arm 0.
        lmin = spectro.lmin[ARM_BLUE]
        dl = spectro.dl[ARM_BLUE]
        npix = int(spectro.npix[ARM_BLUE])
        lam_far = lmin + dl * (npix + 1000)
        signal_window, iref = snr.get_signal(
            params, spectro, ARM_BLUE, lam_far, 1e-16, 70.0, 0.3
        )
        np.testing.assert_array_equal(signal_window, 0.0)
        # iref is still clipped into a valid range for safe indexing.
        assert 0 <= int(iref[0]) <= npix - NP_WIN

        ref_signal, ref_iref = _scalar_get_signal(
            spectro, params, ARM_BLUE, lam_far, 1e-16, 70.0, 0.3
        )
        assert ref_iref is None
        assert all(v == 0.0 for v in ref_signal)

    def test_vectorized_over_multiple_lines(self, spectro, params):
        lam = np.array([500.0, 550.0, 600.0])
        signal_window, iref = snr.get_signal(
            params, spectro, ARM_BLUE, lam, 1e-16, 70.0, 0.3
        )
        assert signal_window.shape == (3, NP_WIN)
        assert iref.shape == (3,)
        for i, lam_i in enumerate(lam):
            scalar_window, scalar_iref = snr.get_signal(
                params, spectro, ARM_BLUE, lam_i, 1e-16, 70.0, 0.3
            )
            assert int(iref[i]) == int(scalar_iref[0])
            np.testing.assert_allclose(signal_window[i], scalar_window[0], rtol=1e-12)


class TestSnrLine:
    def test_optimal_at_least_uniform(self, spectro, params, arm_noise):
        lam = 559.0
        opt = snr.snr_line(
            params, spectro, ARM_BLUE, lam, 1e-16, 70.0, 0.3, arm_noise[0].noise, 0
        )
        uni = snr.snr_line(
            params, spectro, ARM_BLUE, lam, 1e-16, 70.0, 0.3, arm_noise[0].noise, 1
        )
        assert opt >= uni > 0.0

    def test_snr_proportional_to_flux(self, spectro, params, arm_noise):
        lam = 559.0
        snr_1 = snr.snr_line(
            params, spectro, ARM_BLUE, lam, 1e-16, 70.0, 0.3, arm_noise[0].noise, 0
        )
        snr_2 = snr.snr_line(
            params, spectro, ARM_BLUE, lam, 2e-16, 70.0, 0.3, arm_noise[0].noise, 0
        )
        assert snr_2 == pytest.approx(2.0 * snr_1, rel=1e-10)

    def test_scalar_input_returns_scalar_float(self, spectro, params, arm_noise):
        result = snr.snr_line(
            params, spectro, ARM_BLUE, 559.0, 1e-16, 70.0, 0.3, arm_noise[0].noise
        )
        assert isinstance(result, float)

    def test_vector_input_returns_array(self, spectro, params, arm_noise):
        result = snr.snr_line(
            params,
            spectro,
            ARM_BLUE,
            [500.0, 559.0],
            1e-16,
            70.0,
            0.3,
            arm_noise[0].noise,
        )
        assert isinstance(result, np.ndarray)
        assert result.shape == (2,)


class TestSnrSingle:
    def test_quirk_transmission_independent_of_sigma_v(
        self, spectro, params, arm_noise, monkeypatch
    ):
        # gsetc.c:1213-1218 QUIRK: the object-continuum transmission
        # average is evaluated at the fixed line wavelength for every
        # quadrature point, so it must be identical for any sigma_v --
        # unlike get_signal's velocity-smeared trans, which does vary with
        # sigma_v. Guards against a future "fix" silently reintroducing
        # velocity smearing here.
        #
        # Spy on snr.atm_transmission (the exact name snr_single calls) to
        # capture what it's actually invoked with/returns for each sigma_v,
        # rather than comparing a fresh atm_transmission(lam, ...) call
        # against itself (which would hold regardless of what snr_single
        # does internally). snr_single's own get_signal call also goes
        # through this same module-level name for its (genuinely
        # sigma_v-dependent) 41-point velocity-smeared quadrature, passing a
        # (Nz, 41)-shaped grid -- filter that out by keeping only the
        # single-point (Nz==1 line, one wavelength each) calls, which is
        # what snr_single's own continuum-transmission line does.
        captured: list[np.ndarray] = []
        real_atm_transmission = snr.atm_transmission

        def _spy(lam_arr, zenith_ang, sky_type):
            result = real_atm_transmission(lam_arr, zenith_ang, sky_type)
            arr = np.atleast_1d(np.asarray(lam_arr, dtype=np.float64))
            if arr.shape == (1,):
                captured.append(np.array(result, dtype=np.float64, copy=True))
            return result

        monkeypatch.setattr(snr, "atm_transmission", _spy)

        lam = 559.0
        snr_narrow = snr.snr_single(
            params, spectro, ARM_BLUE, 22.5, lam, 1e-16, 10.0, 0.3, arm_noise[0].noise
        )
        snr_wide = snr.snr_single(
            params, spectro, ARM_BLUE, 22.5, lam, 1e-16, 300.0, 0.3, arm_noise[0].noise
        )

        assert len(captured) == 2
        # The continuum transmission snr_single actually used is identical
        # across sigma_v (the quirk) ...
        np.testing.assert_array_equal(captured[0], captured[1])
        # ... and equals the plain pointwise transmission at the fixed line
        # wavelength, not some sigma_v-smeared average.
        expected = real_atm_transmission(
            np.asarray([lam]), params.zenith_ang, params.sky_type
        )
        np.testing.assert_allclose(captured[0], expected)

        # sigma_v still affects the *line* signal (get_signal), so the two
        # SNRs need not be equal.
        assert snr_narrow != snr_wide  # sanity: sigma_v does change something

    def test_fainter_host_increases_line_snr(self, spectro, params, arm_noise):
        # A brighter host continuum adds more shot noise to the line
        # measurement, so a fainter host (larger mag) must give an equal or
        # higher line SNR at fixed line flux.
        lam = 559.0
        snr_bright_host = snr.snr_single(
            params, spectro, ARM_BLUE, 18.0, lam, 1e-16, 70.0, 0.3, arm_noise[0].noise
        )
        snr_faint_host = snr.snr_single(
            params, spectro, ARM_BLUE, 26.0, lam, 1e-16, 70.0, 0.3, arm_noise[0].noise
        )
        assert snr_faint_host >= snr_bright_host

    def test_mag_array_matches_scalar_per_line(self, spectro, params, arm_noise):
        lam = np.array([500.0, 559.0])
        mag = np.array([20.0, 24.0])
        vec = snr.snr_single(
            params, spectro, ARM_BLUE, mag, lam, 1e-16, 70.0, 0.3, arm_noise[0].noise
        )
        for i in range(2):
            scalar = snr.snr_single(
                params,
                spectro,
                ARM_BLUE,
                mag[i],
                lam[i],
                1e-16,
                70.0,
                0.3,
                arm_noise[0].noise,
            )
            assert vec[i] == pytest.approx(scalar, rel=1e-12)


class TestSnrOii:
    def test_roii_clamped(self, spectro, params, arm_noise):
        z = 0.5
        below = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            z,
            1e-16,
            70.0,
            0.3,
            0.0,
            0.1,
            arm_noise[0].noise,
            2,
        )
        at_floor = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            z,
            1e-16,
            70.0,
            0.3,
            0.0,
            0.667,
            arm_noise[0].noise,
            2,
        )
        assert below == pytest.approx(at_floor, rel=1e-12)

        above = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            z,
            1e-16,
            70.0,
            0.3,
            0.0,
            10.0,
            arm_noise[0].noise,
            2,
        )
        at_ceiling = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            z,
            1e-16,
            70.0,
            0.3,
            0.0,
            3.87,
            arm_noise[0].noise,
            2,
        )
        assert above == pytest.approx(at_ceiling, rel=1e-12)

    def test_optimal_at_least_uniform(self, spectro, params, arm_noise):
        z = 0.5
        opt = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            z,
            1e-16,
            70.0,
            0.3,
            0.0,
            1.0,
            arm_noise[0].noise,
            0,
        )
        uni = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            z,
            1e-16,
            70.0,
            0.3,
            0.0,
            1.0,
            arm_noise[0].noise,
            1,
        )
        assert opt >= uni > 0.0

    def test_snr_proportional_to_flux(self, spectro, params, arm_noise):
        z = 0.5
        snr_1 = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            z,
            1e-16,
            70.0,
            0.3,
            0.0,
            1.0,
            arm_noise[0].noise,
            2,
        )
        snr_2 = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            z,
            2e-16,
            70.0,
            0.3,
            0.0,
            1.0,
            arm_noise[0].noise,
            2,
        )
        assert snr_2 == pytest.approx(2.0 * snr_1, rel=1e-10)

    def test_out_of_arm_range_gives_zero(self, spectro, params, arm_noise):
        # z=1.0 pushes both doublet lines (372.71/372.98nm rest-frame) to
        # ~745nm, well redward of arm 0's [380, 650]nm range, but still a
        # perfectly ordinary wavelength for the atmosphere/extinction
        # models -- unlike a much larger z, which would push the doublet
        # past those models' valid domain and isn't a case any real caller
        # would hit (the arm-gate check, gsetc.c:2036, restricts z per arm
        # before ever calling gsGetSNR_OII; that gating is the T10 engine's
        # responsibility, not this function's -- see module docstring).
        result = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            1.0,
            1e-16,
            70.0,
            0.3,
            0.0,
            1.0,
            arm_noise[0].noise,
            2,
        )
        assert result == pytest.approx(0.0, abs=1e-8)

    def test_vectorized_matches_scalar_per_z(self, spectro, params, arm_noise):
        z = np.array([0.3, 0.8, 1.6])
        vec = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            z,
            1e-16,
            70.0,
            0.3,
            1e-29,
            1.3,
            arm_noise[0].noise,
            2,
        )
        for i in range(z.size):
            scalar = snr.snr_oii(
                params,
                spectro,
                ARM_BLUE,
                z[i],
                1e-16,
                70.0,
                0.3,
                1e-29,
                1.3,
                arm_noise[0].noise,
                2,
            )
            assert vec[i] == pytest.approx(scalar, rel=1e-10)

    @pytest.mark.parametrize(
        "z,src_cont,roii", [(0.3, 0.0, 1.0), (0.9, 1e-29, 0.8), (1.7, 5e-30, 2.5)]
    )
    def test_matches_scalar_transcription(
        self, spectro, params, arm_noise, z, src_cont, roii
    ):
        vec = snr.snr_oii(
            params,
            spectro,
            ARM_BLUE,
            z,
            1e-16,
            70.0,
            0.3,
            src_cont,
            roii,
            arm_noise[0].noise,
            2,
        )
        ref = _scalar_snr_oii_type2(
            spectro,
            params,
            ARM_BLUE,
            z,
            1e-16,
            70.0,
            0.3,
            src_cont,
            roii,
            arm_noise[0].noise,
        )
        assert vec == pytest.approx(ref, rel=1e-10)


class TestSnrContinuum:
    def test_shape(self, spectro, params, arm_noise):
        result = snr.snr_continuum(
            params, spectro, ARM_BLUE, 22.5, 0.3, arm_noise[0].noise
        )
        npix = int(spectro.npix[ARM_BLUE])
        for arr in (
            result.snr,
            result.counts,
            result.noise,
            result.mag,
            result.trans,
            result.sample_factor,
        ):
            assert arr.shape == (npix,)

    def test_mag_99_9_gives_zero_snr(self, spectro, params, arm_noise):
        result = snr.snr_continuum(
            params, spectro, ARM_BLUE, 99.9, 0.3, arm_noise[0].noise
        )
        assert np.all(result.snr < 1e-15)

    def test_snr_increases_with_brightness(self, spectro, params, arm_noise):
        faint = snr.snr_continuum(
            params, spectro, ARM_BLUE, 24.0, 0.3, arm_noise[0].noise
        )
        bright = snr.snr_continuum(
            params, spectro, ARM_BLUE, 20.0, 0.3, arm_noise[0].noise
        )
        assert np.all(bright.snr >= faint.snr)

    def test_left_edge_wavelength_quirk(self, spectro, params, arm_noise):
        # gsGetSNR_Continuum evaluates each pixel's wavelength at the left
        # edge (lmin + dl*ipix), not the center used by
        # noise.compute_noise_arm's sky continuum grid -- verified here by
        # checking pixel 0's dimensionless throughput curve (`.trans`,
        # `counts/src_cont`, so independent of the chosen `mag`) against a
        # direct evaluation at `lam=lmin` (the left edge).
        result = snr.snr_continuum(
            params, spectro, ARM_BLUE, 22.5, 0.3, arm_noise[0].noise
        )
        lmin = spectro.lmin[ARM_BLUE]
        dl = spectro.dl[ARM_BLUE]

        # Use noise.smoothed_transmission (the same PSF-smoothed average
        # snr_continuum itself uses) rather than the raw pointwise
        # atmosphere.transmission, so this isolates the left-edge-vs-center
        # wavelength convention rather than conflating it with PSF-smoothing
        # differences.
        from pfsspecsim.etc.noise import smoothed_transmission

        expected_trans_pix0 = snr._continuum_counts_at(
            params,
            spectro,
            ARM_BLUE,
            np.array([lmin]),
            0.3,
            np.array([1.0]),
            smoothed_transmission(
                spectro, ARM_BLUE, np.array([lmin]), params.zenith_ang, params.sky_type
            ),
            apply_sample_factor=False,
        )
        assert result.trans[0] == pytest.approx(expected_trans_pix0[0], rel=1e-8)
        # And explicitly *not* the pixel-center convention.
        lam_center0 = lmin + 0.5 * dl
        expected_trans_center0 = snr._continuum_counts_at(
            params,
            spectro,
            ARM_BLUE,
            np.array([lam_center0]),
            0.3,
            np.array([1.0]),
            smoothed_transmission(
                spectro,
                ARM_BLUE,
                np.array([lam_center0]),
                params.zenith_ang,
                params.sky_type,
            ),
            apply_sample_factor=False,
        )
        assert result.trans[0] != pytest.approx(expected_trans_center0[0], rel=1e-6)

    def test_sample_factor_matches_helper(self, spectro, params, arm_noise):
        from pfsspecsim.etc.noise import sample_factor_for_arm

        result = snr.snr_continuum(
            params, spectro, ARM_NIR, 22.5, 0.3, arm_noise[2].noise
        )
        expected = sample_factor_for_arm(spectro, ARM_NIR, params.hgcdte_sutr)
        np.testing.assert_array_equal(result.sample_factor, expected)


class TestOiiUnionWindow:
    def test_max_offset_bound_holds_for_packaged_configs(self, spectro):
        # Empirical check backing the `_OII_MAX_OFFSET` derivation in
        # snr.py: over the full [OII]-curve/prescan z sweep (gsetc.c:2031,
        # 2129), the two doublet lines' clipped `iref`s never differ by
        # more than `_OII_MAX_OFFSET` pixels, for every arm.
        z = np.arange(0.1, 2.5001, 0.0002)
        lam0 = OII_LAMBDA[0] * (1.0 + z)
        lam1 = OII_LAMBDA[1] * (1.0 + z)
        for ia in range(spectro.N_arms):
            dl = spectro.dl[ia]
            npix = int(spectro.npix[ia])
            pos0 = (lam0 - spectro.lmin[ia]) / dl
            pos1 = (lam1 - spectro.lmin[ia]) / dl
            iref0 = np.clip(
                np.floor(pos0 - NP_WIN / 2.0).astype(np.int64), 0, npix - NP_WIN
            )
            iref1 = np.clip(
                np.floor(pos1 - NP_WIN / 2.0).astype(np.int64), 0, npix - NP_WIN
            )
            assert np.max(np.abs(iref0 - iref1)) <= snr._OII_MAX_OFFSET
