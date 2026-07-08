"""Tests for pfsspecsim.etc.psf (task T6).

Ports of gsGeometricThroughput (gsetc.c:504-570), gsAeff (gsetc.c:576-615),
gsSpectroMTF (gsetc.c:621-667), gsSpectroDist (gsetc.c:673-690), and
gsFracTrace (gsetc.c:696-710).

The vectorized `spectro_mtf`/`spectro_dist` are the hardest port in the
whole plan (matrix-product decomposition of a triple loop), so in addition
to physical-sanity checks this file includes a literal, unvectorized,
loop-based scalar transcription of the C triple loops and cross-checks the
vectorized implementation against it to ~1e-12 -- the strongest available
guard against a silent transcription error (see task brief).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.special import j1 as _scipy_j1

from pfsspecsim.etc import psf
from pfsspecsim.etc.config import find_config_file, load_spectrograph_config
from pfsspecsim.etc.materials import si_abslength, si_index_real


@pytest.fixture(scope="module")
def spectro():
    path = find_config_file("20240714")
    return load_spectrograph_config(path, degrade=1.0)


# Arm 0 (blue, Dtype=0/Si) and arm 2 (NIR, Dtype=1/HgCdTe) exercise both
# branches of the Si-defocus block in gsSpectroMTF.
ARM_SI = 0
ARM_HGCDTE = 2


# --- Scalar (loop-based) transcription of gsSpectroMTF/gsSpectroDist ------
#
# Deliberately *not* sharing code with psf.py's vectorized implementation
# beyond the same scipy Bessel-function substitution (getJ1 -> scipy.special.j1,
# per the task brief) and materials.si_abslength/si_index_real (already unit
# tested independently in T3) -- everything else is a line-by-line, scalar,
# Python transcription of gsetc.c:621-690.


def _scalar_spectro_mtf(spectro_, i_arm, lam, u):
    """Scalar transcription of gsSpectroMTF (gsetc.c:621-667)."""
    d_spot = spectro_.diam[i_arm] / spectro_.pix[i_arm]
    mtf = 1.0
    mtf *= (
        2.0 * _scipy_j1(math.pi * d_spot * u) / (math.pi * d_spot * u)
        if abs(d_spot * u) > 1e-6
        else 1.0
    )
    mtf *= math.sin(math.pi * u) / (math.pi * u) if abs(u) > 1e-9 else 1.0

    sigma = spectro_.rms_cam[i_arm] / spectro_.pix[i_arm]
    mtf *= math.exp(-2.0 * math.pi**2 * sigma**2 * u**2)

    if spectro_.Dtype[i_arm] == 0:
        n = 3
        ddepth = spectro_.thick[i_arm] / n
        numer = denom = 0.0
        mfp = float(si_abslength(lam, spectro_.temperature[i_arm]))
        nsi = float(si_index_real(lam))
        if mfp < 1e4 * ddepth:
            d0 = (
                mfp
                * (1 - (1 + ddepth / mfp) * math.exp(-ddepth / mfp))
                / (1 - math.exp(-ddepth / mfp))
            )
        else:
            d0 = 0.5 * ddepth
        for i in range(2 * n):
            depth = d0 + ddepth * i
            contrib = math.exp(-ddepth * i / mfp) * (
                0.3 if depth > spectro_.thick[i_arm] else 1.0
            )
            rspot = depth / nsi / 2.0 / spectro_.fratio[i_arm] / spectro_.pix[i_arm]
            arg = 2 * math.pi * rspot * u
            denom += contrib
            numer += contrib * (
                0.1666666667 * math.cos(0.866025404 * arg)
                + 0.5 * math.cos(0.5 * arg)
                + 0.3333333333
            )
        mtf *= numer / denom

    mtf *= math.exp(-lam / spectro_.dl[i_arm] / spectro_.nline[i_arm] * abs(u))
    return mtf


def _scalar_spectro_dist(spectro_, i_arm, lam, pos, sigma, N):
    """Scalar transcription of gsSpectroDist (gsetc.c:673-690)."""
    nu, du = 1000, 0.005
    fr = [0.0] * N
    for iu in range(nu):
        u = du * (iu + 0.5)
        mtf1d = _scalar_spectro_mtf(spectro_, i_arm, lam, u) * math.exp(
            -2.0 * math.pi**2 * sigma**2 * u**2
        )
        for ip in range(N):
            fr[ip] += 2.0 * du * math.cos(2.0 * math.pi * u * (pos - ip)) * mtf1d
    return fr


class TestSpectroMTFScalarTranscription:
    @pytest.mark.parametrize(
        "i_arm,lam,u",
        [
            (ARM_SI, 500.0, 0.01),
            (ARM_SI, 500.0, 0.3),
            (ARM_SI, 620.0, 1.5),
            (ARM_HGCDTE, 1000.0, 0.01),
            (ARM_HGCDTE, 1200.0, 0.7),
        ],
    )
    def test_matches_scalar_reference(self, spectro, i_arm, lam, u):
        vec = psf.spectro_mtf(spectro, i_arm, [lam], [u])
        assert vec.shape == (1, 1)
        scalar = _scalar_spectro_mtf(spectro, i_arm, lam, u)
        assert vec[0, 0] == pytest.approx(scalar, rel=1e-12, abs=1e-15)


class TestSpectroDistScalarTranscription:
    @pytest.mark.parametrize("i_arm", [ARM_SI, ARM_HGCDTE])
    def test_matches_scalar_reference(self, spectro, i_arm):
        lam = 500.0 if i_arm == ARM_SI else 1000.0
        pos, sigma, N = 2.5, 0.3, 6
        vec = psf.spectro_dist(spectro, i_arm, [lam], pos, sigma, N)
        assert vec.shape == (1, N)
        scalar = _scalar_spectro_dist(spectro, i_arm, lam, pos, sigma, N)
        np.testing.assert_allclose(vec[0], scalar, rtol=1e-12, atol=1e-15)


class TestSpectroMTF:
    def test_u_to_zero_is_one(self, spectro):
        # At u=0 every factor is exactly 1: fiber/pixel guards trigger,
        # Gaussian spot -> 1, grating -> 1, and the Si defocus numer/denom
        # ratio -> 1 since arg=0 makes every cos term 1 (1/6+1/2+1/3=1).
        mtf = psf.spectro_mtf(spectro, ARM_SI, [600.0], [0.0])
        assert mtf[0, 0] == pytest.approx(1.0, abs=1e-12)

    def test_u_near_zero_approaches_one(self, spectro):
        mtf = psf.spectro_mtf(spectro, ARM_SI, [600.0], [1e-6])
        assert mtf[0, 0] == pytest.approx(1.0, abs=1e-6)

    def test_hgcdte_arm_skips_si_defocus(self, spectro):
        # Dtype=1 (HgCdTe) -- still must return a finite, sane MTF.
        mtf = psf.spectro_mtf(spectro, ARM_HGCDTE, [1000.0], psf.U_GRID)
        assert np.all(np.isfinite(mtf))
        assert mtf[0, 0] < mtf[0, -1] * 0 + 1.0  # MTF <= 1 identically true
        assert np.all(mtf <= 1.0 + 1e-12)

    def test_vectorized_over_wavelength_matches_single_calls(self, spectro):
        lams = [420.0, 500.0, 600.0]
        u = psf.U_GRID[:20]
        vec = psf.spectro_mtf(spectro, ARM_SI, lams, u)
        for i, lam in enumerate(lams):
            single = psf.spectro_mtf(spectro, ARM_SI, [lam], u)
            np.testing.assert_array_equal(vec[i], single[0])

    def test_mtf_bounded_by_one(self, spectro):
        mtf = psf.spectro_mtf(spectro, ARM_SI, [420.0, 500.0, 600.0], psf.U_GRID)
        assert np.all(mtf <= 1.0 + 1e-12)
        assert np.all(mtf >= -1.0 - 1e-12)


class TestSpectroDist:
    def test_symmetric_and_normalized(self, spectro):
        # N=32 (matching SP_PSF_LEN), pos=0.5*(N-1)=15.5 exactly centers the
        # window; sigma=0 (no extra smearing).
        N = 32
        pos = 0.5 * (N - 1)
        fr = psf.spectro_dist(spectro, ARM_SI, [550.0], pos, 0.0, N)[0]
        np.testing.assert_allclose(fr, fr[::-1], atol=1e-10)
        total = fr.sum()
        assert 0.9 < total < 1.001

    def test_vectorized_over_wavelength_matches_single_calls(self, spectro):
        lams = [420.0, 500.0, 600.0]
        N = 16
        pos = 0.5 * (N - 1)
        vec = psf.spectro_dist(spectro, ARM_SI, lams, pos, 0.5, N)
        for i, lam in enumerate(lams):
            single = psf.spectro_dist(spectro, ARM_SI, [lam], pos, 0.5, N)
            np.testing.assert_array_equal(vec[i], single[0])

    def test_precomputed_mtf_matches_internal(self, spectro):
        N = 16
        pos = 3.0
        lam = [500.0, 600.0]
        mtf = psf.spectro_mtf(spectro, ARM_SI, lam, psf.U_GRID)
        via_mtf = psf.spectro_dist(spectro, ARM_SI, lam, pos, 0.2, N, mtf=mtf)
        internal = psf.spectro_dist(spectro, ARM_SI, lam, pos, 0.2, N)
        np.testing.assert_allclose(via_mtf, internal, rtol=1e-12)


class TestFracTrace:
    def test_including_adjacent_traces_captures_more_flux(self, spectro):
        lam = np.array([450.0, 550.0, 620.0])
        tr0 = psf.frac_trace(spectro, ARM_SI, lam, tr=0)
        tr1 = psf.frac_trace(spectro, ARM_SI, lam, tr=1)
        assert np.all(tr1 >= tr0 - 1e-12)

    def test_shape_matches_lambda(self, spectro):
        lam = np.array([450.0, 550.0, 620.0, 640.0])
        frac = psf.frac_trace(spectro, ARM_SI, lam, tr=1)
        assert frac.shape == lam.shape


class TestGeometricThroughput:
    def test_monotonic_decreasing_with_r_eff(self, spectro):
        r_effs = [0.0, 0.2, 0.5, 1.0, 2.0]
        ee = [
            psf.geometric_throughput(spectro, 600.0, r_eff, 0.0, 0.0, 0.80)
            for r_eff in r_effs
        ]
        assert np.all(np.diff(ee) < 0)

    def test_monotonic_decreasing_with_decenter(self, spectro):
        decents = [0.0, 0.1, 0.2, 0.4, 0.8]
        ee = [
            psf.geometric_throughput(spectro, 600.0, 0.0, decent, 0.0, 0.80)
            for decent in decents
        ]
        assert np.all(np.diff(ee) < 0)

    def test_scalar_returns_float(self, spectro):
        ee = psf.geometric_throughput(spectro, 600.0, 0.0, 0.0, 0.0, 0.80)
        assert isinstance(ee, float)

    def test_vectorized_over_wavelength_matches_single_calls(self, spectro):
        lams = np.array([420.0, 600.0, 900.0])
        vec = psf.geometric_throughput(spectro, lams, 0.3, 0.1, 0.2, 0.80)
        for i, lam in enumerate(lams):
            single = psf.geometric_throughput(spectro, float(lam), 0.3, 0.1, 0.2, 0.80)
            assert vec[i] == pytest.approx(single, rel=1e-12)

    def test_bounded_in_unit_interval(self, spectro):
        ee = psf.geometric_throughput(spectro, 600.0, 0.3, 0.1, 0.2, 0.80)
        assert 0.0 < ee < 1.0


class TestEffectiveArea:
    def test_positive_and_bounded_by_geometric_area(self, spectro):
        geometric_area = np.pi / 4.0 * spectro.D_outer**2 * (1.0 - spectro.centobs**2)
        for i_arm in range(spectro.N_arms):
            lam = 0.5 * (spectro.lmin[i_arm] + spectro.lmax[i_arm])
            aeff = psf.effective_area(spectro, i_arm, lam, 0.0)
            assert 0.0 < aeff < geometric_area

    def test_vectorized_over_wavelength_matches_single_calls(self, spectro):
        lams = np.array([400.0, 500.0, 600.0])
        vec = psf.effective_area(spectro, ARM_SI, lams, 0.1)
        for i, lam in enumerate(lams):
            single = psf.effective_area(spectro, ARM_SI, float(lam), 0.1)
            assert vec[i] == pytest.approx(single, rel=1e-12)

    def test_clamped_outside_throughput_grid(self, spectro):
        imin = int(spectro.istart[ARM_SI])
        lo = spectro.l[imin]
        below = psf.effective_area(spectro, ARM_SI, lo - 100.0, 0.0)
        at_edge = psf.effective_area(spectro, ARM_SI, lo, 0.0)
        assert below == pytest.approx(at_edge)
