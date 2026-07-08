"""Fast unit tests for pfsspecsim.etc.noise (task T8).

Port of `gsGetNoise` (gsetc.c:730-1059). The full numerical regression
against the C engine lives in test_noise_reference.py (marker `slow`);
this file covers the fast, structural properties: smoothed-transmission
shape/normalization, the `count<0 -> 0` clamp, the deposit-window `iref`
clip to `[0, Npix - SP_PSF_LEN]`, the HgCdTe SUTR sample factor, and the
empirical adjust-factor indexing by the *internal* arm index `ia`.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pfsspecsim.etc import noise
from pfsspecsim.etc.atmosphere import airmass, transmission
from pfsspecsim.etc.config import load_spectrograph_config
from pfsspecsim.etc.constants import ADJUST_NOISE_LR, SP_PSF_LEN
from pfsspecsim.etc.params import EtcParams

# Protected fixture; never modify (see task brief).
CONFIG_FIXTURE = Path(__file__).resolve().parents[1] / "PFS.20211220.dat"

SKY_TYPE = "11005"


@pytest.fixture(scope="module")
def spectro():
    return load_spectrograph_config(CONFIG_FIXTURE, degrade=1.0)


class TestSampleFactor:
    def test_hgcdte_sutr_arm(self, spectro):
        # Arm 2 is the HgCdTe (Dtype=1) arm in the LR config.
        assert noise.sample_factor_for_arm(spectro, 2, hgcdte_sutr=True) == 1.2

    def test_hgcdte_without_sutr(self, spectro):
        assert noise.sample_factor_for_arm(spectro, 2, hgcdte_sutr=False) == 1.0

    def test_si_arm_unaffected(self, spectro):
        assert noise.sample_factor_for_arm(spectro, 0, hgcdte_sutr=True) == 1.0
        assert noise.sample_factor_for_arm(spectro, 1, hgcdte_sutr=True) == 1.0


class TestSmoothedTransmission:
    def test_shape_and_bounds(self, spectro):
        lam = np.linspace(500.0, 640.0, 25)
        trans = noise.smoothed_transmission(spectro, 0, lam, 45.0, SKY_TYPE)
        assert trans.shape == lam.shape
        assert np.all(trans > 0.0)
        assert np.all(trans <= 1.0)

    def test_matches_pointwise_where_transmission_is_smooth(self, spectro):
        # Below the Kitt Peak table start (500nm) the line-absorption factor
        # is identically 1 and only the slowly varying continuum opacity
        # remains, so PSF smoothing is a near no-op: the weighted average
        # must collapse to (approximately) the pointwise transmission --
        # this is the normalization check (sum FR[j/5]*trans / sum FR[j/5]).
        lam = np.array([420.0, 450.0, 480.0])
        smoothed = noise.smoothed_transmission(spectro, 0, lam, 45.0, SKY_TYPE)
        pointwise = transmission(lam, 45.0, SKY_TYPE)
        np.testing.assert_allclose(smoothed, pointwise, rtol=1e-4)

    def test_smooths_sharp_absorption_features(self, spectro):
        # Around the 760nm O2 A-band the raw transmission has deep, narrow
        # lines; the PSF-smoothed version must lie strictly between the
        # local min and max of the raw transmission (a true average).
        lam = np.array([760.5])
        dl = spectro.dl[1]
        smoothed = noise.smoothed_transmission(spectro, 1, lam, 45.0, SKY_TYPE)
        window = np.linspace(lam[0] - 16 * dl, lam[0] + 16 * dl, 500)
        raw = transmission(window, 45.0, SKY_TYPE)
        assert raw.min() < smoothed[0] < raw.max()


class TestLineCountClipping:
    def test_negative_intensity_clipped_to_zero(self, spectro):
        # A couple of UVES atlas entries have (slightly) negative intensity;
        # gsetc.c clips only the derived count (`if (count<0) count=0;`,
        # gsetc.c:808), not the intensity itself.
        lam = np.array([500.0, 510.0])
        intensity = np.array([-3.0, 5.0])
        count, _frac = noise._line_counts(
            spectro,
            0,
            lam,
            intensity,
            fieldang=0.675,
            rad=0.5,
            t_exp=900.0,
            am=airmass(45.0),
            sky_type=SKY_TYPE,
            ref_airmass=1.1,
            brightness_scale=1.0,
        )
        assert count[0] == 0.0
        assert count[1] > 0.0


class TestDepositBounds:
    def test_iref_clipped_to_valid_window(self):
        npix = 100
        # Positions: far left (clips to 0), center, far right (clips to
        # npix - SP_PSF_LEN).
        pos = np.array([-15.9, 50.0, float(npix) + 14.9])
        iref = noise._line_deposit_positions(pos, npix)
        assert iref[0] == 0
        assert iref[1] == int(np.floor(50.0 - (SP_PSF_LEN / 2 - 0.5)))
        assert iref[2] == npix - SP_PSF_LEN
        # Every deposit window [iref, iref+32) stays in bounds.
        assert np.all(iref >= 0)
        assert np.all(iref + SP_PSF_LEN <= npix)

    def test_edge_lines_deposit_in_bounds(self, spectro):
        # Lines just inside the selection window at both arm edges must
        # deposit without an out-of-bounds index (np.add.at would raise).
        npix = int(spectro.npix[0])
        lmin, dl = spectro.lmin[0], spectro.dl[0]
        lam = np.array(
            [lmin - 15.0 * dl, lmin + (npix + 14.0) * dl]  # pos=-15, npix+14
        )
        acc = np.zeros(npix)
        noise._deposit_lines(
            acc,
            spectro,
            0,
            lam,
            np.array([10.0, 10.0]),
            fieldang=0.675,
            rad=0.5,
            t_exp=900.0,
            am=airmass(45.0),
            sky_type=SKY_TYPE,
            ref_airmass=1.1,
            brightness_scale=1.0,
            sample_factor=1.0,
        )
        # Deposits land in the clipped 32-pixel windows at each edge.
        assert acc[:SP_PSF_LEN].sum() != 0.0
        assert acc[npix - SP_PSF_LEN :].sum() != 0.0
        assert np.all(acc[SP_PSF_LEN : npix - SP_PSF_LEN] == 0.0)


class TestAdjustFactorIndexing:
    def test_indexed_by_internal_arm_index(self, spectro):
        # With sky lines and continuum disabled (sky_type '00000') and the
        # Moon below the horizon, the only noise contributions are dark
        # current + read noise, so per arm ia:
        #   noise = ADJUST_NOISE_LR[ia] * (dark*t*samp + read^2) * width
        # -- a closed form that fails if the adjust factor were indexed by
        # the output arm id instead of ia (the factors differ per arm).
        params = EtcParams(
            seeing=0.8,
            zenith_ang=45.0,
            field_ang=0.675,
            fiber_offset=0.0,
            moon_zenith_ang=135.0,
            exp_time=900.0,
            exp_num=8,
            mag=22.5,
            sky_type="00000",
            sky_sub_floor=0.01,
            diffuse_stray=0.02,
            degrade=1.0,
            obsc_fov_dep=False,
            hgcdte_sutr=True,
        )
        for ia in range(spectro.N_arms):
            arm = noise.compute_noise_arm(params, spectro, ia)
            samp = noise.sample_factor_for_arm(spectro, ia, True)
            expected = (
                ADJUST_NOISE_LR[ia]
                * (spectro.dark[ia] * 900.0 * samp + spectro.read[ia] ** 2)
                * spectro.width[ia]
            )
            np.testing.assert_allclose(arm.noise, expected, rtol=1e-12)
            np.testing.assert_allclose(arm.sky, 0.0, atol=1e-30)
            assert arm.sample_factor == samp


class TestComputeNoise:
    def test_returns_all_arms(self, spectro):
        params = EtcParams(
            moon_zenith_ang=135.0,
            sky_type="00000",
            mag=22.5,
            obsc_fov_dep=False,
        )
        result = noise.compute_noise(params, spectro)
        assert len(result) == spectro.N_arms
        for ia in range(spectro.N_arms):
            assert result[ia].noise.shape == (int(spectro.npix[ia]),)
            assert result[ia].sky.shape == (int(spectro.npix[ia]),)
            assert result[ia].lam.shape == (int(spectro.npix[ia]),)
