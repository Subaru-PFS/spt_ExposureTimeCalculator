"""Tests for pfsspecsim.etc.sky (task T7).

Ports of the sky-line loading (gsetc.c:786-866), sky-continuum models
0x0-0x6 (gsetc.c:889-953), and Krisciunas & Schaefer (1991) moonlight
model (gsetc.c:956-1004) inside `gsGetNoise`.

In addition to physical-sanity checks, this file includes:

* a closed-form hand computation for sky continuum model 0x6 at
  lambda=700nm (the C formula transcribed independently, as a literal, in
  the test body/comment -- not by calling `sky.sky_continuum`);
* a scalar, line-by-line C transcription of `sky_continuum` model 0x6 and
  `moon_continuum`, cross-checked against the vectorized implementation to
  ~1e-12 (the same style used in test_psf.py for gsSpectroMTF/gsSpectroDist).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from pfsspecsim.etc.atmosphere import cont_opacity
from pfsspecsim.etc.constants import DEG_TO_RAD, PHOTONS_PER_ERG_1NM
from pfsspecsim.etc.sky import load_sky_lines, moon_continuum, sky_continuum

# Default sky_type (gsetc.c's SKYMODELS = '11006'): continuum model 0x6
# (bits 0-3), moonlight model 0x0 (bits 4-7), opacity model 0x0 (bits 8-11).
SKY_TYPE_DEFAULT = "11006"


# ---------------------------------------------------------------------------
# load_sky_lines
# ---------------------------------------------------------------------------


class TestLoadSkyLines:
    def test_shapes(self):
        sl = load_sky_lines()
        assert sl.uves_lambda_vac_nm.shape == sl.uves_intensity.shape
        assert sl.oh_lambda_vac_nm.shape == sl.oh_intensity.shape
        assert sl.uves_lambda_vac_nm.shape == (2816,)
        assert sl.oh_lambda_vac_nm.shape == (698,)

    def test_uves_air_to_vac_increases_wavelength(self):
        # Index of refraction of air > 1 in the optical/near-IR, so the
        # vacuum wavelength must be strictly greater than the (raw, air)
        # UVES table wavelength at every line.
        from pfsspecsim.etc._modeldata import load_modeldata

        raw_air = load_modeldata().uves_lambda
        sl = load_sky_lines()
        assert np.all(sl.uves_lambda_vac_nm > raw_air)
        # ...but only by a small (sub-nm-scale) amount.
        assert np.all(sl.uves_lambda_vac_nm - raw_air < 1.0)

    def test_oh_unmodified_from_modeldata(self):
        from pfsspecsim.etc._modeldata import load_modeldata

        oh = load_modeldata().oh_data
        sl = load_sky_lines()
        np.testing.assert_array_equal(sl.oh_lambda_vac_nm, oh[:, 0])
        np.testing.assert_array_equal(sl.oh_intensity, oh[:, 1])

    def test_intensities_mostly_nonnegative(self):
        # The raw UVES table has a couple of lines with slightly negative
        # tabulated intensity (noise-floor artifacts in the atlas); gsetc.c
        # does *not* clip these at load time -- it clips the derived
        # per-pixel `count` downstream instead (`if (count<0) count=0;`,
        # gsetc.c:807), which is a T8/noise.py responsibility. So this
        # loader intentionally passes the small negative values through
        # unmodified; just check they're rare and small in magnitude.
        sl = load_sky_lines()
        assert np.sum(sl.uves_intensity < 0.0) <= 5
        assert np.all(sl.uves_intensity[sl.uves_intensity < 0.0] > -1.0)
        assert np.all(sl.oh_intensity >= 0.0)

    def test_cached_returns_same_arrays(self):
        # functools.lru_cache: repeated calls should return the identical
        # (read-only) arrays, not recompute/reallocate.
        sl1 = load_sky_lines()
        sl2 = load_sky_lines()
        assert sl1.uves_lambda_vac_nm is sl2.uves_lambda_vac_nm
        assert not sl1.uves_lambda_vac_nm.flags.writeable


# ---------------------------------------------------------------------------
# sky_continuum
# ---------------------------------------------------------------------------


class TestSkyContinuumModelZero:
    def test_no_continuum(self):
        lam = np.array([400.0, 700.0, 1000.0])
        # Replace only the last hex digit (continuum-model nibble) with 0.
        sky_type_model0 = (int(SKY_TYPE_DEFAULT, 16) & ~0xF) | 0x0
        result = sky_continuum(lam, 1.2, 0.9, sky_type_model0)
        np.testing.assert_array_equal(result, np.zeros_like(lam))


class TestSkyContinuumModel6ClosedForm:
    def test_lambda_700_matches_hand_computation(self):
        # Closed-form hand computation of gsetc.c:940-947 (case 0x6) at
        # lambda=700nm, airmass=1, trans=1, independently transcribed here
        # (not by calling sky.sky_continuum):
        #
        #   mag = 24.316 - 5.199e-3*700 + 1.465e-6*700^2
        #         - 0.55*exp(-0.005*(700-594)^2) - 6.14656e9/700^4
        #       = 24.316 - 3.6393 + 0.71785 - 0.55*exp(-59.58) - 0.0255998...
        #       = 21.36895 (exp(-59.58) ~ 0, negligible)
        #   continuum = 0.01089 * 10**(0.4*(22.5-mag)) * 1e6/700^2
        #   continuum *= (-1.02278215e-3*700 + 1.77400498)   [lambda<=800]
        #   continuum *= 1e-11 * PHOTONS_PER_ERG_1NM * 700
        #   continuum *= airmass(=1) * trans(=1)
        lam = 700.0
        mag = (
            24.316
            - 5.199e-03 * lam
            + 1.465e-06 * lam**2
            - 0.55 * math.exp(-0.005 * (lam - 594.0) ** 2)
            - 6.14656e9 / lam**4
        )
        continuum = 0.01089 * 10.0 ** (0.4 * (22.5 - mag)) * 1e6 / lam**2
        continuum *= -1.02278215e-03 * lam + 1.77400498e00
        continuum *= 1e-11 * PHOTONS_PER_ERG_1NM * lam
        continuum *= 1.0 * 1.0

        expected = 0.23484556720877373  # computed independently via the above
        assert continuum == pytest.approx(expected, rel=1e-12)

        result = sky_continuum(lam, 1.0, 1.0, SKY_TYPE_DEFAULT)
        assert result == pytest.approx(expected, rel=1e-9)
        assert result == pytest.approx(continuum, rel=1e-12)


class TestSkyContinuumBoundaries:
    """Model 0x6's mag formula switches fit branches at 600nm, and its
    low-lambda correction switches from a linear ramp to a flat 1.0 at
    800nm (gsetc.c:927-947). Both are independent piecewise-fit pieces
    (not algebraically forced to agree at the switch point), so exact
    continuity is not guaranteed by construction; these tests bound how
    large the jump actually is, and otherwise check positivity/no blow-up
    on both sides.
    """

    def test_positive_across_600nm_boundary(self):
        lam = np.array([590.0, 599.0, 599.99, 600.0, 600.01, 601.0, 610.0])
        result = sky_continuum(lam, 1.1, 0.95, SKY_TYPE_DEFAULT)
        assert np.all(result > 0.0)

    def test_positive_across_800nm_boundary(self):
        lam = np.array([790.0, 799.0, 799.99, 800.0, 800.01, 801.0, 810.0])
        result = sky_continuum(lam, 1.1, 0.95, SKY_TYPE_DEFAULT)
        assert np.all(result > 0.0)

    def test_continuity_600nm_boundary(self):
        below = sky_continuum(599.999, 1.0, 1.0, SKY_TYPE_DEFAULT)
        above = sky_continuum(600.001, 1.0, 1.0, SKY_TYPE_DEFAULT)
        # The two independent polynomial fits (lambda>600 vs lambda<=600)
        # very nearly agree at the switch point: <1% relative jump.
        assert above == pytest.approx(below, rel=1e-2)

    def test_no_large_discontinuity_at_800nm_boundary(self):
        below = sky_continuum(799.999, 1.0, 1.0, SKY_TYPE_DEFAULT)
        above = sky_continuum(800.001, 1.0, 1.0, SKY_TYPE_DEFAULT)
        # The lambda<=800 linear correction does not evaluate to exactly
        # 1.0 as lambda->800 (it is ~0.9558 there), so there is a genuine
        # ~4-5% jump in the C formula itself at this boundary -- verified
        # by direct transcription, not a porting bug. Bound it generously
        # so a real regression (a much larger jump) would still be caught.
        assert above == pytest.approx(below, rel=0.1)


class TestSkyContinuumOtherModels:
    @pytest.mark.parametrize("model", [0x1, 0x2, 0x3, 0x4, 0x5, 0x6])
    def test_positive_over_wide_range(self, model):
        sky_type = (int(SKY_TYPE_DEFAULT, 16) & ~0xF) | model
        lam = np.linspace(320.0, 1300.0, 200)
        result = sky_continuum(lam, 1.2, 0.9, sky_type)
        assert np.all(result > 0.0)

    def test_model1_ir_correction_above_1040(self):
        sky_type = (int(SKY_TYPE_DEFAULT, 16) & ~0xF) | 0x1
        below = sky_continuum(1039.0, 1.0, 1.0, sky_type)
        above = sky_continuum(1041.0, 1.0, 1.0, sky_type)
        assert below > 0.0 and above > 0.0

    def test_model2_has_no_ir_correction(self):
        # Model 0x2 shares model 0x1's formula except for the lambda>1040
        # IR correction (gsetc.c:904); at 1100nm this makes a measurable
        # difference between the two models.
        sky_type1 = (int(SKY_TYPE_DEFAULT, 16) & ~0xF) | 0x1
        sky_type2 = (int(SKY_TYPE_DEFAULT, 16) & ~0xF) | 0x2
        val1 = sky_continuum(1100.0, 1.0, 1.0, sky_type1)
        val2 = sky_continuum(1100.0, 1.0, 1.0, sky_type2)
        assert val1 != pytest.approx(val2, rel=1e-6)

    def test_illegal_model_raises(self):
        sky_type = (int(SKY_TYPE_DEFAULT, 16) & ~0xF) | 0x7
        with pytest.raises(ValueError):
            sky_continuum(700.0, 1.0, 1.0, sky_type)


class TestSkyContinuumVectorScalarConsistency:
    def test_vector_matches_scalar(self):
        lam = np.linspace(320.0, 1300.0, 151)
        vec = sky_continuum(lam, 1.15, 0.92, SKY_TYPE_DEFAULT)
        scalars = np.array(
            [sky_continuum(float(x), 1.15, 0.92, SKY_TYPE_DEFAULT) for x in lam]
        )
        np.testing.assert_array_equal(vec, scalars)

    def test_scalar_returns_float(self):
        assert isinstance(sky_continuum(700.0, 1.1, 0.9, SKY_TYPE_DEFAULT), float)


# ---------------------------------------------------------------------------
# moon_continuum
# ---------------------------------------------------------------------------


class TestMoonContinuumBelowHorizon:
    def test_zenith_ang_90_is_zero(self):
        assert moon_continuum(600.0, 90.0, 60.0, 0.5, 45.0, SKY_TYPE_DEFAULT) == 0.0

    def test_zenith_ang_above_90_is_zero(self):
        assert moon_continuum(600.0, 135.0, 90.0, 0.25, 45.0, SKY_TYPE_DEFAULT) == 0.0

    def test_vectorized_zero(self):
        lam = np.array([500.0, 700.0, 900.0])
        result = moon_continuum(lam, 91.0, 60.0, 0.4, 45.0, SKY_TYPE_DEFAULT)
        np.testing.assert_array_equal(result, np.zeros_like(lam))

    def test_just_below_horizon_is_positive(self):
        assert moon_continuum(600.0, 89.9, 60.0, 0.5, 45.0, SKY_TYPE_DEFAULT) > 0.0


class TestMoonContinuumPhase:
    def test_increases_toward_full_moon(self):
        phases = np.array([0.02, 0.1, 0.2, 0.3, 0.4, 0.5])
        vals = np.array(
            [
                moon_continuum(600.0, 30.0, 60.0, float(p), 45.0, SKY_TYPE_DEFAULT)
                for p in phases
            ]
        )
        assert np.all(np.diff(vals) > 0.0)

    def test_symmetric_about_full_moon(self):
        for delta in (0.05, 0.1, 0.2, 0.4):
            before = moon_continuum(
                600.0, 30.0, 60.0, 0.5 - delta, 45.0, SKY_TYPE_DEFAULT
            )
            after = moon_continuum(
                600.0, 30.0, 60.0, 0.5 + delta, 45.0, SKY_TYPE_DEFAULT
            )
            assert before == pytest.approx(after, rel=1e-9)

    def test_phase_uses_fractional_part_only(self):
        # gsetc.c:959: lunarphase = obs->lunarphase - floor(obs->lunarphase).
        a = moon_continuum(600.0, 30.0, 60.0, 0.25, 45.0, SKY_TYPE_DEFAULT)
        b = moon_continuum(600.0, 30.0, 60.0, 3.25, 45.0, SKY_TYPE_DEFAULT)
        assert a == pytest.approx(b, rel=1e-12)


class TestMoonContinuumKVConsistency:
    def test_kv_is_cont_opacity_at_550(self):
        # The Krisciunas & Schaefer model uses kV = gsAtmContOp(obs,550,flags)
        # (gsetc.c:983); reproduce the full formula with kV taken directly
        # from atmosphere.cont_opacity(550, sky_type) and check it matches
        # moon_continuum's output exactly.
        lam = 650.0
        moon_za, moon_ang, phase, za = 40.0, 70.0, 0.35, 30.0
        kV = cont_opacity(550.0, SKY_TYPE_DEFAULT)
        assert kV == 0.12  # Mauna Kea median extinction at 550nm (T5)

        expected = _scalar_moon_continuum(
            lam, moon_za, moon_ang, phase, za, kV
        )
        result = moon_continuum(lam, moon_za, moon_ang, phase, za, SKY_TYPE_DEFAULT)
        assert result == pytest.approx(expected, rel=1e-12)


class TestMoonContinuumVectorScalarConsistency:
    def test_vector_matches_scalar(self):
        lam = np.linspace(400.0, 1200.0, 81)
        vec = moon_continuum(lam, 40.0, 70.0, 0.35, 30.0, SKY_TYPE_DEFAULT)
        scalars = np.array(
            [
                moon_continuum(float(x), 40.0, 70.0, 0.35, 30.0, SKY_TYPE_DEFAULT)
                for x in lam
            ]
        )
        np.testing.assert_array_equal(vec, scalars)

    def test_scalar_returns_float(self):
        assert isinstance(
            moon_continuum(600.0, 40.0, 70.0, 0.35, 30.0, SKY_TYPE_DEFAULT), float
        )


# ---------------------------------------------------------------------------
# Scalar, line-by-line C transcription oracles (gsetc.c:927-947, 956-1004)
# ---------------------------------------------------------------------------
#
# Deliberately not sharing code with sky.py's vectorized implementation
# beyond the constants/atmosphere functions already unit-tested in T2/T5
# -- everything else below is a scalar, Python transcription of the C
# statements, using math.* instead of numpy.


def _scalar_sky_continuum_model6(lam, airmass, trans):
    """Scalar transcription of gsetc.c:940-947 (sky continuum model 0x6)."""
    mag = (
        (24.316 if lam > 600 else 27.166)
        + (-5.199e-03 if lam > 600 else -1.419e-02) * lam
        + (1.465e-06 if lam > 600 else 8.541e-06) * lam * lam
        - 0.55 * math.exp(-0.005 * (lam - 594) * (lam - 594))
        - 6.14656e9 / lam / lam / lam / lam
    )
    continuum = 0.01089 * math.pow(10, 0.4 * (22.5 - mag)) * 1e6 / lam / lam
    continuum *= 1.0 if lam > 800 else (-1.02278215e-03 * lam + 1.77400498e00)
    continuum *= 1e-11 * PHOTONS_PER_ERG_1NM * lam
    continuum *= airmass * trans
    return continuum


def _scalar_moon_continuum(lam, moon_za_deg, moon_ang_deg, phase, za_deg, kV):
    """Scalar transcription of gsetc.c:956-1004 (Krisciunas & Schaefer 1991).

    `kV` is passed in directly (rather than recomputed via cont_opacity)
    so this oracle can be exercised independent of atmosphere.py.
    """
    if moon_za_deg >= 90:
        return 0.0

    lunarphase = phase - math.floor(phase)

    scale_rs = (math.exp(2480.0 / 550.0) - 1.0) / (
        math.exp(2480.0 / lam) - 1.0
    ) * math.pow(lam / 550.0, -7.0)
    scale_ms = (math.exp(2480.0 / 550.0) - 1.0) / (
        math.exp(2480.0 / lam) - 1.0
    ) * math.pow(lam / 550.0, -4.3)

    alpha = 360 * abs(lunarphase - 0.5)
    Istar = math.pow(10.0, -0.4 * (3.84 + 0.026 * alpha + 4e-9 * alpha**4))
    if alpha < 7:
        Istar *= 1.35 - 0.05 * alpha

    moon_ang_rad = moon_ang_deg * DEG_TO_RAD
    moon_za_rad = moon_za_deg * DEG_TO_RAD
    za_rad = za_deg * DEG_TO_RAD

    f1 = 2.29e5 * (1.06 + math.cos(moon_ang_rad) * math.cos(moon_ang_rad))
    f2 = math.pow(10.0, 6.15 - moon_ang_deg / 40.0)
    Bmoon = (
        (f1 * scale_rs + f2 * scale_ms)
        * Istar
        * math.pow(10.0, -0.4 * kV / math.sqrt(1.0 - 0.96 * math.sin(moon_za_rad) ** 2))
        * (
            1.0
            - math.pow(
                10.0, -0.4 * kV / math.sqrt(1.0 - 0.96 * math.sin(za_rad) ** 2)
            )
        )
    )
    return Bmoon / 3.408e10 * 5.48e10 / 550.0


class TestScalarTranscriptionOracle:
    @pytest.mark.parametrize(
        "lam,airmass_,trans_",
        [
            (400.0, 1.0, 1.0),
            (600.0, 1.1, 0.95),
            (700.0, 1.2, 0.9),
            (800.0, 1.05, 0.99),
            (1000.0, 1.3, 0.8),
            (1250.0, 1.15, 0.85),
        ],
    )
    def test_sky_continuum_model6(self, lam, airmass_, trans_):
        expected = _scalar_sky_continuum_model6(lam, airmass_, trans_)
        result = sky_continuum(lam, airmass_, trans_, SKY_TYPE_DEFAULT)
        assert result == pytest.approx(expected, rel=1e-12)

    @pytest.mark.parametrize(
        "lam,moon_za,moon_ang,phase,za",
        [
            (450.0, 20.0, 45.0, 0.5, 30.0),
            (600.0, 40.0, 70.0, 0.35, 30.0),
            (700.0, 60.0, 100.0, 0.05, 60.0),
            (900.0, 89.0, 150.0, 0.9, 45.0),
            (1100.0, 10.0, 30.0, 0.6, 20.0),
        ],
    )
    def test_moon_continuum(self, lam, moon_za, moon_ang, phase, za):
        kV = cont_opacity(550.0, SKY_TYPE_DEFAULT)
        expected = _scalar_moon_continuum(lam, moon_za, moon_ang, phase, za, kV)
        result = moon_continuum(lam, moon_za, moon_ang, phase, za, SKY_TYPE_DEFAULT)
        assert result == pytest.approx(expected, rel=1e-12)
