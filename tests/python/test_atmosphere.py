"""Tests for pfsspecsim.etc.atmosphere (task T5).

Ports of gs_n_air / gs_vac2air / gs_air2vac (gsetc.c:209-228), the airmass
approximation used at gsetc.c:759, gsAtmContOp (gsetc.c:390-425, opacity
model 0x0 only), and gsAtmTrans (gsetc.c:431-491).
"""

import numpy as np
import pytest

from pfsspecsim.etc.atmosphere import (
    air2vac,
    airmass,
    cont_opacity,
    n_air,
    transmission,
    vac2air,
)

# Default sky_type (gsetc.c's SKYMODELS = '11006'): opacity model 0x0
# (bits 8-11), line absorption model 0x1 -- KP below 900nm, MK 3mm at and
# above (bits 12-15).
SKY_TYPE_DEFAULT = "11006"
# Same as SKY_TYPE_DEFAULT but line absorption model 0x0 -- Kitt Peak table
# only (bits 12-15 = 0 instead of 1).
SKY_TYPE_KP_ONLY = "10006"


class TestVacAirConversion:
    def test_round_trip(self):
        lam_air = np.array([400.0, 550.0, 656.3, 900.0, 1200.0])
        lam_vac = air2vac(lam_air)
        back = vac2air(lam_vac)
        np.testing.assert_allclose(back, lam_air, atol=1e-9)

    def test_air_lambda_shorter_than_vacuum(self):
        # n_air > 1 everywhere in the optical -> air wavelength < vacuum.
        assert vac2air(600.0) < 600.0
        assert air2vac(600.0) > 600.0

    def test_scalar_returns_float(self):
        assert isinstance(n_air(600.0), float)
        assert isinstance(vac2air(600.0), float)
        assert isinstance(air2vac(600.0), float)

    def test_vectorized_matches_scalar(self):
        lam = np.linspace(350.0, 1300.0, 97)
        vec = air2vac(lam)
        scalars = np.array([air2vac(float(x)) for x in lam])
        np.testing.assert_array_equal(vec, scalars)


class TestAirmass:
    def test_reference_value_45deg(self):
        assert airmass(45.0) == pytest.approx(1.386750, abs=1e-6)

    def test_zenith_is_unity(self):
        assert airmass(0.0) == pytest.approx(1.0, abs=1e-12)

    def test_increases_with_zenith_angle(self):
        za = np.array([0.0, 20.0, 40.0, 60.0, 75.0])
        am = airmass(za)
        assert np.all(np.diff(am) > 0.0)

    def test_scalar_returns_float(self):
        assert isinstance(airmass(30.0), float)


class TestContOpacity:
    def test_reference_values(self):
        assert cont_opacity(550.0, SKY_TYPE_DEFAULT) == 0.12
        assert cont_opacity(305.0, SKY_TYPE_DEFAULT) == 1.37

    def test_baseline_beyond_900nm(self):
        assert cont_opacity(1000.0, SKY_TYPE_DEFAULT) == 0.05

    def test_scalar_returns_float(self):
        assert isinstance(cont_opacity(600.0, SKY_TYPE_DEFAULT), float)

    def test_vectorized_matches_scalar(self):
        lam = np.linspace(305.0, 950.0, 101)
        vec = cont_opacity(lam, SKY_TYPE_DEFAULT)
        scalars = np.array([cont_opacity(float(x), SKY_TYPE_DEFAULT) for x in lam])
        np.testing.assert_array_equal(vec, scalars)

    def test_unimplemented_opacity_model_raises(self):
        with pytest.raises(NotImplementedError):
            cont_opacity(550.0, "11106")


class TestTransmission:
    def test_monotonically_decreasing_with_zenith_angle(self):
        lam = np.array([450.0, 550.0, 700.0, 950.0, 1200.0])
        za_grid = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0]
        trans = np.array([transmission(lam, za, SKY_TYPE_DEFAULT) for za in za_grid])
        # For each fixed wavelength, transmission strictly decreases as
        # zenith angle increases.
        for i in range(lam.size):
            assert np.all(np.diff(trans[:, i]) < 0.0)

    def test_bounded_in_unit_interval(self):
        # Deep telluric bands in the KP/MK tables can saturate to ~0
        # transmission, so only a non-negative lower bound is guaranteed.
        lam = np.linspace(310.0, 1400.0, 200)
        trans = transmission(lam, 45.0, SKY_TYPE_DEFAULT)
        assert np.all(trans >= 0.0)
        assert np.all(trans <= 1.0)

    def test_below_500nm_has_no_line_absorption(self):
        # Verbatim quirk (gsetc.c:447, 479): lambda<500nm skips the KP
        # table lookup entirely (factor 1), leaving only the continuum
        # term -- for both line-absorption models 0x0 and 0x1.
        for sky_type in (SKY_TYPE_KP_ONLY, SKY_TYPE_DEFAULT):
            k = cont_opacity(450.0, sky_type)
            expected = 10.0 ** (-0.4 * k / np.cos(np.radians(20.0)))
            assert transmission(450.0, 20.0, sky_type) == pytest.approx(expected)

    def test_900nm_boundary_switches_table(self):
        # Model 0x1 switches from the Kitt Peak table to the Mauna Kea 3mm
        # table strictly above 900nm (gsetc.c:462); at/below 900nm it's KP.
        below = transmission(900.0, 30.0, SKY_TYPE_DEFAULT)
        kp_only = transmission(900.0, 30.0, SKY_TYPE_KP_ONLY)
        assert below == pytest.approx(kp_only)

    def test_kp_only_model_matches_default_below_900nm(self):
        lam = np.linspace(500.0, 899.0, 50)
        np.testing.assert_allclose(
            transmission(lam, 40.0, SKY_TYPE_KP_ONLY),
            transmission(lam, 40.0, SKY_TYPE_DEFAULT),
        )

    def test_kp_only_model_errors_above_1500nm(self):
        with pytest.raises(ValueError):
            transmission(1600.0, 40.0, SKY_TYPE_KP_ONLY)

    def test_scalar_returns_float(self):
        assert isinstance(transmission(600.0, 30.0, SKY_TYPE_DEFAULT), float)

    def test_vectorized_matches_scalar(self):
        lam = np.linspace(400.0, 1300.0, 61)
        vec = transmission(lam, 25.0, SKY_TYPE_DEFAULT)
        scalars = np.array(
            [transmission(float(x), 25.0, SKY_TYPE_DEFAULT) for x in lam]
        )
        np.testing.assert_array_equal(vec, scalars)

    def test_unrecognized_line_absorption_model_raises(self):
        # bits 12-15 = 0x2: opacity model stays implemented (0x0) so this
        # exercises the line-absorption model dispatch specifically.
        with pytest.raises(ValueError):
            transmission(600.0, 30.0, "12006")
