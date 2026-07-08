"""Tests for pfsspecsim.etc.materials (task T3).

Ports of gsOP_Si_abslength (gsetc.c:239-275) and gsOP_Si_indexreal
(gsetc.c:293-330). Reference values come straight from the C lookup table
(si_index_real is exact at table nodes) and from physical sanity of the
Rajkanan et al. (1979) absorption model.
"""

import numpy as np
import pytest

from pfsspecsim.etc.materials import si_abslength, si_index_real

# Representative PFS detector temperatures (K); the red/NIR CCDs run cold,
# room temperature included as an extra regime for the phonon occupancies.
TEMPERATURES = [140.0, 173.0, 293.0]


class TestSiIndexReal:
    def test_table_nodes_exact(self):
        # xfrac == 0 at 5 nm grid nodes -> the Lagrange cubic returns the
        # table entry exactly (si_table[60] and si_table[160] in gsetc.c).
        assert si_index_real(500.0) == 4.2975
        assert si_index_real(1000.0) == 3.6106

    def test_boundary_nodes_exact(self):
        # xint clipping to [1, 178] still reproduces the end nodes exactly
        # (xfrac = -1 at 200 nm, +2 at 1100 nm zero out the other weights).
        assert si_index_real(200.0) == pytest.approx(0.9317, rel=1e-12)
        assert si_index_real(1100.0) == pytest.approx(3.5797, rel=1e-12)

    def test_interpolation_between_nodes(self):
        # Off-node values stay bracketed by neighboring nodes in this
        # monotonically decreasing part of the table.
        n = si_index_real(502.5)
        assert 4.2716 < n < 4.2975

    def test_scalar_returns_float(self):
        assert isinstance(si_index_real(600.0), float)

    def test_vectorized_matches_scalar(self):
        lam = np.linspace(210.0, 1090.0, 173)
        vec = si_index_real(lam)
        assert isinstance(vec, np.ndarray)
        assert vec.shape == lam.shape
        scalars = np.array([si_index_real(float(x)) for x in lam])
        np.testing.assert_array_equal(vec, scalars)

    @pytest.mark.parametrize("lam", [199.9, 1100.1, 0.0, 5000.0])
    def test_out_of_range_raises(self, lam):
        with pytest.raises(ValueError):
            si_index_real(lam)

    def test_out_of_range_in_array_raises(self):
        with pytest.raises(ValueError):
            si_index_real(np.array([500.0, 1200.0]))


class TestSiAbslength:
    @pytest.mark.parametrize("temperature", TEMPERATURES)
    def test_positive_and_monotonic_400_1000(self, temperature):
        # Absorption length must be positive and monotonically increasing
        # with wavelength over 400-1000 nm (silicon absorbs strongly in the
        # blue and becomes transparent toward the NIR, i.e. the absorption
        # coefficient alpha decreases monotonically with wavelength).
        lam = np.linspace(400.0, 1000.0, 601)
        abslen = si_abslength(lam, temperature)
        assert np.all(abslen > 0.0)
        assert np.all(np.diff(abslen) > 0.0)

    def test_colder_silicon_is_more_transparent(self):
        # The band gap widens and phonon occupancies drop as T decreases,
        # so the absorption length at fixed wavelength grows.
        assert si_abslength(1000.0, 140.0) > si_abslength(1000.0, 293.0)

    def test_reference_value_room_temperature(self):
        # Rajkanan et al. model at 300 K, 1000 nm: alpha ~ 8.8e1 cm^-1
        # -> ~ 1.1e2 um absorption length (value locked from this port,
        # which matches a line-by-line transcription of gsetc.c:239-275 to
        # ~2e-4 relative; difference is astropy-derived hc, k_B vs the
        # truncated C literals).
        assert si_abslength(1000.0, 293.0) == pytest.approx(113.3, rel=1e-2)

    def test_scalar_returns_float(self):
        assert isinstance(si_abslength(600.0, 173.0), float)

    def test_vectorized_matches_scalar(self):
        lam = np.linspace(210.0, 1090.0, 173)
        vec = si_abslength(lam, 173.0)
        assert isinstance(vec, np.ndarray)
        assert vec.shape == lam.shape
        scalars = np.array([si_abslength(float(x), 173.0) for x in lam])
        np.testing.assert_array_equal(vec, scalars)

    @pytest.mark.parametrize("lam", [199.9, 1100.1])
    def test_out_of_range_raises(self, lam):
        with pytest.raises(ValueError):
            si_abslength(lam, 173.0)

    def test_out_of_range_in_array_raises(self):
        with pytest.raises(ValueError):
            si_abslength(np.array([500.0, 150.0]), 173.0)

    def test_range_endpoints_allowed(self):
        # C checks are strict inequalities: 200 and 1100 nm are legal.
        assert si_abslength(200.0, 173.0) > 0.0
        assert si_abslength(1100.0, 173.0) > 0.0
