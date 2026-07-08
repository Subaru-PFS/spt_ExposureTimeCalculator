"""Tests for pfsspecsim.etc.extinction (task T3).

Port of gsGalactic_Alambda__EBV (gsetc.c:343-385): linear interpolation of
the Weingartner & Draine Milky Way R_V=3.1 dust curve on a log-spaced
wavelength grid (100 points per dex from 10 um downward).
"""

import numpy as np
import pytest

from pfsspecsim.etc._modeldata import load_modeldata
from pfsspecsim.etc.extinction import alambda_over_ebv


def test_reference_value_1000nm():
    # lambda = 1 um -> lset = 100 exactly -> norm[100] = 1.31336
    # (dust_norm[100] in modeldata.npz, norm[] at gsetc.c:358).
    assert alambda_over_ebv(1000.0) == pytest.approx(1.31336, rel=1e-6)


def test_table_endpoints():
    # 10 um -> lset = 0 -> norm[0]; 0.1 um -> lset = 200, il clipped to 199,
    # lf = 1 -> norm[200]. Both endpoints are legal (strict C inequalities).
    assert alambda_over_ebv(10000.0) == pytest.approx(0.24174, rel=1e-6)
    assert alambda_over_ebv(100.0) == pytest.approx(13.92363, rel=1e-6)


def test_v_band_magnitude_scale():
    # A_V / E(B-V) = R_V ~ 3.1 for the Milky Way diffuse ISM curve.
    assert alambda_over_ebv(550.0) == pytest.approx(3.1, abs=0.15)


def test_interpolation_bracketed_by_nodes():
    norm = load_modeldata().dust_norm
    lam = 1234.5  # nm
    lset = 100.0 / np.log(10.0) * np.log(10.0 / (lam / 1e3))
    il = int(np.floor(lset))
    value = alambda_over_ebv(lam)
    lo, hi = sorted((norm[il], norm[il + 1]))
    assert lo <= value <= hi


def test_monotonic_decreasing_in_optical_nir():
    # Redward of the 2175 A bump the extinction curve declines steadily.
    lam = np.linspace(400.0, 2500.0, 500)
    curve = alambda_over_ebv(lam)
    assert np.all(np.diff(curve) < 0.0)


def test_scalar_returns_float():
    assert isinstance(alambda_over_ebv(800.0), float)


def test_vectorized_matches_scalar():
    lam = np.geomspace(105.0, 9900.0, 211)
    vec = alambda_over_ebv(lam)
    assert isinstance(vec, np.ndarray)
    assert vec.shape == lam.shape
    scalars = np.array([alambda_over_ebv(float(x)) for x in lam])
    np.testing.assert_array_equal(vec, scalars)


@pytest.mark.parametrize("lam", [99.9, 10000.1, 1.0, 1e6])
def test_out_of_range_raises(lam):
    with pytest.raises(ValueError):
        alambda_over_ebv(lam)


def test_out_of_range_in_array_raises():
    with pytest.raises(ValueError):
        alambda_over_ebv(np.array([500.0, 50.0]))
