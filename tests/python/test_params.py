"""Tests for `pfsspecsim.etc.params` (task T2).

`tests/mag_18.dat` is a shared acceptance fixture (do not modify): a
2-column [wavelength_nm, mag] file that is flat at AB mag 18.0 over its
whole 300-1300nm span.
"""

import pickle

import numpy as np
import pytest

from pfsspecsim.etc.params import (
    EtcParams,
    MagSpec,
    calc_obscuration,
    load_params,
    resolve_degrade,
)

from pathlib import Path

MAG_18_DAT = Path(__file__).resolve().parent.parent / "mag_18.dat"


# --- EtcParams / validate() -------------------------------------------------


def test_defaults_validate():
    EtcParams().validate()


def test_mag_xor_mag_file_both_none_raises():
    with pytest.raises(ValueError):
        EtcParams(mag=None, mag_file=None).validate()


def test_mag_xor_mag_file_both_set_raises():
    with pytest.raises(ValueError):
        EtcParams(mag=22.5, mag_file=MAG_18_DAT).validate()


def test_mag_file_alone_is_valid():
    EtcParams(mag=None, mag_file=MAG_18_DAT).validate()


def test_invalid_sky_type_raises():
    with pytest.raises(ValueError):
        EtcParams(sky_type="not_hex").validate()


def test_invalid_spectrograph_raises():
    with pytest.raises(ValueError):
        EtcParams(spectrograph="bogus").validate()


@pytest.mark.parametrize(
    "field,value",
    [
        ("seeing", -0.1),
        ("zenith_ang", 90.0),
        ("exp_time", 0.0),
        ("exp_num", 0),
        ("line_width", 0.0),
        ("degrade", 0.0),
        ("min_snr", 0.0),
        ("n_workers", 0),
    ],
)
def test_out_of_range_raises(field, value):
    params = EtcParams(**{field: value})
    with pytest.raises(ValueError):
        params.validate()


def test_etcparams_picklable():
    params = EtcParams()
    restored = pickle.loads(pickle.dumps(params))
    assert restored == params


# --- calc_obscuration / resolve_degrade -------------------------------------


def test_calc_obscuration_hand_computed():
    obsc, corr = calc_obscuration(0.45)
    # obsc = (0.36-0.27)/0.675*0.45 + 0.27
    assert obsc == pytest.approx(0.33, abs=1e-12)
    # corr = (1-obsc)/(1-0.19)
    assert corr == pytest.approx((1 - 0.33) / (1 - 0.19), abs=1e-12)


def test_resolve_degrade_matches_hand_computed():
    params = EtcParams(field_ang=0.45, degrade=1.0, obsc_fov_dep=True)
    _, corr = calc_obscuration(0.45)
    assert resolve_degrade(params) == pytest.approx(1.0 * corr, abs=1e-12)


def test_resolve_degrade_disabled_returns_degrade_unchanged():
    params = EtcParams(field_ang=0.45, degrade=2.0, obsc_fov_dep=False)
    assert resolve_degrade(params) == 2.0


# --- load_params: TOML + override priority ----------------------------------


def test_load_params_defaults_only():
    params = load_params()
    assert params == EtcParams()


def test_load_params_toml_overrides_defaults(tmp_path):
    toml_path = tmp_path / "params.toml"
    toml_path.write_text("seeing = 0.65\nexp_num = 6\n")
    params = load_params(toml_path)
    assert params.seeing == 0.65
    assert params.exp_num == 6
    assert params.zenith_ang == EtcParams().zenith_ang  # untouched default


def test_load_params_cli_overrides_beat_toml(tmp_path):
    toml_path = tmp_path / "params.toml"
    toml_path.write_text("seeing = 0.65\nexp_num = 6\n")
    params = load_params(toml_path, overrides={"exp_num": 10})
    assert params.seeing == 0.65  # from TOML
    assert params.exp_num == 10  # CLI beats TOML


def test_load_params_unknown_toml_key_raises(tmp_path):
    toml_path = tmp_path / "bad.toml"
    toml_path.write_text("not_a_real_field = 1\n")
    with pytest.raises(ValueError):
        load_params(toml_path)


def test_load_params_unknown_override_key_raises():
    with pytest.raises(ValueError):
        load_params(overrides={"not_a_real_field": 1})


def test_load_params_path_fields_coerced(tmp_path):
    toml_path = tmp_path / "params.toml"
    toml_path.write_text(f'mag_file = "{MAG_18_DAT}"\n')
    # Default `mag` is non-None too, so it must be explicitly cleared via an
    # override to satisfy the mag/mag_file XOR.
    params = load_params(toml_path, overrides={"mag": None})
    assert isinstance(params.mag_file, Path)
    assert params.mag_file == MAG_18_DAT


# --- MagSpec -----------------------------------------------------------------


def test_magspec_scalar():
    spec = MagSpec(mag=22.5)
    lam = np.array([400.0, 800.0, 1200.0])
    np.testing.assert_array_equal(spec(lam), np.full(3, 22.5))


def test_magspec_xor_raises():
    with pytest.raises(ValueError):
        MagSpec(mag=22.5, mag_file=MAG_18_DAT)
    with pytest.raises(ValueError):
        MagSpec()


def test_magspec_file_flat_18_over_full_range():
    spec = MagSpec(mag_file=MAG_18_DAT)
    lam = np.linspace(380.0, 1260.0, 200)
    np.testing.assert_allclose(spec(lam), 18.0)


def test_magspec_padding_interior_nonpositive_to_99p9(tmp_path):
    # gsetc.c:1932: any interior mag <= 0 is replaced with 99.9.
    mag_file = tmp_path / "mag_pad.dat"
    mag_file.write_text("350.0 20.0\n400.0 -1.0\n1200.0 19.0\n")
    spec = MagSpec(mag_file=mag_file)
    np.testing.assert_allclose(spec._lam, [300.0, 350.0, 400.0, 1200.0, 1300.0])
    np.testing.assert_allclose(spec._mag, [99.9, 20.0, 99.9, 19.0, 99.9])


def test_magspec_padding_endpoint_reuse_quirk(tmp_path):
    # gsetc.c:1918-1926, 1940-1947: when the data already extends past
    # 300/1300nm, the added pad point reuses the *raw* endpoint magnitude
    # (even if <= 0) rather than re-applying the interior <=0 -> 99.9 rule.
    mag_file = tmp_path / "mag_pad2.dat"
    mag_file.write_text("250.0 -2.0\n400.0 20.0\n1350.0 -3.0\n")
    spec = MagSpec(mag_file=mag_file)
    np.testing.assert_allclose(spec._lam, [249.0, 250.0, 400.0, 1350.0, 1351.0])
    np.testing.assert_allclose(spec._mag, [-2.0, 99.9, 20.0, 99.9, -3.0])


def test_magspec_file_requires_increasing_wavelength(tmp_path):
    mag_file = tmp_path / "mag_unsorted.dat"
    mag_file.write_text("400.0 20.0\n350.0 19.0\n")
    with pytest.raises(ValueError):
        MagSpec(mag_file=mag_file)
