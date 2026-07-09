"""Tests for pfsspecsim.simspec: the modern snake_case Python API for the
spectral simulator (`SimSpecParams`/`load_params`/`run_sim_spec`), mirroring
`pfsspecsim.etc.params`'s architecture.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import pytest

from pfsspecsim.simspec import SimSpecParams, load_params

# --- SimSpecParams / validate() ---------------------------------------------


def test_defaults_validate():
    SimSpecParams().validate()


def test_mag_xor_mag_file_both_none_raises():
    with pytest.raises(ValueError, match="mutually exclusive"):
        SimSpecParams(mag=None, mag_file=None).validate()


def test_mag_xor_mag_file_both_set_raises():
    with pytest.raises(ValueError, match="mutually exclusive"):
        SimSpecParams(mag=22.5, mag_file=Path("mag.dat")).validate()


def test_mag_file_alone_is_valid():
    SimSpecParams(mag=None, mag_file=Path("mag.dat")).validate()


def test_simspecparams_picklable():
    params = SimSpecParams()
    restored = pickle.loads(pickle.dumps(params))
    assert restored == params


# --- load_params: TOML + override priority ----------------------------------


def test_load_params_defaults_only():
    assert load_params() == SimSpecParams()


def test_load_params_toml_overrides_defaults(tmp_path):
    toml_path = tmp_path / "params.toml"
    toml_path.write_text('nrealize = 2\nascii_table = "fromtoml"\n')
    params = load_params(toml_path)
    assert params.nrealize == 2
    assert params.ascii_table == "fromtoml"
    assert params.ra == SimSpecParams().ra  # untouched default


def test_overrides_beat_toml(tmp_path):
    toml_path = tmp_path / "params.toml"
    toml_path.write_text("nrealize = 2\n")
    params = load_params(toml_path, overrides={"nrealize": 7})
    assert params.nrealize == 7


def test_path_fields_coerced_to_path(tmp_path):
    toml_path = tmp_path / "params.toml"
    toml_path.write_text('out_dir = "somewhere"\n')
    params = load_params(toml_path)
    assert params.out_dir == Path("somewhere")


def test_unknown_toml_key_raises(tmp_path):
    toml_path = tmp_path / "bad.toml"
    toml_path.write_text("not_a_real_field = 1\n")
    with pytest.raises(ValueError, match="not_a_real_field"):
        load_params(toml_path)


def test_unknown_override_key_raises():
    with pytest.raises(ValueError, match="bogus_field"):
        load_params(overrides={"bogus_field": 1})
