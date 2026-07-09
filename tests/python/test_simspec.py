"""Tests for pfsspecsim.simspec: the modern snake_case Python API for the
spectral simulator (`SimSpecParams`/`load_params`/`run_sim_spec`), mirroring
`pfsspecsim.etc.params`'s architecture.
"""

from __future__ import annotations

import pickle
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from pfsspecsim import pfsspec
from pfsspecsim.simspec import SimSpecParams, load_params, run_sim_spec

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


# --- run_sim_spec: translation + end-to-end tests ----------------------------

_N = 128


def _snc_arrays() -> dict[str, np.ndarray]:
    """The 11 SNC columns, in the ETC's ECSV schema order (same fixture
    shape as `tests/python/test_sim_spec.py::_snc_arrays`)."""
    wav = np.linspace(400.0, 1100.0, _N)
    return {
        "arm": np.zeros(_N, dtype=int),
        "pixel": np.arange(_N),
        "wavelength": wav,
        "snr": np.full(_N, 5.0),
        "signal": np.full(_N, 100.0),
        "noise_variance": np.linspace(50.0, 80.0, _N),
        "noise_variance_tot": np.linspace(60.0, 90.0, _N),
        "input_mag": np.full(_N, 22.5),
        "conversion_factor": np.linspace(2.0e26, 3.0e26, _N),
        "sampling_factor": np.ones(_N),
        "sky": np.linspace(80.0, 120.0, _N),
    }


@pytest.fixture(scope="module")
def snc_ecsv(tmp_path_factory) -> Path:
    path = tmp_path_factory.mktemp("simspec_run_fixtures") / "snc.ecsv"
    cols = _snc_arrays()
    tbl = Table(list(cols.values()), names=list(cols.keys()))
    tbl.meta["params"] = {"exp_num": 5}
    tbl.write(path, format="ascii.ecsv")
    return path


class TestRunSimSpec:
    def test_writes_ascii_table(self, snc_ecsv, tmp_path):
        np.random.seed(1234)
        params = SimSpecParams(
            etc_file=snc_ecsv,
            out_dir=tmp_path,
            write_fits=False,
            ascii_table="sim",
            mag=21.0,
            exp_num=3,
        )
        sim = run_sim_spec(params)
        assert (tmp_path / "sim.dat").is_file()
        assert Path(sim.outdir) == tmp_path

    def test_mag_file_variant_is_accepted(self, snc_ecsv, tmp_path):
        mag_file = tmp_path / "mag.dat"
        mag_file.write_text("400.0 20.0\n1100.0 20.0\n")
        np.random.seed(1234)
        params = SimSpecParams(
            etc_file=snc_ecsv,
            out_dir=tmp_path / "out",
            write_fits=False,
            ascii_table="sim",
            mag=None,
            mag_file=mag_file,
            exp_num=3,
        )
        sim = run_sim_spec(params)
        assert (tmp_path / "out" / "sim.dat").is_file()

    def test_matches_legacy_pfsspec_byte_for_byte(self, snc_ecsv, tmp_path):
        """`run_sim_spec` must produce output identical to the equivalent
        legacy `Pfsspec().set_param(...)` call chain -- proves the
        snake_case -> legacy-dict translation is a pure relabeling."""
        np.random.seed(1234)
        params = SimSpecParams(
            etc_file=snc_ecsv,
            out_dir=tmp_path / "modern",
            write_fits=False,
            ascii_table="sim",
            mag=21.0,
            exp_num=3,
        )
        run_sim_spec(params)
        modern_out = (tmp_path / "modern" / "sim.dat").read_bytes()

        np.random.seed(1234)
        legacy = pfsspec.Pfsspec()
        legacy.set_param("etcFile", str(snc_ecsv))
        legacy.set_param("outDir", str(tmp_path / "legacy"))
        legacy.set_param("writeFits", "f")
        legacy.set_param("asciiTable", "sim")
        legacy.set_param("MAG_FILE", "21.0")
        legacy.set_param("EXP_NUM", "3")
        assert legacy.make_sim_spec() == 0
        legacy_out = (tmp_path / "legacy" / "sim.dat").read_bytes()

        assert modern_out == legacy_out
