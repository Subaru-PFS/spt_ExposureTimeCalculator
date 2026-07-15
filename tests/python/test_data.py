"""Verify the extracted modeldata.npz archive (task T1).

The npz is now also exposed as installed package data, loaded via
`importlib.resources` by `pfsspecsim.etc._modeldata.load_modeldata` (T15).
This test instead locates the archive by a path relative to the repository
root, so it can check the raw extracted file directly (independent of the
installed-package loader and its caching).
"""

from pathlib import Path

import numpy as np
import pytest

MODELDATA_NPZ = (
    Path(__file__).resolve().parent.parent.parent
    / "src"
    / "pfsspecsim"
    / "etc"
    / "data"
    / "modeldata.npz"
)


@pytest.fixture(scope="module")
def modeldata():
    assert MODELDATA_NPZ.exists(), f"missing {MODELDATA_NPZ}"
    with np.load(MODELDATA_NPZ, allow_pickle=False) as data:
        yield {key: data[key] for key in data.files}


def test_uves_shapes(modeldata):
    assert modeldata["uves_lambda"].shape == (2816,)
    assert modeldata["uves_int"].shape == (2816,)


def test_oh_data_shape(modeldata):
    assert modeldata["oh_data"].shape == (698, 2)


def test_atm_trans_kp_length(modeldata):
    assert modeldata["atm_trans_kp"].shape[0] >= 40000


def test_mk_trans_3mm_length(modeldata):
    assert modeldata["mk_trans_3mm"].shape[0] >= 30001


def test_dust_norm(modeldata):
    dust_norm = modeldata["dust_norm"]
    assert dust_norm.shape == (201,)
    assert dust_norm[0] == 0.24174
    assert dust_norm[200] == 13.92363


def test_si_index_table(modeldata):
    si_index_table = modeldata["si_index_table"]
    assert si_index_table.shape == (181,)
    assert si_index_table[60] == 4.2975
