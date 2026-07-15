"""Tests for pfsspecsim.etc.io (task T10).

ECSV writer/reader for `EtcResults` tables: the `table.meta` block (task
brief's "ECSV 出力スキーマ" section) and the `overwrite`/pre-existing-file
behavior split between this module (delegates straight to
`astropy.table.Table.write`) and `engine.run_etc_files` (the actual
before-computation existence guard, tested in test_engine.py).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.table import Table

from pfsspecsim.etc import io
from pfsspecsim.etc.params import EtcParams, resolve_degrade


@pytest.fixture
def sample_table() -> Table:
    return Table(
        {
            "arm": np.array([0, 0, 1], dtype=np.int64),
            "pixel": np.array([0, 1, 0], dtype=np.int64),
            "wavelength": np.array([380.0, 380.5, 630.0]) * u.nm,
            "variance": np.array([1.0 / 3.0, 2.0 / 7.0, np.pi]) * u.electron**2,
            "sky": np.array([0.1, 0.2, 0.3]) * u.electron,
        }
    )


class TestBuildMeta:
    def test_contains_required_keys(self):
        params = EtcParams()
        meta = io.build_meta(params, "/some/path/PFS.20240714.dat")
        assert meta["etc_version"]
        assert meta["created"]
        assert meta["instr_config"] == "PFS.20240714.dat"
        assert isinstance(meta["params"], dict)
        # Top-level resolved degrade (mirroring the resolved `instr_config`
        # precedent); `meta["params"]["degrade"]` stays the raw input.
        assert meta["degrade_resolved"] == pytest.approx(resolve_degrade(params))
        assert meta["params"]["degrade"] == params.degrade

    def test_degrade_resolved_tracks_obsc_fov_dep(self):
        corrected = io.build_meta(EtcParams(degrade=2.0, obsc_fov_dep=True), "PFS.dat")
        uncorrected = io.build_meta(
            EtcParams(degrade=2.0, obsc_fov_dep=False), "PFS.dat"
        )
        assert corrected["degrade_resolved"] < 2.0
        assert uncorrected["degrade_resolved"] == 2.0
        assert corrected["params"]["degrade"] == 2.0  # raw input unchanged

    def test_path_fields_stringified(self):
        params = EtcParams(
            mag=None, mag_file=Path("/tmp/spec.dat"), outdir=Path("/tmp/out")
        )
        meta = io.build_meta(params, "PFS.dat")
        assert meta["params"]["mag_file"] == "/tmp/spec.dat"
        assert meta["params"]["outdir"] == "/tmp/out"
        assert isinstance(meta["params"]["mag_file"], str)

    def test_none_fields_pass_through(self):
        params = EtcParams()
        meta = io.build_meta(params, "PFS.dat")
        assert meta["params"]["mag_file"] is None
        assert meta["params"]["instr_config"] is None
        assert meta["params"]["oii_cat_in"] is None

    def test_all_dataclass_fields_present(self):
        import dataclasses

        params = EtcParams()
        meta = io.build_meta(params, "PFS.dat")
        field_names = {f.name for f in dataclasses.fields(EtcParams)}
        assert set(meta["params"]) == field_names


class TestWriteReadRoundTrip:
    def test_values_and_units_preserved(self, tmp_path, sample_table):
        params = EtcParams()
        path = tmp_path / "sub" / "noise.ecsv"
        io.write_table(sample_table, path, params, "PFS.20240714.dat")

        assert path.is_file()  # parent dir was created
        reloaded = io.read_table(path)

        np.testing.assert_array_equal(
            np.asarray(reloaded["variance"]), np.asarray(sample_table["variance"])
        )
        assert reloaded["wavelength"].unit == u.nm
        assert reloaded["variance"].unit == u.electron**2
        assert reloaded["sky"].unit == u.electron
        np.testing.assert_array_equal(np.asarray(reloaded["arm"]), [0, 0, 1])

    def test_meta_round_trips(self, tmp_path, sample_table):
        params = EtcParams(mag=None, mag_file=Path("spec.dat"))
        path = tmp_path / "noise.ecsv"
        io.write_table(sample_table, path, params, "PFS.20240714.dat")

        reloaded = io.read_table(path)
        assert reloaded.meta["instr_config"] == "PFS.20240714.dat"
        assert reloaded.meta["params"]["mag_file"] == "spec.dat"
        assert reloaded.meta["params"]["mag"] is None
        assert reloaded.meta["params"]["exp_num"] == 4

    def test_write_does_not_mutate_caller_table(self, tmp_path, sample_table):
        params = EtcParams()
        assert sample_table.meta == {}
        io.write_table(
            sample_table, tmp_path / "noise.ecsv", params, "PFS.20240714.dat"
        )
        assert sample_table.meta == {}

    def test_extra_meta_preserved_alongside_standard_block(
        self, tmp_path, sample_table
    ):
        # engine._compute_oii_catalog attaches a "redshift_histogram" key to
        # a table's meta *before* handing it to write_table; write_table
        # must merge its standard block in without clobbering that.
        sample_table.meta["redshift_histogram"] = {"counts": [1, 2, 3]}
        params = EtcParams()
        path = tmp_path / "oii_cat.ecsv"
        io.write_table(sample_table, path, params, "PFS.20240714.dat")

        reloaded = io.read_table(path)
        assert reloaded.meta["redshift_histogram"] == {"counts": [1, 2, 3]}
        assert reloaded.meta["etc_version"]


class TestOverwrite:
    def test_overwrite_false_raises_if_exists(self, tmp_path, sample_table):
        params = EtcParams()
        path = tmp_path / "noise.ecsv"
        io.write_table(sample_table, path, params, "PFS.20240714.dat", overwrite=True)
        with pytest.raises(OSError):
            io.write_table(
                sample_table, path, params, "PFS.20240714.dat", overwrite=False
            )

    def test_overwrite_true_replaces_existing(self, tmp_path, sample_table):
        params = EtcParams()
        path = tmp_path / "noise.ecsv"
        io.write_table(sample_table, path, params, "PFS.20240714.dat", overwrite=True)
        io.write_table(sample_table, path, params, "PFS.20240714.dat", overwrite=True)
        assert path.is_file()
