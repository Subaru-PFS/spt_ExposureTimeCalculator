"""Tests for pfsspecsim.cli.etc (task T12).

`typer.testing.CliRunner`-driven tests for the `pfs-spec etc` subcommand:
`--help`, TOML-config + CLI-option override precedence (landing in the
written ECSV's `table.meta["params"]`, per `io.build_meta`), the
`--mag`/`--mag-file` mutual-exclusivity CLI-usage error, and a fast
end-to-end run producing the three default ECSV outputs. The z-grid
sweeps (`engine._oii_curve_z_grid` etc.) are monkeypatched to a handful of
points module-wide -- the same "private hooks" pattern `test_engine.py`
uses -- so these CLI-level runs stay fast while still exercising the real
`engine.run_etc_files`/`io.write_table` code paths.
"""

from __future__ import annotations

import dataclasses
import inspect

import numpy as np
import pytest
from astropy.table import Table
from typer.testing import CliRunner

from pfsspecsim.cli import app
from pfsspecsim.cli import etc as cli_etc
from pfsspecsim.etc import engine
from pfsspecsim.etc.params import EtcParams, load_params, resolve_degrade

_TINY_Z = np.array([0.3, 0.8, 1.6])

runner = CliRunner()


@pytest.fixture(autouse=True)
def _tiny_z_grids(monkeypatch):
    """Shrink every z-sweep grid to 3 points for the whole module, exactly
    like `test_engine.py`'s fixture of the same name -- these are the
    private hooks the task brief calls out for a fast CLI smoke run.
    """
    monkeypatch.setattr(engine, "_oii_curve_z_grid", lambda: _TINY_Z)
    monkeypatch.setattr(engine, "_snl_z_grid", lambda: _TINY_Z)
    monkeypatch.setattr(engine, "_oii_prescan_z_grid", lambda: _TINY_Z)


class TestCliContractMatchesEtcParams:
    def test_every_etcparams_field_has_a_cli_option(self):
        # `etc_command`'s `overrides` dict is built by indexing `locals()`
        # with every `EtcParams` field name (see cli/etc.py's docstring/
        # `etc_command` body); if a new field were ever added to
        # `EtcParams` without a matching parameter on `etc_command`, that
        # indexing raises a `KeyError` at CLI-run time rather than at
        # import/definition time. Guard the contract directly by
        # introspection instead.
        field_names = {f.name for f in dataclasses.fields(EtcParams)}
        cli_param_names = set(inspect.signature(cli_etc.etc_command).parameters)
        missing = field_names - cli_param_names
        assert not missing, f"EtcParams field(s) with no CLI option: {missing}"


class TestHelpAndVersion:
    def test_help_exits_zero_and_lists_options(self):
        result = runner.invoke(app, ["etc", "--help"])
        assert result.exit_code == 0
        assert "--seeing" in result.output
        assert "--mag-file" in result.output
        assert "--config" in result.output

    def test_version_exits_zero(self):
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "pfs-spec" in result.output


class TestMagMutualExclusivity:
    def test_mag_and_mag_file_together_is_an_error(self, tmp_path):
        mag_file = tmp_path / "mag.dat"
        mag_file.write_text("400.0 20.0\n800.0 20.0\n")
        result = runner.invoke(
            app,
            [
                "etc",
                "--mag",
                "20.0",
                "--mag-file",
                str(mag_file),
                "--outdir",
                str(tmp_path),
            ],
        )
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_mag_file_alone_nulls_default_mag(self, tmp_path):
        mag_file = tmp_path / "mag.dat"
        mag_file.write_text("400.0 20.0\n800.0 20.0\n")
        outdir = tmp_path / "out"
        result = runner.invoke(
            app, ["etc", "--mag-file", str(mag_file), "--outdir", str(outdir)]
        )
        assert result.exit_code == 0, result.output
        table = Table.read(outdir / "ref.noise.ecsv", format="ascii.ecsv")
        assert table.meta["params"]["mag"] is None
        assert table.meta["params"]["mag_file"] == str(mag_file)


class TestTomlAndCliOverridePriority:
    def test_cli_seeing_override_beats_toml_others_kept(self, tmp_path):
        toml_path = tmp_path / "params.toml"
        toml_path.write_text("seeing = 0.65\nexp_num = 6\nmag = 20.0\n")
        outdir = tmp_path / "out"

        result = runner.invoke(
            app,
            [
                "etc",
                "--config",
                str(toml_path),
                "--seeing",
                "1.2",
                "--outdir",
                str(outdir),
            ],
        )
        assert result.exit_code == 0, result.output

        table = Table.read(outdir / "ref.noise.ecsv", format="ascii.ecsv")
        params_meta = table.meta["params"]
        assert params_meta["seeing"] == 1.2  # CLI beat TOML
        assert params_meta["exp_num"] == 6  # TOML kept (no CLI override)
        assert params_meta["mag"] == 20.0  # TOML kept (no CLI override)

    def test_cli_degrade_override_lands_in_meta_unresolved(self, tmp_path):
        """Both degrade values are recorded in the written meta block:
        `meta["params"]["degrade"]` is the raw (CLI > TOML > default)
        *input* value, while the top-level `meta["degrade_resolved"]` is
        `params.resolve_degrade`'s output -- the obscuration-corrected
        value actually applied to the throughput model (mirroring the
        resolved top-level `instr_config` precedent).
        """
        outdir = tmp_path / "out"
        result = runner.invoke(
            app, ["etc", "--degrade", "2.0", "--outdir", str(outdir)]
        )
        assert result.exit_code == 0, result.output
        table = Table.read(outdir / "ref.noise.ecsv", format="ascii.ecsv")
        assert table.meta["params"]["degrade"] == 2.0
        assert table.meta["params"]["obsc_fov_dep"] is True
        # Defaults: field_ang=0.45, obsc_fov_dep=True -> resolved value is
        # raw * calc_obscuration(0.45)'s correction factor.
        expected = resolve_degrade(load_params(overrides={"degrade": 2.0}))
        assert table.meta["degrade_resolved"] == pytest.approx(expected)
        assert table.meta["degrade_resolved"] < 2.0  # correction shrinks it

    def test_degrade_resolved_equals_raw_when_obsc_fov_dep_off(self, tmp_path):
        outdir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "etc",
                "--degrade",
                "2.0",
                "--no-obsc-fov-dep",
                "--outdir",
                str(outdir),
            ],
        )
        assert result.exit_code == 0, result.output
        table = Table.read(outdir / "ref.noise.ecsv", format="ascii.ecsv")
        assert table.meta["params"]["degrade"] == 2.0
        assert table.meta["params"]["obsc_fov_dep"] is False
        assert table.meta["degrade_resolved"] == 2.0

    def test_unknown_toml_key_is_a_clean_cli_error(self, tmp_path):
        toml_path = tmp_path / "bad.toml"
        toml_path.write_text("not_a_real_field = 1\n")
        result = runner.invoke(
            app, ["etc", "--config", str(toml_path), "--outdir", str(tmp_path)]
        )
        assert result.exit_code != 0
        assert "not_a_real_field" in result.output


class TestFastEndToEndRun:
    def test_default_run_writes_three_readable_ecsv_files(self, tmp_path):
        result = runner.invoke(app, ["etc", "--outdir", str(tmp_path)])
        assert result.exit_code == 0, result.output

        for name in ("ref.noise.ecsv", "ref.snc.ecsv", "ref.snl.ecsv"):
            path = tmp_path / name
            assert path.is_file()
            table = Table.read(path, format="ascii.ecsv")
            assert len(table) > 0

        # outfile_oii defaults to None -- no [OII]-curve file should exist.
        assert not (tmp_path / "ref.oii.ecsv").exists()

        assert "aperture_factor_800_target" in result.output
        assert "aperture_factor_800_point" in result.output
        assert str(tmp_path / "ref.noise.ecsv") in result.output
