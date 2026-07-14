"""Task 7: smoke tests for the deprecated argparse/`@file` console scripts,
`src/pfsspecsim/scripts/run_etc.py` (`pfs-run-etc`) and
`src/pfsspecsim/scripts/gen_sim_spec.py` (`pfs-gen-sim-spec`), neither of
which had any test coverage.

Both scripts are thin argparse shims: `main()` builds a legacy `Etc`/
`Pfsspec` wrapper instance, feeds every parsed argument through
`set_param`, and drives the (already-tested-elsewhere) legacy engines --
see `test_compat.py` for `pfsetc.Etc` and `test_sim_spec.py` for
`pfsspec.Pfsspec`. So the goal here is coverage of the argparse/`@file`
plumbing (`convert_arg_line_to_args`, `--help`, argument -> `set_param`
wiring) rather than re-testing the engines themselves; both `main()`
functions are invoked in-process with a monkeypatched `sys.argv`, per
`test_cli.py`/`test_sim_spec.py`'s existing in-process-CLI pattern (no
`subprocess`/`uv run` needed).

`pfs-run-etc`'s minimal end-to-end run uses the `@file` interface directly
(the historical "Hirata-style parameter file" it exists for); the SNL
z-grid is monkeypatched to a handful of points, exactly like
`test_engine.py`/`test_cli.py`'s `_tiny_z_grids` fixture, so the run stays
fast (`OUTFILE_OII` is left at its script default of "-" -> `None`, so the
much larger [OII]-curve sweep never runs at all). `pfs-gen-sim-spec`'s
end-to-end run also uses `@file`, feeding a small synthetic SNC ECSV
fixture (same shape as `test_sim_spec.py`'s `_snc_arrays`).
"""

from __future__ import annotations

import sys

import numpy as np
import pytest
from astropy.table import Table

from pfsspecsim.etc import engine
from pfsspecsim.scripts import gen_sim_spec, run_etc

_TINY_Z = np.array([0.3, 0.8, 1.6])

_N = 128


def _snc_ecsv(path):
    wav = np.linspace(400.0, 1100.0, _N)
    cols = {
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
    tbl = Table(list(cols.values()), names=list(cols.keys()))
    tbl.meta["params"] = {"exp_num": 5}
    tbl.write(path, format="ascii.ecsv")
    return path


class TestHelpExitsZero:
    def test_run_etc_help_exits_zero(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["pfs-run-etc", "--help"])
        with pytest.raises(SystemExit) as exc_info:
            run_etc.main()
        assert exc_info.value.code == 0
        out = capsys.readouterr().out
        assert "--SEEING" in out
        assert "--OUTFILE_OII" in out

    def test_gen_sim_spec_help_exits_zero(self, monkeypatch, capsys):
        monkeypatch.setattr(sys, "argv", ["pfs-gen-sim-spec", "--help"])
        with pytest.raises(SystemExit) as exc_info:
            gen_sim_spec.main()
        assert exc_info.value.code == 0
        out = capsys.readouterr().out
        assert "--etcFile" in out
        assert "--nrealize" in out


class TestEndToEndViaAtFile:
    """Both scripts' `fromfile_prefix_chars="@"` + custom
    `convert_arg_line_to_args` (the "Hirata-style parameter file" -- one
    `PARAM value` pair per line, `#`-prefixed lines skipped) is the
    historical interface these scripts exist for; drive each through that
    interface rather than plain `--flag value` args.
    """

    def test_run_etc_at_file_drives_the_engine(self, tmp_path, monkeypatch):
        monkeypatch.setattr(engine, "_snl_z_grid", lambda: _TINY_Z)

        noise_path = tmp_path / "noise.ecsv"
        snc_path = tmp_path / "snc.ecsv"
        snl_path = tmp_path / "snl.ecsv"
        param_file = tmp_path / "run_etc.par"
        param_file.write_text(
            f"OUTFILE_NOISE {noise_path}\n"
            f"OUTFILE_SNC {snc_path}\n"
            f"OUTFILE_SNL {snl_path}\n"
            "EXP_NUM 2\n"
        )
        monkeypatch.setattr(sys, "argv", ["pfs-run-etc", f"@{param_file}"])

        assert run_etc.main() == 0

        assert noise_path.is_file()
        assert snc_path.is_file()
        assert snl_path.is_file()

    def test_gen_sim_spec_at_file_drives_the_engine(self, tmp_path, monkeypatch):
        etc_file = _snc_ecsv(tmp_path / "snc.ecsv")
        out_dir = tmp_path / "out"
        param_file = tmp_path / "gen_sim_spec.par"
        param_file.write_text(
            f"etcFile {etc_file}\n"
            f"outDir {out_dir}\n"
            "asciiTable sim\n"
            "writeFits False\n"
            "MAG_FILE 21.0\n"
            "EXP_NUM 3\n"
        )
        monkeypatch.setattr(sys, "argv", ["pfs-gen-sim-spec", f"@{param_file}"])

        np.random.seed(1234)
        assert gen_sim_spec.main() == 0

        assert (out_dir / "sim.dat").is_file()
        data = np.loadtxt(out_dir / "sim.dat")
        assert data.shape == (_N, 6)
        assert np.all(np.isfinite(data[:, 1]))  # flux
