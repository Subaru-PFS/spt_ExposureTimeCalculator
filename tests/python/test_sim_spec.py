"""Tests for task T14: `Pfsspec.make_sim_spec`'s ECSV input path (+ legacy
plain-text fallback) and the `pfs-spec sim` typer CLI
(`pfsspecsim.cli`/`pfsspecsim.cli.sim`).

`make_sim_spec` now reads its `etcFile` as the Astropy ECSV table the
pure-Python ETC writes for `outfile_snc` (columns `arm`/`wavelength`/
`noise_variance`/`conversion_factor`/`sampling_factor`/`sky` standing in
for the old whitespace-format usecols (0, 2, 5, 8, 9, 10), and
`table.meta["params"]["exp_num"]` standing in for the old `#  EXP_NUM: n`
header grep); on ECSV parse failure it falls back to the old
header-grep + `np.loadtxt` path, as insurance for old-format files.

The equivalence test seeds `np.random` identically and feeds the *same*
numbers through both formats (the legacy file is written with `%.17e` so
`np.loadtxt` round-trips the exact float64 values the ECSV carries): the
two ascii outputs must be byte-identical, which pins the ECSV column
mapping AND the fallback in one shot. A meta-sensitivity test then varies
*only* `meta["params"]["exp_num"]` to prove that value (not the params
dict's own `EXP_NUM`, not anything in the columns) is what feeds the
`nexp_etc` sky-systematics rescaling.

FITS output (`writeFits='t'`, the pfs.datamodel path) is deliberately not
exercised -- T14's brief says the datamodel output path is untouched, and
`writeFits='f'` + `asciiTable` covers everything this task changed.
"""

from __future__ import annotations

import dataclasses
import inspect
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table
from typer.testing import CliRunner

from conftest import CONFIG_FIXTURE, reference_params
from pfsspecsim import sim as simspec
from pfsspecsim.sim import pfsspec
from pfsspecsim.cli import app
from pfsspecsim.cli import sim as cli_sim
from pfsspecsim.etc import engine

runner = CliRunner()

#: One synthetic "arm 0" SNC curve; wide enough (400-1100nm) that
#: `calculateFiberMagnitude`'s g/r/i/z/y bandpass integrals all see a
#: non-empty bandpass (an empty one would 0/0 into RuntimeWarnings).
_N = 128


def _snc_arrays() -> dict[str, np.ndarray]:
    """The 11 SNC columns, in the ETC's ECSV schema order, as plain
    float64 arrays shared by both fixture formats below.
    """
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


def _write_snc_ecsv(path: Path, exp_num: int = 5) -> Path:
    cols = _snc_arrays()
    tbl = Table(list(cols.values()), names=list(cols.keys()))
    tbl.meta["params"] = {"exp_num": exp_num}
    tbl.write(path, format="ascii.ecsv")
    return path


def _write_snc_legacy(path: Path, exp_num: int = 5) -> Path:
    """Same numbers in the old C engine's plain-text format: `#  KEY: value`
    header lines (`EXP_NUM` grepped by the fallback) + 11 whitespace
    columns. `%.17e` guarantees `np.loadtxt` reproduces the ECSV's float64
    values bit-for-bit, so the equivalence test below can require
    byte-identical output.
    """
    cols = _snc_arrays()
    with open(path, "w") as fd:
        fd.write("#  SEEING: 0.80\n")
        fd.write(f"#  EXP_NUM: {exp_num}\n")
        for row in zip(*cols.values()):
            arm_v, pix_v, *floats = row
            fd.write(
                f"{int(arm_v)} {int(pix_v)} "
                + " ".join("%.17e" % v for v in floats)
                + "\n"
            )
    return path


@pytest.fixture(scope="module")
def synthetic_dir(tmp_path_factory) -> Path:
    path = tmp_path_factory.mktemp("snc_fixtures")
    _write_snc_ecsv(path / "snc.ecsv")
    _write_snc_legacy(path / "snc_legacy.dat")
    return path


@pytest.fixture(scope="module")
def engine_snc_ecsv(tmp_path_factory) -> Path:
    """A *real* `pfs-spec etc` SNC ECSV: one full-resolution noise+SNC engine
    run (a few seconds; SNL/[OII] z-sweeps disabled via their outfile
    fields, so no z-grid monkeypatching is needed) over the protected
    `tests/PFS.20211220.dat` config, shared by the whole module.
    """
    outdir = tmp_path_factory.mktemp("engine_out")
    params = reference_params(
        instr_config=CONFIG_FIXTURE,
        mag=20.0,
        exp_num=4,
        outdir=outdir,
        outfile_noise=Path("noise.ecsv"),
        outfile_snc=Path("snc.ecsv"),
        outfile_snl=None,
        outfile_oii=None,
    )
    engine.run_etc_files(params)
    return outdir / "snc.ecsv"


def _run_sim(etc_file: Path, outdir: Path, seed: int = 1234, **params) -> Path:
    """Seeded `Pfsspec` run with `writeFits='f'` + `asciiTable='sim'`;
    returns the ascii output path (which `write_ascii` names `<table>.dat`).
    """
    np.random.seed(seed)
    sim = pfsspec.Pfsspec()
    sim.set_param("etcFile", str(etc_file))
    sim.set_param("outDir", str(outdir))
    sim.set_param("writeFits", "f")
    sim.set_param("asciiTable", "sim")
    sim.set_param("MAG_FILE", "21.0")
    sim.set_param("EXP_NUM", "3")
    for key, value in params.items():
        sim.set_param(key, value)
    assert sim.make_sim_spec() == 0
    return outdir / "sim.dat"


class TestMakeSimSpecEcsv:
    def test_real_engine_snc_ecsv_generates_ascii(self, engine_snc_ecsv, tmp_path):
        # Precondition: the engine really embedded exp_num in the meta the
        # new read path consumes.
        tbl = Table.read(engine_snc_ecsv, format="ascii.ecsv")
        assert tbl.meta["params"]["exp_num"] == 4

        out = _run_sim(engine_snc_ecsv, tmp_path)
        assert out.is_file()
        data = np.loadtxt(out)
        assert data.shape == (len(tbl), 6)  # one output row per SNC pixel
        np.testing.assert_allclose(  # column 1 = wavelength passthrough
            data[:, 0], np.round(np.asarray(tbl["wavelength"]), 3)
        )
        assert np.all(np.isfinite(data[:, 1]))  # flux
        assert np.all(data[:, 2] > 0)  # error

    def test_ecsv_and_legacy_fallback_produce_identical_output(
        self, synthetic_dir, tmp_path
    ):
        """Same numbers, both formats, same RNG seed -> byte-identical
        ascii output. Pins the ECSV column mapping to the old usecols
        (0, 2, 5, 8, 9, 10) AND proves the old-format fallback still
        works end-to-end.
        """
        out_ecsv = _run_sim(synthetic_dir / "snc.ecsv", tmp_path / "ecsv")
        out_legacy = _run_sim(synthetic_dir / "snc_legacy.dat", tmp_path / "legacy")
        assert out_ecsv.read_bytes() == out_legacy.read_bytes()

    def test_nexp_etc_is_read_from_ecsv_meta(self, tmp_path):
        """Two ECSVs identical except `meta["params"]["exp_num"]` must
        produce different output: `nexp_etc` scales the sky-subtraction
        systematic variance (`nsv_sys * nexp / nexp_etc`), so if the meta
        value weren't being read the outputs would coincide.
        """
        # exp_num=20 keeps nsv_sys = (floor*sqrt(nexp_etc)*sky)^2 well below
        # the noise_variance column, so nsv_rnd stays positive (no NaNs) --
        # the two outputs differ purely through the nexp/nexp_etc rescaling.
        ecsv_a = _write_snc_ecsv(tmp_path / "a.ecsv", exp_num=5)
        ecsv_b = _write_snc_ecsv(tmp_path / "b.ecsv", exp_num=20)
        out_a = _run_sim(ecsv_a, tmp_path / "out_a")
        out_b = _run_sim(ecsv_b, tmp_path / "out_b")
        assert out_a.read_bytes() != out_b.read_bytes()


class TestCliContractMatchesSimSpecParams:
    def test_every_field_has_a_cli_parameter(self):
        # `sim_command`'s `overrides` dict comprehension indexes `locals()`
        # by `_FIELD_NAMES` (== `SimSpecParams`'s field names); if a field
        # were ever added to `SimSpecParams` without a matching CLI
        # parameter, `sim_command`'s comprehension would silently skip it
        # whenever that field is actually passed. Guard the contract
        # directly by introspection.
        cli_param_names = set(inspect.signature(cli_sim.sim_command).parameters)
        field_names = {f.name for f in dataclasses.fields(simspec.SimSpecParams)}
        missing = field_names - cli_param_names
        assert not missing, f"SimSpecParams field(s) with no CLI parameter: {missing}"


class TestCliHelpAndVersion:
    def test_help_exits_zero_and_lists_options(self):
        result = runner.invoke(app, ["sim", "--help"])
        assert result.exit_code == 0
        assert "--etc-file" in result.output
        assert "--mag-file" in result.output
        assert "--config" in result.output
        assert "--ascii-table" in result.output

    def test_version_exits_zero(self):
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert "pfs-spec" in result.output


class TestCliErrors:
    def test_mag_and_mag_file_together_is_an_error(self, tmp_path):
        mag_file = tmp_path / "mag.dat"
        mag_file.write_text("400.0 20.0\n1100.0 20.0\n")
        result = runner.invoke(
            app, ["sim", "--mag", "20.0", "--mag-file", str(mag_file)]
        )
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_mag_and_mag_file_both_in_toml_is_an_error(self, tmp_path):
        toml_path = tmp_path / "params.toml"
        toml_path.write_text('mag = 20.0\nmag_file = "mag.dat"\n')
        result = runner.invoke(app, ["sim", "--config", str(toml_path)])
        assert result.exit_code != 0
        assert "mutually exclusive" in result.output

    def test_unknown_toml_key_is_a_clean_cli_error(self, tmp_path):
        toml_path = tmp_path / "bad.toml"
        toml_path.write_text("not_a_real_field = 1\n")
        result = runner.invoke(app, ["sim", "--config", str(toml_path)])
        assert result.exit_code != 0
        assert "not_a_real_field" in result.output

    def test_missing_etc_file_is_a_clean_cli_error(self, tmp_path):
        result = runner.invoke(
            app,
            [
                "sim",
                "--etc-file",
                str(tmp_path / "nope.ecsv"),
                "--no-write-fits",
                "--ascii-table",
                "sim",
                "--out-dir",
                str(tmp_path / "out"),
            ],
        )
        assert result.exit_code != 0
        assert "Error" in result.output


class TestCliRunAndPriority:
    def test_options_only_run_writes_ascii(self, synthetic_dir, tmp_path):
        outdir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "sim",
                "--etc-file",
                str(synthetic_dir / "snc.ecsv"),
                "--no-write-fits",
                "--ascii-table",
                "simcli",
                "--out-dir",
                str(outdir),
                "--mag",
                "21.0",
                "--exp-num",
                "3",
            ],
        )
        assert result.exit_code == 0, result.output
        assert (outdir / "simcli.dat").is_file()
        assert "ascii_table:" in result.output

    def test_cli_option_beats_toml_and_toml_beats_default(
        self, synthetic_dir, tmp_path
    ):
        """`--ascii-table` on the command line must beat the TOML value;
        the TOML's other keys (`etc_file`, `write_fits`, `mag`, `exp_num`)
        must survive the merge and drive the run (the built-in defaults
        would try `out/ref.snc.dat` + FITS output and fail).
        """
        toml_path = tmp_path / "params.toml"
        toml_path.write_text(
            f'etc_file = "{synthetic_dir / "snc.ecsv"}"\n'
            'ascii_table = "fromtoml"\n'
            "write_fits = false\n"
            "mag = 21.0\n"
            "exp_num = 3\n"
        )
        outdir = tmp_path / "out"
        result = runner.invoke(
            app,
            [
                "sim",
                "--config",
                str(toml_path),
                "--ascii-table",
                "fromcli",
                "--out-dir",
                str(outdir),
            ],
        )
        assert result.exit_code == 0, result.output
        assert (outdir / "fromcli.dat").is_file()  # CLI beat TOML
        assert not (outdir / "fromtoml.dat").exists()

    def test_toml_only_mag_file_no_cli_flags_works(self, synthetic_dir, tmp_path):
        """Regression: --config pointing at a TOML that sets only
        `mag_file` (no `mag`, no --mag/--mag-file CLI flags) must not
        fail with a spurious mutual-exclusivity error."""
        mag_file = tmp_path / "mag.dat"
        mag_file.write_text("400.0 20.0\n1100.0 20.0\n")
        toml_path = tmp_path / "params.toml"
        outdir = tmp_path / "out"
        toml_path.write_text(
            f'etc_file = "{synthetic_dir / "snc.ecsv"}"\n'
            f'mag_file = "{mag_file}"\n'
            'ascii_table = "fromtoml"\n'
            "write_fits = false\n"
        )
        result = runner.invoke(
            app, ["sim", "--config", str(toml_path), "--out-dir", str(outdir)]
        )
        assert result.exit_code == 0, result.output
        assert (outdir / "fromtoml.dat").is_file()
