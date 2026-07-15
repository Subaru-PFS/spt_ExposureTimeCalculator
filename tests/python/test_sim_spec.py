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

FITS output (`writeFits='t'`, the pfs.datamodel path) was originally
deliberately not exercised here -- T14's brief said the datamodel output
path was untouched, and `writeFits='f'` + `asciiTable` covered everything
that task changed. Task 6 (below, `TestFitsWritePath`) closes that gap with
a dedicated end-to-end FITS-write test.

Task 6 also adds coverage for the three previously-untested `sky_sub_mode`
branches (`TestSkySubModes`) and the multi-object `mag_file` path
(`TestMultiObjectPath`), both in `pfsspec.py`'s `make_sim_spec`.
"""

from __future__ import annotations

import dataclasses
import inspect
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from typer.testing import CliRunner

from conftest import CONFIG_FIXTURE, reference_params, strip_ansi
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


class TestMakeSimSpecEcsvSchemaDrift:
    """`make_sim_spec` dispatches on the file's `# %ECSV` magic line, not on
    whether the ECSV read happens to raise (src/pfsspecsim/sim/pfsspec.py):
    a real ECSV file with a drifted schema (missing meta key or column)
    must propagate the parse error, not get silently rerouted to the
    legacy plain-text parser.
    """

    def test_ecsv_missing_exp_num_meta_raises(self, tmp_path):
        cols = _snc_arrays()
        tbl = Table(list(cols.values()), names=list(cols.keys()))
        tbl.meta["params"] = {}  # no "exp_num" key
        path = tmp_path / "bad_meta.ecsv"
        tbl.write(path, format="ascii.ecsv")

        with pytest.raises(KeyError):
            _run_sim(path, tmp_path / "out")

    def test_ecsv_missing_required_column_raises(self, tmp_path):
        cols = _snc_arrays()
        del cols["sky"]  # one of the 6 columns make_sim_spec extracts
        tbl = Table(list(cols.values()), names=list(cols.keys()))
        tbl.meta["params"] = {"exp_num": 5}
        path = tmp_path / "bad_column.ecsv"
        tbl.write(path, format="ascii.ecsv")

        with pytest.raises(KeyError):
            _run_sim(path, tmp_path / "out")


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
        output = strip_ansi(result.output)
        assert "--etc-file" in output
        assert "--mag-file" in output
        assert "--config" in output
        assert "--ascii-table" in output

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
        would try `out/ref.snc.ecsv` + FITS output and fail).
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


# --- Task 6: sky_sub_mode branches, FITS write path, multi-object path -----


class TestSkySubModes:
    """`Pfsspec.make_sim_spec` has *two* independent `sky_sub_mode`
    dispatches, and only "random" (the default) had any test coverage:

    1. `snr1`/`sigma1` (pfsspec.py:452-464): the sky-subtraction systematic
       term (`nsv_sys_mtrx`) is folded into the noise used to draw the
       simulated FLUX only when `sky_sub_mode == "random"` *exactly* --
       every other string (including "systematic"/"wavecalib"/"psfvar",
       and any unrecognized value) gets a smaller, systematic-free
       `sigma1`. `snr2`/`sigma2` (the reported ERROR column) always folds
       the systematic term in, regardless of `sky_sub_mode` -- so the
       ERROR column is mode-*independent* by construction, a cheap,
       deterministic invariant to pin.
    2. The per-spectrum residual-injection loop (pfsspec.py:596-691):
       "systematic"/"wavecalib"/"psfvar" each add their own structured
       residual (drawn from `nexp` synthetic exposures, then averaged) on
       top of that smaller `sigma1`'s noise draw; any other string
       (including "random" itself) adds no residual.
    """

    @pytest.mark.parametrize("mode", ["random", "systematic", "wavecalib", "psfvar"])
    def test_mode_runs_and_produces_finite_output(self, synthetic_dir, tmp_path, mode):
        out = _run_sim(synthetic_dir / "snc.ecsv", tmp_path / mode, SKY_SUB_MODE=mode)
        data = np.loadtxt(out)
        assert data.shape == (_N, 6)
        assert np.all(np.isfinite(data[:, 1]))  # flux
        assert np.all(data[:, 2] > 0)  # error

    @pytest.mark.parametrize("mode", ["systematic", "wavecalib", "psfvar"])
    def test_mode_flux_differs_from_random_but_error_column_matches(
        self, synthetic_dir, tmp_path, mode
    ):
        out_mode = _run_sim(
            synthetic_dir / "snc.ecsv", tmp_path / f"{mode}_out", SKY_SUB_MODE=mode
        )
        out_random = _run_sim(synthetic_dir / "snc.ecsv", tmp_path / f"{mode}_random")
        data_mode = np.loadtxt(out_mode)
        data_random = np.loadtxt(out_random)

        # Different residual-injection code path -> different flux draw.
        assert not np.array_equal(data_mode[:, 1], data_random[:, 1])
        # The reported error (sigma2) never depends on sky_sub_mode.
        np.testing.assert_array_equal(data_mode[:, 2], data_random[:, 2])

    def test_unrecognized_mode_does_not_raise_but_silently_understates_noise(
        self, synthetic_dir, tmp_path
    ):
        """QUIRK, preserved on the legacy `Pfsspec` surface only: a
        typo'd/unrecognized `sky_sub_mode` raises nothing there (pre-2.0
        behavior -- legacy/python_wrapper/pfsspecsim/pfsspec.py:380,452-473
        has the identical double dispatch -- kept verbatim). The modern API
        is stricter: `SimSpecParams.validate()` rejects unknown modes with
        `ValueError` (test_sim.py::test_unknown_sky_sub_mode_raises), so
        this pin covers only the legacy path driven via `Pfsspec.set_param`.

        It is NOT a silent alias for "random": dispatch #1 above keys off
        `sky_sub_mode == "random"` by exact string match, so any other
        value -- including a typo -- gets the smaller, systematic-free
        `sigma1`, same as "systematic"/"wavecalib"/"psfvar" (hence its FLUX
        differs from "random"'s). But it also never reaches any of the
        three named branches in dispatch #2, so it gets none of their
        compensating residual injection either -- it silently falls into
        the same code as "random" there, just with the smaller `sigma1`.

        Net effect: the simulated FLUX for an unrecognized mode is
        scattered using a noise budget that omits the sky-subtraction
        systematic term entirely, while the reported ERROR column still
        reports the full (systematic-inclusive) uncertainty -- flux and
        error silently go out of sync, with no warning or error raised.
        """
        out_bogus = _run_sim(
            synthetic_dir / "snc.ecsv",
            tmp_path / "bogus",
            SKY_SUB_MODE="not_a_real_mode",
        )
        out_random = _run_sim(synthetic_dir / "snc.ecsv", tmp_path / "random")
        data_bogus = np.loadtxt(out_bogus)
        data_random = np.loadtxt(out_random)

        # Not a silent alias for "random": flux differs (smaller sigma1)...
        assert not np.array_equal(data_bogus[:, 1], data_random[:, 1])
        # ...yet the reported error is identical (mode-independent by
        # construction) -- the silent flux/error mismatch this test pins.
        np.testing.assert_array_equal(data_bogus[:, 2], data_random[:, 2])


class TestFitsWritePath:
    """`writeFits='t'` (the default) drives `pfsDesign`/`pfsConfig`/
    `pfsArm`/`pfsObject` `.write()` calls (pfsspec.py:710-736) that no test
    exercised. Uses glob patterns matching each pfs.datamodel class's own
    `fileNameFormat`/`filenameFormat` (verified interactively against the
    installed `pfs.datamodel` package) rather than hardcoding a filename,
    since e.g. `pfsDesignId`/`pfsVisitHash` are content hashes.
    """

    def test_writes_readable_nonempty_fits_outputs(self, synthetic_dir, tmp_path):
        np.random.seed(1234)
        sim = pfsspec.Pfsspec()
        sim.set_param("etcFile", str(synthetic_dir / "snc.ecsv"))
        sim.set_param("outDir", str(tmp_path))
        sim.set_param("MAG_FILE", "21.0")
        sim.set_param("EXP_NUM", "3")
        # writeFits/writePfsArm both default to "t"; pfsObject is written
        # unconditionally alongside them whenever writeFits is on.
        assert sim.make_sim_spec() == 0

        design_files = sorted(tmp_path.glob("pfsDesign-0x*.fits"))
        config_files = sorted(tmp_path.glob("pfsConfig-0x*-*.fits"))
        arm_files = sorted(tmp_path.glob("pfsArm-*.fits"))
        object_files = sorted(tmp_path.glob("pfsObject-*.fits"))
        assert len(design_files) == 1
        assert len(config_files) == 1
        assert len(arm_files) == 1  # only the "b" arm is present in the fixture
        assert len(object_files) == 1

        for path in design_files + config_files + arm_files + object_files:
            with fits.open(path) as hdul:
                hdul.verify("exception")  # re-readable, well-formed FITS

        for path in arm_files + object_files:
            with fits.open(path) as hdul:
                flux = hdul["FLUX"].data
                assert flux.size > 0
                assert np.any(flux != 0)


class TestMultiObjectPath:
    """`make_sim_spec`'s `nobj > 1` branch (pfsspec.py:244-330), reached via
    a multi-column `mag_file` table, was untested: `_replicate_ids`
    (pfsspec.py:30-47) must give each of the `nobj` objects a distinct
    `objId`, and objects with different input magnitudes must end up with
    correspondingly different flux. `writeFits='f'` skips FITS I/O but
    `make_sim_spec` still populates `self.pfsObjects` (one `PfsObject` per
    input object) unconditionally, so that list is enough to check both
    properties without needing any file output.
    """

    def test_multiple_objects_get_distinct_ids_and_differing_flux(
        self, synthetic_dir, tmp_path
    ):
        # Magnitudes chosen bright enough (5/8/11, ~15.8x flux ratio per
        # step) that the true per-object flux difference swamps the
        # synthetic fixture's noise floor -- `sigma1`'s background/
        # read-noise term is ~constant regardless of source magnitude
        # (dominated by `noise_variance`, not object counts, whenever the
        # source is faint relative to the sky/detector background), so
        # fainter mags in this fixture (e.g. 18/20/22) would be swamped by
        # noise and not reliably ordered.
        mag_file = tmp_path / "multi_mag.dat"
        mag_file.write_text("400.0 5.0 8.0 11.0\n1100.0 5.0 8.0 11.0\n")

        np.random.seed(1234)
        sim = pfsspec.Pfsspec()
        sim.set_param("etcFile", str(synthetic_dir / "snc.ecsv"))
        sim.set_param("outDir", str(tmp_path / "out"))
        sim.set_param("writeFits", "f")
        sim.set_param("MAG_FILE", str(mag_file))
        sim.set_param("EXP_NUM", "3")
        assert sim.make_sim_spec() == 0

        assert len(sim.pfsObjects) == 3  # one PfsObject per input object
        obj_ids = [obj.target.objId for obj in sim.pfsObjects]
        assert len(set(obj_ids)) == 3  # _replicate_ids: distinct objIds

        # Brighter (lower mag) objects must have higher mean flux: mag
        # columns are 5/8/11 in objId order (1, 2, 3).
        mean_flux_by_obj_id = {
            obj.target.objId: float(np.mean(obj.flux)) for obj in sim.pfsObjects
        }
        assert mean_flux_by_obj_id[1] > mean_flux_by_obj_id[2] > mean_flux_by_obj_id[3]
