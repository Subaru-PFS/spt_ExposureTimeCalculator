"""Tests for pfsspecsim.pfsetc (task T13): the deprecated ALL_CAPS
`Etc().set_param()`/`run()` compatibility layer wrapping
`pfsspecsim.etc.engine` for old scripts/notebooks.

The z-grid sweeps (`engine._oii_curve_z_grid` etc.) are monkeypatched to a
handful of points, exactly as `test_engine.py`/`test_cli.py` do -- `Etc`'s
own defaults leave `OUTFILE_OII` (the [OII]-curve output) enabled, unlike
`EtcParams`'s own default, so every `run()` here would otherwise walk the
full ~12000-point [OII]-curve/single-line z grids.

Every `Etc()` construction now emits a `DeprecationWarning` (see
`pfsetc.Etc.__init__`); every construction below is wrapped in
`pytest.warns` so that warning is consumed rather than leaking into the
test session's warnings summary.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from pfsspecsim import pfsetc
from pfsspecsim.etc import engine

_TINY_Z = np.array([0.3, 0.8, 1.6])


@pytest.fixture(autouse=True)
def _tiny_z_grids(monkeypatch):
    """Shrink every z-sweep grid to 3 points for the whole module -- same
    "private hook" pattern as test_engine.py/test_cli.py.
    """
    monkeypatch.setattr(engine, "_oii_curve_z_grid", lambda: _TINY_Z)
    monkeypatch.setattr(engine, "_snl_z_grid", lambda: _TINY_Z)
    monkeypatch.setattr(engine, "_oii_prescan_z_grid", lambda: _TINY_Z)


def _new_etc() -> pfsetc.Etc:
    with pytest.warns(DeprecationWarning, match="deprecated"):
        return pfsetc.Etc()


@pytest.fixture
def etc_instance() -> pfsetc.Etc:
    return _new_etc()


class TestDeprecation:
    def test_construction_warns(self):
        with pytest.warns(DeprecationWarning, match="deprecated"):
            pfsetc.Etc()

    @pytest.mark.parametrize(
        "omp_num_threads,expected_n_workers",
        [
            (2, 2),  # within [1, 3]: passed through
            (16, 3),  # legacy default: capped at the 3 spectrograph arms
            (1, 1),  # serial
            (0, 1),  # subprocess-era harmless value: floored, not a crash
        ],
    )
    def test_omp_num_threads_maps_to_n_workers(
        self, omp_num_threads, expected_n_workers
    ):
        with pytest.warns(DeprecationWarning):
            etc = pfsetc.Etc(omp_num_threads=omp_num_threads)
        assert etc.omp_num_threads == omp_num_threads
        assert etc._to_new_params().n_workers == expected_n_workers
        # ... and no ETC_SRC/HOME_DIR subprocess-era attribute is created.
        assert not hasattr(etc, "ETC_SRC")


class TestCalcObscurationReExport:
    def test_reexported_from_etc_params(self):
        from pfsspecsim.etc.params import calc_obscuration as new_calc_obscuration

        assert pfsetc.calc_obscuration is new_calc_obscuration
        obsc, corr = pfsetc.calc_obscuration(0.45)
        assert 0.0 < obsc < 1.0
        assert corr > 0.0


class TestUnknownParam:
    def test_unknown_param_prints_message_and_does_not_raise(
        self, etc_instance, capsys
    ):
        result = etc_instance.set_param("NOT_A_REAL_PARAM", 1.0)
        assert result == 0
        captured = capsys.readouterr()
        assert "NOT_A_REAL_PARAM" in captured.out
        assert "can not be recognized" in captured.out
        assert "NOT_A_REAL_PARAM" not in etc_instance.params


class TestLoadParamFile:
    def test_old_style_keys_land_in_params(self, etc_instance, tmp_path):
        param_file = tmp_path / "old_style.par"
        param_file.write_text(
            "# a comment line, ignored\n" "SEEING 1.20\n" "EXP_TIME 1200\n"
        )
        result = etc_instance.load_param_file(str(param_file))
        assert result == 0
        assert etc_instance.params["SEEING"] == "1.20"
        assert etc_instance.params["EXP_TIME"] == "1200"

        # And the values actually flow through to the resolved EtcParams.
        params = etc_instance._to_new_params()
        assert params.seeing == 1.20
        assert params.exp_time == 1200.0


class TestMagFileSplit:
    def test_float_string_becomes_mag(self, etc_instance):
        etc_instance.set_param("MAG_FILE", 20.0)
        params = etc_instance._to_new_params()
        assert params.mag == 20.0
        assert params.mag_file is None

    def test_path_string_becomes_mag_file(self, etc_instance):
        etc_instance.set_param("MAG_FILE", "tests/mag_18.dat")
        params = etc_instance._to_new_params()
        assert params.mag is None
        assert params.mag_file == Path("tests/mag_18.dat")

    def test_default_mag_file_is_flat_mag(self, etc_instance):
        params = etc_instance._to_new_params()
        assert params.mag == 22.5
        assert params.mag_file is None


class TestBoolAndDashConversion:
    def test_y_n_convert_to_bool(self, etc_instance):
        etc_instance.set_param("MR_MODE", "Y")
        etc_instance.set_param("OVERWRITE", "N")
        params = etc_instance._to_new_params()
        assert params.mr_mode is True
        assert params.overwrite is False

    def test_dash_converts_to_none(self, etc_instance):
        params = etc_instance._to_new_params()
        assert params.oii_cat_in is None
        assert params.oii_cat_out is None


class TestRunProducesEcsvAtDefaultPath:
    def test_run_writes_ecsv_with_mag_in_meta(
        self, tmp_path, monkeypatch, etc_instance
    ):
        monkeypatch.chdir(tmp_path)
        etc_instance.set_param("MAG_FILE", 20.0)

        result = etc_instance.run()
        assert result == 0

        out_path = tmp_path / "out" / "ref.noise.ecsv"
        assert out_path.is_file()
        table = Table.read(out_path, format="ascii.ecsv")
        assert table.meta["params"]["mag"] == 20.0

    def test_old_attributes_restored_from_tables(
        self, tmp_path, monkeypatch, etc_instance
    ):
        monkeypatch.chdir(tmp_path)
        etc_instance.run()

        noise_tbl = Table.read(tmp_path / "out" / "ref.noise.ecsv", format="ascii.ecsv")
        np.testing.assert_array_equal(
            etc_instance.nsm_arms, np.asarray(noise_tbl["arm"])
        )
        np.testing.assert_array_equal(
            etc_instance.nsm_pixs, np.asarray(noise_tbl["pixel"])
        )
        np.testing.assert_array_equal(
            etc_instance.nsm_lams, np.asarray(noise_tbl["wavelength"])
        )
        np.testing.assert_array_equal(
            etc_instance.nsm_nois, np.asarray(noise_tbl["variance"])
        )
        np.testing.assert_array_equal(
            etc_instance.nsm_skys, np.asarray(noise_tbl["sky"])
        )

        snc_tbl = Table.read(tmp_path / "out" / "ref.snc.ecsv", format="ascii.ecsv")
        np.testing.assert_array_equal(etc_instance.snc_sncs, np.asarray(snc_tbl["snr"]))
        np.testing.assert_array_equal(
            etc_instance.snc_sigs, np.asarray(snc_tbl["signal"])
        )
        np.testing.assert_array_equal(
            etc_instance.snc_nois_mobj, np.asarray(snc_tbl["noise_variance"])
        )
        np.testing.assert_array_equal(
            etc_instance.snc_nois, np.asarray(snc_tbl["noise_variance_tot"])
        )
        np.testing.assert_array_equal(
            etc_instance.snc_conv, np.asarray(snc_tbl["conversion_factor"])
        )
        np.testing.assert_array_equal(
            etc_instance.snc_samp, np.asarray(snc_tbl["sampling_factor"])
        )

        snl_tbl = Table.read(tmp_path / "out" / "ref.snl.ecsv", format="ascii.ecsv")
        np.testing.assert_array_equal(
            etc_instance.snl_snls, np.asarray(snl_tbl["snr_tot"])
        )
        np.testing.assert_array_equal(
            etc_instance.snl_fcov, np.asarray(snl_tbl["fiber_aperture_factor"])
        )

        sno2_tbl = Table.read(tmp_path / "out" / "ref.sno2.ecsv", format="ascii.ecsv")
        np.testing.assert_array_equal(
            etc_instance.sno2_sno2, np.asarray(sno2_tbl["snr_tot"])
        )
        np.testing.assert_array_equal(etc_instance.sno2_zsps, np.asarray(sno2_tbl["z"]))

    def test_get_accessors_match_restored_attributes(
        self, tmp_path, monkeypatch, etc_instance
    ):
        monkeypatch.chdir(tmp_path)
        etc_instance.run()

        lams, nois = etc_instance.get_noise()
        assert lams is etc_instance.nsm_lams
        assert nois is etc_instance.nsm_nois

        lams, sncs = etc_instance.get_snc()
        assert lams is etc_instance.snc_lams
        assert sncs is etc_instance.snc_sncs

        lams, snls = etc_instance.get_snl()
        assert lams is etc_instance.snl_lams
        assert snls is etc_instance.snl_snls

        zsps, sno2 = etc_instance.get_sno2()
        assert zsps is etc_instance.sno2_zsps
        assert sno2 is etc_instance.sno2_sno2


class TestMakeStepMethods:
    def test_make_noise_then_make_snc_reuses_noise(
        self, tmp_path, monkeypatch, etc_instance
    ):
        monkeypatch.chdir(tmp_path)
        etc_instance.set_param("OUTFILE_SNL", "-")
        etc_instance.set_param("OUTFILE_OII", "-")

        etc_instance.make_noise_model()
        assert (tmp_path / "out" / "ref.noise.ecsv").is_file()
        assert not (tmp_path / "out" / "ref.snc.ecsv").exists()

        etc_instance.make_snc()
        assert (tmp_path / "out" / "ref.snc.ecsv").is_file()
        assert not (tmp_path / "out" / "ref.snl.ecsv").exists()
        assert not (tmp_path / "out" / "ref.sno2.ecsv").exists()

    def test_make_snc_reloads_noise_without_rewriting_file(
        self, tmp_path, monkeypatch, etc_instance
    ):
        """The C engine's reload branch (gsetc.c:1980-2001) only ever
        *reads* the noise file. So in the legacy chained pattern
        `make_noise_model(); set_param(...); make_snc()`, make_snc must
        (a) reload the noise vector from OUTFILE_NOISE rather than
        recomputing it, and (b) leave the file itself byte-for-byte (and
        mtime-) untouched -- its meta must keep describing the run that
        actually produced it. Proven by tampering with the file's variance
        column in between: the tampered values must flow into
        `snc_nois_mobj`, and the tampered file must survive unchanged.
        """
        monkeypatch.chdir(tmp_path)
        etc_instance.set_param("OUTFILE_SNL", "-")
        etc_instance.set_param("OUTFILE_OII", "-")

        etc_instance.make_noise_model()
        noise_path = tmp_path / "out" / "ref.noise.ecsv"
        assert noise_path.is_file()

        # Tamper: scale the variance column so a reload is distinguishable
        # from a recompute.
        tampered = Table.read(noise_path, format="ascii.ecsv")
        tampered["variance"] = np.asarray(tampered["variance"]) * 2.0
        tampered.write(noise_path, format="ascii.ecsv", overwrite=True)

        bytes_before = noise_path.read_bytes()
        mtime_before = os.stat(noise_path).st_mtime_ns

        # A param change between the two steps must NOT leak into the
        # noise file (this is exactly the misdescribing-meta hazard).
        etc_instance.set_param("MAG_FILE", 20.0)
        etc_instance.make_snc()

        # (a) reloaded, not recomputed: the tampered variance came through.
        np.testing.assert_array_equal(
            etc_instance.snc_nois_mobj, np.asarray(tampered["variance"])
        )
        # (b) not rewritten: content and mtime unchanged.
        assert noise_path.read_bytes() == bytes_before
        assert os.stat(noise_path).st_mtime_ns == mtime_before


class TestRunMulti:
    def test_run_multi_writes_proc_multi_named_outputs(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        etc = _new_etc()
        # Disable the two z-sweep curves: run_multi forks worker processes
        # (multiprocessing.Pool) that do not inherit this test's in-process
        # z-grid monkeypatch, so keep the nproc=2 smoke test fast by
        # skipping the ~12000-point sweeps rather than relying on that.
        etc.set_param("OUTFILE_SNL", "-")
        etc.set_param("OUTFILE_OII", "-")

        result = etc.run_multi(2, "EXP_TIME", ["450", "900"])
        assert result == 0

        for value in ("450", "900"):
            noise_path = tmp_path / "out" / f"ref.noise.ecsv.EXP_TIME.{value}"
            snc_path = tmp_path / "out" / f"ref.snc.ecsv.EXP_TIME.{value}"
            snl_path = tmp_path / "out" / f"ref.snl.ecsv.EXP_TIME.{value}"
            assert noise_path.is_file()
            assert snc_path.is_file()
            assert not snl_path.exists()

            table = Table.read(noise_path, format="ascii.ecsv")
            assert table.meta["params"]["exp_time"] == float(value)


class TestOldImportPathAliases:
    """Guard `pfsspecsim.__init__`'s `sys.modules` alias block: the
    pre-v2.0 public import paths (`pfsspecsim.pfsetc`/`.pfsspec`/
    `.dm_utils`) must keep resolving to the real modules now living in
    `pfsspecsim.legacy`/`pfsspecsim.sim`, for both `import x.y` and
    `from x.y import z` forms.
    """

    def test_import_pfsspecsim_pfsspec_module_form(self):
        import pfsspecsim.pfsspec

        import pfsspecsim.sim.pfsspec

        assert pfsspecsim.pfsspec is pfsspecsim.sim.pfsspec

    def test_import_pfsspecsim_dm_utils_module_form(self):
        import pfsspecsim.dm_utils
        import pfsspecsim.sim.dm_utils

        assert pfsspecsim.dm_utils is pfsspecsim.sim.dm_utils

    def test_import_pfsspecsim_pfsetc_module_form(self):
        import pfsspecsim.legacy.pfsetc
        import pfsspecsim.pfsetc

        assert pfsspecsim.pfsetc is pfsspecsim.legacy.pfsetc

    def test_from_pfsspecsim_pfsspec_import_pfsspec_class(self):
        from pfsspecsim.pfsspec import Pfsspec
        from pfsspecsim.sim.pfsspec import Pfsspec as NewPfsspec

        assert Pfsspec is NewPfsspec

    def test_from_pfsspecsim_pfsetc_import_etc_class(self):
        from pfsspecsim.legacy.pfsetc import Etc as NewEtc
        from pfsspecsim.pfsetc import Etc

        assert Etc is NewEtc
