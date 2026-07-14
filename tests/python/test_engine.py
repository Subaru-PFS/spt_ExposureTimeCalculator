"""Tests for pfsspecsim.etc.engine (task T10).

Port of `main`'s control flow (gsetc.c:1690-2177). The z-grid sweeps that
drive the [OII]-curve/single-line/[OII]-catalog-prescan outputs
(gsetc.c:2031, 2067, 2129 -- thousands of points each) are monkeypatched
to a handful of points here (`engine._oii_curve_z_grid` etc. are exactly
the "private hooks" the task brief calls for) so the full pipeline -- all
tables, all quirks -- can be exercised in well under a second of z-sweep
work; the per-arm noise computation itself (over the real, full-size
packaged 4096-pixel-per-arm config) is left untouched, since T8/T9's own
test suites already exercise that at full size without a runtime problem.
"""

from __future__ import annotations

import dataclasses
import warnings
from pathlib import Path

import numpy as np
import pytest

from pfsspecsim.etc import engine, io, psf
from pfsspecsim.etc.config import (
    find_config_file,
    load_spectrograph_config,
    spectro_arm,
)
from pfsspecsim.etc.params import EtcParams

from conftest import CONFIG_FIXTURE, MAG_FILE_FIXTURE, reference_params

_TINY_Z = np.array([0.3, 0.8, 1.6])


@pytest.fixture(autouse=True)
def _tiny_z_grids(monkeypatch):
    """Shrink every z-sweep grid to 3 points for the whole module -- the
    "private hook" the task brief asks for, exercising every column but
    none of the (physically already-tested-elsewhere) full-resolution
    sweep cost.
    """
    monkeypatch.setattr(engine, "_oii_curve_z_grid", lambda: _TINY_Z)
    monkeypatch.setattr(engine, "_snl_z_grid", lambda: _TINY_Z)
    monkeypatch.setattr(engine, "_oii_prescan_z_grid", lambda: _TINY_Z)


@pytest.fixture
def oii_cat_in(tmp_path) -> Path:
    """Tiny synthetic [OII] catalog: C column order `id z r_eff ROII FOII
    contOII sigma` (gsetc.c:2147). Row 1 is a bright, easily detected
    source; row 2 has `r_eff=0` (point source, still eligible); row 3's
    flux is far below any plausible MDLF and must be dropped.
    """
    path = tmp_path / "oii_cat.txt"
    path.write_text(
        "1 0.5 0.30 1.0 5.0e-16 0.0 70\n"
        "2 0.9 0.00 1.2 1.0e-16 0.0 70\n"
        "3 1.6 0.30 1.0 1.0e-20 0.0 70\n"
    )
    return path


def _base_params(tmp_path, **overrides) -> EtcParams:
    values = dict(
        exp_time=900.0,
        exp_num=4,
        outdir=tmp_path,
        outfile_noise=Path("noise.ecsv"),
        outfile_snc=Path("snc.ecsv"),
        outfile_snl=Path("snl.ecsv"),
        outfile_oii=Path("oii_curve.ecsv"),
    )
    values.update(overrides)
    return EtcParams(**values)


class TestRunEtcSmoke:
    def test_all_tables_and_meta_produced(self, tmp_path, oii_cat_in):
        params = _base_params(
            tmp_path,
            oii_cat_in=oii_cat_in,
            oii_cat_out=tmp_path / "oii_cat_out.ecsv",
        )
        results = engine.run_etc_files(params)

        assert results.noise is not None and len(results.noise) > 0
        assert results.snc is not None and len(results.snc) > 0
        assert results.snl is not None and len(results.snl) == _TINY_Z.size
        assert results.oii_curve is not None and len(results.oii_curve) == _TINY_Z.size
        assert results.oii_catalog is not None
        assert results.oii_histogram is not None
        assert results.oii_histogram.shape == (24,)
        assert results.oii_n_targets == len(results.oii_catalog)

        for name in ("noise.ecsv", "snc.ecsv", "snl.ecsv", "oii_curve.ecsv"):
            assert (tmp_path / name).is_file()
        assert (tmp_path / "oii_cat_out.ecsv").is_file()

        # Row 3 (flux 1e-20) must be dropped; rows 1/2 (bright / point
        # source r_eff=0) must survive.
        assert set(results.oii_catalog["obj_id"]) == {1, 2}

        # `io.write_table` writes a *copy* with meta attached (it must not
        # mutate the in-memory `EtcResults` tables, see test_io.py), so the
        # meta block is checked on the round-tripped, on-disk files.
        for name in ("noise.ecsv", "snc.ecsv", "snl.ecsv", "oii_curve.ecsv"):
            table = io.read_table(tmp_path / name)
            assert table.meta["etc_version"]
            assert table.meta["created"]
            assert table.meta["instr_config"] == "PFS.20240714.dat"
            assert table.meta["params"]["exp_num"] == 4
            assert table.meta["params"]["mag_file"] is None

        oii_cat_table = io.read_table(tmp_path / "oii_cat_out.ecsv")
        assert "redshift_histogram" in oii_cat_table.meta

    def test_skips_outputs_with_none_paths(self, tmp_path):
        params = _base_params(
            tmp_path, outfile_snc=None, outfile_snl=None, outfile_oii=None
        )
        results = engine.run_etc(params)
        assert results.snc is None
        assert results.snl is None
        assert results.oii_curve is None
        assert results.oii_catalog is None
        assert results.oii_histogram is None
        assert results.noise is not None

    def test_aperture_factors_match_geometric_throughput(self, tmp_path):
        params = _base_params(tmp_path)
        results = engine.run_etc(params)

        spectro = load_spectrograph_config(find_config_file("20240714"), degrade=1.0)
        expected_target = psf.geometric_throughput(
            spectro,
            800.0,
            params.reff,
            params.fiber_offset,
            params.field_ang,
            params.seeing,
        )
        expected_point = psf.geometric_throughput(
            spectro, 800.0, 0.0, params.fiber_offset, params.field_ang, params.seeing
        )
        assert results.aperture_factor_800_target == pytest.approx(expected_target)
        assert results.aperture_factor_800_point == pytest.approx(expected_point)
        # A larger effective radius must not increase encircled energy.
        assert results.aperture_factor_800_target <= results.aperture_factor_800_point


class TestOiiCatalogValidationWarning:
    """The [OII] catalog pathway (`oii_cat_in`) has unit/oracle test
    coverage but no end-to-end C-reference regression gate (the C fixture
    at `tests/gsetc_params.txt` was generated with `InFileOII='-'`, i.e.
    no catalog). `run_etc` must warn
    once whenever the catalog pathway is actually entered, and must stay
    silent for ordinary runs (including ones that enable the separately-
    gated sno2/`outfile_oii` curve, which is *not* what this warning is
    about).
    """

    def test_warns_when_catalog_pathway_used(self, tmp_path, oii_cat_in):
        params = _base_params(
            tmp_path,
            oii_cat_in=oii_cat_in,
            oii_cat_out=tmp_path / "oii_cat_out.ecsv",
        )
        with pytest.warns(UserWarning, match="\\[OII\\] emitter-catalog pathway"):
            engine.run_etc(params)

    def test_no_warning_for_ordinary_run_without_catalog(self, tmp_path):
        # sno2 curve (outfile_oii) is enabled by _base_params's defaults --
        # it has its own (separately gated) C-reference test coverage and
        # must not trip the catalog-only warning.
        params = _base_params(tmp_path)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            engine.run_etc(params)
        catalog_warnings = [
            w
            for w in caught
            if issubclass(w.category, UserWarning)
            and "[OII] emitter-catalog pathway" in str(w.message)
        ]
        assert not catalog_warnings


class TestNoiseComputeOrReload:
    def test_noise_reused_round_trip(self, tmp_path):
        params_compute = _base_params(
            tmp_path, outfile_snc=None, outfile_snl=None, outfile_oii=None
        )
        computed = engine.run_etc_files(params_compute)

        params_reload = _base_params(
            tmp_path,
            outfile_snc=None,
            outfile_snl=None,
            outfile_oii=None,
            noise_reused=True,
        )
        reloaded = engine.run_etc(params_reload)

        np.testing.assert_array_equal(
            np.asarray(computed.noise["variance"]),
            np.asarray(reloaded.noise["variance"]),
        )
        np.testing.assert_array_equal(
            np.asarray(computed.noise["sky"]), np.asarray(reloaded.noise["sky"])
        )
        np.testing.assert_array_equal(
            np.asarray(computed.noise["arm"]), np.asarray(reloaded.noise["arm"])
        )
        np.testing.assert_array_equal(
            np.asarray(computed.noise["pixel"]), np.asarray(reloaded.noise["pixel"])
        )

    def test_noise_reused_without_outfile_noise_raises(self, tmp_path):
        params = _base_params(
            tmp_path,
            outfile_noise=None,
            outfile_snc=None,
            outfile_snl=None,
            outfile_oii=None,
            noise_reused=True,
        )
        with pytest.raises(ValueError, match="noise_reused"):
            engine.run_etc(params)

    def test_reload_rejects_row_count_mismatch(self, tmp_path):
        spectro = load_spectrograph_config(find_config_file("20240714"), degrade=1.0)
        from astropy.table import Table
        from astropy import units as u

        bad_table = Table(
            {
                "arm": np.zeros(10, dtype=np.int64),
                "pixel": np.arange(10, dtype=np.int64),
                "wavelength": np.linspace(380.0, 400.0, 10) * u.nm,
                "variance": np.ones(10) * u.electron**2,
                "sky": np.ones(10) * u.electron,
            }
        )
        with pytest.raises(ValueError, match="noise_reused"):
            engine._noise_arrays_from_table(spectro, bad_table)


class TestOverwriteGuard:
    def test_overwrite_false_raises_before_compute(self, tmp_path, monkeypatch):
        # Pre-create only the noise output file; the guard must fire from
        # that alone (gsetc.c-era wrapper quirk: only noise/snc/snl are
        # even checked -- see engine.run_etc_files's docstring).
        (tmp_path).mkdir(exist_ok=True)
        (tmp_path / "noise.ecsv").write_text("placeholder")

        def _boom(_params):
            raise AssertionError("run_etc must not be called when the guard fires")

        monkeypatch.setattr(engine, "run_etc", _boom)

        params = _base_params(tmp_path, overwrite=False)
        with pytest.raises(FileExistsError, match="noise.ecsv"):
            engine.run_etc_files(params)

    def test_overwrite_false_ignores_oii_outputs(self, tmp_path, oii_cat_in):
        # Verbatim quirk (pfsetc.py:236-249): the OII curve/catalog outputs
        # are never checked by the overwrite guard, even though this
        # function does go on to (over)write them.
        oii_curve_path = tmp_path / "oii_curve.ecsv"
        oii_curve_path.parent.mkdir(parents=True, exist_ok=True)
        oii_curve_path.write_text("placeholder")

        params = _base_params(
            tmp_path,
            outfile_snc=None,
            outfile_snl=None,
            oii_cat_in=oii_cat_in,
            oii_cat_out=tmp_path / "oii_cat_out.ecsv",
            overwrite=False,
        )
        # Must NOT raise, and must clobber the placeholder file.
        results = engine.run_etc_files(params)
        assert results.oii_curve is not None
        assert oii_curve_path.read_text() != "placeholder"

    def test_overwrite_true_allows_rerun(self, tmp_path):
        params = _base_params(
            tmp_path, outfile_snc=None, outfile_snl=None, outfile_oii=None
        )
        engine.run_etc_files(params)
        engine.run_etc_files(params)  # must not raise


class TestQuirks:
    def test_oii_curve_aperture_factor_uses_fieldang_zero(self, tmp_path):
        # gsetc.c:2049 QUIRK: the [OII] curve's aperture-factor column is
        # computed at fieldang=0 regardless of params.field_ang, unlike the
        # SNL curve's (gsetc.c:2084, real fieldang) -- see snr.py's module
        # docstring and the task brief.
        params = _base_params(tmp_path, field_ang=0.675)
        results = engine.run_etc(params)

        spectro = load_spectrograph_config(find_config_file("20240714"), degrade=1.0)
        z = _TINY_Z
        lam_mid = engine._OII_MID * (1.0 + z)
        expected_fieldang_zero = psf.geometric_throughput(
            spectro, lam_mid, params.reff, params.fiber_offset, 0.0, params.seeing
        )
        expected_real_fieldang = psf.geometric_throughput(
            spectro,
            lam_mid,
            params.reff,
            params.fiber_offset,
            params.field_ang,
            params.seeing,
        )
        np.testing.assert_allclose(
            results.oii_curve["fiber_aperture_factor"], expected_fieldang_zero
        )
        assert not np.allclose(
            np.asarray(results.oii_curve["fiber_aperture_factor"]),
            expected_real_fieldang,
        )

    def test_snl_uses_real_fieldang(self, tmp_path):
        params = _base_params(tmp_path, field_ang=0.675)
        results = engine.run_etc(params)

        spectro = load_spectrograph_config(find_config_file("20240714"), degrade=1.0)
        lam = engine._SNL_LAMBDA_REST * (1.0 + _TINY_Z)
        expected = psf.geometric_throughput(
            spectro,
            lam,
            params.reff,
            params.fiber_offset,
            params.field_ang,
            params.seeing,
        )
        np.testing.assert_allclose(results.snl["fiber_aperture_factor"], expected)

    def test_snc_wavelength_is_pixel_left_edge(self, tmp_path):
        params = _base_params(tmp_path, outfile_snl=None, outfile_oii=None)
        results = engine.run_etc(params)
        spectro = load_spectrograph_config(find_config_file("20240714"), degrade=1.0)
        arm0 = results.snc[results.snc["arm"] == spectro_arm(0, spectro.MR)]
        assert arm0["wavelength"][0] == pytest.approx(spectro.lmin[0])

    def test_noise_wavelength_is_pixel_center(self, tmp_path):
        params = _base_params(
            tmp_path, outfile_snc=None, outfile_snl=None, outfile_oii=None
        )
        results = engine.run_etc(params)
        spectro = load_spectrograph_config(find_config_file("20240714"), degrade=1.0)
        arm0 = results.noise[results.noise["arm"] == spectro_arm(0, spectro.MR)]
        expected_first = spectro.lmin[0] + 0.5 * spectro.dl[0]
        assert arm0["wavelength"][0] == pytest.approx(expected_first)

    def test_snr_scaled_by_sqrt_n_exp(self):
        # snr.snr_oii itself returns a single-exposure SNR (per its own
        # docstring); the engine must multiply by sqrt(exp_num)
        # (gsetc.c:2038). Isolated from noise.compute_noise_arm's own
        # (unrelated) exp_num dependence -- the sky-subtraction systematic
        # floor's sqrt(n_exp) rescaling, gsetc.c:1834 -- by supplying a
        # fixed, hand-built noise vector directly to `_compute_oii_curve`
        # rather than going through the full `run_etc` noise computation.
        import math

        spectro = load_spectrograph_config(find_config_file("20240714"), degrade=1.0)
        noise_arrays = [np.full(int(n), 100.0) for n in spectro.npix]

        params1 = EtcParams(exp_num=1, outfile_oii=Path("x.ecsv"))
        params4 = EtcParams(exp_num=4, outfile_oii=Path("x.ecsv"))

        curve1 = engine._compute_oii_curve(params1, spectro, noise_arrays)
        curve4 = engine._compute_oii_curve(params4, spectro, noise_arrays)

        ratio = np.asarray(curve4["snr_tot"]) / np.asarray(curve1["snr_tot"])
        finite = np.isfinite(ratio) & (np.asarray(curve1["snr_tot"]) > 0)
        assert finite.any()
        np.testing.assert_allclose(ratio[finite], math.sqrt(4.0), rtol=1e-8)


class TestMRMode:
    def test_snl_column_named_snr_m_in_mr_mode(self, tmp_path):
        params = _base_params(
            tmp_path, mr_mode=True, outfile_snc=None, outfile_oii=None
        )
        results = engine.run_etc(params)
        assert "snr_m" in results.snl.colnames
        assert "snr_r" not in results.snl.colnames

    def test_noise_arm_ids_remap_red_arm_to_3(self, tmp_path):
        params = _base_params(
            tmp_path, mr_mode=True, outfile_snc=None, outfile_snl=None, outfile_oii=None
        )
        results = engine.run_etc(params)
        assert set(np.unique(results.noise["arm"])) == {0, 2, 3}


class TestOiiNgtotQuirk:
    def test_detected_object_outside_histogram_range_not_counted(self, tmp_path):
        # gsetc.c:2160-2165 QUIRK: `ngtot++` (and `ngal[j]++`) sit *inside*
        # the `if (j>=0 && j<NZ_OII)` bin-range check, so a detected object
        # with z outside [0.1, 2.5) is written to the output catalog but
        # counted in neither the histogram nor ngtot. Reachable with the
        # packaged LR config: at z=0.05 the doublet sits at ~391nm, inside
        # the blue arm's SNR gate (373.8*1.05 = 392.5 > lmin = 380), yet
        # bins to floor((0.05-0.1)/0.1) = -1.
        #
        # The catalog also deliberately mixes rows sharing one r_eff (rows
        # 1/2: a multi-row vectorized group within one arm gate -- the case
        # that exposed the array-r_eff broadcast bug against snr.snr_oii's
        # scalar-r_eff contract) with a row at a different r_eff (row 3,
        # exercising the unique-r_eff grouping loop).
        cat = tmp_path / "oii_cat_lowz.txt"
        cat.write_text(
            "1 0.05 0.30 1.0 5.0e-15 0.0 70\n"  # detected, z below ZMIN_OII
            "2 0.50 0.30 1.0 5.0e-15 0.0 70\n"  # detected, in histogram range
            "3 0.45 0.15 1.0 5.0e-15 0.0 70\n"  # detected, different r_eff
        )
        params = _base_params(
            tmp_path,
            outfile_snc=None,
            outfile_snl=None,
            outfile_oii=None,
            oii_cat_in=cat,
        )
        results = engine.run_etc(params)

        # All three objects pass the min_snr cut and land in the catalog...
        assert set(results.oii_catalog["obj_id"]) == {1, 2, 3}
        # ...but only the in-range ones are counted (ngtot == histogram
        # sum, strictly smaller than the catalog length here).
        assert results.oii_n_targets == 2
        assert results.oii_histogram.sum() == 2
        assert results.oii_n_targets < len(results.oii_catalog)
        # And the in-range objects landed in their bins:
        # z=0.45 -> bin 3, z=0.5 -> bin 4.
        assert results.oii_histogram[3] == 1
        assert results.oii_histogram[4] == 1


class TestMapMasked:
    def test_only_masked_positions_evaluated(self):
        n = 6
        mask = np.array([False, True, False, True, True, False])
        seen = []

        def fn(chunk):
            seen.append(chunk.copy())
            return chunk * 2.0

        out = engine._map_masked(n, mask, 2, fn, np.arange(n, dtype=np.float64))
        expected = np.zeros(n)
        expected[mask] = np.arange(n, dtype=np.float64)[mask] * 2.0
        np.testing.assert_array_equal(out, expected)
        # Only masked positions were ever passed to fn.
        seen_concat = np.concatenate(seen) if seen else np.array([])
        assert set(seen_concat.tolist()) == set(np.arange(n)[mask].tolist())

    def test_empty_mask_returns_all_zero(self):
        out = engine._map_masked(
            5, np.zeros(5, dtype=bool), 2, lambda c: c, np.arange(5, dtype=np.float64)
        )
        np.testing.assert_array_equal(out, np.zeros(5))


class TestArmParallelism:
    @pytest.mark.parametrize("n_workers", [1, 2, 8])
    def test_serial_and_parallel_are_bit_identical(self, monkeypatch, n_workers):
        # `base` below sets outfile_oii on top of reference_params' default
        # outfile_snc/outfile_snl, so all three of oii_curve/snl/snc are
        # requested -- with n_workers>1 this exercises every parallelism
        # level at once: run_etc's run_products (product-level, see
        # `_parallel.py`) submitting all three concurrently, each of which
        # in turn uses map_arms (arm-level) and, inside that,
        # engine._map_masked / psf.spectro_mtf's row chunks (chunk-level)
        # internally. n_workers=8 (the new default's cap) additionally
        # exercises the arm x chunk oversubscription case.
        #
        # Small, fixed z grids (module-local, distinct from the
        # module-autouse `_tiny_z_grids` fixture's `_TINY_Z`, per the task
        # brief) -- just enough points to exercise every column-assembly
        # code path without the full-resolution sweep cost.
        tiny_z = np.arange(0.4, 0.5, 0.01)
        monkeypatch.setattr(engine, "_oii_curve_z_grid", lambda: tiny_z)
        monkeypatch.setattr(engine, "_snl_z_grid", lambda: tiny_z)

        base = reference_params(
            instr_config=CONFIG_FIXTURE,
            mag=None,
            mag_file=MAG_FILE_FIXTURE,
            reff=0.3,
            outfile_oii=Path("oii_curve.ecsv"),
        )

        serial = engine.run_etc(dataclasses.replace(base, n_workers=1))
        parallel = engine.run_etc(dataclasses.replace(base, n_workers=n_workers))

        for name in ("noise", "snc", "snl", "oii_curve"):
            serial_table = getattr(serial, name)
            parallel_table = getattr(parallel, name)
            assert serial_table.colnames == parallel_table.colnames
            for col in serial_table.colnames:
                # Bit-identical, not merely close: arm-level and
                # chunk-level parallelism must not perturb floating-point
                # aggregation order (see `_parallel.py`'s module
                # docstring).
                np.testing.assert_array_equal(
                    np.asarray(serial_table[col]),
                    np.asarray(parallel_table[col]),
                    err_msg=f"table={name!r} column={col!r} n_workers={n_workers}",
                )
