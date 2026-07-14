"""Full-pipeline acceptance gate against the C reference (task T11).

Runs `engine.run_etc` once per frozen reference set (module-scoped fixture,
parametrized over `conftest.REFERENCE_SETS` -- see `tests/README.md`) with
every output enabled -- noise, continuum SNR ("snc"), single-line SNR
("snl"), and the [OII]-doublet SNR-vs-redshift curve ("sno2") -- for the
same observing conditions as T8's noise-only gate (`conftest.reference_params`)
plus the `tests/gsetc_params.txt` fields T8 didn't need:
`mag_file=tests/mag_18.dat` (`mag=None`), `reff=0.3`, `line_flux=1e-17`,
`line_width=70` (`sigma_v` km/s). The four reference sets are the original
LR/sky_type=11005 run (config `tests/PFS.20211220.dat`) plus MR/11005,
LR/11006, and MR/11006 siblings (configs `tests/PFS.redMR.20211220.dat`,
`tests/PFS.20240714.dat`, `tests/PFS.redMR.20240714.dat`). The [OII]
*catalog* output is not exercised here: `InFileOII='-'` in every
`gsetc_params*.txt`, so no C reference run ever produced one, and there is
no `tests/master_results/*.dat` to compare against.

Each table is compared column-by-column against its frozen
`tests/master_results/{noise,snc,snl,sno2}_<suffix>.dat` (protected
fixtures -- never regenerate or edit them) under a combined tolerance,
``|python - C| <= atol + rtol * |C|`` (`rtol=1.5e-3` throughout, matching
T8), with `atol` derived from each column's C `printf` precision
(gsetc.c:2010 [noise], 2109 [snc], 2084 [snl], 2048 [sno2]): half the last
printed decimal place for fixed-point columns (e.g. a `%.4lf` column gets
`atol=5e-5`), or half the last mantissa digit *scaled by the printed
value's own magnitude* for `%...le` scientific-notation columns (a 5-digit
mantissa is accurate to 1 part in 2e5). The `atol` term matters only for
columns that legitimately print as (near-)zero -- e.g. `snr` far from any
line/continuum flux -- where a bare relative comparison would be
meaningless; every column passes at 100% of rows (>= the required 99.9%)
well before `atol` needs to do any work for the bulk of "real signal" rows.

Achieved max relative deviation per table, as measured during development
against the original `lr_11005` reference set (on rows whose C value is not
itself swept into the `atol` floor, i.e. `|C| > 0.05` for the SNR columns --
below that, the print-quantization `atol` dominates and a plain relative
deviation stops being meaningful): noise 1.4e-05, snc 3.2e-04, snl 4.3e-04,
sno2 7.3e-04 (all comfortably under the `rtol=1.5e-3` gate; the `snX`/`snoY`
curves' SNR columns carry the largest deviations of any table, consistent
with T8's own noise-only margin of 1.4e-05 -- the extra ~20x is the
`snr.py` SNR-formula chain built on top of that noise vector, not the noise
vector itself). See CLAUDE.md for the achieved max deviations of the other
three reference sets (`mr_11005`, `lr_11006`, `mr_11006`). Every column
below is currently at 100% of rows within tolerance (the print loop -- run
with `-s`/`-v` -- shows the live max/pass-fraction per column, per
reference set, if that ever regresses).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pfsspecsim.etc import engine
from pfsspecsim.etc.config import load_spectrograph_config, spectro_arm

from conftest import (
    MAG_FILE_FIXTURE,
    MASTER_RESULTS_DIR,
    REFERENCE_SETS,
    reference_params,
)

RTOL = 1.5e-3
MIN_PASS_FRACTION = 0.999


def _atol_fixed(decimals: int) -> float:
    """Half the last printed decimal place of a fixed-point `%.<decimals>f`
    column -- the C engine's own print quantization."""
    return 0.5 * 10 ** (-decimals)


def _atol_sci(theirs: np.ndarray, decimals: int = 5) -> np.ndarray:
    """Print quantization for a `%.<decimals>le`/`%.<decimals>lE`
    scientific-notation column: the mantissa carries `decimals` digits
    after the point, so the absolute quantization scales with the printed
    value's own magnitude."""
    return 0.5 * 10 ** (-decimals) * np.abs(theirs)


def _check_column(
    ours: np.ndarray, theirs: np.ndarray, atol
) -> tuple[float, np.ndarray, np.ndarray]:
    """Combined-tolerance pass fraction (``|a-b| <= atol + rtol*|b|``) plus
    a plain relative deviation (using `atol` as a zero-guard denominator
    floor) for worst-offender reporting.
    """
    ours = np.asarray(ours, dtype=np.float64)
    theirs = np.asarray(theirs, dtype=np.float64)
    diff = np.abs(ours - theirs)
    atol_arr = np.broadcast_to(np.atleast_1d(atol), theirs.shape).astype(np.float64)
    limit = atol_arr + RTOL * np.abs(theirs)
    passed = diff <= limit
    denom = np.maximum(np.abs(theirs), atol_arr)
    denom = np.where(denom > 0, denom, 1.0)
    rel = diff / denom
    return float(np.mean(passed)), rel, diff


def _check_table(
    table_name: str,
    columns: list[tuple[str, np.ndarray, np.ndarray, object]],
    row_labels: np.ndarray,
    failures: list[str],
) -> None:
    """Check every `(col_name, ours, theirs, atol)` column of one table,
    appending a worst-offender message per failing column to `failures`.
    `row_labels` identifies each row for the failure message (e.g. a
    per-arm pixel index, or a z-grid index).
    """
    for col_name, ours, theirs, atol in columns:
        pass_fraction, rel, diff = _check_column(ours, theirs, atol)
        if pass_fraction < MIN_PASS_FRACTION:
            worst = int(np.argmax(rel))
            failures.append(
                f"{table_name}.{col_name}: only {pass_fraction:.5%} of rows "
                f"within tolerance (need >= {MIN_PASS_FRACTION:.1%}); worst "
                f"offender row={row_labels[worst]}, "
                f"python={np.asarray(ours, dtype=np.float64)[worst]:.6e}, "
                f"C={np.asarray(theirs, dtype=np.float64)[worst]:.6e}, "
                f"rel_dev={rel[worst]:.3e}"
            )
        else:
            print(
                f"{table_name}.{col_name}: max rel deviation {rel.max():.3e}, "
                f"pass fraction {pass_fraction:.5%}"
            )


@pytest.fixture(scope="module", params=REFERENCE_SETS, ids=lambda rs: rs.id)
def reference_results(request):
    """Run the full pipeline once per reference set, with every
    C-comparable output enabled.

    `outfile_oii` (default `None`) is set to enable the [OII]-curve
    ("sno2") computation; `outfile_noise`/`outfile_snc`/`outfile_snl`
    already default to non-`None` paths. `run_etc` (not `run_etc_files`)
    is used, so nothing is written to disk -- the output paths only gate
    which tables get computed. Parametrized over `conftest.REFERENCE_SETS`
    (see `tests/README.md`): LR/11005, MR/11005, LR/11006, MR/11006.
    """
    ref_set = request.param
    params = reference_params(
        mag=None,
        mag_file=MAG_FILE_FIXTURE,
        reff=0.3,
        line_flux=1.0e-17,
        line_width=70.0,
        min_snr=9.0,
        sky_type=ref_set.sky_type,
        instr_config=ref_set.config,
        outfile_oii=Path("reference_oii_curve.ecsv"),
    )
    spectro = load_spectrograph_config(ref_set.config, degrade=1.0)
    results = engine.run_etc(params)
    return results, spectro, ref_set


@pytest.mark.slow
def test_noise_table_matches_c_reference(reference_results):
    results, spectro, ref_set = reference_results
    ref = np.loadtxt(MASTER_RESULTS_DIR / f"noise_{ref_set.suffix}.dat")
    arm_col = ref[:, 0].astype(int)
    table_arm = np.asarray(results.noise["arm"])
    table_pixel = np.asarray(results.noise["pixel"])

    failures: list[str] = []
    for ia in range(spectro.N_arms):
        arm_id = spectro_arm(ia, spectro.MR)
        mask = arm_col == arm_id
        sel = table_arm == arm_id
        assert mask.sum() == spectro.npix[ia]
        assert sel.sum() == spectro.npix[ia]

        columns = [
            (
                "wavelength",
                np.asarray(results.noise["wavelength"])[sel],
                ref[mask, 2],
                _atol_fixed(4),
            ),
            (
                "variance",
                np.asarray(results.noise["variance"])[sel],
                ref[mask, 3],
                _atol_sci(ref[mask, 3]),
            ),
            (
                "sky",
                np.asarray(results.noise["sky"])[sel],
                ref[mask, 4],
                _atol_sci(ref[mask, 4]),
            ),
        ]
        _check_table(f"noise[arm={ia}]", columns, table_pixel[sel], failures)

    assert not failures, "\n".join(failures)


@pytest.mark.slow
def test_snc_table_matches_c_reference(reference_results):
    results, spectro, ref_set = reference_results
    assert results.snc is not None
    ref = np.loadtxt(MASTER_RESULTS_DIR / f"snc_{ref_set.suffix}.dat")
    arm_col = ref[:, 0].astype(int)
    table_arm = np.asarray(results.snc["arm"])
    table_pixel = np.asarray(results.snc["pixel"])

    # Column order mirrors the C `fprintf` (gsetc.c:2109-2110) and
    # `engine._compute_snc`'s Table column order (both 11 columns: arm,
    # pixel, then these 9).
    fixed_point_cols = [
        ("wavelength", 2, 3),
        ("snr", 3, 4),
    ]
    sci_cols = [
        ("signal", 4),
        ("noise_variance", 5),
        ("noise_variance_tot", 6),
        ("input_mag", 7),
        ("conversion_factor", 8),
        ("sampling_factor", 9),
        ("sky", 10),
    ]

    failures: list[str] = []
    for ia in range(spectro.N_arms):
        arm_id = spectro_arm(ia, spectro.MR)
        mask = arm_col == arm_id
        sel = table_arm == arm_id
        assert mask.sum() == spectro.npix[ia]
        assert sel.sum() == spectro.npix[ia]

        columns = [
            (name, np.asarray(results.snc[name])[sel], ref[mask, col], _atol_fixed(dec))
            for name, col, dec in fixed_point_cols
        ] + [
            (
                name,
                np.asarray(results.snc[name])[sel],
                ref[mask, col],
                _atol_sci(ref[mask, col]),
            )
            for name, col in sci_cols
        ]
        _check_table(f"snc[arm={ia}]", columns, table_pixel[sel], failures)

    assert not failures, "\n".join(failures)


@pytest.mark.slow
def test_snl_table_matches_c_reference(reference_results):
    results, spectro, ref_set = reference_results
    assert results.snl is not None
    ref = np.loadtxt(MASTER_RESULTS_DIR / f"snl_{ref_set.suffix}.dat")
    assert len(results.snl) == ref.shape[0]
    row_idx = np.arange(ref.shape[0])

    mid_name = "snr_m" if spectro.MR else "snr_r"
    # Column order mirrors the C `fprintf` (gsetc.c:2083-2084) and
    # `engine._compute_snl`'s Table column order (7 columns total).
    columns = [
        ("wavelength", ref[:, 0], _atol_fixed(2)),
        ("fiber_aperture_factor", ref[:, 1], _atol_fixed(6)),
        ("effective_area", ref[:, 2], _atol_fixed(5)),
        ("snr_b", ref[:, 3], _atol_fixed(4)),
        (mid_name, ref[:, 4], _atol_fixed(4)),
        ("snr_n", ref[:, 5], _atol_fixed(4)),
        ("snr_tot", ref[:, 6], _atol_fixed(4)),
    ]
    columns = [
        (name, np.asarray(results.snl[name]), theirs, atol)
        for name, theirs, atol in columns
    ]

    failures: list[str] = []
    _check_table("snl", columns, row_idx, failures)
    assert not failures, "\n".join(failures)


@pytest.mark.slow
def test_oii_curve_table_matches_c_reference(reference_results):
    results, spectro, ref_set = reference_results
    assert results.oii_curve is not None
    ref = np.loadtxt(MASTER_RESULTS_DIR / f"sno2_{ref_set.suffix}.dat")
    assert len(results.oii_curve) == ref.shape[0]
    row_idx = np.arange(ref.shape[0])

    # Column order mirrors the C `fprintf` (gsetc.c:2047-2050) and
    # `engine._compute_oii_curve`'s Table column order (9 columns total).
    columns = [
        ("z", ref[:, 0], _atol_fixed(4)),
        ("wavelength0", ref[:, 1], _atol_fixed(2)),
        ("wavelength1", ref[:, 2], _atol_fixed(2)),
        ("fiber_aperture_factor", ref[:, 3], _atol_fixed(6)),
        ("effective_area", ref[:, 4], _atol_fixed(5)),
        ("snr_b", ref[:, 5], _atol_fixed(4)),
        ("snr_r", ref[:, 6], _atol_fixed(4)),
        ("snr_n", ref[:, 7], _atol_fixed(4)),
        ("snr_tot", ref[:, 8], _atol_fixed(4)),
    ]
    columns = [
        (name, np.asarray(results.oii_curve[name]), theirs, atol)
        for name, theirs, atol in columns
    ]

    failures: list[str] = []
    _check_table("oii_curve", columns, row_idx, failures)
    assert not failures, "\n".join(failures)
