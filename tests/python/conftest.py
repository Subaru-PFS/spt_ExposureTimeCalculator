"""Shared fixtures/helpers for the `tests/python` suite.

`reference_params` builds the `EtcParams` matching the C+OpenMP reference
run recorded in `tests/gsetc_params.txt` (config `tests/PFS.20211220.dat`);
`test_noise_reference.py` (T8) and `test_reference_outputs.py` (T11) both
call it rather than duplicating the field-by-field construction, tweaking
only what each gate actually needs (T8 runs `compute_noise` directly with a
scalar `mag`; T11 runs the full `run_etc` pipeline and needs `mag_file`,
`reff`, `line_flux`, `line_width`, `outfile_oii` on top).
"""

from __future__ import annotations

from pathlib import Path

from pfsspecsim.etc.params import EtcParams

#: Protected fixtures (see the task briefs); never modify.
TESTS_DIR = Path(__file__).resolve().parents[1]
CONFIG_FIXTURE = TESTS_DIR / "PFS.20211220.dat"
MAG_FILE_FIXTURE = TESTS_DIR / "mag_18.dat"
MASTER_RESULTS_DIR = TESTS_DIR / "master_results"


def reference_params(**overrides) -> EtcParams:
    """`EtcParams` matching the stdin values in `tests/gsetc_params.txt`.

    Common ground for both the noise-only gate (T8) and the full-pipeline
    gate (T11): config `tests/PFS.20211220.dat` (via `instr_config`,
    supplied separately -- not an `EtcParams` field), degrade=1.0,
    sky_type='11005', seeing=0.8, zenith_ang=45.0, galactic_ext=0.0,
    field_ang=0.675, fiber_offset=0.0, moon=(30, 60, 0.0), exp_time=900,
    exp_num=8, sky_sub_floor=0.01, diffuse_stray=0.02. `obsc_fov_dep=False`
    because the C engine has no `calc_obscuration` correction (that was
    applied by the old Python wrapper, not by gsetc.x, and the fixture was
    produced by driving gsetc.x directly).

    Callers must still set `instr_config` (both gates point it at
    `CONFIG_FIXTURE`) and whichever of `mag`/`mag_file` they need.
    """
    values = dict(
        seeing=0.8,
        zenith_ang=45.0,
        galactic_ext=0.0,
        field_ang=0.675,
        fiber_offset=0.0,
        moon_zenith_ang=30.0,
        moon_target_ang=60.0,
        moon_phase=0.0,
        exp_time=900.0,
        exp_num=8,
        mag=22.5,
        sky_type="11005",
        sky_sub_floor=0.01,
        diffuse_stray=0.02,
        degrade=1.0,
        obsc_fov_dep=False,
        hgcdte_sutr=True,
    )
    values.update(overrides)
    return EtcParams(**values)
