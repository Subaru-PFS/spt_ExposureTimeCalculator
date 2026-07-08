"""Slow regression gate for pfsspecsim.etc.noise (task T8).

Compares `compute_noise` against the C engine's frozen reference output
`tests/master_results/noise_omp.dat` (columns: output arm id, pixel,
lambda[nm], Noise variance, SkyMod), produced by the OpenMP C binary with
the stdin parameters recorded in `tests/gsetc_params.txt`. Both files are
protected fixtures -- never regenerate or edit them.

Acceptance criterion (task brief): rtol=1.5e-3 agreement on >= 99.9% of
pixels, per column (col 3 = noise variance, col 4 = sky), across all arms.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from pfsspecsim.etc.config import load_spectrograph_config
from pfsspecsim.etc.noise import compute_noise
from pfsspecsim.etc.params import EtcParams

TESTS_DIR = Path(__file__).resolve().parents[1]
# Protected fixtures; never modify (see task brief).
CONFIG_FIXTURE = TESTS_DIR / "PFS.20211220.dat"
REFERENCE_FIXTURE = TESTS_DIR / "master_results" / "noise_omp.dat"

RTOL = 1.5e-3
MIN_PASS_FRACTION = 0.999


def _reference_params() -> EtcParams:
    """Parameters matching tests/gsetc_params.txt (see task brief).

    The C run used config PFS.20211220.dat, degrade=1.0, sky_type='11005',
    seeing=0.8, zenith_ang=45.0, galactic_ext=0.0, field_ang=0.675,
    fiber_offset=0.0, moon=(30, 60, 0.0), exp_time=900, exp_num=8,
    sky_sub_floor=0.01, diffuse_stray=0.02. obsc_fov_dep=False because the
    C engine has no calc_obscuration correction (that was applied by the
    old Python wrapper, not by gsetc.x, and the fixture was produced by
    driving gsetc.x directly).
    """
    return EtcParams(
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


@pytest.mark.slow
def test_noise_matches_c_reference():
    ref = np.loadtxt(REFERENCE_FIXTURE)
    arm_col = ref[:, 0].astype(int)

    spectro = load_spectrograph_config(CONFIG_FIXTURE, degrade=1.0)
    assert spectro.MR is False  # LR config: output arm id == internal ia.
    result = compute_noise(_reference_params(), spectro)

    failures = []
    for ia in range(spectro.N_arms):
        mask = arm_col == ia
        assert mask.sum() == spectro.npix[ia]
        arm = result[ia]

        np.testing.assert_allclose(arm.lam, ref[mask, 2], atol=5.5e-5)

        for name, ours, theirs in (
            ("noise", arm.noise, ref[mask, 3]),
            ("sky", arm.sky, ref[mask, 4]),
        ):
            rel = np.abs(ours - theirs) / np.abs(theirs)
            pass_fraction = np.mean(rel <= RTOL)
            if pass_fraction < MIN_PASS_FRACTION:
                failures.append(
                    f"arm {ia} {name}: only {pass_fraction:.5%} of pixels "
                    f"within rtol={RTOL} (need >= {MIN_PASS_FRACTION:.1%}); "
                    f"max rel deviation {rel.max():.3e} at pixel "
                    f"{int(rel.argmax())} "
                    f"(lambda={arm.lam[rel.argmax()]:.3f} nm, "
                    f"ours={ours[rel.argmax()]:.6e}, "
                    f"ref={theirs[rel.argmax()]:.6e})"
                )
            else:
                # Keep the achieved margin visible in verbose runs.
                print(
                    f"arm {ia} {name}: max rel deviation {rel.max():.3e} "
                    f"(gate {RTOL}), pass fraction {pass_fraction:.5%}"
                )

    assert not failures, "\n".join(failures)
