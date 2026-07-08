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

import numpy as np
import pytest

from pfsspecsim.etc.config import load_spectrograph_config
from pfsspecsim.etc.noise import compute_noise

from conftest import CONFIG_FIXTURE, MASTER_RESULTS_DIR, reference_params

REFERENCE_FIXTURE = MASTER_RESULTS_DIR / "noise_omp.dat"

RTOL = 1.5e-3
MIN_PASS_FRACTION = 0.999


@pytest.mark.slow
def test_noise_matches_c_reference():
    ref = np.loadtxt(REFERENCE_FIXTURE)
    arm_col = ref[:, 0].astype(int)

    spectro = load_spectrograph_config(CONFIG_FIXTURE, degrade=1.0)
    assert spectro.MR is False  # LR config: output arm id == internal ia.
    result = compute_noise(reference_params(), spectro)

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
