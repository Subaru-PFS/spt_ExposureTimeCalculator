"""Verify astropy-derived constants match the C literals (task T2).

The C source hard-codes these as `#define`s in `src/gsetc.c`; `constants.py`
re-derives them from `astropy.constants` / `astropy.units` so that the
package no longer needs to hand-type unit-conversion factors. This test
locks down that the two never drift apart by more than the ~1e-3 physical
regression tolerance used elsewhere in the port (here we hold constants to
a much tighter 1e-6 relative tolerance, since these are exact unit
conversions / well-known physical constants, not numerical approximations).
"""

import math

from pfsspecsim.etc import constants


def test_photons_per_erg_1nm():
    # C literal: PHOTONS_PER_ERG_1NM = 5.03411747e8 (gsetc.c:19)
    assert math.isclose(
        constants.PHOTONS_PER_ERG_1NM, 5.03411747e8, rel_tol=1e-6
    )


def test_ab_zeropoint_cgs():
    # C literal: 3.631e-20 erg/s/cm2/Hz (gsetc.c:1221, 1428)
    assert math.isclose(constants.AB_ZEROPOINT_CGS, 3.631e-20, rel_tol=1e-6)


def test_arcsec_per_urad():
    # C literal: ARCSEC_PER_URAD = 0.206264806247097 (gsetc.c:17)
    assert math.isclose(
        constants.ARCSEC_PER_URAD, 0.206264806247097, rel_tol=1e-6
    )


def test_c_nm_per_s():
    # C literal: 2.99792458e17 nm/s (e.g. gsetc.c:1231)
    assert math.isclose(constants.C_NM_PER_S, 2.99792458e17, rel_tol=1e-6)


def test_c_km_per_s():
    assert math.isclose(constants.C_KM_PER_S, 299792.458, rel_tol=1e-9)


def test_deg_to_rad():
    assert math.isclose(constants.DEG_TO_RAD, math.pi / 180.0, rel_tol=1e-12)


def test_algorithm_literals_unchanged():
    # Not physical constants -- must be kept verbatim from gsetc.c.
    assert constants.RAT_HL_SL_EXP == 1.67834  # gsetc.c:22
    assert constants.SP_PSF_LEN == 32  # gsetc.c:34
    assert constants.NP_WIN == 32  # gsetc.c:1076
    assert constants.REF_SIZE == 0.30  # gsetc.c:2025
    assert (constants.ZMIN_OII, constants.NZ_OII, constants.DZ_OII) == (
        0.10,
        24,
        0.10,
    )  # gsetc.c:9-11
    assert constants.ADJUST_NOISE_LR == (0.752, 0.722, 2.694)  # gsetc.c:722-723
    assert constants.ADJUST_NOISE_MR == (0.752, 1.015, 2.694)  # gsetc.c:722-723
    assert constants.ADJUST_THROUGHPUT_LR == (1.08, 1.05, 1.24)  # gsetc.c:1458-1459
    assert constants.ADJUST_THROUGHPUT_MR == (1.08, 1.20, 1.24)  # gsetc.c:1458-1459
    assert constants.OII_LAMBDA == (372.71, 372.98)  # gsetc.c:1286-1287
    assert constants.L_CHUNK == 4096
    assert constants.Z_CHUNK == 2048
