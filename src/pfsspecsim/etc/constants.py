"""Physical constants, unit conversions, and algorithm-specific literals.

Physical constants and unit conversions are derived once, at import time,
from ``astropy.constants`` / ``astropy.units`` and exposed as plain
module-level ``float``s. Hot loops elsewhere in the ``etc`` package must
never hold on to ``astropy.units.Quantity`` objects (too slow); they use
these plain floats with bare ``numpy`` arrays instead. ``Quantity`` usage is
confined to this module and to column-unit tagging in ``io.py``.

Do not add hand-typed conversion factors (``0.206264...``, ``5.034e8``,
``pi/180``, ...) anywhere else in the package -- derive them here from
astropy and import the resulting float.

Values that are *not* physical constants (empirical fit coefficients,
window/chunk sizes, catalog-defining literals, ...) are kept verbatim from
the C source `src/gsetc.c`, with a line-number reference, since they carry
no derivable physical meaning.
"""

from __future__ import annotations

from astropy import constants as const
from astropy import units as u

# --- Physical constants / unit conversions (astropy-derived) --------------

#: Photons per erg for a 1 nm photon (E_photon = h*c/lambda).
#: Matches the C literal ``PHOTONS_PER_ERG_1NM = 5.03411747e8`` (gsetc.c:19).
PHOTONS_PER_ERG_1NM = (1.0 * u.nm / (const.h * const.c)).to_value(1 / u.erg)

#: arcsec per microradian.
#: Matches the C literal ``ARCSEC_PER_URAD = 0.206264806247097`` (gsetc.c:17).
ARCSEC_PER_URAD = (1.0 * u.urad).to_value(u.arcsec)

#: Speed of light in nm/s.
#: Matches the C literal ``2.99792458e17`` used e.g. at gsetc.c:1231.
C_NM_PER_S = const.c.to_value(u.nm / u.s)

#: Speed of light in km/s; used for the line_width (sigma_v) -> sigma_pix
#: conversion (gsetc.c:1115).
C_KM_PER_S = const.c.to_value(u.km / u.s)

#: Degrees -> radians conversion factor (replaces hand-typed ``pi/180``).
DEG_TO_RAD = u.deg.to(u.rad)

#: AB-magnitude zero-point flux density, in erg/s/cm2/Hz.
#:
#: The C source uses the literal ``3.631e-20`` (gsetc.c:1221, 1428), which is
#: the commonly quoted "3631 Jy" AB zero point. Deriving this the "obvious"
#: astropy way, ``(0.0 * u.ABmag).to_value(u.erg/u.s/u.cm**2/u.Hz)``, instead
#: yields ``3.6307805...e-20`` erg/s/cm2/Hz: astropy's ``ABmag`` equivalency
#: is defined from the Oke & Gunn (1983) zero point of -48.60 mag
#: (``10**(-48.6/2.5)`` erg/s/cm2/Hz = 3630.78 Jy), which differs from the
#: rounded "3631 Jy" convention baked into gsetc.c by a relative ~6e-5 --
#: two orders of magnitude looser than the ~1e-6 tolerance required for T2's
#: constants regression test. Anchoring on the 3631 Jy literal (itself an
#: algorithm-defining constant, not something with more "true" precision to
#: derive) and letting astropy perform the exact Jy -> cgs unit conversion
#: keeps the promise of "no hand-typed conversion factors" while matching
#: the C literal to full double precision.
AB_ZEROPOINT_CGS = (3631.0 * u.Jy).to_value(u.erg / u.s / u.cm**2 / u.Hz)

#: Photon energy conversion h*c in eV nm; ``hnu[eV] = HC_EV_NM / lambda[nm]``.
#: Matches the C literal ``1.239842e-6`` eV m used in gsOP_Si_abslength
#: (gsetc.c:262) to within ~1e-8 relative.
HC_EV_NM = (const.h * const.c).to_value(u.eV * u.nm)

#: Boltzmann constant in eV/K.
#: Matches the C literal ``k = 8.617e-5`` eV/K (gsetc.c:245) to within
#: ~4e-5 relative (well inside the ~1e-3 regression tolerance).
K_B_EV_PER_K = const.k_B.to_value(u.eV / u.K)

# --- Algorithm-specific literals (verbatim from gsetc.c; not physical) -----

#: Ratio of half-light to scale radius for an exponential-profile galaxy
#: (gsetc.c:22).
RAT_HL_SL_EXP = 1.67834

#: Spectrograph PSF length in pixels; must be even (gsetc.c:34).
SP_PSF_LEN = 32

#: Window (in pixels) used by gsGetSignal around a spectral feature
#: (gsetc.c:1076).
NP_WIN = 32

#: Fiducial effective radius (arcsec) used for the reference [OII] SNR curve
#: (gsetc.c:2025).
REF_SIZE = 0.30

#: [OII] recovery histogram grid (gsetc.c:9-11).
ZMIN_OII = 0.10
NZ_OII = 24
DZ_OII = 0.10

#: Magnitude-to-nepers conversion, ``2.5/ln(10)`` rounded to 3 decimals
#: (gsetc.c:813, 851-852: ``exp(-k*airmass/1.086)``). The C source itself
#: hardcodes the rounded ``1.086``, not the exact ``2.5/np.log(10) =
#: 1.085736...``, so this is an algorithm-defining literal (kept verbatim
#: for bit-for-bit regression fidelity) rather than an astropy-derived
#: conversion factor.
MAG_TO_NEP_APPROX = 1.086

#: Empirical noise-variance adjustment factors per spectrograph arm, indexed
#: by the internal arm index `ia` (0=Blue, 1=Red[LR/MR], 2=NIR), derived
#: from past observations (gsetc.c:722-723).
ADJUST_NOISE_LR = (0.752, 0.722, 2.694)
ADJUST_NOISE_MR = (0.752, 1.015, 2.694)

#: Empirical throughput adjustment factors per spectrograph arm, same `ia`
#: indexing as above (gsetc.c:1458-1459).
ADJUST_THROUGHPUT_LR = (1.08, 1.05, 1.24)
ADJUST_THROUGHPUT_MR = (1.08, 1.20, 1.24)

#: Rest-frame [OII] doublet wavelengths in nm (gsetc.c:1286-1287).
OII_LAMBDA = (372.71, 372.98)

#: Chunk sizes for vectorized computation over (redshift, wavelength) grids,
#: chosen to keep peak memory for any materialized (Nlines, Nu, Npix)-like
#: array bounded (see common context note on memory discipline).
L_CHUNK = 4096
Z_CHUNK = 2048
