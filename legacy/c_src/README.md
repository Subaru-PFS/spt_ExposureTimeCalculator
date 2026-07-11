# legacy/c_src -- original C ETC source (frozen reference)

This directory vendors the original C exposure-time-calculator program that
`pfsspecsim.etc` (`src/pfsspecsim/etc/`) is a pure-Python port of. It is kept
as a permanent, frozen reference; it is **not** built or imported by anything
in this repository.

Files:

- `gsetc.c` -- the ETC engine. Header comments identify it as `VERSION 5`,
  `MODIFICATION 6`: Christopher M. Hirata's original `gsetc`, with
  PFS-specific modifications by Y. Moritani, K. Yabe, and others (dated
  inline comments throughout the file, e.g. "Added by Y.Moritani for input
  mag. file: 20150422", "Added by K.Yabe for input mag. file: 20160205").
- `modeldata.h` -- the data header `gsetc.c` includes: UVES/OH sky line
  atlases, Kitt Peak / Mauna Kea atmospheric transmission grids, and other
  lookup tables baked in as C array literals.
- `Makefile` -- builds `gsetc.x` from `gsetc.c`. The reference binary used to
  produce the acceptance fixtures (see below) was compiled with
  `CFLAGS = -O3 -DHGCDTE_SUTR -DMOONLIGHT_`, exactly as this Makefile
  specifies.

## Why this is kept

This is the reference implementation for the pure-Python port living at
`src/pfsspecsim/etc/`. The C-reference acceptance fixtures under
`tests/master_results/` (`noise_omp.dat`, `snc_omp.dat`, `snl_omp.dat`,
`sno2_omp.dat`) were produced by this exact code, compiled with the flags
above, driven by the stdin parameters in `tests/gsetc_params.txt`. The slow
regression gates (`tests/python/test_noise_reference.py`,
`tests/python/test_reference_outputs.py`) compare the Python engine's output
against those fixtures and are the arbiter of numeric correctness for the
port (see the repository's `CLAUDE.md`, "Non-negotiable invariants"). The
underlying algorithm and physics are documented in `docs/Manual_v5.pdf`
(Hirata 2012).

## Line-number contract -- do not edit this copy

Docstrings and comments throughout `src/pfsspecsim/etc/` cite `gsetc.c:<line>`
to point at the exact C statement a given Python function or constant was
ported from. Some older docstrings spell the path as `src/gsetc.c` -- that
was this file's location before the pure-Python port removed the C sources
from `src/` (see `docs/review-c-to-python-port-2026-07-10.md`); it is the
same file, byte-identical, now living here instead.

Every one of those line-number citations resolves against **this exact
copy** of `gsetc.c`. Because of that, this copy must never be edited,
reformatted, or otherwise modified -- treat it exactly like the fixtures
under `tests/master_results/`: a frozen input, not a file to improve. If the
upstream C source is ever revisited, do so in a separate location and leave
this copy untouched.

## Known bugs and quirks in this C source

The 2026-07 port review (`docs/review-c-to-python-port-2026-07-10.md`)
catalogued the following. They are documented here so that future readers of
`gsetc.c` do not mistake them for unknowns -- but per the line-number
contract above, **none of them may be fixed in this copy**.

### Intentional quirks, preserved verbatim in the Python port

Each of these is reproduced exactly by `src/pfsspecsim/etc/` and pinned by a
dedicated regression test; "fixing" any of them on the Python side breaks
the acceptance gates.

- `gsetc.c:1213-1218` -- `gsGetSNR_Single`'s 41-point Gaussian transmission
  average evaluates `gsAtmTrans` at the *same, unshifted* wavelength at
  every quadrature point (unlike `gsGetSignal`'s genuinely velocity-smeared
  average at 1100-1105). Mathematically equal to the plain pointwise
  transmission, which is how the port implements it.
- `gsetc.c:2049` vs `2084` -- the [OII] curve's diagnostic aperture-factor
  column calls `gsGeometricThroughput` with `fieldang=0`, while the
  single-line curve uses the real field angle.
- `gsetc.c:1401` vs `882` -- `gsGetSNR_Continuum` evaluates each pixel's
  wavelength at the pixel *left edge* (`lmin + dl*ipix`), whereas
  `gsGetNoise`'s sky-continuum grid uses the pixel *center*
  (`lmin + (ipix+0.5)*dl`).
- `gsetc.c:805-813` vs `842-852` -- the UVES sky lines are rescaled to a
  reference airmass of 1.1, the OH lines to 1.0 (plus an
  `exp((14.8-15.8)/1.086)` brightness-level factor).
- `gsetc.c:2160-2165` -- `ngtot++` (and `ngal[j]++`) sit *inside* the
  `j>=0 && j<NZ_OII` histogram-bin check, so a detected object with z
  outside [0.1, 2.5) is written to the output catalog but counted in
  neither the histogram nor `ngtot`.
- `gsetc.c:1918-1947` -- the magnitude-file padding reuses the *raw*,
  unthresholded first/last data values when extending past the data range
  (the `mag <= 0 -> 99.9` substitution is not re-applied to those points).
- `gsetc.c:813, 851` -- the magnitude-to-nepers conversion hardcodes the
  rounded `1.086` rather than the exact `2.5/ln(10) = 1.0857...`.

### Latent bug, deliberately NOT replicated by the Python port

- `gsetc.c:1980-2000` (`flag_reused==1` noise-reload branch) -- rows are
  assigned to arms by fixed offsets against a flat row counter `k`, with
  `else if (k>=spectro.npix[1])` where `npix[0]` was evidently intended
  (the index used is `k-spectro.npix[0]`). This only works because arm 0
  and arm 1 happen to share a pixel count in the PFS configs; any config
  where they differ mis-assigns rows. The port
  (`engine._noise_arrays_from_table`) instead groups rows explicitly by the
  file's own `arm` column, which is identical for all shipped configs and
  correct in general.

### Dead code (not ported)

- `geterf` (`gsetc.c:161-204`) is defined but never called.
- The moonlight block computes `k = gsAtmContOp(...)` (`gsetc.c:960`) and
  never uses it (only `kV`, evaluated at 550 nm, feeds the Krisciunas &
  Schaefer formula).
- `REF_SIZE` (`gsetc.c:2025`) is defined but unused by the live code.

### Documentation discrepancy (manual typo, code is correct)

- `docs/Manual_v5.pdf` §4.I prints the aperture perimeter-to-area ratio as
  `4(1-upsilon)/D_outer`. The code implements the diffraction scale as
  `theta_D = lambda/(D_outer*(1-centobs))` (`gsetc.c:556`), i.e.
  `sigma = 4/[D(1-upsilon)]` -- the correct ratio for an annular aperture
  (outer rim `pi*D` + obscuration rim `pi*upsilon*D` over area
  `(pi/4)*D^2*(1-upsilon^2)`). The manual's version is also unphysical: it
  would make diffraction losses *decrease* as the obscuration grows,
  whereas a thinner annulus has more edge per unit area and diffracts
  more. Confirmed against the printed PDF 2026-07-10 (the exponent
  `exp(-(4/pi)*u*theta_D)` does match the manual): **the manual's (1-upsilon)
  placement is a typo; this code is right.** Effect on the aperture factor
  is sub-percent for PFS parameters either way.

Not a bug in this C code, but related: the pre-2.0 *Python wrapper* had a
degrade-compounding bug (the field-angle obscuration correction was
multiplied into the stored `degrade` parameter on every `run()`/`make_*()`
call), documented in detail in
`docs/review-c-to-python-port-2026-07-10.md` §6.1. The wrapper drove this C
binary but the bug lived entirely in the Python layer.

## `modeldata.h` provenance

`src/pfsspecsim/etc/data/modeldata.npz` was extracted from this `modeldata.h`
(and from the inline tables in this `gsetc.c`) by `tools/extract_modeldata.py`
-- see `src/pfsspecsim/etc/data/README.md` for the full table-by-table
provenance. Regenerate `modeldata.npz` only by re-running that script against
this directory; do not hand-edit the archive or this header.
