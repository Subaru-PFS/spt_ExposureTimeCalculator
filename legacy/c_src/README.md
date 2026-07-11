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

## `modeldata.h` provenance

`src/pfsspecsim/etc/data/modeldata.npz` was extracted from this `modeldata.h`
(and from the inline tables in this `gsetc.c`) by `tools/extract_modeldata.py`
-- see `src/pfsspecsim/etc/data/README.md` for the full table-by-table
provenance. Regenerate `modeldata.npz` only by re-running that script against
this directory; do not hand-edit the archive or this header.
