# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`pfsspecsim` is the PFS (Prime Focus Spectrograph) Exposure Time Calculator and
spectrum simulator. As of v2.0 the ETC engine — formerly the C program
`gsetc.c`/`gsetc_omp.c` driven as a subprocess — is a **pure-Python package**
(`pfsspecsim.etc`) built on numpy/scipy/astropy. There is no C compiler,
Makefile, or OpenMP dependency. README.md covers install and end-user usage;
this file covers what you need to work *inside* the code.

## Commands

Use `uv` for everything (project targets Python **>=3.11**; `tomllib` is required).

```bash
uv sync                              # install package + dev group (pytest, black, ruff, ty)
uv run pytest tests/python -q                     # fast suite (slow gates deselected)
uv run pytest tests/python -m slow -q             # C-reference acceptance gates (~10s)
uv run pytest tests/python -q -k "test_name"      # a single test / pattern
uv run pytest tests/python/test_noise.py -q       # a single file
uv run ruff check src tests/python                # lint (config in pyproject.toml)
uv run black src tests/python                      # format
uv build                             # sdist + wheel
```

CLIs (installed as entry points): `pfs-etc` (new typer CLI), `pfs-sim-spec`
(simulator), and the deprecated `pfs-run-etc` / `pfs-gen-sim-spec` (argparse
`@file` interface over the legacy `Etc` wrapper). Example:
`uv run pfs-etc --config examples/pfs_etc_example.toml`.

Note: editing a `.py` file auto-runs `ruff --fix` + `black` on it via a
PostToolUse hook, so you don't need to format manually after edits.

## Non-negotiable invariants

**Never modify these C-reference fixtures** — they are the acceptance criteria,
produced by the original C+OpenMP binary and used by the regression gates:
`tests/master_results/`, `tests/gsetc_params.txt`, `tests/PFS.20211220.dat`,
`tests/mag_18.dat`, `tests/analyze_diff.py`. (A PreToolUse hook blocks edits to
them.)

**The slow gates are the arbiter of correctness.** `test_noise_reference.py`
(noise only) and `test_reference_outputs.py` (full pipeline) compare every
output column against the C reference at `rtol=1.5e-3` and require ≥99.9% of
values to match. Any change to a numeric kernel
(`constants/config/materials/extinction/atmosphere/psf/sky/noise/snr/engine`)
**must** keep `uv run pytest tests/python -m slow` green. Current achieved max
relative deviation: noise 1.4e-05, snc 3.2e-04, snl 4.3e-04, sno2 7.3e-04.

**Do not "fix" the ported quirks.** `pfsspecsim.etc` is a faithful port of
`gsetc.c` (deleted from the tree; readable in git history). Several intentional
C quirks/bugs are preserved verbatim, each marked with a `gsetc.c:<line>`
comment and pinned by a dedicated regression test — e.g. `snr_single`'s
fixed-λ transmission, the OII-curve aperture factor computed at `fieldang=0`,
SNC wavelength = pixel *left* edge vs noise = pixel *center*, the OH/UVES
airmass reference split (1.0 vs 1.1), and `oii_n_targets` counting only
redshift-in-range detections. If a quirk looks wrong, it is deliberate; changing
it will break a gate. Ported functions carry `gsetc.c` line references in their
docstrings — use those when cross-checking behavior against the original.

**Unit/constant discipline.** Physical constants and unit conversions are
derived from `astropy.constants`/`astropy.units` in `etc/constants.py` and
converted once to module-level plain floats; do not hand-type values like
`5.034e8`, `3.631e-20`, or `π/180` in the kernels. Algorithm-specific literals
(grid sizes, empirical fit coefficients, `ADJUST_*` factors) stay verbatim from
the C source and live in `constants.py` too.

**The `numpy<2.5` pin is intentional** and forced by `pfs.datamodel @ w.2024.06`
(uses a legacy dtype alias removed in numpy 2.5), *not* by the ETC. Don't try to
lift it here.

## Architecture

Package is **src-layout** (`src/pfsspecsim/`). The engine lives in
`src/pfsspecsim/etc/`, built in strict dependency order (each module only
imports from earlier ones):

```
constants + _modeldata (npz loader)      ← physical constants, C data tables
  → materials, extinction, atmosphere    ← Si/dust/atmosphere physics
  → config (Spectrograph dataclass)      ← instrument .dat parser
  → psf                                  ← geometric throughput, spectrograph MTF/PSF (hardest vectorization)
  → sky                                  ← sky lines, continuum, moonlight
  → noise                                ← per-arm noise/sky vectors (gsGetNoise port)
  → snr                                  ← signal + line/OII/continuum SNR
  → engine + io                          ← main() orchestrator + ECSV writer/reader
  → cli                                  ← typer app `pfs-etc`
```

Public API (`from pfsspecsim.etc import ...`): `EtcParams` (dataclass of all
inputs; `load_params(toml, overrides)` with priority CLI > TOML > defaults),
`run_etc(params) -> EtcResults` (no file I/O), `run_etc_files(params)` (runs +
writes ECSV), plus `MagSpec`, `resolve_degrade`, `calc_obscuration`.

**Data tables** (UVES/OH sky lines, atmospheric transmission, Si index, dust
norm) were extracted once from the C `modeldata.h` by `tools/extract_modeldata.py`
into `src/pfsspecsim/etc/data/modeldata.npz` (loaded via `importlib.resources`
+ `lru_cache`; `etc/data/README.md` records provenance). Regenerate only via
that script — the C source is gone from the tree.

**Outputs are Astropy ECSV** (`.ecsv`). Every table carries `table.meta` with
the resolved `EtcParams`, `etc_version`, `instr_config`, and top-level
`degrade_resolved` (the obscuration-corrected value the engine actually
applied — distinct from the raw `params.degrade` under `meta["params"]`).

**Backward-compatibility layer.** `src/pfsspecsim/pfsetc.py` keeps the legacy
`Etc` class surface (uppercase param dict, `set_param`/`run`/`make_*`/`run_multi`,
`__init__(omp_num_threads=…)` accepted+ignored with a `DeprecationWarning`). It
translates old params to `EtcParams`, drives the new engine, and restores the
old `nsm_*`/`snc_*`/`snl_*`/`sno2_*` attributes from the result tables so old
scripts and notebooks keep working. One deliberate behavioral change: the old
`degrade`-compounding bug (obscuration correction accumulating across repeated
`run()` calls) is **not** replicated, so a chained
`make_noise_model(); make_snc()` differs from pre-2.0 output by an obscuration
factor (~0.83 at `field_ang=0.45`). `src/pfsspecsim/pfsspec.py` reads the ETC's
SNC output as ECSV (with a legacy `np.loadtxt` fallback for old-format files);
its datamodel FITS output path is unchanged.

## Testing conventions

`tests/python/` mirrors the module layout (`test_<module>.py`). The `slow`
marker (see `pyproject.toml`) gates the two C-reference regression tests; the
default `-m "not slow"` run stays fast. `tests/python/conftest.py::reference_params`
builds the `EtcParams` matching `tests/gsetc_params.txt` and is shared by both
gates. The hardest kernels (`psf`, `sky`, `snr`) additionally carry independent
**scalar-transcription oracle** tests: a direct Python transliteration of the C
loop, compared against the vectorized implementation to ~1e-10 — when editing
those modules, keep the oracle tests meaningful (they are the guard against
silent vectorization errors). When porting or fixing kernel code, cite the
`gsetc.c` line in the docstring/comment as the existing modules do.
