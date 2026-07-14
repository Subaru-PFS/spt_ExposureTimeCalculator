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

CLIs (installed as entry points): `pfs-spec` (umbrella typer app, subcommands
`etc` and `sim`), and the deprecated `pfs-run-etc` / `pfs-gen-sim-spec`
(argparse `@file` interface over the legacy `Etc` wrapper). Example:
`uv run pfs-spec etc --config examples/pfs_etc_example.toml`.

Note: editing a `.py` file auto-runs `ruff --fix` + `black` on it via a
PostToolUse hook, so you don't need to format manually after edits.

Note: editing a `.md` file similarly auto-runs `prettier` (general formatting,
incl. table-column alignment) then `markdownlint-cli2 --fix` (project style
rules from `.markdownlint.jsonc`) via a PostToolUse hook. These tools come
from `package.json` (optional — run `npm install` once to get fast local
binaries via `node_modules/.bin`); without it, the hook transparently falls
back to `npx`, which resolves the tool online on first use and from npx's
local cache thereafter.

## Non-negotiable invariants

**Never modify these C-reference fixtures** — they are the acceptance criteria,
produced by the original C+OpenMP binary and used by the regression gates:
`tests/master_results/`, `tests/gsetc_params.txt`, `tests/PFS.20211220.dat`,
`tests/mag_18.dat`, `tests/analyze_diff.py`, `tests/gsetc_params_mr.txt`,
`tests/gsetc_params_lr11006.txt`, `tests/gsetc_params_mr11006.txt`,
`tests/PFS.redMR.20211220.dat`, `tests/PFS.20240714.dat`,
`tests/PFS.redMR.20240714.dat`. (A PreToolUse hook blocks edits to them.) The
original C source vendored at `legacy/c_src/` (see Architecture) is likewise
frozen and must never be edited. See `tests/README.md` for full fixture
provenance.

**The slow gates are the arbiter of correctness.** `test_noise_reference.py`
(noise only) and `test_reference_outputs.py` (full pipeline) each run
parametrized across four frozen C-reference sets (`conftest.REFERENCE_SETS`:
LR/sky_type=11005, MR/11005, LR/11006, MR/11006 — see `tests/README.md`), for
20 slow tests total. Every output column is compared against the matching C
reference at `rtol=1.5e-3` and requires ≥99.9% of values to match. Any change
to a numeric kernel
(`constants/config/materials/extinction/atmosphere/psf/sky/noise/snr/engine`)
**must** keep `uv run pytest tests/python -m slow` green. Current achieved max
relative deviation per reference set: `lr_11005` (original) noise 1.4e-05,
snc 3.2e-04, snl 4.3e-04, sno2 7.3e-04; `mr_11005` noise 5.1e-06,
snc 3.2e-04, snl 4.3e-04, sno2 7.3e-04; `lr_11006` noise 1.7e-05,
snc 3.3e-04, snl 3.0e-04, sno2 8.2e-04; `mr_11006` noise 5.7e-06,
snc 3.3e-04, snl 3.0e-04, sno2 8.2e-04 (all comfortably under the
`rtol=1.5e-3` gate).

**Do not "fix" the ported quirks.** `pfsspecsim.etc` is a faithful port of
`gsetc.c` (vendored, frozen, at `legacy/c_src/gsetc.c`). Several intentional
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

Package is **src-layout** (`src/pfsspecsim/`). The ETC engine lives in
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
```

The simulator engine lives in `src/pfsspecsim/sim/`, mirroring `etc/`'s
package shape: `params.py` (`SimSpecParams`, `load_params`), `engine.py`
(`run_sim_spec`, translating a `SimSpecParams` into `Pfsspec.set_param` calls),
`pfsspec.py` (the `Pfsspec` engine class — datamodel FITS/pfsArm output),
`dm_utils.py` (the `pfs.datamodel` object builders `pfsspec.py` calls), and
`__init__.py` re-exporting `SimSpecParams`/`load_params`/`run_sim_spec`. The
pure back-compat wrapper, `pfsetc.Etc`, lives separately in
`src/pfsspecsim/legacy/` (see Backward-compatibility below). The CLI is a
separate top-level package, `src/pfsspecsim/cli/`, which imports both engines
but nothing imports it back — `app.py` (the `pfs-spec` umbrella typer app +
shared `--version`), `etc.py` (the `etc` subcommand, drives `pfsspecsim.etc`),
`sim.py` (the `sim` subcommand, drives `pfsspecsim.sim`).

Public API (`from pfsspecsim.etc import ...`): `EtcParams` (dataclass of all
inputs; `load_params(toml, overrides)` with priority CLI > TOML > defaults),
`run_etc(params) -> EtcResults` (no file I/O), `run_etc_files(params)` (runs +
writes ECSV), plus `MagSpec`, `resolve_degrade`, `calc_obscuration`.

**Data tables** (UVES/OH sky lines, atmospheric transmission, Si index, dust
norm) were extracted once from the C `modeldata.h` by `tools/extract_modeldata.py`
into `src/pfsspecsim/etc/data/modeldata.npz` (loaded via `importlib.resources`

- `lru_cache`; `etc/data/README.md` records provenance). Regenerate only via
  that script — the C source is kept frozen at `legacy/c_src/` (see
  `legacy/c_src/README.md`), not re-derived elsewhere.

**Outputs are Astropy ECSV** (`.ecsv`). Every table carries `table.meta` with
the resolved `EtcParams`, `etc_version`, `instr_config`, and top-level
`degrade_resolved` (the obscuration-corrected value the engine actually
applied — distinct from the raw `params.degrade` under `meta["params"]`).

**Backward-compatibility layer.** `src/pfsspecsim/legacy/pfsetc.py` keeps the
legacy `Etc` class surface (uppercase param dict, `set_param`/`run`/`make_*`/
`run_multi`, `__init__(omp_num_threads=…)` accepted+ignored with a
`DeprecationWarning`). It translates old params to `EtcParams`, drives the new
engine, and restores the old `nsm_*`/`snc_*`/`snl_*`/`sno2_*` attributes from
the result tables so old scripts and notebooks keep working. One deliberate
behavioral change: the old `degrade`-compounding bug (obscuration correction
accumulating across repeated `run()` calls) is **not** replicated, so a
chained `make_noise_model(); make_snc()` differs from pre-2.0 output by an
obscuration factor (~0.83 at `field_ang=0.45`). `src/pfsspecsim/sim/pfsspec.py`
reads the ETC's SNC output as ECSV (with a legacy `np.loadtxt` fallback for
old-format files); its datamodel FITS output path is unchanged. `legacy/` is
for this pure-wrapper-only module; `pfsspec.py`/`dm_utils.py` are current,
non-deprecated implementation code and live in `sim/` instead, not `legacy/`.
The pre-2.0 top-level import paths (`pfsspecsim.pfsetc`/`.pfsspec`/
`.dm_utils`) are kept working via a `sys.modules` alias block in
`src/pfsspecsim/__init__.py`, so `import pfsspecsim.pfsspec` and
`from pfsspecsim.pfsetc import Etc` still resolve to the real modules above.

## Testing conventions

`tests/python/` mirrors the module layout (`test_<module>.py`). The `slow`
marker (see `pyproject.toml`) gates the two C-reference regression test files
(`test_noise_reference.py`, `test_reference_outputs.py`), each parametrized
over `conftest.REFERENCE_SETS`' four fixture sets for 20 slow tests total; the
default `-m "not slow"` run stays fast. `tests/python/conftest.py::reference_params`
builds the `EtcParams` matching `tests/gsetc_params.txt` (overridden per
reference set with `instr_config`/`sky_type`) and is shared by both gate
files. The hardest kernels (`psf`, `sky`, `snr`) additionally carry independent
**scalar-transcription oracle** tests: a direct Python transliteration of the C
loop, compared against the vectorized implementation to ~1e-10 — when editing
those modules, keep the oracle tests meaningful (they are the guard against
silent vectorization errors). When porting or fixing kernel code, cite the
`gsetc.c` line in the docstring/comment as the existing modules do.
