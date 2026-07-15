# legacy/python_wrapper -- original Python wrapper (frozen reference)

This directory vendors the pre-2.0 Python wrapper that drove `gsetc.x` (built
from `legacy/c_src/`) as a subprocess. It is kept as a permanent, frozen
reference; it is **not** built, packaged, or imported by anything in this
repository. It was vendored byte-identical from the `master` branch at commit
`19245a8bdd1311283bef47d601d59dd3051e0a3b`.

Files (subpaths preserved from their `master`-branch locations):

- `pfsspecsim/__init__.py`, `pfsspecsim/pfsetc.py`, `pfsspecsim/pfsspec.py`,
  `pfsspecsim/dm_utils.py` -- vendored from `python/pfsspecsim/` on `master`.
  `pfsetc.py` defines the `Etc` class: the subprocess wrapper around
  `gsetc.x` (uppercase param dict, `set_param`/`run`/`make_noise_model`/
  `make_snc`/`make_snl`/`make_sno2`/`run_multi`). `pfsspec.py`/`dm_utils.py`
  are the pre-2.0 simulator engine and `pfs.datamodel` object builders.
  `python/pfsspecsim/config/` (instrument `.dat` files, already packaged
  under `src/pfsspecsim/config/`) and `python/pfsspecsim/bin/` (compiled
  binaries; the build recipe is already vendored at `legacy/c_src/Makefile`)
  were intentionally **not** vendored here.
- `scripts/__init__.py`, `scripts/run_etc.py`, `scripts/run_etc.defaults`,
  `scripts/gen_sim_spec.py`, `scripts/gen_sim_spec.defaults` -- vendored from
  `scripts/` on `master`. These are the argparse `@file`-interface CLI
  entry points (`run_etc.py`/`gen_sim_spec.py`) and their default-parameter
  files, superseded by `pfs-run-etc`/`pfs-gen-sim-spec` in this repo.

## Why this is kept

`src/pfsspecsim/etc/` and `src/pfsspecsim/legacy/pfsetc.py` are a pure-Python
port of, and compatibility shim for, this wrapper. Several docstrings across
the 2.0 codebase cite `pfsetc.py:<line>` to point at the exact statement in
the legacy `Etc` class a given behavior was ported from or is deliberately
*not* replicating. This copy is what those citations resolve against.

It also documents the origin of two 2.0 defaults and one deliberate
behavioral change:

- `EtcParams.sky_type` defaults to `"11006"`, inherited from
  `SKYMODELS = '11006'` at `pfsspecsim/pfsetc.py:17`.
- `EtcParams.fiber_offset` defaults to `0.10`, inherited from
  `OFFSET_FIB = '0.10'` at `pfsspecsim/pfsetc.py:18`.
- `scripts/run_etc.defaults` and `scripts/gen_sim_spec.defaults` record the
  legacy `@file`-interface default parameters for `run_etc.py`/
  `gen_sim_spec.py`, the argparse ancestors of `pfs-run-etc`/
  `pfs-gen-sim-spec`.

## Line-number contract -- do not edit this copy

Docstrings and comments throughout `src/pfsspecsim/etc/` and
`src/pfsspecsim/legacy/pfsetc.py` cite `pfsetc.py:<line>` (occasionally
spelled `pfsspecsim/pfsetc.py:<line>`) to point at the exact statement in the
legacy wrapper a given Python function, constant, or deliberately-not-ported
quirk corresponds to. Every one of those citations resolves against **this
exact copy** of `pfsspecsim/pfsetc.py`.

Because of that, this copy must never be edited, reformatted, or otherwise
modified -- treat it exactly like `legacy/c_src/gsetc.c` and the fixtures
under `tests/master_results/`: a frozen input, not a file to improve. This
`README.md` is the one file in this directory that may be edited.

## Known bug: the degrade-compounding bug

`pfsspecsim/pfsetc.py`'s `Etc` class recomputes the field-angle obscuration
correction and multiplies it into the *stored* `self.params['degrade']` at
each of its five public entry points, rather than applying the correction to
a local/derived value:

- `run()` -- `degrade = float(self.params['degrade']) * corr` at line 254
  (guarded by the `obscFoVDep` check starting at line 252).
- `make_noise_model()` -- line 450.
- `make_snc()` -- line 541.
- `make_snl()` -- line 620.
- `make_sno2()` -- line 697.

Because each of these mutates `self.params['degrade']` in place, calling more
than one of them on the same `Etc` instance (e.g. the common
`make_noise_model(); make_snc()` pattern) compounds the correction factor
`corr` across calls: after *n* calls the effective degrade carries `corr**n`
instead of `corr`. `pfsspecsim` 2.0 deliberately does **not** reproduce this:
`src/pfsspecsim/legacy/pfsetc.py` applies the obscuration correction once per
result, not cumulatively, so chained legacy-style calls produce output that
differs from pre-2.0 by a factor of `corr` (~0.83 at `field_ang=0.45`) per
extra call beyond the first. See
`docs/review-c-to-python-port-2026-07-10.md` §6.1 for the full analysis.

## Not importable, not packaged

This directory lives outside `src/` and is never imported by this package --
it is a historical reference only, kept for the line-number contract and bug
documentation above. If you need a working, supported drop-in replacement for
the legacy `Etc` class, use `src/pfsspecsim/legacy/pfsetc.py` (see the
top-level `CLAUDE.md`, "Backward-compatibility layer").
