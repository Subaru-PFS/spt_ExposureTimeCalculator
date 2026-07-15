# C-reference acceptance fixtures

This directory holds the frozen C+OpenMP reference runs the `slow`
regression gates (`tests/python/test_noise_reference.py`,
`tests/python/test_reference_outputs.py`) compare `pfsspecsim.etc` against.
**Everything listed below is a protected fixture — never edit or
regenerate it by hand.** (`.claude/hooks/protect-files.sh` blocks edits to
all of it via its first `case` pattern.) See `CLAUDE.md`'s "Non-negotiable
invariants" for the project-wide policy.

## Fixture sets

Four reference sets, each a `gsetc_params_*.txt` stdin script paired with
its `master_results/{noise,snc,snl,sno2}_<suffix>.dat` outputs:

| id (pytest) | params file                      | config                         | sky_type | output suffix  |
| ----------- | -------------------------------- | ------------------------------ | -------- | -------------- |
| `lr_11005`  | `tests/gsetc_params.txt`         | `tests/PFS.20211220.dat`       | 11005    | `_omp.dat`     |
| `mr_11005`  | `tests/gsetc_params_mr.txt`      | `tests/PFS.redMR.20211220.dat` | 11005    | `_mr.dat`      |
| `lr_11006`  | `tests/gsetc_params_lr11006.txt` | `tests/PFS.20240714.dat`       | 11006    | `_lr11006.dat` |
| `mr_11006`  | `tests/gsetc_params_mr11006.txt` | `tests/PFS.redMR.20240714.dat` | 11006    | `_mr11006.dat` |

`lr_11005` is the original reference set. The other three were generated
**2026-07-13** by re-running the same frozen `legacy/c_src/gsetc.c`, built
with the Makefile's exact flags:

```sh
cc -O3 -DHGCDTE_SUTR -DMOONLIGHT_ -o gsetc.x gsetc.c -lm
```

(Apple clang, macOS arm64). Every other observing-condition parameter
(seeing, zenith angle, field angle, moon geometry, exposure time/count,
sky-subtraction floor, diffuse stray light, degrade, mag file, effective
radius, line flux/width, min SNR) is identical to `tests/gsetc_params.txt`
— only the instrument config and, for the `*11006` sets, `sky_type` differ.
`tests/PFS.redMR.20211220.dat`, `tests/PFS.20240714.dat`, and
`tests/PFS.redMR.20240714.dat` are frozen config copies, verified
byte-identical to their `src/pfsspecsim/config/` counterparts at the time
they were captured.

MR sets have `snc`/`noise` output arm ids `{0, 3, 2}` (the MR red arm
remaps to output id 3, per `spectro_arm`/gsetc.c:76-83); LR sets have
`{0, 1, 2}`.

### Sanity check (2026-07-13)

Before generating the three new sets, the binary was re-run against the
existing `tests/gsetc_params.txt` to confirm the rebuilt binary reproduces
the original fixtures: `noise_omp.dat` and `snc_omp.dat` reproduced
byte-identically; `snl_omp.dat` differed only in its wavelength column
(87/93198 cells, max rel deviation 2.6e-5) and `sno2_omp.dat` only in one
wavelength column (3/108009 cells, max rel deviation 1.5e-5) — consistent
with last-printf-digit differences from serial-vs-OpenMP z-grid
accumulation order, not a behavioral change.

### Deriving the new params files

Each new `gsetc_params_*.txt` was derived from `tests/gsetc_params.txt`
with `sed`:

```sh
# mr (tests/gsetc_params_mr.txt)
sed -e '1s|.*|../PFS.redMR.20211220.dat|' -e 's/_omp\.dat/_mr.dat/' \
    tests/gsetc_params.txt > tests/gsetc_params_mr.txt

# lr11006 (tests/gsetc_params_lr11006.txt)
sed -e '1s|.*|../PFS.20240714.dat|' -e '3s/11005/11006/' \
    -e 's/_omp\.dat/_lr11006.dat/' \
    tests/gsetc_params.txt > tests/gsetc_params_lr11006.txt

# mr11006 (tests/gsetc_params_mr11006.txt)
sed -e '1s|.*|../PFS.redMR.20240714.dat|' -e '3s/11005/11006/' \
    -e 's/_omp\.dat/_mr11006.dat/' \
    tests/gsetc_params.txt > tests/gsetc_params_mr11006.txt
```

`gsetc.x` was then run from within `tests/` (so the `../PFS.*.dat` config
path resolves relative to `tests/`, and the relative output filenames land
in `tests/master_results/`) with each params file piped to stdin, e.g.:

```sh
cd tests && /path/to/gsetc.x < gsetc_params_mr.txt
```

## Other protected fixtures

- `tests/mag_18.dat` — magnitude-file input for the full-pipeline gate.
- `tests/analyze_diff.py` — reference diffing helper used during initial
  C-vs-Python validation.
- `legacy/c_src/gsetc.c` (and `modeldata.h`, `Makefile`) — the frozen
  original C source; see `legacy/c_src/README.md`.

## Achieved deviations

See `CLAUDE.md`'s "Non-negotiable invariants" section for the current
achieved max relative deviation per table, per reference set.
