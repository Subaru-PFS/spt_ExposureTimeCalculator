# PFS Exposure Time Calculator and Spectrum Simulator

## This package is developed by the following people

The original ETC was developed and written by Christopher Hirata (Ohio State University), which is based on the package developed for WFIRST (C. Hirata; arXiv:1204.5151) and altered for use in PFS project.

The code modification, the python wrapping, and the development of the spectral simulator were mainly done by Kiyoto Yabe, Yuki Moritani, Masato Onodera (Subaru Telescope), Atsushi Shimono (Kavli IPMU) and Robert Lupton (Princeton University).

## Release Note

- Version 1.0 Feb. 26, 2016
- Version 1.1 Apr. 27, 2016
- Version 1.2 Feb. 05, 2021
- Version 1.3 Jul. 31, 2024
- Version 2.0 The ETC engine (formerly the C program `gsetc.c`/`gsetc_omp.c`, driven as a subprocess) has been rewritten as a pure-Python package (`pfsspecsim.etc`), built on numpy/scipy/astropy. There is no C compiler, Makefile, or OpenMP dependency any more. See "Migrating from the old C-backed ETC" below if you have scripts built around the pre-2.0 interface.

## Requirements

- Python 3.11 or higher (needed for `tomllib`)
- numpy (2.0 and higher, but below 2.5 -- see the note below)
- scipy
- matplotlib (for plotting options)
- astropy (for ECSV/FITS file handling and the PFS datamodel package)

All required dependencies are installed automatically by `pip install .` or `uv sync` (see Installation below); there is no C compiler or OpenMP requirement.

**Note on the numpy upper bound**: `numpy<2.5` is pinned in `pyproject.toml`. This is *not* a requirement of the ETC engine itself (which works fine with newer numpy) -- it is forced by the pinned `pfs.datamodel @ w.2024.06` dependency, whose FITS I/O code uses a legacy dtype alias (`'a'`/`'|S'`) that numpy 2.5 removed. Once the `pfs.datamodel` pin is bumped past that fix, this upper bound can be lifted; that is out of scope for this package.

## Installation

Clone the repository and install with `pip`:

```bash
git clone https://github.com/Subaru-PFS/spt_ExposureTimeCalculator.git
cd spt_ExposureTimeCalculator
pip install .
```

To update:

```bash
git pull
pip install . --upgrade
```

### Alternative installation with `uv` (recommended for development)

```bash
git clone https://github.com/Subaru-PFS/spt_ExposureTimeCalculator.git
cd spt_ExposureTimeCalculator
uv sync
# or, to install into an existing environment instead of the project .venv:
uv pip install .
```

`uv sync` also installs the `dev` dependency group (pytest, black, ruff, ty) by default. To update:

```bash
git pull
uv sync
# or
uv pip install . --upgrade
```

Note: `pfs.datamodel` is installed automatically from GitHub as part of either install method.

Optional: `npm install` installs `prettier`/`markdownlint-cli2` locally so Claude
Code's Markdown auto-format hook runs fast; without it, the hook falls back to
`npx` (resolves the tool online on first use, then from npx's local cache).

Before you use the package, please read the sections below carefully.

## Description

This package includes two parts: the exposure time calculator (`pfs-spec etc`) and the spectral simulator (`pfs-spec sim`). You can get S/N information of an object in a given exposure time and various conditions by using the ETC, which is based on the "Chris Hirata's simulator" (now a pure-Python port of it). By using the results from the ETC, you can get the simulated spectra in the format of the current PFS datamodel with the spectral simulator.

Both are subcommands of the single `typer`-based `pfs-spec` command-line application. Every parameter can be given either as an individual `--option`, gathered into a TOML file passed via `--config`, or left at its built-in default. When the same parameter is given in more than one place, the priority is:

```text
CLI option  >  TOML file (--config)  >  built-in default
```

## Exposure Time Calculator (ETC)

Run `pfs-spec etc --help` for the full list of options. A minimal run using all defaults:

```bash
pfs-spec etc
```

writes `out/ref.noise.ecsv`, `out/ref.snc.ecsv`, and `out/ref.snl.ecsv` (the [OII]-doublet curve, `--outfile-oii`, is off by default). Passing options directly:

```bash
pfs-spec etc --mag 23.0 --line-flux 5.0e-17 --line-width 100 --exp-time 900 --exp-num 4
```

### TOML configuration files

Instead of (or in addition to) individual options, you can collect parameters in a TOML file and pass it with `--config`. TOML keys are the same snake_case names as the CLI options (with `-` replaced by `_`). A ready-to-run copy of the example below is committed at `examples/pfs_etc_example.toml`:

```toml
# examples/pfs_etc_example.toml
seeing = 0.80
zenith_ang = 35.00
galactic_ext = 0.00
moon_zenith_ang = 30.0
moon_target_ang = 60.0
moon_phase = 0.125
exp_time = 900
exp_num = 4
field_ang = 0.45

mag = 22.5          # or use mag_file = "path/to/spectrum.dat" instead (mutually exclusive)
reff = 0.3
line_flux = 1.0e-17
line_width = 70

mr_mode = false
outdir = "out"
outfile_noise = "ref.noise.ecsv"
outfile_snc = "ref.snc.ecsv"
outfile_snl = "ref.snl.ecsv"
```

```bash
pfs-spec etc --config examples/pfs_etc_example.toml
# override just one field on top of the TOML file:
pfs-spec etc --config examples/pfs_etc_example.toml --mag 23.0
```

`mag` (a flat AB magnitude) and `mag_file` (a 2-column wavelength[nm]/magnitude file) are mutually exclusive -- passing both (from any combination of CLI/TOML) is an error. Please note that the ETC does not convolve the instrument resolution into a `mag_file` spectrum, so the resolution should already be included in the input spectrum; also note that the wavelength is resampled with a sampling of 0.5 Angstrom (slightly larger than the PFS detector's own pixel sampling of ~0.7A/~0.8A/~0.8A for blue/red/NIR). Do not include "NaN" or other non-numeric values in the input file.

### Full parameter reference

| Parameter          | Default value            | Description                                                                                | Unit                          |
| :----------------- | :----------------------- | :----------------------------------------------------------------------------------------- | :---------------------------- |
| `seeing`           | 0.80                     | Seeing FWHM @800nm                                                                         | [arcsec]                      |
| `zenith_ang`       | 35.00                    | Zenith angle                                                                               | [deg.]                        |
| `galactic_ext`     | 0.00                     | Galactic extinction E(B-V)                                                                 | [ABmag.]                      |
| `field_ang`        | 0.45                     | Field angle                                                                                | [deg.; center=0, edge=0.675]  |
| `fiber_offset`     | 0.10                     | Fiber centering offset                                                                     | [arcsec]                      |
| `moon_zenith_ang`  | 30.0                     | Moon zenith angle (>=90 disables moonlight)                                                | [deg.]                        |
| `moon_target_ang`  | 60.0                     | Moon-target separation                                                                     | [deg.]                        |
| `moon_phase`       | 0.125                    | Moon phase                                                                                 | [0=New,0.25=quarter,0.5=full] |
| `exp_time`         | 900                      | Single exposure time                                                                       | [sec.]                        |
| `exp_num`          | 4                        | The number of exposures                                                                    |                               |
| `mag`              | 22.5                     | Flat AB magnitude; mutually exclusive with `mag_file`                                      | [ABmag]                       |
| `mag_file`         | (none)                   | 2-column wavelength[nm]/mag file; mutually exclusive with `mag`                            | *filename*                    |
| `reff`             | 0.3                      | Effective radius                                                                           | [arcsec]                      |
| `line_flux`        | 1.0e-17                  | Emission line flux                                                                         | [erg/s/cm^2]                  |
| `line_width`       | 70                       | Emission line width sigma                                                                  | [km/s]                        |
| `mr_mode`          | off                      | Medium resolution mode in the red arm                                                      | [bool]                        |
| `throughput_model` | "20240714"               | Throughput model tag                                                                       |                               |
| `spectrograph`     | "ave"                    | Spectrograph unit: ave\|sm1\|sm2\|sm3\|sm4                                                 |                               |
| `instr_config`     | (none)                   | Explicit spectrograph config file (overrides `throughput_model`/`spectrograph`/`mr_mode`)  | *filename*                    |
| `degrade`          | 1.0                      | Throughput degradation factor                                                              |                               |
| `obsc_fov_dep`     | on                       | Apply field-angle-dependent obscuration correction to `degrade`                            | [bool]                        |
| `sky_type`         | "11006"                  | Sky-model hex bitmask                                                                      |                               |
| `hgcdte_sutr`      | on                       | HgCdTe up-the-ramp sampling noise model                                                    | [bool]                        |
| `sky_sub_floor`    | 0.01                     | Sky-subtraction systematic floor fraction                                                  |                               |
| `diffuse_stray`    | 0.02                     | Diffuse stray-light fraction                                                               |                               |
| `oii_cat_in`       | (none)                   | Input catalog for [OII] emitters (enables the catalog output)                              | *filename*                    |
| `oii_cat_out`      | (none)                   | Output catalog for [OII] emitters                                                          | *filename*                    |
| `min_snr`          | 9.0                      | Minimum SNR for [OII] emission catalog detection                                           |                               |
| `n_workers`        | `min(8, os.cpu_count())` | Thread count for parallel stages; 1=serial. Results are identical regardless of this value |                               |
| `noise_reused`     | off                      | Reload the noise vector from `outfile_noise` instead of recomputing                        | [bool]                        |
| `overwrite`        | on                       | Allow overwriting existing output files                                                    | [bool]                        |
| `outdir`           | "out"                    | Output directory                                                                           |                               |
| `outfile_noise`    | "ref.noise.ecsv"         | Output file for the noise vector, relative to `outdir`                                     | *filename*                    |
| `outfile_snc`      | "ref.snc.ecsv"           | Output file for continuum S/N, relative to `outdir`                                        | *filename*                    |
| `outfile_snl`      | "ref.snl.ecsv"           | Output file for emission line S/N, relative to `outdir`                                    | *filename*                    |
| `outfile_oii`      | (none)                   | Output file for [OII] doublet S/N curve, relative to `outdir`                              | *filename*                    |

If you set `noise_reused` to true, the ETC will skip the process of generating the noise vector, which is a time-consuming step. The process time can be reduced to roughly half by this mode -- useful if you want to calculate the S/N of objects with various magnitudes and line fluxes under the same noise assumption (zenith angle, field angle, lunar condition, and exposure time). If you use this mode, `outfile_noise` must point at a previously written noise ECSV file.

### Output files (Astropy ECSV)

All ETC outputs are Astropy ECSV tables (plain text with a YAML header) -- readable with `astropy.table.Table.read(path, format="ascii.ecsv")` or as plain text. Every table's `table.meta["params"]` records the fully-resolved input `EtcParams` used to produce it (so an output file is self-describing), and `table.meta["degrade_resolved"]` records the actual (obscuration-corrected, if `obsc_fov_dep`) throughput degrade factor applied.

- `outfile_noise` (noise variance per pixel): `arm`, `pixel`, `wavelength` [nm, vacuum, pixel center], `variance` [e^2/pix, per exposure, without object continuum], `sky` [e/pix, per exposure]
- `outfile_snc` (continuum S/N per pixel): `arm`, `pixel`, `wavelength` [nm, vacuum, pixel left edge], `snr`, `signal` [e, per exposure], `noise_variance` [e^2/pix, without object continuum], `noise_variance_tot` [e^2/pix, with object continuum], `input_mag` [ABmag], `conversion_factor`, `sampling_factor`, `sky` [e]
- `outfile_snl` (single emission line S/N): `wavelength` [nm], `fiber_aperture_factor`, `effective_area` [m^2], `snr_b`, `snr_r` (or `snr_m` in `mr_mode`), `snr_n`, `snr_tot`
- `outfile_oii` (the [OII] doublet S/N curve vs. redshift, assuming a 1:1 line ratio): `z`, `wavelength0` ([OII]3726), `wavelength1` ([OII]3729), `fiber_aperture_factor`, `effective_area` [m^2], `snr_b`, `snr_r`, `snr_n`, `snr_tot`
- `oii_cat_out` (only written when `oii_cat_in` is given, driving an input catalog through the [OII] S/N model): `obj_id`, `z`, `reff`, `flux_oii`, `snr`

Some caution again:

1. `noise_reused` is only valid for the same observing condition (zenith angle, field angle, lunar condition, and exposure time). Please be careful when you use this option and do not reuse a noise file produced under different conditions.
2. You can input your own spectra by using `mag_file`, but please note that the ETC does not convolve the instrument resolution into the input spectrum internally, so the resolution should already be considered in the input spectrum. Also note that the wavelength sampling should be finer than the pixel sampling of the PFS detector (~0.7A, ~0.8A, and ~0.8A for blue, red, and NIR arm, respectively).

## Spectral Simulator

We can generate simulated outputs for use by the 1-D pipeline using a subset of the outputs from the ETC above and then running `pfs-spec sim`. The only input we use from the ETC is the continuum S/N ECSV (`outfile_snc`), and only a subset of its columns. In particular, the ETC's own `exp_num`/`mag` settings have *no* effect on the simulator's output.

Parameters that do matter (from the `pfs-spec etc` run) are `mr_mode` and anything affecting the observing conditions (exposure time, seeing, moon, ...): the settings used when running `pfs-spec etc` are the ones that determine the simulated spectrum.

So, before running `pfs-spec sim`, prepare an ETC output:

```bash
pfs-spec etc --no-mr-mode --exp-time 450 --outdir out --outfile-snc etc-t450-lr.ecsv
pfs-spec etc --mr-mode    --exp-time 450 --outdir out --outfile-snc etc-t450-mr.ecsv
```

You are now ready to generate simulated spectra. The input spectrum is given by `--mag`/`--mag-file` (same convention as the ETC) and the simulator adds appropriate noise, taking into account the number of integrations requested (`--nrealize`):

```bash
pfs-spec sim --etc-file out/etc-t450-lr.ecsv
```

Use `pfs-spec sim --help` to list all options, or pass `--config some.toml` for a TOML file of overrides (same CLI > TOML > default priority as `pfs-spec etc`). A more complete example, writing an ASCII table instead of FITS files:

```bash
pfs-spec sim --config sim_spec.toml --out-dir out --etc-file out/etc-t450-lr.ecsv --ascii-table test.sim --no-write-fits --mag 20.0
```

to write the ASCII file (`out/test.sim.dat`) instead of the FITS files, simulating a 20th-magnitude flat-spectrum source. If you specify `--nrealize` > 1 then multiple realizations of the noise are generated, each with its own `obj_id`, and the corresponding number of `pfsObject` files are written.

Key `pfs-spec sim` options:

| Parameter                                  | Default value      | Description                                                                                                               |
| :----------------------------------------- | :----------------- | :------------------------------------------------------------------------------------------------------------------------ |
| `etc_file`                                 | "out/ref.snc.ecsv" | Input noise file for the simulator (an ECSV `outfile_snc` from `pfs-spec etc`; legacy plain-text files are also accepted) |
| `nrealize`                                 | 1                  | The number of realizations                                                                                                |
| `ascii_table`                              | (none)             | Output ASCII table name without extension                                                                                 |
| `tract`/`patch`/`visit0`/`obj_id`/`cat_id` | 0/"0,0"/1/1/0      | Datamodel filename fields                                                                                                 |

The ASCII table columns are:

- (1) `WAVELENGTH` [nm]
- (2) `FLUX` [nJy]
- (3) `ERROR` [nJy]
- (4) `MASK` [1=masked]
- (5) `SKY` [nJy]
- (6) `ARM` [0=blue,1=red,2=NIR,3=redMR]

FITS output (`pfsArm`/`pfsObject`/`pfsConfig` files, following the [PFS datamodel](https://github.com/Subaru-PFS/datamodel/blob/master/datamodel.txt)) is written unless `--no-write-fits` is given.

### Realization of multiple spectra

If you have many spectra you want to realize, you can do that in a single run using an input magnitude file (`--mag-file`) containing each spectrum's information (columns: wavelength magnitude1 magnitude2 ... magnitudeN). Then you get one output file per spectrum. Please note `--nrealize 1` is required in this mode.

## Old -> new parameter name migration table

Pre-2.0 releases drove the C engine via `pfsspecsim.pfsetc.Etc` with ALL_CAPS `params` dict keys (or the `scripts/run_etc.py @file` argparse interface). Those interfaces still work (see "Migrating from the old C-backed ETC" below), but new code and TOML files should use the snake_case names below:

| Old (ALL_CAPS)                                                                                                                                         | New (snake_case)                                                                                                | Notes                                                     |
| :----------------------------------------------------------------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------------- | :-------------------------------------------------------- |
| `SEEING`/`ZENITH_ANG`/`GALACTIC_EXT`/`MOON_ZENITH_ANG`/`MOON_TARGET_ANG`/`MOON_PHASE`/`EXP_TIME`/`EXP_NUM`/`FIELD_ANG`/`REFF`/`LINE_FLUX`/`LINE_WIDTH` | same name, snake_case (`seeing`, `zenith_ang`, ... `line_width`)                                                | `exp_num` is now `int`                                    |
| `MAG_FILE` (a number)                                                                                                                                  | `mag`                                                                                                           | split into two mutually-exclusive fields                  |
| `MAG_FILE` (a path)                                                                                                                                    | `mag_file`                                                                                                      | split into two mutually-exclusive fields                  |
| `NOISE_REUSED`, `MR_MODE`, `OVERWRITE` (`'Y'`/`'N'`)                                                                                                   | `noise_reused`, `mr_mode`, `overwrite`                                                                          | now real `bool`                                           |
| `INFILE_OIICat` / `OUTFILE_OIICat`                                                                                                                     | `oii_cat_in` / `oii_cat_out`                                                                                    | `'-'` -> `None`                                           |
| `minSNR`                                                                                                                                               | `min_snr`                                                                                                       |                                                           |
| `degrade` / `SKY_SUB_FLOOR` / `DIFFUSE_STRAY` / `throughput_model` / `spectrograph` / `OUTDIR` / `obscFoVDep`                                          | `degrade` / `sky_sub_floor` / `diffuse_stray` / `throughput_model` / `spectrograph` / `outdir` / `obsc_fov_dep` |                                                           |
| `OUTFILE_NOISE`/`OUTFILE_SNC`/`OUTFILE_SNL`/`OUTFILE_OII`                                                                                              | `outfile_noise`/`outfile_snc`/`outfile_snl`/`outfile_oii`                                                       | `'-'` -> `None`; content is now ECSV                      |
| (hardcoded `SKYMODELS='11006'`)                                                                                                                        | `sky_type`                                                                                                      | now a configurable parameter                              |
| (hardcoded `OFFSET_FIB='0.10'`)                                                                                                                        | `fiber_offset`                                                                                                  | now a configurable parameter                              |
| (compile flag `-DHGCDTE_SUTR`)                                                                                                                         | `hgcdte_sutr`                                                                                                   | default `True`                                            |
| (compile flag `-DMOONLIGHT_`)                                                                                                                          | --                                                                                                              | always active; disable by setting `moon_zenith_ang >= 90` |
| `TMPDIR`, `BINDIR`                                                                                                                                     | --                                                                                                              | removed; accepted as no-ops by the compatibility layer    |

## Migrating from the old C-backed ETC

The pre-2.0 subprocess/C engine has been fully replaced by the pure-Python `pfsspecsim.etc` package; there is nothing left to compile. Existing code and notebooks built around the old interfaces keep working through two deprecated compatibility layers, each of which now internally drives the new pure-Python engine and emits a `DeprecationWarning`:

- **`pfsspecsim.pfsetc.Etc`** (`set_param`/`run`/`make_noise_model`/`make_snc`/`make_snl`/`make_sno2`/`run_multi`/...): same ALL_CAPS `params` dict as before, but `run()` now calls the pure-Python engine and every output *file* is now Astropy ECSV rather than the old whitespace-delimited plain text (the file name you configure via `OUTFILE_*` is unchanged; only the content format changed).
- **`pfs-run-etc`/`pfs-gen-sim-spec`** console scripts (and the equivalent `python scripts/run_etc.py`/`python scripts/gen_sim_spec.py @file` invocations): unchanged argparse interfaces on top of the same `Etc`/`Pfsspec` compatibility layers.

New code should use `pfs-spec etc`/`pfs-spec sim` (or `pfsspecsim.etc.load_params` + `run_etc_files`, and `pfsspecsim.sim.load_params` + `run_sim_spec`, directly -- see "Calling the modern API directly from Python" below) instead.

**Known behavioral difference in the chained `make_noise_model(); make_snc()` pattern**: the old C engine had no field-angle-dependent obscuration correction at all (`calc_obscuration` was applied only by the old *Python* wrapper on top of it, inconsistently, across the two calls in some scripts). The pure-Python engine resolves `degrade` (via `resolve_degrade`/`obsc_fov_dep`) exactly once per run and does not replicate that old double-application/compounding bug. Concretely, at the default `field_ang=0.45`, the obscuration correction factor is `calc_obscuration(0.45)[1] ~= 0.83`, so outputs produced through the compatibility layer's chained calls can differ from old, buggy, doubly-corrected runs by around that factor -- this is intentional; it is the C-parity path (single application of `resolve_degrade`) that was verified against the frozen C reference outputs (see "Numerical accuracy" below), not the old wrapper's compounding behavior.

```python
from pfsspecsim import pfsetc

etc = pfsetc.Etc()  # emits DeprecationWarning
etc.set_param('EXP_TIME', 1200)
etc.set_param('EXP_NUM', 3)
etc.set_param('OUTFILE_NOISE', 'out/ref.noise.ecsv')
etc.set_param('OUTFILE_SNC', 'out/ref.snc.ecsv')
etc.set_param('OUTFILE_SNL', 'out/ref.snl.ecsv')
etc.set_param('OUTFILE_OII', 'out/ref.snoii.dat')
etc.run()
```

```python
from pfsspecsim import pfsspec

sim = pfsspec.Pfsspec()
sim.set_param('ra', 150.0)
sim.set_param('dec', 2.0)
sim.set_param('etcFile', 'out/ref.snc.ecsv')  # legacy plain-text ETC files are still readable too
sim.set_param('MAG_FILE', 19.0)
sim.set_param('EXP_NUM', 16)
sim.set_param('asciiTable', 'test')
sim.set_param('nrealize', 1)
sim.set_param('plotObject', 't')
sim.set_param('plotArmSet', 'f')
sim.make_sim_spec()
```

See `example/notebooks/ETC Example.ipynb` for a fuller example (the `omp_num_threads=` constructor argument on `Etc` is still accepted for backward compatibility but has no effect -- there is no OpenMP thread pool any more).

### Calling the modern API directly from Python

New scripts do not need the compatibility layers above. Both engines expose a
snake_case dataclass of parameters plus a `run_*` function that returns the
results in memory (and, for the ETC, a `run_*_files` variant that also writes
the ECSV outputs):

```python
from pfsspecsim.etc import load_params, run_etc_files

params = load_params(overrides={
    "mag": 23.0,
    "exp_time": 900,
    "exp_num": 4,
    "outdir": "out",
})
results = run_etc_files(params)  # writes out/ref.{noise,snc,snl}.ecsv
print(results.snc)  # astropy.table.Table
```

```python
from pfsspecsim.sim import SimSpecParams, run_sim_spec

params = SimSpecParams(
    etc_file="out/ref.snc.ecsv",  # an outfile_snc ECSV written by pfs-spec etc
    mag=20.0,
    exp_num=4,
    nrealize=1,
    ascii_table="test",
    write_fits=False,  # skip pfsArm/pfsObject FITS output
)
sim = run_sim_spec(params)
print(sim.outdir, sim.asciiTable)  # 'out' 'test'
```

`load_params(toml_path=..., overrides=...)` (ETC) and
`SimSpecParams(...)`/`sim.load_params(toml_path=..., overrides=...)`
(simulator) both accept the same snake_case TOML files and CLI-priority
merging as the `pfs-spec etc`/`pfs-spec sim` command-line tools -- see the
"Full parameter reference" table above for the ETC's fields and
`pfsspecsim/sim/params.py`'s `SimSpecParams` dataclass for the simulator's.

## Multi processing

Multi-core processing for `pfsetc` is available using `run_multi`, where the number of cores is specified by `nproc`, the parameter name and the values that you want to calculate in parallel are specified by `param_name` and `param_values`. An example is shown below.

```python
from pfsspecsim import pfsetc

etc = pfsetc.Etc()
etc.set_param('EXP_TIME', 1200)
etc.set_param('EXP_NUM', 3)
etc.set_param('OUTFILE_NOISE', 'out/test.noise.dat')
etc.set_param('OUTFILE_SNC', 'out/test.snc.dat')
etc.set_param('OUTFILE_SNL', '-')
etc.set_param('OUTFILE_OII', '-')
etc.run_multi(nproc=3, param_name='MAG_FILE', param_values=['20.0', '21.0', '22.0', '23.0', '24.0'])
```

Multi-core processing for `pfsspec` is available using `make_sim_spec_multi`, where the number of cores is specified by `nproc`, all parameters are specified as a list of dictionaries `params`. An example is shown below.

```python
from pfsspecsim import pfsspec

params = [{'objId': 1000,
          'catId': 0.0,
          'fiberId': 71,
          'ra': 150.77132,
          'dec': 2.34267,
          'tract': 1,
          'patch': '0,0',
          'fiberMag': [19.5, 19.5, 19.5, 19.5, 19.5],
          'filterName': ['g', 'r', 'i', 'z', 'y'],
          'MAG_FILE': './tmp/tmp0.dat',
          'visit0': 100},
         {'objId': 1001,
          'catId': 0.0,
          'fiberId': 65,
          'ra': 150.02075,
          'dec': 2.47669,
          'tract': 2,
          'patch': '1,1',
          'fiberMag': [20.5, 20.5, 20.5, 20.5, 20.5],
          'filterName': ['g', 'r', 'i', 'z', 'y'],
          'MAG_FILE': './tmp/tmp1.dat',
          'visit0': 101},
         {'objId': 1002,
          'catId': 1.0,
          'fiberId': 42,
          'ra': 150.63364,
          'dec': 2.00197,
          'tract': 3,
          'patch': '2,2',
          'fiberMag': [21.5, 21.5, 21.5, 21.5, 21.5],
          'filterName': ['g', 'r', 'i', 'z', 'y'],
          'MAG_FILE': './tmp/tmp2.dat',
          'visit0': 102},
         {'objId': 1003,
          'catId': 1.0,
          'fiberId': 67,
          'ra': 150.74880,
          'dec': 2.25609,
          'tract': 4,
          'patch': '3,3',
          'fiberMag': [22.5, 22.5, 22.5, 22.5, 22.5],
          'filterName': ['g', 'r', 'i', 'z', 'y'],
          'MAG_FILE': './tmp/tmp3.dat',
          'visit0': 103}]

sim = pfsspec.Pfsspec()
sim.set_param('etcFile', 'out/test.snc.dat')
sim.set_param('EXP_NUM', 3)
sim.set_param('asciiTable', 'test')
sim.set_param('nrealize', 1)
sim.set_param('plotObject', 'f')
sim.set_param('plotArmSet', 'f')
sim.make_sim_spec_multi(nproc=3, params=params)
```

## Numerical accuracy vs. the old C engine

The pure-Python engine is verified against a frozen set of reference outputs produced by the old C engine (`gsetc_omp.x`) for a fixed set of observing conditions and inputs (see `tests/python/test_reference_outputs.py` and `tests/python/test_noise_reference.py`). The acceptance criterion is `rtol=1.5e-3` agreement on at least 99.9% of rows in every output table; the actual maximum relative deviation achieved during development is well within that margin for every table:

| Table                | Max relative deviation from C reference |
| :------------------- | :-------------------------------------- |
| noise                | 1.4e-05                                 |
| snc                  | 3.2e-04                                 |
| snl                  | 4.3e-04                                 |
| sno2 (`[OII]` curve) | 7.3e-04                                 |

### Known C-side quirks preserved in the port

A small number of C-side quirks and known bugs were intentionally *not* fixed in the port (they are ported verbatim, each with a comment pointing at the corresponding `gsetc.c` line numbers), to keep the numerical agreement above meaningful as a regression gate rather than silently changing the physics. If a result looks odd in one of the ways below, this is why -- it is a preserved property of the original C engine, not a porting bug. In practice these mostly matter for interpreting output columns correctly; none of them require changing how you invoke the ETC.

- **`outfile_snc`'s wavelength is the pixel *left edge*; `outfile_noise`'s is the pixel *center*.** The two tables' `wavelength` columns are computed on a half-pixel-shifted grid from each other (`gsetc.c:1401` vs `gsetc.c:882`) -- a genuine cross-function inconsistency in the original C code. If you cross-match rows between the two tables by index, remember the wavelengths are not the same grid.
- **`outfile_oii`'s `fiber_aperture_factor` diagnostic column ignores `field_ang`.** It is always computed as if the target were at the center of the field (`fieldang=0`, `gsetc.c:2049`), regardless of the actual `--field-ang` value used for the rest of the run. `outfile_snl`'s equivalent column does *not* have this quirk -- it correctly uses the real `field_ang` (`gsetc.c:2084`).
- **The `[OII]` catalog mode (`oii_cat_in`/`oii_cat_out`) ignores per-object line width.** The input catalog's velocity-dispersion (`sigma`) column is read but never used; every catalog object is scored with a hardcoded `sigma_v = 70 km/s` regardless of what the catalog file says (`gsetc.c:2154`).
- **The `[OII]` catalog's target count can undercount the output rows.** `oii_n_targets` (the redshift-recovery histogram total) only counts detected objects whose recovered redshift falls inside the histogram's fixed `[ZMIN_OII, ZMIN_OII + NZ_OII*DZ_OII)` range (`gsetc.c:2160-2165`). A detected object outside that range is still written to `oii_cat_out`, just not counted in the histogram/total -- so `n_targets` can be less than the number of rows in the output catalog.
- **The UVES and OH sky-line atlases are normalized to different reference airmasses.** UVES (optical) lines scale from a reference airmass of 1.1 (`gsetc.c:805-813`); OH (NIR) airglow lines scale from 1.0 (`gsetc.c:842-852`). This split is a property of the original model/data, not a bug introduced by the port.
- **`snr_single`'s (single-emission-line SNR) object-continuum transmission uses a fixed wavelength.** It is nominally a 41-point Gaussian average over the line's velocity-smeared wavelength range, but every quadrature point is evaluated at the same, unshifted line wavelength (`gsetc.c:1213-1218`), so the "average" is mathematically just the plain pointwise transmission. `get_signal`'s own line-shape average (used elsewhere) is genuinely velocity-smeared; only this particular transmission factor is not.

See the docstrings in `pfsspecsim/etc/*.py` (search for `QUIRK`) for the full technical detail behind each item above.

## Planned reorganization (not yet implemented)

The following restructuring is planned but not yet done -- noted here so the
intended direction is clear, not as a description of current behavior:

- **`example/`** (singular -- pre-2.0 sample noise/output files, spectra, and
  the `ETC Example.ipynb` notebook, all built around the old C-driven
  workflow) is expected to move under `legacy/` alongside `legacy/c_src` and
  `legacy/python_wrapper`, consistent with those being frozen references to
  the pre-2.0 toolchain.
- **`scripts/`** (top-level -- the thin `python scripts/run_etc.py`/
  `gen_sim_spec.py` shims described in "Migrating from the old C-backed ETC"
  above) is likewise expected to move under `legacy/` once the
  `pfs-run-etc`/`pfs-gen-sim-spec` console scripts are considered the sole
  supported path, or to be dropped entirely.
- **`examples/`** (plural -- currently just `pfs_etc_example.toml`) is
  expected to grow into the primary set of runnable examples for the modern
  `pfs-spec etc`/`pfs-spec sim` CLI, replacing what `example/` covers today:
  more TOML configs (e.g. an `mr_mode` comparison, an `oii_cat_in`/
  `oii_cat_out` catalog run, a multi-object `mag_file` realization) and
  updated example spectra/notebooks, drawn from `example/` and the snippets
  already in this README.

## Cautions and known issues

- The ETC has been partly validated using observed data for a limited type of objects taken during the commissionig, but the validation is still on-going.

- The throughput model is based on data taken during commissioning. However, there is a slight chromatic difference between the model predictions based on lab measurements and the actual observations. We are investigating the cause of this discrepancy.

- We have not verified yet a long integration for faint targets (i.e. whether the signal-to-noise ratio is proportional to the square root of the integration time).

- Data Reduction Pipeline is still under development, so the processed data itself has also not yet been fully validated. In particular, the validation of the ETC in the NIR regime remains uncertain.

- The nominal (default) values are mostly chosen based on statistics and typical values. However, the quality of individual spectra depends on various factors such as fiber configuration, guiding quality, and the location within the focal plane. Therefore, please note that actual observing results may differ from the ETC predictions.

## Other notes

- The other parameters that are implemented in the original Chris Hirata's ETC are indicated in `pfsspecsim/scripts/run_etc.py` (the argparse-based legacy entry point) including instrument setups and sky subtraction floor. Users can change these parameters on your own responsibility.

- The fraction of light covered by fiber aperture may be overestimated by up to 10% depending on the field angle compared to the ray-tracing calculation by using the PFS optical model including the telescope and the wide-field corrector.

- We assumed no fiber central position offset, additional 1% systematic sky subtraction error (on each 1D pixel), and additional 2% diffuse stray light (the entire 2D detector surface) for the noise model calculation.

- As a fiducial sky continuum model, we use a sky model generated based on the observations in the SDSS/BOSS in optical and MOSFIRE data in NIR and adjusted slightly with PFS commissioning data at <800 nm. We use the sky emission line model taken from UVES visible line atlas and theoretical model (Rousselot et al. 2000, A&A, 354, 1134) in NIR. As a fiducial atmospheric transmission model, we use Kitt Peak model for <900nm and a simulated ATRAN model with 3 mm PWV at longer wavelengths.

- The ETC physics (originally the C program `gsetc.c`, developed by Chris Hirata and partially modified by some people listed at the top of this page) is now implemented in the `pfsspecsim.etc` pure-Python package. The manual for the original code (`docs/Manual_v5.pdf`) may still be useful for understanding some assumptions in the noise calculation; the model data tables it documents are now bundled as `pfsspecsim/etc/data/modeldata.npz` (see `pfsspecsim/etc/data/README.md` and `tools/extract_modeldata.py` for provenance).

- The latest throughput model includes observed results of flux standard stars in the commissioning runs, which is avarage over all spectrographs.

Any feedback is highly appreciated!

## Contact

Kiyoto Yabe (Subaru Telescope, NAOJ)
e-mail: <kiyoyabe@naoj.org>
