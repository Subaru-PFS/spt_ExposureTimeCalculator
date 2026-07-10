"""The `etc` subcommand of `pfs-spec`, running the pure-Python ETC engine.

CLI options take priority over a TOML config file, which in turn takes
priority over `EtcParams` defaults (see `params.load_params`). Every
`EtcParams` field has a corresponding snake_case option (typer renders
`some_field` as `--some-field` automatically); every option defaults to
`None`, a sentinel meaning "not given on the command line" so that
`load_params`'s CLI > TOML > dataclass-default priority can tell "the user
typed this" apart from "this happens to equal the dataclass default".

`--mag`/`--mag-file` are handled as a genuinely mutually-exclusive pair
(gsetc.c `MAG_FILE` -> `EtcParams.mag`/`mag_file` split, see `params.py`):
passing both on the command line is a CLI-usage error (exit code != 0)
*before* `load_params`/`EtcParams.validate` even run, and passing either one
alone also clears the other in the `overrides` dict passed to
`load_params` -- otherwise a TOML file that sets the opposite one of the
pair would survive into the merged `EtcParams` and trip its own XOR check
downstream. Passing neither leaves both keys out of `overrides` entirely,
so a TOML file (or the dataclass defaults) fully control them as usual.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import typer

from ..etc import engine
from ..etc.params import EtcParams, load_params

_FIELD_NAMES = frozenset(f.name for f in dataclasses.fields(EtcParams))


def etc_command(
    config: Path | None = typer.Option(
        None, "--config", help="TOML file of EtcParams overrides (snake_case keys)."
    ),
    # --- Observing conditions ------------------------------------------
    seeing: float | None = typer.Option(
        None, help="Seeing FWHM @800nm, arcsec (default: 0.80)."
    ),
    zenith_ang: float | None = typer.Option(
        None, help="Zenith angle, deg (default: 35.0)."
    ),
    galactic_ext: float | None = typer.Option(
        None, help="Galactic extinction E(B-V), mag (default: 0.00)."
    ),
    field_ang: float | None = typer.Option(
        None, help="Field angle, deg (default: 0.45)."
    ),
    fiber_offset: float | None = typer.Option(
        None, help="Fiber centering offset, arcsec (default: 0.10)."
    ),
    moon_zenith_ang: float | None = typer.Option(
        None,
        help="Moon zenith angle, deg; >=90 disables moonlight (default: 30.0).",
    ),
    moon_target_ang: float | None = typer.Option(
        None, help="Moon-target angular separation, deg (default: 60.0)."
    ),
    moon_phase: float | None = typer.Option(
        None, help="Moon phase, 0=new .. 0.5=full (default: 0.125)."
    ),
    # --- Exposure --------------------------------------------------------
    exp_time: float | None = typer.Option(
        None, help="Exposure time per exposure, s (default: 900.0)."
    ),
    exp_num: int | None = typer.Option(
        None, help="Number of exposures to coadd (default: 4)."
    ),
    # --- Target ------------------------------------------------------------
    mag: float | None = typer.Option(
        None, help="Flat AB magnitude; mutually exclusive with --mag-file."
    ),
    mag_file: Path | None = typer.Option(
        None,
        help="2-column wavelength[nm]/mag file; mutually exclusive with --mag.",
    ),
    reff: float | None = typer.Option(
        None, help="Target effective radius, arcsec (default: 0.3)."
    ),
    line_flux: float | None = typer.Option(
        None, help="Emission line flux, erg/cm2/s (default: 1.0e-17)."
    ),
    line_width: float | None = typer.Option(
        None, help="Emission line width (sigma), km/s (default: 70.0)."
    ),
    # --- Instrument / model selection --------------------------------------
    mr_mode: bool | None = typer.Option(
        None, "--mr-mode/--no-mr-mode", help="Medium-resolution mode (default: off)."
    ),
    throughput_model: str | None = typer.Option(
        None, help="Throughput model tag (default: '20240714')."
    ),
    spectrograph: str | None = typer.Option(
        None, help="Spectrograph unit: ave|sm1|sm2|sm3|sm4 (default: 'ave')."
    ),
    instr_config: Path | None = typer.Option(
        None,
        help="Explicit spectrograph config file (overrides --throughput-model/--spectrograph/--mr-mode).",
    ),
    degrade: float | None = typer.Option(
        None, help="Throughput degrade factor (default: 1.0)."
    ),
    obsc_fov_dep: bool | None = typer.Option(
        None,
        "--obsc-fov-dep/--no-obsc-fov-dep",
        help="Apply the field-angle-dependent obscuration correction to --degrade (default: on).",
    ),
    sky_type: str | None = typer.Option(
        None, help="Sky-model hex bitmask (default: '11006')."
    ),
    hgcdte_sutr: bool | None = typer.Option(
        None,
        "--hgcdte-sutr/--no-hgcdte-sutr",
        help="HgCdTe up-the-ramp sampling model (default: on).",
    ),
    # --- Systematics ------------------------------------------------------
    sky_sub_floor: float | None = typer.Option(
        None, help="Sky-subtraction systematic floor fraction (default: 0.01)."
    ),
    diffuse_stray: float | None = typer.Option(
        None, help="Diffuse stray-light fraction (default: 0.02)."
    ),
    # --- [OII] catalog -----------------------------------------------------
    oii_cat_in: Path | None = typer.Option(
        None, help="Input [OII]-emitter catalog file (enables the catalog output)."
    ),
    oii_cat_out: Path | None = typer.Option(
        None, help="Output [OII]-emitter catalog ECSV path (not outdir-relative)."
    ),
    min_snr: float | None = typer.Option(
        None, help="Minimum SNR for [OII] catalog detection (default: 9.0)."
    ),
    # --- Execution control / output ----------------------------------------
    n_workers: int | None = typer.Option(
        None,
        help="Thread count for parallel stages; 1=serial (default: 3). "
        "Results are identical regardless of this value.",
    ),
    noise_reused: bool | None = typer.Option(
        None,
        "--noise-reused/--no-noise-reused",
        help="Reload the noise vector from --outfile-noise instead of recomputing (default: off).",
    ),
    overwrite: bool | None = typer.Option(
        None,
        "--overwrite/--no-overwrite",
        help="Allow overwriting existing noise/snc/snl output files (default: on).",
    ),
    outdir: Path | None = typer.Option(None, help="Output directory (default: 'out')."),
    outfile_noise: Path | None = typer.Option(
        None,
        help="Noise ECSV filename, relative to --outdir (default: 'ref.noise.ecsv').",
    ),
    outfile_snc: Path | None = typer.Option(
        None,
        help="Continuum SNR ECSV filename, relative to --outdir (default: 'ref.snc.ecsv').",
    ),
    outfile_snl: Path | None = typer.Option(
        None,
        help="Single-line SNR ECSV filename, relative to --outdir (default: 'ref.snl.ecsv').",
    ),
    outfile_oii: Path | None = typer.Option(
        None,
        help="[OII]-doublet SNR curve ECSV filename, relative to --outdir (default: disabled).",
    ),
) -> None:
    """Run the PFS exposure time calculator and write its ECSV outputs."""
    if mag is not None and mag_file is not None:
        typer.echo(
            "Error: --mag and --mag-file are mutually exclusive; give at most one.",
            err=True,
        )
        raise typer.Exit(code=1)

    local_values = locals()
    overrides = {
        name: local_values[name]
        for name in _FIELD_NAMES
        if local_values[name] is not None
    }
    if mag is not None:
        overrides["mag_file"] = None
    if mag_file is not None:
        overrides["mag"] = None

    try:
        params = load_params(config, overrides=overrides)
    except (ValueError, OSError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    try:
        results = engine.run_etc_files(params)
    except (
        Exception
    ) as exc:  # noqa: BLE001 -- surface any engine failure as a clean CLI error
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    outdir_resolved = Path(params.outdir)
    for label, outfile in (
        ("noise", params.outfile_noise),
        ("snc", params.outfile_snc),
        ("snl", params.outfile_snl),
        ("oii_curve", params.outfile_oii),
    ):
        if outfile is not None:
            typer.echo(f"{label}: {outdir_resolved / outfile}")
    if params.oii_cat_out is not None:
        typer.echo(f"oii_catalog: {params.oii_cat_out}")

    typer.echo(f"aperture_factor_800_target: {results.aperture_factor_800_target:.6f}")
    typer.echo(f"aperture_factor_800_point: {results.aperture_factor_800_point:.6f}")
