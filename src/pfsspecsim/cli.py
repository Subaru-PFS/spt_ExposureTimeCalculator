"""`typer` CLI app (`pfs-sim-spec`) wrapping `pfsspecsim.simspec.run_sim_spec`.

Mirrors `pfsspecsim.etc.cli`'s pattern: every `SimSpecParams` field has a
corresponding snake_case CLI option, every option defaults to `None` -- a
sentinel meaning "not given on the command line" -- and an optional
`--config` TOML file can supply overrides too, with priority
CLI > TOML > `SimSpecParams`'s own dataclass defaults (`simspec.load_params`).

`--mag`/`--mag-file` are a genuinely mutually-exclusive pair (both fields on
`SimSpecParams`, mirroring `EtcParams.mag`/`mag_file`): passing both on the
command line is a CLI-usage error (exit code != 0) before `load_params` even
runs; passing either one alone also clears the other in the `overrides` dict
so a TOML file's opposite setting cannot leak through into the merged value.

`--fiber-mag`/`--filter-name` are comma-separated strings on the command
line (there is no repeated-option or list syntax here), parsed into
`SimSpecParams.fiber_mag`/`filter_name` (`list[float]`/`list[str]`) before
being added to `overrides`; a TOML file or a direct `SimSpecParams(...)` call
supplies these as native TOML/Python arrays instead.

`--etc-file` (`SimSpecParams.etc_file`) is read by `Pfsspec.make_sim_spec`
as an Astropy ECSV file (T14), typically the `outfile_snc` ECSV written by
`pfs-etc`; old-format plain-text files are still accepted via that method's
fallback path.
"""

from __future__ import annotations

import dataclasses
from importlib import metadata
from pathlib import Path
from typing import Any

import typer

from . import simspec
from .simspec import SimSpecParams

#: Distribution name to look up for `--version` (see `pyproject.toml`).
_PKG_NAME = "pfsspecsim"

app = typer.Typer(
    name="pfs-sim-spec",
    add_completion=False,
    no_args_is_help=False,
    help="PFS spectral simulator (pure-Python engine).",
)

_FIELD_NAMES = frozenset(f.name for f in dataclasses.fields(SimSpecParams))


def _package_version() -> str:
    try:
        return metadata.version(_PKG_NAME)
    except metadata.PackageNotFoundError:
        return "0.0.0+unknown"


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"pfs-sim-spec {_package_version()}")
        raise typer.Exit()


@app.command()
def main(
    config: Path | None = typer.Option(
        None,
        "--config",
        help="TOML file of SimSpecParams overrides (snake_case keys).",
    ),
    etc_file: Path | None = typer.Option(
        None,
        help="ETC continuum-SNR ECSV input file (default: 'out/ref.snc.dat').",
    ),
    exp_num: int | None = typer.Option(
        None, help="Number of exposures to coadd (default: 4)."
    ),
    mag: float | None = typer.Option(
        None, help="Flat AB magnitude; mutually exclusive with --mag-file."
    ),
    mag_file: Path | None = typer.Option(
        None,
        help="2-column wavelength[nm]/mag file (or multi-object magnitude "
        "table); mutually exclusive with --mag.",
    ),
    counts_min: float | None = typer.Option(
        None, help="Minimum counts per pixel used for the noise floor (default: 0.1)."
    ),
    nrealize: int | None = typer.Option(
        None, help="Number of realizations (default: 1)."
    ),
    out_dir: Path | None = typer.Option(
        None, help="Output directory (default: 'out')."
    ),
    ascii_table: str | None = typer.Option(
        None,
        help="Ascii table basename to also write, relative to --out-dir "
        "(default: disabled).",
    ),
    ra: float | None = typer.Option(None, help="Fiber RA, deg (default: 150.0)."),
    dec: float | None = typer.Option(None, help="Fiber Dec, deg (default: 2.0)."),
    tract: int | None = typer.Option(None, help="Tract (default: 0)."),
    patch: str | None = typer.Option(None, help="Patch, 'x,y' (default: '0,0')."),
    visit0: int | None = typer.Option(None, help="First visit number (default: 1)."),
    cat_id: int | None = typer.Option(None, help="Catalog id (default: 0)."),
    obj_id: int | None = typer.Option(None, help="Object id (default: 1)."),
    fiber_id: int | None = typer.Option(None, help="Fiber id (default: 1)."),
    fiber_mag: str | None = typer.Option(
        None,
        help="Comma-separated fiducial fiber magnitudes for g,r,i,z,y "
        "(default: '22.5,22.5,22.5,22.5,22.5').",
    ),
    filter_name: str | None = typer.Option(
        None,
        help="Comma-separated filter names for g,r,i,z,y "
        "(default: 'hcs_g,hcs_r,hcs_i,hcs_z,hcs_y').",
    ),
    spectrograph: int | None = typer.Option(
        None, help="Spectrograph number (default: 1)."
    ),
    pfs_config_full: bool | None = typer.Option(
        None,
        "--pfs-config-full/--no-pfs-config-full",
        help="Write a full PfsConfig (default: off).",
    ),
    write_fits: bool | None = typer.Option(
        None,
        "--write-fits/--no-write-fits",
        help="Write datamodel FITS files (default: on).",
    ),
    write_pfs_arm: bool | None = typer.Option(
        None,
        "--write-pfs-arm/--no-write-pfs-arm",
        help="Write pfsArm FITS files, requires --write-fits (default: on).",
    ),
    plot_arm_set: bool | None = typer.Option(
        None,
        "--plot-arm-set/--no-plot-arm-set",
        help="Plot the pfsArmSet data (default: off).",
    ),
    plot_object: bool | None = typer.Option(
        None,
        "--plot-object/--no-plot-object",
        help="Plot the pfsObject data (default: off).",
    ),
    sky_sub_floor: float | None = typer.Option(
        None, help="Sky-subtraction systematic floor fraction (default: 0.01)."
    ),
    sky_sub_mode: str | None = typer.Option(
        None,
        help="Sky-subtraction residual mode: random|systematic|wavecalib|psfvar "
        "(default: 'random').",
    ),
    sky_sub_seed: int | None = typer.Option(
        None, help="Sky-subtraction residual RNG seed (default: 0)."
    ),
    version: bool | None = typer.Option(
        None,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show the pfs-sim-spec version and exit.",
    ),
) -> None:
    """Simulate PFS spectra from an ETC continuum-SNR ECSV file."""
    if mag is not None and mag_file is not None:
        typer.echo(
            "Error: --mag and --mag-file are mutually exclusive; give at most one.",
            err=True,
        )
        raise typer.Exit(code=1)

    _special = frozenset({"fiber_mag", "filter_name"})
    local_values = locals()
    overrides: dict[str, Any] = {
        name: local_values[name]
        for name in _FIELD_NAMES
        if name not in _special and local_values.get(name) is not None
    }
    if mag is not None:
        overrides["mag_file"] = None
    if mag_file is not None:
        overrides["mag"] = None
    if fiber_mag is not None:
        overrides["fiber_mag"] = [float(v) for v in fiber_mag.split(",")]
    if filter_name is not None:
        overrides["filter_name"] = [v.strip() for v in filter_name.split(",")]

    try:
        params = simspec.load_params(config, overrides=overrides)
    except (ValueError, OSError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    try:
        sim = simspec.run_sim_spec(params)
    except Exception as exc:  # noqa: BLE001 -- surface any failure as a clean CLI error
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    typer.echo(f"outdir: {sim.outdir}")
    if sim.asciiTable != "None":
        typer.echo(f"ascii_table: {sim.outdir}/{sim.asciiTable}")
    if sim.writeFits:
        typer.echo(f"n_pfsObject: {len(sim.pfsObjects)}")


if __name__ == "__main__":
    app()
