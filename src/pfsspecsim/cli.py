"""`typer` CLI app (`pfs-sim-spec`) wrapping `pfsspecsim.pfsspec.Pfsspec`.

Mirrors `pfsspecsim.etc.cli`'s pattern (see that module's docstring for the
full rationale): every `Pfsspec.params` entry has a corresponding
snake_case CLI option, every option defaults to `None` -- a sentinel
meaning "not given on the command line" -- and an optional `--config` TOML
file can supply overrides too, with priority CLI > TOML > `Pfsspec`'s own
built-in defaults (`Pfsspec.__init__`). Unlike `EtcParams`, `Pfsspec` has no
dataclass of its own (legacy mixed-case `params` dict, e.g. `EXP_NUM`,
`MAG_FILE`, `countsMin`, `outDir`, ...); `_KEY_MAP` below is this module's
own snake_case-option -> legacy-dict-key migration table, used only here.

`--mag`/`--mag-file` are a mutually exclusive pair mapping onto the single
legacy `MAG_FILE` parameter (`Pfsspec.make_sim_spec` tells them apart by
whether `float()` succeeds on it -- see pfsspec.py:175-178): passing both
on the command line is a CLI-usage error (exit code != 0); passing either
one alone clears the other from `overrides` so a TOML file's opposite
setting cannot leak through into the merged value.

`--etc-file` (the legacy `etcFile` parameter) is read by
`Pfsspec.make_sim_spec` as an Astropy ECSV file (T14), typically the
`outfile_snc` ECSV written by `pfs-etc`; old-format plain-text files are
still accepted via that method's fallback path.
"""

from __future__ import annotations

import tomllib
from importlib import metadata
from pathlib import Path
from typing import Any

import typer

from .pfsspec import Pfsspec

#: Distribution name to look up for `--version` (see `pyproject.toml`).
_PKG_NAME = "pfsspecsim"

app = typer.Typer(
    name="pfs-sim-spec",
    add_completion=False,
    no_args_is_help=False,
    help="PFS spectral simulator (pure-Python engine).",
)

#: CLI/TOML snake_case option name -> legacy `Pfsspec.params` dict key.
#: `mag`/`mag_file` both map onto the single `MAG_FILE` legacy key; which
#: one wins is resolved by the mutual-exclusivity handling in `main`.
_KEY_MAP: dict[str, str] = {
    "exp_num": "EXP_NUM",
    "mag": "MAG_FILE",
    "mag_file": "MAG_FILE",
    "counts_min": "countsMin",
    "etc_file": "etcFile",
    "nrealize": "nrealize",
    "out_dir": "outDir",
    "ascii_table": "asciiTable",
    "ra": "ra",
    "dec": "dec",
    "tract": "tract",
    "patch": "patch",
    "visit0": "visit0",
    "cat_id": "catId",
    "obj_id": "objId",
    "fiber_id": "fiberId",
    "fiber_mag": "fiberMag",
    "filter_name": "filterName",
    "spectrograph": "spectrograph",
    "pfs_config_full": "pfsConfigFull",
    "write_fits": "writeFits",
    "write_pfs_arm": "writePfsArm",
    "plot_arm_set": "plotArmSet",
    "plot_object": "plotObject",
    "sky_sub_floor": "SKY_SUB_FLOOR",
    "sky_sub_mode": "SKY_SUB_MODE",
    "sky_sub_seed": "SKY_SUB_SEED",
}

#: Legacy keys whose value `Pfsspec.make_sim_spec` parses with `strToBool`
#: (pfsspec.py's own "1"/"t"/"true" vs "0"/"f"/"false" string convention,
#: case-insensitive) rather than treating it as a native Python bool.
_STRTOBOOL_KEYS = frozenset({"writeFits", "writePfsArm", "pfsConfigFull"})

#: Every recognized CLI/TOML option name (the keys of `_KEY_MAP`).
_FIELD_NAMES = frozenset(_KEY_MAP)


def _package_version() -> str:
    try:
        return metadata.version(_PKG_NAME)
    except metadata.PackageNotFoundError:
        return "0.0.0+unknown"


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"pfs-sim-spec {_package_version()}")
        raise typer.Exit()


def _load_toml(path: Path) -> dict[str, Any]:
    with open(path, "rb") as fh:
        data = tomllib.load(fh)
    unknown = set(data) - _FIELD_NAMES
    if unknown:
        raise ValueError(f"Unknown key(s) in TOML file {path}: {sorted(unknown)}")
    return data


def _apply_overrides(sim: Pfsspec, overrides: dict[str, Any]) -> None:
    """Push each override through `Pfsspec.set_param`, translating the
    snake_case name to its legacy dict key and adapting the value to the
    string convention `make_sim_spec` expects for `_STRTOBOOL_KEYS`
    (`plotArmSet`/`plotObject` are consumed as native Python bools by
    `make_sim_spec`, so those two pass through unchanged).
    """
    for name, value in overrides.items():
        key = _KEY_MAP[name]
        if key in _STRTOBOOL_KEYS and isinstance(value, bool):
            value = "true" if value else "false"
        sim.set_param(key, value)


@app.command()
def main(
    config: Path | None = typer.Option(
        None,
        "--config",
        help="TOML file of Pfsspec parameter overrides (snake_case keys).",
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
    out_dir: Path | None = typer.Option(None, help="Output directory (default: 'out')."),
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

    _special = frozenset({"mag", "mag_file", "fiber_mag", "filter_name"})
    local_values = locals()
    overrides: dict[str, Any] = {
        name: local_values[name]
        for name in _FIELD_NAMES
        if name not in _special and local_values.get(name) is not None
    }
    if mag is not None:
        overrides["mag"] = mag
    if mag_file is not None:
        overrides["mag_file"] = mag_file
    if fiber_mag is not None:
        overrides["fiber_mag"] = [float(v) for v in fiber_mag.split(",")]
    if filter_name is not None:
        overrides["filter_name"] = [v.strip() for v in filter_name.split(",")]

    try:
        toml_overrides = _load_toml(config) if config is not None else {}
    except (ValueError, OSError) as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    # Both `mag` and `mag_file` map onto the single legacy `MAG_FILE` key,
    # so the merge must never end up holding both: passing either on the
    # command line clears the other from the TOML layer (same rationale as
    # pfs-etc's handling), and a TOML file supplying both is an error
    # rather than a silent dict-order race.
    if mag is not None:
        toml_overrides.pop("mag_file", None)
    if mag_file is not None:
        toml_overrides.pop("mag", None)
    if "mag" in toml_overrides and "mag_file" in toml_overrides:
        typer.echo(
            "Error: 'mag' and 'mag_file' are mutually exclusive in the TOML "
            "config; give at most one.",
            err=True,
        )
        raise typer.Exit(code=1)

    merged = dict(toml_overrides)
    merged.update(overrides)

    sim = Pfsspec()
    try:
        _apply_overrides(sim, merged)
        result = sim.make_sim_spec()
    except SystemExit as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1) from exc
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
