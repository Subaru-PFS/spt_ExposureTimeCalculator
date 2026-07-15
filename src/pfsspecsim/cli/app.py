"""The `pfs-spec` umbrella `typer` app, hosting the `etc` and `sim` subcommands.

`pfs-spec etc` runs the pure-Python exposure time calculator engine
(`pfsspecsim.etc`); `pfs-spec sim` runs the spectral simulator
(`pfsspecsim.sim`). Each subcommand's option surface and behavior is
defined in the sibling `etc.py`/`sim.py` modules -- this module only wires
them into a single discoverable entry point and hosts the one shared
`--version` flag.
"""

from __future__ import annotations

from importlib import metadata

import typer

from .etc import etc_command
from .sim import sim_command

#: Distribution name to look up for `--version` (see `pyproject.toml`).
_PKG_NAME = "pfsspecsim"

app = typer.Typer(
    name="pfs-spec",
    add_completion=False,
    no_args_is_help=True,
    help="PFS spectral tools: exposure time calculator (etc) and simulator (sim).",
)

app.command("etc", help="Run the exposure time calculator and write its ECSV outputs.")(
    etc_command
)
app.command("sim", help="Simulate PFS spectra from an ETC continuum-SNR ECSV file.")(
    sim_command
)


def _package_version() -> str:
    try:
        return metadata.version(_PKG_NAME)
    except metadata.PackageNotFoundError:
        return "0.0.0+unknown"


def _version_callback(value: bool) -> None:
    if value:
        typer.echo(f"pfs-spec {_package_version()}")
        raise typer.Exit()


@app.callback(invoke_without_command=True)
def _root(
    version: bool | None = typer.Option(
        None,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show the pfs-spec version and exit.",
    ),
) -> None:
    """PFS spectral tools: exposure time calculator (etc) and simulator (sim)."""


if __name__ == "__main__":
    app()
