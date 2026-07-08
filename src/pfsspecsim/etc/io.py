"""Astropy ECSV writer/reader for `EtcResults` tables.

Implements the "ECSV 出力スキーマ" section of the task brief: every table
`run_etc_files` writes carries a common `table.meta` block (package
version, UTC creation timestamp, resolved spectrograph config filename,
and every resolved `EtcParams` field -- `Path` values stringified, `None`
values passed through as-is since ECSV's YAML-header meta block already
renders a Python `None` as `null`). Column units (`nm`, `m**2`,
`electron`/`electron**2`, ...) are attached here (and nowhere else --
`engine.py`'s computational kernels stay in plain `ndarray`s, per
`constants.py`'s unit discipline note) immediately before writing.

`overwrite=False` protection: this module's `write_table` passes
`overwrite` straight through to `astropy.table.Table.write`, which raises
`OSError` at write time if the file already exists. The task brief's
"raise before running the (possibly expensive) computation" behavior
(a port of `pfsspecsim/pfsetc.py:236-249`'s existence pre-check) is
implemented by `engine.run_etc_files`, one level up -- by the time this
module's `write_table` is called, that decision has already been made.
"""

from __future__ import annotations

import dataclasses
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

from astropy.table import Table

from .params import EtcParams

#: Distribution name to look up for `etc_version` (see `pyproject.toml`).
_PKG_NAME = "pfsspecsim"


def _package_version() -> str:
    """Installed `pfsspecsim` version, or a fallback if run from an
    uninstalled checkout (e.g. `python -m pytest` without `pip install -e`,
    which should not happen under `uv run` but costs nothing to guard).
    """
    try:
        return metadata.version(_PKG_NAME)
    except metadata.PackageNotFoundError:
        return "0.0.0+unknown"


def _params_to_meta(params: EtcParams) -> dict:
    """`EtcParams` as a plain, YAML/ECSV-meta-safe dict: every
    `dataclasses.fields(EtcParams)` entry, with `Path` values stringified.
    `None` values are left as `None` (ECSV's YAML meta header already
    renders that as `null` natively; no special-casing needed).
    """
    out: dict = {}
    for f in dataclasses.fields(EtcParams):
        value = getattr(params, f.name)
        out[f.name] = str(value) if isinstance(value, Path) else value
    return out


def build_meta(params: EtcParams, instr_config: str | Path) -> dict:
    """Build the common `table.meta` block attached to every output ECSV
    table (task brief's ECSV schema section).

    Parameters
    ----------
    params : EtcParams
        The resolved input parameters for this run.
    instr_config : str or Path
        The spectrograph config file actually used (only its filename is
        recorded, not the full path -- packaged config files live under an
        installation-specific `pfsspecsim/config/` directory that isn't
        meaningful to persist).

    Returns
    -------
    dict
        ``{"etc_version": ..., "created": ..., "instr_config": ...,
        "params": {...}}``.
    """
    return {
        "etc_version": _package_version(),
        "created": datetime.now(timezone.utc).isoformat(),
        "instr_config": Path(instr_config).name,
        "params": _params_to_meta(params),
    }


def write_table(
    table: Table,
    path: str | Path,
    params: EtcParams,
    instr_config: str | Path,
    overwrite: bool = True,
) -> None:
    """Write one `EtcResults` table to `path` as ECSV, with `build_meta`'s
    standard meta block attached.

    Creates `path`'s parent directory if it doesn't already exist. Does
    not mutate the caller's `table` (writes a shallow copy with the meta
    block merged in, so `table.meta` -- typically empty, from the plain
    `Table(...)` constructor calls in `engine.py` -- is left untouched).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = table.copy()
    out.meta.update(build_meta(params, instr_config))
    out.write(path, format="ascii.ecsv", overwrite=overwrite)


def read_table(path: str | Path) -> Table:
    """Read back an ECSV table written by :func:`write_table` (also used
    by `engine.run_etc` to reload the noise vector when
    `EtcParams.noise_reused` is set).
    """
    return Table.read(Path(path), format="ascii.ecsv")
