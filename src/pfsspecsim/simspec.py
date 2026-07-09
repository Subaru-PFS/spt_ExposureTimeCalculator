"""Modern snake_case Python API for the spectral simulator.

Mirrors `pfsspecsim.etc.params`'s architecture: `SimSpecParams` is a flat,
picklable dataclass of every `pfsspecsim.pfsspec.Pfsspec.params` entry
(replacing the legacy mixed-case `EXP_NUM`/`MAG_FILE`/`outDir`/... dict keys
with predictable snake_case field names), and `load_params` builds one with
CLI-overrides > TOML > dataclass-default priority (same contract as
`pfsspecsim.etc.params.load_params`). See `run_sim_spec` (below) for how a
`SimSpecParams` drives the (unchanged) `Pfsspec` class.

`mag`/`mag_file` are split into two mutually-exclusive fields, exactly as
`EtcParams.mag`/`mag_file` split the single legacy `MAG_FILE` parameter (see
`pfsspecsim.etc.params.EtcParams`): both still multiplex onto the single
legacy `Pfsspec.params["MAG_FILE"]` key underneath.
"""

from __future__ import annotations

import dataclasses
import tomllib
from pathlib import Path
from typing import Any

from .pfsspec import Pfsspec

#: Fields whose values are filesystem paths; coerced to `Path` by
#: `load_params` (mirrors `pfsspecsim.etc.params._PATH_FIELDS`).
_PATH_FIELDS = frozenset({"etc_file", "mag_file", "out_dir"})


@dataclasses.dataclass
class SimSpecParams:
    """Resolved spectral-simulator input parameters (one flat, picklable,
    snake_case set), the modern counterpart to `Pfsspec.params`.

    Defaults reproduce `Pfsspec.__init__`'s built-in defaults.
    """

    etc_file: Path = Path("out/ref.snc.dat")
    exp_num: int = 4
    mag: float | None = 22.5  # AB mag (flat); XOR with mag_file
    mag_file: Path | None = None  # wavelength[nm]/mag file, or multi-object table
    counts_min: float = 0.1
    nrealize: int = 1
    out_dir: Path = Path("out")
    ascii_table: str | None = None
    ra: float = 150.0
    dec: float = 2.0
    tract: int = 0
    patch: str = "0,0"
    visit0: int = 1
    cat_id: int = 0
    obj_id: int = 1
    fiber_id: int = 1
    fiber_mag: list[float] = dataclasses.field(
        default_factory=lambda: [22.5, 22.5, 22.5, 22.5, 22.5]
    )
    filter_name: list[str] = dataclasses.field(
        default_factory=lambda: ["hcs_g", "hcs_r", "hcs_i", "hcs_z", "hcs_y"]
    )
    # Spectrograph unit number (1-4); unrelated to `EtcParams.spectrograph`,
    # which is an 'ave'/'sm1'..'sm4' throughput-model tag.
    spectrograph: int = 1
    pfs_config_full: bool = False
    write_fits: bool = True
    write_pfs_arm: bool = True
    plot_arm_set: bool = False
    plot_object: bool = False
    sky_sub_floor: float = 0.01
    sky_sub_mode: str = "random"
    sky_sub_seed: int = 0

    def validate(self) -> None:
        """Cross-field checks; raises `ValueError` on failure."""
        if (self.mag is None) == (self.mag_file is None):
            raise ValueError(
                "Exactly one of `mag` (flat AB magnitude) or `mag_file` "
                "(wavelength-dependent magnitude spectrum, or multi-object "
                "magnitude table) must be set -- they are mutually "
                f"exclusive, got mag={self.mag!r} mag_file={self.mag_file!r}"
            )


def load_params(
    toml_path: str | Path | None = None,
    overrides: dict[str, Any] | None = None,
) -> SimSpecParams:
    """Build a `SimSpecParams`, merging defaults, an optional TOML file, and
    explicit overrides, with priority CLI (`overrides`) > TOML > defaults.

    `toml_path` is a flat TOML file with the same snake_case keys as
    `SimSpecParams` fields. Unknown keys, in either the TOML file or
    `overrides`, raise `ValueError`.
    """
    field_names = {f.name for f in dataclasses.fields(SimSpecParams)}
    values: dict[str, Any] = {}

    if toml_path is not None:
        with open(toml_path, "rb") as fh:
            toml_data = tomllib.load(fh)
        unknown = set(toml_data) - field_names
        if unknown:
            raise ValueError(
                f"Unknown key(s) in TOML file {toml_path}: {sorted(unknown)}"
            )
        values.update(toml_data)

    if overrides:
        unknown_overrides = set(overrides) - field_names
        if unknown_overrides:
            raise ValueError(f"Unknown override key(s): {sorted(unknown_overrides)}")
        values.update(overrides)

    # If exactly one of `mag`/`mag_file` was actually supplied (by the TOML
    # file and/or overrides) and the other was not mentioned at all, force
    # the unmentioned one to `None` -- otherwise it would silently fall back
    # to `SimSpecParams`'s own default (`mag=22.5`), spuriously colliding
    # with the one the caller *did* specify and tripping `validate()`'s XOR
    # check. This mirrors what `pfsspecsim.cli.main` already does for
    # CLI-flag-level mag/mag_file exclusivity, but applies it uniformly to
    # TOML-only settings too (which `main`'s CLI-flag-only logic misses).
    if ("mag" in values) != ("mag_file" in values):
        if "mag" in values:
            values["mag_file"] = None
        else:
            values["mag"] = None

    for name in _PATH_FIELDS:
        if values.get(name) is not None:
            values[name] = Path(values[name])

    params = SimSpecParams(**values)
    params.validate()
    return params


#: snake_case `SimSpecParams` field name -> legacy `Pfsspec.params` dict
#: key. `mag`/`mag_file` both target the single legacy `MAG_FILE` key --
#: `run_sim_spec` sends only whichever one of the pair is not `None`.
_LEGACY_KEY_MAP: dict[str, str] = {
    "etc_file": "etcFile",
    "exp_num": "EXP_NUM",
    "mag": "MAG_FILE",
    "mag_file": "MAG_FILE",
    "counts_min": "countsMin",
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
#: (a `"1"/"t"/"true"` vs `"0"/"f"/"false"` string convention, case
#: insensitive) rather than a native Python bool.
_STRTOBOOL_FIELDS = frozenset({"write_fits", "write_pfs_arm", "pfs_config_full"})


def run_sim_spec(params: SimSpecParams) -> Pfsspec:
    """Run the spectral simulator for `params` and return the populated
    `Pfsspec` instance (`sim.pfsObjects`, `sim.outdir`, `sim.asciiTable`,
    ... -- the same attributes `Pfsspec.make_sim_spec` has always set).

    Translates `params` to `Pfsspec.set_param` calls and invokes
    `Pfsspec.make_sim_spec()`; `Pfsspec`'s own internals (FITS/datamodel
    output path included) are untouched by this function.
    """
    params.validate()
    sim = Pfsspec()
    for f in dataclasses.fields(params):
        name = f.name
        if name == "mag" and params.mag is None:
            continue
        if name == "mag_file" and params.mag_file is None:
            continue
        value = getattr(params, name)
        if name == "ascii_table":
            value = "None" if value is None else value
        elif name in _STRTOBOOL_FIELDS:
            value = "true" if value else "false"
        sim.set_param(_LEGACY_KEY_MAP[name], value)
    sim.make_sim_spec()
    return sim
