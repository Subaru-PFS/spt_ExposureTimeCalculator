"""`run_sim_spec`: translates a `SimSpecParams` into a `Pfsspec` run.

Maps `SimSpecParams`'s snake_case fields onto the legacy `Pfsspec.params`
ALL_CAPS/camelCase dict keys (`_LEGACY_KEY_MAP`) and legacy string
conventions (`_STRTOBOOL_FIELDS`), then drives the (unchanged)
`pfsspecsim.sim.pfsspec.Pfsspec` class -- see `pfsspecsim.sim.params` for
the `SimSpecParams` dataclass and `load_params` itself.
"""

from __future__ import annotations

import dataclasses

from .params import SimSpecParams
from .pfsspec import Pfsspec

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
