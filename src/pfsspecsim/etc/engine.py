"""Orchestrator: ports `main()` from `src/gsetc.c`.

Will define `run_etc(params) -> EtcResults` and `run_etc_files(params) ->
EtcResults`, and the `EtcResults` dataclass (noise/snc/snl/oii_curve/
oii_catalog tables plus the two `aperture_factor_800_*` scalars from
gsetc.c:1959-1960). Not yet implemented.
"""
