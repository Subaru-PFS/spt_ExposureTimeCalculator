"""Deprecated pre-v2.0 public API compatibility layer.

Reproduces the old subprocess-driven ETC wrapper's surface (`Etc`, its
ALL_CAPS `params` dict, `set_param`/`run`/`make_*`/`run_multi`) on top of
the pure-Python `pfsspecsim.etc` engine, so existing scripts/notebooks keep
working unmodified. May be removed in a future release; new code should use
`pfsspecsim.etc` or `pfsspecsim.sim` directly.
"""
