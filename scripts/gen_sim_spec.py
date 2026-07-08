#!/usr/bin/env python

"""Thin shim for the historical ``python scripts/gen_sim_spec.py ...``
invocation.

The real implementation lives in `pfsspecsim.scripts.gen_sim_spec` (also
the ``pfs-gen-sim-spec`` console script installed by the package); this
file exists only so that the old repo-root invocation keeps working,
together with the sibling ``gen_sim_spec.defaults`` file. See the README
for the current recommended usage (the ``pfs-sim-spec`` CLI).
"""

from pfsspecsim.scripts.gen_sim_spec import main

if __name__ == "__main__":
    raise SystemExit(main())
