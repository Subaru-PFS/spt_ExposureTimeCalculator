#!/usr/bin/env python

"""Thin shim for the historical ``python scripts/run_etc.py ...`` invocation.

The real implementation lives in `pfsspecsim.scripts.run_etc` (also the
``pfs-run-etc`` console script installed by the package); this file exists
only so that the old repo-root invocation keeps working, together with the
sibling ``run_etc.defaults`` file. See the README for the current
recommended usage (the ``pfs-spec etc`` CLI).
"""

from pfsspecsim.scripts.run_etc import main

if __name__ == "__main__":
    raise SystemExit(main())
