"""Arm-level thread parallelism helper.

Not a port of any `gsetc.c` code -- new infrastructure added on top of the
otherwise-faithful port, so it carries no `gsetc.c:<line>` references.

The three spectrograph arms (b/r/n) are independent, read-only computations
over shared `Spectrograph`/`EtcParams` state, and their dominant cost is
large numpy ufunc/matmul work that releases the GIL -- so a plain
`concurrent.futures.ThreadPoolExecutor` gets real wall-clock speedup without
any multiprocessing pickling/import overhead (see the project plan's package
survey). No new dependency: stdlib only.

**Bit-identity guarantee.** `map_arms` returns results in ``ia=0..n_arms-1``
order regardless of `n_workers` (`ThreadPoolExecutor.map` preserves input
order even though workers may finish out of order). Callers that aggregate
per-arm results (summing into an accumulator, `vstack`-ing per-arm tables,
assigning into a fixed-order `snr_arms[ia]` row) must consume this list in
that same ``ia`` order -- never reduce it with an order-sensitive operation
performed in completion order -- so that floating-point addition order,
and therefore every bit of the output, is identical to the serial
(`n_workers=1`) code path. `n_workers<=1` or `n_arms<=1` short-circuits to a
plain serial list comprehension, which is not just an optimization: it is
the same code path structure as the threaded branch (one `fn(ia)` call per
arm, collected in order), so there is no separate "serial algorithm" to
drift out of sync with the parallel one.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Callable, TypeVar

T = TypeVar("T")


def map_arms(fn: Callable[[int], T], n_arms: int, n_workers: int) -> list[T]:
    """Apply `fn(ia)` for `ia` in `range(n_arms)`, returning results in
    `ia` order.

    Serial (`[fn(ia) for ia in range(n_arms)]`) when `n_workers <= 1` or
    `n_arms <= 1`; otherwise threaded via `ThreadPoolExecutor.map`, which
    preserves input order in its returned iterator regardless of worker
    completion order -- see the module docstring's bit-identity guarantee.
    """
    if n_workers <= 1 or n_arms <= 1:
        return [fn(ia) for ia in range(n_arms)]
    with ThreadPoolExecutor(max_workers=min(n_workers, n_arms)) as ex:
        return list(ex.map(fn, range(n_arms)))
