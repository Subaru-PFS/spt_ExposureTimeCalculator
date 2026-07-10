"""Arm-level and chunk-level thread parallelism helpers.

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


def map_index_chunks(
    fn: Callable[[int, int], T], n: int, chunk_size: int, n_workers: int
) -> list[T]:
    """Apply `fn(start, stop)` for each chunk of `range(0, n, chunk_size)`,
    returning results in chunk order.

    Serial (`[fn(start, stop) for ...]`) when `n_workers <= 1` or there is
    only one chunk; otherwise threaded via `ThreadPoolExecutor.map`, which
    -- like `map_arms` -- preserves input order in its returned iterator
    regardless of worker completion order.

    **Bit-identity guarantee.** Each chunk covers a disjoint half-open index
    range `[start, stop)`, so `fn` is expected to compute (and a caller
    assembling the chunks, e.g. via `out[start:stop] = ...` or
    `np.concatenate`, to place) each chunk's output into a region that does
    not overlap any other chunk's. There is no order-sensitive reduction
    across chunks -- the same structure as `map_arms`'s per-arm results --
    so which chunk a worker happens to finish first can never change any
    element of the final result: it is bit-identical to the serial
    (`n_workers=1`) code path for any `n_workers`. `n_workers<=1` or a
    single chunk short-circuits to a plain serial list comprehension, the
    same code path structure as the threaded branch, so there is no
    separate "serial algorithm" to drift out of sync with the parallel one.
    """
    starts = list(range(0, n, chunk_size))
    stops = [min(start + chunk_size, n) for start in starts]
    if n_workers <= 1 or len(starts) <= 1:
        return [fn(start, stop) for start, stop in zip(starts, stops)]
    with ThreadPoolExecutor(max_workers=min(n_workers, len(starts))) as ex:
        return list(ex.map(fn, starts, stops))


def run_products(
    tasks: list[tuple[bool, Callable[[], T]]], n_workers: int
) -> list[T | None]:
    """Run each `fn` in `tasks` for which `enabled` is True, returning a
    list the same length as `tasks` (`None` for entries whose `enabled` is
    False).

    Serial and in `tasks` order (`[fn() if enabled else None for enabled,
    fn in tasks]`) when `n_workers <= 1` or fewer than two entries are
    enabled; otherwise each enabled `fn` is submitted to a
    `ThreadPoolExecutor` (one worker per enabled entry, capped at
    `n_workers`) and collected with `.result()`, which re-raises any
    exception raised inside `fn`.

    **Bit-identity guarantee.** Unlike `map_arms`, nothing here accumulates
    across `tasks` -- each entry computes and returns its own independent
    result (e.g. one whole output `Table`), so worker completion order can
    never affect any individual result's floating-point value. Each `fn`
    closes over whatever `map_arms`/numpy calls it already made when run
    serially, unchanged by this helper -- so every result is bit-identical
    to the `n_workers=1` code path for any `n_workers`.
    """
    if n_workers <= 1 or sum(1 for enabled, _ in tasks if enabled) <= 1:
        return [fn() if enabled else None for enabled, fn in tasks]
    results: list[T | None] = [None] * len(tasks)
    enabled_idx = [i for i, (enabled, _) in enumerate(tasks) if enabled]
    with ThreadPoolExecutor(max_workers=min(n_workers, len(enabled_idx))) as ex:
        futures = [(i, ex.submit(tasks[i][1])) for i in enabled_idx]
        for i, fut in futures:
            results[i] = fut.result()
    return results
