"""Tests for pfsspecsim.etc._parallel.

Not a port of any `gsetc.c` code -- see `_parallel.py`'s module docstring.
`map_arms` is already exercised indirectly (bit-identity, via
`test_engine.py::TestArmParallelism`); this file covers `map_index_chunks`
directly: serial-vs-threaded equality and its chunk-boundary edge cases
(n divisible by chunk_size, remainder, single chunk, n_workers=1).
"""

from __future__ import annotations

import threading

import pytest

from pfsspecsim.etc._parallel import map_index_chunks


def _double(start: int, stop: int) -> list[int]:
    return [i * 2 for i in range(start, stop)]


class TestMapIndexChunks:
    @pytest.mark.parametrize(
        "n,chunk_size",
        [
            (10, 5),  # n divisible by chunk_size: exactly 2 chunks
            (10, 3),  # remainder: chunks of 3,3,3,1
            (5, 100),  # chunk_size > n: a single chunk
            (0, 5),  # empty range: no chunks at all
        ],
    )
    @pytest.mark.parametrize("n_workers", [1, 2, 4])
    def test_matches_serial_reference(self, n, chunk_size, n_workers):
        got = map_index_chunks(_double, n, chunk_size, n_workers)
        expected = [
            _double(start, min(start + chunk_size, n))
            for start in range(0, n, chunk_size)
        ]
        assert got == expected

    def test_single_chunk_is_serial_regardless_of_n_workers(self):
        # n <= chunk_size -> exactly one chunk -> always the serial code
        # path (see the docstring's "single chunk" short-circuit), so the
        # calling thread must be the one that ran fn.
        main_thread = threading.current_thread()
        seen = []

        def _record(start: int, stop: int) -> None:
            seen.append(threading.current_thread())

        map_index_chunks(_record, 5, 100, n_workers=4)
        assert seen == [main_thread]

    def test_n_workers_le_1_is_serial(self):
        main_thread = threading.current_thread()
        seen = []

        def _record(start: int, stop: int) -> None:
            seen.append(threading.current_thread())

        map_index_chunks(_record, 20, 5, n_workers=1)
        assert seen == [main_thread] * 4

    def test_threaded_path_actually_uses_worker_threads(self):
        # Sanity check that n_workers>1 with >1 chunk really does dispatch
        # off the calling thread for at least one chunk (not a strict
        # requirement of the bit-identity guarantee, but confirms the
        # threaded branch is reachable/exercised).
        main_thread = threading.current_thread()
        seen = []

        def _record(start: int, stop: int) -> None:
            seen.append(threading.current_thread())

        map_index_chunks(_record, 20, 5, n_workers=4)
        assert len(seen) == 4
        assert any(t is not main_thread for t in seen)

    def test_chunks_are_disjoint_and_cover_the_full_range(self):
        n, chunk_size = 23, 4
        chunks = map_index_chunks(
            lambda start, stop: list(range(start, stop)), n, chunk_size, n_workers=4
        )
        flat = [i for chunk in chunks for i in chunk]
        assert flat == list(range(n))
