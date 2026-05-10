"""HPC safety helpers: RAM monitoring, peak-memory tracking, and a
joblib wrapper that respects per-host limits.

Designed for *shared* lab servers where polite resource use matters more
than absolute throughput. Default policy:

  * Workers: at most `os.cpu_count() // 2` parallel jobs.
  * Memory: pause dispatch if available RAM drops below 20%.
  * Per-job memory: each job runs in its own process (joblib's "loky"
    backend), so a leaking job does not poison siblings.

All functions are pure utilities --- they do not mutate global state.
"""

from __future__ import annotations

import os
import threading
import time
from collections.abc import Callable, Iterable
from contextlib import contextmanager
from dataclasses import dataclass
from typing import TypeVar

import psutil

_T = TypeVar("_T")
_R = TypeVar("_R")

# ---------- RAM monitoring -------------------------------------------------


def available_ram_fraction() -> float:
    """Return `available / total` system RAM (in [0, 1])."""
    vm = psutil.virtual_memory()
    return float(vm.available) / float(vm.total)


def wait_for_ram(
    min_fraction: float = 0.20,
    poll_seconds: float = 10.0,
    max_wait_seconds: float = 600.0,
    on_wait: Callable[[float, float], None] | None = None,
) -> None:
    """Block until available RAM rises above `min_fraction`.

    Raises `TimeoutError` after `max_wait_seconds`. Calls `on_wait(frac,
    elapsed)` once per poll if supplied (handy for logging).
    """
    if not 0.0 < min_fraction < 1.0:
        raise ValueError("min_fraction must be in (0, 1)")
    start = time.time()
    while True:
        frac = available_ram_fraction()
        if frac >= min_fraction:
            return
        elapsed = time.time() - start
        if on_wait is not None:
            on_wait(frac, elapsed)
        if elapsed > max_wait_seconds:
            raise TimeoutError(
                f"available RAM stayed below {min_fraction:.0%} for "
                f"{max_wait_seconds:.0f}s (last seen {frac:.0%})"
            )
        time.sleep(poll_seconds)


def default_max_workers(reserve_ratio: float = 0.5) -> int:
    """`os.cpu_count() * reserve_ratio`, clamped to >= 1."""
    cpu = os.cpu_count() or 1
    return max(1, int(cpu * reserve_ratio))


# ---------- Peak memory tracking ------------------------------------------


@dataclass
class MemorySample:
    rss_bytes: int
    elapsed: float


@contextmanager
def peak_memory(poll_seconds: float = 0.25):
    """Context manager that records peak RSS of the current process.

    Usage::

        with peak_memory() as track:
            ... # do work
        print(track.peak_bytes)
    """
    proc = psutil.Process(os.getpid())
    start = time.time()
    state = {"peak": proc.memory_info().rss, "stop": False}

    def _poll() -> None:
        while not state["stop"]:
            try:
                rss = proc.memory_info().rss
            except psutil.NoSuchProcess:
                return
            if rss > state["peak"]:
                state["peak"] = rss
            time.sleep(poll_seconds)

    t = threading.Thread(target=_poll, daemon=True)
    t.start()

    class _Track:
        @property
        def peak_bytes(self) -> int:
            return int(state["peak"])

        @property
        def elapsed_seconds(self) -> float:
            return float(time.time() - start)

    track = _Track()
    try:
        yield track
    finally:
        state["stop"] = True
        t.join(timeout=1.0)


# ---------- joblib wrapper -------------------------------------------------


def safe_parallel_map(
    func: Callable[..., _R],
    args_list: Iterable[tuple],
    *,
    max_workers: int | None = None,
    ram_threshold: float = 0.20,
    backend: str = "loky",
    verbose: int = 0,
) -> list[_R]:
    """Run `func(*args) for args in args_list` across processes safely.

    * `max_workers`: defaults to `default_max_workers()`.
    * `ram_threshold`: do not start the dispatch if available RAM is
      below this fraction.
    * `backend`: joblib backend; `"loky"` keeps each job in its own
      process (memory safety) at modest startup cost.
    """
    from joblib import Parallel, delayed

    if max_workers is None:
        max_workers = default_max_workers()
    if max_workers < 1:
        raise ValueError("max_workers must be >= 1")

    # One-time RAM check before launch.
    if available_ram_fraction() < ram_threshold:
        wait_for_ram(min_fraction=ram_threshold, poll_seconds=15.0)

    items = list(args_list)
    return Parallel(n_jobs=max_workers, backend=backend, verbose=verbose)(
        delayed(func)(*args) for args in items
    )


# ---------- formatting helpers --------------------------------------------


def format_bytes(b: int) -> str:
    if b < 1024:
        return f"{b} B"
    units = ("KiB", "MiB", "GiB", "TiB")
    v = float(b) / 1024
    for u in units:
        if v < 1024:
            return f"{v:.1f} {u}"
        v /= 1024
    return f"{v:.1f} PiB"
