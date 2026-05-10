"""Append-only CSV that resumes after interruption.

Usage:

    log = ResumableCSV(
        path=Path("experiments/results/raw/foo.csv"),
        key_cols=["n", "t_count", "instance"],
        all_cols=["n", "t_count", "instance", "mae", "wall_seconds"],
    )
    if not log.is_done(key):
        result = run_one(key)
        log.append(result)

The header row is written exactly once. Subsequent runs read the existing
file, populate the `done` set, and skip already-completed keys.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


class ResumableCSV:
    def __init__(self, path: Path, key_cols: list[str], all_cols: list[str]) -> None:
        if not set(key_cols).issubset(all_cols):
            raise ValueError("key_cols must be a subset of all_cols")
        self.path = Path(path)
        self.key_cols = list(key_cols)
        self.all_cols = list(all_cols)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.done: set[tuple[str, ...]] = set()
        if self.path.exists():
            with self.path.open() as f:
                reader = csv.DictReader(f)
                for row in reader:
                    self.done.add(tuple(row[k] for k in self.key_cols))
        else:
            with self.path.open("w", newline="") as f:
                csv.writer(f).writerow(self.all_cols)

    def _key_tuple(self, key: dict[str, Any]) -> tuple[str, ...]:
        return tuple(str(key[k]) for k in self.key_cols)

    def is_done(self, key: dict[str, Any]) -> bool:
        return self._key_tuple(key) in self.done

    def append(self, row: dict[str, Any]) -> None:
        for c in self.all_cols:
            if c not in row:
                raise ValueError(f"missing column {c}")
        with self.path.open("a", newline="") as f:
            csv.writer(f).writerow([row[c] for c in self.all_cols])
        self.done.add(self._key_tuple(row))

    @classmethod
    def read_all(cls, path: Path) -> list[dict[str, str]]:
        if not Path(path).exists():
            return []
        with Path(path).open() as f:
            return list(csv.DictReader(f))
