"""Shared helpers for the experiment scripts."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "experiments" / "results"
FIGS_DIR = REPO_ROOT / "report" / "figs"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGS_DIR.mkdir(parents=True, exist_ok=True)


def parse_seed_only() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def make_rng(seed: int) -> np.random.Generator:
    return np.random.default_rng(seed)
