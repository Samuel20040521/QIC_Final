"""Sprint Task 3: phase-diagram heatmap of Algorithm 1 prediction MAE
over (n, T-count).

X-axis: number of qubits `n`. Y-axis: T-count of the underlying
Clifford+T circuit. Color: mean prediction MAE at fixed query budget
`N`. The figure is annotated with the MAE value in each cell, and the
"safe region" (MAE <= SAFE_THRESHOLD) is outlined.

Resumable: per-cell rows in `experiments/results/raw/11_*.csv`.
"""

from __future__ import annotations

import argparse
import os
import time
from dataclasses import asdict, dataclass

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import psutil
from joblib import Parallel, delayed
from qiskit.quantum_info import Statevector

from qpsq.algorithm import (
    epsilon_tilde,
    gather_and_learn_streaming,
    k_from_epsilon,
)
from qpsq.designs import (
    clifford_plus_t_circuit,
    random_stabilizer_product_state,
)
from qpsq.observables import pauli_z_first
from qpsq.oracle import GaussianQPStat, UnitaryProcess

from _common import FIGS_DIR, RESULTS_DIR
from _hpc import default_max_workers, wait_for_memory_under_cap, wait_for_ram
from _resumable import ResumableCSV

EPS = 0.3
TAU = 0.05
N_BUDGET = 5000
N_TEST = 40
DEFAULT_NS = (4, 6, 8, 10, 12)
DEFAULT_TS = (0, 2, 4, 6, 8, 12, 16)
DEFAULT_INSTANCES = 2
DEFAULT_DEPTH = 3
SAFE_THRESHOLD = 0.10
MEM_CAP_BYTES = 70 * 1024**3

KEY_COLS = ["n", "t_count", "instance"]
ALL_COLS = KEY_COLS + ["mae", "wall_seconds"]


@dataclass
class CellResult:
    n: int
    t_count: int
    instance: int
    mae: float
    wall_seconds: float


def _run_one(n: int, t_count: int, instance: int, seed: int, depth: int) -> CellResult:
    rng = np.random.default_rng(seed)
    U = clifford_plus_t_circuit(n, depth, t_count, rng)
    process = UnitaryProcess(U)
    obs = pauli_z_first(n)
    k = k_from_epsilon(EPS)
    et = epsilon_tilde(EPS, n, k)
    oracle = GaussianQPStat(process, rng=rng)

    test_states = [
        random_stabilizer_product_state(n, rng).to_statevector()
        for _ in range(N_TEST)
    ]
    truths = np.array([process.expectation(s, obs) for s in test_states])

    t0 = time.time()
    model = gather_and_learn_streaming(
        oracle=oracle, n_qubits=n, observable=obs, tau=TAU,
        n_total=N_BUDGET, k=k, eps_tilde=et, rng=rng, chunk_size=500,
    )
    wall = time.time() - t0

    op = model.as_sparse_pauli_op()
    preds = np.array(
        [float(np.real(s.expectation_value(op))) for s in test_states]
    )
    mae = float(np.mean(np.abs(preds - truths)))
    return CellResult(n=n, t_count=t_count, instance=instance,
                      mae=mae, wall_seconds=wall)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--ns", type=str, default=",".join(str(n) for n in DEFAULT_NS))
    p.add_argument("--ts", type=str, default=",".join(str(t) for t in DEFAULT_TS))
    p.add_argument("--instances", type=int, default=DEFAULT_INSTANCES)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--depth", type=int, default=DEFAULT_DEPTH)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    ns = tuple(int(x) for x in args.ns.split(","))
    ts = tuple(int(x) for x in args.ts.split(","))

    out = ResumableCSV(
        path=RESULTS_DIR / "raw" / "11_innovation_heatmap.csv",
        key_cols=KEY_COLS,
        all_cols=ALL_COLS,
    )
    print(f"Resuming with {len(out.done)} cells.")

    pending: list[tuple] = []
    sc = args.seed * 991
    for n in ns:
        for t in ts:
            for inst in range(args.instances):
                sc += 1
                if not out.is_done({"n": n, "t_count": t, "instance": inst}):
                    pending.append((n, t, inst, sc, args.depth))

    print(f"Pending: {len(pending)} cells. Workers: {args.workers}.")
    if pending:
        wait_for_ram(min_fraction=0.20)
        wait_for_memory_under_cap(MEM_CAP_BYTES)
        gen = Parallel(n_jobs=args.workers, return_as="generator", backend="loky", verbose=5)(
            delayed(_run_one)(*a) for a in pending
        )
        for i, result in enumerate(gen):
            out.append(asdict(result))
            if (i + 1) % 4 == 0:
                used = psutil.virtual_memory().used
                if used > MEM_CAP_BYTES:
                    print(f"  RAM at {used / 1024**3:.1f} GiB; pausing 60s")
                    time.sleep(60)

    # ----- aggregate -----
    rows = ResumableCSV.read_all(out.path)
    by_cell: dict[tuple[int, int], list[float]] = {}
    for r in rows:
        nv, tv = int(r["n"]), int(r["t_count"])
        if nv not in ns or tv not in ts:
            continue
        by_cell.setdefault((nv, tv), []).append(float(r["mae"]))

    grid = np.full((len(ts), len(ns)), np.nan)
    for j, n in enumerate(ns):
        for i, t in enumerate(ts):
            vals = by_cell.get((n, t))
            if vals:
                grid[i, j] = float(np.mean(vals))

    # ----- plot -----
    fig, ax = plt.subplots(figsize=(8.5, 6.0), dpi=150)
    vmin = float(np.nanmin(grid)) if np.isfinite(np.nanmin(grid)) else 0.0
    vmax = max(float(np.nanmax(grid)) if np.isfinite(np.nanmax(grid)) else 1.0,
               2 * SAFE_THRESHOLD)
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=SAFE_THRESHOLD, vmax=vmax)
    im = ax.imshow(grid, cmap="RdYlGn_r", norm=norm, aspect="auto", origin="lower")
    ax.set_xticks(range(len(ns)), labels=[str(n) for n in ns])
    ax.set_yticks(range(len(ts)), labels=[str(t) for t in ts])
    ax.set_xlabel("number of qubits n")
    ax.set_ylabel("T-count")
    ax.set_title(
        f"Algorithm 1 hardness diagram   (N = {N_BUDGET}, depth = {args.depth}, "
        f"{args.instances} instance{'s' if args.instances > 1 else ''} per cell)"
    )

    # Annotate
    for i in range(len(ts)):
        for j in range(len(ns)):
            v = grid[i, j]
            if np.isnan(v):
                txt = "--"
            else:
                txt = f"{v:.2f}"
            ax.text(j, i, txt, ha="center", va="center",
                    color="black" if grid[i, j] < SAFE_THRESHOLD * 1.5 else "white",
                    fontsize=9)

    # Outline the "safe region" (MAE <= SAFE_THRESHOLD)
    safe = grid <= SAFE_THRESHOLD
    for i in range(len(ts)):
        for j in range(len(ns)):
            if safe[i, j]:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           fill=False, edgecolor="navy",
                                           linewidth=1.6))

    cbar = fig.colorbar(im, ax=ax, label="prediction MAE")
    cbar.ax.axhline(SAFE_THRESHOLD, color="navy", linewidth=1.5)
    fig.tight_layout()
    out_pdf = FIGS_DIR / "11_innovation_heatmap.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out.path} and {out_pdf}")


if __name__ == "__main__":
    main()
