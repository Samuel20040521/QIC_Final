"""Sprint Task 1: T-count hardness boundary at n in {8, 10, 12}.

For each (n, t_count, instance, N) cell we run Algorithm 1 against a fresh
Clifford+T circuit and record MAE on a stabilizer test set. We then
report, per (n, t_count), the smallest N at which the average MAE drops
to or below `MAE_TARGET`.

Resumable: per-cell rows live in `experiments/results/raw/12_*.csv` and
are skipped on restart. Memory-safe: max-workers cap, 70 GiB system-RAM
guard, and the streaming `gather_and_learn_streaming` learner.
"""

from __future__ import annotations

import argparse
import os
import time
from dataclasses import asdict, dataclass

import matplotlib.pyplot as plt
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
from _hpc import (
    available_ram_fraction,
    default_max_workers,
    wait_for_memory_under_cap,
    wait_for_ram,
)
from _resumable import ResumableCSV

EPS = 0.3
TAU = 0.05
N_TEST = 40
DEFAULT_NS = (8, 10, 12)
DEFAULT_TS = (0, 2, 4, 6, 8, 10, 12, 14, 16)
DEFAULT_BUDGETS = (1000, 2000, 4000, 8000)
DEFAULT_INSTANCES = 3
DEFAULT_DEPTH = 3
MAE_TARGET = 0.08
MEM_CAP_BYTES = 70 * 1024**3  # 70 GiB

KEY_COLS = ["n", "t_count", "instance", "n_query"]
ALL_COLS = KEY_COLS + ["mae", "wall_seconds", "n_nonzero"]


@dataclass
class CellResult:
    n: int
    t_count: int
    instance: int
    n_query: int
    mae: float
    wall_seconds: float
    n_nonzero: int


def _run_one(
    n: int,
    t_count: int,
    instance: int,
    n_query: int,
    seed: int,
    depth: int,
) -> CellResult:
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
        oracle=oracle,
        n_qubits=n,
        observable=obs,
        tau=TAU,
        n_total=n_query,
        k=k,
        eps_tilde=et,
        rng=rng,
        chunk_size=500,
    )
    wall = time.time() - t0

    op = model.as_sparse_pauli_op()
    preds = np.array(
        [float(np.real(s.expectation_value(op))) for s in test_states]
    )
    mae = float(np.mean(np.abs(preds - truths)))
    n_nz = sum(1 for v in model.coefficients.values() if v != 0.0)
    return CellResult(n=n, t_count=t_count, instance=instance,
                      n_query=n_query, mae=mae, wall_seconds=wall, n_nonzero=n_nz)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--ns", type=str, default=",".join(str(n) for n in DEFAULT_NS))
    p.add_argument("--ts", type=str, default=",".join(str(t) for t in DEFAULT_TS))
    p.add_argument("--budgets", type=str,
                   default=",".join(str(b) for b in DEFAULT_BUDGETS))
    p.add_argument("--instances", type=int, default=DEFAULT_INSTANCES)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--depth", type=int, default=DEFAULT_DEPTH)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    ns = tuple(int(x) for x in args.ns.split(","))
    ts = tuple(int(x) for x in args.ts.split(","))
    budgets = tuple(int(x) for x in args.budgets.split(","))

    out = ResumableCSV(
        path=RESULTS_DIR / "raw" / "12_t_count_hardness.csv",
        key_cols=KEY_COLS,
        all_cols=ALL_COLS,
    )
    print(f"Resuming with {len(out.done)} cells already done.")

    seed_counter = args.seed * 31_337
    pending: list[tuple] = []
    for n in ns:
        for t in ts:
            for inst in range(args.instances):
                for N in budgets:
                    seed_counter += 1
                    key = {"n": n, "t_count": t, "instance": inst, "n_query": N}
                    if out.is_done(key):
                        continue
                    pending.append((n, t, inst, N, seed_counter, args.depth))

    print(f"Pending: {len(pending)} cells. Workers: {args.workers}.")
    print(f"Available RAM: {available_ram_fraction():.1%}, "
          f"system used: {psutil.virtual_memory().used / 1024**3:.1f} GiB.")

    if pending:
        wait_for_ram(min_fraction=0.20)
        wait_for_memory_under_cap(MEM_CAP_BYTES)

        t_start = time.time()
        gen = Parallel(n_jobs=args.workers, return_as="generator", backend="loky", verbose=5)(
            delayed(_run_one)(*args_) for args_ in pending
        )
        for i, result in enumerate(gen):
            out.append(asdict(result))
            if (i + 1) % 6 == 0:
                # Periodic memory check between batches.
                used = psutil.virtual_memory().used
                if used > MEM_CAP_BYTES:
                    print(f"  RAM at {used / 1024**3:.1f} GiB > 70 GiB; pausing 60s")
                    time.sleep(60)
        print(f"Completed in {time.time() - t_start:.1f}s.")

    # ----- aggregate + plot -----
    rows = ResumableCSV.read_all(out.path)
    by_cell: dict[tuple[int, int, int, int], float] = {}
    for r in rows:
        key = (int(r["n"]), int(r["t_count"]), int(r["instance"]), int(r["n_query"]))
        by_cell[key] = float(r["mae"])

    # Average MAE over instances at each (n, t, N).
    avg_mae: dict[tuple[int, int, int], float] = {}
    for n in ns:
        for t in ts:
            for N in budgets:
                vals = [by_cell.get((n, t, i, N))
                        for i in range(args.instances)
                        if by_cell.get((n, t, i, N)) is not None]
                if vals:
                    avg_mae[(n, t, N)] = float(np.mean(vals))

    # min N to reach MAE_TARGET.
    min_N: dict[tuple[int, int], float | None] = {}
    for n in ns:
        for t in ts:
            ns_pairs = sorted(
                ((N, avg_mae.get((n, t, N))) for N in budgets), key=lambda x: x[0]
            )
            hit = next(
                (N for N, m in ns_pairs if m is not None and m <= MAE_TARGET),
                None,
            )
            min_N[(n, t)] = float(hit) if hit is not None else None

    # ----- plot -----
    fig, ax = plt.subplots(figsize=(8.5, 5.0), dpi=150)
    cmap = plt.get_cmap("plasma")
    for i, n in enumerate(ns):
        xs, ys = [], []
        for t in ts:
            v = min_N.get((n, t))
            if v is not None:
                xs.append(t); ys.append(v)
        color = cmap(0.15 + 0.7 * i / max(1, len(ns) - 1))
        ax.plot(xs, ys, marker="o", color=color, linewidth=1.6,
                label=f"n = {n}")
    # log_2(n) reference markers
    ymin = min((v for v in min_N.values() if v is not None), default=1)
    ymax = max((v for v in min_N.values() if v is not None), default=10)
    for n in ns:
        ax.axvline(np.log2(n), color="k", linestyle=":", linewidth=0.5, alpha=0.4)
    ax.set_xlabel("T-count")
    ax.set_ylabel(f"min N to reach mean MAE <= {MAE_TARGET}")
    ax.set_yscale("log")
    ax.set_title(
        f"T-count hardness boundary  (depth={args.depth}, "
        f"{args.instances} instances per cell, eps={EPS}, tau={TAU})"
    )
    ax.grid(True, which="both", linestyle=":", alpha=0.4)
    ax.legend(loc="lower right", fontsize=9)
    fig.tight_layout()
    out_pdf = FIGS_DIR / "12_t_count_hardness.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out.path} and {out_pdf}")


if __name__ == "__main__":
    main()
