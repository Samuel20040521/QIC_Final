"""Extension: empirical scaling of Algorithm 1 with system size.

Goal: demonstrate the *quasipolynomial* nature of Algorithm 1's query
complexity by sweeping the number of qubits `n` and recording, for each
`n`:

  * mean prediction MAE on a stabilizer test set as a function of query
    budget `N`,
  * wall-clock training time,
  * peak resident-set memory (process RSS),
  * the size (number of non-zero coefficients) of the learned model.

Memory-safety protocols (see `experiments/_hpc.py`):
  * Workers capped at `os.cpu_count() // 2` (or the user's `--workers`).
  * Available system RAM is checked before parallel dispatch; if below
    20%, the launcher waits.
  * Each job runs in its own process (joblib `loky`); memory is reclaimed
    on completion regardless of any leak in user code.
  * The training step uses `gather_and_learn_streaming` with chunked
    samples and periodic dictionary pruning, so peak per-job memory is
    bounded by the size of `x_hat` --- not by `N x 2^n`.

For evaluation we use Qiskit's *sparse* `Statevector.expectation_value`
path: we never materialise a `2^n x 2^n` matrix. At `n = 14` that would
have been a 4 GiB allocation per evaluation.

CLI:
  uv run python experiments/10_ext_scaling_analysis.py \
      --seed 0 --max-n 10 --workers 4

Outputs:
  experiments/results/ext_scaling_analysis.csv
  report/figs/ext_scaling_analysis.pdf
"""

from __future__ import annotations

import argparse
import csv
import os
import time
from dataclasses import asdict, dataclass

import matplotlib.pyplot as plt
import numpy as np
from qiskit.quantum_info import Statevector

from qpsq.algorithm import (
    epsilon_tilde,
    gather_and_learn_streaming,
    k_from_epsilon,
)
from qpsq.designs import haar_unitary, random_stabilizer_product_state
from qpsq.observables import pauli_z_first
from qpsq.oracle import GaussianQPStat, UnitaryProcess

from _common import FIGS_DIR, RESULTS_DIR
from _hpc import (
    available_ram_fraction,
    default_max_workers,
    format_bytes,
    peak_memory,
    safe_parallel_map,
)

EPS = 0.3
TAU = 0.05
N_TEST = 40
DEFAULT_NS = (4, 6, 8, 10)
DEFAULT_BUDGETS = (500, 1000, 2000, 4000)
N_INSTANCES = 3
CHUNK_SIZE = 500
EPS_TARGET = 0.10  # threshold used for the "min N" extraction


@dataclass
class CellResult:
    n_qubits: int
    instance: int
    n_query: int
    chunk_size: int
    mae: float
    wall_seconds: float
    peak_rss_bytes: int
    n_nonzero_coefficients: int


def _run_single(
    n_qubits: int,
    instance: int,
    n_query: int,
    seed: int,
    chunk_size: int,
) -> CellResult:
    """One worker job: build the process, train Algorithm 1 on `n_query`
    streamed samples, evaluate on a fresh stabilizer test set, return
    a `CellResult`. Runs in its own joblib process for memory isolation.
    """
    rng = np.random.default_rng(seed)
    U = haar_unitary(n_qubits, rng)
    process = UnitaryProcess(U)
    obs = pauli_z_first(n_qubits)
    k = k_from_epsilon(EPS)
    et = epsilon_tilde(EPS, n_qubits, k)
    oracle = GaussianQPStat(process, rng=rng)

    test_states = [
        random_stabilizer_product_state(n_qubits, rng).to_statevector()
        for _ in range(N_TEST)
    ]
    truths = np.array([process.expectation(s, obs) for s in test_states])

    with peak_memory() as track:
        t0 = time.time()
        model = gather_and_learn_streaming(
            oracle=oracle,
            n_qubits=n_qubits,
            observable=obs,
            tau=TAU,
            n_total=n_query,
            k=k,
            eps_tilde=et,
            rng=rng,
            chunk_size=chunk_size,
        )
        wall = time.time() - t0
        peak = track.peak_bytes

    op = model.as_sparse_pauli_op()
    preds = np.array(
        [float(np.real(s.expectation_value(op))) for s in test_states]
    )
    mae = float(np.mean(np.abs(preds - truths)))
    n_nonzero = sum(1 for v in model.coefficients.values() if v != 0.0)

    return CellResult(
        n_qubits=n_qubits,
        instance=instance,
        n_query=n_query,
        chunk_size=chunk_size,
        mae=mae,
        wall_seconds=wall,
        peak_rss_bytes=peak,
        n_nonzero_coefficients=n_nonzero,
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--max-n", type=int, default=10,
                   help="Largest n to include in the sweep (default 10).")
    p.add_argument("--ns", type=str, default=None,
                   help="Comma-separated n values (overrides --max-n).")
    p.add_argument("--budgets", type=str, default=None,
                   help="Comma-separated query budgets (default 500,1000,2000,4000).")
    p.add_argument("--instances", type=int, default=N_INSTANCES)
    p.add_argument("--chunk-size", type=int, default=CHUNK_SIZE)
    p.add_argument("--workers", type=int, default=None,
                   help="Parallel workers (default os.cpu_count() // 2).")
    p.add_argument("--ram-threshold", type=float, default=0.20)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.ns:
        ns = tuple(int(x) for x in args.ns.split(","))
    else:
        ns = tuple(n for n in DEFAULT_NS if n <= args.max_n)
        if 12 <= args.max_n:
            ns = ns + (12,) if 12 not in ns else ns
        if 14 <= args.max_n:
            ns = ns + (14,) if 14 not in ns else ns

    if args.budgets:
        budgets = tuple(int(x) for x in args.budgets.split(","))
    else:
        budgets = DEFAULT_BUDGETS

    workers = args.workers or default_max_workers()

    print("=" * 72)
    print(f"Scaling sweep: n in {ns}, budgets {budgets}, instances={args.instances}")
    print(f"Workers: {workers} (cpu_count={os.cpu_count()})")
    print(f"Available RAM: {available_ram_fraction():.1%}; threshold={args.ram_threshold:.0%}")
    print(f"Streaming chunk size: {args.chunk_size}")
    print("=" * 72)

    # Build the (n, instance, N, seed) job list.
    job_args: list[tuple[int, int, int, int, int]] = []
    seed_counter = args.seed * 1_000_003
    for n in ns:
        for inst in range(args.instances):
            for N in budgets:
                seed_counter += 1
                job_args.append((n, inst, N, seed_counter, args.chunk_size))

    # We submit jobs sorted by (n_qubits ascending, n_query ascending).
    # joblib will then pick them up in roughly that order.
    job_args.sort(key=lambda t: (t[0], t[2], t[1]))

    t_overall = time.time()
    results: list[CellResult] = safe_parallel_map(
        _run_single,
        job_args,
        max_workers=workers,
        ram_threshold=args.ram_threshold,
        verbose=5,
    )
    t_overall = time.time() - t_overall
    print(f"\nFinished {len(results)} jobs in {t_overall:.1f}s.")

    # CSV
    csv_path = RESULTS_DIR / "ext_scaling_analysis.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        if results:
            w.writerow(list(asdict(results[0]).keys()))
            for r in results:
                w.writerow(list(asdict(r).values()))

    # ----- Aggregations -----
    by_nN: dict[tuple[int, int], list[CellResult]] = {}
    for r in results:
        by_nN.setdefault((r.n_qubits, r.n_query), []).append(r)

    def _agg(field: str) -> dict[tuple[int, int], tuple[float, float]]:
        out: dict[tuple[int, int], tuple[float, float]] = {}
        for key, rs in by_nN.items():
            vals = np.array([getattr(r, field) for r in rs], dtype=float)
            out[key] = (float(vals.mean()), float(vals.std()))
        return out

    mae_agg = _agg("mae")
    wall_agg = _agg("wall_seconds")
    peak_agg = _agg("peak_rss_bytes")
    nz_agg = _agg("n_nonzero_coefficients")

    # ----- Plot -----
    fig = plt.figure(figsize=(13, 9))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.28)

    # (1) MAE vs N for each n.
    ax1 = fig.add_subplot(gs[0, 0])
    cmap = plt.get_cmap("viridis")
    for i, n in enumerate(ns):
        xs, ys, ystd = [], [], []
        for N in budgets:
            if (n, N) in mae_agg:
                m, s = mae_agg[(n, N)]
                xs.append(N); ys.append(m); ystd.append(s)
        ax1.errorbar(xs, ys, yerr=ystd, marker="o",
                     color=cmap(0.1 + 0.8 * i / max(1, len(ns) - 1)),
                     label=f"n={n}", capsize=3)
    ax1.axhline(EPS_TARGET, color="k", linestyle=":", linewidth=1)
    ax1.set_xscale("log")
    ax1.set_xlabel("query budget N")
    ax1.set_ylabel("MAE")
    ax1.set_title("Convergence: MAE vs N for each n")
    ax1.grid(True, which="both", linestyle=":", alpha=0.4)
    ax1.legend(loc="upper right", fontsize=8)

    # (2) Min N to reach EPS_TARGET vs n.
    ax2 = fig.add_subplot(gs[0, 1])
    min_N_per_n: list[float | None] = []
    for n in ns:
        ns_pairs = sorted(((N, mae_agg.get((n, N), (np.nan, 0))[0]) for N in budgets),
                           key=lambda x: x[0])
        hit = next((N for N, mae in ns_pairs if mae <= EPS_TARGET and not np.isnan(mae)), None)
        min_N_per_n.append(float(hit) if hit is not None else np.nan)
    valid = [(n, N) for n, N in zip(ns, min_N_per_n) if N is not None and not np.isnan(N)]
    has_legend_2 = False
    if valid:
        xs = [v[0] for v in valid]; ys = [v[1] for v in valid]
        ax2.plot(xs, ys, marker="o", linewidth=1.5, color="C0", label="empirical")
        has_legend_2 = True
        if len(xs) >= 2:
            n_min = min(xs)
            ref_const = ys[0]
            ax2.plot(xs, [ref_const * (xx / n_min) ** 2 for xx in xs],
                     "k--", linewidth=0.8, label=r"$n^2$ reference")
            ax2.plot(xs, [ref_const * (xx / n_min) ** 3 for xx in xs],
                     "k:", linewidth=0.8, label=r"$n^3$ reference")
    ax2.set_yscale("log")
    ax2.set_xlabel("n")
    ax2.set_ylabel(f"min N to reach MAE <= {EPS_TARGET}")
    ax2.set_title("Sample complexity vs n")
    ax2.grid(True, which="both", linestyle=":", alpha=0.4)
    if has_legend_2:
        ax2.legend(loc="upper left", fontsize=8)

    # (3) Wall time vs N for each n.
    ax3 = fig.add_subplot(gs[1, 0])
    for i, n in enumerate(ns):
        xs, ys, ystd = [], [], []
        for N in budgets:
            if (n, N) in wall_agg:
                m, s = wall_agg[(n, N)]
                xs.append(N); ys.append(m); ystd.append(s)
        ax3.errorbar(xs, ys, yerr=ystd, marker="s",
                     color=cmap(0.1 + 0.8 * i / max(1, len(ns) - 1)),
                     label=f"n={n}", capsize=3)
    ax3.set_xscale("log"); ax3.set_yscale("log")
    ax3.set_xlabel("query budget N")
    ax3.set_ylabel("wall-clock time (s)")
    ax3.set_title("Wall-clock time per training run")
    ax3.grid(True, which="both", linestyle=":", alpha=0.4)
    ax3.legend(loc="lower right", fontsize=8)

    # (4) Peak memory vs n at the largest budget.
    ax4 = fig.add_subplot(gs[1, 1])
    largest_N = max(budgets)
    xs, ys, ystd = [], [], []
    for n in ns:
        if (n, largest_N) in peak_agg:
            m, s = peak_agg[(n, largest_N)]
            xs.append(n); ys.append(m); ystd.append(s)
    ax4.errorbar(xs, [y / 1e6 for y in ys],
                 yerr=[s / 1e6 for s in ystd],
                 marker="^", color="C3", capsize=3)
    ax4.set_xlabel("n")
    ax4.set_ylabel(f"peak RSS at N = {largest_N} (MB)")
    ax4.set_title("Per-job peak memory")
    ax4.grid(True, which="both", linestyle=":", alpha=0.4)

    fig.suptitle(
        f"Algorithm 1 scaling sweep "
        f"(eps={EPS}, tau={TAU}, instances={args.instances}, chunk={args.chunk_size})"
    )
    out = FIGS_DIR / "ext_scaling_analysis.pdf"
    fig.savefig(out, bbox_inches="tight")
    print(f"wrote {csv_path} and {out}")

    # Print summary.
    print("\n--- Summary table (mean over instances) ---")
    print(f"{'n':>3}  {'N':>5}  {'MAE':>7}  {'wall':>7}  {'peak_RSS':>10}  {'#alpha':>6}")
    for n in ns:
        for N in budgets:
            if (n, N) in mae_agg:
                m_mae, _ = mae_agg[(n, N)]
                m_wall, _ = wall_agg[(n, N)]
                m_peak, _ = peak_agg[(n, N)]
                m_nz, _ = nz_agg[(n, N)]
                print(f"{n:>3}  {N:>5}  {m_mae:>7.4f}  {m_wall:>7.2f}  "
                      f"{format_bytes(int(m_peak)):>10}  {int(m_nz):>6}")


if __name__ == "__main__":
    main()
