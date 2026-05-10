"""Extension: empirical sample complexity of Algorithm 1 vs. T-count.

The Wadhwa-Doosti diamond-distance lower bounds (Theorems 2 and 4)
characterize two extreme ensembles --- exact 2-designs (Cliffords are a
3-design and inherit the bound) and Haar-random unitaries --- but say
nothing about the regime in between, where T-count parameterises
non-stabilizer structure.

We run Algorithm 1 on Clifford+T circuits for a sweep of T-count values
and measure prediction MAE at growing query budgets. We compare to a
classical-shadow baseline that takes `Q / T_test` randomized Pauli
snapshots per test state. The total query budget `Q` is matched between
the two methods.

Outputs:
  experiments/results/ext_t_count_complexity.csv
  report/figs/ext_t_count_complexity.pdf
"""

from __future__ import annotations

import csv
import time

import matplotlib.pyplot as plt
import numpy as np
from qiskit.quantum_info import DensityMatrix

from qpsq.algorithm import (
    Algorithm1Result,
    epsilon_tilde,
    gather,
    k_from_epsilon,
    learn,
)
from qpsq.designs import (
    clifford_plus_t_circuit,
    random_stabilizer_product_state,
)
from qpsq.observables import pauli_l1_norm, pauli_z_first
from qpsq.oracle import GaussianQPStat, UnitaryProcess
from qpsq.shadows import shadow_evolution_estimator

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

N_QUBITS = 5
N_INSTANCES = 4
DEPTH = 3
T_COUNTS = (0, 1, 2, 4, 8, 16)
QUERY_BUDGETS = (500, 1000, 2000, 4000, 8000)
EPS = 0.3
TAU = 0.05
N_TEST = 30
EPS_TARGET = 0.10  # phase-transition threshold


def qpsq_mae(
    process: UnitaryProcess,
    oracle: GaussianQPStat,
    obs,
    n_qubits: int,
    n_query: int,
    test_states,
    truths: np.ndarray,
    rng: np.random.Generator,
) -> float:
    k = k_from_epsilon(EPS)
    et = epsilon_tilde(EPS, n_qubits, k)
    samples = gather(oracle, n_qubits, obs, TAU, n_query, rng)
    model = learn(samples, n_qubits, k, et, pauli_l1_norm(obs))
    M = model.as_sparse_pauli_op().to_matrix()
    preds = np.array(
        [float(np.real(s.data.conj() @ M @ s.data)) for s in test_states]
    )
    return float(np.mean(np.abs(preds - truths)))


def shadow_mae(
    process: UnitaryProcess,
    obs,
    n_query: int,
    test_states,
    truths: np.ndarray,
    rng: np.random.Generator,
) -> float:
    """Each test state gets `n_query / T_test` randomized Pauli snapshots
    of `E(rho_test)`."""
    n_test = len(test_states)
    n_shadows_per = max(1, n_query // n_test)
    preds = []
    for psi in test_states:
        evolved = process.evolve_density(DensityMatrix(psi))
        preds.append(shadow_evolution_estimator(evolved, obs, n_shadows_per, rng))
    return float(np.mean(np.abs(np.array(preds) - truths)))


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)
    obs = pauli_z_first(N_QUBITS)

    rows: list[tuple] = []
    t_start = time.time()
    for t_count in T_COUNTS:
        for inst in range(N_INSTANCES):
            U = clifford_plus_t_circuit(N_QUBITS, DEPTH, t_count, rng)
            process = UnitaryProcess(U)
            test_states = [
                random_stabilizer_product_state(N_QUBITS, rng).to_statevector()
                for _ in range(N_TEST)
            ]
            truths = np.array([process.expectation(s, obs) for s in test_states])

            for n_query in QUERY_BUDGETS:
                qpsq_oracle = GaussianQPStat(process, rng=rng)
                mae_qpsq = qpsq_mae(
                    process, qpsq_oracle, obs, N_QUBITS, n_query,
                    test_states, truths, rng,
                )
                mae_shadow = shadow_mae(
                    process, obs, n_query, test_states, truths, rng,
                )
                rows.append((t_count, inst, n_query, "qpsq", mae_qpsq))
                rows.append((t_count, inst, n_query, "shadow", mae_shadow))
        print(f"  t_count={t_count} done at t={time.time() - t_start:.1f}s")

    csv_path = RESULTS_DIR / "ext_t_count_complexity.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_count", "instance", "n_query", "method", "mae"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], f"{r[4]:.6f}"])

    # Plot 1: MAE vs query budget for each t_count (QPSQ + shadow)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
    for ax, method in zip(axes, ("qpsq", "shadow"), strict=True):
        for t_count in T_COUNTS:
            xs, ys, ystds = [], [], []
            for n_query in QUERY_BUDGETS:
                vals = [r[4] for r in rows
                        if r[0] == t_count and r[2] == n_query and r[3] == method]
                xs.append(n_query)
                ys.append(float(np.mean(vals)))
                ystds.append(float(np.std(vals)))
            ax.errorbar(xs, ys, yerr=ystds, marker="o",
                        label=f"T-count={t_count}", capsize=3, linewidth=1)
        ax.axhline(EPS_TARGET, color="k", linestyle=":", linewidth=1,
                   label=f"target MAE = {EPS_TARGET}")
        ax.set_xscale("log")
        ax.set_xlabel("total query budget Q")
        ax.set_title(f"{method.upper()} (Algorithm 1)" if method == "qpsq"
                     else f"{method.upper()} (per-test classical shadows)")
        ax.grid(True, which="both", linestyle=":", alpha=0.5)
    axes[0].set_ylabel("mean absolute prediction error")
    axes[1].legend(loc="upper right", fontsize=7, ncol=2)
    fig.suptitle(
        f"Sample complexity vs T-count "
        f"(n={N_QUBITS}, depth={DEPTH}, {N_INSTANCES} instances each)"
    )
    fig.tight_layout()
    out = FIGS_DIR / "ext_t_count_complexity.pdf"
    fig.savefig(out)

    # Plot 2: minimum Q to reach target MAE vs T-count (the phase-transition view)
    fig2, ax2 = plt.subplots(figsize=(7, 4.5))
    for method, color, marker in (("qpsq", "C0", "o"), ("shadow", "C3", "s")):
        ys, ystds = [], []
        for t_count in T_COUNTS:
            min_qs: list[float] = []
            for inst in range(N_INSTANCES):
                inst_rows = [r for r in rows
                             if r[0] == t_count and r[1] == inst and r[3] == method]
                inst_rows.sort(key=lambda r: r[2])
                hit = next((r[2] for r in inst_rows if r[4] <= EPS_TARGET), None)
                if hit is not None:
                    min_qs.append(float(hit))
            if min_qs:
                ys.append(float(np.mean(min_qs)))
                ystds.append(float(np.std(min_qs)))
            else:
                ys.append(np.nan)
                ystds.append(0.0)
        ax2.errorbar(T_COUNTS, ys, yerr=ystds, marker=marker, color=color,
                     label=method.upper(), capsize=4, linewidth=1.5)
    ax2.set_xlabel("T-count")
    ax2.set_ylabel(f"minimum Q to reach MAE <= {EPS_TARGET}")
    ax2.set_yscale("log")
    ax2.set_title(
        f"Phase transition in sample complexity vs T-count (n={N_QUBITS})"
    )
    ax2.grid(True, which="both", linestyle=":", alpha=0.5)
    ax2.legend()
    fig2.tight_layout()
    out2 = FIGS_DIR / "ext_t_count_phase_transition.pdf"
    fig2.savefig(out2)
    print(f"wrote {csv_path}, {out}, and {out2}")


if __name__ == "__main__":
    main()
