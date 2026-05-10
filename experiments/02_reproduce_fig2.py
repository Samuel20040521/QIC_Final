"""Reproduce paper Figure 2.

For each of `N_UNITARIES` Haar-random 6-qubit unitaries, run Algorithm 1 with
the QPStat-tolerance values `TAUS` and growing query budgets `NS`. Predict
`tr(Z_0 E(rho))` on 50 target states drawn from each of three distributions:

  * uniform over computational basis states
  * uniform over stabilizer product states
  * Haar-random pure states

Plot mean absolute prediction error vs. query budget `N`, one subplot per
target distribution, one curve per tau.

Outputs:
  experiments/results/fig2.csv      — long-form results
  report/figs/fig2_reproduction.pdf — three-panel plot
"""

from __future__ import annotations

import csv
import time

import matplotlib.pyplot as plt
import numpy as np
from qiskit.quantum_info import Statevector, random_statevector

from qpsq.algorithm import gather, k_from_epsilon, learn, epsilon_tilde
from qpsq.designs import haar_unitary, random_stabilizer_product_state
from qpsq.observables import pauli_l1_norm, pauli_z_first
from qpsq.oracle import GaussianQPStat, UnitaryProcess

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

N_QUBITS = 6
N_UNITARIES = 10
EPS = 0.3                 # target prediction error
TAUS = (0.05, 0.1, 0.2)   # QPStat tolerances
NS = (1000, 2000, 4000, 8000)  # query budgets
N_TEST = 50               # test states per distribution


def _make_test_states(
    distribution: str, n: int, n_states: int, rng: np.random.Generator
) -> list[Statevector]:
    if distribution == "computational":
        out: list[Statevector] = []
        for _ in range(n_states):
            idx = int(rng.integers(0, 2**n))
            vec = np.zeros(2**n, dtype=complex)
            vec[idx] = 1.0
            out.append(Statevector(vec))
        return out
    if distribution == "stabilizer":
        return [
            random_stabilizer_product_state(n, rng).to_statevector()
            for _ in range(n_states)
        ]
    if distribution == "haar":
        seeds = [int(rng.integers(0, 2**32 - 1)) for _ in range(n_states)]
        return [Statevector(random_statevector(2**n, seed=s)) for s in seeds]
    raise ValueError(distribution)


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)
    obs = pauli_z_first(N_QUBITS)
    obs_l1 = pauli_l1_norm(obs)
    k = k_from_epsilon(EPS)
    eps_tilde = epsilon_tilde(EPS, N_QUBITS, k)

    csv_path = RESULTS_DIR / "fig2.csv"
    fcsv = csv_path.open("w", newline="")
    writer = csv.writer(fcsv)
    writer.writerow(["unitary_idx", "tau", "N", "distribution", "mae", "rmse"])

    distributions = ("computational", "stabilizer", "haar")
    # Aggregate: results[distribution][tau][N] -> list of mean-abs-errors per unitary
    results: dict[str, dict[float, dict[int, list[float]]]] = {
        d: {t: {n: [] for n in NS} for t in TAUS} for d in distributions
    }

    t_start = time.time()
    for u_idx in range(N_UNITARIES):
        U = haar_unitary(N_QUBITS, rng)
        process = UnitaryProcess(U)

        # Pre-compute true expectation values for a fixed test set so all
        # (tau, N) cells share the same evaluation set per unitary.
        test_sets = {
            d: _make_test_states(d, N_QUBITS, N_TEST, rng) for d in distributions
        }
        truths = {
            d: np.array([process.expectation(s, obs) for s in test_sets[d]])
            for d in distributions
        }

        for tau in TAUS:
            oracle = GaussianQPStat(process, rng=rng)
            # Gather max(NS) once, then slice for smaller budgets.
            samples = gather(oracle, N_QUBITS, obs, tau, max(NS), rng)
            for n_query in NS:
                model = learn(samples[:n_query], N_QUBITS, k, eps_tilde, obs_l1)
                op_matrix = model.as_sparse_pauli_op().to_matrix()
                for d in distributions:
                    preds = np.array(
                        [
                            float(np.real(s.data.conj() @ op_matrix @ s.data))
                            for s in test_sets[d]
                        ]
                    )
                    err = preds - truths[d]
                    mae = float(np.mean(np.abs(err)))
                    rmse = float(np.sqrt(np.mean(err**2)))
                    results[d][tau][n_query].append(mae)
                    writer.writerow([u_idx, tau, n_query, d, f"{mae:.6f}", f"{rmse:.6f}"])
        print(f"  unitary {u_idx + 1}/{N_UNITARIES} done at "
              f"t={time.time() - t_start:.1f}s")
    fcsv.close()

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    distribution_titles = {
        "computational": "computational basis",
        "stabilizer": "stabilizer product",
        "haar": "Haar-random",
    }
    for ax, d in zip(axes, distributions, strict=True):
        for tau, marker in zip(TAUS, ("o", "s", "^"), strict=True):
            xs = list(NS)
            ys = [float(np.mean(results[d][tau][n_query])) for n_query in NS]
            yerr = [float(np.std(results[d][tau][n_query])) for n_query in NS]
            ax.errorbar(xs, ys, yerr=yerr, marker=marker, label=f"tau={tau}",
                        capsize=3, linewidth=1)
        ax.set_xscale("log")
        ax.set_xlabel("N (number of QPSQ queries)")
        ax.set_title(distribution_titles[d])
        ax.grid(True, which="both", linestyle=":", alpha=0.4)
    axes[0].set_ylabel("mean abs. prediction error")
    axes[-1].legend(loc="upper right", fontsize=8)
    fig.suptitle(
        f"Fig. 2 reproduction (n={N_QUBITS}, eps={EPS}, "
        f"{N_UNITARIES} Haar unitaries, observable Z_0)"
    )
    fig.tight_layout()
    out = FIGS_DIR / "fig2_reproduction.pdf"
    fig.savefig(out)
    print(f"wrote {csv_path} and {out} (total runtime {time.time() - t_start:.1f}s)")


if __name__ == "__main__":
    main()
