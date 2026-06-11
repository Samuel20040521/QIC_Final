"""Reproduce paper Figure 2.

For each of `N_UNITARIES` Haar-random 6-qubit unitaries, run Algorithm 1 with
the QPStat-tolerance values `TAUS` and growing query budgets `NS`. Predict
`tr(Z_0 E(rho))` on target states drawn from three distributions:

  * all 2^n computational basis states
  * uniform over stabilizer product states (100 states)
  * Haar-random pure states (100 states)

Plot mean SQUARED prediction error vs. query budget `N`, one panel per
target distribution, one curve per tau — mirroring the paper's Figure 2.

Hyperparameter note (reproducibility gap we document in the report):
the paper's Eq. (42) states `k = ceil(log_1.5(2 / eps^2))`, but the
authors' published simulation code (qpsq-learning, `coeff.py:hyperparams`)
uses `k = ceil(log_1.5(2 / eps))` with `eps = 0.9`, giving k = 2. Their
code also computes an effective threshold `eps_tilde ~ 1.5e-11`, i.e. the
two-part threshold of Algorithm 1 is inactive in the published figure.
We match the authors' code here, because that is what produced the paper's
Figure 2: k = 2, thresholds off (EPS_TILDE ~ 0), error metric = MSE.

Choosing instead eps = 0.3 under the Eq. (42) formula gives k = 8 >= n,
i.e. *no truncation*: every Pauli coefficient up to weight n is estimated
with a 3^|P| variance amplification, the prediction error explodes, and
the Haar-state panel inverts from best to worst. See report Section 4.

Outputs:
  experiments/results/fig2.csv      — long-form results (mse + mae)
  report/figs/fig2_reproduction.pdf — three-panel plot
"""

from __future__ import annotations

import csv
import math
import time

import matplotlib.pyplot as plt
import numpy as np
from qiskit.quantum_info import Statevector, random_statevector

from qpsq.algorithm import gather, learn
from qpsq.designs import haar_unitary, random_stabilizer_product_state
from qpsq.observables import pauli_l1_norm, pauli_z_first
from qpsq.oracle import GaussianQPStat, UnitaryProcess

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

N_QUBITS = 6
N_UNITARIES = 10
EPS = 0.9                      # authors' choice in sim_haar.py
# Authors' code formula (NOT Eq. 42 — see module docstring): k = 2 for eps=0.9.
K_TRUNC = math.ceil(math.log(2.0 / EPS) / math.log(1.5))
EPS_TILDE = 1e-12              # thresholds inactive, as in the authors' code
TAUS = (0.05, 0.15, 0.25)      # tolerances plotted in the paper's Figure 2
NS = tuple(range(50, 501, 50))  # query budgets, as in the paper
N_TEST = 100                   # test states for stabilizer / haar distributions


def _make_test_states(
    distribution: str, n: int, rng: np.random.Generator
) -> list[Statevector]:
    if distribution == "computational":
        # The paper evaluates on ALL 2^n computational basis states.
        out: list[Statevector] = []
        for idx in range(2**n):
            vec = np.zeros(2**n, dtype=complex)
            vec[idx] = 1.0
            out.append(Statevector(vec))
        return out
    if distribution == "stabilizer":
        return [
            random_stabilizer_product_state(n, rng).to_statevector()
            for _ in range(N_TEST)
        ]
    if distribution == "haar":
        seeds = [int(rng.integers(0, 2**32 - 1)) for _ in range(N_TEST)]
        return [Statevector(random_statevector(2**n, seed=s)) for s in seeds]
    raise ValueError(distribution)


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)
    obs = pauli_z_first(N_QUBITS)
    obs_l1 = pauli_l1_norm(obs)
    print(f"eps={EPS} -> k={K_TRUNC} (authors' formula), eps_tilde={EPS_TILDE}")

    csv_path = RESULTS_DIR / "fig2.csv"
    fcsv = csv_path.open("w", newline="")
    writer = csv.writer(fcsv)
    writer.writerow(["unitary_idx", "tau", "N", "distribution", "mse", "mae"])

    distributions = ("computational", "stabilizer", "haar")
    # As in the authors' code: one fixed set of unitaries and one fixed set of
    # test states, shared across all (tau, N) cells.
    unitaries = [haar_unitary(N_QUBITS, rng) for _ in range(N_UNITARIES)]
    test_sets = {d: _make_test_states(d, N_QUBITS, rng) for d in distributions}

    # results[distribution][tau][N] -> list of per-unitary MSEs
    results: dict[str, dict[float, dict[int, list[float]]]] = {
        d: {t: {n: [] for n in NS} for t in TAUS} for d in distributions
    }

    t_start = time.time()
    for u_idx, U in enumerate(unitaries):
        process = UnitaryProcess(U)
        truths = {
            d: np.array([process.expectation(s, obs) for s in test_sets[d]])
            for d in distributions
        }
        for tau in TAUS:
            for n_query in NS:
                # Fresh dataset per (unitary, tau, N) — matches the authors'
                # code, which calls learn() with a new dataset for each cell.
                oracle = GaussianQPStat(process, rng=rng)
                samples = gather(oracle, N_QUBITS, obs, tau, n_query, rng)
                model = learn(samples, N_QUBITS, K_TRUNC, EPS_TILDE, obs_l1)
                op_matrix = model.as_sparse_pauli_op().to_matrix()
                for d in distributions:
                    preds = np.array(
                        [
                            float(np.real(s.data.conj() @ op_matrix @ s.data))
                            for s in test_sets[d]
                        ]
                    )
                    err = preds - truths[d]
                    mse = float(np.mean(err**2))
                    mae = float(np.mean(np.abs(err)))
                    results[d][tau][n_query].append(mse)
                    writer.writerow(
                        [u_idx, tau, n_query, d, f"{mse:.6f}", f"{mae:.6f}"]
                    )
        print(f"  unitary {u_idx + 1}/{N_UNITARIES} done at "
              f"t={time.time() - t_start:.1f}s")
    fcsv.close()

    # Plot: paper layout — comp / stab on the top row, haar centered below.
    # Independent y-axes per panel (the haar panel's range is ~10x smaller).
    distribution_titles = {
        "computational": "comp",
        "stabilizer": "stab",
        "haar": "haar",
    }
    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(2, 4)
    axes = {
        "computational": fig.add_subplot(gs[0, 0:2]),
        "stabilizer": fig.add_subplot(gs[0, 2:4]),
        "haar": fig.add_subplot(gs[1, 1:3]),
    }
    for d, ax in axes.items():
        for tau in TAUS:
            ys = [float(np.mean(results[d][tau][n_query])) for n_query in NS]
            ax.plot(NS, ys, linestyle="--", label=f"Tolerance: {tau}")
        ax.set_xlabel("Number of queries")
        ax.set_ylabel("Error")
        ax.set_title(
            f"Learning error for haar-random unitaries, "
            f"Distr = {distribution_titles[d]}, n = {N_QUBITS}",
            fontsize=10,
        )
        ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    out = FIGS_DIR / "fig2_reproduction.pdf"
    fig.savefig(out)
    print(f"wrote {csv_path} and {out} (total runtime {time.time() - t_start:.1f}s)")


if __name__ == "__main__":
    main()
