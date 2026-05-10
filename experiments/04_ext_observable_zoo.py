"""Extension: empirical query complexity vs. observable weight.

The paper bounds the number of QPSQ queries as
`N = tau^2 M log(M n / delta) * 2^O(log n log 1/eps)` for `M` observables
each given as a sum of few-body terms (Corollary 1). Each `O_i` is assumed
to have low Pauli weight. This experiment replaces the paper's `O = Z_0`
with random Pauli strings of varying weight `w in {1, ..., n}` and measures
prediction MAE under a fixed query budget.

Outputs:
  experiments/results/ext_observable_zoo.csv
  report/figs/ext_observable_zoo.pdf
"""

from __future__ import annotations

import csv
import time

import matplotlib.pyplot as plt
import numpy as np
from qiskit.quantum_info import SparsePauliOp

from qpsq.algorithm import gather, k_from_epsilon, learn, epsilon_tilde
from qpsq.designs import haar_unitary, random_stabilizer_product_state
from qpsq.observables import pauli_l1_norm
from qpsq.oracle import GaussianQPStat, UnitaryProcess

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

N_QUBITS = 5
N_UNITARIES = 5
N_OBSERVABLES_PER_WEIGHT = 4
EPS = 0.3
TAU = 0.05
NS = (500, 1000, 2000, 4000)
N_TEST = 50


def random_pauli_with_weight(n: int, weight: int, rng: np.random.Generator) -> SparsePauliOp:
    """Build a Pauli string with exactly `weight` non-identity sites in Qiskit form."""
    positions = rng.choice(n, size=weight, replace=False)
    chars = ["I"] * n
    for pos in positions:
        chars[int(pos)] = ("X", "Y", "Z")[int(rng.integers(0, 3))]
    label = "".join(reversed(chars))
    return SparsePauliOp.from_list([(label, 1.0)])


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)
    k = k_from_epsilon(EPS)
    eps_tilde = epsilon_tilde(EPS, N_QUBITS, k)

    rows: list[tuple] = []
    t_start = time.time()
    for u_idx in range(N_UNITARIES):
        U = haar_unitary(N_QUBITS, rng)
        process = UnitaryProcess(U)
        test_states = [
            random_stabilizer_product_state(N_QUBITS, rng).to_statevector()
            for _ in range(N_TEST)
        ]
        for weight in range(1, N_QUBITS + 1):
            for obs_idx in range(N_OBSERVABLES_PER_WEIGHT):
                obs = random_pauli_with_weight(N_QUBITS, weight, rng)
                obs_l1 = pauli_l1_norm(obs)
                truths = np.array([process.expectation(s, obs) for s in test_states])
                oracle = GaussianQPStat(process, rng=rng)
                samples = gather(oracle, N_QUBITS, obs, TAU, max(NS), rng)
                for n_query in NS:
                    model = learn(samples[:n_query], N_QUBITS, k, eps_tilde, obs_l1)
                    M = model.as_sparse_pauli_op().to_matrix()
                    preds = np.array(
                        [float(np.real(s.data.conj() @ M @ s.data)) for s in test_states]
                    )
                    mae = float(np.mean(np.abs(preds - truths)))
                    rows.append((u_idx, weight, obs_idx, n_query, mae))
        print(f"  unitary {u_idx + 1}/{N_UNITARIES} done at "
              f"t={time.time() - t_start:.1f}s")

    csv_path = RESULTS_DIR / "ext_observable_zoo.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["unitary_idx", "weight", "obs_idx", "N", "mae"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], f"{r[4]:.6f}"])

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for weight in range(1, N_QUBITS + 1):
        xs, ys, ystd = [], [], []
        for n_query in NS:
            vals = [r[4] for r in rows if r[1] == weight and r[3] == n_query]
            xs.append(n_query)
            ys.append(float(np.mean(vals)))
            ystd.append(float(np.std(vals)))
        ax.errorbar(xs, ys, yerr=ystd, marker="o", label=f"weight {weight}", capsize=3)
    ax.set_xscale("log")
    ax.set_xlabel("N (queries)")
    ax.set_ylabel("Algorithm 1 prediction MAE")
    ax.set_title(
        f"MAE vs. query budget for varying Pauli weight (n={N_QUBITS}, tau={TAU})"
    )
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    out = FIGS_DIR / "ext_observable_zoo.pdf"
    fig.savefig(out)
    print(f"wrote {csv_path} and {out}")


if __name__ == "__main__":
    main()
