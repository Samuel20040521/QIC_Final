"""Extension: how does Algorithm 1 perform on different circuit ensembles?

The paper's lower bound (Theorem 2) proves an exponential separation for
unitary 2-designs *with respect to diamond distance*. Algorithm 1 solves an
average-case prediction problem (a much weaker task) and is meant to work on
any quantum process. This experiment checks empirically: holding query
budget and tau fixed, does the algorithm's MAE depend on the structure of
the underlying unitary?

We compare:
  * Haar-random unitaries (the paper's setting),
  * uniformly random Cliffords (a 3-design),
  * brick-wall random circuits at depths 1, 2, 4, 8 (approaching a 2-design
    as depth grows).

Outputs:
  experiments/results/ext_circuit_classes.csv
  report/figs/ext_circuit_classes.pdf
"""

from __future__ import annotations

import csv
import time

import matplotlib.pyplot as plt
import numpy as np

from qpsq.algorithm import gather, k_from_epsilon, learn, epsilon_tilde
from qpsq.designs import (
    brickwall_random_circuit,
    haar_unitary,
    random_clifford_unitary,
    random_stabilizer_product_state,
)
from qpsq.observables import pauli_l1_norm, pauli_z_first
from qpsq.oracle import GaussianQPStat, UnitaryProcess

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

N_QUBITS = 6
N_INSTANCES = 8
EPS = 0.3
TAU = 0.05
N_QUERIES = 4000
N_TEST = 60

CLASSES: tuple[tuple[str, str], ...] = (
    ("haar", "Haar"),
    ("clifford", "random Clifford"),
    ("brickwall_d1", "brick-wall depth 1"),
    ("brickwall_d2", "brick-wall depth 2"),
    ("brickwall_d4", "brick-wall depth 4"),
    ("brickwall_d8", "brick-wall depth 8"),
)


def sample_unitary(class_id: str, n: int, rng: np.random.Generator):
    if class_id == "haar":
        return haar_unitary(n, rng)
    if class_id == "clifford":
        return random_clifford_unitary(n, rng)
    if class_id.startswith("brickwall_d"):
        depth = int(class_id.split("d")[-1])
        return brickwall_random_circuit(n, depth, rng)
    raise ValueError(class_id)


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)
    obs = pauli_z_first(N_QUBITS)
    obs_l1 = pauli_l1_norm(obs)
    k = k_from_epsilon(EPS)
    eps_tilde = epsilon_tilde(EPS, N_QUBITS, k)

    rows: list[tuple] = []
    t_start = time.time()
    for class_id, _label in CLASSES:
        for inst_idx in range(N_INSTANCES):
            U = sample_unitary(class_id, N_QUBITS, rng)
            process = UnitaryProcess(U)
            test_states = [
                random_stabilizer_product_state(N_QUBITS, rng).to_statevector()
                for _ in range(N_TEST)
            ]
            truths = np.array([process.expectation(s, obs) for s in test_states])
            oracle = GaussianQPStat(process, rng=rng)
            samples = gather(oracle, N_QUBITS, obs, TAU, N_QUERIES, rng)
            model = learn(samples, N_QUBITS, k, eps_tilde, obs_l1)
            M = model.as_sparse_pauli_op().to_matrix()
            preds = np.array(
                [float(np.real(s.data.conj() @ M @ s.data)) for s in test_states]
            )
            mae = float(np.mean(np.abs(preds - truths)))
            rows.append((class_id, inst_idx, mae, float(np.std(preds - truths))))
        print(f"  class {class_id} done at t={time.time() - t_start:.1f}s")

    csv_path = RESULTS_DIR / "ext_circuit_classes.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["class", "instance", "mae", "err_std"])
        for r in rows:
            w.writerow([r[0], r[1], f"{r[2]:.6f}", f"{r[3]:.6f}"])

    fig, ax = plt.subplots(figsize=(8, 4.5))
    labels: list[str] = []
    means: list[float] = []
    stds: list[float] = []
    for class_id, label in CLASSES:
        vals = [r[2] for r in rows if r[0] == class_id]
        labels.append(label)
        means.append(float(np.mean(vals)))
        stds.append(float(np.std(vals)))
    xs = np.arange(len(CLASSES))
    ax.bar(xs, means, yerr=stds, capsize=4, color="C0", alpha=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Algorithm 1 prediction MAE")
    ax.set_title(
        f"Prediction MAE across unitary ensembles "
        f"(n={N_QUBITS}, tau={TAU}, N={N_QUERIES}, {N_INSTANCES} instances each)"
    )
    ax.grid(True, axis="y", linestyle=":", alpha=0.5)
    fig.tight_layout()
    out = FIGS_DIR / "ext_circuit_classes.pdf"
    fig.savefig(out)
    print(f"wrote {csv_path} and {out}")


if __name__ == "__main__":
    main()
