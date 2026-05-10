"""Extension: how does Algorithm 1 degrade as the underlying process gets noisier?

The paper's analysis assumes the QPStat oracle returns within tolerance tau,
but says nothing about whether the *true* process being sampled is unitary.
Real devices implement noisy unitaries — depolarizing or amplitude damping
on each qubit after the gate. This experiment fixes a Haar-random 6-qubit
unitary, sweeps a per-qubit noise rate `p`, and measures Algorithm 1's
prediction MAE.

Two noise models are compared:

  * symmetric depolarizing: rho -> (1 - p) rho + (p / (2^n - 1)) * (I - rho)
    on each qubit (single-qubit form).
  * amplitude damping with damping parameter `gamma = p`.

Outputs:
  experiments/results/ext_noisy_channel.csv
  report/figs/ext_noisy_channel.pdf
"""

from __future__ import annotations

import csv
import time

import matplotlib.pyplot as plt
import numpy as np
from qiskit.quantum_info import Kraus

from qpsq.algorithm import gather, k_from_epsilon, learn, epsilon_tilde
from qpsq.designs import haar_unitary, random_stabilizer_product_state
from qpsq.observables import pauli_l1_norm, pauli_z_first
from qpsq.oracle import GaussianQPStat, _NoisyUnitary

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

N_QUBITS = 6
N_UNITARIES = 5
EPS = 0.3
TAU = 0.1
N_QUERIES = 4000
N_TEST = 60
NOISE_RATES = (0.0, 0.02, 0.05, 0.1, 0.2)


def depolarizing_kraus(p: float) -> Kraus:
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return Kraus(
        [
            np.sqrt(1 - 3 * p / 4) * I,
            np.sqrt(p / 4) * X,
            np.sqrt(p / 4) * Y,
            np.sqrt(p / 4) * Z,
        ]
    )


def amplitude_damping_kraus(gamma: float) -> Kraus:
    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype=complex)
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return Kraus([K0, K1])


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)
    obs = pauli_z_first(N_QUBITS)
    obs_l1 = pauli_l1_norm(obs)
    k = k_from_epsilon(EPS)
    eps_tilde = epsilon_tilde(EPS, N_QUBITS, k)

    rows: list[tuple] = []
    t_start = time.time()
    for u_idx in range(N_UNITARIES):
        U = haar_unitary(N_QUBITS, rng)
        # Pre-fix one stabilizer test set per unitary so results across noise
        # levels are directly comparable.
        test_states = [
            random_stabilizer_product_state(N_QUBITS, rng).to_statevector()
            for _ in range(N_TEST)
        ]

        for noise_kind, builder in (
            ("depolarizing", depolarizing_kraus),
            ("amplitude_damping", amplitude_damping_kraus),
        ):
            for p in NOISE_RATES:
                if p == 0.0:
                    # No-noise baseline: still go through the noisy wrapper for
                    # an apples-to-apples comparison (Kraus has float roundoff).
                    proc = _NoisyUnitary(unitary=U, per_qubit_kraus=builder(0.0))
                else:
                    proc = _NoisyUnitary(unitary=U, per_qubit_kraus=builder(p))

                truths = np.array([proc.expectation(s, obs) for s in test_states])
                oracle = GaussianQPStat(proc, rng=rng)
                samples = gather(oracle, N_QUBITS, obs, TAU, N_QUERIES, rng)
                model = learn(samples, N_QUBITS, k, eps_tilde, obs_l1)
                M = model.as_sparse_pauli_op().to_matrix()
                preds = np.array(
                    [float(np.real(s.data.conj() @ M @ s.data)) for s in test_states]
                )
                err = preds - truths
                mae = float(np.mean(np.abs(err)))
                rows.append((u_idx, noise_kind, p, mae, float(np.std(err))))
        print(f"  unitary {u_idx + 1}/{N_UNITARIES} done at "
              f"t={time.time() - t_start:.1f}s")

    csv_path = RESULTS_DIR / "ext_noisy_channel.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["unitary_idx", "noise_kind", "p", "mae", "err_std"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], f"{r[3]:.6f}", f"{r[4]:.6f}"])

    fig, ax = plt.subplots(figsize=(7, 4))
    for noise_kind, marker in (("depolarizing", "o"), ("amplitude_damping", "s")):
        xs, ys, yerr = [], [], []
        for p in NOISE_RATES:
            mae_list = [r[3] for r in rows if r[1] == noise_kind and r[2] == p]
            xs.append(p)
            ys.append(float(np.mean(mae_list)))
            yerr.append(float(np.std(mae_list)))
        ax.errorbar(xs, ys, yerr=yerr, marker=marker, label=noise_kind, capsize=3)
    ax.set_xlabel("per-qubit noise rate p")
    ax.set_ylabel("Algorithm 1 prediction MAE")
    ax.set_title(
        f"Algorithm 1 robustness vs. noise (n={N_QUBITS}, eps={EPS}, tau={TAU}, N={N_QUERIES})"
    )
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.legend()
    fig.tight_layout()
    out = FIGS_DIR / "ext_noisy_channel.pdf"
    fig.savefig(out)
    print(f"wrote {csv_path} and {out}")


if __name__ == "__main__":
    main()
