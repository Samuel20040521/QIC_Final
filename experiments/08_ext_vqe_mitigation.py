"""Extension: surrogate VQE mitigation using a one-shot QPSQ model.

We pick a 4-qubit transverse-field Ising Hamiltonian
    H = -J sum_i Z_i Z_{i+1} - h sum_i X_i
and a hardware-efficient ansatz `U(theta)`. The "noisy device" applies an
ideal `U(theta)` followed by an effective decoherence channel `Lambda`
(per-qubit depolarising). This factorisation is the standard
"noise-after-circuit" approximation for fixed-depth ansaetze on a
calibrated device. We additionally check the assumption against a more
realistic Aer `NoiseModel` (hardware realism point in the review).

We compare three estimators of `E_noisy(theta) = <H>_{Lambda * U(theta)}`:

  * Noiseless reference: <H>_{U(theta)} (the VQE optimum target).
  * Raw noisy: tr(H Lambda(U(theta)|0><0|U(theta)^dagger)). Computed
    every iteration; would cost one full circuit run per VQE step.
  * Surrogate QPSQ mitigation: train Algorithm 1 ONCE on `Lambda`
    (as a stand-alone CPTP map) using N samples; for any `theta`,
    evaluate <psi_theta|Lambda^dagger(H)|psi_theta> in CLOSED FORM ---
    no further oracle queries.
  * Classical-shadow baseline: K = N / (number of evaluation points)
    randomized Pauli snapshots per evaluation against the noisy circuit.

The factorisation: the QPSQ oracle's process is `Lambda` alone. The
learned operator `M_Lambda` ~= Lambda^dagger(H) does not contain any
ansatz unitary, so it transfers across `theta` correctly. This avoids
the bug where a unitary baked into the surrogate at theta_* would
mis-evaluate at theta != theta_*.

Outputs:
  experiments/results/ext_vqe_mitigation.csv
  report/figs/ext_vqe_mitigation.pdf
"""

from __future__ import annotations

import csv
import time

import matplotlib.pyplot as plt
import numpy as np
from qiskit.circuit import ParameterVector, QuantumCircuit
from qiskit.quantum_info import (
    DensityMatrix,
    Kraus,
    Operator,
    SparsePauliOp,
    Statevector,
)

from qpsq.algorithm import (
    epsilon_tilde,
    gather,
    k_from_epsilon,
    learn,
)
from qpsq.observables import pauli_l1_norm
from qpsq.oracle import (
    CircuitNoiseProcess,
    GaussianQPStat,
    NoisyUnitaryProcess,
    default_hardware_noise_model,
)
from qpsq.shadows import shadow_evolution_estimator

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

N_QUBITS = 4
J = 1.0
H_FIELD = 0.5
ANSATZ_LAYERS = 2
N_QPSQ_SAMPLES = 4000
TAU = 0.05
EPS = 0.3
N_THETA_NEIGHBORS = 32
THETA_RADIUS = 1.0
SHADOW_BUDGET_PER_EVAL = 200
P_DEPOL = 0.04  # effective per-qubit depolarising rate of `Lambda`


def tfim_hamiltonian(n: int = N_QUBITS, j: float = J, h: float = H_FIELD) -> SparsePauliOp:
    terms: list[tuple[str, float]] = []
    for i in range(n - 1):
        label = ["I"] * n
        label[n - 1 - i] = "Z"
        label[n - 1 - (i + 1)] = "Z"
        terms.append(("".join(label), -j))
    for i in range(n):
        label = ["I"] * n
        label[n - 1 - i] = "X"
        terms.append(("".join(label), -h))
    return SparsePauliOp.from_list(terms)


def make_ansatz(n: int = N_QUBITS, layers: int = ANSATZ_LAYERS):
    n_params = layers * 2 * n
    theta = ParameterVector("theta", n_params)
    qc = QuantumCircuit(n)
    p_idx = 0
    for _ in range(layers):
        for q in range(n):
            qc.ry(theta[p_idx], q); p_idx += 1
        for q in range(n):
            qc.rz(theta[p_idx], q); p_idx += 1
        for q in range(n - 1):
            qc.cx(q, q + 1)
    return qc, theta


def bind_ansatz(qc, theta, values: np.ndarray):
    return qc.assign_parameters(dict(zip(theta, values, strict=True)))


def depolarizing_kraus(p: float) -> Kraus:
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    return Kraus([
        np.sqrt(1 - 3 * p / 4) * I,
        np.sqrt(p / 4) * X,
        np.sqrt(p / 4) * Y,
        np.sqrt(p / 4) * Z,
    ])


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)
    H = tfim_hamiltonian()
    obs_l1 = pauli_l1_norm(H)
    qc_template, theta_param = make_ansatz()
    print(f"  H l1-norm = {obs_l1:.2f}, ansatz parameters = {qc_template.num_parameters}")

    # --- Define the noise channel `Lambda` (factorized model) -----------
    lam_kraus = depolarizing_kraus(P_DEPOL)
    lambda_proc = NoisyUnitaryProcess(
        unitary=Operator(np.eye(2**N_QUBITS)),
        per_qubit_kraus=lam_kraus,
    )

    # --- Step 1: coarse VQE on the noisy device -------------------------
    eval_count = {"raw_calls": 0}

    def noisy_energy_factored(theta_vals: np.ndarray) -> float:
        eval_count["raw_calls"] += 1
        bound = bind_ansatz(qc_template, theta_param, theta_vals)
        proc = NoisyUnitaryProcess(unitary=Operator(bound), per_qubit_kraus=lam_kraus)
        return proc.expectation(Statevector.from_label("0" * N_QUBITS), H)

    from scipy.optimize import minimize

    rng_init = np.random.default_rng(args.seed + 1)
    theta0 = rng_init.uniform(-0.5, 0.5, size=qc_template.num_parameters)
    print("  coarse VQE on factored noisy device...")
    t_start = time.time()
    res = minimize(noisy_energy_factored, theta0, method="COBYLA",
                   options={"maxiter": 80, "rhobeg": 0.4, "disp": False})
    theta_star = np.asarray(res.x)
    print(f"  coarse VQE done in {time.time() - t_start:.1f}s, "
          f"{eval_count['raw_calls']} noisy evals; E_noisy(theta*) = {res.fun:.4f}")

    # --- Step 2: train QPSQ on `Lambda` ALONE (factored from U) --------
    oracle_lambda = GaussianQPStat(lambda_proc, rng=rng)
    k = k_from_epsilon(EPS)
    et = epsilon_tilde(EPS, N_QUBITS, k)
    print(f"  training QPSQ on Lambda (N={N_QPSQ_SAMPLES})...")
    t0 = time.time()
    samples = gather(oracle_lambda, N_QUBITS, H, TAU, N_QPSQ_SAMPLES, rng)
    model = learn(samples, N_QUBITS, k, et, obs_l1)
    surrogate_M = model.as_sparse_pauli_op().to_matrix()
    print(f"  surrogate trained in {time.time() - t0:.1f}s")

    # Sanity check: at any rho, h(rho) ~= tr(H Lambda(rho)).
    sanity_state = Statevector.from_label("0" * N_QUBITS)
    sanity_truth = lambda_proc.expectation(sanity_state, H)
    sanity_pred = float(np.real(sanity_state.data.conj() @ surrogate_M @ sanity_state.data))
    print(f"  sanity at |0>: truth={sanity_truth:.4f}, surrogate={sanity_pred:.4f}, "
          f"diff={abs(sanity_truth - sanity_pred):.4f}")

    # --- Step 3: sweep theta around theta_*; evaluate energy 4 ways ----
    rng_n = np.random.default_rng(args.seed + 2)
    deltas = rng_n.uniform(-1, 1, size=(N_THETA_NEIGHBORS, qc_template.num_parameters))
    norms = np.linalg.norm(deltas, axis=1, keepdims=True)
    deltas /= np.maximum(norms, 1e-12)
    deltas *= THETA_RADIUS * rng_n.uniform(0.0, 1.0, size=(N_THETA_NEIGHBORS, 1))
    theta_grid = theta_star[None, :] + deltas

    # Hardware-realistic noise check: an Aer NoiseModel with per-gate errors,
    # used only at evaluation time (not in QPSQ training).
    nm_aer = default_hardware_noise_model(p1=0.001, p2=0.01, readout_p=0.0)

    rows: list[tuple] = []
    print(f"  evaluating {N_THETA_NEIGHBORS} neighbors...")
    t0 = time.time()
    for j, theta_vals in enumerate(theta_grid):
        bound = bind_ansatz(qc_template, theta_param, theta_vals)
        U_theta = Operator(bound)
        psi_theta = Statevector.from_label("0" * N_QUBITS).evolve(U_theta)

        # Noiseless reference.
        e_noiseless = float(np.real(psi_theta.expectation_value(H)))

        # Factored noisy: Lambda after U(theta).
        proc_noisy = NoisyUnitaryProcess(unitary=U_theta, per_qubit_kraus=lam_kraus)
        e_noisy_factored = proc_noisy.expectation(
            Statevector.from_label("0" * N_QUBITS), H
        )

        # Surrogate: <psi_theta | Lambda^dagger(H) | psi_theta>, closed form.
        e_surrogate = float(np.real(psi_theta.data.conj() @ surrogate_M @ psi_theta.data))

        # Hardware-realistic Aer NoiseModel for cross-check.
        proc_aer = CircuitNoiseProcess(circuit=bound, noise_model=nm_aer)
        e_noisy_aer = proc_aer.expectation(Statevector.from_label("0" * N_QUBITS), H)

        # Shadow baseline: snapshots against the factored noisy circuit.
        rho = proc_noisy.evolve_density(DensityMatrix(Statevector.from_label("0" * N_QUBITS)))
        e_shadow = shadow_evolution_estimator(rho, H, SHADOW_BUDGET_PER_EVAL, rng)

        rows.append((j, "noiseless", e_noiseless))
        rows.append((j, "noisy_factored", e_noisy_factored))
        rows.append((j, "noisy_aer", e_noisy_aer))
        rows.append((j, "surrogate", e_surrogate))
        rows.append((j, "shadow", e_shadow))
    print(f"  neighbor sweep done in {time.time() - t0:.1f}s")

    csv_path = RESULTS_DIR / "ext_vqe_mitigation.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["theta_idx", "method", "energy"])
        for r in rows:
            w.writerow([r[0], r[1], f"{r[2]:.6f}"])

    by: dict[str, np.ndarray] = {
        m: np.array([r[2] for r in rows if r[1] == m])
        for m in ("noiseless", "noisy_factored", "noisy_aer", "surrogate", "shadow")
    }
    print("\n  Energy MAE relative to factored noisy reference:")
    for m in ("surrogate", "shadow", "noisy_aer"):
        mae = float(np.mean(np.abs(by[m] - by["noisy_factored"])))
        print(f"    {m:>14}  MAE = {mae:.4f}")
    print("\n  Energy MAE relative to noiseless reference:")
    for m in ("noisy_factored", "noisy_aer", "surrogate", "shadow"):
        mae = float(np.mean(np.abs(by[m] - by["noiseless"])))
        print(f"    {m:>14}  MAE = {mae:.4f}")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    truth = by["noisy_factored"]
    lo, hi = min(truth.min(), by["noiseless"].min()) - 0.4, max(truth.max(), by["noiseless"].max()) + 0.4
    for ax, method, color in zip(
        axes, ("surrogate", "shadow", "noisy_aer"),
        ("C0", "C2", "C3"), strict=True,
    ):
        preds = by[method]
        ax.scatter(truth, preds, color=color, s=22, alpha=0.85, label=method)
        ax.plot([lo, hi], [lo, hi], "k--", linewidth=0.7)
        mae = float(np.mean(np.abs(preds - truth)))
        ax.set_title(f"{method} vs noisy_factored\n(MAE = {mae:.3f})")
        ax.set_xlabel("noisy_factored energy (target)")
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.grid(True, linestyle=":", alpha=0.4)
    axes[0].set_ylabel("estimated energy")
    fig.suptitle(
        f"Surrogate QPSQ mitigation for 4-qubit TFIM (N_train={N_QPSQ_SAMPLES}, "
        f"{N_THETA_NEIGHBORS} parameter neighbors of theta*)"
    )
    fig.tight_layout()
    out = FIGS_DIR / "ext_vqe_mitigation.pdf"
    fig.savefig(out)
    print(f"\n  wrote {csv_path} and {out}")


if __name__ == "__main__":
    main()
