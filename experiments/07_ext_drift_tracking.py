"""Extension: drift-aware (online) Algorithm 1 for non-Markovian QPStat.

The QPStat oracle in Wadhwa-Doosti is implicitly memoryless: every query
is against the same channel `E`. Real devices drift over the timescale
of `N` queries (1/f noise, calibration wander). This experiment models a
linear drift in per-qubit depolarising rate from `p_start` to `p_end`
over the full query budget, then compares three estimators predicting
properties of the *end* channel `E_{t = N}`:

  * vanilla `learn` (uniform 1/N weighting --- the paper's Algorithm 1),
  * `online_learn` with several `half_life` settings (exponential
    weighting; smaller half-life = more recency bias),
  * a classical-shadow baseline that takes K = N / T_test snapshots
    against the current end-channel for each test state (no training).

Outputs:
  experiments/results/ext_drift_tracking.csv
  report/figs/ext_drift_tracking.pdf
"""

from __future__ import annotations

import csv
import time

import matplotlib.pyplot as plt
import numpy as np
from qiskit.quantum_info import DensityMatrix, Kraus

from qpsq.algorithm import (
    GatheredSample,
    epsilon_tilde,
    k_from_epsilon,
    learn,
    online_learn,
)
from qpsq.designs import (
    haar_unitary,
    random_stabilizer_product_state,
)
from qpsq.observables import pauli_l1_norm, pauli_z_first
from qpsq.oracle import (
    DriftingProcess,
    GaussianQPStat,
    NoisyUnitaryProcess,
    linear_depolarizing_drift,
)
from qpsq.shadows import shadow_evolution_estimator

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

N_QUBITS = 4
N_INSTANCES = 4
N_QUERIES = 4000
P_START = 0.0
P_END_VALUES = (0.05, 0.10, 0.20, 0.40)
HALF_LIVES = (200, 500, 1000, 2000)  # finite + we add "vanilla" as a label
EPS = 0.3
TAU = 0.05
N_TEST = 50


def end_process(base_unitary, p_end: float) -> NoisyUnitaryProcess:
    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)
    p = max(0.0, min(1.0, p_end))
    kraus = Kraus([
        np.sqrt(1 - 3 * p / 4) * I,
        np.sqrt(p / 4) * X,
        np.sqrt(p / 4) * Y,
        np.sqrt(p / 4) * Z,
    ])
    return NoisyUnitaryProcess(unitary=base_unitary, per_qubit_kraus=kraus)


def predict_mae(model, test_states, truths) -> float:
    M = model.as_sparse_pauli_op().to_matrix()
    preds = np.array(
        [float(np.real(s.data.conj() @ M @ s.data)) for s in test_states]
    )
    return float(np.mean(np.abs(preds - truths)))


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)
    obs = pauli_z_first(N_QUBITS)
    obs_l1 = pauli_l1_norm(obs)
    k = k_from_epsilon(EPS)
    et = epsilon_tilde(EPS, N_QUBITS, k)

    rows: list[tuple] = []
    t_start = time.time()
    for inst in range(N_INSTANCES):
        base = haar_unitary(N_QUBITS, rng)
        test_states = [
            random_stabilizer_product_state(N_QUBITS, rng).to_statevector()
            for _ in range(N_TEST)
        ]
        for p_end in P_END_VALUES:
            current = end_process(base, p_end)
            truths = np.array([current.expectation(s, obs) for s in test_states])

            # 1) Run drift-aware gather: drift rate sweeps p across queries.
            factory = linear_depolarizing_drift(
                base, p_start=P_START, p_end=p_end, n_total_queries=N_QUERIES
            )
            drifting = DriftingProcess(factory)
            oracle = GaussianQPStat(drifting, rng=rng)
            samples: list[GatheredSample] = []
            for _ in range(N_QUERIES):
                state = random_stabilizer_product_state(N_QUBITS, rng)
                y = oracle.query(state.to_statevector(), obs, TAU)
                samples.append(GatheredSample(state=state, y=y))

            # 2) Vanilla learn.
            model_flat = learn(samples, N_QUBITS, k, et, obs_l1)
            mae_flat = predict_mae(model_flat, test_states, truths)
            rows.append((inst, p_end, "vanilla", float("inf"), mae_flat))

            # 3) Online learn at several half_life values.
            for hl in HALF_LIVES:
                model_online = online_learn(
                    samples, N_QUBITS, k, et, obs_l1, half_life=float(hl)
                )
                mae_online = predict_mae(model_online, test_states, truths)
                rows.append((inst, p_end, "online", float(hl), mae_online))

            # 4) Shadow baseline against the current channel.
            n_shadows_per = max(1, N_QUERIES // N_TEST)
            shadow_preds = []
            for s in test_states:
                evolved = current.evolve_density(DensityMatrix(s))
                shadow_preds.append(
                    shadow_evolution_estimator(evolved, obs, n_shadows_per, rng)
                )
            mae_shadow = float(np.mean(np.abs(np.array(shadow_preds) - truths)))
            rows.append((inst, p_end, "shadow", 0.0, mae_shadow))

        print(f"  instance {inst + 1}/{N_INSTANCES} done at "
              f"t={time.time() - t_start:.1f}s")

    csv_path = RESULTS_DIR / "ext_drift_tracking.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["instance", "p_end", "method", "half_life", "mae"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], f"{r[4]:.6f}"])

    fig, ax = plt.subplots(figsize=(8, 4.5))
    # x-axis: p_end (drift magnitude). One curve per method (vanilla, shadow,
    # and one per half_life).
    for label_method, marker, color in (
        ("vanilla", "o", "C0"),
        ("shadow", "s", "C3"),
    ):
        ys, ystd = [], []
        for p_end in P_END_VALUES:
            vals = [r[4] for r in rows if r[1] == p_end and r[2] == label_method]
            ys.append(float(np.mean(vals)))
            ystd.append(float(np.std(vals)))
        ax.errorbar(P_END_VALUES, ys, yerr=ystd, marker=marker, color=color,
                    label=("Algorithm 1 (uniform)" if label_method == "vanilla"
                           else "shadow baseline"),
                    capsize=4, linewidth=1.6)
    cmap = plt.get_cmap("viridis")
    for i, hl in enumerate(HALF_LIVES):
        ys, ystd = [], []
        for p_end in P_END_VALUES:
            vals = [r[4] for r in rows
                    if r[1] == p_end and r[2] == "online" and r[3] == float(hl)]
            ys.append(float(np.mean(vals)))
            ystd.append(float(np.std(vals)))
        ax.errorbar(P_END_VALUES, ys, yerr=ystd, marker="^",
                    color=cmap(0.15 + 0.7 * i / max(1, len(HALF_LIVES) - 1)),
                    label=f"online (half_life={hl})", capsize=3, linewidth=1)
    ax.set_xlabel("drift end-rate p_end (per-qubit depolarising)")
    ax.set_ylabel("MAE on end-channel test set")
    ax.set_title(
        f"Drift-aware Algorithm 1 (n={N_QUBITS}, N={N_QUERIES}, "
        f"{N_INSTANCES} instances)"
    )
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    out = FIGS_DIR / "ext_drift_tracking.pdf"
    fig.savefig(out)
    print(f"wrote {csv_path} and {out}")


if __name__ == "__main__":
    main()
