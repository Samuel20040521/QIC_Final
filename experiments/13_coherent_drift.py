"""Sprint Task 2: Coherent unitary drift, vanilla vs online learn.

For a sequence of drift rates `eta`, run a single 4000-step QPStat
trajectory against a `CoherentDriftProcess`. Apply *both* `learn`
(uniform weights) and `online_learn` at several half-life settings to
the same trajectory, then evaluate prediction MAE on the *current*
unitary `U(t = N - 1)`.

Resumable per (eta, instance) trajectory: each instance saves a JSON
sidecar with its trajectory ID; on restart we skip done IDs.

Goal of the experiment: at moderate drift rates, online with a
well-chosen half-life should achieve >= 30% lower tracking MAE than
vanilla.
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
    CoherentDriftProcess,
    GaussianQPStat,
    UnitaryProcess,
    random_pauli_drift_hamiltonian,
)

from _common import FIGS_DIR, RESULTS_DIR
from _hpc import default_max_workers, wait_for_memory_under_cap, wait_for_ram
from _resumable import ResumableCSV

N_QUBITS = 4
N_QUERIES = 4000
EPS = 0.3
TAU = 0.05
N_TEST = 60
DRIFT_RATES = (0.0, 0.0025, 0.005, 0.01, 0.02, 0.04)
HALF_LIVES = (250, 500, 1000, 2000)
DEFAULT_INSTANCES = 4
DRIFT_HAMILTONIAN_TERMS = 12
MEM_CAP_BYTES = 70 * 1024**3

KEY_COLS = ["eta", "instance", "method", "half_life"]
ALL_COLS = KEY_COLS + ["mae", "wall_seconds"]


@dataclass
class CellRow:
    eta: float
    instance: int
    method: str  # "vanilla" | "online"
    half_life: float  # 0.0 for vanilla
    mae: float
    wall_seconds: float


def _trajectory_results(
    eta: float,
    instance: int,
    seed: int,
) -> list[CellRow]:
    """One trajectory: gather 4000 samples through CoherentDriftProcess,
    then apply vanilla learn + online_learn at every half_life. Return one
    row per (method, half_life)."""
    rng = np.random.default_rng(seed)
    U0 = haar_unitary(N_QUBITS, rng)
    H_drift = random_pauli_drift_hamiltonian(
        N_QUBITS, n_terms=DRIFT_HAMILTONIAN_TERMS, rng=rng
    )
    drift_proc = CoherentDriftProcess(U0, H_drift, drift_rate=eta)
    obs = pauli_z_first(N_QUBITS)
    obs_l1 = pauli_l1_norm(obs)
    k = k_from_epsilon(EPS)
    et = epsilon_tilde(EPS, N_QUBITS, k)
    oracle = GaussianQPStat(drift_proc, rng=rng)

    # 1) Gather a single 4000-step trajectory.
    samples: list[GatheredSample] = []
    t_gather = time.time()
    for _ in range(N_QUERIES):
        s = random_stabilizer_product_state(N_QUBITS, rng)
        y = oracle.query(s.to_statevector(), obs, TAU)
        samples.append(GatheredSample(state=s, y=y))
    gather_time = time.time() - t_gather

    # 2) End-channel ground truth: the unitary at t = N_QUERIES - 1.
    end_unitary = drift_proc._U_at(N_QUERIES - 1)
    from qiskit.quantum_info import Operator
    end_proc = UnitaryProcess(Operator(end_unitary))
    test_states = [
        random_stabilizer_product_state(N_QUBITS, rng).to_statevector()
        for _ in range(N_TEST)
    ]
    truths = np.array([end_proc.expectation(s, obs) for s in test_states])

    rows: list[CellRow] = []

    # Vanilla learn
    t_l = time.time()
    model = learn(samples, N_QUBITS, k, et, obs_l1)
    op = model.as_sparse_pauli_op()
    preds = np.array([float(np.real(s.expectation_value(op))) for s in test_states])
    mae = float(np.mean(np.abs(preds - truths)))
    rows.append(CellRow(eta=eta, instance=instance, method="vanilla",
                        half_life=0.0, mae=mae, wall_seconds=time.time() - t_l + gather_time))

    # Online learn at each half_life
    for hl in HALF_LIVES:
        t_l = time.time()
        m = online_learn(samples, N_QUBITS, k, et, obs_l1, half_life=float(hl))
        op = m.as_sparse_pauli_op()
        preds = np.array([float(np.real(s.expectation_value(op))) for s in test_states])
        mae = float(np.mean(np.abs(preds - truths)))
        rows.append(CellRow(eta=eta, instance=instance, method="online",
                            half_life=float(hl), mae=mae,
                            wall_seconds=time.time() - t_l + gather_time))

    return rows


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--instances", type=int, default=DEFAULT_INSTANCES)
    p.add_argument("--workers", type=int, default=6)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out = ResumableCSV(
        path=RESULTS_DIR / "raw" / "13_coherent_drift.csv",
        key_cols=KEY_COLS,
        all_cols=ALL_COLS,
    )
    print(f"Resuming with {len(out.done)} rows already logged.")

    # Trajectory key: (eta, instance). Each trajectory yields multiple rows.
    pending_traj: list[tuple[float, int, int]] = []
    seed_counter = args.seed * 7919
    for eta in DRIFT_RATES:
        for inst in range(args.instances):
            seed_counter += 1
            # Skip if every (method, hl) row for this trajectory exists.
            already = all(
                out.is_done({"eta": eta, "instance": inst, "method": m, "half_life": hl})
                for (m, hl) in [("vanilla", 0.0), *(("online", float(h)) for h in HALF_LIVES)]
            )
            if not already:
                pending_traj.append((eta, inst, seed_counter))

    print(f"Pending trajectories: {len(pending_traj)}")
    if pending_traj:
        wait_for_ram(min_fraction=0.20)
        wait_for_memory_under_cap(MEM_CAP_BYTES)

        gen = Parallel(n_jobs=args.workers, return_as="generator", backend="loky", verbose=5)(
            delayed(_trajectory_results)(eta, inst, sd) for (eta, inst, sd) in pending_traj
        )
        for batch in gen:
            for row in batch:
                key = {"eta": row.eta, "instance": row.instance,
                       "method": row.method, "half_life": row.half_life}
                if not out.is_done(key):
                    out.append(asdict(row))

    # ----- aggregate -----
    rows = ResumableCSV.read_all(out.path)

    def collect(method: str, half_life: float | None = None) -> dict[float, list[float]]:
        by_eta: dict[float, list[float]] = {}
        for r in rows:
            if r["method"] != method:
                continue
            if half_life is not None and float(r["half_life"]) != float(half_life):
                continue
            by_eta.setdefault(float(r["eta"]), []).append(float(r["mae"]))
        return by_eta

    vanilla = collect("vanilla")
    online_by_hl = {hl: collect("online", hl) for hl in HALF_LIVES}

    # Find best half-life per eta (lowest MAE).
    best_online_mae: dict[float, float] = {}
    best_hl: dict[float, float] = {}
    for eta in DRIFT_RATES:
        best = (np.inf, None)
        for hl in HALF_LIVES:
            vals = online_by_hl[hl].get(eta, [])
            if not vals:
                continue
            m = float(np.mean(vals))
            if m < best[0]:
                best = (m, hl)
        best_online_mae[eta] = best[0]
        best_hl[eta] = best[1] if best[1] is not None else 0.0

    # ----- plot -----
    fig, ax = plt.subplots(figsize=(8.5, 5.0), dpi=150)
    eta_vals = list(DRIFT_RATES)
    vanilla_means = [float(np.mean(vanilla.get(e, [np.nan]))) for e in eta_vals]
    vanilla_stds = [float(np.std(vanilla.get(e, [0.0]))) for e in eta_vals]
    ax.errorbar(eta_vals, vanilla_means, yerr=vanilla_stds, marker="o",
                color="C3", capsize=4, linewidth=1.6, label="vanilla learn")
    cmap = plt.get_cmap("viridis")
    for i, hl in enumerate(HALF_LIVES):
        ys = [float(np.mean(online_by_hl[hl].get(e, [np.nan]))) for e in eta_vals]
        ystd = [float(np.std(online_by_hl[hl].get(e, [0.0]))) for e in eta_vals]
        ax.errorbar(eta_vals, ys, yerr=ystd, marker="^",
                    color=cmap(0.15 + 0.7 * i / max(1, len(HALF_LIVES) - 1)),
                    capsize=3, linewidth=1, label=f"online (hl = {hl})")
    # Best-of-online overlay
    ys_best = [best_online_mae[e] for e in eta_vals]
    ax.plot(eta_vals, ys_best, marker="*", color="C0", linewidth=2,
            markersize=12, label="online (best half_life)")
    ax.set_xlabel("coherent drift rate eta")
    ax.set_ylabel("MAE on the end-channel test set")
    ax.set_title(
        f"Coherent drift tracking (n={N_QUBITS}, N={N_QUERIES} queries, "
        f"{args.instances} instances)"
    )
    ax.grid(True, linestyle=":", alpha=0.4)
    ax.legend(loc="upper left", fontsize=9, ncol=2)
    fig.tight_layout()
    out_pdf = FIGS_DIR / "13_coherent_drift.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"wrote {out.path} and {out_pdf}")
    print("\nMean MAE summary:")
    for e in eta_vals:
        v = vanilla_means[eta_vals.index(e)]
        b = best_online_mae[e]
        rel = (v - b) / v * 100 if v > 0 else 0.0
        print(f"  eta={e:.4f}: vanilla={v:.4f}, online_best={b:.4f} "
              f"(hl={best_hl[e]:.0f}), improvement={rel:+.1f}%")


if __name__ == "__main__":
    main()
