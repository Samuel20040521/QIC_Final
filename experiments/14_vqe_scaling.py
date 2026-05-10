"""Sprint Task 4: VQE surrogate mitigation at n in {6, 8}.

Same factored-noise design as `08_ext_vqe_mitigation.py`, but at larger
n. We keep the QPSQ training one-shot and the surrogate fully closed-form
at evaluation time.

We compare:
  * noiseless reference (truth),
  * noisy_factored (`Lambda` after `U(theta)`),
  * surrogate (closed form),
  * shadow baseline (K = `SHADOW_BUDGET` snapshots per evaluation).

Resumable per-(n, theta_idx, method) row.
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
from qiskit.circuit import ParameterVector, QuantumCircuit
from qiskit.quantum_info import (
    DensityMatrix,
    Kraus,
    Operator,
    SparsePauliOp,
    Statevector,
)
from scipy.optimize import minimize

from qpsq.algorithm import (
    epsilon_tilde,
    gather,
    k_from_epsilon,
    learn,
)
from qpsq.observables import pauli_l1_norm
from qpsq.oracle import (
    GaussianQPStat,
    NoisyUnitaryProcess,
)
from qpsq.shadows import shadow_evolution_estimator

from _common import FIGS_DIR, RESULTS_DIR
from _hpc import wait_for_memory_under_cap, wait_for_ram
from _resumable import ResumableCSV

J, H_FIELD = 1.0, 0.5
ANSATZ_LAYERS = 2
N_QPSQ_SAMPLES = 4000
TAU = 0.05
EPS = 0.3
N_THETA_NEIGHBORS = 24
THETA_RADIUS = 0.8
SHADOW_BUDGET = 200
P_DEPOL = 0.04
DEFAULT_NS = (6, 8)
MEM_CAP_BYTES = 70 * 1024**3

KEY_COLS = ["n", "theta_idx", "method"]
ALL_COLS = KEY_COLS + ["energy", "wall_seconds"]


def tfim_hamiltonian(n: int) -> SparsePauliOp:
    terms: list[tuple[str, float]] = []
    for i in range(n - 1):
        label = ["I"] * n
        label[n - 1 - i] = "Z"
        label[n - 1 - (i + 1)] = "Z"
        terms.append(("".join(label), -J))
    for i in range(n):
        label = ["I"] * n
        label[n - 1 - i] = "X"
        terms.append(("".join(label), -H_FIELD))
    return SparsePauliOp.from_list(terms)


def make_ansatz(n: int):
    n_params = ANSATZ_LAYERS * 2 * n
    theta = ParameterVector("theta", n_params)
    qc = QuantumCircuit(n)
    p_idx = 0
    for _ in range(ANSATZ_LAYERS):
        for q in range(n):
            qc.ry(theta[p_idx], q); p_idx += 1
        for q in range(n):
            qc.rz(theta[p_idx], q); p_idx += 1
        for q in range(n - 1):
            qc.cx(q, q + 1)
    return qc, theta


def bind(qc, theta, vals):
    return qc.assign_parameters(dict(zip(theta, vals, strict=True)))


def depol_kraus(p: float) -> Kraus:
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


@dataclass
class Row:
    n: int
    theta_idx: int
    method: str
    energy: float
    wall_seconds: float


def _scenario_for_n(n: int, seed: int) -> dict:
    """Coarse-VQE + surrogate-train pipeline for one `n`."""
    rng = np.random.default_rng(seed)
    H = tfim_hamiltonian(n)
    obs_l1 = pauli_l1_norm(H)
    qc_template, theta_param = make_ansatz(n)
    lam_kraus = depol_kraus(P_DEPOL)
    lambda_proc = NoisyUnitaryProcess(
        unitary=Operator(np.eye(2**n)), per_qubit_kraus=lam_kraus,
    )

    rng_init = np.random.default_rng(seed + 1)
    theta0 = rng_init.uniform(-0.5, 0.5, size=qc_template.num_parameters)

    def noisy_E(theta_vals):
        bound = bind(qc_template, theta_param, theta_vals)
        proc = NoisyUnitaryProcess(unitary=Operator(bound), per_qubit_kraus=lam_kraus)
        return proc.expectation(Statevector.from_label("0" * n), H)

    res = minimize(noisy_E, theta0, method="COBYLA",
                   options={"maxiter": 60, "rhobeg": 0.4})
    theta_star = np.asarray(res.x)

    # Train QPSQ on Lambda.
    oracle = GaussianQPStat(lambda_proc, rng=rng)
    k = k_from_epsilon(EPS)
    et = epsilon_tilde(EPS, n, k)
    samples = gather(oracle, n, H, TAU, N_QPSQ_SAMPLES, rng)
    model = learn(samples, n, k, et, obs_l1)
    return {
        "n": n,
        "H": H,
        "qc_template": qc_template,
        "theta_param": theta_param,
        "lam_kraus": lam_kraus,
        "theta_star": theta_star,
        "surrogate_op": model.as_sparse_pauli_op(),
        "rng": rng,
    }


def _evaluate_one(scenario: dict, theta_idx: int, theta_vals: np.ndarray, seed: int) -> list[Row]:
    n = scenario["n"]
    H = scenario["H"]
    qc_template = scenario["qc_template"]
    theta_param = scenario["theta_param"]
    lam_kraus = scenario["lam_kraus"]
    surrogate_op = scenario["surrogate_op"]
    rng = np.random.default_rng(seed)

    rows: list[Row] = []
    bound = bind(qc_template, theta_param, theta_vals)
    U_theta = Operator(bound)
    psi = Statevector.from_label("0" * n).evolve(U_theta)

    t0 = time.time()
    e_noiseless = float(np.real(psi.expectation_value(H)))
    rows.append(Row(n=n, theta_idx=theta_idx, method="noiseless",
                    energy=e_noiseless, wall_seconds=time.time() - t0))

    t0 = time.time()
    proc_noisy = NoisyUnitaryProcess(unitary=U_theta, per_qubit_kraus=lam_kraus)
    e_noisy = proc_noisy.expectation(Statevector.from_label("0" * n), H)
    rows.append(Row(n=n, theta_idx=theta_idx, method="noisy_factored",
                    energy=e_noisy, wall_seconds=time.time() - t0))

    t0 = time.time()
    e_surr = float(np.real(psi.expectation_value(surrogate_op)))
    rows.append(Row(n=n, theta_idx=theta_idx, method="surrogate",
                    energy=e_surr, wall_seconds=time.time() - t0))

    t0 = time.time()
    rho_noisy = proc_noisy.evolve_density(DensityMatrix(Statevector.from_label("0" * n)))
    e_shadow = shadow_evolution_estimator(rho_noisy, H, SHADOW_BUDGET, rng)
    rows.append(Row(n=n, theta_idx=theta_idx, method="shadow",
                    energy=e_shadow, wall_seconds=time.time() - t0))

    return rows


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--ns", type=str, default=",".join(str(n) for n in DEFAULT_NS))
    p.add_argument("--workers", type=int, default=4)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    ns = tuple(int(x) for x in args.ns.split(","))
    out = ResumableCSV(
        path=RESULTS_DIR / "raw" / "14_vqe_scaling.csv",
        key_cols=KEY_COLS,
        all_cols=ALL_COLS,
    )
    print(f"Resuming with {len(out.done)} rows.")

    for n in ns:
        # Build scenario (coarse VQE + QPSQ surrogate). Not parallelised:
        # the gather step is the dominant cost and small N is fast.
        wait_for_ram(min_fraction=0.2)
        wait_for_memory_under_cap(MEM_CAP_BYTES)
        scenario_seed = args.seed * 12345 + n
        print(f"\n--- n = {n}: building surrogate ---")
        t_s = time.time()
        scenario = _scenario_for_n(n, scenario_seed)
        print(f"  scenario built in {time.time() - t_s:.1f}s; "
              f"surrogate has {scenario['surrogate_op'].size} Pauli terms")

        rng_n = np.random.default_rng(scenario_seed + 999)
        n_params = scenario["qc_template"].num_parameters
        deltas = rng_n.uniform(-1, 1, size=(N_THETA_NEIGHBORS, n_params))
        norms = np.linalg.norm(deltas, axis=1, keepdims=True)
        deltas = deltas / np.maximum(norms, 1e-12)
        deltas *= THETA_RADIUS * rng_n.uniform(0.0, 1.0, size=(N_THETA_NEIGHBORS, 1))
        theta_grid = scenario["theta_star"][None, :] + deltas

        # Pending = (n, theta_idx) where any of the 4 methods is missing.
        pending: list[tuple[int, int, np.ndarray, int]] = []
        for j, theta_vals in enumerate(theta_grid):
            keys_missing = [
                {"n": n, "theta_idx": j, "method": m}
                for m in ("noiseless", "noisy_factored", "surrogate", "shadow")
            ]
            if any(not out.is_done(k) for k in keys_missing):
                pending.append((n, j, theta_vals, scenario_seed + 7000 + j))

        print(f"  pending {len(pending)} theta-grid points")
        if not pending:
            continue

        gen = Parallel(n_jobs=args.workers, return_as="generator", backend="loky", verbose=5)(
            delayed(_evaluate_one)(scenario, j, tv, sd) for (_, j, tv, sd) in pending
        )
        for batch in gen:
            for row in batch:
                if not out.is_done({"n": row.n, "theta_idx": row.theta_idx, "method": row.method}):
                    out.append(asdict(row))

    # ----- plot -----
    rows = ResumableCSV.read_all(out.path)
    by_n: dict[int, dict[int, dict[str, float]]] = {}
    for r in rows:
        nv = int(r["n"])
        if nv not in ns:
            continue
        idx = int(r["theta_idx"])
        by_n.setdefault(nv, {}).setdefault(idx, {})[r["method"]] = float(r["energy"])

    fig, axes = plt.subplots(1, len(ns), figsize=(6.5 * len(ns), 5.0), dpi=150,
                             sharey=False)
    if len(ns) == 1:
        axes = [axes]
    for ax, n in zip(axes, ns, strict=True):
        per = by_n.get(n, {})
        if not per:
            ax.set_title(f"n={n} (no data)")
            continue
        idxs = sorted(per.keys())
        truth = np.array([per[i].get("noiseless", np.nan) for i in idxs])
        for method, color, marker in (
            ("noisy_factored", "C3", "o"),
            ("surrogate", "C0", "s"),
            ("shadow", "C2", "^"),
        ):
            preds = np.array([per[i].get(method, np.nan) for i in idxs])
            valid = ~(np.isnan(preds) | np.isnan(truth))
            if not valid.any():
                continue
            mae = float(np.mean(np.abs(preds[valid] - truth[valid])))
            ax.scatter(truth[valid], preds[valid], color=color, marker=marker,
                       s=24, alpha=0.85, label=f"{method} (MAE={mae:.3f})")
        valid = ~np.isnan(truth)
        if valid.any():
            lo, hi = truth[valid].min() - 0.4, truth[valid].max() + 0.4
            ax.plot([lo, hi], [lo, hi], "k--", linewidth=0.7)
            ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel("noiseless reference energy")
        ax.set_ylabel("estimated energy")
        ax.set_title(f"n = {n}")
        ax.grid(True, linestyle=":", alpha=0.4)
        ax.legend(loc="upper left", fontsize=8)
    fig.suptitle(
        f"Surrogate-VQE scaling (TFIM, J={J}, h={H_FIELD}, p_depol={P_DEPOL}, "
        f"N_train={N_QPSQ_SAMPLES}, {N_THETA_NEIGHBORS} theta-neighbors)"
    )
    fig.tight_layout()
    out_pdf = FIGS_DIR / "14_vqe_scaling.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"\nwrote {out.path} and {out_pdf}")


if __name__ == "__main__":
    main()
