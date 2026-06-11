"""Reproduce paper Figure 1.

Compare two simulators of the QPStat oracle on the same task:

  * `GaussianQPStat`: adds normal noise with sigma calibrated so the
    deviation lies within tau with probability 1 - delta.
  * Classical-shadow estimator (Huang-Kueng-Preskill): finite-shot Pauli
    measurement averages of `tr(O E(rho))`.

Setup matches the paper: a single-qubit fixed Haar unitary, observable
`O = Z`, tau = 0.2, delta = 0.0455, evaluated on random single-qubit
stabilizer states. The paper reports that the shadow method tends to
produce smaller deviations than the calibrated Gaussian, so the paper's
emulator is a worst-case proxy.

The shadow shot budget is NOT a free parameter: it is derived from the
classical-shadow sample-complexity bound (Chebyshev) so that the shadow
estimate satisfies the same `(tau, delta)` guarantee as the Gaussian
oracle:

    n_shots = ceil(3^w * (sum_P |c_P|)^2 / (delta * tau^2)) = 1649

for the weight-1 observable `Z`. Because this bound is loose, the actual
shadow deviations come out visibly narrower than the calibrated Gaussian
— which is exactly the point of the paper's Figure 1. (Choosing an ad-hoc
shot count breaks the comparison: e.g. 256 shots gives a shadow std of
~0.108 ~= the Gaussian's sigma = 0.1, and the histograms overlap.)

Outputs:
  experiments/results/fig1.csv     — per-trial deviations
  report/figs/fig1_reproduction.pdf — overlaid histograms
"""

from __future__ import annotations

import csv
import math

import matplotlib.pyplot as plt
import numpy as np
from qiskit.quantum_info import DensityMatrix, Statevector

from qpsq.designs import haar_unitary, random_stabilizer_product_state
from qpsq.observables import pauli_z_first
from qpsq.oracle import GaussianQPStat, UnitaryProcess
from qpsq.shadows import shadow_evolution_estimator

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

TAU = 0.2
DELTA = 0.0455
N_TRIALS = 1000  # matches the paper's N = 1000
N_QUBITS = 1
# Weight of O and Pauli-1-norm of O (O = Z: weight 1, sum of |coeffs| = 1).
OBS_WEIGHT = 1
OBS_COEFF_L1 = 1.0
# Shadow shot budget from the (tau, delta) guarantee — see module docstring.
N_SHADOWS = math.ceil(3**OBS_WEIGHT * OBS_COEFF_L1**2 / (DELTA * TAU**2))


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)

    # Fix one Haar-random 1-qubit unitary for the whole experiment.
    U = haar_unitary(N_QUBITS, rng)
    process = UnitaryProcess(U)
    obs = pauli_z_first(N_QUBITS)
    oracle = GaussianQPStat(process, delta=DELTA, rng=rng)

    gauss_devs: list[float] = []
    shadow_devs: list[float] = []

    for _ in range(N_TRIALS):
        state = random_stabilizer_product_state(N_QUBITS, rng)
        psi: Statevector = state.to_statevector()
        truth = process.expectation(psi, obs)

        gauss_est = oracle.query(psi, obs, TAU)
        gauss_devs.append(gauss_est - truth)

        evolved = process.evolve_density(DensityMatrix(psi))
        shadow_est = shadow_evolution_estimator(evolved, obs, N_SHADOWS, rng)
        shadow_devs.append(shadow_est - truth)

    gauss_devs_arr = np.asarray(gauss_devs)
    shadow_devs_arr = np.asarray(shadow_devs)

    csv_path = RESULTS_DIR / "fig1.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["trial", "gauss_deviation", "shadow_deviation"])
        for i, (g, s) in enumerate(zip(gauss_devs_arr, shadow_devs_arr, strict=True)):
            w.writerow([i, f"{g:.6f}", f"{s:.6f}"])

    print(
        f"Gaussian: |deviation| < tau in "
        f"{(np.abs(gauss_devs_arr) < TAU).mean():.3%} of trials. "
        f"std = {gauss_devs_arr.std():.4f}"
    )
    print(
        f"Shadow ({N_SHADOWS} shots): |deviation| < tau in "
        f"{(np.abs(shadow_devs_arr) < TAU).mean():.3%} of trials. "
        f"std = {shadow_devs_arr.std():.4f}"
    )

    # Histogram style mirrors the paper's Figure 1: 35 bins over (-0.35, 0.35),
    # classical shadow drawn first (narrow), normal-distribution oracle second.
    fig, ax = plt.subplots(figsize=(6, 4.5))
    bins = np.linspace(-0.35, 0.35, 36)
    ax.hist(
        shadow_devs_arr,
        bins=bins,
        alpha=0.75,
        label="Classical Shadow",
        color="C0",
        edgecolor="black",
        linewidth=0.4,
    )
    ax.hist(
        gauss_devs_arr,
        bins=bins,
        alpha=0.75,
        label="Normal Distr",
        color="C1",
        edgecolor="black",
        linewidth=0.4,
    )
    ax.set_xlabel("error")
    ax.set_ylabel("Count")
    ax.set_title("Comparing deviations using two methods")
    ax.legend(loc="upper right")
    fig.tight_layout()
    out = FIGS_DIR / "fig1_reproduction.pdf"
    fig.savefig(out)
    print(f"wrote {csv_path} and {out}")


if __name__ == "__main__":
    main()
