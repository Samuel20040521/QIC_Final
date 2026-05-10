"""End-to-end checks on Algorithm 1."""

from __future__ import annotations

import numpy as np
from qiskit.quantum_info import Operator

from qpsq.algorithm import run_algorithm1
from qpsq.designs import random_stabilizer_product_state
from qpsq.observables import pauli_z_first
from qpsq.oracle import GaussianQPStat, UnitaryProcess


def test_identity_process_dominant_coefficient() -> None:
    """For the identity channel and `O = Z_0`, learn `alpha_{IIZ} ~= 1` and
    predict `tr(O rho)` accurately on test stabilizer states."""
    rng = np.random.default_rng(0)
    n = 3
    process = UnitaryProcess(Operator(np.eye(2**n)))
    oracle = GaussianQPStat(process, rng=rng)
    res = run_algorithm1(
        oracle=oracle,
        n_qubits=n,
        observable=pauli_z_first(n),
        eps=0.3,
        tau=0.05,
        n_samples=2000,
        rng=rng,
    )
    target_label = "I" * (n - 1) + "Z"
    assert abs(res.model.coefficients.get(target_label, 0.0) - 1.0) < 0.2

    # Prediction error on stabilizer test states should be small.
    truths, preds = [], []
    for _ in range(100):
        state = random_stabilizer_product_state(n, rng)
        psi = state.to_statevector()
        truths.append(process.expectation(psi, pauli_z_first(n)))
        preds.append(res.model.predict(psi))
    truths = np.asarray(truths)
    preds = np.asarray(preds)
    assert float(np.mean((truths - preds) ** 2)) < 0.15


def test_x_gate_process_recovers_minus_z() -> None:
    """For `E(rho) = X_0 rho X_0` and `O = Z_0`, `E^*(O) = -Z_0`."""
    rng = np.random.default_rng(1)
    n = 3
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    I = np.eye(2, dtype=complex)
    # Qiskit little-endian: qubit 0 is the rightmost factor.
    U = I
    for _ in range(n - 1):
        U = np.kron(I, U)
    U = np.kron(np.eye(2**(n - 1)), X)  # X on qubit 0
    process = UnitaryProcess(Operator(U))
    oracle = GaussianQPStat(process, rng=rng)
    res = run_algorithm1(
        oracle=oracle,
        n_qubits=n,
        observable=pauli_z_first(n),
        eps=0.3,
        tau=0.05,
        n_samples=2000,
        rng=rng,
    )
    target_label = "I" * (n - 1) + "Z"
    assert abs(res.model.coefficients.get(target_label, 0.0) + 1.0) < 0.2
