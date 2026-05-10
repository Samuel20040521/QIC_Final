"""Tests for `gather_and_learn_streaming`."""

from __future__ import annotations

import numpy as np
from qiskit.quantum_info import Operator

from qpsq.algorithm import (
    epsilon_tilde,
    gather,
    gather_and_learn_streaming,
    k_from_epsilon,
    learn,
)
from qpsq.designs import random_stabilizer_product_state
from qpsq.observables import pauli_l1_norm, pauli_z_first
from qpsq.oracle import GaussianQPStat, UnitaryProcess


def test_streaming_matches_batch_when_pruning_disabled() -> None:
    """With prune_threshold=0 the streaming code reproduces gather + learn
    up to floating-point order-of-operations differences."""
    rng_a = np.random.default_rng(0)
    rng_b = np.random.default_rng(0)
    n = 3
    process = UnitaryProcess(Operator(np.eye(2**n)))
    obs = pauli_z_first(n)
    obs_l1 = pauli_l1_norm(obs)
    k = k_from_epsilon(0.3)
    et = epsilon_tilde(0.3, n, k)

    # Use the same RNG seed for the oracle to make the two paths comparable.
    oracle_a = GaussianQPStat(process, rng=np.random.default_rng(7))
    oracle_b = GaussianQPStat(process, rng=np.random.default_rng(7))
    samples = gather(oracle_a, n, obs, tau=0.05, n_samples=600, rng=rng_a)
    model_batch = learn(samples, n, k, et, obs_l1)

    model_stream = gather_and_learn_streaming(
        oracle=oracle_b, n_qubits=n, observable=obs, tau=0.05, n_total=600,
        k=k, eps_tilde=et, rng=rng_b, chunk_size=200, prune_threshold=0.0,
    )
    # Same set of non-zero labels.
    nonzero_batch = {l for l, v in model_batch.coefficients.items() if v != 0}
    nonzero_stream = {l for l, v in model_stream.coefficients.items() if v != 0}
    assert nonzero_batch == nonzero_stream
    for label in nonzero_batch:
        assert abs(model_batch.coefficients[label] - model_stream.coefficients[label]) < 1e-9


def test_streaming_recovers_x_gate_dominant_coefficient() -> None:
    """For E(rho) = X_0 rho X_0 and O = Z_0, the dominant alpha is for label
    `'I' * (n-1) + 'Z'` with value -1. Streaming should recover this."""
    rng = np.random.default_rng(0)
    n = 3
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    U = np.kron(np.eye(2 ** (n - 1)), X)  # X on qubit 0 in Qiskit little-endian
    process = UnitaryProcess(Operator(U))
    obs = pauli_z_first(n)
    k = k_from_epsilon(0.3)
    et = epsilon_tilde(0.3, n, k)
    oracle = GaussianQPStat(process, rng=rng)
    model = gather_and_learn_streaming(
        oracle=oracle, n_qubits=n, observable=obs, tau=0.05, n_total=2000,
        k=k, eps_tilde=et, rng=rng, chunk_size=400,
    )
    target_label = "I" * (n - 1) + "Z"
    assert abs(model.coefficients.get(target_label, 0.0) + 1.0) < 0.2
