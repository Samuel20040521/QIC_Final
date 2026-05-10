"""QPStat tolerance and reproducibility."""

from __future__ import annotations

import numpy as np
from qiskit.quantum_info import Operator

from qpsq.designs import random_stabilizer_product_state
from qpsq.observables import pauli_z_first
from qpsq.oracle import GaussianQPStat, ShotQPStat, UnitaryProcess


def test_gaussian_oracle_within_tolerance_at_high_confidence() -> None:
    """Empirical |error| < tau in >= (1 - delta) fraction of queries."""
    rng = np.random.default_rng(0)
    n = 4
    process = UnitaryProcess(Operator(np.eye(2**n)))
    delta = 0.0455
    oracle = GaussianQPStat(process, delta=delta, rng=rng)
    tau = 0.1
    obs = pauli_z_first(n)
    n_trials = 5000
    deviations = []
    for _ in range(n_trials):
        state = random_stabilizer_product_state(n, rng)
        truth = process.expectation(state.to_statevector(), obs)
        est = oracle.query(state.to_statevector(), obs, tau)
        deviations.append(abs(est - truth))
    deviations = np.asarray(deviations)
    fraction_within = float(np.mean(deviations <= tau))
    # Should comfortably exceed (1 - delta).
    assert fraction_within >= (1 - delta) - 0.01, (fraction_within, 1 - delta)


def test_gaussian_oracle_reproducible_with_seed() -> None:
    rng_a = np.random.default_rng(123)
    rng_b = np.random.default_rng(123)
    n = 3
    process = UnitaryProcess(Operator(np.eye(2**n)))
    obs = pauli_z_first(n)
    oa = GaussianQPStat(process, rng=rng_a)
    ob = GaussianQPStat(process, rng=rng_b)
    rng_state = np.random.default_rng(0)
    state = random_stabilizer_product_state(n, rng_state)
    out_a = [oa.query(state.to_statevector(), obs, 0.1) for _ in range(20)]
    out_b = [ob.query(state.to_statevector(), obs, 0.1) for _ in range(20)]
    assert out_a == out_b


def test_shot_oracle_recovers_expectation() -> None:
    rng = np.random.default_rng(0)
    n = 3
    process = UnitaryProcess(Operator(np.eye(2**n)))
    oracle = ShotQPStat(process, shots=8000, rng=rng)
    state = random_stabilizer_product_state(n, rng)
    truth = process.expectation(state.to_statevector(), pauli_z_first(n))
    est = oracle.query(state.to_statevector(), pauli_z_first(n), tau=0.05)
    assert abs(est - truth) < 0.05
