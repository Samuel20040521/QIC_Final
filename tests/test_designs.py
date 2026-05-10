"""Endian conventions and stabilizer-state utilities."""

from __future__ import annotations

import numpy as np
import pytest
from qiskit.quantum_info import DensityMatrix, Pauli

from qpsq.designs import (
    pauli_expectation_on_stabilizer,
    random_stabilizer_product_state,
)


@pytest.mark.parametrize("seed", [0, 1, 7, 42])
def test_fast_expectation_matches_qiskit(seed: int) -> None:
    """`pauli_expectation_on_stabilizer` agrees with `DensityMatrix.expectation_value`."""
    rng = np.random.default_rng(seed)
    for _ in range(10):
        n = int(rng.integers(2, 6))
        state = random_stabilizer_product_state(n, rng)
        rho = DensityMatrix(state.to_statevector())
        for _ in range(8):
            label = "".join(rng.choice(list("IXYZ"), size=n))
            fast = pauli_expectation_on_stabilizer(label, state)
            slow = float(np.real(rho.expectation_value(Pauli(label))))
            assert abs(fast - slow) < 1e-9, (label, state, fast, slow)


def test_stabilizer_factor_ranges() -> None:
    """`pauli_expectation_on_stabilizer` returns -1, 0, or +1."""
    rng = np.random.default_rng(0)
    seen = set()
    for _ in range(500):
        n = 3
        state = random_stabilizer_product_state(n, rng)
        label = "".join(rng.choice(list("IXYZ"), size=n))
        seen.add(pauli_expectation_on_stabilizer(label, state))
    assert seen.issubset({-1, 0, 1})
    # Over 500 random trials we expect to see all three values.
    assert seen == {-1, 0, 1}
