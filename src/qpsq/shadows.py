"""Classical shadow tomography with randomized Pauli measurements.

Used by Figure 1 to compare the Gaussian-noise emulator of `QPStat` against
a finite-shot classical-shadow estimator. Conventions match Qiskit
(little-endian: qubit 0 is the LSB; Pauli labels read right-to-left).

Single-qubit inverse-channel formula `3 |s><s| - I` follows Huang, Kueng &
Preskill (Nat. Phys. 2020).
"""

from __future__ import annotations

import numpy as np
from qiskit.quantum_info import DensityMatrix, Operator, SparsePauliOp

from .designs import _AXES

# Single-qubit basis-rotation matrices: map the chosen Pauli basis onto Z.
_H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
_SDG = np.array([[1, 0], [0, -1j]], dtype=complex)
_BASIS_ROT = {"X": _H, "Y": _H @ _SDG, "Z": np.eye(2, dtype=complex)}


def _basis_change_unitary(axes: tuple[str, ...]) -> Operator:
    """U such that measuring `U rho U^dag` in the Z basis = measuring `rho` in
    the per-qubit Pauli basis described by `axes`.

    Qiskit qubit ordering: qubit 0 is LSB, so qubit 0's matrix is the
    rightmost factor in the tensor product.
    """
    n = len(axes)
    mat = _BASIS_ROT[axes[n - 1]]
    for i in range(n - 2, -1, -1):
        mat = np.kron(mat, _BASIS_ROT[axes[i]])
    return Operator(mat)


def random_pauli_shadow(
    rho: DensityMatrix, n_shadows: int, rng: np.random.Generator
) -> list[tuple[tuple[str, ...], tuple[int, ...]]]:
    """Collect `n_shadows` randomized Pauli snapshots.

    Each snapshot is `(axes, outcomes)` where `axes[i]` is the measurement
    basis on qubit `i` and `outcomes[i]` is the eigenvalue (+/- 1) measured
    on qubit `i`.
    """
    n = rho.num_qubits
    shadows: list[tuple[tuple[str, ...], tuple[int, ...]]] = []
    for _ in range(n_shadows):
        axes = tuple(_AXES[rng.integers(0, 3)] for _ in range(n))
        rotated = rho.evolve(_basis_change_unitary(axes))
        probs = np.real(np.diag(rotated.data))
        probs = np.clip(probs, 0.0, None)
        if probs.sum() <= 0:
            continue
        probs /= probs.sum()
        idx = int(rng.choice(2**n, p=probs))
        outcomes = tuple(+1 if ((idx >> i) & 1) == 0 else -1 for i in range(n))
        shadows.append((axes, outcomes))
    return shadows


def shadow_expectation(
    shadows: list[tuple[tuple[str, ...], tuple[int, ...]]],
    pauli_label: str,
) -> float:
    """Estimator for `tr(P rho)` from a classical-shadow snapshot set.

    `pauli_label` is in Qiskit's little-endian form (rightmost = qubit 0).
    A snapshot contributes `prod_i 3 * outcome_i` over the qubits where `P`
    is non-identity, but only if the snapshot's basis on those qubits matches
    `P`. Otherwise it contributes zero.
    """
    if not shadows:
        return 0.0
    samples: list[float] = []
    for axes, outcomes in shadows:
        contribution = 1.0
        ok = True
        for p_char, axis, out in zip(reversed(pauli_label), axes, outcomes, strict=True):
            if p_char == "I":
                continue
            if p_char != axis:
                ok = False
                break
            contribution *= 3.0 * out
        samples.append(contribution if ok else 0.0)
    return float(np.mean(samples))


def shadow_evolution_estimator(
    process_evolved_density: DensityMatrix,
    observable: SparsePauliOp,
    n_shadows: int,
    rng: np.random.Generator,
) -> float:
    """Classical-shadow estimator of `tr(O E(|psi><psi|))` from `n_shadows` snapshots."""
    shadows = random_pauli_shadow(process_evolved_density, n_shadows, rng)
    total = 0.0
    for label, coeff in zip(observable.paulis.to_labels(), observable.coeffs, strict=True):
        total += float(np.real(coeff)) * shadow_expectation(shadows, label)
    return total
