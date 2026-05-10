"""Random state and unitary samplers (paper Section 2.2 + 4.3).

A `StabilizerProductState` is a tensor product of single-qubit Pauli
eigenstates, i.e. an element of `stab_1^{otimes n}`. We keep the per-qubit
(axis, sign) representation rather than a 2^n state vector, because:

  * The sampling distribution in Algorithm 1 is uniform over this set.
  * `tr(P |psi><psi|)` for any Pauli string `P` reduces to a product of
    single-qubit factors in `{-1, 0, +1}` and never requires 2^n storage.

This matters: Algorithm 1 may issue tens of thousands of queries; computing
expectation values via dense state vectors would be wasteful.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import IGate
from qiskit.quantum_info import Operator, Statevector, random_clifford, random_unitary

# (axis, sign) -> 2-vector
_SINGLE_QUBIT = {
    ("Z", +1): np.array([1.0, 0.0], dtype=complex),                       # |0>
    ("Z", -1): np.array([0.0, 1.0], dtype=complex),                       # |1>
    ("X", +1): np.array([1.0, 1.0], dtype=complex) / np.sqrt(2),          # |+>
    ("X", -1): np.array([1.0, -1.0], dtype=complex) / np.sqrt(2),         # |->
    ("Y", +1): np.array([1.0, 1j], dtype=complex) / np.sqrt(2),           # |+i>
    ("Y", -1): np.array([1.0, -1j], dtype=complex) / np.sqrt(2),          # |-i>
}

_AXES = ("X", "Y", "Z")
_SIGNS = (+1, -1)


@dataclass(frozen=True)
class StabilizerProductState:
    """`|psi> = otimes_i |s_i>`, with each `|s_i>` an eigenstate of a single Pauli.

    `axes[i]` and `signs[i]` describe qubit `i`. Qiskit uses little-endian
    state-vector indexing (qubit 0 = LSB), so `to_statevector` builds the
    tensor product with high-index qubits first.
    """

    axes: tuple[str, ...]   # one of "X", "Y", "Z" per qubit
    signs: tuple[int, ...]  # +1 or -1 per qubit

    @property
    def n(self) -> int:
        return len(self.axes)

    def to_statevector(self) -> Statevector:
        vec = np.array([1.0 + 0j])
        for i in range(self.n - 1, -1, -1):
            vec = np.kron(vec, _SINGLE_QUBIT[(self.axes[i], self.signs[i])])
        return Statevector(vec)


def random_stabilizer_product_state(n: int, rng: np.random.Generator) -> StabilizerProductState:
    """Sample uniformly from `stab_1^{otimes n}` (6^n states)."""
    axes = tuple(_AXES[rng.integers(0, 3)] for _ in range(n))
    signs = tuple(int(_SIGNS[rng.integers(0, 2)]) for _ in range(n))
    return StabilizerProductState(axes=axes, signs=signs)


def pauli_expectation_on_stabilizer(label: str, state: StabilizerProductState) -> int:
    """`tr(P |psi><psi|)` for a Qiskit-form Pauli label `P` and a stabilizer
    product state.

    Returns an integer in `{-1, 0, +1}`: a product of single-qubit factors. A
    factor is 0 iff the Pauli on that qubit disagrees with the qubit's
    stabilizer axis (and is non-identity).

    Qiskit Pauli labels are little-endian: position 0 of the string is qubit
    `n-1`, position `n-1` is qubit 0. We zip the reversed label with `axes`
    (which is indexed by qubit number) so the right qubit is paired with the
    right Pauli.
    """
    if len(label) != state.n:
        raise ValueError(f"Pauli label length {len(label)} != n_qubits {state.n}")
    factor = 1
    for p, axis, sign in zip(reversed(label), state.axes, state.signs, strict=True):
        if p == "I":
            continue
        if p == axis:
            factor *= sign
        else:
            return 0
    return factor


def haar_unitary(n: int, rng: np.random.Generator) -> Operator:
    """Wrap qiskit's Haar sampler with a numpy `Generator` seed."""
    seed = int(rng.integers(0, 2**32 - 1))
    return Operator(random_unitary(2**n, seed=seed))


def random_clifford_unitary(n: int, rng: np.random.Generator) -> Operator:
    seed = int(rng.integers(0, 2**32 - 1))
    return Operator(random_clifford(n, seed=seed))


def clifford_plus_t_circuit(
    n: int,
    depth: int,
    t_count: int,
    rng: np.random.Generator,
) -> Operator:
    """A Clifford+T circuit with `depth` random Cliffords interleaved with
    `t_count` T-gates distributed at random slots.

    The circuit has the schematic form

        C_0  T...  C_1  T...  C_2  ...  C_{depth-1}  T...

    where T-gates are inserted at uniformly chosen "slots" (`depth + 1`
    positions: before each Clifford and after the last). The T-gate's
    qubit is sampled uniformly per gate. With `t_count = 0` the circuit
    is a pure Clifford (a unitary 3-design when `depth >= 1`). As
    `t_count` grows, magic accumulates and the operator distribution
    departs from any unitary 2-design.

    Reproducible: all randomness flows through `rng`.
    """
    if depth < 0 or t_count < 0:
        raise ValueError("depth and t_count must be non-negative")
    qc = QuantumCircuit(n)
    if depth == 0 and t_count == 0:
        return Operator(qc)

    # Sample T-gate placements once.
    if t_count > 0:
        t_slots = rng.integers(0, depth + 1, size=t_count)
        t_qubits = rng.integers(0, n, size=t_count)
    else:
        t_slots = np.array([], dtype=int)
        t_qubits = np.array([], dtype=int)
    slot_t_indices: list[list[int]] = [[] for _ in range(depth + 1)]
    for i, s in enumerate(t_slots):
        slot_t_indices[int(s)].append(i)

    for slot in range(depth + 1):
        for i in slot_t_indices[slot]:
            qc.t(int(t_qubits[i]))
        if slot < depth:
            seed = int(rng.integers(0, 2**32 - 1))
            cliff = random_clifford(n, seed=seed)
            qc.compose(cliff.to_circuit(), inplace=True)

    return Operator(qc)


def brickwall_random_circuit(n: int, depth: int, rng: np.random.Generator) -> Operator:
    """A brick-wall random circuit: alternating layers of Haar-random 2-qubit gates.

    Even layers couple (0,1), (2,3), ...; odd layers couple (1,2), (3,4), ....
    `depth` counts total layers. Useful as a tunable approximate t-design knob.
    """
    qc = QuantumCircuit(n)
    for layer in range(depth):
        offset = layer % 2
        for i in range(offset, n - 1, 2):
            seed = int(rng.integers(0, 2**32 - 1))
            two_qubit = random_unitary(4, seed=seed)
            qc.unitary(two_qubit, [i, i + 1])
        # Apply identity to leftover qubit so the circuit still spans n wires
        if offset == 1 and n % 2 == 0:
            qc.append(IGate(), [n - 1])
    return Operator(qc)
