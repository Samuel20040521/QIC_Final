"""Pauli-string observables used throughout the project.

Conventions match Section 2.1 of Wadhwa & Doosti: a Pauli string `P` over `n`
qubits is an element of `{I, X, Y, Z}^{otimes n}`, encoded as a length-`n`
string with the leftmost character acting on qubit 0. The `degree` `|P|` is
the number of non-identity tensor factors.
"""

from __future__ import annotations

from collections.abc import Iterator
from itertools import product

from qiskit.quantum_info import SparsePauliOp


def pauli_z_first(n: int) -> SparsePauliOp:
    """The observable used in paper Fig. 2: Z on qubit 0, identity elsewhere.

    Qiskit Pauli labels are little-endian: the rightmost character is qubit 0.
    So Z on qubit 0 with identity elsewhere is the label `"I...IZ"`.
    """
    label = "I" * (n - 1) + "Z"
    return SparsePauliOp.from_list([(label, 1.0)])


def pauli_string(label: str) -> SparsePauliOp:
    return SparsePauliOp.from_list([(label, 1.0)])


def degree(label: str) -> int:
    return sum(1 for c in label if c != "I")


def enumerate_low_weight_paulis(n: int, k: int) -> Iterator[str]:
    """Yield Pauli labels P over n qubits with degree |P| <= k.

    Total count: sum_{j=0}^{k} C(n, j) * 3^j  (the (3n)^k upper bound the paper
    uses is a loose envelope of this).
    """
    alphabet = ("I", "X", "Y", "Z")
    for label in product(alphabet, repeat=n):
        s = "".join(label)
        if degree(s) <= k:
            yield s


def pauli_l1_norm(observable: SparsePauliOp) -> float:
    """`||O||_{Pauli,1} = sum_P |a_P|`, used in Algorithm 1's threshold."""
    return float(sum(abs(c) for c in observable.coeffs))
