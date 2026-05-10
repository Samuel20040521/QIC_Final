"""Quantum processes and the QPStat oracle (paper Section 3 + 4.3).

The oracle interface follows Definition 17: `QPStat_E(rho, O, tau)` returns a
real number `alpha` with `|alpha - tr(O E(rho))| <= tau`. We provide two
emulators:

  * `GaussianQPStat`: computes the exact expectation analytically and adds
    Gaussian noise calibrated so that the deviation stays within tau with
    user-specified confidence (1 - delta). This matches the paper's own
    simulation methodology (Section 4.3, Fig 1).
  * `ShotQPStat`: a finite-shot Pauli-basis estimator; closer to a real
    device. Used for sanity checks.

Both emulators take a `QuantumProcess`, which abstracts over a unitary, a
unitary plus per-qubit noise, or any Kraus channel.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np
from qiskit.quantum_info import DensityMatrix, Kraus, Operator, SparsePauliOp, Statevector
from scipy.stats import norm

# ----- Quantum processes -----------------------------------------------------


class QuantumProcess(ABC):
    """A CPTP map `E: rho -> E(rho)` we can query for expectation values."""

    @abstractmethod
    def expectation(self, state: Statevector, observable: SparsePauliOp) -> float:
        """Return the real part of `tr(O E(rho))` where `rho = |psi><psi|`."""

    @abstractmethod
    def evolve_density(self, rho: DensityMatrix) -> DensityMatrix:
        """Return `E(rho)` as a `DensityMatrix` (used by `ShotQPStat`)."""


@dataclass
class UnitaryProcess(QuantumProcess):
    unitary: Operator

    def expectation(self, state: Statevector, observable: SparsePauliOp) -> float:
        evolved = state.evolve(self.unitary)
        return float(np.real(evolved.expectation_value(observable)))

    def evolve_density(self, rho: DensityMatrix) -> DensityMatrix:
        return rho.evolve(self.unitary)


@dataclass
class KrausProcess(QuantumProcess):
    """A general CPTP channel given by Kraus operators acting on all `n` qubits."""

    kraus: Kraus

    def expectation(self, state: Statevector, observable: SparsePauliOp) -> float:
        rho = DensityMatrix(state).evolve(self.kraus)
        return float(np.real(rho.expectation_value(observable)))

    def evolve_density(self, rho: DensityMatrix) -> DensityMatrix:
        return rho.evolve(self.kraus)


def compose_with_per_qubit_kraus(unitary: Operator, per_qubit_kraus: Kraus) -> "_NoisyUnitary":
    """A unitary followed by independent identical noise on each qubit."""
    return _NoisyUnitary(unitary=unitary, per_qubit_kraus=per_qubit_kraus)


@dataclass
class _NoisyUnitary(QuantumProcess):
    unitary: Operator
    per_qubit_kraus: Kraus

    def expectation(self, state: Statevector, observable: SparsePauliOp) -> float:
        rho = DensityMatrix(state).evolve(self.unitary)
        n = state.num_qubits
        for q in range(n):
            rho = rho.evolve(self.per_qubit_kraus, qargs=[q])
        return float(np.real(rho.expectation_value(observable)))

    def evolve_density(self, rho: DensityMatrix) -> DensityMatrix:
        rho = rho.evolve(self.unitary)
        n = rho.num_qubits
        for q in range(n):
            rho = rho.evolve(self.per_qubit_kraus, qargs=[q])
        return rho


# ----- QPStat oracle ---------------------------------------------------------


class QPStat(ABC):
    """`QPStat_E(rho, O, tau) -> alpha` with `|alpha - tr(O E(rho))| <= tau` w.h.p."""

    @abstractmethod
    def query(self, state: Statevector, observable: SparsePauliOp, tau: float) -> float: ...


@dataclass
class GaussianQPStat(QPStat):
    """Add Gaussian noise calibrated so |error| <= tau w.p. 1 - delta.

    Standard deviation: `sigma = tau / Phi^{-1}(1 - delta/2)`.
    For delta = 0.0455 (the paper's choice in Fig 1), `Phi^{-1}` ~= 2.0,
    so sigma ~= tau / 2.
    """

    process: QuantumProcess
    delta: float = 0.0455
    rng: np.random.Generator | None = None

    def __post_init__(self) -> None:
        if not 0 < self.delta < 1:
            raise ValueError("delta must lie in (0, 1)")
        self._z = float(norm.ppf(1.0 - self.delta / 2.0))

    def query(self, state: Statevector, observable: SparsePauliOp, tau: float) -> float:
        true_val = self.process.expectation(state, observable)
        sigma = tau / self._z
        rng = self.rng if self.rng is not None else np.random.default_rng()
        return true_val + float(rng.normal(0.0, sigma))


@dataclass
class ShotQPStat(QPStat):
    """Finite-shot oracle: decompose `O` into Paulis, sample each in its basis."""

    process: QuantumProcess
    shots: int = 1024
    rng: np.random.Generator | None = None

    def query(self, state: Statevector, observable: SparsePauliOp, tau: float) -> float:
        # The `tau` argument is informational here: the oracle's accuracy is
        # whatever the shot count gives us. Callers should set `shots` large
        # enough that the empirical deviation is below tau with high prob.
        del tau
        rho = self.process.evolve_density(DensityMatrix(state))
        rng = self.rng if self.rng is not None else np.random.default_rng()
        accumulator = 0.0
        for label, coeff in zip(observable.paulis.to_labels(), observable.coeffs, strict=True):
            mean = float(np.real(rho.expectation_value(SparsePauliOp.from_list([(label, 1.0)]))))
            mean = max(-1.0, min(1.0, mean))
            p_plus = (1.0 + mean) / 2.0
            outcomes = rng.binomial(1, p_plus, size=self.shots)  # +1 if 1, -1 if 0
            est = (2.0 * outcomes.mean() - 1.0)
            accumulator += float(np.real(coeff)) * est
        return accumulator
