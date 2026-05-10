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
from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np
from qiskit.circuit import QuantumCircuit
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


@dataclass
class NoisyUnitaryProcess(QuantumProcess):
    """A unitary followed by independent identical Kraus noise on each qubit."""

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


# Backwards-compatible alias for older code paths (experiment 03).
_NoisyUnitary = NoisyUnitaryProcess


def compose_with_per_qubit_kraus(
    unitary: Operator, per_qubit_kraus: Kraus
) -> NoisyUnitaryProcess:
    """A unitary followed by independent identical noise on each qubit."""
    return NoisyUnitaryProcess(unitary=unitary, per_qubit_kraus=per_qubit_kraus)


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


# ----- Drifting (non-Markovian) process ------------------------------------


class DriftingProcess(QuantumProcess):
    """Wrap a base process whose parameters drift with query index `t`.

    The user supplies a `process_at(t: int) -> QuantumProcess` factory
    callable; each call to `expectation` or `evolve_density` advances `t`
    by one and constructs the process for that index. This deliberately
    keeps the drift model out of the wrapper: any time-dependent channel
    (drifting unitary, drifting depolarizing rate, etc.) plugs in as a
    factory.

    Reproducibility: the factory is responsible for any randomness it
    needs; the wrapper only tracks `t`.
    """

    def __init__(self, process_at: Callable[[int], QuantumProcess], start_t: int = 0) -> None:
        if not callable(process_at):
            raise TypeError("process_at must be callable: int -> QuantumProcess")
        self._process_at = process_at
        self._t = int(start_t)

    @property
    def t(self) -> int:
        return self._t

    def reset(self, t: int = 0) -> None:
        self._t = int(t)

    def expectation(self, state: Statevector, observable: SparsePauliOp) -> float:
        proc = self._process_at(self._t)
        out = proc.expectation(state, observable)
        self._t += 1
        return out

    def evolve_density(self, rho: DensityMatrix) -> DensityMatrix:
        proc = self._process_at(self._t)
        out = proc.evolve_density(rho)
        self._t += 1
        return out


def linear_depolarizing_drift(
    base_unitary: Operator,
    p_start: float,
    p_end: float,
    n_total_queries: int,
) -> Callable[[int], QuantumProcess]:
    """Factory: per-qubit depolarizing rate sweeps linearly from `p_start`
    to `p_end` over `n_total_queries` calls.

    Use case: a slow drift in 1Q decoherence rate over the experiment.
    """
    if n_total_queries <= 0:
        raise ValueError("n_total_queries must be positive")

    I = np.eye(2, dtype=complex)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    Z = np.array([[1, 0], [0, -1]], dtype=complex)

    def _kraus_at(p: float) -> Kraus:
        p = max(0.0, min(1.0, p))
        return Kraus([
            np.sqrt(1 - 3 * p / 4) * I,
            np.sqrt(p / 4) * X,
            np.sqrt(p / 4) * Y,
            np.sqrt(p / 4) * Z,
        ])

    def factory(t: int) -> QuantumProcess:
        frac = min(1.0, max(0.0, t / max(1, n_total_queries - 1)))
        p_t = p_start + frac * (p_end - p_start)
        return NoisyUnitaryProcess(unitary=base_unitary, per_qubit_kraus=_kraus_at(p_t))

    return factory


# ----- Hardware-realistic circuit + Aer NoiseModel -------------------------


@dataclass
class CircuitNoiseProcess(QuantumProcess):
    """Run a `QuantumCircuit` through Qiskit Aer's density-matrix simulator
    with an attached `NoiseModel`. This is the most hardware-realistic
    process emulator we expose: 1Q/2Q gate errors, readout errors, and any
    other Aer-supported noise propagate end-to-end.

    The circuit must be fully bound (no free parameters). The number of
    qubits is taken from the circuit; input states must match.

    Reuses one `AerSimulator` instance per process to amortize setup cost.
    """

    circuit: QuantumCircuit
    noise_model: object = None  # qiskit_aer.noise.NoiseModel; optional
    _sim: object = field(init=False, default=None, repr=False)

    def __post_init__(self) -> None:
        # Lazy-import to avoid a hard dep at module-import time.
        from qiskit_aer import AerSimulator

        kwargs = {"method": "density_matrix"}
        if self.noise_model is not None:
            kwargs["noise_model"] = self.noise_model
        self._sim = AerSimulator(**kwargs)
        if self.circuit.num_parameters != 0:
            raise ValueError(
                f"circuit must be bound (got {self.circuit.num_parameters} free params)"
            )

    def _run(self, init_dm: DensityMatrix) -> DensityMatrix:
        n = self.circuit.num_qubits
        if init_dm.num_qubits != n:
            raise ValueError(
                f"input density matrix has {init_dm.num_qubits} qubits, circuit has {n}"
            )
        qc = QuantumCircuit(n)
        qc.set_density_matrix(init_dm)
        qc.compose(self.circuit, inplace=True)
        qc.save_density_matrix()
        result = self._sim.run(qc).result()
        return result.data(0).get("density_matrix")

    def expectation(self, state: Statevector, observable: SparsePauliOp) -> float:
        rho = self._run(DensityMatrix(state))
        return float(np.real(rho.expectation_value(observable)))

    def evolve_density(self, rho: DensityMatrix) -> DensityMatrix:
        return self._run(rho)


def default_hardware_noise_model(
    p1: float = 0.001,
    p2: float = 0.01,
    readout_p: float = 0.02,
):
    """A pragmatic default: 1Q depolarizing on `u`, `rz`, `sx`, `x`, `h`,
    `t`, `tdg`; 2Q depolarizing on `cx`, `cz`; symmetric readout error.

    Returns a `qiskit_aer.noise.NoiseModel` ready to attach to
    `CircuitNoiseProcess`.
    """
    from qiskit_aer.noise import (
        NoiseModel,
        ReadoutError,
        depolarizing_error,
    )

    nm = NoiseModel()
    one_q = depolarizing_error(p1, 1)
    two_q = depolarizing_error(p2, 2)
    nm.add_all_qubit_quantum_error(one_q, ["u", "u1", "u2", "u3", "rz", "sx", "x", "h", "t", "tdg", "ry", "rx"])
    nm.add_all_qubit_quantum_error(two_q, ["cx", "cz"])
    if readout_p > 0:
        nm.add_all_qubit_readout_error(
            ReadoutError([[1 - readout_p, readout_p], [readout_p, 1 - readout_p]])
        )
    return nm
