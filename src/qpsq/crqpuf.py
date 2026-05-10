"""Classical-Readout Quantum PUFs and the QPSQ-based attack (paper Section 6).

A CR-QPUF wraps a fixed (possibly noisy) quantum process `E`. The
authentication protocol of Section 6.2 has the form:

  * Verifier sends a classical challenge string `c`, which is interpreted as
    a description of an input state `rho_c` (here: a stabilizer product
    state) and an observable `O_c` (here: a fixed Pauli-Z on the first qubit).
  * Prover prepares `rho_c`, applies `E`, measures `O_c`, and returns the
    estimated expectation value `r`.
  * Verifier accepts if `|r - r_enrolled(c)| <= threshold`.

The attack instantiates `E` as the QPSQ oracle on the adversary's side and
runs Algorithm 1 to learn a predictor `h` for `tr(O_c E(rho))`. The adversary
then forges responses by computing `h(rho_c)` for any future challenge.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from qiskit.quantum_info import SparsePauliOp, Statevector

from .algorithm import Algorithm1Result, run_algorithm1
from .designs import StabilizerProductState, random_stabilizer_product_state
from .observables import pauli_z_first
from .oracle import GaussianQPStat, QuantumProcess, UnitaryProcess


@dataclass
class CRQPUF:
    """A device that, on challenge `rho`, produces `tr(O E(rho))` (noisy)."""

    process: QuantumProcess
    observable: SparsePauliOp
    n_qubits: int

    @classmethod
    def from_haar_unitary(cls, n: int, rng: np.random.Generator) -> "CRQPUF":
        from .designs import haar_unitary

        return cls(
            process=UnitaryProcess(haar_unitary(n, rng)),
            observable=pauli_z_first(n),
            n_qubits=n,
        )

    def respond(self, challenge: StabilizerProductState, tau: float, rng: np.random.Generator) -> float:
        """Honest prover response — exact expectation plus query-tolerance noise."""
        oracle = GaussianQPStat(self.process, rng=rng)
        return oracle.query(challenge.to_statevector(), self.observable, tau)


def random_challenge_set(
    n_qubits: int, n_challenges: int, rng: np.random.Generator
) -> list[StabilizerProductState]:
    return [random_stabilizer_product_state(n_qubits, rng) for _ in range(n_challenges)]


@dataclass
class AuthenticationProtocol:
    crpuf: CRQPUF
    accept_threshold: float = 0.1
    response_tau: float = 0.05

    def enroll(
        self, challenges: list[StabilizerProductState], rng: np.random.Generator
    ) -> list[float]:
        """Server collects ground-truth responses to a set of challenges."""
        return [self.crpuf.respond(c, self.response_tau, rng) for c in challenges]

    def verify(
        self,
        challenge: StabilizerProductState,
        claimed_response: float,
        enrolled_response: float,
    ) -> bool:
        return abs(claimed_response - enrolled_response) <= self.accept_threshold


# ----- the attack -----------------------------------------------------------


@dataclass
class AttackResult:
    learned: Algorithm1Result
    n_queries: int


def qpsq_attack(
    crpuf: CRQPUF,
    n_queries: int,
    eps: float,
    tau: float,
    rng: np.random.Generator,
) -> AttackResult:
    """Adversary runs Algorithm 1 against the CR-QPUF's QPStat oracle."""
    oracle = GaussianQPStat(crpuf.process, rng=rng)
    learned = run_algorithm1(
        oracle=oracle,
        n_qubits=crpuf.n_qubits,
        observable=crpuf.observable,
        eps=eps,
        tau=tau,
        n_samples=n_queries,
        rng=rng,
    )
    return AttackResult(learned=learned, n_queries=n_queries)


def attack_success_rate(
    crpuf: CRQPUF,
    protocol: AuthenticationProtocol,
    attack: AttackResult,
    test_challenges: list[StabilizerProductState],
    rng: np.random.Generator,
) -> float:
    """Fraction of test challenges for which `h(rho)` passes verification."""
    truthful = protocol.enroll(test_challenges, rng)
    success = 0
    for challenge, enrolled in zip(test_challenges, truthful, strict=True):
        forged: float = attack.learned.model.predict(challenge.to_statevector())
        if protocol.verify(challenge, forged, enrolled):
            success += 1
    return success / len(test_challenges)
