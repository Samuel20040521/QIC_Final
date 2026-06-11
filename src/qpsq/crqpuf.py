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

Defense (Remark 4 of the paper, made concrete): the device adds a secret
structured bias `b(c) = A * prod_{i in M} sign_i(c)` to every response, for a
secret qubit subset `M` of size `w`. This deliberately violates Assumption 1
(Eq. 53) in a way that is invisible to any degree-`< w` Fourier estimate
(see `compute_structured_bias`), so a `k < w` truncated Algorithm 1 converges
to the *unbiased* low-degree predictor and its forgeries miss `b(c)` entirely.
Because `b(c)` is deterministic per challenge, it cancels between enrollment
and any honest re-query — legitimate authentication is unaffected.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from qiskit.quantum_info import SparsePauliOp, Statevector

from .algorithm import Algorithm1Result, run_algorithm1
from .designs import StabilizerProductState, random_stabilizer_product_state
from .observables import pauli_z_first
from .oracle import GaussianQPStat, QPStat, QuantumProcess, UnitaryProcess, gaussian_noise_z


# ----- the defense: secret high-weight sign-parity bias ----------------------


def compute_structured_bias(
    signs: tuple[int, ...], secret_mask: tuple[int, ...], amplitude: float
) -> float:
    """`b(x) = A * prod_{i in M} sign_i(x)` for secret mask `M`, `|M| = w`.

    Why this hides from a degree-`k < w` learner: over the uniform stabilizer
    product ensemble, the bias's Fourier coefficient against a Pauli `P` is
    `b_hat(P) = E_x[b(x) * tr(P rho_x)]`, which factorizes per qubit. For a
    qubit `i in M`, the factor is `E[sign_i * tr_i]`: zero if `P_i = I`
    (`E[sign_i] = 0`) and `1/3` if `P_i != I` (axis matches `P_i` w.p. `1/3`,
    then `sign_i^2 = 1`). For `i not in M`, the factor is `E[tr_i]`: one if
    `P_i = I`, zero otherwise (`E[sign_i] = 0`). Hence `b_hat(P)` is non-zero
    **only** for `P` non-identity on exactly `M`, where it equals
    `A * (1/3)^w` — every Pauli of weight `< w` (and any weight-`w` Pauli not
    supported on `M`) carries exactly zero bias mass, so Algorithm 1's
    `x_hat_P` estimates for `|P| <= k < w` are unchanged in expectation.
    """
    return amplitude * math.prod(signs[i] for i in secret_mask)


def structured_bias_for_challenge(
    challenge: StabilizerProductState, secret_mask: tuple[int, ...], amplitude: float
) -> float:
    """Thin wrapper over `compute_structured_bias` reading `challenge.signs`."""
    return compute_structured_bias(challenge.signs, secret_mask, amplitude)


@dataclass
class BiasedQPStat(QPStat):
    """A QPStat oracle that *deliberately violates Assumption 1* (Eq. 53):

        query(rho_x, O, tau) = tr(O E(rho_x)) + b(x) + N(0, (tau/z)^2)

    with `b(x)` the secret sign-parity bias above. The Gaussian part reuses
    the `GaussianQPStat` calibration (`sigma = tau / Phi^{-1}(1 - delta/2)`).
    The bias is keyed off the keyword-only `challenge`; when no challenge
    identity is supplied the oracle degrades to the unbiased response.
    """

    process: QuantumProcess
    secret_mask: tuple[int, ...]
    amplitude: float
    delta: float = 0.0455
    rng: np.random.Generator | None = None

    def __post_init__(self) -> None:
        self._z = gaussian_noise_z(self.delta)

    def query(
        self,
        state: Statevector,
        observable: SparsePauliOp,
        tau: float,
        *,
        challenge: object | None = None,
    ) -> float:
        true_val = self.process.expectation(state, observable)
        bias = 0.0
        if challenge is not None:
            bias = structured_bias_for_challenge(challenge, self.secret_mask, self.amplitude)
        sigma = tau / self._z
        rng = self.rng if self.rng is not None else np.random.default_rng()
        return true_val + bias + float(rng.normal(0.0, sigma))


@dataclass
class CRQPUF:
    """A device that, on challenge `rho`, produces `tr(O E(rho))` (noisy).

    When `bias_mask` is set, the device carries the sign-parity defense: every
    response — to verifier and adversary alike — includes the secret bias
    `b(c)`. The bias is a property of the *device*, not of who is asking.
    """

    process: QuantumProcess
    observable: SparsePauliOp
    n_qubits: int
    bias_mask: tuple[int, ...] | None = None
    bias_amplitude: float = 0.0

    @classmethod
    def from_haar_unitary(
        cls,
        n: int,
        rng: np.random.Generator,
        bias_mask: tuple[int, ...] | None = None,
        bias_amplitude: float = 0.0,
    ) -> "CRQPUF":
        from .designs import haar_unitary

        return cls(
            process=UnitaryProcess(haar_unitary(n, rng)),
            observable=pauli_z_first(n),
            n_qubits=n,
            bias_mask=bias_mask,
            bias_amplitude=bias_amplitude,
        )

    def make_oracle(self, rng: np.random.Generator) -> QPStat:
        """The device's QPStat interface, honoring the defense config."""
        if self.bias_mask is not None:
            return BiasedQPStat(
                process=self.process,
                secret_mask=self.bias_mask,
                amplitude=self.bias_amplitude,
                rng=rng,
            )
        return GaussianQPStat(self.process, rng=rng)

    def respond(self, challenge: StabilizerProductState, tau: float, rng: np.random.Generator) -> float:
        """Honest prover response — exact expectation plus query-tolerance noise
        (plus the secret bias when the defense is on)."""
        oracle = self.make_oracle(rng)
        return oracle.query(
            challenge.to_statevector(), self.observable, tau, challenge=challenge
        )


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
    k_override: int | None = None,
    eps_tilde_override: float | None = None,
) -> AttackResult:
    """Adversary runs Algorithm 1 against the CR-QPUF's QPStat oracle.

    The oracle honors the device's defense config: a defended CR-QPUF answers
    the adversary's training queries through the *same* biased channel it
    answers the verifier through. `k_override` pins the adversary's truncation
    degree explicitly (the degree-bounded threat model); pass a tiny
    `eps_tilde_override` (~1e-12) to keep Algorithm 1's thresholds inactive
    so the sweep over `k` is the only knob.
    """
    oracle = crpuf.make_oracle(rng)
    learned = run_algorithm1(
        oracle=oracle,
        n_qubits=crpuf.n_qubits,
        observable=crpuf.observable,
        eps=eps,
        tau=tau,
        n_samples=n_queries,
        rng=rng,
        k_override=k_override,
        eps_tilde_override=eps_tilde_override,
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
