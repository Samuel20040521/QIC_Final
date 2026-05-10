"""Algorithm 1 — average-case shadow tomography of a quantum process.

Faithful reimplementation of the pseudocode on p. 12 of Wadhwa & Doosti.

The three stages:

  1. **Gather**: for `l = 1..N`, sample a stabilizer product state
     `|psi_l>` uniformly from `stab_1^{otimes n}` and query
     `y_l = QPStat(rho_l, O, tau)`.
  2. **Learn**: for every Pauli string `P` with `|P| <= k`,
     `x_hat_P = (1/N) sum_l y_l * tr(P |psi_l><psi_l|)`.
     Apply the paper's two-part threshold:
       * `(1/3)^|P| > 2 * eps_tilde`           (degree filter)
       * `|x_hat_P| > 2 * 3^{|P|/2} * sqrt(eps_tilde) * ||O||_{Pauli,1}`
                                              (magnitude filter)
     If both hold, set `alpha_hat_P = 3^|P| * x_hat_P`. Else `alpha_hat_P = 0`.
  3. **Predict**: `h(rho) = sum_{|P| <= k} alpha_hat_P * tr(P rho)`.

Hyperparameters (Equation 42):
    `k = ceil(log_{1.5}(2 / eps^2))`
    `eps_tilde = c * eps^2 / (2n)^k`   (the paper's `Theta(.)`; we expose `c`)
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from qiskit.quantum_info import DensityMatrix, SparsePauliOp, Statevector

from .designs import (
    StabilizerProductState,
    random_stabilizer_product_state,
)
from .observables import enumerate_low_weight_paulis, pauli_l1_norm
from .oracle import QPStat


# ----- hyperparameters -------------------------------------------------------


def k_from_epsilon(eps: float) -> int:
    """Algorithm 1 sets `k = ceil(log_{1.5}(2 / eps^2))`."""
    return math.ceil(math.log(2.0 / (eps**2)) / math.log(1.5))


def epsilon_tilde(eps: float, n: int, k: int, c: float = 1.0) -> float:
    """The paper's `eps_tilde = Theta(eps^2 / (2n)^k)`. We expose the constant."""
    return c * (eps**2) / ((2 * n) ** k)


# ----- gather ----------------------------------------------------------------


@dataclass
class GatheredSample:
    state: StabilizerProductState
    y: float


def gather(
    oracle: QPStat,
    n_qubits: int,
    observable: SparsePauliOp,
    tau: float,
    n_samples: int,
    rng: np.random.Generator,
) -> list[GatheredSample]:
    samples: list[GatheredSample] = []
    for _ in range(n_samples):
        s = random_stabilizer_product_state(n_qubits, rng)
        y = oracle.query(s.to_statevector(), observable, tau)
        samples.append(GatheredSample(state=s, y=y))
    return samples


# ----- learn -----------------------------------------------------------------


@dataclass
class LearnedModel:
    """Sparse Pauli expansion of `E^*(O)` truncated to degree `<= k`."""

    n_qubits: int
    k: int
    coefficients: dict[str, float]  # label -> alpha_hat_P

    def as_sparse_pauli_op(self) -> SparsePauliOp:
        """Bake non-zero coefficients into a single `SparsePauliOp`.

        Reused across many predictions: one `DensityMatrix.expectation_value`
        call replaces a Python loop over thousands of Pauli labels.
        """
        items = [(label, alpha) for label, alpha in self.coefficients.items() if alpha != 0.0]
        if not items:
            zero_label = "I" * self.n_qubits
            return SparsePauliOp.from_list([(zero_label, 0.0)])
        return SparsePauliOp.from_list(items)

    def predict(self, target: Statevector | DensityMatrix) -> float:
        op = self.as_sparse_pauli_op()
        if isinstance(target, Statevector):
            return float(np.real(target.expectation_value(op)))
        return float(np.real(target.expectation_value(op)))


def learn(
    samples: list[GatheredSample],
    n_qubits: int,
    k: int,
    eps_tilde: float,
    observable_pauli_l1: float,
) -> LearnedModel:
    """Estimate `x_hat_P` for each Pauli with `|P| <= k`, then threshold.

    Speedup over a naive 4^n enumeration: for any stabilizer product state,
    only 2^n Pauli strings have non-zero `tr(P |psi><psi|)` (those whose
    non-identity factors match the stabilizer axes on each qubit). We
    enumerate exactly those per sample and accumulate.
    """
    n_samples = len(samples)
    if n_samples == 0:
        raise ValueError("Need at least one gathered sample")
    x_hat: dict[str, float] = {}
    inv_n = 1.0 / n_samples
    for sample in samples:
        axes = sample.state.axes
        signs = sample.state.signs
        for subset in range(1 << n_qubits):
            deg = subset.bit_count()
            if deg > k:
                continue
            label_chars = ["I"] * n_qubits
            factor = 1
            for q in range(n_qubits):
                if (subset >> q) & 1:
                    label_chars[q] = axes[q]
                    factor *= signs[q]
            # Qiskit little-endian: position 0 of label = qubit n-1.
            label = "".join(reversed(label_chars))
            x_hat[label] = x_hat.get(label, 0.0) + (factor * sample.y) * inv_n

    coefficients: dict[str, float] = {}
    for label in enumerate_low_weight_paulis(n_qubits, k):
        deg = sum(1 for c in label if c != "I")
        xv = x_hat.get(label, 0.0)
        if (1.0 / 3.0) ** deg <= 2.0 * eps_tilde:
            coefficients[label] = 0.0
            continue
        threshold = 2.0 * (3.0 ** (deg / 2.0)) * math.sqrt(eps_tilde) * observable_pauli_l1
        coefficients[label] = (3.0**deg) * xv if abs(xv) > threshold else 0.0
    return LearnedModel(n_qubits=n_qubits, k=k, coefficients=coefficients)


# ----- end-to-end convenience wrapper ---------------------------------------


@dataclass
class Algorithm1Result:
    model: LearnedModel
    n_samples: int
    k: int
    eps_tilde: float


def run_algorithm1(
    *,
    oracle: QPStat,
    n_qubits: int,
    observable: SparsePauliOp,
    eps: float,
    tau: float,
    n_samples: int,
    rng: np.random.Generator,
    eps_tilde_constant: float = 1.0,
    eps_tilde_override: float | None = None,
) -> Algorithm1Result:
    """One-call convenience wrapper around gather -> learn.

    `eps_tilde_override` lets experiments calibrate the implicit constant in
    the paper's `Theta(eps^2 / (2n)^k)`. When `None`, we use
    `eps_tilde_constant * eps^2 / (2n)^k`.
    """
    k = k_from_epsilon(eps)
    et = (
        eps_tilde_override
        if eps_tilde_override is not None
        else epsilon_tilde(eps, n_qubits, k, c=eps_tilde_constant)
    )
    samples = gather(oracle, n_qubits, observable, tau, n_samples, rng)
    model = learn(samples, n_qubits, k, et, pauli_l1_norm(observable))
    return Algorithm1Result(model=model, n_samples=n_samples, k=k, eps_tilde=et)


# ----- prediction helper ----------------------------------------------------


def predict_callable(model: LearnedModel) -> Callable[[Statevector], float]:
    return lambda state: model.predict(state)


# ----- online (drift-aware) learning ----------------------------------------


def online_learn(
    samples: list[GatheredSample],
    n_qubits: int,
    k: int,
    eps_tilde: float,
    observable_pauli_l1: float,
    half_life: float | None = None,
) -> LearnedModel:
    """Drift-aware variant of `learn`.

    Each sample at index `l` (0-based, oldest first) gets weight
    `gamma^(N - 1 - l)` with `gamma = 2^(-1 / half_life)`, so the most
    recent sample has weight 1 and the oldest has weight `gamma^(N-1)`.
    Setting `half_life=None` (or `+inf`) recovers uniform weighting,
    i.e. exactly the behavior of `learn`.

    The same threshold rule as the original Algorithm 1 is applied at
    the end, with `x_hat` interpreted as the weighted average.
    """
    n_samples = len(samples)
    if n_samples == 0:
        raise ValueError("Need at least one gathered sample")

    if half_life is None or not np.isfinite(half_life) or half_life <= 0:
        gamma = 1.0
    else:
        gamma = float(2.0 ** (-1.0 / half_life))

    # Precompute weights with the most recent sample weight = 1.
    if gamma == 1.0:
        weights = np.ones(n_samples, dtype=float)
    else:
        weights = np.array(
            [gamma ** (n_samples - 1 - l) for l in range(n_samples)], dtype=float
        )
    w_sum = float(weights.sum())
    if w_sum <= 0:
        raise ValueError("weights sum to zero")

    x_hat: dict[str, float] = {}
    for sample, w in zip(samples, weights, strict=True):
        if w == 0:
            continue
        axes = sample.state.axes
        signs = sample.state.signs
        contribution_factor = (w * sample.y) / w_sum
        for subset in range(1 << n_qubits):
            deg = subset.bit_count()
            if deg > k:
                continue
            label_chars = ["I"] * n_qubits
            factor = 1
            for q in range(n_qubits):
                if (subset >> q) & 1:
                    label_chars[q] = axes[q]
                    factor *= signs[q]
            label = "".join(reversed(label_chars))
            x_hat[label] = x_hat.get(label, 0.0) + factor * contribution_factor

    coefficients: dict[str, float] = {}
    for label in enumerate_low_weight_paulis(n_qubits, k):
        deg = sum(1 for c in label if c != "I")
        xv = x_hat.get(label, 0.0)
        if (1.0 / 3.0) ** deg <= 2.0 * eps_tilde:
            coefficients[label] = 0.0
            continue
        threshold = 2.0 * (3.0 ** (deg / 2.0)) * math.sqrt(eps_tilde) * observable_pauli_l1
        coefficients[label] = (3.0**deg) * xv if abs(xv) > threshold else 0.0
    return LearnedModel(n_qubits=n_qubits, k=k, coefficients=coefficients)
