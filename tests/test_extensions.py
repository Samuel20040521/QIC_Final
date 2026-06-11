"""Tests for the Phase-G extensions: drift, online learn, Clifford+T,
hardware-noise circuit process."""

from __future__ import annotations

import numpy as np
import pytest
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import HGate
from qiskit.quantum_info import Clifford, Operator

from qpsq.algorithm import (
    epsilon_tilde,
    gather,
    k_from_epsilon,
    learn,
    online_learn,
)
from qpsq.designs import (
    clifford_plus_t_circuit,
    haar_unitary,
    random_stabilizer_product_state,
)
from qpsq.observables import pauli_l1_norm, pauli_z_first
from qpsq.oracle import (
    CircuitNoiseProcess,
    DriftingProcess,
    GaussianQPStat,
    NoisyUnitaryProcess,
    UnitaryProcess,
    default_hardware_noise_model,
    linear_depolarizing_drift,
)


# --- Clifford + T circuit ---------------------------------------------------


def test_clifford_plus_t_is_unitary() -> None:
    rng = np.random.default_rng(0)
    for n in (2, 3, 4):
        for t_count in (0, 1, 4):
            U = clifford_plus_t_circuit(n=n, depth=2, t_count=t_count, rng=rng)
            mat = U.data
            err = np.linalg.norm(mat @ mat.conj().T - np.eye(2**n))
            assert err < 1e-9, (n, t_count, err)


def test_clifford_plus_t_zero_T_normalizes_paulis() -> None:
    """With `t_count = 0` the circuit is in the Clifford group, so for any
    Pauli `P`, `U^dag P U` is again a Pauli (up to a global sign), which
    means each row of the matrix `U^dag P U` has exactly one non-zero
    entry.
    """
    from qiskit.quantum_info import Pauli

    rng = np.random.default_rng(1)
    for trial in range(3):
        U = clifford_plus_t_circuit(n=2, depth=3, t_count=0, rng=rng)
        for label in ("IX", "IZ", "XI", "ZZ", "XY"):
            P = Operator(Pauli(label))
            evolved = U.adjoint() @ P @ U
            mat = evolved.data
            non_zero = np.abs(mat) > 1e-9
            row_counts = non_zero.sum(axis=1)
            assert (row_counts == 1).all(), (trial, label, row_counts)


def test_clifford_plus_t_with_T_breaks_pauli_normalization() -> None:
    """With `t_count > 0` we generally leave the Clifford group, in which
    case `U^dag P U` for some Pauli `P` has more than one non-zero entry
    per row. We try several seeds and several Paulis; at least one should
    show non-Pauli structure.
    """
    from qiskit.quantum_info import Pauli

    found_break = False
    for seed in range(20):
        rng = np.random.default_rng(seed)
        U = clifford_plus_t_circuit(n=2, depth=2, t_count=1, rng=rng)
        for label in ("XI", "IX", "XY", "YI"):
            evolved = U.adjoint() @ Operator(Pauli(label)) @ U
            if (np.abs(evolved.data) > 1e-9).sum(axis=1).max() > 1:
                found_break = True
                break
        if found_break:
            break
    assert found_break, "a single T-gate should break Pauli normalisation"


# --- DriftingProcess --------------------------------------------------------


def test_drifting_process_advances_t_per_query() -> None:
    rng = np.random.default_rng(0)
    base = haar_unitary(2, rng)
    factory = linear_depolarizing_drift(base, p_start=0.0, p_end=0.5, n_total_queries=10)
    drift = DriftingProcess(factory)
    assert drift.t == 0
    state = random_stabilizer_product_state(2, rng).to_statevector()
    obs = pauli_z_first(2)
    drift.expectation(state, obs)
    drift.expectation(state, obs)
    assert drift.t == 2


def test_coherent_drift_advances_t_and_changes_output() -> None:
    """`CoherentDriftProcess` increments `t` per query and yields different
    expectation values at successive `t`."""
    from qpsq.oracle import CoherentDriftProcess, random_pauli_drift_hamiltonian

    rng = np.random.default_rng(0)
    n = 3
    U0 = haar_unitary(n, rng)
    H = random_pauli_drift_hamiltonian(n, n_terms=6, rng=rng)
    proc = CoherentDriftProcess(U0, H, drift_rate=0.1)
    obs = pauli_z_first(n)
    state = random_stabilizer_product_state(n, rng).to_statevector()
    e0 = proc.expectation(state, obs)
    e1 = proc.expectation(state, obs)
    assert proc.t == 2
    assert abs(e0 - e1) > 1e-6
    proc.reset()
    assert proc.t == 0


def test_random_pauli_drift_hamiltonian_is_hermitian_and_normalized() -> None:
    from qpsq.oracle import random_pauli_drift_hamiltonian

    rng = np.random.default_rng(0)
    H = random_pauli_drift_hamiltonian(4, n_terms=10, rng=rng, spectral_norm=2.0)
    # Hermitian
    assert np.allclose(H, H.conj().T, atol=1e-9)
    # Spectral norm matches request
    s = float(np.linalg.norm(H, ord=2))
    assert abs(s - 2.0) < 1e-9


def test_drifting_process_outputs_change_with_t() -> None:
    rng = np.random.default_rng(0)
    n = 2
    base = haar_unitary(n, rng)
    factory = linear_depolarizing_drift(base, p_start=0.0, p_end=0.6, n_total_queries=2)
    drift = DriftingProcess(factory)
    state = random_stabilizer_product_state(n, rng).to_statevector()
    obs = pauli_z_first(n)
    e0 = drift.expectation(state, obs)
    e1 = drift.expectation(state, obs)
    assert abs(e0 - e1) > 1e-6, "depolarizing rate should differ between queries"


# --- online_learn -----------------------------------------------------------


def test_online_learn_matches_learn_when_half_life_infinite() -> None:
    """With infinite half-life, online_learn should exactly recover learn()."""
    rng = np.random.default_rng(0)
    n = 3
    process = UnitaryProcess(Operator(np.eye(2**n)))
    oracle = GaussianQPStat(process, rng=np.random.default_rng(7))
    obs = pauli_z_first(n)
    samples = gather(oracle, n, obs, tau=0.05, n_samples=300, rng=rng)
    k = k_from_epsilon(0.3)
    et = epsilon_tilde(0.3, n, k)

    a = learn(samples, n, k, et, pauli_l1_norm(obs))
    b = online_learn(samples, n, k, et, pauli_l1_norm(obs), half_life=None)
    assert set(a.coefficients) == set(b.coefficients)
    for label in a.coefficients:
        assert abs(a.coefficients[label] - b.coefficients[label]) < 1e-9


def test_online_learn_short_half_life_emphasizes_recent() -> None:
    """With small half-life, predictions follow the recent process,
    not the average of an old + new mix.
    """
    rng = np.random.default_rng(0)
    n = 2
    base = haar_unitary(n, rng)
    # First half of samples come from the "old" process, second half "new".
    old = UnitaryProcess(base)
    X = np.array([[0, 1], [1, 0]], dtype=complex)
    new_unitary = Operator(np.kron(np.eye(2), X)) @ base  # post-X on qubit 0
    new = UnitaryProcess(new_unitary)
    obs = pauli_z_first(n)

    rng_state = np.random.default_rng(11)
    rng_oracle_old = np.random.default_rng(101)
    rng_oracle_new = np.random.default_rng(202)
    oracle_old = GaussianQPStat(old, rng=rng_oracle_old)
    oracle_new = GaussianQPStat(new, rng=rng_oracle_new)
    samples = []
    for _ in range(400):
        s = random_stabilizer_product_state(n, rng_state)
        y = oracle_old.query(s.to_statevector(), obs, 0.02)
        from qpsq.algorithm import GatheredSample
        samples.append(GatheredSample(state=s, y=y))
    for _ in range(400):
        s = random_stabilizer_product_state(n, rng_state)
        y = oracle_new.query(s.to_statevector(), obs, 0.02)
        from qpsq.algorithm import GatheredSample
        samples.append(GatheredSample(state=s, y=y))

    k = k_from_epsilon(0.3)
    et = epsilon_tilde(0.3, n, k)
    obs_l1 = pauli_l1_norm(obs)
    flat = online_learn(samples, n, k, et, obs_l1, half_life=None)
    drift_aware = online_learn(samples, n, k, et, obs_l1, half_life=50.0)

    # Test on stabilizer states under the *new* process: drift-aware should
    # have lower MAE than uniform.
    test_states = [random_stabilizer_product_state(n, rng_state).to_statevector()
                   for _ in range(60)]
    truths = np.array([new.expectation(s, obs) for s in test_states])
    flat_pred = np.array([flat.predict(s) for s in test_states])
    drift_pred = np.array([drift_aware.predict(s) for s in test_states])
    flat_mae = float(np.mean(np.abs(flat_pred - truths)))
    drift_mae = float(np.mean(np.abs(drift_pred - truths)))
    assert drift_mae < flat_mae, (drift_mae, flat_mae)


# --- CircuitNoiseProcess (smoke) -------------------------------------------


def test_circuit_noise_process_no_noise_matches_unitary() -> None:
    n = 2
    qc = QuantumCircuit(n)
    qc.h(0)
    qc.cx(0, 1)
    proc_noisy = CircuitNoiseProcess(circuit=qc, noise_model=None)
    proc_unitary = UnitaryProcess(Operator(qc))
    rng = np.random.default_rng(0)
    state = random_stabilizer_product_state(n, rng).to_statevector()
    obs = pauli_z_first(n)
    a = proc_noisy.expectation(state, obs)
    b = proc_unitary.expectation(state, obs)
    assert abs(a - b) < 1e-8


def test_circuit_noise_process_mps_matches_density_matrix() -> None:
    """For a moderately entangling unitary, MPS and density-matrix backends
    agree on observable expectation values within MPS truncation tolerance."""
    from qiskit.quantum_info import Operator, SparsePauliOp, Statevector

    n = 4
    qc = QuantumCircuit(n)
    qc.h(0); qc.cx(0, 1); qc.cx(1, 2); qc.cx(2, 3); qc.rz(0.7, 0); qc.ry(0.3, 2)
    proc_dm = CircuitNoiseProcess(circuit=qc, noise_model=None, method="density_matrix")
    proc_mps = CircuitNoiseProcess(circuit=qc, noise_model=None, method="matrix_product_state")
    proc_auto_small = CircuitNoiseProcess(circuit=qc, noise_model=None, method="auto")
    assert proc_auto_small._resolved_method == "density_matrix"
    psi = Statevector.from_label("0" * n)
    obs = SparsePauliOp.from_list([("ZZZZ", 1.0), ("IIIX", 0.5)])
    a = proc_dm.expectation(psi, obs)
    b = proc_mps.expectation(psi, obs)
    assert abs(a - b) < 1e-6, (a, b)


def test_circuit_noise_process_auto_method_for_large_n() -> None:
    """`method="auto"` picks MPS for n > 10."""
    qc = QuantumCircuit(11)
    qc.h(0)
    proc = CircuitNoiseProcess(circuit=qc, noise_model=None, method="auto")
    assert proc._resolved_method == "matrix_product_state"


def test_circuit_noise_process_with_default_noise_changes_output() -> None:
    """H + CX prepares a Bell state from |00>; tr(ZZ |Bell><Bell|) = 1.
    A noisy version should give strictly less."""
    from qiskit.quantum_info import SparsePauliOp, Statevector

    n = 2
    qc = QuantumCircuit(n)
    qc.h(0)
    qc.cx(0, 1)
    nm = default_hardware_noise_model(p1=0.05, p2=0.1, readout_p=0.0)
    proc_noisy = CircuitNoiseProcess(circuit=qc, noise_model=nm)
    proc_clean = UnitaryProcess(Operator(qc))
    state = Statevector.from_label("0" * n)
    obs = SparsePauliOp.from_list([("ZZ", 1.0)])
    a = proc_noisy.expectation(state, obs)
    b = proc_clean.expectation(state, obs)
    assert b > 0.99 and a < b - 1e-2, (a, b)


# --- CR-QPUF sign-parity defense (Remark 4 instantiation) --------------------


def test_structured_bias_low_degree_fourier_mass_is_exactly_zero() -> None:
    """Over the *full* stabilizer product ensemble (6^n states), the bias
    `b(x) = A * prod_{i in M} sign_i(x)` has exactly zero Fourier mass on
    every Pauli of weight < w, and mass `A / 3^w` on every Pauli that is
    non-identity exactly on `M`."""
    from itertools import product

    from qpsq.crqpuf import compute_structured_bias
    from qpsq.designs import (
        StabilizerProductState,
        pauli_expectation_on_stabilizer,
    )

    n = 4
    mask = (0, 2, 3)  # w = 3
    w = len(mask)
    amplitude = 0.7

    ensemble = [
        StabilizerProductState(axes=axes, signs=signs)
        for axes in product("XYZ", repeat=n)
        for signs in product((1, -1), repeat=n)
    ]
    assert len(ensemble) == 6**n

    def fourier_mass(label: str) -> float:
        total = 0.0
        for state in ensemble:
            b = compute_structured_bias(state.signs, mask, amplitude)
            total += b * pauli_expectation_on_stabilizer(label, state)
        return total / len(ensemble)

    # Every Pauli of weight <= w - 1: exactly zero mass.
    for chars in product("IXYZ", repeat=n):
        weight = sum(c != "I" for c in chars)
        if weight > w - 1:
            continue
        label = "".join(reversed(chars))  # chars indexed by qubit -> Qiskit label
        assert abs(fourier_mass(label)) < 1e-12, (label, fourier_mass(label))

    # Weight-w Paulis supported exactly on M: mass A / 3^w, for any axes.
    for on_mask_axes in (("X", "X", "X"), ("Z", "Y", "X")):
        chars = ["I"] * n
        for q, ax in zip(mask, on_mask_axes, strict=True):
            chars[q] = ax
        label = "".join(reversed(chars))
        assert abs(fourier_mass(label) - amplitude / 3**w) < 1e-12, label


def test_defended_attack_fails_below_parity_weight() -> None:
    """At truncation degree k = 2 < w = 3, the defended device's responses
    are unforgeable: attack success collapses versus the undefended device."""
    from qpsq.crqpuf import (
        AuthenticationProtocol,
        CRQPUF,
        attack_success_rate,
        qpsq_attack,
        random_challenge_set,
    )
    from qpsq.observables import pauli_z_first

    rng = np.random.default_rng(7)
    n = 4
    mask = (0, 1, 3)  # w = 3
    amplitude = 0.4
    unitary = haar_unitary(n, rng)
    challenges = random_challenge_set(n, 60, rng)

    rates = {}
    for defended in (False, True):
        crpuf = CRQPUF(
            process=UnitaryProcess(unitary),
            observable=pauli_z_first(n),
            n_qubits=n,
            bias_mask=mask if defended else None,
            bias_amplitude=amplitude if defended else 0.0,
        )
        protocol = AuthenticationProtocol(
            crpuf=crpuf, accept_threshold=0.1, response_tau=0.02
        )
        attack = qpsq_attack(
            crpuf=crpuf,
            n_queries=4000,
            eps=0.3,
            tau=0.05,
            rng=rng,
            k_override=2,
            eps_tilde_override=1e-12,
        )
        rates[defended] = attack_success_rate(crpuf, protocol, attack, challenges, rng)

    assert rates[True] < rates[False] - 0.2, rates


def test_defended_device_keeps_honest_completeness() -> None:
    """The bias is deterministic per challenge, so it cancels between
    enrollment and an honest re-query: completeness stays ~1."""
    from qpsq.crqpuf import AuthenticationProtocol, CRQPUF, random_challenge_set

    rng = np.random.default_rng(11)
    n = 4
    crpuf = CRQPUF.from_haar_unitary(n, rng, bias_mask=(0, 1, 2), bias_amplitude=0.4)
    protocol = AuthenticationProtocol(crpuf=crpuf, accept_threshold=0.1, response_tau=0.02)
    challenges = random_challenge_set(n, 50, rng)
    enrolled = protocol.enroll(challenges, rng)
    passes = sum(
        protocol.verify(c, crpuf.respond(c, 0.02, rng), e)
        for c, e in zip(challenges, enrolled, strict=True)
    )
    assert passes / len(challenges) >= 0.95
