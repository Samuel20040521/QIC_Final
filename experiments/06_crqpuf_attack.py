"""Extension: QPSQ attack on a CR-QPUF, with and without a sign-parity defense.

Undefended baseline (Section 6 of Wadhwa & Doosti): the adversary runs
Algorithm 1 against the device's QPStat oracle and forges responses with the
learned predictor `h(rho)`. Theorem 5's guarantee rests on Assumption 1
(Eq. 53) — the oracle is *unbiased*. Remark 4 observes that biased noise
breaks the attack and "the resulting protocol may be secure", and leaves the
construction open. We instantiate it:

  Defense: the device adds a secret bias `b(c) = A * prod_{i in M} sign_i(c)`
  to every response, for a secret qubit subset `M` with `|M| = w`. The bias's
  Fourier mass sits entirely on weight-`w` Paulis supported on `M`
  (`qpsq.crqpuf.compute_structured_bias`), so any degree-`k < w` truncated
  Algorithm 1 converges to the *unbiased* low-degree predictor — its forgery
  is missing `b(c)` and fails verification by `~A`. Since `b(c)` is
  deterministic per challenge, it cancels between enrollment and honest
  re-queries: completeness is untouched.

We sweep the adversary's truncation degree `k in {2, 3, 4}` at `n = 6` with
`w = 4` to exhibit the boundary directly: the defense holds for `k < w` and
collapses at `k >= w` (the adversary's reach extends to the hiding Paulis).

Outputs:
  experiments/results/ext_crqpuf_attack.csv
  report/figs/ext_crqpuf_attack.pdf
"""

from __future__ import annotations

import csv
import time

import matplotlib.pyplot as plt
import numpy as np

from qpsq.algorithm import k_from_epsilon
from qpsq.crqpuf import (
    AuthenticationProtocol,
    CRQPUF,
    attack_success_rate,
    qpsq_attack,
    random_challenge_set,
    structured_bias_for_challenge,
)
from qpsq.designs import haar_unitary, pauli_expectation_on_stabilizer
from qpsq.observables import pauli_z_first
from qpsq.oracle import UnitaryProcess

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

N_QUBITS = 6
SECRET_MASK = (0, 2, 3, 5)  # |M| = w = 4
BIAS_AMPLITUDE = 0.3        # A > ACCEPT_THRESHOLD: a missing b(c) always rejects
N_QPUF_INSTANCES = 4
N_TEST_CHALLENGES = 80
ACCEPT_THRESHOLD = 0.10
RESPONSE_TAU = 0.02
EPS_FOR_ATTACK = 0.3        # bookkeeping only: k is pinned by K_SWEEP below
ATTACK_TAU = 0.05
K_SWEEP = (2, 3, 4)
NS = (200, 500, 1000, 2000, 4000, 8000, 16000)
# Keep Algorithm 1's two-part thresholds inactive so the truncation degree k
# is the only knob (same trick as the fig2 calibration).
EPS_TILDE_OVERRIDE = 1e-12

W = len(SECRET_MASK)
assert W > min(K_SWEEP), "defense is vacuous unless some adversary has k < w"
assert max(K_SWEEP) >= W, "sweep must include k >= w to exhibit the collapse"


def empirical_bias_fourier_mass(rng: np.random.Generator, n_mc: int = 30_000) -> None:
    """Monte-Carlo check of the hiding property: `E_x[b(x) tr(P rho_x)]` is
    ~0 for every `|P| < w`, and `A * (1/3)^w` for `P` non-identity exactly
    on `M`. Logged as a sanity line; the exact (full-ensemble) version is in
    tests/test_extensions.py."""
    challenges = random_challenge_set(N_QUBITS, n_mc, rng)
    biases = np.array(
        [structured_bias_for_challenge(c, SECRET_MASK, BIAS_AMPLITUDE) for c in challenges]
    )

    def mass(label: str) -> float:
        tr = np.array([pauli_expectation_on_stabilizer(label, c) for c in challenges])
        return float(np.mean(biases * tr))

    # Random low-weight Paulis (|P| < w), 25 per weight.
    worst_low = 0.0
    axes = "XYZ"
    for weight in range(1, W):
        for _ in range(25):
            qubits = rng.choice(N_QUBITS, size=weight, replace=False)
            chars = ["I"] * N_QUBITS
            for q in qubits:
                chars[q] = axes[rng.integers(0, 3)]
            worst_low = max(worst_low, abs(mass("".join(reversed(chars)))))

    # One weight-w Pauli supported exactly on M.
    chars = ["I"] * N_QUBITS
    for q in SECRET_MASK:
        chars[q] = "X"
    on_mask = mass("".join(reversed(chars)))

    mc_err = BIAS_AMPLITUDE / np.sqrt(n_mc)
    print(
        f"  [sanity] bias Fourier mass ({n_mc} MC samples): "
        f"max |b_hat(P)| over |P|<{W} = {worst_low:.5f} (MC err ~{mc_err:.5f}); "
        f"b_hat(P on M) = {on_mask:.5f} (theory A/3^w = {BIAS_AMPLITUDE / 3**W:.5f})"
    )


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)

    print(
        f"[threat model] n={N_QUBITS}, secret mask M={SECRET_MASK} (w={W}), "
        f"A={BIAS_AMPLITUDE}, accept_threshold={ACCEPT_THRESHOLD}, "
        f"k sweep={K_SWEEP} (paper's Eq. 42 would give "
        f"k={k_from_epsilon(EPS_FOR_ATTACK)} at eps={EPS_FOR_ATTACK}), "
        f"eps_tilde_override={EPS_TILDE_OVERRIDE} (thresholds inactive)"
    )
    empirical_bias_fourier_mass(rng)

    rows: list[tuple] = []
    t_start = time.time()
    for inst in range(N_QPUF_INSTANCES):
        # Same Haar unitary for the defended and undefended device, so the
        # defense toggle is the only difference within an instance.
        unitary = haar_unitary(N_QUBITS, rng)
        test_challenges = random_challenge_set(N_QUBITS, N_TEST_CHALLENGES, rng)

        for defense in ("off", "on"):
            crpuf = CRQPUF(
                process=UnitaryProcess(unitary),
                observable=pauli_z_first(N_QUBITS),
                n_qubits=N_QUBITS,
                bias_mask=SECRET_MASK if defense == "on" else None,
                bias_amplitude=BIAS_AMPLITUDE if defense == "on" else 0.0,
            )
            protocol = AuthenticationProtocol(
                crpuf=crpuf,
                accept_threshold=ACCEPT_THRESHOLD,
                response_tau=RESPONSE_TAU,
            )

            # Honest completeness: a fresh device response vs the enrolled one.
            # The deterministic bias b(c) appears in both and cancels.
            enrolled = protocol.enroll(test_challenges, rng)
            honest_pass = sum(
                protocol.verify(c, crpuf.respond(c, RESPONSE_TAU, rng), e)
                for c, e in zip(test_challenges, enrolled, strict=True)
            )
            honest_rate = honest_pass / N_TEST_CHALLENGES

            # Random-guess baseline: adversary outputs a uniform value in [-1, 1].
            random_pass = sum(
                abs(float(rng.uniform(-1, 1)) - e) <= ACCEPT_THRESHOLD for e in enrolled
            )
            random_rate = random_pass / N_TEST_CHALLENGES

            for k in K_SWEEP:
                for n_query in NS:
                    attack = qpsq_attack(
                        crpuf=crpuf,
                        n_queries=n_query,
                        eps=EPS_FOR_ATTACK,
                        tau=ATTACK_TAU,
                        rng=rng,
                        k_override=k,
                        eps_tilde_override=EPS_TILDE_OVERRIDE,
                    )
                    rate = attack_success_rate(
                        crpuf, protocol, attack, test_challenges, rng
                    )
                    rows.append(
                        (inst, defense, k, n_query, rate, honest_rate, random_rate)
                    )
            print(
                f"  qpuf {inst + 1}/{N_QPUF_INSTANCES} defense={defense}: "
                f"honest={honest_rate:.2f} random={random_rate:.2f} "
                f"t={time.time() - t_start:.1f}s"
            )

    csv_path = RESULTS_DIR / "ext_crqpuf_attack.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            ["qpuf_instance", "defense", "k", "n_queries",
             "attack_success", "honest_success", "random_success"]
        )
        for r in rows:
            w.writerow([r[0], r[1], r[2], r[3], f"{r[4]:.4f}", f"{r[5]:.4f}", f"{r[6]:.4f}"])

    # ----- figure: two panels, undefended vs defended -----
    fig, axes_panels = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, defense, title in zip(
        axes_panels,
        ("off", "on"),
        (
            f"defense off (n={N_QUBITS})",
            f"defense on: secret parity w={W}, A={BIAS_AMPLITUDE}",
        ),
        strict=True,
    ):
        sub = [r for r in rows if r[1] == defense]
        for ki, k in enumerate(K_SWEEP):
            means, stds = [], []
            for n_query in NS:
                vals = [r[4] for r in sub if r[2] == k and r[3] == n_query]
                means.append(float(np.mean(vals)))
                stds.append(float(np.std(vals)))
            style = "-o" if k < W else "--s"
            label = f"k={k}" + (" (= w)" if k == W else "")
            ax.errorbar(NS, means, yerr=stds, fmt=style, color=["C0", "C1", "C3"][ki],
                        capsize=3, label=label)
        honest = float(np.mean([r[5] for r in sub]))
        rand = float(np.mean([r[6] for r in sub]))
        ax.axhline(honest, color="black", linestyle="--", alpha=0.7,
                   label=f"honest prover ({honest:.2f})")
        ax.axhline(rand, color="grey", linestyle=":",
                   label=f"random guess ({rand:.2f})")
        ax.set_xscale("log")
        ax.set_xlabel("N (adversary's QPStat queries)")
        ax.set_title(title)
        ax.set_ylim(-0.02, 1.05)
        ax.grid(True, which="both", linestyle=":", alpha=0.5)
        ax.legend(loc="center left", fontsize=8)
    axes_panels[0].set_ylabel("verification acceptance rate")
    fig.suptitle(
        f"CR-QPUF attack vs sign-parity defense "
        f"(threshold={ACCEPT_THRESHOLD}, {N_QPUF_INSTANCES} devices)"
    )
    fig.tight_layout()
    out = FIGS_DIR / "ext_crqpuf_attack.pdf"
    fig.savefig(out)
    print(f"wrote {csv_path} and {out}")


if __name__ == "__main__":
    main()
