"""Extension: end-to-end demo of the QPSQ-based attack on a CR-QPUF authentication.

We follow the protocol sketch in Section 6 of Wadhwa & Doosti:

  * Server holds a fixed Haar-random unitary `U` (the CR-QPUF). Authentication
    queries are stabilizer product states `rho`. The honest device returns
    `tr(O E(rho))` (with measurement noise of size `tau`) for `O = Z_0`.
  * Verifier accepts a claimed response if it is within `accept_threshold` of
    the enrolled response.
  * Adversary uses Algorithm 1 against the device's QPStat oracle to learn a
    predictor `h(rho)`. With `N` queries it gets a quasipolynomial-time
    approximation of `tr(O E(.))`. We measure attack success rate (fraction
    of test challenges accepted) as a function of `N`.

This is a *toy* demonstration: we do not formally break CR-QPUF security
(quasipolynomial complexity), but we show that for small n the success rate
climbs steeply with N, exactly as the paper's analysis predicts.

Outputs:
  experiments/results/ext_crqpuf_attack.csv
  report/figs/ext_crqpuf_attack.pdf
"""

from __future__ import annotations

import csv
import time

import matplotlib.pyplot as plt
import numpy as np

from qpsq.crqpuf import (
    AuthenticationProtocol,
    CRQPUF,
    attack_success_rate,
    qpsq_attack,
    random_challenge_set,
)

from _common import FIGS_DIR, RESULTS_DIR, make_rng, parse_seed_only

N_QUBITS = 4
N_QPUF_INSTANCES = 6
N_TEST_CHALLENGES = 80
ACCEPT_THRESHOLD = 0.10
RESPONSE_TAU = 0.02
EPS_FOR_ATTACK = 0.3
ATTACK_TAU = 0.05
NS = (200, 500, 1000, 2000, 4000, 8000)


def main() -> None:
    args = parse_seed_only()
    rng = make_rng(args.seed)

    rows: list[tuple] = []
    t_start = time.time()
    for inst in range(N_QPUF_INSTANCES):
        crpuf = CRQPUF.from_haar_unitary(N_QUBITS, rng)
        protocol = AuthenticationProtocol(
            crpuf=crpuf,
            accept_threshold=ACCEPT_THRESHOLD,
            response_tau=RESPONSE_TAU,
        )
        test_challenges = random_challenge_set(N_QUBITS, N_TEST_CHALLENGES, rng)

        # Honest baseline: how often does an *honest* prover pass verification?
        # Bounded by RESPONSE_TAU << ACCEPT_THRESHOLD, so should be ~100%.
        honest_pass = 0
        for c in test_challenges:
            r1 = crpuf.respond(c, RESPONSE_TAU, rng)
            r2 = crpuf.respond(c, RESPONSE_TAU, rng)
            if abs(r1 - r2) <= ACCEPT_THRESHOLD:
                honest_pass += 1
        honest_rate = honest_pass / N_TEST_CHALLENGES

        # Random-guess baseline: the adversary outputs a random number in [-1, 1].
        random_pass = 0
        truthful = protocol.enroll(test_challenges, rng)
        for r in truthful:
            guess = float(rng.uniform(-1, 1))
            if abs(guess - r) <= ACCEPT_THRESHOLD:
                random_pass += 1
        random_rate = random_pass / N_TEST_CHALLENGES

        rows.append((inst, "honest", 0, honest_rate))
        rows.append((inst, "random_guess", 0, random_rate))

        for n_query in NS:
            attack = qpsq_attack(
                crpuf=crpuf,
                n_queries=n_query,
                eps=EPS_FOR_ATTACK,
                tau=ATTACK_TAU,
                rng=rng,
            )
            rate = attack_success_rate(crpuf, protocol, attack, test_challenges, rng)
            rows.append((inst, "qpsq_attack", n_query, rate))
        print(f"  qpuf {inst + 1}/{N_QPUF_INSTANCES} done at "
              f"t={time.time() - t_start:.1f}s")

    csv_path = RESULTS_DIR / "ext_crqpuf_attack.csv"
    with csv_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["qpuf_instance", "strategy", "n_queries", "success_rate"])
        for r in rows:
            w.writerow([r[0], r[1], r[2], f"{r[3]:.4f}"])

    fig, ax = plt.subplots(figsize=(7, 4.5))
    qpsq_means = []
    qpsq_stds = []
    for n_query in NS:
        vals = [r[3] for r in rows if r[1] == "qpsq_attack" and r[2] == n_query]
        qpsq_means.append(float(np.mean(vals)))
        qpsq_stds.append(float(np.std(vals)))
    ax.errorbar(
        NS, qpsq_means, yerr=qpsq_stds, marker="o",
        color="C3", capsize=4, label="QPSQ attack (Algorithm 1)",
    )
    honest = float(np.mean([r[3] for r in rows if r[1] == "honest"]))
    rand = float(np.mean([r[3] for r in rows if r[1] == "random_guess"]))
    ax.axhline(honest, color="C2", linestyle="--", label=f"honest prover ({honest:.2f})")
    ax.axhline(rand, color="grey", linestyle=":", label=f"random guess ({rand:.2f})")
    ax.set_xscale("log")
    ax.set_xlabel("N (adversary's QPStat queries)")
    ax.set_ylabel("verification acceptance rate")
    ax.set_title(
        f"CR-QPUF attack via Algorithm 1 (n={N_QUBITS}, threshold={ACCEPT_THRESHOLD}, "
        f"{N_QPUF_INSTANCES} devices)"
    )
    ax.set_ylim(-0.02, 1.05)
    ax.grid(True, which="both", linestyle=":", alpha=0.5)
    ax.legend(loc="lower right")
    fig.tight_layout()
    out = FIGS_DIR / "ext_crqpuf_attack.pdf"
    fig.savefig(out)
    print(f"wrote {csv_path} and {out}")


if __name__ == "__main__":
    main()
