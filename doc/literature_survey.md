# Literature Survey — Learning Quantum Processes with Quantum Statistical Queries

Focused survey (2021–2026) of work closely related to the project's focal paper:

> **Wadhwa & Doosti**, *Learning Quantum Processes with Quantum Statistical
> Queries*, **Quantum 9, 1739 (2025)**, DOI `10.22331/q-2025-05-12-1739`,
> arXiv:2310.02075. (University of Edinburgh.)

The paper initiates the study of learning quantum **processes** (channels /
unitaries) under a **Quantum Process Statistical Query (QPSQ)** access model — a
single-copy, no-entangled-ancilla oracle `QPStat_E(ρ, O, τ)` returning a
`τ`-accurate estimate of `tr(O·E(ρ))`. Its four pillars: (i) an efficient
average-case shadow-process-tomography algorithm (Algorithm 1) with a nearly
matching `Ω(M)` observable-count lower bound; (ii) an exponential QPSQ lower
bound for unitary 2-designs and (iii) a **doubly** exponential lower bound for
Haar-random unitaries in diamond distance (via a many-vs-one distinguishing
reduction); (iv) a cryptanalytic attack on Classical-Readout Quantum PUF
(CR-QPUF) authentication.

**Verification key:** ✓✓ = adversarially verified (3-0 unanimous) against a
primary source during the deep-research pass; ✓ = primary source fetched and
author/title/venue confident, not in the final adversarial set. Every entry
below is a primary arXiv/journal source.

---

## Bibliographic corrections (surfaced by the survey)

1. **Focal paper year/venue:** it is **Quantum 9, 1739 (2025)**, DOI
   `10.22331/q-2025-05-12-1739` — *not* 2024. (`refs.bib` key `wadhwa2024qpsq`
   corrected in place; the citation key string is retained to avoid breaking
   `\cite` references in `main.tex`.)
2. **arXiv:2203.03591 authorship:** the authors are **Armando Angrisani &
   Elham Kashefi** (IEEE Trans. Inf. Theory 71(5), 2025) — **Mina Doosti is
   not** an author of this paper.
3. **Chen–Cotler–Huang–Li Pauli memory bound** must be written `2^{(n-k)/3}`
   (division *inside* the exponent), not `(2^{n-k})/3`.
4. **Two distinct Nietner papers** — do not conflate: arXiv:2305.05765 (output
   distributions, existing key `nietner2023qsqshadow`) vs. arXiv:2310.17716
   (unifying SQ/QSQ/parametrized, new key `nietner2023unifying`).

---

## Category 1 — Restricted / statistical access models for learning processes

State-QSQ → process-QSQ transitions; single-copy, no-ancilla, noisy-oracle,
limited-input-control regimes.

| bibkey | paper | arXiv | venue | relevance |
| --- | --- | --- | --- | --- |
| `arunachalam2020qsq` ✓✓ | Arunachalam, Grilo, Yuen — *Quantum Statistical Query Learning* | 2002.08240 | — | original **state-level QSQ** model; the lineage W&D lift to processes |
| `angrisani2023learningunitaries` ✓✓ | Angrisani — *Learning Unitaries with Quantum Statistical Queries* | 2310.02254 | Quantum (2025) | **sibling work**: QSQ on Choi states; QSQ Goldreich–Levin (no `U⁻¹`); exp LB for phase-oracle unitaries, double-exp LB for unitarity testing |
| `wadhwa2024noisetolerant` ✓✓ | Wadhwa, Doosti — *Noise-tolerant learnability of shallow quantum circuits … cost of quantum pseudorandomness* | 2405.12085 | — | W&D follow-up: QSQ noise-tolerance; shallow-circuit learning at linear overhead; avg-case diamond LB; constant-depth ⇏ PRU |
| `angrisani2022qldp` ✓✓ | Angrisani, Kashefi — *Quantum Local Differential Privacy and the QSQ Model* | 2203.03591 | IEEE TIT 71(5) (2025) | QSQ ≡ quantum local DP; frames QSQ as limited-resource / private access |
| `nietner2023unifying` ✓✓ | Nietner — *Unifying (Quantum) Statistical and Parametrized (Quantum) Algorithms* | 2310.17716 | — | "evaluation oracle" unifying classical SQ / QSQ / VQA; unconditional LBs |

**Search strings**
1. `"quantum statistical query" (process OR channel OR unitary) learning lower bound`
2. `Choi state "statistical quer*" learning unitary Goldreich-Levin`
3. `quantum statistical query "local differential privacy" OR "noise-tolerant" shallow circuit`
4. `"evaluation oracle" quantum statistical query parametrized lower bound Nietner`

**Synthesis.** This cluster shares the QPSQ weak-access premise — only noisy
bounded-`τ` expectation estimates, single copy, no entangled ancilla. The
objects differ: Arunachalam–Grilo–Yuen learn states; Angrisani (2310.02254)
learns unitaries via Choi states (a parallel, not identical, model); W&D learn
general channels. The noise-tolerance follow-up (2405.12085) extends the same
machinery to robustness and ties learnability to quantum-pseudorandomness lower
bounds — the bridge between §5 hardness and the §6 CR-QPUF application.

---

## Category 2 — Shadow process tomography / classical shadows for predicting process properties

Query/sample-complexity optimization, average-case guarantees, non-local
observables; the Huang–Chen–Preskill line.

| bibkey | paper | arXiv | venue | relevance |
| --- | --- | --- | --- | --- |
| `huang2023predictprocess` ✓ | Huang, Chen, Preskill — *Learning to Predict Arbitrary Quantum Processes* | 2210.14894 | PRX Quantum 4, 040337 (2023) | **ref [40]** — the average-case algorithm Algorithm 1 directly adapts; the category-2 anchor |
| `raza2024online` ✓✓ | Raza, Caro, Eisert, Khatri — *Online Learning of Quantum Processes* | 2406.04250 | — | online/regret + mistake-bound learning of bounded-gate & Pauli channels; sample-efficient shadow tomography for Pauli channels (online contrast to QPSQ) |
| `chen2022singlecopy` ✓✓ | Chen, Cotler, Huang, Li — *Exponential separations between learning with and without quantum memory* | 2111.05881 | FOCS | memory-vs-sample separation; `Ω(2^{(n-k)/3})` Pauli-estimation bound; shadow tomography needs `Ω(min(M,2ⁿ))` without memory |
| `huang2022shadowprocess` ✓ | Huang, Kueng, Torlai, Albert, Preskill — *Provably efficient ML for quantum many-body problems* | — | Science 377 (2022) | classical-shadow property-prediction framework feeding the process setting |
| `kunjummen2023shadowprocess` ✓ | Kunjummen, Tran, Carney, Taylor — *Shadow process tomography of quantum channels* | 2110.03629 | PRA 107 (2023) | Choi/ancilla-assisted shadow process tomography — access-model **contrast** to no-ancilla QPSQ |
| `levy2023classical` ✓ | Levy, Luo, Clark — *Classical shadows for quantum process tomography (ShadowQPT)* | 2110.02965 | PRR 6 (2024) | near-term Choi-matrix shadow reconstruction; baseline vs. restricted-access QPSQ |
| `wadhwa2024agnostic` ✓ | Wadhwa, Lewis, Kashefi, Doosti — *Agnostic Process Tomography* | 2410.11957 | PRX Quantum (2025) | W&D-team extension of channel learning to the agnostic setting |

**Search strings**
1. `classical shadows quantum process OR channel tomography Choi prediction observable`
2. `"shadow tomography" quantum process sample complexity lower bound Huang Chen`
3. `online learning quantum channels Pauli "shadow tomography" Caro regret`
4. `learning to predict properties quantum process average-case local observable`

**Synthesis.** W&D's shadow-process-tomography result is an adaptation of
Huang–Chen–Preskill (2210.14894) with the access *downgraded* from running the
channel to bounded-`τ` statistical queries — paying a linear overhead in the
number of observables (with a near-matching lower bound). The Choi/ancilla
methods (2110.03629, 2110.02965) and the online framework (2406.04250) chart the
neighboring access models that make the near-term, single-copy, no-ancilla QPSQ
trade-offs legible.

---

## Category 3 — Quantum hardware security & cryptanalysis (CR-QPUFs)

ML attacks / unforgeability of (Classical-Readout) QPUFs; learning-theory lower
bounds as cryptographic security proofs.

| bibkey | paper | arXiv | venue | relevance |
| --- | --- | --- | --- | --- |
| `pirnay2022crqpuf` ✓✓ | Pirnay, Pappa, Seifert — *Learning Classical Readout Quantum PUFs based on single-qubit gates* | 2112.06661 | Quantum Mach. Intell. 4:14 (2022) | **SQ-model precursor** to W&D's CR-QPUF attack; low-degree-polynomial modelling attack, validated on IBM Q |
| `arapinis2021qpuf` ✓ | Arapinis, Delavar, Doosti, Kashefi — *Quantum Physical Unclonable Functions: Possibilities and Impossibilities* | 1910.02126 | Quantum 5 (2021) | foundational QPUF unforgeability definitions used to scope §6 security games |
| `doosti2021crqpuf` ✓ | Doosti, Kumar, Delavar, Kashefi — *Client-Server Identification Protocols with Quantum PUF* | — | ACM TQC (2021) | QPUF-based authentication protocol context |

**Search strings**
1. `Classical Readout Quantum PUF statistical query attack unforgeability single-qubit`
2. `quantum physical unclonable function unforgeability learning model selective`
3. `machine learning modelling attack quantum PUF authentication forge signature`
4. `learning theory lower bound cryptographic security proof quantum pseudorandom`

**Synthesis.** Pirnay–Pappa–Seifert (2112.06661) attack CR-QPUFs in the *state*
SQ model; W&D lift this to QPSQ by treating the PUF as an unknown **channel** to
be learned, yielding a quasipolynomial attack that *tightens* the security
requirement rather than formally breaking it (Remark 4: a deliberately **biased**
response may restore security — the open problem the project's Direction-4
defense targets). Arapinis et al. (1910.02126) supply the unforgeability
definitions that bound "how much an adversary must learn to forge."

---

## Cross-cutting candidates fetched but not finalized

These arXiv IDs were retrieved during the search but did not enter the final
adversarial-verification set; confirm metadata before citing:
`2104.06244`, `2110.09469`, `2210.17545` (QPUF / unforgeability), and
`2212.04471` (Caro — Pauli-transfer-matrix process learning, present as
`caro2023verification`).

## Open follow-ups
- Additional Huang–Chen–Preskill shadow / property-prediction works for
  non-local process observables not surfaced here.
- Works proving QSQ/QPSQ lower bounds that *fully* (not partially) resolve
  CR-QPUF unforgeability beyond Pirnay et al. and the focal paper.
- Final peer-reviewed venues/DOIs for the still-preprint anchors (2405.12085,
  2406.04250, 2310.17716) before locking `refs.bib`.

_Provenance: distilled from a fan-out deep-research pass (22 primary sources,
101 extracted claims, 25 adversarially verified). Verbatim per-claim evidence
retained in the run transcript._
