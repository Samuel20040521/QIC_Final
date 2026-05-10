# QPSQ — QIC Final Project

Reimplementation and extension of:

> Chirag Wadhwa and Mina Doosti, *Learning Quantum Processes with Quantum Statistical Queries*, Quantum (2024).

This repository accompanies the final project for *Quantum Information and Computation*. It contains:

- An independent reimplementation of the QPStat oracle and Algorithm 1 (average-case shadow tomography of quantum processes).
- Reproductions of Figures 1 and 2 from the paper.
- Extension experiments: noisy-channel robustness, observable-weight scaling, circuit-class comparison, and a small CR-QPUF authentication-attack demo.
- The written report and recorded video presentation.

## Quickstart

```bash
# Sync the environment (Python 3.11, Qiskit 1.x)
uv sync

# Run the tests
uv run pytest

# Reproduce paper figures
uv run python experiments/01_reproduce_fig1.py --seed 0
uv run python experiments/02_reproduce_fig2.py --seed 0

# Run extensions
uv run python experiments/03_ext_noisy_channel.py --seed 0
uv run python experiments/04_ext_observable_zoo.py --seed 0
uv run python experiments/05_ext_circuit_classes.py --seed 0
uv run python experiments/06_crqpuf_attack.py --seed 0

# Compile the report (PDF lands at report/main.pdf)
cd report && pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

All experiment scripts accept `--seed` and write deterministic CSVs under `experiments/results/` and figures under `report/figs/`.

## Project layout

```
src/qpsq/        # core library
experiments/     # numbered scripts that produce paper figures
tests/           # pytest unit tests
report/          # LaTeX source for the final PDF
notebooks/       # exploratory work
doc/             # the paper PDF, the project guideline PDF, and their
                 # markdown conversions under doc/md/
```

## Status

| phase | deliverable | status |
| --- | --- | --- |
| A | uv project bootstrap | done |
| B | literature survey skeleton (refs.bib, ~30 citations) | done |
| C | `src/qpsq/` core library + 10 unit tests | done |
| D | reproduction of Figs 1 and 2 | done |
| E | extension experiments (noisy / observable zoo / circuit classes / CR-QPUF) | done |
| F | report skeleton (`report/main.tex` + compiles to `main.pdf`) | done; needs team polish + YouTube link |
| F | video presentation | TODO before 2026-06-11 |

See `/home/b11202015/.claude/plans/according-to-out-md-finalproject-guideli-immutable-quilt.md` for the full plan and timeline.

## Reference (do not copy)

The authors' code repository is at <https://github.com/chirag-w/qpsq-learning>. We use it only as a sanity-check reference; everything in `src/qpsq/` is reimplemented from scratch and any inherited mathematical formula is cited inline.
