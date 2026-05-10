"""Sequentially run the four sprint tasks.

Each task is dispatched as a subprocess so its memory footprint is
released cleanly before the next task begins, even if the lab server
becomes contested in between. Each task script has its own resumable
checkpoint, so re-running this orchestrator picks up where it left off.

Usage:
    uv run python run_overnight_sprint.py --seed 0
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parent
LOG_DIR = REPO / "experiments" / "results" / "raw" / "sprint_logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)


# (Script name relative to experiments/, args)
TASKS: list[tuple[str, str, list[str]]] = [
    # (label, script_basename, extra_args)
    ("Task 1: T-count hardness (n=8,10,12)", "12_t_count_hardness.py", []),
    ("Task 2: Coherent drift", "13_coherent_drift.py", []),
    ("Task 3: Innovation heatmap", "11_final_innovation_heatmap.py", []),
    ("Task 4: VQE surrogate scaling (n=6,8)", "14_vqe_scaling.py", []),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--workers", type=int, default=6,
                   help="max-workers passed to each task script")
    p.add_argument(
        "--only", type=str, default=None,
        help="comma-separated 1-based indices of tasks to run; default = all",
    )
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if args.only:
        wanted = {int(x) for x in args.only.split(",")}
        tasks = [t for i, t in enumerate(TASKS, start=1) if i in wanted]
    else:
        tasks = list(TASKS)

    print("=" * 80)
    print(f"Sprint runner   pid={os.getpid()}   tasks={len(tasks)}   workers={args.workers}")
    print("=" * 80)
    sprint_t0 = time.time()
    rc_summary: list[tuple[str, int, float]] = []

    for i, (label, script, extra) in enumerate(tasks, start=1):
        print(f"\n[{i}/{len(tasks)}] {label}")
        cmd = [
            sys.executable, str(REPO / "experiments" / script),
            "--seed", str(args.seed),
            "--workers", str(args.workers),
            *extra,
        ]
        log_path = LOG_DIR / f"{Path(script).stem}.log"
        print(f"      cmd: {' '.join(cmd)}")
        print(f"      log: {log_path}")

        if args.dry_run:
            rc_summary.append((label, 0, 0.0))
            continue

        t0 = time.time()
        with log_path.open("a") as f:
            f.write(f"\n=== {time.strftime('%Y-%m-%d %H:%M:%S')} ===\n")
            f.flush()
            proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT,
                                  cwd=str(REPO), check=False)
        dt = time.time() - t0
        rc_summary.append((label, proc.returncode, dt))
        print(f"      exit={proc.returncode}, wall={dt:.0f}s")
        if proc.returncode != 0:
            print(f"      *** task failed; continuing anyway. ***")

    total = time.time() - sprint_t0
    print("\n" + "=" * 80)
    print(f"Sprint summary  ({total:.0f}s total)")
    print("=" * 80)
    for label, rc, dt in rc_summary:
        status = "OK" if rc == 0 else f"FAIL({rc})"
        print(f"  {status:10s}  {dt:7.0f}s   {label}")
    return 0 if all(rc == 0 for _, rc, _ in rc_summary) else 1


if __name__ == "__main__":
    raise SystemExit(main())
