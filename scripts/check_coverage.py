#!/usr/bin/env python3
"""Coverage check of the chain-sieve production logs (CLAUDE.md rule 1, third bullet).

For every completed production run, run work/gsieve/aggregate.py and require: every c in [0, N] covered
for each of the main, direct and odd paths; no MEMBER lines; and the total point count recorded in
docs/CHAIN_SIEVE.md. Add a line to RUNS when a new shell is completed and documented.
"""
import os, re, subprocess, sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GS = os.path.join(REPO, 'work', 'gsieve')
RUNS = [  # (N, N0, total points as recorded in docs/CHAIN_SIEVE.md)
    (36000, -1, 46_710_375_040_091),
    (40000, 36000, 17_356_762_084_000),
]
ok = True
for N, N0, total in RUNS:
    out = subprocess.run([sys.executable, 'aggregate.py', str(N), str(N0), 'logs'], cwd=GS,
                         capture_output=True, text=True, check=True).stdout
    cov = re.findall(r"^(main|direct|odd): c covered (\d+)/(\d+), missing \[\]$", out, re.M)
    full = len(cov) == 3 and all(int(a) == int(b) == N + 1 for _, a, b in cov)
    pts = int(re.search(r"^points total: (\d+)$", out, re.M).group(1))
    none = re.search(r"^members: none$", out, re.M) is not None
    good = full and pts == total and none
    ok &= good
    print(f"N={N} N0={N0}: coverage {'complete' if full else 'INCOMPLETE'}, points {pts:,} "
          f"({'matches' if pts == total else 'DIFFERS FROM'} {total:,}), members {'none' if none else 'PRESENT'}"
          f" -> {'OK' if good else 'FAIL'}")
    if not good: print(out)
sys.exit(0 if ok else 1)
