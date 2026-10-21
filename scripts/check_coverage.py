#!/usr/bin/env python3
"""Coverage check of the chain-sieve production logs (CLAUDE.md rule 1, third bullet).

For every completed production run, run work/gsieve/aggregate.py and require: every c in [0, N] covered
for each of the main, direct and odd paths; no MEMBER lines; and the total point count recorded in
docs/CHAIN_SIEVE.md. Runs whose logs are kept as archives (work/gsieve/logs/prod_N<N>_n0<N0>.tar.xz, see
scripts/archive_logs.sh) are unpacked into a temporary directory first, and every unpacked file must match the
SHA-256 list stored next to the archive. Add a line to RUNS when a new run is completed and documented.
"""
import os, re, subprocess, sys, tarfile, tempfile, hashlib

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GS = os.path.join(REPO, 'work', 'gsieve')
RUNS = [  # (N, N0, total points as recorded in docs/CHAIN_SIEVE.md)
    (36000, -1, 46_710_375_040_091),
    (40000, 36000, 17_356_762_084_000),
    (50000, 40000, 61_037_785_210_000),        # session 3 shells (CUDA variant), logs archived
    (60000, 50000, 91_046_185_210_000),
    (70000, 60000, 127_054_585_210_000),
    (80000, 70000, 169_062_985_210_000),
    (90000, 80000, 217_071_385_210_000),
    (100000, 90000, 271_079_785_210_000),
]

def unpack(N, N0, tmp):
    """unpack logs/prod_N<N>_n0<N0>.tar.xz into tmp, check every file against the .sha256 list; return the log dir"""
    base = f"prod_N{N}_n0{N0}"
    arc, lst = os.path.join(GS, 'logs', base + '.tar.xz'), os.path.join(GS, 'logs', base + '.sha256')
    with tarfile.open(arc) as t: t.extractall(tmp, filter='data')
    expect = dict(reversed(l.split(maxsplit=1)) for l in open(lst).read().splitlines())
    got = {}
    for root, _, files in os.walk(os.path.join(tmp, base)):
        for f in files:
            p = os.path.join(root, f)
            got[os.path.relpath(p, tmp).replace(os.sep, '/')] = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    if got != expect: sys.exit(f"{arc}: unpacked files do not match {lst}")
    return os.path.join(tmp, base)

ok = True
for N, N0, total in RUNS:
    with tempfile.TemporaryDirectory() as tmp:
        archived = os.path.exists(os.path.join(GS, 'logs', f"prod_N{N}_n0{N0}.tar.xz"))
        logdir = unpack(N, N0, tmp) if archived else 'logs'
        out = subprocess.run([sys.executable, 'aggregate.py', str(N), str(N0), logdir], cwd=GS,
                             capture_output=True, text=True, check=True).stdout
    cov = re.findall(r"^(main|direct|odd): c covered (\d+)/(\d+), missing \[\]$", out, re.M)
    full = len(cov) == 3 and all(int(a) == int(b) == N + 1 for _, a, b in cov)
    pts = int(re.search(r"^points total: (\d+)$", out, re.M).group(1))
    none = re.search(r"^members: none$", out, re.M) is not None
    good = full and pts == total and none
    ok &= good
    print(f"N={N} N0={N0}{' (archive)' if archived else ''}: coverage {'complete' if full else 'INCOMPLETE'}, "
          f"points {pts:,} ({'matches' if pts == total else 'DIFFERS FROM'} {total:,}), members {'none' if none else 'PRESENT'}"
          f" -> {'OK' if good else 'FAIL'}")
    if not good: print(out)
sys.exit(0 if ok else 1)
