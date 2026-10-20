#!/usr/bin/env python3
"""Base-10 smoke test (ROADMAP P0.3): re-derive tree nodes from scratch and compare with the ledger.

In a fresh scratch copy of work/ whose bagres*/ and ckpt*/ directories are empty, and with any
proofC_<node>.json of the tested nodes moved aside, run `python3 prove_cum.py <nodes>`. Then check:
  1. the table regenerated from the committed records alone (results_table.py) is byte-identical to
     docs/results_base10.md (records and ledger agree);
  2. every re-derived node has exactly the preimage set listed in docs/results_base10.md;
  3. with the re-derived proofC files added, the regenerated table differs from the ledger at most in
     the 'proved length bound' column of the re-derived rows, and every bound is <= 5000 (the length
     up to which Z.json is complete). For 98 and 147 the committed records are gap-chain proofs
     (proofJ_*, bound 187); the cumulative prover proves a smaller bound.
The repository itself is never written to. Usage: scripts/smoke_base10.py [node ...]  (default 98 147)
"""
import json, os, shutil, subprocess, sys, tempfile, time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEDGER = os.path.join(REPO, 'docs', 'results_base10.md')
nodes = [int(a) for a in sys.argv[1:]] or [98, 147]

def table_rows(text):
    """node -> list of the 7 cells of its table row"""
    rows = {}
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        if len(cells) == 7 and cells[2].isdigit():
            rows[int(cells[2])] = cells
    return rows

def preimages(cell):
    return [] if cell == 'none' else [int(x) for x in cell.split(',')]

ledger_text = open(LEDGER).read()
ledger = table_rows(ledger_text)
for s in nodes:
    assert s in ledger, f"node {s} not in {LEDGER}"
for b in ('csrc/topc', 'csrc/topc3', 'csrc/bottomc', 'csrc/bottomc2', 'csrc/bottomc3'):
    assert os.path.exists(os.path.join(REPO, 'work', b)), f"work/{b} missing: run `make helpers` first"

tmp = tempfile.mkdtemp(prefix='smoke_base10_')
work = os.path.join(tmp, 'work')
shutil.copytree(os.path.join(REPO, 'work'), work, symlinks=True,
                ignore=shutil.ignore_patterns('gsieve', '__pycache__', '*.pyc'))
for d in os.listdir(work):
    if (d.startswith('bagres') or d.startswith('ckpt')) and os.path.isdir(os.path.join(work, d)):
        shutil.rmtree(os.path.join(work, d)); os.mkdir(os.path.join(work, d))
os.mkdir(os.path.join(tmp, 'aside'))
for s in nodes:
    f = os.path.join(work, f'proofC_{s}.json')
    if os.path.exists(f): shutil.move(f, os.path.join(tmp, 'aside'))
print(f"scratch copy: {work} (bagres*/ckpt* emptied)")

def regenerate():
    subprocess.run([sys.executable, 'results_table.py'], cwd=work, check=True, capture_output=True, text=True)
    return open(os.path.join(work, 'results.md')).read()

# 1. records -> ledger
assert regenerate() == ledger_text, "table regenerated from the committed records differs from docs/results_base10.md"
print("1. table regenerated from committed records == docs/results_base10.md")

# 2. re-derive
t0 = time.time()
r = subprocess.run([sys.executable, '-u', 'prove_cum.py', *map(str, nodes)], cwd=work, capture_output=True, text=True)
open(os.path.join(tmp, 'prove_cum.log'), 'w').write(r.stdout + r.stderr)
print(r.stdout.rstrip())
assert r.returncode == 0, f"prove_cum.py failed (exit {r.returncode}); log in {tmp}/prove_cum.log"
print(f"   prove_cum.py {' '.join(map(str, nodes))}: {time.time() - t0:.1f}s")
for s in nodes:
    rec = json.load(open(os.path.join(work, f'proofC_{s}.json')))[str(s)]
    want = preimages(ledger[s][6])
    assert rec['preimages'] == want, (s, rec['preimages'], want)
    assert rec['Lbound'] <= 5000, (s, rec['Lbound'])
    print(f"2. node {s}: preimages {rec['preimages'] or 'none'} == ledger; cumulative bound L <= {rec['Lbound']} "
          f"(ledger bound {ledger[s][5]})")

# 3. ledger with the re-derived records
new = table_rows(regenerate())
assert set(new) == set(ledger), "node set changed"
for s, cells in ledger.items():
    if new[s] == cells: continue
    assert s in nodes, f"row of node {s} changed although it was not re-derived"
    diff = [i for i in range(7) if new[s][i] != cells[i]]
    assert diff == [5], f"node {s}: columns {diff} changed (only the bound column may change)"
    assert int(new[s][5]) <= 5000
print("3. regenerated table with re-derived records differs from the ledger only in the bound column of "
      f"re-derived rows: {[(s, ledger[s][5], new[s][5]) for s in nodes if new[s] != ledger[s]]}")
shutil.rmtree(tmp)
print("SMOKE TEST PASSED")
