"""Main-path totals of the logged chain-sieve production for m <= 100000, for comparison with one tree-sieve run.

The production covered the main path in pieces: m <= 36000 (two invocations, the first of them killed and later re-run in
session 3), 36000 < m <= 40000, and the six shells 40000 < m <= 100000 (archives). The n0 split partitions the region
exactly and every counter is additive, so the sums must equal the DONE counters of a single run over m <= 100000
(`tsieve_gpu tp13.txt G M 100000`, or `gsieve_k13 main G 100000 0 100001`).
For c in [0, 7540) of the m <= 36000 run the production log has no DONE line (the invocation was killed); its counters
come from the session-3 re-run (validation_pc/task3_cuda/rerun_cuda_final.tar.xz), whose points, s4 and members are
checked here against the production log's last CDONE line of that invocation.
usage (from work/tsieve): python3 main_totals.py [TREE_DONE_LINE_FILE]   (prints the totals; with a file holding a DONE
line, compares every counter)"""
import sys, os, re, tarfile, hashlib

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
LOGS = os.path.join(R, 'work', 'gsieve', 'logs')
KEYS = ('points', 's1', 's2', 's3', 's4', 'gmp', 'full', 'members')

def done_fields(line):
    return {k: int(v) for k, v in re.findall(r'(\w+)=(\d+)(?=\s|$)', line) if k in KEYS}

def done_lines(text):
    return [l for l in text.splitlines() if l.startswith('DONE ')]

def main():
    tot = dict.fromkeys(KEYS, 0); pieces = []
    # m <= 36000: [7540, 7898) and [7898, 36001) from production, [0, 7540) from the re-run (checked against the CDONE)
    t = open(os.path.join(LOGS, 'prod_N36000_n0-1_main_0_7898.log')).read()
    inv = t.split('# invocation start c=')
    killed = [x for x in inv if x.startswith('0 ')][0]
    last = [l for l in killed.splitlines() if l.startswith('CDONE')][-1]
    assert last.split()[1] == '7539', last
    cd = done_fields(last)
    with tarfile.open(os.path.join(R, 'validation_pc', 'task3_cuda', 'rerun_cuda_final.tar.xz')) as tf:
        rr = tf.extractfile('rerun_cuda_final/prod_N36000_n0-1_main_0_7540.log').read().decode()
    d = done_fields(done_lines(rr)[-1])
    for k in ('points', 's4', 'members'): assert d[k] == cd[k], (k, d[k], cd[k])
    pieces.append(('N=36000 c=[0,7540) (re-run; points/s4/members = production CDONE 7539)', d))
    for x in done_lines(t): pieces.append(('N=36000 ' + re.search(r'c=\[\S+\)', x).group(0), done_fields(x)))
    for name in ('prod_N36000_n0-1_main_7898_36001.log', 'prod_N40000_n036000_main_0_11824.log', 'prod_N40000_n036000_main_11824_40001.log'):
        for x in done_lines(open(os.path.join(LOGS, name)).read()):
            pieces.append((name, done_fields(x)))
    for n in range(50000, 100001, 10000):
        base = f'prod_N{n}_n0{n - 10000}'
        sums = dict(l.split()[::-1] for l in open(os.path.join(LOGS, base + '.sha256')).read().splitlines())
        with tarfile.open(os.path.join(LOGS, base + '.tar.xz')) as tf:
            member = f'{base}/{base}_main_0_{n + 1}.log'
            data = tf.extractfile(member).read()
        assert hashlib.sha256(data).hexdigest() == sums[member], member
        ds = done_lines(data.decode())
        assert len(ds) == 1, (member, len(ds))
        pieces.append((member, done_fields(ds[0])))
    for name, d in pieces:
        for k in KEYS: tot[k] += d[k]
        print(f"  {name}: " + ' '.join(f"{k}={d[k]}" for k in KEYS))
    print("TOTAL m <= 100000 main path: " + ' '.join(f"{k}={tot[k]}" for k in KEYS))
    if len(sys.argv) > 1:
        line = [l for l in open(sys.argv[1]).read().splitlines() if l.startswith('DONE ')][-1]
        d = done_fields(line); bad = [k for k in KEYS if d.get(k) != tot[k]]
        for k in KEYS: print(f"  {k}: production {tot[k]}  tree {d.get(k)}  {'ok' if k not in bad else 'MISMATCH'}")
        print("MATCH" if not bad else "MISMATCH")
        return 0 if not bad else 1
    return 0

if __name__ == '__main__':
    sys.exit(main())
