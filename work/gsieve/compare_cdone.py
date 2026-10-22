# Compare per-c counts (points, s4, members) of two sets of chain-sieve logs, e.g. the production logs and a
# re-run with another binary. CDONE lines are cumulative within an invocation; every pair of consecutive CDONE
# lines gives the counts of an interval of c values. Both log sets are reduced to cumulative totals G(c) at their
# CDONE boundaries, and G is compared at every boundary the two sets share (all production boundaries when the
# re-run has a CDONE line per c). Killed invocations (no DONE line) contribute up to their last CDONE line.
# usage: compare_cdone.py N N0 DIR_A DIR_B [--modes main,direct,odd]   (--modes: compare only these paths)
import collections, glob, sys
args = sys.argv[1:]; MODES = ("main", "direct", "odd")
if '--modes' in args: i = args.index('--modes'); MODES = tuple(args[i + 1].split(',')); del args[i:i + 2]
N, N0, DA, DB = int(args[0]), int(args[1]), args[2], args[3]
KEYS = ("points", "s4", "members")
def intervals(d):
    out = collections.defaultdict(list)              # mode -> [(lo, hi, deltas)], interval = lo <= c <= hi
    for path in sorted(glob.glob(f"{d}/prod_N{N}_n0{N0}_*.log")):
        mode = path.split("/")[-1].split("_")[3]
        prev_c = prev = None
        for line in open(path):
            if line.startswith("# invocation start c="):
                prev_c = int(line.split("c=")[1].split()[0]) - 1; prev = (0, 0, 0)
            elif line.startswith("CDONE"):
                p = line.split(); c = int(p[1]); kv = dict(x.split("=") for x in p[2:])
                cum = tuple(int(kv[k]) for k in KEYS)
                if c > prev_c: out[mode].append((prev_c + 1, c, tuple(a - b for a, b in zip(cum, prev))))
                elif cum != prev: sys.exit(f"{path}: repeated CDONE {c} with different counts")
                prev_c, prev = c, cum
    return out
def cumulative(ivs):
    ivs = sorted(ivs); G = {}; tot = (0, 0, 0); nxt = 0
    for lo, hi, dl in ivs:
        if lo != nxt: sys.exit(f"gap or overlap at c={lo} (expected {nxt})")
        tot = tuple(a + b for a, b in zip(tot, dl)); G[hi] = tot; nxt = hi + 1
    if nxt != N + 1: sys.exit(f"coverage ends at c={nxt - 1}, not N={N}")
    return G
A, B = intervals(DA), intervals(DB)
ok = True
for mode in MODES:
    GA, GB = cumulative(A[mode]), cumulative(B[mode])
    common = sorted(set(GA) & set(GB))
    bad = [c for c in common if GA[c] != GB[c]]
    ok &= not bad and GA[N] == GB[N]
    print(f"{mode}: boundaries {len(GA)} vs {len(GB)}, compared at {len(common)}; totals {dict(zip(KEYS, GA[N]))} vs "
          f"{dict(zip(KEYS, GB[N]))}; " + ("all equal" if not bad else f"{len(bad)} DIFFER, first at c={bad[0]}: {GA[bad[0]]} vs {GB[bad[0]]}"))
print("per-c counts identical" if ok else "MISMATCH")
sys.exit(0 if ok else 1)
