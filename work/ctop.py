"""Python driver for the C top-side loop (csrc/top).  Builds the job from a FastTop lattice, runs the
binary, and re-checks every level-2 survivor exactly with valid_points (cached by target m)."""
import subprocess, os, tempfile, itertools
from fractions import Fraction
import flint
from solver import valid_points, _gs, close_points, _attach

HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, 'csrc', 'top')

def lattice_job(top):
    """Exact data for the level-1 and level-2 tests of a FastTop lattice."""
    lat = top.lat
    if not hasattr(lat, '_c_data'):
        n = lat.dim
        gs2, mu = _gs(lat.B)
        assert [Fraction(x) for x in gs2] == [Fraction(x) for x in lat.gs2]
        Bi = flint.fmpq_mat(flint.fmpz_mat(lat.B)).inv()
        row = [Fraction(int(Bi[n - 1, i].p), int(Bi[n - 1, i].q)) for i in range(n)]
        D = 1
        for x in row: D = D * x.denominator // __import__('math').gcd(D, x.denominator)
        a = [int(x * D) for x in row]
        assert Fraction(a[n - 1], D) == Fraction(top.qa, top.qb)
        mu32 = mu[n - 1][n - 2]
        R2 = Fraction(lat.R2)
        lat._c_data = dict(C=lat.C, qa=top.qa, qb=top.qb, K1=top.K1, K2=top.K2, D=D, a2=a[n - 2], a3=a[n - 1],
                           P=mu32.numerator, Q=mu32.denominator,
                           g2n=gs2[n - 2].numerator, g2d=gs2[n - 2].denominator,
                           g3n=gs2[n - 1].numerator, g3d=gs2[n - 1].denominator,
                           Rn=R2.numerator, Rd=R2.denominator)
    return lat._c_data

def run_top(top, j, ranges, orders, skip=0, progress=None, maxsurv=10**7, debug=False):
    """ranges[i] = known[i] (gap g_i in [1, ranges[i])); orders = list of digit tuples (length j).
    Yields parsed output lines."""
    d = lattice_job(top)
    prec = max(top.prec, d['C'].bit_length() + 64)
    lines = [f"prec {prec}"]
    for kname in ('C', 'qa', 'qb', 'K1', 'K2', 'D', 'a2', 'a3', 'P', 'Q', 'g2n', 'g2d', 'g3n', 'g3d', 'Rn', 'Rd'):
        lines.append(f"{kname} {d[kname]}")
    lines.append(f"j {j}")
    lines.append("ranges " + " ".join(map(str, ranges)))
    lines.append(f"nord {len(orders)}")
    for o in orders: lines.append(" ".join(str(9 * (x - 1)) for x in o))
    lines.append(f"skip {skip}"); lines.append(f"maxsurv {maxsurv}"); lines.append(f"debug {int(debug)}")
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    try:
        p = subprocess.Popen([BIN, path], stdout=subprocess.PIPE, text=True, bufsize=1)
        for line in p.stdout:
            yield line.split()
        rc = p.wait()
        if rc: raise RuntimeError(f"C top failed rc={rc}")
    finally:
        os.unlink(path)

def config_N(digs, gaps):
    s = sum(gaps); N = 10**s; acc = 0
    for x, g in zip(digs, gaps):
        acc += g; N += 9 * (x - 1) * 10**(s - acc)
    return N, s

def candidates(lat, m):
    """Lattice points within R of (0,..,0,-m) that pass the integrality and box filters of valid_points
    (these filters do not depend on 1+z).  If empty, valid_points(lat, m, opz) == [] for every opz."""
    _attach(lat)
    n = lat.dim
    pts = close_points(lat.B, lat.gs2, lat._mu, lat._Binv, [0] * (n - 1) + [-m], lat.R2)
    good = []
    for v in pts:
        if n == 4:
            if v[0] % lat._lt or v[2] % lat._lw: continue
            t, up, w, last = v[0] // lat._lt, v[1], v[2] // lat._lw, v[3]
        else:
            if v[1] % lat._lw: continue
            t, up, w, last = 0, v[0], v[1] // lat._lw, v[2]
        num = t * lat.rC[2] + up * lat.rC[3] + w * lat.rC[7] - last
        if num % lat.rC[10]: continue
        L = num // lat.rC[10]
        if not (0 <= t <= lat.T and 2 <= up <= lat._U and 0 <= w <= lat._W and 1 <= L <= lat.XL): continue
        good.append((t, up, w, L))
    return good

def prove_top_c(j, bag, G, top, known, stats, orders, skip=0, progress=lambda d: None):
    """Same statement as prove_iter2.prove_top for j >= 1; returns True iff every configuration is killed.
    C kills configurations by exact level-1/level-2 lattice tests; every surviving run of configurations
    (consecutive last gaps sharing one target m) is re-checked here: its target m is recomputed in Python
    (interval arithmetic) at both ends of the run, and the run's lattice points within R are enumerated
    exactly (Fincke-Pohst); configurations with a candidate point get the full valid_points check."""
    lat = top.lat
    ranges = [known[i] for i in range(j)]
    cache = {}
    ok = True
    for f in run_top(top, j, ranges, orders, skip=skip):
        if f[0] == 'SURV':
            o = int(f[1]); pre = list(map(int, f[2:2 + j - 1])); lo, hi, mC = int(f[-3]), int(f[-2]), int(f[-1])
            stats['enum'] += 1
            for gl in (lo, hi):
                N, s = config_N(orders[o], pre + [gl])
                if top.target(N, s) != mC: raise RuntimeError(f"C/Python target mismatch at {o} {pre} {gl}")
            if mC not in cache: cache[mC] = candidates(lat, mC)
            if cache[mC]:
                for gl in range(lo, hi + 1):
                    N, s = config_N(orders[o], pre + [gl])
                    if valid_points(lat, top.target(N, s), Fraction(N, 10**s)): ok = False
        elif f[0] == 'PROG':
            progress(int(f[1]))
        elif f[0] == 'DONE':
            stats['top'] += int(f[1]); stats['l1fail'] = stats.get('l1fail', 0) + int(f[2])
    return ok
