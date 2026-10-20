"""Driver for csrc/bottomb (cumulative bottom test, one-free-prime bases)."""
import subprocess, os, tempfile
HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, 'csrc', os.environ.get('BOTTOMB', 'bottomb'))

def run(base, pb, Xc, Gb, k, j, heff, orders, skip=0, aklo=0, akhi=10**6):
    b = base.b
    lines = [f"b {b}", f"q {base.q}", f"nd {len(base.bdiv)}"]
    for p in base.bdiv:
        e = 0; x = b
        while x % p == 0: x //= p; e += 1
        lines += [f"p {p}", f"e {e}", f"pb {pb.get(p, 0)}"]
    lines += [f"Xc {int(Xc)}", f"Gb {Gb}", f"k {k}", f"j {j}", "heff " + " ".join(map(str, heff)), f"nord {len(orders)}"]
    lines += [" ".join(str(d - 1) for d in o) for o in orders] + [f"skip {skip}", f"aklo {aklo}", f"akhi {akhi}"]
    fd, path = tempfile.mkstemp(suffix='.job'); os.write(fd, ("\n".join(lines) + "\n").encode()); os.close(fd)
    try:
        p = subprocess.Popen([BIN, path], stdout=subprocess.PIPE, text=True, bufsize=1)
        for line in p.stdout: yield line.split()
        if p.wait(): raise RuntimeError("C bottomb failed")
    finally:
        os.unlink(path)

def candidate_possible(base, bag, cfg, M, c0, betas, K=40):
    """The bottom congruence pinned the exponent: the only y with this configuration (digits cfg below
    height M, 1s up to M) is y = B q^c0.  Return False if y cannot be a solution: one of its digits at
    heights M..M+K-1 is 0 or not among the digits left after removing cfg from the bag (with multiplicity).
    A short y is checked completely (if it were a genuine solution the case stays open)."""
    import math
    from collections import Counter
    import flint
    b, q = base.b, base.q
    B = 1
    for p, e in zip(base.bdiv, betas): B *= p**e
    if c0 * math.log(q) / math.log(b) < M + K + 12:
        y = B * q**c0
        ds = flint.fmpz(y).str(b)
        return '0' not in ds and sorted(int(ch, b) for ch in ds if ch != '1') == sorted(bag)
    rem = Counter(bag) - Counter(cfg)
    yy = B * pow(q, c0, b**(M + K)) % b**(M + K) // b**M
    for _ in range(K):
        d = yy % b; yy //= b
        if d == 0: return False
        if d != 1:
            rem[d] -= 1
            if rem[d] < 0: return False
    return True

def bottom_test_c(j, bag, base, pb, Xc, Gb, bot, st, orders, heff, skip=0, progress=lambda d: None, aklo=0, akhi=10**6):
    ok = None
    possible = False
    for f in run(base, pb, Xc, Gb, len(bag), j, heff, orders, skip, aklo, akhi):
        if f[0] == 'PROG': progress(int(f[1]))
        elif f[0] == 'FAIL': ok = False; st['bottom_fail'] = f
        elif f[0] == 'CAND':
            o = int(f[1]); nb = len(bag) - j
            M = int(f[3 + nb + 1]); c0 = int(f[3 + nb + 3]); betas = list(map(int, f[3 + nb + 5:]))
            if candidate_possible(base, bag, orders[o], M, c0, betas):
                possible = True; st['bottom_fail'] = f
            else:
                st['resolved_bottom'] = st.get('resolved_bottom', 0) + 1
        elif f[0] == 'WEAK':
            # congruence c = c0 (mod N) with N <= Xc: every exponent of the class up to Xc is checked
            o = int(f[1]); nb = len(bag) - j
            M = int(f[3 + nb + 1]); c0 = int(f[3 + nb + 3]); N = int(f[3 + nb + 5]); betas = list(map(int, f[3 + nb + 7:]))
            count = (Xc - c0) // N + 1 if c0 <= Xc else 0
            if count > 200000:
                possible = True; st['bottom_fail'] = f
            else:
                for i in range(count):
                    if candidate_possible(base, bag, orders[o], M, c0 + i * N, betas):
                        possible = True; st['bottom_fail'] = f; break
                else:
                    st['resolved_weak'] = st.get('resolved_weak', 0) + 1
        elif f[0] == 'DONE': ok = True; st['bottom'] += int(f[1])
    if ok is None: raise RuntimeError("no verdict")
    return ok and not possible
