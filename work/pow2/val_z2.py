#!/usr/bin/env python3
"""Exact check of z2sieve on a region of exponents (docs/ZEROLESS_POW2.md, "Validation").

Runs EXE N_LO N_HI dump threads=T and checks everything it reports against exact arithmetic that shares no code with
it (Python integers, python-flint ball arithmetic instead of MPFR):
  - S2 lines: exactly the n in [N_LO, N_HI) whose last 2K digits contain no 0 (stages 1-2, in both directions);
  - R3: the first TD digits of 2^n contain a 0;
  - R4 / S4: digits 2K+1..9 NL of 2^n contain a 0 / contain none;
  - after S4: R5 needs a 0 among the last 1000 digits, the first 40, or 2^n in full; MEMBER needs 2^n without a 0;
    CAND needs no 0 in those windows;
  - CDONE lines cover [N_LO, N_HI) exactly; the DONE counters equal the exact counts (s1 = n whose last K digits
    contain no 0).
So every n of the region is accounted for: either a 0 of 2^n is confirmed here, or 2^n is a MEMBER confirmed here.
usage: val_z2.py EXE N_LO N_HI [THREADS]      (exit status 0 only if everything matches)
"""
import re, subprocess, sys
import flint
sys.set_int_max_str_digits(0)
FULLMAX = 10_000_000


def lead(n, m):
    """The first m digits of 2^n, by ball arithmetic (None if undecided at 1024 bits)."""
    for prec in (256, 1024):
        flint.ctx.prec = prec
        x = flint.arb(n) * flint.arb(2).log() / flint.arb(10).log()
        k = x.floor().unique_fmpz()
        if k is None:
            continue
        v = (flint.arb(10) ** (x - k + (m - 1))).floor().unique_fmpz()
        if v is not None and len(str(v)) == m:
            return str(v)
    return None


def full_has_zero(n):
    return '0' in str(flint.fmpz(2) ** n)


def main():
    exe, nlo, nhi = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    T = sys.argv[4] if len(sys.argv) > 4 else "4"
    p = subprocess.run([exe, str(nlo), str(nhi), "dump", f"threads={T}"], capture_output=True, text=True)
    if p.returncode:
        sys.exit(f"{exe} failed (exit {p.returncode}): {p.stderr.strip()[-500:]}")
    K, NL, TB, TD, NMIN = map(int, re.search(r"K=(\d+) NL=(\d+) TB=(\d+) TD=(\d+) NMIN=(\d+)", p.stderr).groups())
    ev, cov, done, errors = {}, [], None, []
    for line in p.stdout.splitlines():
        tag, _, rest = line.partition(" ")
        if tag in ("S2", "R3", "R4", "S4", "R5", "MEMBER", "CAND"):
            ev.setdefault(int(rest), []).append(tag)
        elif tag == "CDONE":
            a, b = map(int, re.search(r"n=\[(\d+),(\d+)\)", rest).groups())
            cov.append((a, b))
        elif tag == "DONE":
            done = {k: int(v) for k, v in re.findall(r"(\w+)=(\d+)(?=\s|$)", rest)}
        else:
            errors.append(f"unexpected line: {line}")
    # coverage
    cov.sort()
    if not cov or cov[0][0] != nlo or cov[-1][1] != nhi or any(cov[i][1] != cov[i + 1][0] for i in range(len(cov) - 1)):
        errors.append(f"CDONE lines do not cover [{nlo},{nhi}) exactly: {cov[:3]}...")
    # stages 1-2, exactly, by doubling 2^n mod 10^(2K)
    m2k = 10 ** (2 * K)
    w = pow(2, nlo, m2k)
    s1, s2set = 0, set()
    for n in range(nlo, nhi):
        sw = str(w).zfill(2 * K)
        if '0' not in sw[K:]:
            s1 += 1
            if '0' not in sw[:K]:
                s2set.add(n)
        w = 2 * w % m2k
    got2 = {n for n, t in ev.items() if t[0] == "S2"}
    if got2 != s2set:
        miss, extra = sorted(s2set - got2), sorted(got2 - s2set)
        errors.append(f"S2 set differs: {len(miss)} missing (UNSOUND) e.g. {miss[:5]}, {len(extra)} extra e.g. {extra[:5]}")
    # every reported decision
    nr3 = nr4 = ns4 = nr5 = nmem = ncand = ns5 = nr3_loose = 0
    mod4 = 10 ** (9 * NL)
    for n in sorted(ev):
        t = ev[n]
        if t == ["S2", "R3"]:
            nr3 += 1
            d = lead(n, TD)
            if d is None or '0' not in d:
                errors.append(f"R3 {n}: first {TD} digits {d} have no 0 (UNSOUND)")
            continue
        d = lead(n, TD)
        if d is not None and '0' in d:
            nr3_loose += 1                   # passed stage 3 although its leading digits have a 0: allowed (margins)
        win = str(pow(2, n, mod4)).zfill(9 * NL)[:9 * NL - 2 * K]      # digits 2K+1..9 NL
        if t == ["S2", "R4"]:
            nr4 += 1
            if '0' not in win:
                errors.append(f"R4 {n}: digits {2 * K + 1}..{9 * NL} have no 0 (UNSOUND)")
            continue
        if len(t) != 3 or t[:2] != ["S2", "S4"] or t[2] not in ("R5", "MEMBER", "CAND"):
            errors.append(f"{n}: bad tag sequence {t}")
            continue
        ns4 += 1
        if '0' in win:
            errors.append(f"S4 {n}: digits {2 * K + 1}..{9 * NL} contain a 0 (stage 4 passed it)")
        if n <= 4000:
            ns5 += 1
            z = full_has_zero(n)
        else:
            low = str(pow(2, n, 10 ** 1000)).zfill(1000)
            top = lead(n, 40) or ""
            z = '0' in low or '0' in top
            if not z:
                ns5 += 1
                if n <= FULLMAX:
                    z = full_has_zero(n)
        if t[2] == "R5":
            nr5 += 1
            if not z:
                errors.append(f"R5 {n}: no 0 found in the final windows / in full (UNSOUND)")
        elif t[2] == "MEMBER":
            nmem += 1
            if z or n > FULLMAX:
                errors.append(f"MEMBER {n}: 2^n contains a 0")
        else:
            ncand += 1
            if z or n <= FULLMAX:
                errors.append(f"CAND {n}: should have been decided")
    exp = dict(s1=s1, s2=len(s2set), s3=len(s2set) - nr3, s4=ns4, s5=ns5, cand=ncand, members=nmem)
    if done is None:
        errors.append("no DONE line")
    else:
        for k, v in exp.items():
            if done.get(k) != v:
                errors.append(f"DONE {k}={done.get(k)} but exact {v}")
    members = sorted(n for n, t in ev.items() if t[-1] == "MEMBER")
    print(f"val_z2 {exe} n=[{nlo},{nhi}) K={K} NL={NL} TB={TB} TD={TD}: s1={s1} s2={len(s2set)} R3={nr3} R4={nr4} "
          f"S4={ns4} R5={nr5} members={members} cand={ncand}; stage 3 passed {nr3_loose} with a 0 in the first {TD} "
          f"digits: {'OK' if not errors else 'FAILED'}")
    for e in errors[:20]:
        print("  ERROR", e)
    sys.exit(1 if errors else 0)


main()
