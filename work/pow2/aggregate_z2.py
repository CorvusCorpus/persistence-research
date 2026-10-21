#!/usr/bin/env python3
"""Coverage and totals of z2sieve production logs: every n in [N_LO, N_HI) must lie in a CDONE range.

Prints the coverage, the summed counters, the configurations seen (K, NL, TB, TD, binary) and any CAND/MEMBER line.
usage: aggregate_z2.py N_LO N_HI LOG...      (exit status 0 only if the coverage is complete and nothing was found)
"""
import collections, re, sys

nlo, nhi = int(sys.argv[1]), int(sys.argv[2])
iv, tot, conf, found = [], collections.Counter(), set(), []
for path in sys.argv[3:]:
    for line in open(path):
        m = re.match(r"CDONE j=\[\d+,\d+\) n=\[(\d+),(\d+)\) (.*)", line)
        if m:
            iv.append((int(m.group(1)), int(m.group(2))))
            for k, v in re.findall(r"(\w+)=(\d+)", m.group(3)):
                tot[k] += int(v)
        elif line.startswith("z2sieve K="):
            conf.add(line.split(":")[0])
        elif line.startswith("# invocation"):
            conf.add(re.sub(r".*exe=(\S+).*", r"exe=\1", line.strip()))
        elif line.startswith(("CAND", "MEMBER")):
            found.append(line.strip())
iv.sort()
gaps, end, overlap = [], nlo, 0
for a, b in iv:
    if a > end:
        gaps.append((end, a))
    if a < end:
        overlap += min(b, end) - a
    end = max(end, b)
if end < nhi:
    gaps.append((end, nhi))
gaps = [(max(a, nlo), min(b, nhi)) for a, b in gaps if min(b, nhi) > max(a, nlo)]
print(f"coverage of [{nlo},{nhi}): {len(iv)} chunks, gaps {gaps[:5]}{'...' if len(gaps) > 5 else ''}, overlap {overlap}")
print("totals:", " ".join(f"{k}={tot[k]}" for k in ("s1", "s2", "s3", "s4", "s5", "cand", "members")))
print("configurations:", sorted(conf))
print("CAND/MEMBER lines:", found if found else "none")
sys.exit(0 if not gaps and not found else 1)
