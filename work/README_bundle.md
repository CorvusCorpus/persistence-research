# Multiplicative persistence: computer-assisted proofs for nonzero final digits

Contents
- `results.md`     base 10: every 7-smooth preimage of each proved node of the even-final-digit trees
- `smallbases.md`  bases 3, 4, 6: complete antecedent trees (theorems) for nonzero final digits
- `METHOD.md`      lemmas, method, validations, caveats (sections 1-4 describe this session's additions)
- `certificate.json`, `proofJ_*.json`, `proofC_*.json`, `base{3,4,6}_cum.json`  per-node / per-digit-set proof records

Pipeline (base 10): `prove_cum.py <node>` (cumulative method, C loops `csrc/topc`, `csrc/bottomc2`);
earlier nodes: `prove_iter2.py` / `prove_iter3.py` (gap chain). Small bases: `prove_cum_b.py <b> <D>`
(C loops `csrc/topb`, `csrc/bottomb`). Build: `gcc -O2 -march=native -o csrc/X csrc/X.c -lmpfr -lgmp -lm`
(sources: topc <- top.c, bottomc2 <- bottom2.c, bottomc <- bottom.c, topb <- topb.c, bottomb <- bottomb.c,
top <- top_v1.c, bottom <- bottom_v1.c). Validation scripts: `val_*.py`, `planted_*.py`, `audit_M.py`.

Inputs: `Z.json` (all zero-free 7-smooth numbers below 10^5000, from the exhaustive sieve) is not included
(size); `trees_b6_D800.json` holds the base-6 analogue (all zero-free 5-smooth numbers below 6^800).

Caveats: computer-assisted, not peer reviewed; relies on Matveev's theorem, the exhaustive lists, the
2-adic/p-adic valuation recursion, and the correctness of the code. Final digit 0 (where the longest
persistences live, e.g. 277777788888899 with persistence 11) is not addressed.
