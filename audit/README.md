# Independent audit kit for the chain sieve

`ref_sieve.py` re-implements the digit tests of the chain sieve (`work/gsieve/gsieve.c`) in plain Python,
from their definitions, so that the fast program's decisions can be audited without trusting its code. It
uses Python integers only, no floating point, and imports nothing from `work/gsieve`. It is slow (about 10⁶
points per second per process) and meant for spot regions.

## What it checks independently

| fast-program stage | reference decision | how |
|---|---|---|
| region bounds (`main`, `direct`, `odd`, shells with `n0`) | point counts | exact integer bounds, including x ≥ 10¹³⁰ (the fast program uses floating-point bounds with a margin) |
| stage 1, bitmap of the last 13 digits | `s1` | last 13 digits of x = 2ᵃ3ᵇ7ᶜ, computed directly |
| stage 2, affine lift, digits 14–26 | `S2` set | last 26 digits of x, computed directly |
| direct path, 54-digit filter | `D4 0` set | last 54 digits of x |
| odd path, b ≤ 39, 54-digit filter | `D4 1` set | last 54 digits of x = 3ᵃ5ᵇ7ᶜ |
| odd path, b ≥ 40, "odd thread" (lemma L4) | `T4` set, thread steps `s1` | thread derived independently (below); every candidate checked digit by digit |
| optional (`--full`): stages 3–5 | members | every survivor's full digits |

The fast program does not test odd candidates with b ≥ 40 one by one. Lemma L4 says that their last 40 digits
can all be odd only on one residue class of b modulo 2³⁸, the odd thread. The kit derives that class
independently of the fast program (`gconst.py`), and the derivation is also the proof that the class is unique:

- Let r be a residue mod 10ᵏ with only odd digits and 5ᵏ | r.
- The lifts r + j·10ᵏ divisible by 5ᵏ⁺¹ are j₀ and j₀ + 5: one odd digit and one even digit.
- So exactly one lift survives at every step, starting from r = 5. There is therefore exactly one 40-digit such
  residue, …9313975971735999337397515773773193359375.
- Writing it as −5ᵗ (mod 2⁴⁰), and using 3 = −5^F3 and 7 = −5^F7, gives the condition: a + c odd, and
  b ≡ t − a·F3 − c·F7 (mod 2³⁸).

`lemmas` also checks every constant in `gconst.h` from its defining property. It checks lemmas L1 (suffix
collapse) and L3 (affine lift) against direct computation, and lemma L2 (exactly 4ᵏ admissible endings)
exhaustively for k ≤ 9.

## Usage (from the repository root, after `make`)

```
python3 audit/ref_sieve.py lemmas
python3 audit/ref_sieve.py main G 36000 35500 35501 --exe work/gsieve/gsieve_k13 --full
python3 audit/ref_sieve.py thread G 36000 0 400 --exe work/gsieve/gsieve_k13 --jobs 4
audit/run_audit.sh            # the standard audit set
```

Each region command prints the reference counts. With `--exe` it also runs the fast binary on the same region
and prints `REFERENCE AGREES` only if the point counts, the survivor sets and the stage counters are equal.
The exit status is 0 only in that case.

## Scope and limits

- The fast program's stages 3–5 (leading-digit table, limb engine, GMP) are not re-implemented; `--full`
  checks their end result instead: whether any survivor is a member. `work/gsieve/val_dump.py` and
  `val_prodrow.py` check the individual stage-3/4 decisions against exact arithmetic.
- Numbers below 10¹³⁰ are covered by `work/gsieve/tiny.py`, an exhaustive enumeration that is already
  independent of the sieve.
- The kit is independent of the fast program's code, but not of the mathematics in `docs/CHAIN_SIEVE.md`:
  the region definitions and the reduction to chain numbers are taken from there.
