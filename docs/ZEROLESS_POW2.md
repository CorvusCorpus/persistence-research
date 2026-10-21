# Zeroless powers of two: extending the search (OEIS A007377)

## Question and status
- **Question.** A007377 lists the n for which 2^n has no digit 0:
  0, 1, …, 9, 13, 14, 15, 16, 18, 19, 24, 25, 27, 28, 31, …, 37, 39, 49, 51, 67, 72, 76, 77, 81, 86.
  "It is an open problem of long standing to show that 86 is the last term."
- **Previous searches** (re-checked 21 Oct 2026, `docs/LITERATURE.md`):
  - OEIS A007377: "Checked up to k = 10^10" (D. Radcliffe, Aug 21 2022). This is the newest bound found.
  - T. Khovanova's blog (2011) cites a check to 4.6·10^7. A blog comment claims 2.5·10^9, unconfirmed.
- **Not the same problem:** Chaffin and Tooth's searches for runs of zeros in 2^n (A006889).
- **Status of this extension:** the code is built and validated; the production bound is still to come. A result is recorded
  in `docs/RESULTS.md` only after a validated production run with full coverage.

## Reduction
Write n = r + jP with P = 4·5^(K−1) and 0 ≤ r < P. P is the order of 2 modulo 5^K: 2 is a primitive root modulo every
5^k, because 2^4 = 16 ≢ 1 (mod 25).

**Lemma 1 (suffix).**
- *Statement.* For n ≥ K, 2^n mod 10^K = w(r), where w(r) = 2^(r+P) mod 10^K.
- *Proof.* Both 2^n and 2^(r+P) are ≡ 0 (mod 2^K). Modulo 5^K they are equal, because n ≡ r + P (mod P). Apply the
  Chinese remainder theorem. ∎

**Lemma 2 (affine lift).**
- *Notation.*
  - s = 2^r mod 5^(2K), and z₀ = s·2^(−2K) mod 5^(2K), written as z₀ = z_lo + 5^K·z_hi with z_lo, z_hi < 5^K.
  - a₀ = ⌊2^K·z_lo/5^K⌋.
  - u is defined by 2^P ≡ 1 + 5^K·u (mod 5^(2K)), and d = z_lo·u mod 5^K.
- *Statement.* For n ≥ 2K, digits K+1..2K of 2^n (zero padded) form the K-digit number
  U(j) = a₀ + 2^K·t(j), with t(j) = (z_hi + j·d) mod 5^K.
- *Proof.*
  - 2^(jP) = (1 + 5^K u)^j ≡ 1 + j·5^K·u (mod 5^(2K)). Hence y = 2^n mod 5^(2K) ≡ s + j·5^K·(s u).
  - By the Chinese remainder theorem, x = 2^n mod 10^(2K) = 2^(2K)·Z with Z = y·2^(−2K) mod 5^(2K) = z_lo + 5^K·t(j).
  - So x = 2^(2K)·z_lo + 10^K·2^K·t(j), where 2^(2K)·z_lo = w(r) + 10^K·a₀ and a₀ < 2^K.
  - Since a₀ + 2^K·t < 2^K·5^K = 10^K, the upper K digits of x are U(j). ∎
- *Origin.* This is the chain sieve's lemma L3 on the row 3^0·7^0. Along j, U steps by the constant 2^K·d modulo 10^K.

**Lemma 3 (deeper digits).** For any L, 2^n mod 10^L = (2^r mod 10^L)·(2^(jP) mod 10^L) mod 10^L.

**Leading digits.**
- The first m digits of 2^n are those of 10^f, with f = frac(n·log₁₀2).
- The code computes f in 128-bit fixed point: n·⌊log₁₀2·2^128⌋ mod 2^128. The error is below n·2^−128 + 2^−64, so
  below 2^−63 for n < 2^62.

## Pipeline (`work/pow2/z2sieve.c`, K = 13, NL = 18)
| stage | digits tested | how | survivors (fraction of all n) |
|---|---|---|---|
| 1 | last 13 | w(r), doubled along r (no table) | 28.2% |
| 2 | 14..26 | U(j), three-part decimal addition along j, bit tables | 7.1% |
| 3 | first 7 | 2^24-bin table of f (45% of bins reject) | about 3.9% |
| 4 | 27..162 | limb engine on R(r)·Q(j), one column at a time, early exit | about 2·10^−8 |
| 5 | last 1000 and first 40 | GMP and MPFR (rigorous enclosure); 2^n in full if n ≤ 10^7 | 0 expected |

- **Stage 3 table.** A bin rejects only if every 7-digit prefix that the bin can hold contains a 0. The bin's interval is
  widened by 2^−40 in f and by 10^−15 relative in the mantissa (as in the chain sieve's TOP table). The two bins next to
  the wrap-around never reject.
- **Work.** Work is linear in the search limit N (the persistence problems are cubic).
  - Stage 2 runs about 0.28·N steps of a few integer operations.
  - Stage 1 runs P steps per chunk of 1024·P ≈ 10^12 exponents.
- **Exact output.** Every n of [N_LO, N_HI) is rejected by an exact test, or by the stage-3 table with its margins, or
  it is reported:
  - `MEMBER`: 2^n computed in full and found to have no 0;
  - `CAND`: no 0 in the 1162 digits examined, and 2^n too large to compute in full.
  - Neither is expected above 86. Any such line stops production.

## Integrity on non-ECC hardware
A chunk's `CDONE` line is written only after all of the following checks pass; otherwise the program stops and writes
nothing of the chunk. So any `CDONE` line is trustworthy, even if the run is later killed, and the log resumes from
its last one.
- **Tables.** The stage-3 table, the digit tables and the constants are re-hashed, and the Q table (2^(jP) mod
  10^162 for the chunk) as well.
- **Running state.** Every block's running state (w, 2^r mod 5^26 and 2^r mod 10^162) is re-derived by modular
  exponentiation at the end of the block.
- **Stage 2 against stage 1.** At every residue, 2^26·z_lo mod 10^13 must equal w(r).
- **Stage 4 against stages 1–2.** On every candidate whose n is a multiple of 64 (on all candidates in dump mode), the
  limb engine's low 26 digits must equal those of stages 1–2.

## Validation (`work/pow2/val_z2.sh`)
- **K: kernel tests.**
  - At 10^6 random n < 2^56, each build checks stage 1, stage 2, the full limb product, the limb engine's verdict, the
    fixed-point f (within 2^−63) and every stage-3 rejection against GMP and MPFR.
- **E: exact regions** (`val_z2.py`).
  - The dump is checked against Python integers and python-flint ball arithmetic. These share no code with the sieve,
    and arb replaces MPFR.
  - The S2 set must be exact in both directions. Every R3, R4, S4, R5, MEMBER and CAND decision is re-derived.
  - The CDONE lines must cover the region exactly, and the DONE counters must equal the exact counts.
  - Regions:
    - the start of the range;
    - the first j boundary;
    - the chunk boundaries 10^12 = 1024·P and 10^15 = 1024000·P;
    - 3.7·10^13;
    - just below 2^61.
- **C: positive control.**
  - The build K = 4, NL = 2, TD = 3 (NMIN = 60) runs the same code with windows short enough for 2^60..2^86.
  - Over [60, 2·10^6) it must report exactly the known terms 67, 72, 76, 77, 81, 86, with every other n's 0 confirmed.
- **T: threads.** Dumps at 1 and 8 threads are identical.
- **Results:** `validation_pc/pow2/`. See `docs/NOTEBOOK.md`, section P.

## Production
To be done. The CPU build first, then a CUDA port: each step of stage 2 is a few integer additions and three bit lookups
in tables that fit in shared memory.
