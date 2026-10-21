# Zero-free 7-smooth numbers: extending the search (OEIS A238985)

## Question and status
- **Question.** A238985 lists the zeroless 7-smooth numbers: 7-smooth numbers whose decimal expansion has no digit 0.
  It "seems finite and contains 12615 terms", the largest being 2^25·3^227·7^28 (140 digits). Its keywords fini/full were
  removed in 2023 because finiteness is only conjectured. The entry states no search bound.
- **Previous bounds.**
  - B. Chaffin (OEIS A003001, 29 Sep 2016): every 7-smooth number in (10^140, 10^20000) contains a 0.
  - This repository: `work/Z.json` is complete below 10^5000 (session 1, `zf7.c` + `verify.py`).
- **Status of the extension:** see "Production" below. The claim is recorded in `docs/RESULTS.md` only once the run's
  coverage, members and validation all check out.

## Reduction
A 7-smooth x without a digit 0 is not divisible by 10, so x = 2^a 3^b 7^c (form i) or x = 3^a 5^b 7^c with b ≥ 1
(form iii).

With the digit length of the chain sieve, m(x) = a/3 + b/2 + c (form i) or a/2 + b + c (form iii):
- form i: m ≤ log10(x)/log10(7), since 1/(3 log10 2) = 1.107, 1/(2 log10 3) = 1.048, 1/log10 7 = 1.183;
- form iii: m ≤ log10(x)/log10(5), since 1/(2 log10 3) = 1.048, 1/log10 5 = 1.431, 1/log10 7 = 1.183.

So for x < 10^D it suffices to search form i with m ≤ N1 = ⌈D/log10 7⌉ + 1 and form iii with m ≤ N3 = ⌈D/log10 5⌉ + 1.
These are the regions of the chain sieve's main + direct paths and its odd path, respectively, in zero-free mode
(`E`). The sieve covers x ≥ 10^130; smaller x are in `work/Z.json`.

## Methods
The chain sieve's zero-free mode was its validation mode: it reproduces the 314 known zero-free numbers of the size
range of `val_e2e.py`. At production scale two of its parts were too slow. About 7% of the candidates pass stage 2,
against 0.38% in chain mode, and the odd path tested every point and called GMP for 0.4% of them.

1. **Stages 3–4 on the GPU (main path, `GSIEVE_GPU34=1`).**
   - The CUDA kernel also runs `later_stages()`: the leading-digit table, then the limb engine `limb_check()`
     (digits 26..125). Only candidates that pass every digit test go to GMP on the CPU.
   - The arithmetic is the CPU code's integer arithmetic: u32 limbs, u64 sums of at most 14 products < 10^18,
     division by the constant 10^9, and digit windows through the same allowed-window mask.
   - The host receives per-segment stage counts. In dump mode it also receives every stage-2 survivor with its
     stage-3/4 class, so dumps are byte-identical to the CPU variants'.
   - About 2.5·10^10 candidates/s in zero-free mode, against 3.2·10^9 with stages 3–4 on the CPU. It also gives
     2.6·10^11/s in chain mode.
2. **2-adic suffix collapse for form iii (`odd2`, `GSIEVE_ODD2=1`).**
   - **Lemma.** Let x = 3^a 5^b 7^c with b ≥ 26, and write 3 ≡ −5^f3, 7 ≡ −5^f7 (mod 2^26). Then
     x mod 10^26 = 5^26·z with z = (−1)^(a+c)·5^e mod 2^26 and e = (b + a·f3 + c·f7 − 26) mod 2^24.
   - *Proof.* x ≡ 0 (mod 5^26), and x ≡ (−1)^(a+c)·5^(b + a f3 + c f7) (mod 2^26), because 3 and 7 are units ≡ 3 (mod 4)
     and (Z/2^26)^× = {±1} × ⟨5⟩. By the Chinese remainder theorem, x mod 10^26 = 5^26·(x·5^−26 mod 2^26). The power of
     5 is taken mod 2^24, the order of 5 modulo 2^26. ∎
   - **Sieve.**
     - A bitmap of 2 × 2^24 bits (4 MB) records the (sign, e) whose 26 digits contain no 0.
     - Along a row (a, c fixed, b increasing), e increases by 1, so the row is a contiguous run of bits; most points
       cost a fraction of a word operation.
     - Survivors then pass the leading-digit table (t = a L3 + b L5 + c L7), the limb engine on 5^b·(3^a 7^c) (table
       T5[b] = 5^b mod 10^126), and GMP.
     - Rows with b ≤ 25 use the unchanged point-by-point code.
   - **Density.** 0.0807 ≈ 0.9^24: a multiple of 25 ends in 25 or 75, so two of the 26 digits are never 0.
3. **Direct path (form i, a < 26):** unchanged.

## Validation (`work/gsieve/val_cuda.sh`, section Z; `work/gsieve/val_odd2.py`)
- **Stages 3–4 on the GPU.**
  - `val_dump.py` K = 13 and 12, G and E modes: every leading-digit and limb-engine decision is verified against exact
    arithmetic.
  - Dumps are byte-identical to the avx2 build, for example 139,741,177 lines for main E N = 1000.
  - A production slice against exact arithmetic, and the logged production invocation (all eight counters MATCH).
- **odd2.** `val_odd2.py` compares with exact big-integer arithmetic on full regions (N = 160, 300, 400, 1000, a shell
  with n0, and NL = 7 builds that contain odd-form members):
  - point counts equal to the old odd path;
  - the 26-digit survivor set exact;
  - the b ≤ 25 rows exact;
  - every leading-digit rejection and every limb-engine decision verified;
  - members exact (82 odd-form zero-free numbers at N = 160 with NL = 7).
- **End to end.** `val_e2e.py` with both paths active (K = 13 and 12) finds exactly the 314 known zero-free numbers
  of that size range.
- **Positive control at production scale.** The only known zero-free 7-smooth numbers with at least 130 digits,
  2^77·3^64·7^98 (137 digits, main path) and 2^25·3^227·7^28 (140 digits, direct path), must be found, and nothing
  else.

## Production
Driver `work/gsieve/run_zerofree.sh RUNDIR D`: both forms run in parallel. The tripwire fires on any MEMBER line other
than the two expected ones, and the logs are archived to `work/gsieve/logs/zf_*.tar.xz`. Results: (to be filled in)
