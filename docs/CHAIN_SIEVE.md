# The chain sieve: a new lower bound for multiplicative persistence 12

## Statement

Let p(n) be the product of the decimal digits of n. Call a 7-smooth number x ≥ 10 a **chain number**
if x has no digit 0 and does not contain both a digit 5 and an even digit. Write G for the set of
chain numbers. For x = 2^a 3^b 5^d 7^c put m(x) = a/3 + b/2 + c + d ("digit length").

**Theorem (computational).** Every chain number x with m(x) ≤ 100,000 is one of 1,117 known numbers,
the largest being 4613824623662469877463718941499166849983276122112 (49 digits).
(m ≤ 40,000: sessions 1–2, re-executed in session 3; 40,000 < m ≤ 100,000: session 3. See "Production runs".)

**Corollary.** No number below 10^100000 has multiplicative persistence ≥ 12 (any final digit). Hence, if
a(12) exists, a(12) > 2.67·10^100000 (the smallest candidate with 100,001 digits is 2677…7; see
`docs/LITERATURE.md`).

Proof of the corollary. Let n have persistence ≥ 12 and c_i = p^(i+1)(n). Since p(x) < x for x ≥ 10,
c_0 > c_1 > c_2 ≥ p^11(n) ≥ 10. c_1 and c_0 contain no 0 (their digit products are nonzero), so c_1 is a
zero-free number ≥ 10, hence not divisible by 10; therefore c_0 cannot contain both a 5 and an even
digit, i.e. c_0 ∈ G. Every digit d of n contributes at most 1 to m(p(n)) (2: 1/3, 3: 1/2, 4: 2/3, 5: 1,
6: 5/6, 7, 8, 9: 1), so m(c_0) ≤ len(n). If len(n) ≤ 100,000, c_0 is one of the 1,117 known chain numbers;
all of them have persistence ≤ 10, so n has persistence ≤ 11. ∎

Previous record: Tim Peters (OEIS A003001, Sep 2023): if a(12) exists it is > 2.67·10^30000.

## Why chain numbers

Every intermediate number c_0, …, c_8 of a persistence-12 orbit is a chain number. G is much thinner
than the zero-free 7-smooth numbers E: E has 12,606 known members up to 140 digits, G has 1,117 up to
49 digits. Per digit a random number survives the chain condition with probability 0.8 (vs 0.9).

Forms: x ∈ G is never divisible by 10, so x = 2^a 3^b 7^c (form i) or x = 3^a 5^b 7^c with b ≥ 1
(form iii). In form iii the last digit is 5, so every digit must be odd.

## Lemmas used by the sieve

**L1 (suffix collapse).** For a ≥ K, 2^a 3^b 7^c ≡ 2^e (mod 10^K) where e ≥ K and
e ≡ a + b·μ3 + c·μ7 (mod 4·5^(K−1)), μp = log₂ p in (Z/5^K)^×.
(x ≡ 0 mod 2^K and x ≡ 2^e mod 5^K since 2 is a primitive root mod 5^K.)
So one bitmap indexed by e mod 4·5^12 (976,562,500 bits) answers "are the last 13 digits allowed" for
every 7-smooth number with a ≥ 13, and along a row (b, c fixed) the index runs through consecutive
residues, 64 exponents per machine word.

**L2 (exact branching).** For even numbers with the chain condition (digits ∉ {0,5}) the set S_K of
admissible residues has exactly 4^K elements. Lifting from K to K+1 digits, the 5 lifts realise the 5
digits of one fixed parity; exactly one of them is 0 or 5. (Build check: 4^13 = 67,108,864 ✓.)

**L3 (affine lift).** Let e = e13 + P13·j (P13 = 4·5^12, j < 5^13). With y0 = 2^(e13−26) mod 5^26 =
A_lo + 5^13·A_hi, w = (2^P13 mod 5^26 − 1)/5^13:
digits 14–26 of x equal h0 + 2^13·((A_hi + A_lo·w·j) mod 5^13), h0 = ⌊2^13·A_lo/5^13⌋.
(Because (1 + 5^13 w)^j ≡ 1 + 5^13 w j mod 5^26.) j is constant along a row, so the next 13 digits cost
one modular multiply; implemented with AVX-512 (Montgomery multiply, IFMA divisions).

**L4 (odd-branch rigidity).** For x = 3^a 5^b 7^c, x mod 2^K = (−1)^(a+c)·5^(b + a·f3 + c·f7).
If b ≥ k ≥ 3, the next digit has two candidates differing by 5 — one odd, one even — so exactly one
lift survives the all-odd condition. Exhaustive lifting from k = 3 shows there is exactly one 2-adic
thread (sign −1). Hence for b ≥ 40: all last 40 digits odd ⟺ a + c odd and
b ≡ t − a·f3 − c·f7 (mod 2^38). Every such number ends in
…9313975971735999337397515773773193359375. This branch costs O(1) per (a, c).

**L5 (leading digits).** t = a·L2 + b·L3 + c·L7 mod 2^64 with L_p = round(2^64 log10 p) approximates
2^64·frac(log10 x) within (a+b+c)/2 units. A 2^18-bin table marks bins whose every point (with margin
2^−40) has a 0 or 5 among the first 5 significant digits; only those are rejected.

## Pipeline (form i, a ≥ 26, x ≥ 10^130)

1. bitmap (last 13 digits) → 6.87% survive
2. affine lift, digits 14–26 → 0.38% survive
3. leading-digit table (5 digits) → ~40% of those survive
4. exact base-10^9 limb product, digits 27–126
5. GMP: last 600 digits, then the full number

Form i with a < 26: three 18-digit limbs stepped along b. Form iii: b ≤ 39 stepped along a with an
18-digit filter; b ≥ 40 via L4. Numbers below 10^131: exhaustive (4,483,467 numbers).

## Implementation variants and threads (session 2, 20 Oct 2026)

**Stage-2 variants.** The main-path stage 2 (bitmap scan, affine lift, check of digits 14–26) exists in four
compile-time variants (`gsieve.c`, macros `S2_*`; the build picks the first one the target supports). All
later stages are shared, and all four make the same decisions.

| variant | needs | ⌊v/10⁴⌋ | survivor compaction |
|---|---|---|---|
| `ifma` (production build of 19–20 Oct) | AVX-512 F/DQ/BW/VL + IFMA | `vpmadd52huq` with ⌈2⁶⁰/10⁴⌉ (exact for v < 3.8·10¹⁴) | `vpcompressq` |
| `avx512` | AVX-512 F/DQ/BW/VL | 32×32-bit multiply-shift on split operands (below) | `vpcompressq` |
| `avx2` | AVX2 | the same | permutation-table left-pack |
| `scalar` | any C compiler | `%` and `/` | — |
| `cuda` (session 3, `make cuda`) | NVIDIA GPU, CUDA 13 | the same, on the GPU | warp ballot + scan, in index order |

The `cuda` variant is described in its own section below.

**Exact division without IFMA.** This uses the standard multiply-shift division by an invariant integer
(Granlund–Montgomery 1994). Let M₄ = ⌈2⁴⁵/10⁴⌉ = 3,518,437,209 < 2³², so M₄·10⁴ − 2⁴⁵ = 1168.

- *Lemma.* ⌊v·M₄/2⁴⁵⌋ = ⌊v/10⁴⌋ for 0 ≤ v < 2⁴⁵/1168 ≈ 3.0·10¹⁰.
- *Proof.* Write v = 10⁴q + r with 0 ≤ r < 10⁴. Then v·M₄/2⁴⁵ = q + r/10⁴ + 1168·v/(10⁴·2⁴⁵). The last two
  terms are ≥ 0 and their sum is < 1 when 1168·v < 2⁴⁵. ∎
- *Applying it to the 13-digit block.* The block xh = h₀ + 2ᴷ·Y can reach 2⁴⁴, so it is never divided
  directly.
  - Write Y = 10⁴·s + r (Y < 5ᴷ < 2³¹) and w = h₀ + 2ᴷ·r < 2ᴷ·10⁴.
  - Then xh = 2ᴷ·10⁴·s + w, so xh mod 10⁴ = w mod 10⁴ and ⌊xh/10⁴⌋ = 2ᴷ·s + ⌊w/10⁴⌋.
  - The quotients that follow (< 10⁹, then < 10⁵) are divided with the same M₄.
  - Every divided value is therefore < 2³¹, and v·M₄ < 2⁶³.
- *Exhaustive check.* `val_kernel.c` checks the lemma for every v < 2³¹, and the split for every Y < 5ᴷ
  (K = 12, 13) together with 2²⁶ random (Y, h₀) pairs.

**Threads.** With `threads=T` (or `GSIEVE_THREADS=T`), T workers share one copy of the tables (1.1 GB).
They take c values from an atomic counter in increasing order. Each worker hands its per-c results
(dump and MEMBER lines, counter deltas) to a reporter, which writes them in increasing c.

- CDONE lines therefore keep their meaning: cumulative counts of the invocation over every c up to the one
  reported. So `prod.py` resume and `aggregate.py` work unchanged.
- The output is byte-identical for every T, apart from the timing fields of the DONE line.
- The tables are also built by T threads, over word-aligned chunks of exponents. They are the same for
  every T: `GSIEVE_TABLE_HASH=1` prints a fingerprint.
- Differences from the single-threaded program of 19–20 Oct:
  - Every c now gets a CDONE line. Before, a c with an empty main or direct region had none, and in E and
    P modes no c on the odd path had one. The final duplicate CDONE line is kept.
  - `time=` is now the CPU time of all threads.
  - `threads=` and `wall=` are appended to the DONE line, and `variant=` is inserted after `NL=`.
- `prod.py` worker `ALL` runs every path in one process (`THREADS`, `EXE` environment variables).

## CUDA variant (session 3, 20 Oct 2026)

`gsieve_k13_cuda` (`make cuda`; `-DS2_CUDA`, device code in `gsieve_cuda.cu`) runs stages 1–2 of the main path
on an NVIDIA GPU. Everything else is the CPU code of the other variants, unchanged: the region bounds, stages 3–5
(`later_stages`), the direct and odd paths, the reporter and the output format.

**Work flow.**
- **Pieces.** A piece is rows b₀ ≤ b < b₁ of one c, about 2²⁹ region points (`GSIEVE_PIECE`). Every c has at
  least one piece.
- **Workers.** The worker threads (`threads=T`) take pieces in increasing order, each with its own CUDA stream and
  buffers. For its piece, a worker computes the row segments with the code of `main_c` (region bounds, λ, j, r, the
  wrap at P1) and runs stages 1–2 of all of them in one launch.
- **What the GPU returns.** Per segment: s1 = r1 − r0, and the bit positions of the stage-2 survivors in
  increasing order. That is exactly the content of `PBUF` in the CPU variants.
- **Later stages.** The worker walks the segments in row order. It makes the same `later_stages()` calls as `scan()`
  in the CPU variants, with the same arguments; tb and M follow the rows as in `main_c` (M at b₀ is
  3^b₀·7^c mod 10^126 by GMP).
- **Output.** Pieces are written in increasing (c, piece) order, so dumps and logs are byte-identical to the CPU
  variants'. Only the DONE line's timing fields and variant name differ.

**Kernel.** One block of 256 threads per row segment.
- **Stage 1.** The block's first thread computes the rank range [r0, r1) from RK and BM, as the CPU variants do.
- **Stage 2.** For each index i in [r0, r1):
  - Y = (A_hi + β·j) mod 5^K, by one Montgomery multiplication (QINV = −5^−K mod 2³²);
  - the K-digit block xh = h0 + 2^K·Y is checked 4 digits at a time, with ⌊v/10⁴⌋ = (v·M₄) >> 45 on the split
    Y = 10⁴·sq + r.
  - This is the integer arithmetic of `S2_AVX512`/`S2_AVX2`: every divided value is < 2³¹, and there is no
    floating point.
- **Compaction.** Survivors are compacted in increasing i by warp ballots and one warp-level scan per tile of 1024
  indices, through shared memory, and written per segment.
- **No silent drops.** A segment with more than 2048 survivors writes earlier chunks as extra records. A launch with
  more survivors than its buffer holds grows the buffer and runs again.
- **Host checks on every launch.** Survivors must be strictly increasing positions inside their segment, at most s1
  of them, and the totals must equal the device counters. Any violation stops the program.
- 40 registers per thread, no spills, 13.7 KB of shared memory per block: full occupancy on sm_86.

**Fail-stop table checks (non-ECC memory).** GeForce memory and the desktop's RAM have no ECC.
- **Upload.** The tables BM, RK, AB, H0, POS and the window mask (1.12 GB for K = 13, mode G) are uploaded once and
  compared byte for byte with the host copies, which stay allocated.
- **Periodic checks.** Every `GSIEVE_VERIFY` seconds (default 60), and once more at the end, the device tables are
  compared with the host copies again. The host tables that stages 3–5 read (TOP, T2, digit-flag tables) are
  re-hashed.
- **Why this suffices.** The tables are never written after the upload, so a memory fault that changes them
  persists until a check finds it.
- **What gets written.** Each piece is stamped with the epoch in which it finished. A c is written only after every
  one of its pieces finished before a check that then passed. A failed check stops the program, so a log never
  contains results computed after the last passed check, and `prod.py` resumes from its last CDONE line.
- **What it does not cover.** Transient faults in registers, shared memory or transfers; an independent second run
  catches those.

**Throughput** (RTX 3060 Ti + i7-12700F, 20 threads, main path G, table build excluded):

| workload | cuda, 10⁹ points/s | avx2 (CPU only) | ratio |
|---|---|---|---|
| N = 36,000, c ∈ [7540, 7740), 4.8·10¹¹ points | 169–209 | 34.8 | 4.9–6.0 |
| shell N = 50,000, n0 = 40,000, c ∈ [0, 200), 5.4·10¹¹ points | 155 | 32.2 | 4.8 |
| shell N = 100,000, n0 = 90,000, c ∈ [0, 100), 5.7·10¹¹ points | 145 | 30.7 | 4.7 |

- **The CPU limits the pipeline.** Stages 1–2 alone (`-DBENCH2`) run at 2.6–2.8·10¹¹ points/s on the GPU.
- **CPU cost on one P-core**, per stage-2 survivor: stage 3 (leading-digit table) 3.3 ns, stages 4–5 14 ns.
- **Stage 3 stays on the CPU.** Moving it to the GPU would gain at most about 1.3×.
- **Estimate for N = 10⁵.** The main path of all shells 40,000 → 100,000 (about 9.4·10¹⁴ points) takes about
  1.7 h. The direct and odd paths take about 0.2 h on the CPU.

## Validation

- Stage decisions (bitmap+lift, leading-digit rejections, limb engine) compared one by one with exact
  big-integer arithmetic on full regions (N = 150; G and zero-free modes; K = 12 and 13 builds).
- Zero-free mode end to end (main + direct + odd paths, N = 160, numbers ≥ 10^67): reproduces exactly
  the 314 zero-free 7-smooth numbers of that size from the independent list (both builds).
- Direct and odd paths: survivor sets equal to brute force (N = 160).
- Region split by digit length (n0) partitions points and members exactly.
- Small range (< 10^131): exhaustive enumeration gives exactly the 1,117 known chain numbers.
- Production scale: complete slices c = 35,500 and c = 34,200 of the N = 36,000 run re-checked against
  exact arithmetic (val_prodrow.py): stage-2 survivor sets identical, all limb-engine decisions and a
  sample of leading-digit rejections verified. Odd-thread candidates verified to end in the 40-digit thread.
- Survival rates over the full runs match theory to 4 significant figures (0.06872, 0.05497).
- **Session 2 (20 Oct 2026): new stage-2 variants and threads**, on a CPU without IFMA (details in
  `docs/NOTEBOOK.md`, task 2).
  - **Production counters.** The avx512, avx2 and scalar variants each reproduce all eight counters of the
    logged production invocation main N = 36,000, c ∈ [7540, 7898): 8.6·10¹¹ points, s1 through members.
  - **Suites.** The validation suite passes for every variant, at T = 1 and T = 4, and for strict
    `-march=haswell` (AVX2-only) and `-march=x86-64` builds, which contain no AVX-512 instructions.
  - **Production slices** (c = 34,200 and 33,000 of N = 36,000; c = 39,500, 38,200 and 35,990 of the shell
    N = 40,000 with n0 = 36,000): the stage-2 survivor sets equal exact arithmetic.
  - **Cross-variant dumps.** Every variant and thread count produces identical dumps and counters on eight
    regions of up to 10⁹ points.
  - **Threads.** T = 4 output is byte-identical to T = 1, and table fingerprints are equal for every T.
- **Session 3 (20 Oct 2026): CUDA variant**, on the desktop (i7-12700F without AVX-512, RTX 3060 Ti). Details are in
  `docs/NOTEBOOK.md`, task 3; logs in `validation_pc/task3_cuda/`.
  - **Kernel tests.** The GPU's stage-1/2 decisions agree with independent host computations: on every table entry
    for 1024 lift indices (6.9·10¹⁰ decisions), and on 2²⁰ random points each against 2^E mod 10²⁶ (GMP).
  - **Suite.** It passes for the CUDA build at T = 1 and T = 8, with E mode end to end at 314/314.
  - **Production slices.** Eight slices (7.6·10⁸ points) have stage-2 survivor sets equal to exact arithmetic.
  - **Dumps.** The dumps equal those of avx2 and scalar on ten regions, and are byte-identical in order (CDONE lines
    included).
  - **The whole logged computation (6.4·10¹³ candidates) was re-executed twice** with the CUDA variant. Every
    counter of every invocation and every per-c count equals the production logs.
  - **Step 1 on the desktop.** On the same commit, the CPU suite (`val_all.sh`), `make test-prod` and the
    independent reference (`audit/run_audit.sh`, now including regions that contain real odd-thread candidates) all
    pass.

## Production runs

Base run N ≤ 36,000, then the shell 36,000 < m ≤ 40,000 (flag n0=36000). Both: every c covered for
every path (aggregate.py), no chain numbers outside the known list.

| run | main | direct | odd | total |
|---|---|---|---|---|
| m ≤ 36,000 | 46,626,195,781,376 | 33,689,253,108 | 50,490,005,607 | 46,710,375,040,091 |
| 36,000 < m ≤ 40,000 | 17,337,008,704,000 | 7,903,308,000 | 11,850,072,000 | 17,356,762,084,000 |
| **m ≤ 40,000** | | | | **64,067,137,124,091** (+ 4,483,467 below 10^131) |

### Base run detail (N = 36,000)

| path | region | candidates | exceptions |
|---|---|---|---|
| main (2^a3^b7^c, a ≥ 26, ≥ 130 digits) | 2a+3b+6c ≤ 216000 | 46,626,195,781,376 | 0 |
| direct (2^a3^b7^c, a < 26, ≥ 130 digits) | same | 33,689,253,108 | 0 |
| odd (3^a5^b7^c, ≥ 130 digits) | a+2b+2c ≤ 72000 | 50,490,005,607 | 0 |
| tiny (all forms, < 10^131) | all | 4,483,467 | 1,117 known |
| **total** | | **46,710,379,523,558** | |

Stage survival on the main path: 0.06872 after the bitmap (theory 0.8^12 = 0.068719),
0.05497 after the affine lift (theory 0.8^13 = 0.054976). Every c in [0, 36000] is covered by
checkpointed logs for each path (aggregate.py). Wall time ≈ 2 h on 2 cores (Xeon, 2.1 GHz, AVX-512).

### Shells 40,000 < m ≤ 100,000 (session 3, 20 Oct 2026)

Six shells of width 10,000 (`run_shells.sh`, flag n0), run in one session on the desktop (RTX 3060 Ti + i7-12700F)
with the CUDA variant `gsieve_k13_cuda` (`7fcff6cf…`, validated by `val_cuda.sh` and `val_all.sh` on the exact
binaries). Fail-stop table checks every 60 s; MEMBER tripwire; logs archived in `work/gsieve/logs/prod_N*_n0*.tar.xz`
(SHA-256 of every file in the `.sha256` lists, all in `MANIFEST.sha256`). Every c covered for every path, no chain
numbers (`scripts/check_coverage.py`).

| shell | main | direct | odd | total | s1/points | s2/s1 | wall |
|---|---|---|---|---|---|---|---|
| 40,000 < m ≤ 50,000 | 60,979,301,760,000 | 23,398,270,000 | 35,085,180,000 | 61,037,785,210,000 | 0.06871948 | 0.05497557 | 6.1 min |
| 50,000 < m ≤ 60,000 | 90,974,701,760,000 | 28,598,270,000 | 42,885,180,000 | 91,046,185,210,000 | 0.06871948 | 0.05497579 | 8.8 min |
| 60,000 < m ≤ 70,000 | 126,970,101,760,000 | 33,798,270,000 | 50,685,180,000 | 127,054,585,210,000 | 0.06871948 | 0.05497557 | 12.2 min |
| 70,000 < m ≤ 80,000 | 168,965,501,760,000 | 38,998,270,000 | 58,485,180,000 | 169,062,985,210,000 | 0.06871948 | 0.05497578 | 15.2 min |
| 80,000 < m ≤ 90,000 | 216,960,901,760,000 | 44,198,270,000 | 66,285,180,000 | 217,071,385,210,000 | 0.06871948 | 0.05497561 | 19.3 min |
| 90,000 < m ≤ 100,000 | 270,956,301,760,000 | 49,398,270,000 | 74,085,180,000 | 271,079,785,210,000 | 0.06871948 | 0.05497557 | 24.7 min |
| **40,000 < m ≤ 100,000** | 935,806,810,560,000 | 218,389,620,000 | 327,511,080,000 | **936,352,711,260,000** | | | 86 min |
| **m ≤ 100,000 (all runs)** | | | | **1,000,419,848,384,091** (+ 4,483,467 below 10^131) | | | |

Theory: 0.8^12 = 0.06871948, 0.8^13 = 0.05497558. On a binomial scale the s2/s1 deviations are within 3σ
(z = −0.11, +2.31, −0.20, +2.98, +0.41, −0.19); the two larger ones are excesses, i.e. more candidates passed on to
the exact later stages, not fewer. The main path ran at about 2.1·10^11 candidates/s.

**Not yet done:** an independent second pass of these six shells (CPU binary, about 8.6 h:
`run_shells.sh shells "40000:50000 … 90000:100000" cpu`). The desktop has no ECC memory; the fail-stop table checks
cover the tables, and the second pass would also cover transient faults.

## Structure of any counterexample above the bound

Let n have persistence ≥ 12 and c_i = p^(i+1)(n). Then c_0, …, c_8 are chain numbers and
m(c_0) > 100,000. Let t be the largest index ≤ 8 with m(c_t) > 100,000. Either t = 8 (nine chain numbers
beyond the bound), or c_t is a **bridge**: a chain number with more than 69,897 digits (m ≤ log10 x / log10 5, so log10 x > 100,000·log10 5 > 69,897)
whose digit product is one of the 1,117 known chain numbers (< 2^162), so at most 161 of its digits
differ from 1 — over 99.77% ones.

## Heuristic (random-digit model)

Model: digits of 7-smooth numbers behave like independent uniform digits. Calibration against the
complete list below 10^5000: zero-free numbers with 20–39 digits predicted 5,622 / actual 5,581;
40–59: 2,059 / 2,111; 120–139: 3.3 / 4; chain numbers with 40–59 digits: 7.3 / 7.

A persistence-12 number exists iff some 7-smooth number has digit product 4996238671872 or
937638166841712. Expected number of such 7-smooth numbers: 0.019, peaking at 20 digits
(86% within 22 digits, 99.99% within 30). Expected number of new chain numbers beyond 25,000
digits: 10^-2428.
