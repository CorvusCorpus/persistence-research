# Lab notebook

Dated entries, newest last. Each entry records the command, binary and build flags, git commit, timing
and outcome. Paths are relative to the repository root unless stated otherwise. Statuses follow
`docs/RESULTS.md` (THEOREM / COMPUTATIONAL / HEURISTIC / OPEN); nothing in this notebook upgrades a status
by itself.

---

## 2026-10-20 · Session 2 (cloud container) · Task 1: set up and reproduce

All runs in this entry are on the baseline commit `f5f65ca` with no source changes. The scaffolding
described in 1.7 (top-level `Makefile`, `scripts/`, `MANIFEST.sha256`, this notebook) is in the commit
that adds this file.

### 1.1 Environment
- **CPU:** Intel Xeon @ 2.80 GHz, family 6 model 85 stepping 7 (Cascade Lake), 4 vCPUs (KVM), 15 GB RAM,
  no swap.
  - `grep -o 'avx512[a-z_]*\|bmi2' /proc/cpuinfo | sort -u` gives avx512f, avx512dq, avx512cd, avx512bw,
    avx512vl, avx512_vnni, bmi2. AVX2 and FMA are also present.
  - **No avx512ifma, avx512vbmi or avx512_vbmi2.** The chain sieve needs IFMA (CLAUDE.md lists F, DQ, BW,
    VL, IFMA, VBMI2), so its production main path does not build here (see 1.2).
- **GPU:** none. `nvidia-smi` and `nvcc` are not installed.
- **Toolchain:** Ubuntu 24.04.5 LTS, kernel 6.18.44, glibc 2.39, gcc 13.3.0, GNU Make 4.3, git-lfs 3.4.1.
  - Installed this session: `apt-get install libgmp-dev libmpfr-dev` (GMP 6.3.0, MPFR 4.2.1) and
    `pip install python-flint mpmath` (python-flint 0.9.0, mpmath 1.4.1).
  - Python 3.13.16, numpy 2.5.3.

### 1.2 Build: `cd work && make`
- All CPU helpers (`csrc/topc topc3 bottomc bottomc2 bottomc3 topb bottomb top bottom`) and `zf7` build
  without errors.
- **gsieve fails to build.** `_mm512_madd52hi_epu64` (AVX-512 IFMA) in phases A and B of the main-path
  stage 2 (`gsieve.c:319,333,337`) gives "inlining failed … target specific option mismatch" under
  `-march=native`.
- Consequences:
  - The production binaries cannot be rebuilt or run on this machine.
  - The claim "the production binary rebuilds byte-identically from gsieve.c" cannot be re-checked here.
- **Workaround for validation:** the existing scalar reference path, unchanged source, built with
  `make -C work/gsieve CFLAGS="-O3 -march=native -Wall -Wno-unused-function -DSCALAR"`.
  - `-DSCALAR` replaces only the stage-1/stage-2 scan of the main path, using the `LO` table and
    per-candidate `%` and `/`. Stages 3–5 and the direct and odd paths are the same code as in production.
  - SHA-256 of the binaries:

    | binary | SHA-256 |
    |---|---|
    | `gsieve_k13` | `e6a17105e652589b51e69c4c995b49f5870e9f661c41dd39f1b3379ad6cfe389` |
    | `gsieve_k12` | `98d6857e128a8f4d1396e4be0aa6084f45009529517cc1043720d667fccc5fa5` |
    | `gsieve_k13_nl7` | `840993e5910b5243f88bf4ffb5a530a5843067f1b97e7bf28923412e94e2ff8e` |
    | `gsieve_k12_nl7` | `9fcf909b82194809019b0d2cba5de808bd231b67f79f2cb9babba49f91899b4d` |

### 1.3 Chain-sieve validation suite (scalar build)
- **Command:** `make -C work/gsieve test CFLAGS="-O3 -march=native -Wall -Wno-unused-function -DSCALAR"`.
- **Timing:** 81 s wall.
- **Outcome:** all nine checks pass.

| check | result |
|---|---|
| `tiny.py 131` | 4,483,467 numbers below 10^131; G-members 1,117 = known list (largest 49 digits) |
| `val_dump.py 13 G 150` | 332,946 points; S2 1,271 = Python; R3 732, R4 539, S4 0 verified |
| `val_dump.py 13 E 150` | S2 23,458 = Python; R3 7,790, R4 15,668 verified |
| `val_dump.py 12 E 150` | 336,158 points; S2 29,314 = Python; R3 9,770, R4 19,543, S4 1 (member 2^25·3^227·7^28, 140 digits) |
| `val_dump.py 12 G 150` | S2 2,024 = Python; R3 1,181, R4 843 verified |
| `val_e2e.py ./gsieve_k13_nl7 160 7` | E-mode end to end: expected 314, found 314, missing 0, extra 0 |
| `val_e2e.py ./gsieve_k12_nl7 160 7` | expected 314, found 314, missing 0, extra 0 |
| `val_direct_odd.py ./gsieve_k13_nl7 13 160 7` | direct 2 = 2; odd b ≤ 39: 0 = 0; thread 0 = 0 |
| `val_prodrow.py 36000 35500` | 725,700 points; S2 2,752 = Python; R3 sample 300/1,652, R4 1,100, S4 0 verified |

**Scope.** This validates the scalar main path and the shared stages on this machine. It does not
validate the AVX-512 IFMA binaries that produced the logs, which cannot run here.

**Audit gap.** At N = 160 the thread check in `val_direct_odd.py` is vacuous: the brute-force set
(all-odd last 40 digits, b ≥ 40) and the thread-candidate set are both empty. The thread formula is
otherwise covered by the planted checks in `gconst.py` and the exhaustive-lifting proof of L4, but the
thread path's region bookkeeping has no non-empty test. To fix in task 5 (independent reference
implementation): compare thread candidates with an independent computation on a slice that contains
them, e.g. N = 36000, c < 400, where `thread_sample.txt` records one.

### 1.4 Production-log coverage (no re-run)
- **Commands:** `cd work/gsieve && python3 aggregate.py 36000 -1 logs` and
  `python3 aggregate.py 40000 36000 logs`.
- **Outcome:** the totals in `docs/CHAIN_SIEVE.md` are reproduced exactly.
  - Every c in [0, N] is covered for main, direct and odd, in both runs. Members: none.
  - N = 36000: 46,710,375,040,091 points (main 46,626,195,781,376; direct 33,689,253,108; odd
    50,490,005,607).
  - Shell 36000–40000: 17,356,762,084,000 points (main 17,337,008,704,000; direct 7,903,308,000; odd
    11,850,072,000).
  - Now automated as `scripts/check_coverage.py`.

### 1.5 Base-10 smoke test
- **Procedure:** fresh copy of `work/` in the scratchpad, every `bagres*/` and `ckpt*/` emptied, then
  `python3 prove_cum.py 98 147`.
- **Deviation from the handoff instructions.** There is no `proofC_98.json` or `proofC_147.json` to move
  aside. These two nodes are recorded by the gap-chain prover (`proofJ_98.json`, `proofJ_147.json`, bound
  187 from bag results in `bagres/`), and `docs/results_base10.md` is generated from those records.
  - With `bagres*/` empty, `prove_cum.py` re-derives both nodes with the cumulative method.
- **Timing:** 1 s.
- **Outcome:** pass. The preimage sets agree; only the bound column differs.
  - Node 98, bag (2,7,7): top {0: 51, 1: 51}, bot {2: 51, 3: 2}, every solution has L ≤ 101; preimages
    none.
  - Node 147, bag (3,7,7): same facts, L ≤ 101; preimages none.
  - Both preimage sets equal the ledger's ("none").
  - The table regenerated from the committed records alone (`results_table.py`) is byte-identical to
    `docs/results_base10.md`.
  - With the two re-derived `proofC` files added, only the length-bound column of rows 98 and 147
    changes (187 → 101). The tree-level summary (92 nodes, same open nodes) is unchanged.
  - The cumulative bound being smaller is expected: L ≤ D_1 + a_2 ≤ 50 + 50.
- Now automated as `scripts/smoke_base10.py` (`make smoke`). It never writes to the repository.

### 1.6 Heavier base-10 spot check: `val_cum_e2e.py`
- **Command:** in the scratch copy, `python3 val_cum_e2e.py 2 2,6 3,4 2,2,3 2,7,7 3,7,7 2,2,6,8,9
  2,3,4,8,9`. This is the cumulative prover from scratch in a temp directory, with gap-chain results
  disabled.
- **Timing:** 2 s.
- **Outcome:** pass. Each bag is proved with L ≤ 101 (bag (2): L ≤ 43). Every genuine solution in `Z.json`
  satisfies every proved fact: 0 violations over 3 + 3 + 1 + 1 + 0 + 0 + 1 + 2 solutions.
  - The gap-chain comparison column prints `?` because the script looks for `runJ_*.log` in `work/`
    rather than `work/logs/`. This is cosmetic.

### 1.7 Repository scaffolding (this commit)
- **Top-level `Makefile`:** `make` builds everything; `make test` runs the gsieve suite, the smoke test
  and the record checks.
  - On CPUs without avx512ifma it builds gsieve with `-DSCALAR` and says so. Override with
    `GSIEVE_VARIANT=`.
  - Pitfall: make does not track CFLAGS. Run `make -C work clean` before switching variants.
- **`work/Makefile`:** new `helpers` target (everything but gsieve).
- **`MANIFEST.sha256`:** SHA-256 of 178 record files (all logs, proof records, `certificate.json`,
  small-base proofs, `Z.json`, `G_small.json`, `trees_b6_D800.json`).
  - `make records-check` verifies the manifest and log coverage.
  - `make manifest` regenerates it after new records are added; the diff must show additions only.
- **Logs.** The production logs were already committed as plain files in the baseline (13 MB in total;
  `HANDOFF_PROMPT.md` says so). They stay as they are, since moving them to LFS would rewrite published
  history.
  - `make logs-tarball` builds a reproducible release tarball,
    `dist/persistence-logs-<commit>.tar.xz` (1.4 MB), with its SHA-256, for a GitHub release.
  - Later shells (≈20 MB of logs each near N = 10^5) should go into release tarballs or LFS, not plain
    blobs. Decide before task 4.
- **Validation of the scaffolding:** `make test` from the repository root passes (see 1.8).

### 1.8 `make test` from the repository root
- **Command:** `make test`, which runs the gsieve suite (scalar build, chosen automatically because the CPU
  lacks avx512ifma), `scripts/smoke_base10.py 98 147`, `sha256sum -c MANIFEST.sha256` and
  `scripts/check_coverage.py`.
- **Timing:** 83 s wall.
- **Outcome:** `ALL TESTS PASSED`. `git status` shows no record file modified. `tiny.py` rewrites
  `G_small.json` byte-identically, and only after its assertion passes.

### 1.9 Status after task 1
- **Unchanged.** Every result in `docs/RESULTS.md` keeps its status.
- **Re-confirmed here:**
  - log coverage and totals (1.4);
  - the base-10 records ↔ ledger consistency and two node re-derivations (1.5, 1.6);
  - the chain-sieve decision logic, on the scalar build (1.3).
- **Not re-confirmed here:** the IFMA production binaries themselves (no IFMA on this CPU) and their
  byte-identical rebuild.
- **Plan for task 2:** re-run logged invocations (e.g. main N = 36000, c ∈ [7540, 7898), whose `DONE`
  line records s1 = 58,998,835,763, s2 = 3,243,476,704, s3 = 1,296,156,297, s4 = 0) on the new code
  paths. That compares stage decisions with the production binary at production scale.

---

## 2026-10-20 · Session 2 · Task 2: portable and faster chain sieve on CPU

Same machine as task 1: Cascade Lake, 4 vCPUs, AVX-512 F/DQ/BW/VL, **no IFMA**.

- **Commits:** code in `f3fc492`, write-up in the commit that adds this entry. The base was `17390a3`.
- **Binaries:** built with `make all variants` (gcc 13.3, `-O3 -march=native`) from `gsieve.c` with
  SHA-256 `3f4827dc…`. Their hashes are in `work/gsieve/validation/2026-10-20_session2/binaries.sha256`.
- **The default build here is avx512:** `gsieve_k13` = `1678ff4a…` = `gsieve_k13_avx512`. The avx2 build is
  `3ce4628d…` and the scalar build `f85dc592…`.

### 2.1 What changed in `work/gsieve/gsieve.c`
1. **Stage-2 variants.** The variants are `S2_IFMA` (unchanged production code path), `S2_AVX512`
   (AVX-512 without IFMA), `S2_AVX2` and `SCALAR`. The default is the best one the target supports.
   - On this CPU the default `gsieve_k13` is byte-identical to `gsieve_k13_avx512`.
   - The IFMA path still compiles (`-march=icelake-server`) but cannot run here.
   - IFMA-free division: ⌊v/10⁴⌋ = ⌊v·M₄/2⁴⁵⌋ with M₄ = ⌈2⁴⁵/10⁴⌉, applied to a split of the 13-digit block
     so that every divided value is < 2³¹ (proof in `docs/CHAIN_SIEVE.md`).
   - AVX2 left-packs survivors with a 16-entry `vpermd` table instead of `vpcompressq`.
2. **pthreads.** With `threads=T` or `GSIEVE_THREADS`, the threads share the tables and pick up c values
   from an atomic counter.
   - A reporter writes per-c output in increasing c, so the logs do not depend on T.
   - The tables are built in parallel; `GSIEVE_TABLE_HASH=1` prints a fingerprint of them.
   - `prod.py` gained worker `ALL` and the `THREADS` and `EXE` environment variables.
3. **Output format.**
   - `variant=` is added to the DONE line and the `tables:` line; `threads=` and `wall=` are added to the
     DONE line.
   - `time=` is now the CPU time of all threads.
   - Every c now gets a CDONE line.
4. **New tools:**
   - `val_kernel.c`: exhaustive check of the division lemma.
   - `val_proddone.py`: re-runs logged invocations and compares their DONE counters.
   - `compare_cdone.py`: per-c comparison of two log sets.
   - `val_xvariant.py`: cross-variant dumps.
   - `bench.py`.
   - `rerun_production.sh`.
   - `val_dump.py` and `val_prodrow.py` now take a binary path, and `val_prodrow.py` also takes `n0`.

### 2.2 Validation (identical decisions)
Each check below compares decisions, not only counts, unless it says otherwise.

1. **Division lemma, exhaustively.**
   - Command: `./val_kernel`, 30 s.
   - Outcome: OK. ⌊v·M₄/2⁴⁵⌋ = ⌊v/10⁴⌋ for all v < 2³¹, and the block split holds for all Y < 5ᴷ
     (K = 12, 13) with h₀ ∈ {0, 1, 2ᴷ−1}, plus 2²⁶ random pairs.
2. **Threads vs single thread vs the pre-thread binary.** Full dumps (`dump prog`) at N = 150, for the
   main, direct and odd paths in G and E modes.
   - T = 4 output is byte-identical to T = 1 in all six cases.
   - Against the pre-thread binary, every dump line and the DONE totals are identical and in the same order.
     The only difference is the expected extra CDONE lines: c with an empty main region (9) or direct
     region (36), and every c of the odd path in E mode (151).
   - With `GSIEVE_TABLE_HASH=1`, the table fingerprint is the same for T = 1, 3, 4: `b375760da296addc`
     (K = 13, avx512).
3. **Production DONE line, re-run on 8.6·10¹¹ points.**
   - Command: `val_proddone.py EXE 4 --only main:7540:7898 logs/prod_N36000_n0-1_main_0_7898.log`.
   - Production line: points 858,547,070,292, s1 58,998,835,763, s2 3,243,476,704, s3 1,296,156,297,
     s4 0, gmp 0, full 0, members 0.
   - **avx512, avx2 and scalar** (4 threads each, pre-parallel-build binaries) reproduce **all eight
     counters exactly**, in 108 s, 150 s and 264 s wall.
   - This is a production-scale check, on a different machine, that every stage keeps exactly as many
     candidates as the IFMA production binary did. Equal counts are strong evidence but not set equality.
     Survivor sets are compared directly in items 5–8.
4. **The whole odd path, re-run against production.**
   - Commands: `gsieve_k13 odd G 36000 0 36001 dump threads=2` (213 s wall) and
     `gsieve_k13 odd G 40000 0 40001 n0=36000 dump threads=2`. Both use the final avx512 binary.
   - N = 36000: points 50,490,005,607, s1 (thread steps) 646,596,760, s4 = gmp = 30, full 0, members 0.
     These equal the sums of the two production invocations.
   - Shell: points 11,850,072,000, s1 798,440,760, s4 = gmp = 11, full 0, members 0. This equals the
     production line.
   - So the odd path of the 10^40000 computation has been re-executed in full, with identical counters.
   - The dumps list the 41 thread candidates (T4 lines). All of them reach the GMP check and fail it.
     They give the independent audit (task 5) a non-empty thread set to check.
5. **Validation suite on every build.**
   - Items 5–8 come from one full validation run, 18:47–19:04 UTC. Its driver was then committed as
     `work/gsieve/val_all.sh` (about 15 min on 4 cores).
   - The committed script's own run on the same binaries is logged in
     `work/gsieve/validation/2026-10-20_session2/` and its `.summary.txt`.
   - Builds tested: avx512, avx2 and scalar at T = 4 (via `GSIEVE_THREADS`, so every validation script also
     exercises the threads), and the default build at T = 1. All pass.
   - What the suite covers: `val_kernel`; `tiny.py` (1,117 = known list); `val_dump` at K = 12, 13 in G and
     E modes; E-mode end to end, 314/314 for both K; direct/odd/thread brute force at N = 160; production row
     c = 35500.
6. **Strict-ISA builds.** `-march=haswell` (AVX2 only) selects avx2, and `-march=x86-64` selects scalar.
   - Both builds pass the whole suite.
   - Neither binary contains an AVX-512 instruction (no zmm or mask registers in `objdump`).
   - So the avx2 and scalar variants run on CPUs without AVX-512, as the desktop may be.
7. **Production slices against exact arithmetic** (`val_prodrow.py`, extended to shells). The stage-2
   survivor sets are identical to Python's, and every limb-engine decision is verified.

   | slice | points | S2 survivors |
   |---|---|---|
   | N = 36000, c = 34200 | 9,632,000 | 36,502 |
   | N = 36000, c = 33000 | 26,853,200 | 101,290 |
   | N = 40000, n0 = 36000, c = 39500 | 725,700 | 2,731 |
   | N = 40000, n0 = 36000, c = 38200 | 9,632,000 | 36,628 |
   | N = 40000, n0 = 36000, c = 35990 | 48,044,000 | 181,475 |

   The c = 35990 slice tests the n0 lower bound of the shell.
8. **Cross-variant dumps** (`val_xvariant.py`). These regions are too large for the Python checkers.
   - Binaries compared: the default build at T = 1, and avx512, avx2 and scalar at T = 4.
   - Regions: main G N = 1000 (about 10⁹ points), main E N = 400, main G N = 36000 c = 20000, main G
     N = 40000 c ∈ [35990, 36000) with n0 = 36000, and direct and odd in G (N = 3000) and E (N = 400).
   - Result: the sorted dump lines (S2/R3/R4/S4, D4, T4, MEMBER) and every DONE counter are identical across
     the four.

### 2.3 Benchmark
- **Command:** `cd work/gsieve && python3 bench.py` (best of 2).
- **Workload:** main path G, N = 36000, c ∈ [7540, 7560), 4.854·10¹⁰ points. Table construction is
  excluded. The machine was otherwise idle.

| variant | T | wall s | 10⁹ points/s | ns/point/thread |
|---|---|---|---|---|
| avx512 | 1 | 20.15 | 2.409 | 0.415 |
| avx512 | 2 | 10.33 | 4.699 | 0.426 |
| avx512 | 4 | 5.35 | 9.073 | 0.441 |
| avx2 | 1 | 28.19 | 1.722 | 0.581 |
| avx2 | 2 | 14.20 | 3.418 | 0.585 |
| avx2 | 4 | 7.27 | 6.676 | 0.599 |
| scalar | 1 | 51.54 | 0.942 | 1.062 |
| scalar | 2 | 26.07 | 1.862 | 1.074 |
| scalar | 4 | 13.32 | 3.644 | 1.098 |

- s1 = 3,335,501,241 and s2 = 183,370,745 in all nine runs.
- **Scaling:** 94% efficiency at 4 threads, so memory bandwidth does not limit at this core count.
- **Comparison:**
  - The IFMA production binary did 0.274 ns/point on one core of the 19–20 Oct machine (2.1 GHz Xeon;
    invocation [7540, 7898)).
  - Per core, this machine's best variant is 1.5× slower; over 4 cores it is 2.5× faster than one
    production core.
  - The table build takes 10 s (wall, 4 threads).
- **Stage costs, single thread, avx512** (`-DBENCH2`/`-DBENCH3` builds, c ∈ [7540, 7542)):
  - stages 1–2: 1.33 s;
  - leading-digit table: about 0.4 s, which is 20 ns per stage-2 survivor;
  - limb engine and GMP: about 0.4 s, which is 53 ns per stage-3 survivor.
- **Relevance to task 3.** If stages 1–2 move to the GPU, the CPU still has to spend about 40 ns per
  stage-2 survivor (0.38% of points).
  - On this 4-core machine that caps the pipeline at about 2.6·10¹⁰ points/s, below the GPU estimate of
    5·10¹⁰–2·10¹¹.
  - Moving stage 3 (a 32 KB bitmap lookup) to the GPU as well would roughly halve the CPU work per
    survivor. It would also shrink the survivor stream by 60%, since stage 3 rejects 60% of its input.
  - Measure this on the desktop first.

### 2.4 Production-scale re-run: not done (deliberately)
- **What it would be.** A re-run of the whole logged computation with the new binary: about 2 h on this
  machine.
  - `rerun_production.sh DIR EXE THREADS shell|base|both` re-runs every complete production invocation and
    compares all its counters.
  - It also re-runs the ranges of the two killed invocations, then checks coverage and per-c counts
    (`aggregate.py`, `compare_cdone.py`).
- **Why not now.** It was skipped at the user's request in this cloud session. What has been re-executed is
  the main-path invocation [7540, 7898) (8.6·10¹¹ points, item 3) and the whole odd path (item 4): about
  1.4% of the 6.4·10¹³ candidates. It is the largest
  remaining confidence gain for the 10^40000 claim, and can run on the desktop.
- **Checking the final binary.** `make test-prod` re-checks item 3 with any binary in about 2 min.
  - Its run with the final avx512 binary (parallel table build) is recorded in 2.6.

### 2.5 Status after task 2
- **Done:** ROADMAP P1.1 (AVX2 fallback; also AVX-512 without IFMA) and P1.2 (pthreads with shared tables).
- **Status:** no result in `docs/RESULTS.md` changes status.
- **For the desktop (task 3):**
  1. `make`, then `make test`, `make test-prod` and `work/gsieve/val_all.sh LOGDIR`.
  2. `python3 work/gsieve/bench.py` for the CPU baseline.
  3. The stage-cost note in 2.3 suggests also moving stage 3 to the GPU if the CPU limits the pipeline.
     Decide after measuring.

### 2.6 Final-binary checks (added after the entry above)
- **`make test-prod`** with the final avx512 binary (parallel table build, 4 threads): **all counters match the production DONE lines**.
  `main N=36000 N0=-1 c=[7540,7898): points=858547070292 s1=58998835763 s2=3243476704 s3=1296156297 s4=0 gmp=0 full=0 members=0 | variant=avx512 threads=4 wall=132s | prod time=235.10s -> MATCH`
- **Committed `val_all.sh`.** Its run on the same binaries was stopped after stage C's fourth slice, to save
  usage.
  - Stages A and B and four of the five C slices passed (summary in
    `work/gsieve/validation/2026-10-20_session2.summary.txt`).
  - The checks it did not reach are the same as in the earlier complete run (items 5–8), which passed.
---

## 2026-10-20 · Session 2 · Task 5 (part): independent reference implementation for audits

This is done ahead of tasks 3–4 at the user's request: it needs no GPU, and it closes the vacuous-thread-test
gap of 1.3. Commit: the one that adds this entry. Audit kit: `audit/ref_sieve.py`, `audit/run_audit.sh`,
`audit/README.md`.

### 5.1 What it is
- **Plain Python integers.** No floating point, and nothing imported from `work/gsieve`.
- **Every decision from its definition.**
  - Region bounds use exact integer comparisons, including x ≥ 10¹³⁰.
  - Stage 1 is the last 13 digits and stage 2 the last 26 digits of x = 2ᵃ3ᵇ7ᶜ, computed directly.
  - The direct and odd paths use the last 54 digits; the odd thread uses the last 40.
  - `--full` adds all digits of every survivor, which is chain membership independent of stages 3–5.
- **The odd thread is derived independently** by digit lifting on the residue itself, starting from r = 5.
  - At every step exactly one lift is divisible by 5ᵏ⁺¹ and has an odd new digit. That proves the
    40-digit residue is unique.
  - Writing it as −5ᵗ mod 2⁴⁰ gives t, and 3 = −5^F3, 7 = −5^F7 give F3 and F7.
  - `gconst.py` instead lifts exponent classes for both signs.
- **`lemmas`.** Checks every constant in `gconst.h` from its defining property; random tests of L1 and L3
  against direct computation; L2 exhaustively for k ≤ 9; the thread derivation; and planted thread points.
- **Comparison with the fast program.** `--exe BIN` runs the fast binary on the same region with `dump` and
  requires equal point counts, survivor sets (S2, D4, T4) and counters (s1, s2, s4), plus equal members with
  `--full`.

### 5.2 Results
- **`python3 audit/ref_sieve.py lemmas`:** all OK, in 1 s.
  - Every `gconst.h` constant checks out from its defining property.
  - The independently derived thread gives t = 235,317,736,100 with sign −1, equal to `OTHREAD`. F3 and F7
    equal `OF3` and `OF7`.
  - 300/300 planted thread points pass, and 0/300 shifted points do.
  - L1 and L3 pass 1000 random tests each for K = 12 and 13. L2 gives exactly 4ᵏ for k = 1..9.
- **Full regions against `work/gsieve/gsieve_k13` (avx512), with `--exe` and `--full`:** main G/E at
  N = 150, and direct G/E and odd G/E at N = 160.
  - Points, survivor sets, the s1/s2/s4 counters and members are all equal.
  - The region bounds are computed exactly and agree with the fast program's floating-point bounds.
  - Examples: main G has 332,946 points, s1 = 22,867 and 1,271 survivors. Direct E finds the one member,
    2^25·3^227·7^28.
- **Not run in this session (usage limit):** the production-scale audit `audit/run_audit.sh`, about 15 min on
  4 cores.
  - It covers main slices of N = 36000 and of the shell, direct slices, and odd slices and thread ranges that
    hold real thread candidates (the full odd re-runs of task 2, item 4, list 41).
  - Until it runs, the reference has checked only regions in which the G-mode direct and odd survivor sets
    are empty. So the vacuous-thread gap of 1.3 is not yet closed.
  - Next step at home: `audit/run_audit.sh`.

---

## 2026-10-20 · Session 3 (desktop PC) · Step 1: environment and validation

All runs in this entry use the committed code of `c471a35` (`gsieve.c` SHA-256 `3f4827dc…`). Logs are in
`validation_pc/` (index in its `README.md`). Times are UTC.

### D.1 Environment
- **Machine:** Windows 11 Home 10.0.26200 with WSL 3.0.3.0 (kernel 6.18.54.1-microsoft-standard-WSL2).
- **CPU:** Intel Core i7-12700F (Alder Lake): 8 P-cores with hyper-threading plus 4 E-cores, so 20 logical CPUs.
  16 GB RAM, of which the WSL VM gets 7.7 GB.
  - `grep -o 'avx512[a-z_]*\|avx2\|bmi2' /proc/cpuinfo | sort -u` gives avx2 and bmi2 (also fma, avx_vnni).
  - There is no AVX-512, so the default chain-sieve build is the `avx2` variant.
- **GPU:** RTX 3060 Ti: compute capability 8.6, 38 SMs, 8 GB GDDR6 on a 256-bit bus (448 GB/s), 3 MB L2,
  PCIe 4.0 x16.
  - Driver 610.88 (CUDA 13.3); `nvidia-smi` works inside WSL.
  - It also drives the Windows desktop, which uses 1–2 GB of its memory.
- **WSL distribution `persistence`** (new, dedicated; the user's own Ubuntu 20.04 distribution is untouched).
  - Image: `ubuntu-24.04.5-wsl-amd64.wsl` from releases.ubuntu.com. Its SHA-256 `bb415d82…` equals the one in
    Microsoft's distribution manifest.
  - Installed with `wsl --install --from-file … --name persistence --no-launch`, then set up by the new
    `scripts/setup_wsl_cuda.sh`, run as root.
  - Toolchain: gcc 13.3.0, GNU Make 4.3, GMP 6.3.0, MPFR 4.2.1.
  - Python 3.12.3 in the venv `/opt/venv`, with python-flint 0.9.0, mpmath 1.4.1 and numpy 2.5.4.
  - CUDA 13.3 (`nvcc` V13.3.73) and Nsight Compute 2026.2.1, from NVIDIA's wsl-ubuntu repository. The keyring
    package was checked against the repository index.
  - Same OS release, compiler, python-flint and mpmath as the cloud container of session 2 (which had Python 3.13).
- **How commands are run.** The repository stays on the Windows file system, and everything runs through a login
  shell: `wsl -d persistence -u root --cd <repo> --exec bash -lc '<command>'`.
  - The login shell puts the venv and CUDA on `PATH`.
  - `--exec` stops an outer shell from expanding `$VARIABLES`, which `wsl -- …` does.
- **Line endings (fixed in this commit).** Git for Windows checks text files out with CRLF (`core.autocrlf=true`).
  That breaks `sha256sum -c MANIFEST.sha256`, since the logs are hashed with LF, and it breaks the bash scripts and
  Makefile recipes under WSL.
  - New `.gitattributes` (`* text=auto eol=lf`). The 195 text files were rewritten from their index blobs.
  - `git diff` is empty and `git add --renormalize .` changes nothing, so the content is identical to the commits.

### D.2 Build: `make`
- All helpers, `zf7` and the four gsieve binaries build. The default variant is `avx2`.
- The only warnings are the pre-existing ones: unused `z` in the AVX2 `scan`, and unused `LO`.
- Binary hashes are in `validation_pc/binaries.sha256`, e.g. `gsieve_k13` = `gsieve_k13_avx2` = `46f8308b…`.

### D.3 `make test`
- **Outcome:** `ALL TESTS PASSED`, 72 s wall (`validation_pc/make_test.log`).
- This covers `val_kernel`, `tiny.py` (1,117 = known list), `val_dump` at K = 12, 13 in G and E modes, and E-mode end
  to end at 314/314 for both K.
- It also covers the direct/odd/thread brute force at N = 160, the production row c = 35500, the base-10 smoke test,
  the manifest (all record files intact) and the coverage of both production runs.

### D.4 `make test-prod`
`main N=36000 N0=-1 c=[7540,7898): points=858547070292 s1=58998835763 s2=3243476704 s3=1296156297 s4=0 gmp=0
full=0 members=0 | variant=avx2 threads=20 wall=25s | prod time=235.10s -> MATCH`.

All eight counters equal the IFMA production binary's DONE line. The cloud machine took 132 s on 4 threads.

### D.5 `cd work/gsieve && ./val_all.sh ../../validation_pc`
**ALL CHECKS PASSED** in 307 s (19:47–19:52). The cloud run took about 15 min.
- **A:** suites for avx2 and scalar at T = 4, and for the default build at T = 1.
- **B:** strict-ISA builds. `-march=haswell` selects avx2 and `-march=x86-64` selects scalar; neither contains an
  AVX-512 instruction.
- **C:** five production slices against exact arithmetic (N = 36000: c = 34200, 33000; shell n0 = 36000:
  c = 39500, 38200, 35990).
- **D:** cross-variant dumps on eight regions.

### D.6 Benchmark: `python3 work/gsieve/bench.py`
Main path G, table build excluded, best of 2 (`validation_pc/bench.txt`).

| workload | variant | T | wall s | 10⁹ points/s | ns/point/thread |
|---|---|---|---|---|---|
| c ∈ [7540, 7560), 4.85·10¹⁰ points | avx2 | 1 | 8.29 | 5.86 | 0.171 |
| | avx2 | 2 | 4.26 | 11.39 | 0.176 |
| | avx2 | 4 | 2.25 | 21.57 | 0.185 |
| | scalar | 1 | 26.27 | 1.85 | 0.541 |
| | scalar | 4 | 7.23 | 6.71 | 0.596 |
| c ∈ [7540, 7740), 4.82·10¹¹ points | avx2 | 8 | 17.49 | 27.58 | 0.290 |
| | avx2 | 12 | 15.25 | 31.63 | 0.379 |
| | avx2 | 16 | 13.88 | 34.75 | 0.460 |
| | avx2 | 20 | 13.85 | 34.82 | 0.574 |
| | scalar | 20 | 38.58 | 12.50 | 1.600 |

- s1 and s2 are identical in every run.
- **Per thread:** one P-core runs avx2 2.4× faster than the cloud's avx512 thread (0.171 vs 0.415 ns/point). It is
  1.6× faster than one core of the IFMA production machine (0.274).
- **Whole machine:** 3.5·10¹⁰ points/s. Beyond 8 threads, hyper-threads and E-cores add little.
- **For task 4:** with the CPU alone, the shells 40,000 → 100,000 (about 9.4·10¹⁴ points) would take about 7.5 h.

### D.7 `audit/run_audit.sh` (independent reference, task 5 kit)
- **Outcome:** `AUDIT: reference agrees with work/gsieve/gsieve_k13 on every region`, in 68 s with 20 jobs
  (`validation_pc/audit.txt`).
- The constants and lemmas check out, and the full small regions agree in G and E modes, with members.
- The production slices agree: main N = 36000 (c = 35500, 34200) and shell (c = 39500, 35990), and direct
  (c = 35000, 20000, 1000).
- Odd slices and thread ranges that contain real thread candidates agree:
  - odd N = 36000, c = 153, 1058, 2164, 29958, and shell c = 1787, 36343: one candidate each;
  - thread ranges c ∈ [0, 400) (1 candidate), [1000, 1400) (2) and shell [1700, 1800) (1).
  - Every candidate was checked digit by digit, and with `--full` its full chain membership (none is a member).
- **This closes the audit gap of 1.3 and 5.2:** the reference has now checked non-empty thread-candidate sets.

### D.8 Status after step 1
- Every check passes on this PC with the committed code.
- No result in `docs/RESULTS.md` changes status. The 10^40000 claim keeps its evidence. The odd-thread bookkeeping
  is now independently audited on regions that contain candidates.

---

## 2026-10-20 · Session 3 (desktop PC) · Task 3: CUDA port of stages 1–2

Same machine and WSL distribution as step 1 (D.1). Base commit `21a0f1f`; the code is in the commit that adds this
entry. Logs are in `validation_pc/task3_cuda/`.

### 3.1 What was built
- **`work/gsieve/gsieve_cuda.cu` and `gsieve_cuda.h`:** stages 1–2 of the main path on the GPU.
  - The kernel takes one block of 256 threads per row segment, and computes the rank range [r0, r1) from RK/BM.
  - Stage-2 arithmetic is the `S2_AVX512`/`S2_AVX2` integer arithmetic: a Montgomery multiply mod 5^K, and
    ⌊v/10⁴⌋ = (v·M₄)>>45 on split operands below 2³¹. No floating point.
  - Survivors are compacted in index order through shared memory.
  - Early chunks are written past 2048 survivors per segment, and buffers grow (with a re-launch) on overflow, so
    nothing is dropped.
  - Compile-time checks of QINV, M₄ and Q < 2³¹. ptxas reports 40 registers, no spills, 13.7 KB of shared memory,
    so occupancy is full.
- **`gsieve.c`, variant `S2_CUDA` (`make cuda` builds `gsieve_k1{3,2}[_nl7]_cuda`).**
  - Pieces: rows b₀ ≤ b < b₁ of one c, about 2²⁹ points each. The row segments are computed with `main_c`'s code, the
    GPU runs one launch per piece, then the same `later_stages()` calls follow.
  - A per-piece reporter writes in (c, piece) order. The direct and odd paths and stages 3–5 are unchanged.
  - Host checks on every launch: survivors strictly increasing inside their segment, at most s1 of them, and totals
    equal to the device counters.
- **Fail-stop table checks.** The device tables are uploaded once and compared with the host copies at upload,
  every `GSIEVE_VERIFY` (60) seconds and at the end. The host tables of stages 3–5 are re-hashed at each check.
  - Output is written only for pieces that finished before a passed check, and a failed check stops the program.
  - The reason is that GeForce memory and the desktop's RAM have no ECC.
  - The first version of the upload check reported a mismatch on the 1.2 KB window mask. A pageable `cudaMemcpy`
    had returned before its DMA landed, and the check read the data back on another stream. Fixed with a device
    synchronization.
- **`kerneltest` mode** of the CUDA build, `val_cuda.sh` (the whole CUDA validation), and a `-DBENCH2` path in
  `main_piece` for benchmarks.
- **CPU variants:** unchanged apart from `main()` choosing the worker function.

### 3.2 Validation of the GPU stages (binary before the fail-stop checks, then the final binary in 3.4)
1. **Kernel tests** (`kerneltest`, `validation_pc/task3_cuda/kerneltest.txt`).
   - Whole table, decisions against the SCALAR formula: y = 2^(e−2K) mod 5^2K recomputed by doubling, with `%` and
     `/`, and neither AB/H0 nor the Montgomery form used. Segments have random lengths (a quarter of them 1–128).
   - Random points against the definition: 2²⁰ random (e, j), each against 2^E mod 10^2K computed with GMP.
   - Zero mismatches in every build:

     | build | stage-2 decisions | segments |
     |---|---|---|
     | K = 13 G | 6.87·10¹⁰ (1024 lift indices × 67,108,864 entries) | 1,271,892 |
     | K = 13 E | 4.4·10⁹ | 19,805 |
     | K = 12 G | 4.3·10⁹ | 63,823 |
     | K = 12 E | 9.8·10⁸ | 3,889 |

   - E mode, with about 25% stage-2 survival, exercises the early chunks and the buffer growth with re-launch.
2. **Validation suite** (`make test V=cuda`) at T = 1 and T = 8: every check passes. E mode end to end finds
   314/314 for K = 13 and K = 12, and every stage decision is verified (`suite_cuda_T*.log`).
3. **Production slices against exact arithmetic** (`val_prodrow.py` with the CUDA binary, `prodrow_*.log`): the
   stage-2 survivor sets are identical and the limb-engine decisions are verified, 7.6·10⁸ points in all.
   - N = 36000: c = 34200, 33000, 30000, 25000.
   - Shell n0 = 36000: c = 39500, 38200, 35990, 30000.
4. **Cross-variant dumps** (`val_xvariant.py`, `xvariant.txt`): CUDA at T = 1 and T = 8, avx2 and scalar give
   identical sorted dumps and counters on ten regions.
   - The regions include the c values 7540 (2.4·10⁹ points, five pieces) and 11824–11825 of the shell.
   - Unsorted, the CUDA dumps at N = 150 and N = 1000 (7,363,950 lines) and every CDONE line are byte-identical to
     avx2.
5. **Full re-run of the logged production** (`rerun_production.sh … ./gsieve_k13_cuda 20 both`, 492 s,
   `rerun_cuda.txt`).
   - Every complete invocation of the base run and the shell reproduces all eight DONE counters; the killed ranges
     main [0, 7540) and direct [0, 11381) were also re-run.
   - `aggregate.py`: full coverage of every path, 46,710,375,040,091 and 17,356,762,084,000 points (equal to the
     ledger), no members.
   - `compare_cdone.py`: per-c counts identical at all 35,993 + 39,993 shared boundaries.
   - **So the whole 10^40000 computation (6.4·10¹³ candidates) has now been re-executed**, with an independent
     implementation of stages 1–2 (GPU), on another machine, and with identical results. Session 2 had re-executed
     about 1.4% of it.

### 3.3 Benchmark (`validation_pc/task3_cuda/bench_cuda.txt`)
Main path G, table build excluded.

| workload | T | cuda 10⁹ points/s | avx2 10⁹ points/s |
|---|---|---|---|
| N = 36000, c ∈ [7540, 7740) (4.8·10¹¹ points) | 1 / 2 / 4 / 8 | 23.6 / 50.0 / 93.8 / 129 | 5.9 / 11.4 / 21.6 / 27.6 |
| same | 12 / 16 / 20 | 154 / 164 / 169 (209 in a single run) | 31.6 / 34.7 / 34.8 |
| shell N = 50000, n0 = 40000, c ∈ [0, 200) | 20 | 155 | 32.2 |
| shell N = 100000, n0 = 90000, c ∈ [0, 100) | 20 | 145 | 30.7 |

- **The pre-port estimate** was 5·10¹⁰–2·10¹¹ points/s. Measured: 1.45–2.1·10¹¹, 4.7–6× the CPU.
- **GPU ceiling.** Stages 1–2 alone (`-DBENCH2`) run at 2.6–2.8·10¹¹ points/s. The full pipeline reaches 60–75% of
  that, so the CPU's stages 3–5 limit it.
- **CPU cost on one P-core** (avx2 `-DBENCH2`/`-DBENCH3` builds), per stage-2 survivor: stage 3 costs 3.3 ns and
  stages 4–5 cost 14 ns. On the cloud machine it was about 40 ns in all.
- **Decision: stage 3 stays on the CPU** (task instruction: ask before deviating if it would help). It is under 20%
  of the CPU work, so moving it would gain at most about 1.3×.
- **Nsight Compute** needs the Windows setting "allow access to GPU performance counters" (`ERR_NVGPUCTRPERM`). It was
  not changed; the benchmarks above use timing builds instead.
- **Task 4 estimate:** about 9.4·10¹⁴ main-path points from 40,000 to 100,000. That is about 1.7 h on the GPU
  pipeline plus about 0.2 h for the direct and odd paths, against about 8.6 h on the CPU alone.

### 3.4 Final validation of the binaries to be used for task 4
- **Sources:** `gsieve.c` `6273ada8…`, `gsieve_cuda.cu` `4d94556a…`, `gsieve_cuda.h` `afd01d38…` (in the commit that
  adds this entry).
- **Binaries** (all hashes in `validation_pc/task3_cuda/val_cuda/binaries.sha256`):
  - `gsieve_k13_cuda` = `7fcff6cf…`, built with nvcc 13.3 `-O3 -arch=sm_86` and gcc 13.3 `-O3 -march=native`, linked
    statically to cudart.
  - The CPU binaries rebuilt from the same source are byte-identical to the step-1 builds (`gsieve_k13` =
    `46f8308b…`). So the CUDA changes leave the CPU variants untouched, and their step-1 validation carries over.
- **Runs** (20:42–21:07 UTC, `validation_pc/task3_cuda/final_validation.out`): all three pass.
  1. `./val_cuda.sh` (12 min): **ALL CHECKS PASSED**.
     - Kernel tests: 0 mismatches in all four builds, part 1 and part 2.
     - Suite V=cuda at T = 1 and T = 8.
     - Eight production slices against exact arithmetic.
     - Ten cross-variant regions.
     - The logged invocation main [7540, 7898): MATCH, 7 s.
  2. `./val_all.sh` (CPU variants, 5 min): **ALL CHECKS PASSED**.
  3. `rerun_production.sh … ./gsieve_k13_cuda 20 both` (479 s): the whole logged 10^40000 computation again.
     - Every complete invocation MATCHes; the killed ranges were re-run.
     - Coverage is complete with 46,710,375,040,091 + 17,356,762,084,000 points and no members, and the per-c counts
       are identical.
     - Logs: `rerun_cuda_final.tar.xz` (SHA-256 `f9d4e0d6…`, with the per-file hashes in `rerun_cuda_final.sha256`).
- **Status:** no result in `docs/RESULTS.md` changes status. Claim 1.1 gains a complete independent re-execution
  (GPU stages 1–2, twice) as evidence; it stays COMPUTATIONAL.
- **For task 4.** Production uses `gsieve_k13_cuda` `7fcff6cf…`. This desktop has no ECC memory: the fail-stop table
  checks guard the tables, and an independent second pass guards against the remaining transient faults.

---

## 2026-10-20 · Session 3 (desktop PC) · Task 4: shells 40,000 → 100,000

### 4.1 Setup
- **Code.** Commit `5b599af`: the binaries of task 3 plus the driver `work/gsieve/run_shells.sh` and
  `work/gsieve/check_member.py`.
- **Binaries.** These are the validated binaries of 3.4:
  - `gsieve_k13_cuda` `7fcff6cf…` (production);
  - `gsieve_k13` `46f8308b…` (avx2, reserved for the second pass).
- **Command.** `work/gsieve/shells/launch_task4.sh`, copied to `validation_pc/task4/`, ran
  `bash run_shells.sh shells "40000:50000 50000:60000 60000:70000 70000:80000 80000:90000 90000:100000" both`.
  - It was started detached from Windows with `Start-Process wsl.exe … --exec bash -l shells/launch_task4.sh`, so that
    it did not depend on the Claude session.
  - For each shell, the driver runs `EXE=gsieve_k13_cuda THREADS=20 prod.py N ALL 86400 N0` in
    `work/gsieve/shells/prod_N<N>_n0<N0>/`. It then checks coverage (every c on all three paths, no members) and the
    survival rates (binomial z within ±6), and archives the logs with `scripts/archive_logs.sh` into `work/gsieve/logs/`.
- **Tripwire.** Any MEMBER line under the run directory kills every process of the step's session, records the line,
  checks it with `check_member.py` (exact recomputation and persistence orbit), and exits with status 42.
  - Tested beforehand: a MEMBER line (for 4996238671872) was planted 8 s into a live production run.
  - The run stopped within seconds, with no process left.
  - A first version killed only the process group, and `prod.py`'s `timeout` had moved the binary into its own group.
    It was fixed to kill the whole session.
- **The PC.** The Claude app's keep-awake held while the session was active, and the GPU also drove the desktop.

### 4.2 Results (21:14–22:42 UTC; `validation_pc/task4/driver.log`)

| shell | main | direct | odd | total | members | s1/points (z) | s2/s1 (z) | wall |
|---|---|---|---|---|---|---|---|---|
| 40,000 < m ≤ 50,000 | 60,979,301,760,000 | 23,398,270,000 | 35,085,180,000 | 61,037,785,210,000 | 0 | 0.06871948 (+0.01) | 0.05497557 (−0.11) | 6.1 min |
| 50,000 < m ≤ 60,000 | 90,974,701,760,000 | 28,598,270,000 | 42,885,180,000 | 91,046,185,210,000 | 0 | 0.06871948 (−0.01) | 0.05497579 (+2.31) | 8.8 min |
| 60,000 < m ≤ 70,000 | 126,970,101,760,000 | 33,798,270,000 | 50,685,180,000 | 127,054,585,210,000 | 0 | 0.06871948 (+0.01) | 0.05497557 (−0.20) | 12.2 min |
| 70,000 < m ≤ 80,000 | 168,965,501,760,000 | 38,998,270,000 | 58,485,180,000 | 169,062,985,210,000 | 0 | 0.06871948 (−0.01) | 0.05497578 (+2.98) | 15.2 min |
| 80,000 < m ≤ 90,000 | 216,960,901,760,000 | 44,198,270,000 | 66,285,180,000 | 217,071,385,210,000 | 0 | 0.06871948 (+0.02) | 0.05497561 (+0.41) | 19.3 min |
| 90,000 < m ≤ 100,000 | 270,956,301,760,000 | 49,398,270,000 | 74,085,180,000 | 271,079,785,210,000 | 0 | 0.06871948 (−0.00) | 0.05497557 (−0.19) | 24.7 min |
| **all six** | 935,806,810,560,000 | 218,389,620,000 | 327,511,080,000 | **936,352,711,260,000** | 0 | | | 86 min |

- **Candidates.** With m ≤ 40,000 the total is 1,000,419,848,384,091 ≈ 1.0004·10^15, plus 4,483,467 numbers below
  10^131.
- **Main path:** about 2.1·10^11 candidates/s.
- **GMP checks.** Candidates that reached the GMP check (s4), per shell:
  - main 14, 42, 46, 55, 78, 88;
  - odd 36 to 165;
  - direct 171,177 to 361,496.

  None passed the last-600-digit check (`full=0`), and there are no members.
- **Survival.** On a binomial scale the stage-2 rates are within 3σ. The two larger deviations are excesses, i.e. more
  candidates passed to the exact later stages. Stage 1 matches 0.8^12 to all printed digits in every shell.

### 4.3 Records
- **Archives.** `work/gsieve/logs/prod_N{50000…100000}_n0{40000…90000}.tar.xz` (0.5–1.0 MB each), each with a
  `.sha256` list of its files. All twelve files were added to `MANIFEST.sha256` by `make manifest`; the diff has
  additions only.
- **Coverage check.** `scripts/check_coverage.py` now unpacks archived runs, checks every file against its list, and
  requires full coverage, no members and the recorded totals.
- **`make records-check`:** "all record files intact", and every run (two old, six new) OK.

### 4.4 Second pass: TODO
- After the last shell, the driver started the CPU second pass at 22:42:09. By decision of the user it was stopped at
  once (driver, launcher and the `val_proddone.py` session were stopped by PID and session), and left as a TODO.
  - Resume command: `run_shells.sh shells "…" cpu`, about 8.6 h. The partial `shells/check_N50000_n040000/` is
    rewritten on resume.
  - It is recommended before an OEIS submission, because the desktop has no ECC memory. The fail-stop table checks
    cover the tables; the second pass would cover transient faults and add an independent implementation of stages 1–2.

### 4.5 Status
- `docs/RESULTS.md` 1.1: every chain number with m ≤ 100,000 is known. COMPUTATIONAL, with evidence and the TODO.
- 1.2: no number below 10^100000 has persistence ≥ 12, so a(12) > 2.67·10^100000. THEOREM given 1.1.
- 1.8: bridges have more than 69,897 digits, over 99.77% ones.
- `docs/CHAIN_SIEVE.md`, README, ROADMAP, the paper and the OEIS draft are updated in the same commit.

---

## 2026-10-21 · Session 3 (desktop PC) · Zero-free extension (option 2: OEIS A238985, Chaffin's 10^20000)

Goal: every 7-smooth x < 10^D without a digit 0 is one of the 12,615 known numbers (`work/Z.json`). Chaffin (2016)
showed that none lies in (10^140, 10^20000). Method and lemma: `docs/ZERO_FREE.md`.

### Z.1 Why new code
Zero-free mode (E) was the chain sieve's validation mode. At production scale:
- **Main path:** 3.2·10^9 points/s, entirely CPU-bound. 7.2% of the points pass stage 2, against 0.38% in mode G, and
  the CPU stages 3–4 cost about 52 ns per survivor. D = 40,000 would take about 9 h.
- **Odd path:** about 3·10^6 points/s per thread. It is GMP-bound: 0.42% of the points reach a 600-digit GMP check at
  about 80 µs each. D = 40,000 would take about 12 days.

### Z.2 What was built (commits `bfe9ea0`, `ef5ebde`)
- **`GSIEVE_GPU34=1`:** stages 3–4 of the main path on the GPU, via a CLASSIFY step at every flush of the kernel.
  - It is an exact port of `later_stages()` and `limb_check()`; TOP and T2 are uploaded and covered by the fail-stop
    checks.
  - Dumps are byte-identical to avx2, and DONE lines say `variant=cuda34`.
  - Zero-free mode now runs at 2.5·10^10 points/s, and chain mode at 2.6·10^11 (1.5×).
- **`GSIEVE_ODD2=1`:** the 2-adic suffix collapse for 3^a 5^b 7^c, b ≥ 26 (mode E), on the CPU.
  - A 4 MB bitmap of the last 26 digits (density 0.0807) is walked along rows in b.
  - Then come the leading-digit table, a limb engine on T5[b]·(3^a 7^c) and GMP. Rows with b ≤ 25 use the old code.
  - Speed: 4.8·10^9 points/s at N = 57,240.
- **New tools:** `val_odd2.py`, `val_cuda.sh` section Z, `prod.py` `DM=`, `run_zerofree.sh`. `check_member.py` now
  reports zero-freeness.
- **The CPU binaries are byte-identical to before** (`gsieve_k13` = `46f8308b…`): all new code is under `S2_CUDA` or
  runtime flags.

### Z.3 Validation on the exact binaries (`gsieve_k13_cuda` `9c8d02af…`; `validation_pc/zerofree/`)
- **`val_cuda.sh`, all sections:** kernel tests (0 mismatches), the suite at T = 1 and 8, eight production slices, ten
  cross-variant regions, and the logged invocation (MATCH).
- **Section Z:**
  - `val_dump` with GPU34 for K = 13/12 and G/E: all decisions verified against exact arithmetic.
  - Dumps are byte-identical to avx2 for main G 1000 (7,363,950 lines), main E 1000 (139,741,177 lines) and the shell
    slice.
  - `val_e2e` with GPU34 + ODD2 for K = 13 and 12: **314/314**.
  - `val_prodrow` and `val_proddone` with GPU34: OK and MATCH.
  - `val_odd2` for N = 400, 400 with n0 = 250, and 300 (NL = 7, with 82 members): exact.
- **One step failed, and it was a harness failure, not a wrong decision.** `val_odd2` at N = 1000 was killed by the
  OOM killer, because the Python sets of about 25 M tuples exceed WSL's 7.8 GB. It was replaced by N = 500 (40.5 M
  points), which passes. Earlier, `odd2` at N = 1000 had shown point counts equal to the old odd path, with no members
  in either.
- **Driver tests** (D = 150): both known members were found (main path via GPU34, direct path) and the run ended with
  ALL PASSED. A planted unexpected member set off the tripwire, which stopped everything (exit 42).

### Z.4 Production run D = 40,000 (started 2026-10-21 01:27:53 UTC)
- **Launch:** `zerofree/launch_zerofree.sh` (detached), which runs `run_zerofree.sh zerofree 40000`.
- **Parameters:** N1 = 47,333 (form i: main with GPU34 plus direct, 4 threads) and N3 = 57,229 (form iii: odd2, 16
  threads), in parallel. Commit `ef5ebde`, binary `9c8d02af…`.
- **Expected members:** exactly 2^77·3^64·7^98 (137 digits) and 2^25·3^227·7^28 (140 digits).
- **Results:** in Z.6, once the run ends.

### Z.5 Interrupted by a Windows Update restart; re-check and resume
- **Cause.** Windows Update restarted the PC at 01:29:12 UTC, about 80 s into the run.
  - System log event 1074: `MoUsoCoreWorker.exe`, then `TrustedInstaller.exe`, "Operating System: Upgrade (Planned)".
    The 2026-10 cumulative update KB5127648 (build 26200.9711) and .NET update KB5127690 were installed.
  - There was no crash and no power loss: no Kernel-Power 41, 6008 or BugCheck events.
  - The PC then slept from 01:49 to 07:00 UTC (the AC sleep timeout is 15 min).
- **Written before the kill.**
  - main: c = 0–215 (1.44·10^12 points), including the member 2^77·3^64·7^98 at c = 98.
  - odd: c = 0–76 (2.52·10^11 points, no members).
  - direct: not started.
- **Are those lines checked?**
  - main: yes. A c is written only after a table check covering it has passed (epochs, `gsieve.c`).
  - odd: no. With odd2, CDONE lines are written as each c finishes, and the odd2 tables (bitmap, T5, TOP) are
    re-hashed only at the end of an invocation. The killed invocation never got there.
- **Re-check of odd c = 0–76.** The production command with c1 = 77 (`DM=E GSIEVE_ODD2=1 THREADS=16 JOBS=odd:0:77
  prod.py 57229 ALL 3600`, in `zerofree/recheck_formiii_N57229_c0_77/`), 57 s.
  - It ended with "odd2: tables unchanged at the end of the run" and a DONE line.
  - Its 77 CDONE lines are identical to the killed invocation's. Both invocations start at c = 0 and the counts are
    cumulative, so the lines must agree exactly.
  - Log: `validation_pc/zerofree/recheck_odd_N57229_c0_77.log`.
  - Rule for the rest of the run: any odd-path invocation that ends without a DONE line gets the same re-check before
    its coverage counts.
- **Unchanged by the update:** NVIDIA driver 610.88; binaries `gsieve_k13_cuda` `9c8d02af…` and `gsieve_k13`
  `46f8308b…`; WSL kernel 6.18.54.1. No reboot was pending afterwards.
- **Resumed at 07:10:40 UTC** with `zerofree/launch_zerofree.sh` (detached). prod.py continued main at c = 216 and odd
  at c = 77.
  - A hidden PowerShell process holds a stay-awake power request (`SetThreadExecutionState`) while the run's `wsl.exe`
    is alive (`zerofree/keepawake.txt`). No settings were changed.

## 2026-10-21 · Session 3 (desktop PC) · A007377: zeroless powers of two (harness)

Goal: verify that 2^n contains a 0 for every 86 < n ≤ N far beyond the published 10^10. Method and lemmas:
`docs/ZEROLESS_POW2.md`. The zero-free run (Z.5) was paused meanwhile; this work used at most 8 threads and no GPU.

### P.1 Literature (`docs/LITERATURE.md`, A007377 section)
- The newest search bound is OEIS A007377's "Checked up to k = 10^10" (Radcliffe, Aug 2022). Older figures:
  4.6·10^7 (Khovanova's blog, 2011) and an unconfirmed claim of 2.5·10^9.
- Residue lifting is standard. Only the bound would be new.

### P.2 Code (`work/pow2/`)
- **`z2sieve.c`.** Stages:
  1. the last 13 digits, along r = n mod P, P = 4·5^12, by doubling;
  2. digits 14–26 by the affine lift along j = n div P;
  3. the first 7 digits by a 2^24-bin table with margins;
  4. digits 27–162 by a limb engine on (2^r mod 10^162)·(2^(jP) mod 10^162);
  5. GMP/MPFR on the last 1000 and first 40 digits.
- **Integrity gates before every `CDONE`:** tables re-hashed, the Q table re-hashed, running state re-derived by
  modular exponentiation, stage 2 against stage 1, and sampled stage 4 against stages 1–2.
- **Builds:** `z2sieve` (K = 13, NL = 18, TB = 24, TD = 7) and the positive-control build `z2sieve_k4` (K = 4, NL = 2,
  TD = 3, NMIN = 60).
- **Tools:** `val_z2.py` (exact check of a dump), `val_z2.sh` (the whole validation), `prod_z2.py` (driver with resume
  and a CAND/MEMBER tripwire), `aggregate_z2.py` (coverage).
- **Prototype check of lemma 2 first:** Python, 20,000 random n < 10^17 for K = 4 and for K = 13, with 0 mismatches.

### P.3 Validation (`validation_pc/pow2/val_z2.out`, source `z2sieve.c` `cdb63fb3…`): ALL CHECKS PASSED
- **K: kernel tests**, 10^6 random n < 2^56 per build, against GMP/MPFR.
  - Every stage had 0 mismatches.
  - K = 13: 450,951 stage-3 rejections were confirmed by rigorous leading digits.
  - K = 4: 209,445 stage-3 rejections, likewise.
- **E: exact regions** (dump against Python integers and python-flint balls; S2 set exact both ways; every decision,
  coverage and counter).
  - K = 13: [1000, 2001000), [P − 10^6, P + 10^6), [10^12, 10^12 + 2·10^6), [10^15 − 999993, 10^15 + 1000007),
    [3.7·10^13, +2·10^6) and [2^61 − 2·10^6, 2^61). All OK, each with about 143,600 stage-2 survivors.
  - K = 4: [10^12, +10^6) and [2^61 − 10^6, 2^61). About 132,000 stage-4 survivors each, every one rejected by
    the final windows and confirmed.
- **C: positive control.** The K = 4 build over [60, 2·10^6) reports MEMBER for exactly 67, 72, 76, 77, 81, 86. Every
  other n has its 0 confirmed independently.
- **T: threads.** Dumps at 1 and 8 threads are identical (K = 13 at 10^15; K = 4 over [60, 2·10^6)).

### P.4 Benchmark (K = 13, 8 threads)
- **Run:** `./z2sieve 1000000000000 1250000000000 threads=8`, 2.5·10^11 exponents in 87.1 s wall (695 CPU-s).
  That is 2.9·10^9 exponents/s, or 2.78 ns per exponent per thread.
- **Survivors:**
  - stage 1: 28.2%;
  - stage 2: 7.17%;
  - stage 3: 3.94%;
  - stage 4: 5,987 (2.4·10^−8, as 0.9^136 predicts);
  - stage 5: 0.
- **Projection at 20 threads:** to be measured, estimated 5–6·10^9/s. So 10^13 in about 0.5 h, 10^14 in about 5 h.
- **Next step: a CUDA port.** Stage 2 is a few additions and three bit lookups in tables that fit in shared memory.
  It should reach 10^15–10^16 in hours, once the GPU is free.

### P.5 CUDA variant (`z2sieve_cuda.cu`, `make cuda`; design in `docs/ZEROLESS_POW2.md`)
- **Design.**
  - One thread per stage-1 survivor, from a host bitmap with ranks and run-start tables every 1024 residues.
  - Lockstep j-walks in blocks of 128 threads, with a shared queue of stage-2 survivors drained every 8 steps.
  - Stage 5 stays on the host. Device tables are compared with their host copies before each `CDONE`.
- **Tuning.**
  - Measured on [10^12, 1.25·10^12), 2.5·10^11 exponents, with JCH = 1024 until the last step.
  - Timing probes: a `-DZ2PROBE` build stops after a given stage. Nsight's counters need a Windows setting that was not
    changed.

  | step | change | wall |
  |---|---|---|
  | 1 | first version | 11.45 s |
  | 2 | per-thread counters instead of shared atomics | 11.45 s (atomics were not the cost) |
  | 3 | probes: setup 0.80 s, +stage 2 2.33 s, +stage 3 2.84 s, +stage 4 11.46 s | the limb engine dominates |
  | 4 | exact multiplier for the 10^9 quotient (m = ⌈2^76/5^9⌉ on acc >> 9, proven, 1.4·10^6 values checked); 32-bit digit windows; self-check sampled 1/64 as on the CPU | 11.45 s (latency-bound) |
  | 5 | limbs of 2^r in L1-cached global memory (higher occupancy), queue drained every 16 steps | 9.69 s |
  | 6 | limbs back in shared memory but limb-major (no strided rows), every 8 steps | 8.16 s |
  | 7 | two passes: stage 3 with compaction; Q rows of the window in shared memory; 16-bit queue entries | 6.02 s |
  | 8 | exact pre-filter of limbs 0–3 with no early exit (whole warps), then the full engine; JCH = 4096 | 3.75 s (10^12 in 15.0 s) |

- **Production speed.** One chunk [4·10^12, 8·10^12) takes 59.1 s, or 6.8·10^10 exponents/s, 23 times the CPU at 8
  threads. The GPU ran at 1920 MHz, 121 W, 68 °C.
  - Survivors: s2 = 2.87·10^11, s3 = 1.57·10^11, s4 = 93,993, s5 = 0.
- **Projection.** 10^14 in about 25 min, 10^15 in about 4.1 h, 10^16 in about 41 h.
- **Safety fix before commit.**
  - At first the probes were switched on by the environment variable `Z2_PROBE` alone.
  - A stray variable would have skipped stages, yet still logged coverage. Probes now exist only in `-DZ2PROBE`
    builds.
  - The production build ignores the variable (checked).

### P.6 Validation of the final binaries (`validation_pc/pow2/val_z2_cuda.out`): ALL CHECKS PASSED (24 of 24)
- **Builds.** Sources `z2sieve.c` `1847c835…`, `z2sieve_cuda.cu` `4cc8b5e1…`, `z2sieve_cuda.h` `07c605c2…` (JCH = 4096).
  - Binaries: `z2sieve` `0fa4226a…`, `z2sieve_k4` `3a57e241…`, `z2sieve_cuda` `7c5cd3b3…`, `z2sieve_k4_cuda` `f0b1fce9…`.
  - GPU: RTX 3060 Ti, driver 610.88.
- **Repeated sections K, E, C, T** (as P.3): identical numbers.
- **G: the GPU builds.**
  - All six exact K = 13 regions are OK, with the same counts as the CPU build.
  - The positive control over [60, 2·10^6) gives MEMBER exactly 67, 72, 76, 77, 81, 86.
  - The dumps at 1 and 8 threads are identical.
- **X: CPU against GPU.**
  - Sorted dumps are identical, both for K = 13 at 10^15 and for K = 4 over [60, 2·10^6).
  - On the production-scale slice [4·10^12, 4.25·10^12) the CDONE counters are identical:
    s1 = 70,512,472,576, s2 = 17,923,428,218, s3 = 9,841,385,431, s4 = 5,935.
- **Record.** `validation_pc/pow2/val_z2.out` remains the record of the first CPU-only version (commit `b65c790`).

## 2026-10-21 · Session 4 (cloud container) · Second pass of the six shells (ROADMAP P2.5)

Goal: the independent CPU second pass of the shells 40,000 < m ≤ 100,000 (NOTEBOOK 4.4), run on cloud server hardware
instead of the desktop. The desktop was unavailable for GPU work, and its memory has no ECC.

### S.1 Machine and build
- **Machine.** Cloud VM (Firecracker): Intel Xeon @ 2.10 GHz, 4 vCPUs, AVX-512 with IFMA, 15 GB RAM. Ubuntu 24.04,
  gcc 13.3.0, GMP 6.3.0 and MPFR 4.2.1 (`apt-get install libgmp-dev libmpfr-dev`), python-flint 0.9.0 (pip).
  - Server hardware of this kind normally has ECC memory, but that cannot be observed from inside the VM (no EDAC, no
    DMI). The same holds for the original m ≤ 40,000 run.
- **Build.** `make` at commit `416796e`: `gsieve.c` `addf4d33…`, `gsieve_k13` `052dc15d…` (variant ifma, `-O3
  -march=native`, NL = 14).

### S.2 Validation on this machine (`validation_pc/task4_check/cloud_*.log`)
- `make -C work/gsieve test`: all passed in 2 min 13 s.
  - Kernel arithmetic: OK.
  - Four `val_dump` regions: every decision verified.
  - E-mode end to end: 314/314 for K = 13 and K = 12.
  - Direct, odd and thread decisions: verified.
  - `val_prodrow` at c = 35,500: OK.
- `make -C work/gsieve test-prod`: MATCH on all counters.
  - main c = [7540, 7898) at N = 36,000: 8.59·10^11 points in 70 s.
  - That is 1.23·10^10 points/s, or 3.1·10^9 per thread.
- **Projection for the second pass:** about 22 h in all.
  - Main paths, in order: 1.4, 2.1, 2.9, 3.8, 4.9 and 6.1 h.
  - Direct and odd paths: about 15 min per shell.

### S.3 Second pass, shells 40,000 → 70,000 (started 2026-10-21 14:01:46 UTC, commit `ac0db01`)
- **Inputs.** The production logs were unpacked from `work/gsieve/logs/prod_N{50000,60000,70000,80000,90000,100000}_n0*.tar.xz`
  into `/root/scratch/shells/`. Every file matches its `.sha256` list.
- **Command.** `validation_pc/task4_check/launch_check.sh "40000:50000 50000:60000 60000:70000" A` (from the repository
  root), which runs `setsid nohup ./run_shells.sh /root/scratch/shells "…" cpu` from `work/gsieve`.
  - It is detached so that it outlives the 2 h limit on background tasks in this environment.
  - For each shell, `val_proddone.py` re-runs each of the three production invocations (main, direct, odd over
    c = [0, N]) with `gsieve_k13` at 4 threads, and requires all eight DONE counters to match.
  - Then `aggregate.py` (full coverage, no members) and `compare_cdone.py` (per-c points, s4 and members) run, and the
    logs are archived to `validation_pc/task4_check/check_N<N>_n0<N0>.tar.xz`.
  - The MEMBER tripwire is active.
- **Scope.** The first three shells take about 7 h. The last three (70,000 → 100,000, about 15 h here) are held for
  the user's approval, because together they exceed the 12 h limit.
- **Results:** in S.4.

### S.4 Results so far
- **Shell 40,000 < m ≤ 50,000: PASSED** (14:01:46–15:27:34 UTC; archive `validation_pc/task4_check/check_N50000_n040000.tar.xz`,
  SHA-256 `ee63d8ec…`).
  - All eight DONE counters match production for each invocation:
    - main: 60,979,301,760,000 points, s1 = 4,190,465,737,806, s2 = 230,373,239,297, s3 = 92,062,786,045, s4 = 14;
    - direct: 23,398,270,000 points, s4 = 171,177;
    - odd: 35,085,180,000 points, s1 = 1,248,050,760, s4 = 36.
  - full = 0 and members = 0 on every path.
  - Coverage is complete (50,001/50,001 c values on all three paths, no members).
  - Per-c points, s4 and members are identical at all 50,001 boundaries on every path.
  - Wall times at 4 threads: main 4,895 s, direct 157 s, odd 89 s. The projection was 1.4 h for main.
- **Shell 50,000 < m ≤ 60,000: PASSED** (15:27:34–17:50:49 UTC; archive `validation_pc/task4_check/check_N60000_n050000.tar.xz`,
  SHA-256 `68c08e15…`).
  - All eight DONE counters match production for each invocation:
    - main: 90,974,701,760,000 points, s1 = 6,251,733,870,317, s2 = 343,694,020,810, s3 = 137,348,477,820, s4 = 42;
    - direct: 28,598,270,000 points, s4 = 209,068;
    - odd: 42,885,180,000 points, s1 = 1,797,660,760, s4 = 55.
  - full = 0 and members = 0 on every path.
  - Coverage is complete (60,001/60,001 c values) with no members, and per-c counts are identical at all 60,001
    boundaries.
  - Main took 8,257 s wall; benchmarks of the sub-cubic prototype shared the CPU part of the time.
- **Shell 60,000 < m ≤ 70,000: stopped by decision of the user at 18:01 UTC** (focus moved to the sieve).
  - Its direct invocation had completed: 33,798,270,000 points, s4 = 247,224, all counters MATCH. The log is kept in
    `validation_pc/task4_check/partial_N70000_n060000/`, and `val_proddone.py --logdir` reuses it on resume.
  - Main had run for 6 minutes. Its partial log was discarded, because `val_proddone.py` re-runs incomplete
    invocations from the start.
  - The driver output up to the stop is in `validation_pc/task4_check/driver_A.out`.
- **Status of P2.5:** 2 of the 6 shells are second-passed (40,000 → 60,000). 60,000 → 100,000 remains: about 17 h
  on this machine, or about 6 h on the desktop CPU.

## 2026-10-21 · Session 4 (cloud container) · Sub-cubic chain sieve: prototype (`work/tsieve/`)

Question from the user: can the N³ cost of the chain sieve be beaten?

### T.1 Idea
- The points (a, b, c) whose last D digits avoid {0, 5} form exactly 4^D cosets of the lattice
  L_D = {a + b·μ3 + c·μ7 ≡ 0 (mod P_D)}, which has index P_D = 4·5^(D−1). This follows from L1 and L2.
- The region has about N³ points. If L_D is well rounded, a coset holds about N³/P_D points.
  - Enumerating cosets instead of points costs about 4^D at the depth where P_D ≈ N³.
  - That is N^(3·log 4/log 5) = N^2.58.
  - LLL check: for K = 13…30, the shortest vector divided by det^(1/3) is 0.34–0.99, so L_K is well rounded.
- Elementary consequence (not yet in RESULTS.md, needs a literature check): the number of chain numbers with
  m ≤ N is ≪ N^(2 + log 4/log 5).

### T.2 Code (prototype, not production)
- **`tparams.py K`.** Constants for depth K:
  - level-2K logarithms;
  - an LLL-reduced basis of L_K in the norm (2a, 3b, 6c), with the walk direction chosen so that its step changes
    J by a unit mod 5;
  - the lifting constants of the tree.
- **`tsieve.c`.** The tree of admissible residues to depth K (4^K leaves, the excluded child found from one digit),
  then every point of each leaf's coset (digits K+1..2K by lemma L3, affine along rows).
  - `-DCHECK` re-derives every leaf state and every digit rule.
- **`tsieve4.c`.** Cosets of depth K, with rows along the walk direction. Within a coset:
  - The lowest M digits of the block h0 + 2^K·Y depend only on Y mod 5^M.
  - The admissible set T has exactly 4^M elements, by the same parity argument as L2. It is built as a tree from h0.
  - Along a row, Y = Ys + n·δ. So the hits are the set bits of U = δ^−1·T in the window [δ^−1·Ys, +len), which is
    read 64 positions per word.
  - Row bounds are exact: inequalities in lattice coordinates, Fourier–Motzkin projections for the outer loops.
- **Validation:** `val_tsieve.py` and `val_tsieve2.py` brute-force every point with exact integers.
  - Regions: N = 60…120, K = 3…6, M = 1…5, including shells with N0.
  - Both versions give an EXACT MATCH on every region:
    - the count of points that pass K + M digits;
    - the exact set of points that pass 2K digits;
    - 4^K leaves, with no duplicates.

### T.3 Measurements
- **Setup.** Single thread, CPU time. The CPU was shared with the 4-thread second pass, and the absolute numbers are
  inflated by about 1.5×. Main path, a ≥ 26.

  | N | points | gsieve (26 digits) | tsieve4, best (K, M) | digits | ratio |
  |---|---|---|---|---|---|
  | 2,000 | 7.9·10⁹ | 5.27 s | 10.5 s (9, 4) | 18 | 2.0 |
  | 4,000 | 6.4·10¹⁰ | 32.5 s | 47.0 s (10, 4) | 20 | 1.45 |
  | 8,000 | 5.1·10¹¹ | 189.8 s | 233 s (11, 4/5) | 22 | 1.23 |

- **Scaling.** gsieve's time per point falls from 0.67 to 0.37 ns over this range: its rows get longer, and its
  asymptote is N³ at about 0.3 ns per point. tsieve4 scales as N^2.2–2.3 here, and its asymptote is N^2.58.
- **Cost breakdown at N = 4,000** (timing builds, noisy):
  - tree nodes about 7 ns each;
  - rows about 20 ns each;
  - each hit about 3.5 ns to iterate, plus 3.5 ns for the digit check.
- **Projection with the present scalar code** (not measured):
  - crossover near N ≈ 1.5–3·10⁴;
  - about 2.5–3× faster than gsieve at 10⁵ and about 6× at 10⁶;
  - more if the hit loop is vectorized (gsieve's stage 2 is AVX-512) or moved to the GPU, since the cosets are
    independent.
- **Not done:**
  - digits beyond 2K (K = 13 gives 26, as gsieve's stage 2 does);
  - the slab 26 ≤ a < 2K;
  - the direct and odd paths;
  - threads;
  - E mode;
  - a production-scale cross-check against gsieve's dumps.

### T.4 `tsieve6`: the full main path, identical to gsieve (session 4, 21 Oct 2026)
- **What it is.** The coset tree with bitmap hit lists (K = 13, so the block covers digits 14–26) feeds gsieve's own
  stages 3–5, included from `gsieve.c`: the TOP table, `limb_check` and `full_check`. M = 3^b 7^c mod 10^126 comes from
  tables (T3, T7).
  - Counters: points and s1 as gsieve counts them. Region points come from gsieve's row loop. For s1, the polytope
    points below 10^XMIN whose last 13 digits pass are subtracted.
  - pthreads over 4^7 subtree roots.
  - Modes G and E. In E, a node has 5 admissible lifts when their digits are odd.
  - Downstream batched per thread: stage 3, then `limb_check`'s first test from the low 4 limbs (the same decision,
    with prefetch), then the full M product and `limb_check` for the rest.
  - `roots=A:B` processes a subset of the 16,384 roots, for benchmarks.
- **Bugs found while building it.**
  - An implicit `exp10l` declaration: `<signal.h>` was included before gsieve.c's `_GNU_SOURCE` (the CLAUDE.md
    pitfall). Fixed by defining `_GNU_SOURCE` first.
  - A sign error in a reciprocal made row bounds slow. They stayed exact.
- **Validation** (`val_tsieve6.sh`: gsieve and tsieve6 on the same region; every DONE counter and the sorted
  S2/R3/R4/S4/MEMBER dump lines):
  - G mode: N = 300, 500, 800; shells (900, n0 = 600) and (1200, n0 = 1100); 2–4 threads. IDENTICAL.
  - E mode: N = 300, 600 and (900, n0 = 700). IDENTICAL. N = 300 and 600 have members=1: the known 137-digit
    2^77·3^64·7^98, found by the tree on the main path.
  - **Positive control** (`val_e2e_tree.py`, NL = 7, N = 160). tsieve6 on the main path (173 members) plus gsieve's
    direct (59) and odd (82) paths gives **314/314** of the known zero-free numbers: none missing, none extra.
  - The E-mode tree has 275,439,346 leaves, which equals gsieve's E-mode bitmap count (`good=`). That is an
    independent check of the tree.
  - Re-validated after every later change (the last three G/E regions and 314/314).
- **Speed** (this machine, 4 threads).
  - gsieve, N = 16,000: 357 s wall (1,263 CPU-s), or 1,405 s with one thread.
  - tsieve6 is slower at that size. With K = 13 there are always 6.7·10⁷ cosets of about N³/P13 points each (520 at
    N = 8,000), and per-coset work dominates.
  - At N = 10⁵, on 1/256 of the cosets (`roots=0:64`, M = 7), extrapolated:

    | version | CPU-s per 1/256 | full run |
    |---|---|---|
    | first | 1,772 | 126 core-h |
    | batched downstream with prefetch | 1,514 | 108 core-h |
    | without stages 3–5 (timing build) | 939 | — |
    | gsieve at 0.326 ns/point | ≈ 1,270 | ≈ 90 core-h |

    With M = 8 (no table, the tree per coset) it took 1,920.
  - Counts per 1/256: s1 2.68·10¹¹, hits 5.6·10¹⁰ (4.8× fewer than s1), rows 5.4·10⁹, s2 1.48·10¹⁰.
  - Costs: about 7 ns per hit, about 90 ns per row (about 10 hits each), and about 39 ns per s2 survivor
    downstream. gsieve's stages 3–5 cost about 41 ns per s2 survivor.
- **Conclusion for the CPU.** The tree sieve is about 1.2× slower than gsieve at N = 10⁵ on this CPU. Its cost grows as
  N^2.58 against N³, so the crossover is near N ≈ 3·10⁵. The constants that matter are per hit and per row, which
  scalar code does not reduce enough.
- **Next step: the GPU.** Per coset, the work is regular integer arithmetic: T(h0) as a tree, U in shared memory
  (2·5^7 bits = 19.5 KB), rows, hits, and the downstream with small tables.
  - Estimate at N = 10⁵: about 10¹⁵ simple instructions, a few minutes on the RTX 3060 Ti. The CUDA gsieve takes
    about 64 min at 2.6·10¹¹/s (G mode with GSIEVE_GPU34).
  - It can only be validated on the desktop. Plan: a CUDA build with the same counters and dumps, and val_tsieve6.sh
    against gsieve as the acceptance test.

### T.5 GPU port (`tsieve_gpu`): written and emulated here, not yet run on a GPU
- **Files.**
  - `work/tsieve/tsieve_gpu.cu`: the kernels. The same source builds for CUDA or, with `-DTG_EMU`, for the CPU.
  - `tsieve_gpu.c`: the host. It includes tsieve6.c with `TSIEVE_NO_MAIN` and a leaf hook.
  - `tsieve_gpu.h`, the Makefile and `val_tsieve_gpu.sh`.
- **Design: one thread block per coset.**
  1. U = δ⁻¹·T(h0) goes into shared memory as a doubled bitmap of 2·5^M bits (M ≤ 7, 19.5 KB). T(h0) is enumerated
     leaf by leaf. In G mode, leaf l takes the (l >> 2i & 3)-th of the 4 admissible lifts at level i. In E mode all
     5^M paths are tried and paths through a 0 are dropped.
  2. For each chunk of up to 1,024 i0 values, the i1 ranges come from the exact projection, followed by prefix sums.
  3. Rows are handed out by index: exact i2 range, then the window of U read 32 positions per word. Each hit's block
     is checked, and gsieve's stage 3 (TOP) and `limb_check` digits 27..9NL run inline, low limbs first, then the full
     product.
  - The host computes each coset's constants and exact i0 range. Empty cosets are skipped.
  - The GPU returns counters and the points that pass stage 4 (for `full_check` on the host). It also returns points
    with b·lg3 + c·lg7 < XMIN + 2, which the host decides with gsieve's region rule. With `dump` it returns every s2
    point with its decision.
- **Emulation.** Blocks run one after another, and every phase between barriers is a loop over the thread ids
  (`TG_EMU_THREADS`, 32 by default). The emulation follows the kernel's control flow, including the strided loops and
  the hand-out of rows.
- **Fail-stop checks for non-ECC memory.** `tg_verify` checks that device tables (T2, T3, T7, TOP) and constants equal
  the host copies, and that the host copies still have their hash from `tg_init`.
  - It runs before the run, every 60 s (`TG_VERIFY`) and at the end.
  - The DONE line is written only after the last check passes. Any failure stops with status 3.
- **Build.** `nvcc 12.0 -arch=sm_86`: 94 registers, a 128-byte stack frame, no spills, about 46 KB of dynamic shared
  memory per block, so 2 blocks of 256 threads per SM. It links. Here it stops at "no CUDA-capable device".
- **Pitfall hit twice.** `<unistd.h>` and `<signal.h>`, included before gsieve.c's `_GNU_SOURCE`, left `exp10l`
  implicitly declared. That would make the TOP table wrong. The tsieve Makefile now builds with
  `-Werror=implicit-function-declaration`, and each host file defines `_GNU_SOURCE` first. gsieve.c is unchanged
  (`addf4d33…`).
- **Validation here (emulated kernel against gsieve).** G N = 300: IDENTICAL (counters and dumps). The other regions
  and 314/314 are in T.6.
- **Acceptance on the desktop.**
  - Build: `make -C work/tsieve cuda tsieve6`, then `./val_tsieve_gpu.sh`. It covers exact regions in G and E,
    314/314, the GPU tree against the CPU tree on a slice of the N = 10⁵ cosets (every counter), and N = 16,000
    against gsieve's DONE line with timings.
  - Expected speed (estimate, not measured): about 10¹⁵ simple instructions for N = 10⁵, a few minutes. gsieve's CUDA
    variant takes about 64–86 min for the same region.

### T.6 Validation of the GPU kernel in emulation (cloud, commit `03acd77`; `validation_pc/tsieve/cloud_emu_validation.txt`)
- **Build.** `tsieve_gpu_emu`: tsieve_gpu.cu built with `-DTG_EMU`, `TG_EMU_THREADS=8`. This is version 2 of the
  kernel: hit offsets from a 33-entry table and u32 leaf arithmetic.
- **Against gsieve** (every DONE counter and the S2/R3/R4/S4/MEMBER dump lines), all IDENTICAL:
  - G, N = 300, M = 4 (this region also exercises the points near 10^XMIN, which are handed to the host);
  - G, shell N = 900 with n0 = 600, M = 3;
  - E, N = 300, M = 3: s4 = 30, members = 1, which exercises the full limb path and the S4 output.
- **Against tsieve6 at production scale** (N = 10⁵, one subtree of 4,096 cosets, 4.19·10⁹ points), every counter
  (s1, s2, s3, s4, gmp, full, members, leaves, rows, hits) MATCHES:
  - roots 0:1, M = 7: 8.79·10⁸ hits, 8.39·10⁷ rows;
  - roots 1:2, M = 6, with a build that forces the i0 chunking (IMAX = 16): 1.10·10⁹ hits.
- **Positive control.** `tsieve_gpu_emu_nl7` (M = 2) on the main path plus gsieve's direct and odd paths:
  **314/314**, none missing, none extra.
- **Not run.** An N = 4·10⁵ slice was stopped for time; its CPU side alone takes more than 16 min here.
- **Earlier attempts, discarded.**
  - An E-mode emulation at N = 600 with M = 5 took more than 2 CPU-hours (2.75·10⁸ cosets × 3,125 paths) and was
    stopped.
  - One shell comparison read DIFFERENT only because a cleanup command killed its gsieve process, so the gsieve
    counters were empty. Re-run: IDENTICAL.
- **Status.** The logic of the GPU kernel is validated. CUDA-specific behaviour (barriers, shared-memory atomics,
  occupancy, speed) still needs the desktop run (`val_tsieve_gpu.sh`).

### T.7 Desktop acceptance, a 3× faster kernel, and the full N = 10⁵ run (session 5, 21–22 Oct 2026)
- **Machine and builds.** The desktop and WSL distribution of D.1 (nvcc 13.3, driver 610.88, RTX 3060 Ti). Branch
  `ccr-33b3eadb-tsaozk` at `f0a6868`: `make -C work/gsieve` had nothing to rebuild (`gsieve_k13` = `46f8308b…`, as in Z.2);
  `make -C work/tsieve cuda tsieve6` built `tsieve_gpu` `f3abeb5b…` and `tsieve6` `4f9c35a8…`.
- **Acceptance of kernel v2** (`./val_tsieve_gpu.sh`, `validation_pc/tsieve/gpu_val/`, 21:49–22:15 UTC): **ALL PASSED**, no
  CUDA-specific failure.
  - A: all 8 exact regions IDENTICAL to gsieve (counters and S2/R3/R4/S4/MEMBER dumps): G N = 300, 500, 2000 and shells
    (900, n0 = 600), (1200, n0 = 1100); E N = 300, 600 (each with the member 2^77·3^64·7^98) and (900, n0 = 700).
  - B: 314/314. C: N = 10⁵, roots 0:16, M = 6 and 7: every counter MATCHes the CPU tree (GPU 8.04 s, CPU tree 30.9 / 29.4 s at
    20 threads). D: N = 16,000 MATCHes gsieve (tree 130.2 s, gsieve avx2 at 20 threads 126.9 s).
- **Where the time went.** Nsight Compute needs the Windows GPU-counter permission (not changed), so timing-only probe builds
  (kernels stopped after each phase; never used for counts) and host experiments were used.
  - The host is not the bottleneck: 4, 10 and 20 host threads give the same wall time for 20 roots, the GPU is at 100%, and the
    host's CPU time is spin-waiting in `cudaStreamSynchronize`.
  - Per 20 roots at N = 10⁵, M = 7: rows, i2 ranges and window scan 0.5 s; with the per-hit block check 1.5 s; with the call of
    `later()` (stages 3–4) 3.5 s; with the TOP lookup and the first limb test 3.7 s; complete 8.05 s. The eager 14-limb product
    after the first limb test, done in local memory by the one to three lanes of a warp that had such a point, took more than
    half. (An FP64 corner test looked costly in a probe with 1 s resolution; removing it changed little.)
- **Kernel v3 (commit `560ddbf`), decisions unchanged:**
  1. Candidate queue: phase 3 runs in rounds of one row per thread; points that pass digits 14–26 are queued in shared memory
     (2,048 entries, a·2^42 + b·2^21 + c), and all threads then drain the queue through `later()`. A full queue falls back to the
     inline call.
  2. Lazy limbs: after the first test, limbs 4 … NL−1 of 3^b 7^c and of 2^a·M are computed in order from that test's carries and
     the loop stops at the first bad limb, which is the same conjunction `limb_check` evaluates.
  3. (a, b, c) of a hit from 32-bit row bases (exact, since a ≤ 3N, b ≤ 2N and c ≤ N are below 2^32 for N ≤ 10⁶; checked).
  4. Points near 10^XMIN go to the host by the integer superset b < 278 and c < 157, which contains every point with
     b·lg3 + c·lg7 < XMIN + 2; the host applies gsieve's region rule and stages 3–5 as before, so no FP64 remains on the GPU.
  5. `floordiv` by a multiply-high with a precomputed reciprocal and an exact correction; `mods_br` without a 64-bit `%`; the
     verifier sleeps in 0.1 s ticks (the final join had waited up to 1 s, which is why every earlier wall time ended in .04).
  6. `TG_NT` defaults to 512 for M = 7 (one block per SM fits in shared memory, so 16 warps) and to 256 otherwise.
  7. Test builds with a tiny queue and 16-value i0 chunks, `tsieve_gpu_emu_q8` (in `make test`) and `tsieve_gpu_q64` (new
     section A2 of `val_tsieve_gpu.sh`), make the inline path and the chunking run often.
- **Tuning** (wall s for roots 0:20 of 16,384 at N = 10⁵; kernel v3 unless marked):

  | | TG_NT 128 | 256 | 512 | 1024 |
  |---|---|---|---|---|
  | M = 5 | — | 3.44 | 3.54 | — |
  | M = 6 | — | 3.24 | 3.44 | — |
  | M = 7 | — | 4.95 | **3.04** | does not launch (96 registers × 1024 threads) |
  | M = 7, kernel v2 | 15.0 | 10.0 | 10.0 | — |

  QCAP (M = 7 / 512, M = 6 / 256): 1,024: 4.34 / 3.34; 2,048: 3.04 / 3.24; 4,096: 3.54 / 5.25 (shared memory then limits the
  blocks per SM). Chosen: M = 7, TG_NT = 512, QCAP = 2,048. The host now rejects an unlaunchable TG_NT with a clear message.
- **Re-validation of v3** (`validation_pc/tsieve/v2/`): `make -C work/tsieve test` passed (12 min; the CPU tree and the emulated
  kernel, including the tiny-queue build, against gsieve; 314/314 twice). `val_tsieve_gpu.sh`: **ALL PASSED**, binary
  `31f95f00…` (A 8/8 and A2 3/3 IDENTICAL; B 314/314; C MATCH with the GPU at 2.64 s against 28.3 / 26.9 s for the CPU tree;
  D MATCH, tree 98.8 s against gsieve 115.6 s).
- **Full run, N = 10⁵ (the request of task 1).** Command (in `work/tsieve`): `./tsieve_gpu tp13.txt G 7 100000 threads=20`, commit
  `560ddbf`, binary `31f95f00…` (nvcc 13.3 `-O3 -arch=sm_86`, host gcc 13.3 `-O3 -march=native`), defaults TG_NT = 512,
  TG_BATCH = 512, QCAP = 2,048; 23:55:32–00:41:15 UTC; log `validation_pc/tsieve/v2/full_N100000.out`.
  - **Wall 2,725.8 s (45 min 26 s)** for the sieve, 45 min 42 s with the table set-up; 1.07 GB host memory; GPU at 1,920 MHz,
    71 °C, 163–193 W; 45 table checks passed (one per minute and at the end).
  - DONE: points = 999,770,015,045,376, s1 = 68,703,672,291,929, s2 = 3,777,027,389,141, s3 = 1,509,389,143,313, s4 = 340,
    gmp = 340, full = 0, members = 0; 4^13 leaves, 1.375·10¹² rows, 1.441·10¹³ hits.
  - **Independent second pass of the main path for m ≤ 100,000.** `work/tsieve/main_totals.py` adds up the logged production
    DONE counters of the main path (m ≤ 36,000 in two invocations, the killed range c < 7,540 from the session-3 re-run, checked
    against the production CDONE line; 36,000 < m ≤ 40,000; the six shell archives, each checked against its `.sha256`). Every
    counter MATCHes (`validation_pc/tsieve/v2/full_N100000_vs_production.txt`). The tree sieve enumerates the points by coset
    instead of by bitmap scan and affine lift, and runs stages 3–4 in its own GPU code, so stages 1–4 of the whole main path
    behind claim 1.1 now have two independent executions (the direct and odd paths are not covered by it).
  - **Speed:** gsieve's CUDA variant took 4,429 s for the main paths of the six shells (their DONE lines) and about 5 min for
    m ≤ 40,000, about 79 min in all; the tree sieve takes 45.4 min, **1.74× faster at N = 10⁵**, and its cost grows as N^2.58
    against N³. Estimate (not measured): the shell 10⁵ < m ≤ 2·10⁵ in about 3.7 h on this GPU.
- **Status.** Task 1 done: the GPU tree sieve passes its acceptance test, is 3× faster than the version written in session 4,
  and has re-executed the main path for m ≤ 10⁵ with every counter equal. `docs/RESULTS.md` 1.1 gains this as evidence (main
  path); its status stays COMPUTATIONAL.
- **Pitfall hit:** a file rewritten with the Bash tool's Windows Python came out with CRLF line endings (Makefile,
  `val_tsieve_gpu.sh`, `tsieve_gpu.c`); converted back before any build. Repository files are now edited only from WSL or the
  editor (also in the session memory).

## 2026-10-22 · Session 5 (desktop PC) · Even-target tree nodes on the GPU (task 2)

Design notes: `docs/GPU_EVEN_TARGETS.md`. No result in `docs/RESULTS.md` changes in this entry.

### E.1 Cost of the open nodes (`work/node_costs.py`, `work/logs/node_costs.out`)
- For every bag of `bags.bags_of(node)`, `prove_cum.prove_bag(bag, dry=True)` replays the bootstrapping plan and adds up its
  cost model (one unit per top case, `BOTTOM_COST = 4` per bottom solve). The split by side is recovered from the facts before
  and after each step and is checked to add up to `res['stats']['cost']`.
- Units, cheapest first: 131712 5.66·10⁹ (14 of 15 bags already proved), 122472 7.88·10⁹ (23 of 26), 24192 8.22·10⁹,
  1229312 6.70·10¹⁰, 314928 1.33·10¹¹, 41472 1.40·10¹¹, 93312 1.72·10¹¹, 221184 2.23·10¹³, 1161216 1.39·10¹⁴,
  11239424 1.32·10¹⁴ + 2 bags, 33718272 1.58·10¹⁵ + 6 bags. The bottom side is 64–74% of the units of every node.
- "+ n bags": for k ≥ 18, `prove_cum.matveev_cum` cannot form the a-priori Matveev bound (its iteration starts at 10³⁰⁰ and
  the bound leaves the float range; `assert nb <= B`). The cost of those bags is unknown until the bound is computed in log
  space. The RESULTS figure "1161216 ≈ 1.4·10¹⁴ units" is reproduced.

### E.2 Profile of topc3 and bottomc3
- `work/bag_plan.py` prints a bag's plan; `work/prof_steps.py BAG N --keep DIR` runs one slice of each of the N costliest steps
  exactly as ctop3/cbottom4 would (one ordering, one group) and saves the job files. Representative bag: 2·2·2·2·2·2·2·3·3·3·7
  (k = 11, the largest bag of 24192).
- Cycle counters in instrumented copies of top3.c and bottom3.c (scratch only; the machine was loaded by the GPU acceptance run
  at the time, so the absolute figures are pessimistic):
  - **topc3, per case:** MPFR logarithm 2,300–5,200 cycles (about two thirds), exact level-1 test 460–2,200, level 2 36–490,
    the prefix 105–2,150 (amortised). Level-1 survival 1.2–2.4%, level-2 0.013–0.05%.
  - **bottomc3, per solve:** about 12,000 cycles: `low(Flow)` 5,400–6,400 amortised (two p-adic logs and two inverses per lower
    configuration, at precision Mmax = 4G + 8 although the solves need only 5^G), the 5-adic series about 2,000–2,500, the
    lattice step 2,150–2,400, mt + β (CRT) 1,050–1,150, the 2-adic part 800–950.
- **Which loop dominates:** the bottom solve, by a little: in time 52–63% of every node (3,000 cycles per bottom unit against
  3,700–8,000 per top case). Both sides have to move for a large gain.

### E.3 The bottom solver on the GPU (`work/csrc/bottomgpu.cu`, `bottomgpu_host.c`, `bottomgpu.h`, `bigfix.h`)
- Drop-in replacement for bottomc3 in `cum` mode (same job file, same CAND/PROG/DONE output). One thread per lower
  configuration: low(Flow) once, then every height of the top digit. Every quantity is the residue bottom3.c computes, with the
  p-adic series rewritten so that no exact division is needed (Σ (−1)^{n+1} p^{nv−e} yⁿ/n′, Horner with tabulated
  coefficients), Newton inverses, Garner's CRT, Montgomery products with one R for 5^{G+s}, 5^G and 5^{G−1}, and a certified
  fixed-point floor for the lattice coordinates (ambiguous within 2^−60 → bottom3.c decides). Rare cases go back to bottom3.c's
  own code (included unchanged); jobs with M varying, G > 200 or another mode go to bottomc3 (exec).
- `work/bottom_ref.py`: Python transcription of bottom3.c next to the restructured form; self-check: the series, plog and β
  agree on random inputs (13,063 (Flow, p) cases at G = 20, 37, 61) and 891/891 planted (t, u′, w) are recovered.
- **Validation of the CPU emulation build** (`work/val_bottomgpu.py`; the GPU build comes next):
  - K: kernel tests at production scale, G = 14 … 200 (all five size classes; boxes up to 10⁹⁴): 2,600 queries, every planted
    point found (1,949/1,949), candidate lists equal to bottom_ref's exact lists, no deferrals; the direct 5-adic (p = 0) and
    2-adic paths are exercised.
  - C: 57 whole steps and slices of the plans of five bags (G = 51 … 383): byte-identical to bottomc3 (47 through the kernel,
    10 handed to bottomc3 for G > 200); `bottom_ref.py job` also reproduces bottomc3 on a slice.
  - P: 60 configurations from the digits of genuine zero-free numbers (small G, so that the configuration space contains the
    genuine digits): byte-identical to bottomc3, the genuine solution found every time. (A first version found no such
    numbers: zero-free numbers have runs of at most five 1s.)
- CUDA build: 255 registers per thread with 40 B to 1.2 KB of spills depending on the size class; speed to be measured.

### E.4 The top filter on the GPU (`work/csrc/topgpu.cu`, `topgpu_host.c`, `topgpu.h`)
- **Observation.** top3.c's decisions depend on m only through m·θ and m·θ₂, and θ = qa/qb, θ₂ = a2/D are tiny (2^−79 at
  G = 51, 2^−197 at G = 123, 2^−318 at G = 197), so m·θ = −Θ·λ + δθ with Θ = C·θ, λ = log(1 + z), |δ| ≤ 1/2 + 2^−50.
  `work/top_filter_ref.py` checks this on 45,000 real cases (G = 51, 123, 197): dist(mθ, ℤ) and dist(ψ₁, ℤ) agree to the double
  precision of the check, and a margin-guarded filter passes exactly the cases topc3 keeps at level 2 (2, 4, 8 of 15,000), with
  no unsound kill.
- **Design.** The kernel computes ψ₁ = frac(−Θλ), ψ₂ = frac(−Θ₂λ) in fixed point (FL ≥ bits(Θ) + 80 fraction bits; the prefix
  logarithms down the prefix tree from host tables, Horner series with intermediates in (0, 1], only the low limbs of Θ·λ) and
  classifies each case: certainly killed at level 1; certain level-1 pass and certainly killed at level 2 (counted in f1); or
  sent to the host, which runs top3.c's own per-case code (copied unchanged in substance: same precision, same operations,
  same run merging). Output and counters are therefore topc3's.
- **First check (CPU emulation build):** on the four saved jobs (G = 51, 87, 123, 197; 989,095 cases) the output is
  byte-identical to topc3: SURV runs and DONE cases, f1, f2, nsurv. The emulated kernel is already 2–4× faster than topc3 on
  one core (fixed-point series instead of MPFR logarithms).
- Drivers: `work/cbottom5.py` and `work/ctop4.py` run a whole step in one process (resumable at ordering boundaries, the same
  checkpoint format as cbottom4/ctop3); `work/prove_cum_gpu.py` runs prove_cum's prove_bag unchanged with these modules and
  proves the bags of a node in parallel; `work/val_prove_gpu.sh` proves a node from scratch with the CPU and the GPU paths and
  compares every bag's facts, bound and counters and the proofC file.

### E.5 Validation and speed on the GPU (`validation_pc/even_gpu/val_gpu.out`)
- **Builds:** `make -C work cuda` (nvcc 13.3 `-O3 -arch=sm_86`): `csrc/bottomgpu` `12de58fd…`, `csrc/topgpu` `c0117bff…`; reference
  binaries `csrc/bottomc3` `f47187db…`, `csrc/topc3` `af706775…`. The top filter's largest size class is L = 22 limbs (Θ up to
  about 2^592, G up to about 180): nvcc 13.3 segfaults while generating PTX for L = 24 (cicc; also with ptxas -O1/-O2), so steps
  with a larger Θ (G = 197 and up) are handed to topc3.
- **bottomgpu:** K: 3,900 production-scale queries (G = 14 … 200), every planted point found (2,927/2,927), candidate lists equal
  to bottom_ref's, no deferrals. C: the 57 jobs, byte-identical to bottomc3. P: 60/60 genuine configurations, identical and
  found. **VALIDATION PASSED.**
- **topgpu:** K (accuracy, dump mode against arb at prec + 128 bits): 24 jobs at G = 51 … 160, worst error 2 units of 2^−64
  (the margins assume 4). C: 79 whole steps and slices of six bags' plans (k up to 13), byte-identical to topc3 (SURV runs and
  DONE cases, f1, f2, nsurv). P: 150/150 planted genuine configurations kept (planted_top_cum3.py's construction, G = 16), outputs
  identical. **VALIDATION PASSED.**
- **The saved slices** of E.2 (four top, four bottom jobs) are byte-identical on the GPU too.
- **Speed** (`work/bench_gpu_steps.py`, whole steps of the representative bag's plan, wall time including start-up):

  | step | work | GPU | rate | CPU, one core (idle machine) |
  |---|---|---|---|---|
  | bottom j=6, G=87 | 4.19·10⁸ solves | 1.6 s | 2.6·10⁸/s | 3.0 µs per solve (bottomc3) → 1,250 s |
  | bottom j=7, G=160 | 6.12·10⁷ | 1.0 s | 5.8·10⁷/s (start-up bound) | 4.0 µs → 245 s |
  | top j=4, G=123 | 4.13·10⁸ cases | 6.3 s | 6.6·10⁷/s | 1.5 µs per case (topc3) → 620 s |
  | top j=4, G=87 | 9.98·10⁷ | 2.1 s | 4.7·10⁷/s | ≈1.3 µs → 130 s |
  | top j=5, G=51 | 2.14·10⁸ | 12.8 s | 1.7·10⁷/s | 1.7 µs → 360 s |
  | top j=3, G=197 | 2.47·10⁷ | 64.8 s | (handed to topc3) | |

  Against one core: about 800× (bottom, G = 87), 100× (top, G = 123), 30× (top, G = 51: prefixes with heights below 5 go to the
  host, and at G = 51 about 3% of prefixes do). The idle-machine CPU rates are about half the loaded-machine profile of E.2.
- **Next:** `work/val_prove_gpu.sh` on a proved node (both paths from scratch, identical facts and proofC), then the three
  cheap open nodes.

### E.6 Four open nodes proved: 131712, 122472, 24192, 1229312 (22 Oct 2026)
- **Command** (in `work/`): `python3 prove_cum_gpu.py --jobs 12 131712 122472 24192` at commit `5c34c2b` (bottomgpu
  `12de58fd…`, topgpu `c0117bff…`: top size classes up to L = 22, so every step with Θ > 2^592 (G ≥ 183) went to topc3;
  bottomc3 `f47187db…`, topc3 `af706775…`); log `work/logs/prove_gpu_131712_122472_24192.log`, 00:58–01:03 UTC.
  - prove_cum's prove_bag unchanged, with ctop4/cbottom5 as its step modules; 12 bags at a time. The two bags left half done
    by the paused CPU runs of session 2 resumed from their checkpoints in `ckptC/` (exact resume, E.4).
  - **131712** (target 8, depth 2): 58 s; **122472** (target 6, depth 3): 135 s. Then 24192 stopped: one worker's atomic
    checkpoint rename (`os.replace` in prove_cum._save) failed with EACCES on the Windows file system. A scanner or indexer
    holding the file is the likely cause; 12 workers rewrite checkpoints often. Nothing was lost (checkpoints are written
    atomically; the pool stopped; the bags resumed). `prove_cum_gpu.py` now retries the rename for up to 30 s (commit
    `cb16cf5`); prove_cum.py is unchanged.
  - **24192** (target 6, depth 3): resumed with `--jobs 12 24192` (log `work/logs/prove_gpu_24192.log`), 126 s.
- **1229312** (target 2, depth 3; k up to 13): `python3 prove_cum_gpu.py --jobs 12 1229312 314928 41472 93312` at commit
  `48fd1c7` (bottomgpu `12de58fd…`, topgpu `8204424a…`, classes up to L = 32), log `work/logs/prove_gpu_k13.log`,
  01:18–01:47 UTC: 1229312 complete in 1,767 s. The run then stopped early in 314928: one topgpu process died without a
  message. The likeliest cause is the 7.9 GB WSL VM running out of memory, because an nvcc build of the next bottomgpu (cicc
  peak 3.7 GB) ran beside the 12 workers. A failed step records nothing, and the bags resume from `ckptC/` (E.7).
- **Results** (records `work/proofC_131712.json`, `proofC_122472.json`, `proofC_24192.json`, `proofC_1229312.json`; every bag in
  `work/bagresC/`):

  | node | target, depth | bags (max k) | length bound | top cases | bottom solves | complete 7-smooth preimage set |
  |---|---|---|---|---|---|---|
  | 131712 | 8, 2 | 15 (11) | 101 | 1.45·10⁹ | 1.04·10⁹ | 112717121716224 |
  | 122472 | 6, 3 | 26 (11) | 105 | 2.06·10⁹ | 1.46·10⁹ | none |
  | 24192 | 6, 3 | 39 (11) | 105 | 2.13·10⁹ | 1.52·10⁹ | 642978, 2667168, 9437184, 13716864, 1162261467, 3281116734 |
  | 1229312 | 2, 3 | 12 (13) | 138 | 2.31·10¹⁰ | 1.10·10¹⁰ | none |

  (Counts include the steps done by the CPU runs of session 2.) Every bag ends with L ≤ 138 < 5000, so the exhaustive list
  `Z.json` gives the complete preimage sets; each listed preimage has the node as its digit product.
- **Independent check:** two bags proved here entirely through the GPU path were proved again from scratch on the CPU path
  (`work/recheck_bag_cpu.py`: cbottom4/bottomc3, ctop3/topc3, empty scratch directories). Logs in `validation_pc/even_gpu/`.
  - 2·3·3·3·3·3·3·3·4·7 (122472): 642 s on the CPU against 21 s; every fact, the bound and every counter IDENTICAL (top cases
    209,293,913, bottom solves 123,163,479, survivor runs 45,333, 55 proof steps).
  - 2·2·2·2·2·3·3·3·4·7 (24192, proved in two GPU runs around the EACCES stop): 1,263 s on the CPU; every fact, the bound and
    every counter IDENTICAL (top cases 347,218,433, bottom solves 256,880,580, survivor runs 85,584, 55 proof steps).
- **A bug caught before use** (commit `48fd1c7`): the next build of topgpu (classes up to L = 32, so that G ≤ 271 runs on the
  GPU) failed val_topgpu.py part C on two jobs at G = 271. The prefix gaps were 8-bit, so a gap of 260 wrapped to 4 and a
  kept case was killed. Fixed (16-bit gaps). Part K now also checks that the dumped cases are exactly the job's kernel cases;
  it had recomputed ψ for the wrapped case and so could not see the error. The proofs of 131712, 122472 and 24192 used the
  earlier build, which handed every step with Θ > 2^592 to topc3. Such steps include every step whose gaps can reach 256, so
  those proofs are unaffected. The fixed build `8204424a` (used for 1229312) passes K (32 jobs, G ≤ 271), C (79 jobs) and P
  (150/150), and val_prove_gpu.sh 49392.
- Nothing in `docs/RESULTS.md` changes yet: 1229312 still needs its CPU recheck, and the docs and tables are updated once the
  remaining k = 13 nodes are done.

### E.7 Where the time went, and the v2 solvers (commit `3951485`)
- **Measurements** (whole top steps of the largest bag of 1229312, 2⁹·7⁴; `TP_TIME=1` prints a breakdown on stderr):
  - The kernel was not the bottleneck. On the step j = 4, G = 160 (4.1·10⁸ cases) it took 1.4 s, the host's exact code 5 s for
    357,379 level-2 survivor runs (about 14 µs each), and the Python side about 110 s: ctop2.candidates runs an exact
    Fincke–Pohst (fractions) on every survivor run, about 310 µs at G = 160 and 113 µs at G = 51. Every run has its own m, so
    its cache never hits, and 0 of 5,500 sampled runs had a non-empty enumeration. Records = f2 = nsurv on every step: the
    kernel's level-2 test already sends exactly the true level-2 survivors.
  - At small G the host's MPFR code took every prefix with a digit below height 5 (3–4% of the cases at G = 51).
  - Bottom steps with G > 200 ran on one core (bottomc3). For node 221184 that is 6·10¹⁰ solves at G = 234, 1.4·10¹⁰ at 271,
    1.2·10⁹ at 308 and 2.3·10⁸ at 346, about 100 core-hours, with single steps of up to 40 h.
  - prove_cum's planner weight (BOTTOM_COST) does not matter: dry-run plans of 221184's six largest bags priced at the GPU
    rates give 19.3 h for weights 4 (current), 2, 1 and 0.5, and 30.7 h for 0.25. Left at 4.
- **Changes:**
  - *topgpu, small-height table:* log A and 1/A for every configuration of digits at heights 1–4 (8⁴ entries, rounded to
    nearest like the level-1 tables); no host fallback remains. 3.3× at G = 87 and 13.8× at G = 51 on whole steps.
  - *topgpu, recurrence:* along a prefix the cases have consecutive heights, so T_n = |Θ| uⁿ/n (96 fraction bits) is carried
    by T_n(e + 1) = ⌊T_n(e)/10ⁿ⌋, an O(L) step instead of O(N·L²) products; level 2 (after a certain level-1 pass, 1–2% of
    cases) uses the direct λ. The error argument is in the topgpu.cu header. Alone it gives about 5%, as the host dominated.
  - *topgpu, survivor runs resolved in C:* ctop4 appends the lattice basis and R² (numerator and denominator: R² is not an
    integer) to the job. A survivor run whose exact enumeration around (0, …, 0, −m) is empty is counted in a fifth PROG field
    instead of printed, and ctop4 adds the counts to stats['enum'], so every record keeps the CPU path's statistics. The
    enumeration is close_points in integers with fixed common denominators (topgpu_host.c, fp_setup). Its only large operations
    are r_i = (−m a_i) mod D. topc3 never reads the section; a run with any lattice vector goes to Python as before.
  - *bottomgpu:* size classes for G ≤ 240, 280, 320, 360 (16-bit positions), and a per-job check that the fixed widths cannot
    overflow. The check bounds |K| + Q|R|, the coordinates ci, cj and up, ww; a failure sends the job to bottomc3. Before, the
    sizes relied on the class sizing alone.
  - *prove_cum_gpu.py:* one pool for the bags of all given nodes, the costliest first (estimates from
    `logs/node_costs.json`), so that a node no longer ends with one bag running alone.
- **Validation** (`validation_pc/even_gpu/*_v2_*`):
  - val_topgpu.py on the installed binary `5ac4892e…`. K: 48 jobs, worst error 2 units of 2⁻⁶⁴ on every path (table, direct,
    recurrence start and step), with samples of the small-height cases. C: 103 jobs byte-identical to topc3, now including the
    first group of first gaps (the small heights). P: 150/150 planted genuine configurations.
  - New part F. F1: the C enumeration agrees with Python's close_points on 3,400/3,400 values. These are random m, and m next to
    genuine lattice vectors of the box at distances inside, on and just beyond R, in dimensions 4 and 3; 1,530 are non-empty.
    F2: on 97 jobs, DONE and PROG equal topc3's, and the printed runs are exactly topc3's runs with a non-empty Python
    enumeration (9,803 runs, 9,355 resolved), with the per-ordering counts matching. P with the section: 150/150 genuine runs
    still printed.
  - val_bottomgpu.py, emulation and GPU builds. K: planted queries at 21 sizes, G = 14 … 360, all found, 0 mismatches, 0
    deferrals. C: 87 jobs byte-identical to bottomc3, including the G = 234 … 346 steps of 221184's largest bag through the
    new classes. P: 60/60.
- **Speed** (installed binaries, the same steps; "before" rebuilt from the measured per-case and per-survivor costs):

  | top step | cases | before | now |
  |---|---|---|---|
  | j = 6, G = 51 | 4.9·10⁸ | ≈ 30 s | 2.2 s |
  | j = 5, G = 87 | 4.9·10⁸ | ≈ 30 s | 3.4 s |
  | j = 4, G = 160 | 4.1·10⁸ | ≈ 115 s | 11.0 s |
  | j = 4, G = 197 | 4.8·10⁸ | ≈ 85 s | 9.3 s |

  Kernel rates are now 2.4–5.7·10⁸ cases/s. The host part (exact m and the enumeration, about 25 µs per survivor run) runs in
  each worker process, in parallel. Estimate for node 221184 (2.2·10¹³ units): about 5 h of top kernel and 5 h of bottom
  kernel on the shared GPU, so roughly 10–12 h, against about a week for the pipeline of E.6 and about 6,400 core-hours on the
  CPU. To be checked against the k = 13 nodes.
- **Pitfalls:** never run a large nvcc build beside production in the WSL VM: memory is 7.9 GB, cicc for bottomgpu peaks at
  3.7 GB, and 12 workers with CUDA contexts use most of the rest. A shell-quoted commit message with an apostrophe committed
  a truncated message (fixed by amending before anything used the hash); messages now go through a file.

### E.8 314928, 41472 and 93312 proved with v2; independent checks; claim 1.1's second pass complete (22 Oct 2026)
- **Run:** `python3 prove_cum_gpu.py --jobs 10 314928 41472 93312` at commit `3951485` (bottomgpu `6d277763…`, topgpu `5ac4892e…`,
  built by `make -C work cuda`; bottomc3 `f47187db…`, topc3 `af706775…`), log `work/logs/prove_gpu_k13b.log`, 07:05:34–07:59:29 UTC
  (53:55 wall), all 219 bags in one pool, the costliest first. Exit status 0; every bag has its `bagresC` record and each proof
  file lists exactly its node's bags.
- **Results:** 314928 (target 2): no 7-smooth preimage; 93312 (target 2): 116169984, 1811939328; 41472 (target 6): 1843968, 9483264,
  46294416, 69441624, 226492416. Length bounds L ≤ 138 for every bag. With 1229312 (E.6), each tree now has a single open node
  at its first incomplete level (11239424, 1161216, 221184, 33718272); the levels in RESULTS 2.2 are unchanged.
- **Timing, and two disturbances.** For the first 30 minutes the machine also ran the direct/odd sieve re-runs (6 threads), the
  installed-bottomgpu validation, a CPU recheck and a dry-run cost analysis (load about 19 of 20 threads). An attempt to pause
  production for a clean benchmark (SIGSTOP on the job processes, then SIGCONT) was interrupted after the STOP, so production,
  the validation and the recheck stood frozen from 07:24 to about 07:33 UTC. Nothing was lost: the processes resumed where they
  stopped. The side jobs were then ended, and production ran alone. The k = 12 bags took about 1,100 s each under that load, the
  three k = 13 bags 2,843–3,049 s (including the freeze), while the k = 11 bags that ran alone took 150–270 s. **So per-bag times
  measured under contention say little; benchmark with nothing else running.**
- **Independent checks** (`validation_pc/even_gpu/checks_after_k13.out`, 08:04–08:46 UTC, nothing else running):
  - val_bottomgpu.py K C P on the installed binary `6d277763…`: VALIDATION PASSED (`val_bottomgpu_v2_installed.out`).
  - One bag of each node re-proved from scratch on the CPU path (`work/recheck_bag_cpu.py`; bags chosen by estimated cost near
    1.2·10⁹ units): 2⁴·4·7⁴·8 (1229312, 1,390 s), 2·3⁹·8 (314928, 1,271 s), 2⁶·3²·8·9 (41472, 1,228 s), 2⁶·3·6·9² (93312, 1,230 s).
    All four IDENTICAL: facts, bound and every counter, including stats['enum'], which on the GPU path counts the survivor runs
    topgpu resolved in C, so the fp bookkeeping matches the CPU path on production bags.
- **Second execution of the direct and odd paths, shells 60,000–100,000** (claim 1.1). `val_proddone.py ./gsieve_k13 6` (avx2,
  `46f8308b…`, commit `3951485`) re-ran the seven missing invocations; direct 60–70k was reused from S.4. All eight DONE lines match
  production. `compare_cdone.py N N0 PROD RERUN --modes direct,odd` (the new option compares only some paths) gives identical
  per-c counts at every boundary in all four shells. Records: `validation_pc/task4_check/do2_val_proddone.out`,
  `do2_compare_cdone.out`, `check_do_N70000-100000.tar.xz` (`.sha256`). With the tree sieve's re-run of the main path (T.7), every
  path for m ≤ 100,000 has now been executed twice. RESULTS 1.1 and 2.1/2.3 updated; `work/certificate.json`,
  `work/results.md` (= `docs/results_base10.md`) and `MANIFEST.sha256` regenerated (99 nodes).

### E.9 Research: two-sided bootstrapping (implemented) and meet-in-the-middle steps (cost model) (22 Oct 2026)
- **Pool throughput** (`work/pool_bench.py`, log `work/logs/pool_bench_221184.out`; idle machine; cost-weighted sample of
  221184's step slices):
  - 10 workers: 9.5·10⁷ units/s; 6 workers: 8.0·10⁷ units/s. The 4-worker run was stopped (slices too large, over an hour).
  - At 10 workers 221184 would have taken about 65 h with Matveev's chain alone.
- **Two chains** (`docs/BOOTSTRAP_TWO_CHAINS.md`). A 5-adic Baker chain bounds each bottom digit's height from the next lower
  one (Yu, Acta Math. 211 (2013), Theorem 2, p = 5, three logarithms 3, 7, 2^t/F). The top chain uses four logarithms
  (3, 7, 10, N_j/2^t, Matveev n = 4).
  - The constants were read on the typeset pages. The PDF's text layer drops ϰ: C₂ = 1.053·10⁸, not the 2.1·10⁷ first
    transcribed. Per level the bottom chain gains about 10^11.5 and the top chain about 10^16.5.
  - Cost model, `PC_BOOTSTRAP=twochain python3 node_costs.py --json logs/node_costs_twochain.json …`, log
    `work/logs/node_costs_twochain.out`. Times at 9.5·10⁷ units/s:
    - 221184: 8.03·10¹¹ units (2.23·10¹³ before, 28×), ≈ 2.3 h;
    - 1161216: 6.27·10¹² (1.39·10¹⁴ before, 22×), ≈ 18 h;
    - 11239424: 7.05·10¹² (before: 1.32·10¹⁴ plus two bags without estimate), ≈ 21 h;
    - 33718272: 1.64·10¹⁴ (before: 1.58·10¹⁵ plus six bags without estimate), ≈ 20 days.
  - The bags with k ≥ 18 had no estimate (a-priori bound beyond the float range); now every bag can be planned.
- **Implementation.** `prove_cum.twochain_cum` with `PC_BOOTSTRAP=twochain`; the default is unchanged.
  - Records gain `"bootstrap": "twochain"` (bagresC) and a fifth element `"twochain"` per bag entry of proofC.
  - `prove_cum_gpu.py` writes the same entries.
  - Facts are independent of the bootstrap, so checkpoints of either kind can be mixed.
- **Validation.**
  - `work/val_twochain.py`: ALL OK. Constants at 60 digits. Default bound equal to the original code on 3,000 random
    states. Two chains never weaker than Matveev's chain. Exact checks on 1,003 genuine numbers (k ≤ 20): 10,040 bottom-chain
    inequalities, every top inequality, and the bound with every partial set of their true facts.
  - The first attempt died silently: `twoadic.tmax` on 20-digit bags used about 5 GB, so WSL ran out of memory. The check now
    uses each number's own t.
  - `work/val_twochain_e2e.sh` re-proves a proved node from scratch with two chains on the GPU (empty scratch directory) and
    checks preimages, bounds, tags and facts against genuine solutions. Outputs in `validation_pc/twochain/`; topgpu v2
    `5ac4892e…`, bottomgpu `6d277763…`, 8 workers, other jobs running:
    - 24192: E2E OK. The 6 committed preimages; 1.18·10⁹ units against 8.2·10⁹ in the committed records (7.0×); 323 s.
    - 1229312: E2E OK. No preimages, as committed; 5.8·10⁹ units against 6.7·10¹⁰ (11.6×); 1,859 s.
- **GPU coverage of the last rounds.** The two-chain plans end at G ≈ 37–38 (XL ≈ 10^12.6), where |θ| ≈ 2⁻⁵¹. topgpu v2 handed
  every job with |θ| ≥ 2⁻⁶⁰ to topc3, so these steps ran on the CPU. That is most of 1229312's wall time above, and 30% of
  221184's two-chain plan (2.4·10¹¹ of 8.0·10¹¹ units; the largest steps are 1–2.6·10¹⁰ cases each).
  - topgpu v3 accepts |θ|, |θ₂| < 2⁻⁴⁰ and widens every margin to 4(ε + η), η = 0.6·max(|θ|, |θ₂|). The margin is unchanged
    (64 units of 2⁻⁶⁴) whenever both are below 2⁻⁶⁰. The kernel is unchanged; only the host's coverage test and margins.
  - Build: `make -C work csrc/topgpu csrc/topgpu_emu` gives topgpu `f0afed4f…` and topgpu_emu `14904d28…`. The previous binary
    `5ac4892e…` is kept at /root/scratch/topgpu_v2_5ac4892e.
  - Validation, `validation_pc/even_gpu/val_topgpu_v3_twochain.out`: `PC_BOOTSTRAP=twochain python3 val_topgpu.py csrc/topgpu
    K C P F`, i.e. the jobs of the two-chain plans, many at G = 37–38. PASSED:
    - K: 45 jobs, worst ψ error 2 units of 2⁻⁶⁴;
    - C: 102 jobs byte-identical to topc3;
    - P: 150/150 planted genuine configurations kept, and the same with fp;
    - F1: 4000/4000 enumerations agree; F2: 96 jobs, 6,217 survivor runs, consistent.
  - `val_topgpu_v3_default_C.out` (today's plans, part C): 103 jobs byte-identical. Every small-G job of two plans now takes
    the GPU path, identical to topc3. A 2.8·10⁶-case job takes 0.23 s instead of 3.3 s, mostly launch overhead.
- **Meet-in-the-middle steps** (`docs/MITM_STEPS.md`, cost model `work/plan_mitm.py`; not implemented):
  - A top step paired with a table of bottom windows through a common lattice (hash join on a 3-torus). Each side needs
    only part of the precision.
  - A deep-digit split inside one step: the digits beyond the linear boundary are matched in two halves.
  - On the largest bag of each open node, either step chosen per round gains another 7.8–9.1× over two chains (347–394× over
    Matveev's chain alone). The middle rounds then dominate.
- **Tried and dropped** (cost model only):
  - a closing step with every digit in one of two windows: never competitive;
  - dropping the bottom threshold's tm + 8 floor: 1.7× on 2¹³·3³ (tmax = 55);
  - a smaller margin than +12: about 1.3–1.5× at tm = 55.
