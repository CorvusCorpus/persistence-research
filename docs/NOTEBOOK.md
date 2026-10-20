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
