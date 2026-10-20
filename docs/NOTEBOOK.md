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
