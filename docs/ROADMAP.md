# Roadmap

Ordered by value per effort. Each item lists the acceptance test that must pass before any claim is made.

## P0 — make the repo solid (first session)
1. **Initialise the repo.**
   - git init with a `.gitignore` covering binaries, `__pycache__` and large logs; Git LFS or a release
     tarball for logs.
   - Add `make all`, `make test` and a quick smoke target.
2. **Check the environment.**
   - Toolchain: gcc ≥ 11, GMP, MPFR, Python ≥ 3.10, `python-flint`, `mpmath`.
   - CPU flags: the chain sieve needs AVX-512 F/DQ/BW/VL/IFMA/VBMI2 plus BMI2.
   - GPU: `nvidia-smi`, `nvcc --version`.
3. **Re-run the validation suites.**
   - `make -C work/gsieve test`.
   - Base-10 smoke: in a scratch copy with empty `bagres*`/`ckpt*`, `python3 prove_cum.py 98 147` must
     reproduce `docs/results_base10.md`.
   - A heavier spot check, e.g. `val_cum_e2e.py`.
   - Acceptance: all scripts print OK / exact match.

## P1 — portability and speed of the chain sieve (CPU)
1. **AVX2 fallback** for CPUs without AVX-512 (most consumer Intel). DONE in session 2: `avx512` (no IFMA)
   and `avx2` variants with exact integer division (see `docs/NOTEBOOK.md`, task 2).
   - Acceptance: identical stage counts and S2/R3/R4 sets on `val_dump.py` and `val_prodrow.py`, against
     the AVX-512 build or brute force.
2. **Multithreading** (pthreads) sharing one copy of the tables, with per-thread c-blocks and checkpoint
   lines. DONE in session 2: dynamic per-c scheduling, per-c output written in c order.
3. **Optional K = 14 tables** (≈4.7 GB shared) for about 20% less work. Re-validate.

## P2 — CUDA port (RTX 3060 Ti) and the 10^5-digit record
1. **Port stage 1 (bitmap scan) and stage 2 (affine lift) to CUDA.**
   - One thread block per row segment; the bitmap and AB/H0/POS tables live in device memory (≈1 GB).
   - Emit stage-2 survivors (positions) to the host, where stages 3–5 run unchanged.
   - DONE in session 3: `gsieve_cuda.cu`, variant `cuda` (`make cuda`), with fail-stop table checks for non-ECC
     memory (`docs/CHAIN_SIEVE.md`, "CUDA variant").
2. **Validation:** the GPU stage-2 survivor sets must equal the CPU sets exactly on test regions (several
   c-slices at N = 36,000, and E-mode at N = 160). Measure throughput.
   - DONE in session 3 (`docs/NOTEBOOK.md`, task 3), all by `val_cuda.sh`:
     - kernel tests on the whole table;
     - the suite (E mode 314/314);
     - eight production slices against exact arithmetic;
     - cross-variant dumps;
     - a full re-run of the logged 10^40000 production, with every counter and per-c count equal.
   - Throughput: 1.45–2.1·10¹¹ points/s, 4.7–6× the CPU. Stages 3–5 on the CPU limit it; stage 3 on the GPU would
     gain at most about 1.3×.
3. **Extend in shells** with `n0`, archiving logs at each step: 40,000 → 50,000 → 60,000 → … → 100,000.
   - Acceptance: `aggregate.py N n0` shows full coverage of every path, no members, and survival rates
     that match theory.
   - DONE in session 3 (20 Oct 2026): six shells, 9.36·10^14 candidates in 86 min, all acceptance checks pass
     (`docs/CHAIN_SIEVE.md`, "Shells"; `docs/NOTEBOOK.md`, task 4).
4. **New claim if successful:** no persistence-12 number below 10^100000.
   - DONE: recorded in `docs/RESULTS.md` (1.1, 1.2) as COMPUTATIONAL / THEOREM given 1.1.
5. **TODO: independent second pass of the six new shells** with the CPU binary (about 8.6 h on the desktop):
   `cd work/gsieve && ./run_shells.sh shells "40000:50000 50000:60000 60000:70000 70000:80000 80000:90000 90000:100000" cpu`.
   It re-runs every logged invocation (all eight DONE counters must match) and compares per-c counts. The desktop has
   no ECC memory; recommended before submitting the OEIS comment.
   - Session 4 (cloud, ECC hardware likely): 40,000 → 60,000 done, every counter and per-c count equal
     (`validation_pc/task4_check/`); stopped by decision of the user. Remaining: `"60000:70000 … 90000:100000" cpu`.
   - Session 5: the main path of all of m ≤ 100,000 re-executed by the GPU tree sieve, every counter equal (NOTEBOOK T.7).
     Only the direct and odd paths of 60,000 → 100,000 (cheap: minutes each on the CPU) still lack a second execution.
6. **Tree sieve (session 4, `work/tsieve/`, NOTEBOOK T).** The points that pass D digits form 4^D cosets of a lattice of
   index 4·5^(D−1). Walking the coset tree costs about N^(3 log 4 / log 5) = N^2.58, against N³ for the bitmap scan.
   - `tsieve6` (CPU): identical counters and dumps to gsieve (G and E); 314/314 positive control. About 1.2× slower
     than gsieve at N = 10⁵ on the cloud CPU; crossover estimated near N ≈ 3·10⁵.
   - `tsieve_gpu` (CUDA): accepted on the desktop and made 3× faster (kernel v3, NOTEBOOK T.7). The main path for
     m ≤ 10⁵ takes 45 min on the RTX 3060 Ti against about 79 min for gsieve's CUDA variant, and its cost grows as
     N^2.58. It is the faster main path for any next shell (estimate: 10⁵ < m ≤ 2·10⁵ in about 3.7 h; direct and odd paths
     by gsieve). A new shell needs the user's go-ahead and a second execution before any claim.

## P3 — write-up and community verification
- **Paper draft** (LaTeX):
  - chain lemma, suffix collapse, exact 4^k branching, affine lift, odd thread, bridge lemma;
  - computation and validation;
  - heuristic and calibration;
  - small-base theorems and base-10 trees as a second part (or a separate note).
- **OEIS A003001 comment draft** with the new bound and a code link.
- **Independent verification kit:** a minimal re-implementation (e.g. a slow, simple Python/numpy
  version of stages 1–2 for spot regions) so that others can audit without trusting the fast code.

## P4 — base-10 even-target trees (CPU-heavy)
- Resume `runner4.py 122472` and `runner4.py 131712 24192`. Checkpoints are in `work/ckptC`, `work/bagresC`.
- Then 24192, 41472, 221184, 93312, 314928, 1229312, 11239424 by cost.
- 1161216 needs about 40,000 core-hours with the current method. It would need a new idea, e.g. stronger
  bottom tests or reuse of facts across bags.

## P5 — research directions (speculative)
- **Zero-free-mode sieve** to extend Chaffin's "no zero-free 7-smooth in (10^140, 10^20000)" (OEIS A238985 has no
  stated search bound; it "seems finite, 12615 terms"). NEXT (chosen in session 3). The form 2^a3^b7^c runs on the
  GPU main path in E mode (region m ≤ 1.183·D covers x < 10^D). The form 3^a5^b7^c needs a 2-adic suffix collapse:
  for 5^K | x the last K digits depend only on (−1)^(a+c)·5^(b + a·f3 + c·f7) mod 2^K, so a 2^(K−1)-bit bitmap
  (4 MB for K = 26) and an affine lift mod 2^K replace point-by-point stepping. About 7% of candidates survive
  stage 2 in this mode, so the CPU later stages set the pace.
- **A007377** (2^n with no digit 0, "checked up to k = 10^10", Radcliffe 2022). IN PROGRESS (session 3): `work/pow2/z2sieve`
  (CPU) and `z2sieve_cuda` (GPU, 6.8·10^10 exponents/s: 10^15 in about 4 h) are built and validated
  (`docs/ZEROLESS_POW2.md`, NOTEBOOK P). NEXT: production on the GPU once the zero-free run is done (runs over 12 h
  need the user's approval), then a second pass and the OEIS comment.
  Original note: the same machinery as a single row
  could reach about 10^14–10^15. The ternary analogues (Erdős's digit-2 problem, zero-free 2^n in base 3) are
  already verified to 2·3^45 ≈ 5.9·10^21 (arXiv 2202.13256); our sieve cannot compete there.
- **Rigorous density theorem:** #chain numbers with m ≤ N is ≪ N^θ with θ < 3, via transfer operators
  plus effective equidistribution.
- **Bridges:** special structure (≥ 99.4% ones), possibly with mixed archimedean/p-adic reductions for
  the D = 9 targets (15–29 non-1 digits).
- **Base 3 final digit 0:** the zero-free ternary powers of 2 question (2^15 = 1122221122₃). Same sieve
  ideas with 3-adic structure.
