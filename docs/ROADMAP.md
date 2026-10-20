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
1. **AVX2 fallback** for CPUs without AVX-512 (most consumer Intel).
   - Acceptance: identical stage counts and S2/R3/R4 sets on `val_dump.py` and `val_prodrow.py`, against
     the AVX-512 build or brute force.
2. **Multithreading** (pthreads) sharing one copy of the tables, with per-thread c-blocks and checkpoint
   lines.
3. **Optional K = 14 tables** (≈4.7 GB shared) for about 20% less work. Re-validate.

## P2 — CUDA port (RTX 3060 Ti) and the 10^5-digit record
1. **Port stage 1 (bitmap scan) and stage 2 (affine lift) to CUDA.**
   - One thread block per row segment; the bitmap and AB/H0/POS tables live in device memory (≈1 GB).
   - Emit stage-2 survivors (positions) to the host, where stages 3–5 run unchanged.
2. **Validation:** the GPU stage-2 survivor sets must equal the CPU sets exactly on test regions (several
   c-slices at N = 36,000, and E-mode at N = 160). Measure throughput.
3. **Extend in shells** with `n0`, archiving logs at each step: 40,000 → 50,000 → 60,000 → … → 100,000.
   - Acceptance: `aggregate.py N n0` shows full coverage of every path, no members, and survival rates
     that match theory.
4. **New claim if successful:** no persistence-12 number below 10^100000.

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
- **Zero-free-mode sieve** to extend Chaffin's "no zero-free 7-smooth in (10^140, 10^20000)". This needs
  vectorised later stages, since about 7% of candidates survive stage 2 in that mode.
- **Rigorous density theorem:** #chain numbers with m ≤ N is ≪ N^θ with θ < 3, via transfer operators
  plus effective equidistribution.
- **Bridges:** special structure (≥ 99.4% ones), possibly with mixed archimedean/p-adic reductions for
  the D = 9 targets (15–29 non-1 digits).
- **Base 3 final digit 0:** the zero-free ternary powers of 2 question (2^15 = 1122221122₃). Same sieve
  ideas with 3-adic structure.
