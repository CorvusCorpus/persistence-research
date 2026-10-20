# Validation logs from the desktop PC

Logs of the validation runs on the desktop (Windows 11, WSL2 distribution `persistence`, i7-12700F, RTX 3060 Ti).
`docs/NOTEBOOK.md` (session 3) describes each run: command, binaries, commit, timing and outcome.

Step 1 (20 Oct 2026, commit `c471a35`, binaries built from `gsieve.c` `3f4827dc…`):

| file | run |
|---|---|
| `make_test.log`, `make_test_prod.log` | `make test`, `make test-prod` from the repository root |
| `val_all.out` and `A_*`, `B_*`, `C_*`, `D_*`, `build.log`, `binaries.sha256` | `cd work/gsieve && ./val_all.sh ../../validation_pc` |
| `bench.txt` | `python3 work/gsieve/bench.py` (two workloads) |
| `audit.txt` | `audit/run_audit.sh` |

Later steps keep their logs in subdirectories.

`task3_cuda/` (task 3, CUDA variant; `docs/NOTEBOOK.md` 3.1–3.4):

| file | run |
|---|---|
| `kerneltest.txt`, `suite_cuda_T{1,8}.log`, `prodrow_*.log`, `xvariant.txt` | first checks of the CUDA build (3.2) |
| `rerun_cuda.txt`, `rerun_cuda.tar.xz` | full re-run of the logged production with the first CUDA build (3.2) |
| `bench_cuda.txt` | benchmarks (3.3) |
| `final_validation.out`, `val_cuda/`, `val_all/` | `val_cuda.sh` and `val_all.sh` on the final binaries (3.4) |
| `rerun_cuda_final.txt`, `rerun_cuda_final.tar.xz` | full re-run of the logged production with the final CUDA binary (3.4) |

The `.tar.xz` archives hold the re-run logs in `aggregate.py` format. Unpack one with `tar -xJf`, then run, for
example, `python3 work/gsieve/aggregate.py 36000 -1 validation_pc/task3_cuda/rerun_cuda_final`.
