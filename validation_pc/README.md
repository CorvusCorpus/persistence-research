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
