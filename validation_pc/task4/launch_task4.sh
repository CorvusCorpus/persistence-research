#!/bin/bash
# Task 4 launch (session 3): shells 40000 -> 100000, production on the GPU then the CPU second pass.
# Started detached from Windows:  Start-Process wsl.exe -ArgumentList -d,persistence,-u,root,--cd,<work/gsieve>,--exec,bash,-l,shells/launch_task4.sh
# Resume after an interruption: run it again (prod.py and val_proddone.py continue from their logs).
cd /mnt/c/Users/xlati/Desktop/persistence/persistence-research/work/gsieve || exit 1
exec >> shells/driver.log 2>&1
echo "== launch $(date -u +%FT%TZ) pid $$ host $(hostname) $(nvidia-smi --query-gpu=name,driver_version --format=csv,noheader)"
bash run_shells.sh shells "40000:50000 50000:60000 60000:70000 70000:80000 80000:90000 90000:100000" both
echo "== driver exit $? at $(date -u +%FT%TZ)"
