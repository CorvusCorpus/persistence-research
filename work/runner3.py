"""Run prove_iter2 node by node from an explicit list; resumable (bag results and proof checkpoints)."""
import sys, os, subprocess
for s in sys.argv[1:]:
    if os.path.exists(f'proofJ_{s}.json'): continue
    subprocess.run(['python3', '-u', 'prove_iter3.py', s])
