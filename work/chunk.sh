#!/bin/sh
# One supervised work chunk: base 6 (node 32768) and the base-10 node list in parallel, both resumable.
cd "$(dirname "$0")"
timeout 560 python3 -u prove_cum_b.py 6 800 32768 >> runB6.log 2>&1 &
P1=$!
timeout 560 python3 -u runner4.py 8192 131712 122472 24192 >> runO_A.log 2>&1 &
P2=$!
wait $P1 $P2
grep -E "^base 6: node|THEOREM|^node .*: all|Traceback|survives" runB6.log runO_A.log | tail -3
ls -t ckptB | head -1; cat ckptB/$(ls -t ckptB | head -1) | python3 -c "import json,sys; d=json.load(sys.stdin); print('top', d['top'], 'bot', d['bot'], 'partial', (d.get('partial') or [None]*5)[:3], 'proofs', d['stats']['proofs'])"
