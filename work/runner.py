import sys, json, os, subprocess
todo = json.load(open('todo_k9.json'))          # (kmax, nbags, node, target), sorted by kmax
part, nparts = int(sys.argv[1]), int(sys.argv[2])
mine = [t for i, t in enumerate(todo) if i % nparts == part]
for km, nb, s, d in mine:
    if os.path.exists(f'proofJ_{s}.json'): continue
    subprocess.run(['python3', 'prove_iter2.py', str(s)])
