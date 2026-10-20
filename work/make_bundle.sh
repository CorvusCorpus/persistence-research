#!/bin/sh
# Package code, certificates and summaries.
cd "$(dirname "$0")"
rm -f persistence_proofs.zip
zip -q -r persistence_proofs.zip README_bundle.md METHOD.md results.md smallbases.md certificate.json \
  base3_cum.json base4_cum.json base6_cum.json trees_b6_D800.json \
  bags.py twoadic.py solver.py fastbottom.py prove_iter.py prove_iter2.py prove_iter3.py prove_cum.py prove_cum_b.py \
  basesolver.py basetrees.py ctop.py ctop2.py cbottom.py cbottom2.py cbottomb.py runner4.py \
  audit_M.py val_ctop.py val_ctop2.py val_ctop_cum.py val_cbottom.py val_cbottom2.py val_cbottom_cum.py val_iter3.py \
  val_cum_e2e.py val_topb.py val_topb_fast.py val_bottomb.py planted_cum.py planted_top_cum.py planted_bottomb.py planted_topb.py \
  results_table.py smallbases_table.py certificate.py cache_trees_fast.py \
  csrc/top.c csrc/top_v1.c csrc/bottom.c csrc/bottom_v1.c csrc/bottom2.c csrc/topb.c csrc/bottomb.c cbottom3.py val_cbottom3.py planted_cum3.py planted_bottomb2.py \
  proofJ_*.json proofC_*.json 2>/dev/null
ls -la persistence_proofs.zip
