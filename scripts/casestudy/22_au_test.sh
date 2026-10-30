#!/bin/bash
# X6 node 2: AU test (RELL, 10,000 replicates) of the constrained topology vs the unconstrained ML tree, LG+G4, MAFFT columns.
R=/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution
CS=$R/data/interim/casestudy
cd $CS
for sf in ${SFS:-actin tubulin}; do
  cat iqtree/${sf}_linsi_aa_ml.treefile constraint/${sf}_cons.treefile > constraint/${sf}_au_input.nwk
  iqtree3 -s aln/${sf}_linsi_aa.fasta -m LG+G4 -z constraint/${sf}_au_input.nwk -n 0 -zb 10000 -au -T ${NT:-2} --seed 20261030 \
    --prefix constraint/${sf}_au -redo --quiet
  grep -A6 "Tree      logL" constraint/${sf}_au.iqtree
done
