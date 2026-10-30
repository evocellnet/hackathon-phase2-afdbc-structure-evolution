#!/bin/bash
# X7 step 7: (a) AA ML tree on the MAFFT L-INS-i columns (IQ-TREE3, LG+G4, UFBoot 1000, seed 20261030, 8 threads);
# (b) branch lengths on the FIXED binary ML topology for each alignment source:
#     AA LG+G4 on the AA columns, 3Di Q.3Di.AF+G4 (data/raw/models/Q.3Di.AF) on the 3Di letters of the SAME columns.
# UFBoot < 95 branches stay in the binary tree as nuisance parameters and are not interpreted (as X6).
# Usage: bash 07_iqtree.sh ml | fixed "linsi famsa foldmason"
set -euo pipefail
R=/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution
K=$R/data/interim/kinase
Q3DI=$R/data/raw/models/Q.3Di.AF
mkdir -p $K/iqtree; cd $K
if [ "$1" == "ml" ]; then
  iqtree3 -s aln/linsi_aa.fasta -m LG+G4 -B 1000 -T 8 --seed 20261030 --prefix iqtree/pkl_linsi_aa_ml -redo --quiet
  exit 0
fi
for src in $2; do
  iqtree3 -s aln/${src}_aa.fasta  -te iqtree/pkl_linsi_aa_ml.treefile -m LG+G4    -T 1 --seed 20261030 --prefix iqtree/pkl_${src}_aa_fixed  -redo --quiet &
  iqtree3 -s aln/${src}_3di.fasta -te iqtree/pkl_linsi_aa_ml.treefile -m $Q3DI+G4 -T 1 --seed 20261030 --prefix iqtree/pkl_${src}_3di_fixed -redo --quiet &
done
wait
grep -H -E "Gamma shape alpha|Log-likelihood of the tree|Total tree length" iqtree/*_fixed.iqtree
