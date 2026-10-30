#!/bin/bash
# X6 first look: branch lengths on the fixed AA ML topology (binary; UFBoot<95 branches kept as nuisance
# parameters and not reported, because IQ-TREE -te resolves polytomies randomly).
# AA: LG+G4 ; 3Di: Q.3Di.AF+G4 (data/raw/models/Q.3Di.AF, md5 8f897ff...) on the same columns.
# Sensitivity: same topology, FoldMason structure-guided columns.
set -euo pipefail
R=/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution
CS=$R/data/interim/casestudy
Q3DI=$R/data/raw/models/Q.3Di.AF
cd $CS
for sf in ${SFS:-actin tubulin}; do
  ln -sf ${sf}_linsi_aa_ml.treefile iqtree/${sf}_foldmason_aa_ml.treefile
  for src in linsi foldmason; do
    iqtree3 -s aln/${sf}_${src}_aa.fasta  -te iqtree/${sf}_linsi_aa_ml.treefile -m LG+G4   -T ${NT:-2} --seed 20261030 --prefix iqtree/${sf}_${src}_aa_fixed_binary  -redo --quiet &
    iqtree3 -s aln/${sf}_${src}_3di.fasta -te iqtree/${sf}_linsi_aa_ml.treefile -m $Q3DI+G4 -T ${NT:-2} --seed 20261030 --prefix iqtree/${sf}_${src}_3di_fixed_binary -redo --quiet &
  done
  wait
done
grep -H -E "Gamma shape alpha|Log-likelihood of the tree|Total tree length" iqtree/*_fixed_binary.iqtree
