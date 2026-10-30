#!/bin/bash
# X6 node 2: after the constrained searches finish, run (i) tubulin column bootstrap on the ML tree, (ii) fixed-topology
# lengths + bootstrap (B=50) on the constrained trees, (iii) LG+C20+G4 AA lengths on the ML tree (saturation/LBA check).
R=/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution
cd $R
until [ -f data/interim/casestudy/constraint/actin_cons.treefile ] && [ -f data/interim/casestudy/constraint/tubulin_cons.treefile ]; do sleep 20; done
python scripts/casestudy/16_iqtree_fixed.py --sf tubulin --src linsi famsa foldmason --tree ml --B 100 --workers 4 > logs/casestudy/iqfix_ml_tub.log 2>&1 &
python scripts/casestudy/16_iqtree_fixed.py --sf actin tubulin --src linsi famsa --tree cons --B 50 --workers 2 --full_threads 2 > logs/casestudy/iqfix_cons.log 2>&1 &
python scripts/casestudy/16_iqtree_fixed.py --sf actin --src linsi --tree ml --B 0 --model_aa LG+C20+G4 --tag _C20 --full_threads 2 > logs/casestudy/iqfix_c20_actin.log 2>&1 &
wait
python scripts/casestudy/16_iqtree_fixed.py --sf tubulin --src linsi --tree ml --B 0 --model_aa LG+C20+G4 --tag _C20 --full_threads 4 > logs/casestudy/iqfix_c20_tubulin.log 2>&1
