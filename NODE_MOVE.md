# NODE_MOVE.md: checkpoint and resume protocol

*Living document; last updated 2026-10-30 ~17:45 on node 1. Final update before 20:30.*

## On the new node: first steps
1. Go to the project directory:
   ```
   cd /scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution
   ```
2. Start the node clock:
   ```
   printf "NODE 2 START: %s\nSTART_EPOCH: %s\n" "$(date -Is)" "$(date +%s)" >> START.txt
   python3 scripts/elapsed.py
   ```
3. Activate the environment and check foldseek and the GPU:
   ```
   source /opt/conda/etc/profile.d/conda.sh && conda activate $PWD/conda_envs/main
   foldseek version && python -c "import torch; print(torch.cuda.is_available())"
   ```
4. Read `RESEARCH_TREE.md` (status at the top), then `DEVIATIONS.md` (all pre-freeze changes), then the job table below.
5. Check that the node-1 jobs really stopped: `ps aux | grep -E "run_pairs|run_tuning|v3_tuned" | grep -v grep` should print nothing.

## Pipeline state and resume rules
Every step is idempotent. A step is complete when its sentinel exists; rerunning skips it.

| # | Step | Sentinel / output | Command | Status at last update |
|---|---|---|---|---|
| 1 | Ancestra euk180 (accurate) | `data/interim/ancestra_euk180/ancestra/hogs.orthoxml` | (done; do not rerun) | DONE 16:20 (5.0 h, 161 GB) |
| 2 | Panel tree (D30, no CI-overlap collapse) | `data/interim/timetree/panel_tree.nwk`, `node_ages.tsv` | `python scripts/hogs/panel_tree.py` | DONE |
| 3 | HOG build (7000 HOGs; split.tsv md5 is in RESULTS_LEDGER) | `data/interim/hogs.parquet`, `split.tsv`, `pairlists/` | **never rebuild** (split must not change) | DONE 16:43 |
| 4 | Pair metrics a–c4 and H1e step d | `logs/pipeline/<step>.done` plus shard `.done` files | `ALLOW_BUILD=0 THREADS=48 scripts/run_after_ancestra.sh` (resumes at the first unfinished shard); progress: `scripts/run_after_ancestra.sh status` | a–c4 DONE 18:11 (1.40M pairs × 3 alignments); d (H1e) running since 18:11 |
| 5 | Tuning part 1 (α from seqonly (C29), floors, boot-sub column dumps, SE calibration) | `logs/tuning/<step>.done` | `PAIRS_THREADS=48 bash scripts/tuning/run_tuning.sh --only "validate_profile_alpha tune_alpha yaml_alpha floors_tuned select_bootsub bootsub_pairs se_calibration yaml_se"` | DONE 17:48 (α AA 1.64, 3Di 1.60) |
| 6 | Tuning part 2 (only after step 4/c4, because apply_tuning rewrites the pair tables in place) | `logs/tuning/*.done` | `PAIRS_THREADS=48 V3_ARGS="--B 100 --pb-K 10" bash scripts/tuning/run_tuning.sh --from apply_tuning` (resumable per step) | launched 18:12 |
| 7 | Fill the §Tuning table, then freeze | `config/tuning.yaml` has no nulls | `python scripts/tuning/fill_tuning_table.py --plan ANALYSIS_PLAN.md --tuning config/tuning.yaml`, then commit, then add a ledger line with the md5 of ANALYSIS_PLAN.md and tuning.yaml | pending |
| 8 | Tune dry run of endpoints (exploratory, `results/tuning/`) | `results/tuning/dryrun/*.json` | `python scripts/endpoints/run_dryrun.py` | pending |
| 9 | Confirmation run (**once**) | `results/endpoints/E*.json` | `CONFIRM_RUN=1 python scripts/endpoints/run_confirm.py && python scripts/endpoints/interpret.py --dir results/endpoints` | pending; needs a clean git tree and the ledger md5 |
| 10 | Exploratory X1/X2/X7 on tune, later on all data | `results/exploratory/` | `python scripts/exploratory/x1_curves.py --aln tdiaa` (also x2_age, x7_pclu) | pending |

## Things not to do
- Do not rerun `build_hogs.py` and do not regenerate `split.tsv`.
- Do not load `split=='confirm'` outside `run_confirm.py`.
- Do not commit anything under `data/interim/timetree/`, or any file that contains ages.
