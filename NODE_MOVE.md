# NODE_MOVE.md — checkpoint / resume protocol (living document; finalised before 20:30 on node 1)

## On the new node: first steps
1. `cd /scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution`
2. Append a new start line: `echo "NODE START: $(date -Is)\nSTART_EPOCH: $(date +%s)" >> START.txt`, then run `python3 scripts/elapsed.py`.
3. Check the envs: `source /opt/conda/etc/profile.d/conda.sh && conda activate $PWD/conda_envs/main && foldseek version && python -c "import torch; print(torch.cuda.is_available())"`.
4. Read `RESEARCH_TREE.md` (status at top), then the job table below.

## Long jobs (state at last update)
| Job | Command | Output / sentinel | Resume rule |
|---|---|---|---|
| Ancestra euk180 (accurate, 40 thr) — launched 2026-10-30 11:18; families done 12:52 (94 min); scores done 15:20 (148 min, 158 GB); merge running since 15:20 | `conda_envs/ancestra/bin/python scripts/ancestra/run_ancestra.py --proteomes data/interim/ancestra_input --taxdump data/interim/ancestra/taxdump --out data/interim/ancestra_euk180 --threads 40 --preset accurate --no-labels` | `data/interim/ancestra_euk180/ancestra/hogs.orthoxml`; stage checkpoints in `.../ancestra/stages/` | Stages are skipped on rerun if complete. A stage interrupted mid-way (esp. S2) restarts from that stage's beginning: rerun the same command first thing on node 2. |
