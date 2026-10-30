# Node-2 process retrospective (fresh reviewer)

Evidence: RESEARCH_TREE.md, LESSONS.md, DEVIATIONS.md, NODE_MOVE.md, git log, agents/*.md headers.

## 1. What should we stop doing?
- **Letting exploratory work crowd the critical path.** At M5 the plan was "no new confirmatory endpoints", yet X4, X5, X6 and X7 (kinase, agents/kinase_report.md) ran while v3_tuned, the freeze, dry run and confirm were still "pending" (NODE_MOVE steps 6-9). Exploration only belongs alongside compute-bound waits.
- **Long jobs without internal checkpoints.** v3_tuned was killed for the node move and restarts from zero (1.5-2 h lost). Add per-step resume.
- **Late discovery of design flaws.** D22 (TM-align primary), D30 (CI-overlap collapse), D31/D32 (estimator bias) were all found after DESIGN v3 and three review rounds. Reviews did not find them; data did. Stop adding review rounds as a substitute for a pilot.
- **Unrecorded process rules colliding.** LESSONS had two rule 13s and a stale M3 count; tidy at every prune.
- **Exploratory analyses touching confirm HOGs** (C36, 118/197 pairs). Logged, but avoidable; use tune-only subsets by default.
- **Hand-copying agent output** when writes are blocked; specify the neutral filename up front.

## 2. What worked well?
- **Pilot on real data before freezing** (euk20 demo exposed TM-align failure, D22).
- **Synthetic-world verification before endpoints** (V2-V6) with honest reporting of both criteria (C34) and refusing to widen margins (E5).
- **DEVIATIONS.md discipline**: every change dated, justified and made pre-endpoint; immediate decisions on agent questions.
- **Parallel specialist agents with written reports** (18 reports in agents/), each with sentinels, seeds and "nothing committed".
- **Resumable, sentinel-based pipeline and NODE_MOVE.md**: node move needed no rebuilding; "things not to do" protected the split.
- **Early launch of Ancestra** while design was in review.
- **Licence/size hygiene** after two early incidents.
