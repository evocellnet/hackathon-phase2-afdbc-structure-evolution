# Process retro M3 (t ≈ 4.5 h, 13:50)

Evidence: PROMPT/TREE/LESSONS/DEVIATIONS/NODE_MOVE, git log, agent report headers. No data read.

## 1. Stop doing

1. **Letting the record drift from reality.** RESEARCH_TREE's top summary is still M2 (t=1.8 h) and states "primary alignment is TM-align", but D22 (commit 459f41d) reversed that. The tree also holds [~] and [x] for the same pair pipeline (4d vs 4h). NODE_MOVE says vacate 20:30; START.txt says 20:52. A context reset would resume from wrong facts. Write an M3 summary now; one source of truth for deadlines.
2. **Over-engineering the statistics ahead of data.** About 20 pre-freeze deviations (D1–D18, C9–C21) were made on synthetic worlds. The result was that E5 acceptance is "practically unreachable" (β CI ±0.24). Stop adding engine features and sensitivity layers. Spend effort only on the endpoints that can decide something.
3. **Front-loading review rounds.** 9 reviewer passes consumed about 1.8 h of the first 4.5 h. The flaw that mattered (TM-align mis-registration, D22) came from running a 56k-pair demo, not from a reviewer. Future reviews should be short and data-triggered.
4. **Building non-critical side work while the critical path runs.** Annotations (4k), the cytoskeleton prep (4l), H1e and the P-clu baseline all ran on a node that Ancestra was using with 40 threads (reports note load 5–42). Side tasks are not free when timing estimates drive the 20:30 decision.
5. **Careless shell and git habits.** Self-matching `pkill`/`pgrep` appears twice in LESSONS, and a `git add` of 6.5 GB cost 10 min. Restricted TimeTree ages leaked into a committed table. Empty "tree" commits (7bb2a13, 01525f2) add noise. Stop using pattern kills and blind `git add`; add a pre-commit grep for restricted columns and large files.
6. **Letting deviation IDs skip and tables break.** D3, D7, D10–12, D15–17 are absent. D22 sits after the "pre-freeze" note, outside the table, and is post-freeze while being filed under pre-freeze. Fix the numbering and the headings so the audit trail is credible.
7. **Vague budgets.** The pair budget is 3.5M in C25 but ~2.5M in planning, and the HOG projection is 1.0–1.8M primary pairs. Settle one number before launching pair metrics.

## 2. What worked well

1. **Scouting before design**, run in parallel (M0). Data access and Ancestra feasibility shaped the design, including the "never use Ancestra duplication ages" decision (D2).
2. **Running slow work during review.** Panel, Ancestra launch and the distance library proceeded while the design was under review, and Ancestra started at 11:18.
3. **Verification before freezing.** 33 + 16 tests, V2–V9 and noise floors found real faults (HC1 under-coverage, node-age leakage, the V5 failure). Reviewers' rules for independent checking were followed.
4. **Pre-registered logic and honest logging.** Equivalence margins were not widened after E5 power was found to be poor. DEVIATIONS records reason and effect. C28 and C26 were fixed before endpoints were computed.
5. **Resumable pipelines.** `.done` sentinels, stage checkpoints and shards, plus NODE_MOVE with resume commands.
6. **Short agent reports with a "bottom line"** and PI decision requests (e.g. the E6 caps) made integration cheap.
7. **Measured scaling.** The 94 min families-stage time, with the finding that the stage is single-threaded, gave an ETA rather than a guess.

## 3. Pruned LESSONS.md (12 rules)

See the final list returned to the PI (identical).

1. Run the whole pipeline on a small real dataset before freezing: it caught D22, which no reviewer predicted.
2. Never `pkill -f`/`pgrep -f` a pattern in your own command; use `.done` sentinels or PID files (2 agents' waits hung; own shell killed).
3. Check `git status` and size before each commit; gitignore run dirs (a 6.5 GB `git add` cost 10 min).
4. Grep committed tables for licence-restricted columns before each commit (TimeTree age leak, M2).
5. Budget wall time, not CPU-h, for Ancestra: families stage is single-threaded, 94 min for 3.46M proteins.
6. Do not plan Ancestra on the GPU; it uses none of it (ancestra_report).
7. Launch slow jobs early and run design review and side work in parallel (Ancestra at 11:18 while review ran).
8. Keep review rounds short and shrinking (5 → 3 → 1); they sharpened design but cannot replace data-driven checks.
9. Every post-freeze change goes into DEVIATIONS with reason and effect before any endpoint is computed; never edit the plan silently.
10. Every long job needs a sentinel, a resume rule in NODE_MOVE and a stated ETA; update them when the plan changes.
11. Ask agents to return content inline if a write is blocked (round-3 reviewer).
12. AFDB downloads: use `_v6`, and `-A "Mozilla/5.0"` for EBI; aria2c -x8 makes bulk S3 cheap.

## 4. Top 5 risks, next 12 h

| # | Risk | Mitigation |
|---|---|---|
| 1 | Ancestra S2 or merge overruns past 20:30, and a mid-stage restart loses hours (S2 restarts from its start) | Check stage progress at 15:00, 17:00. Set a go/no-go at 18:30; if merge is not near done, pre-stage the rerun on node 2 and keep the stage checkpoints and `scores.tmp` copied. |
| 2 | CPU contention slows Ancestra | Freeze side jobs on node 1; run only GPU/light work there. |
| 3 | Pair-metrics throughput (4–5 h) leaves no tuning/confirmation margin | Shard by HOG; run tuning 20% first; set the pair budget now; make the endpoint code ready on euk20 data. |
| 4 | Silent state loss at node move (paths, env, sentinels, restricted data) | Dry-run NODE_MOVE steps on a test shard at 19:00. Back up `data/interim` and `results`; commit and verify. |
| 5 | Record drift and unreviewed late changes | Update the tree and deviation numbering now; require that tuning-only values are changed before the confirmation set is touched. |
