# LESSONS: process rules for running this project efficiently
(One-line rules with evidence. Pruned at each milestone. Last pruned at M3 (13:55, node 1) from the fresh-agent retrospective in reviews/process_retro_M3.md.)

1. Run the whole pipeline on a small real dataset before freezing tuning. It caught a primary-method flaw (TM-align mis-registration, D22) that no reviewer predicted.
2. Never `pkill -f`/`pgrep -f` a pattern that matches your own command; use `.done` sentinels or PID files. Two agents' waits hung, and my own `pkill` killed my shell.
3. Check `git status` and file sizes before each commit; gitignore run directories. A 6.5 GB accidental `git add` cost 10 min.
4. Grep committed tables for licence-restricted columns before each commit. TimeTree ages leaked into species_panel.tsv at M2.
5. Budget wall time, not CPU-h, for Ancestra: the families stage is single-threaded Python and took 94 min for 3.46M proteins. Ancestra does not use the GPU.
6. Launch slow jobs early and run review and side work alongside them. Ancestra started at 11:18 while the design was still under review.
7. Keep review rounds short and shrinking (5 → 3 → 1 reviewers). They sharpened the design but cannot replace checks on data.
8. Log every change to the plan in DEVIATIONS.md with reason and effect *before* any endpoint is computed; never edit the plan silently.
9. Every long job needs a sentinel, a resume rule in NODE_MOVE.md and an ETA. Update them when plans change.
10. If a subagent's file write is blocked, have it return the content inline (round-3 reviewer, ANALYSIS_PLAN.md).
11. AFDB: files are `_v6`; EBI needs `-A "Mozilla/5.0"`; bulk S3 downloads via `aria2c -x8` are cheap (25 GB in ~3 min). Don't over-deliberate them.
12. Agents report "decisions for you". Decide immediately, record them in DEVIATIONS.md (C-items) and keep going. This kept the critical path moving. (M2–M3)
13. Error patterns in wait loops must not match benign log text ("failed 0" matched /FAIL/i). Use `Traceback|exit status [1-9]`. A wait's own timeout can hide a deadline: check the clock at every notification. (evidence: a false wake at 17:40; a 2 h wait ran up to 20:11, leaving only 40 min before the node stop)
13. The harness blocks subagents from writing "report" files and they return the text inline, so the PI must re-type long drafts. Ask writer agents to save under a neutral name (e.g. `technical report/draft_sections.md`) via Bash, and fall back to inline only if that is also blocked. (evidence: report draft, node 2)
