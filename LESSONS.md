# LESSONS — process rules for running this project efficiently
(One-line rules with evidence. Pruned at each milestone.)

1. Scout data access, tooling and the orthology software in parallel before designing — design depends on what is actually obtainable. (evidence: M0 plan)
2. AFDB files are now `_v6` (v4 returns 404) and EBI returns 403 to curl's default user agent — use `-A "Mozilla/5.0"`. (evidence: env agent, M0)
3. foldseek 10 GPU search (`makepaddedseqdb` + `--gpu 1`) works on the A100 in conda_envs/main. (evidence: env agent test)
4. Agents write multi-GB run dirs under results/; gitignore `results/**/runs/` and never `git add -A` without a `git status` check first. Never `pkill -f` a pattern that matches your own shell command. (evidence: a stuck 6.5 GB `git add` cost ~10 min, M0)
5. FoldComp index/lookup + afdb_rep_v4 (25 GB) downloaded from S3 in ~3 min with aria2c -x8 — bulk S3 downloads are cheap; don't over-deliberate them. (evidence: M0)
6. Three short review rounds (5 → 3 → 1 reviewers) turned a plausible design into a much sharper one: per-HOG branch lengths, pre-set equivalence margins, a positive control and power checks. Running data acquisition in parallel with review cost nothing. (evidence: M2)
7. Subagents may be blocked from writing some files; ask them to return content inline as a fallback. (evidence: the round-3 reviewer could not write ANALYSIS_PLAN.md)
8. Agents put licence-restricted values (TimeTree ages) into "small" committed tables. Grep committed tables for restricted columns before every commit. (evidence: `nearest_age_at_selection_mya` in species_panel.tsv, M2)
9. Run the pipeline end to end on a small real dataset before freezing tuning; it caught a primary-method flaw (TM-align mis-registration) that no reviewer predicted. (evidence: D22)
10. Ancestra's families stage is single-threaded Python over all DIAMOND hits: 94 min for 3.46M proteins. Budget wall time, not CPU-h, when extrapolating. (evidence: euk180 run)
11. `pgrep -f`/`pkill -f PATTERN` matches the caller's own command line. Use sentinel files (`.done`) or PID files for waits, never a self-matching pattern. (evidence: 2 agents' wait loops never exited; my own `pkill -f "git add -A"` killed my shell)
