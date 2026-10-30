# LESSONS — process rules for running this project efficiently
(One-line rules with evidence. Pruned at each milestone.)

1. Scout data access, tooling and the orthology software in parallel before designing — design depends on what is actually obtainable. (evidence: M0 plan)
2. AFDB files are now `_v6` (v4 returns 404) and EBI returns 403 to curl's default user agent — use `-A "Mozilla/5.0"`. (evidence: env agent, M0)
3. foldseek 10 GPU search (`makepaddedseqdb` + `--gpu 1`) works on the A100 in conda_envs/main. (evidence: env agent test)
4. Agents write multi-GB run dirs under results/; gitignore `results/**/runs/` and never `git add -A` without a `git status` check first. Never `pkill -f` a pattern that matches your own shell command. (evidence: a stuck 6.5 GB `git add` cost ~10 min, M0)
5. FoldComp index/lookup + afdb_rep_v4 (25 GB) downloaded from S3 in ~3 min with aria2c -x8 — bulk S3 downloads are cheap; don't over-deliberate them. (evidence: M0)
