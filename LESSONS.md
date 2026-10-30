# LESSONS — process rules for running this project efficiently
(One-line rules with evidence. Pruned at each milestone.)

1. Scout data access, tooling and the orthology software in parallel before designing — design depends on what is actually obtainable. (evidence: M0 plan)
2. AFDB files are now `_v6` (v4 returns 404) and EBI returns 403 to curl's default user agent — use `-A "Mozilla/5.0"`. (evidence: env agent, M0)
3. foldseek 10 GPU search (`makepaddedseqdb` + `--gpu 1`) works on the A100 in conda_envs/main. (evidence: env agent test)
