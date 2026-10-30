"""strucevo.stats -- statistics layer of ANALYSIS_PLAN.md (branch engine, pairwise curvature engine,
pigeonhole bootstrap, endpoints E1-E6, Holm, interpretation, synthetic worlds).

Modules
-------
common            seeds (config/seeds.yaml), confirm-split guard, provenance, CIs, Holm, JSON writer
trees             panel tree, pruning to S_h, design matrices, NJ, age draws
branch            branch engine (two-pass WLS, sandwich SEs, cov_aa, tip correction, misfit, hidden paralogy)
pairwise          pairwise curvature engine (profiled a_h, global c, SIMEX) on torch/GPU
bootstrap         pigeonhole weights
reml              fast REML for the E4 mixed model (bootstrap replicates)
endpoints_curve   E1, E2 (+ d*)
endpoints_branch  E3, E4 (cells, merge rule, MixedLM, age redraws)
endpoints_lineage E5 (lineages, floor F, sigma_L, sister contrasts, spike-in power, N0 decision)
endpoints_misc    E6, split halves, Holm, interpretation (§11.3 row ids)
synth             synthetic worlds for V2-V4, V6
"""
from .common import load_pairs, load_branches, stage_rng, ConfirmSplitGuard, holm  # noqa: F401
