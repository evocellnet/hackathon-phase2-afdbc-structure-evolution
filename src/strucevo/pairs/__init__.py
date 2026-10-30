"""strucevo.pairs: per-pair alignment + divergence table (ANALYSIS_PLAN §1.2).

store     FoldComp subset -> Foldseek DB + ragged CA / pLDDT / AA / 3Di arrays (+ sequence md5)
align     Foldseek tmalign / structurealign (3Di+AA) on exhaustive given-pair prefilters; Biopython seqonly
posthoc   core mask, counts, lDDT (lddtfull or strucevo.lddt), core RMSD
ml        ML distances from counts (LG / Q.3Di.AF + Gamma4, d_max 20, SE floor), floor-corrected lddt_D
pipeline  sharded resumable runner with .done sentinels
"""
