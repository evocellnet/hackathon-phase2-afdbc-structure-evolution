# Round 3 review: consistency and executability (DESIGN v3)

2026-10-30. **Verdict: freeze after the MUST items** (all adopted in `ANALYSIS_PLAN.md`).

## MUST

| # | Problem | Fix |
|---|---|---|
| M1 | E4 sign flips. §11.2 defines y = log(3Di/AA), but the §11.3 table is written for k = AA/3Di. | Fit y; read the table on β_k = −β. |
| M2 | Confidence levels clash: E5 uses 90%/95%, the multiplicity rule uses Bonferroni. α and m are undefined, and E1 and E6 each hold two tests. | α = 0.05, m = 5, so 98% CIs and 99% one-sided bounds. E1: "linear" beats "concave". E6: p = 2·min(p_a, p_b). Holm over the five endpoints. |
| M3 | E5 noise floor defined twice (bootstrap variance in §11.2, parametric bootstrap in §11.5). | Floor = parametric bootstrap under N0. CI = pigeonhole bootstrap. |
| M4 | Pruned HOG edges span several panel edges, so E4 ages and E5 lineages are undefined. | Map each HOG edge by the panel MRCAs of its child and parent clades. E5 uses edges that map one-to-one; E3 and E4 use all edges. |
| M5 | The geometric midpoint is 0 for terminal branches, and the bin edges are missing. | Arithmetic midpoint; bin edges 0/30/100/300/1000 Mya. |
| M6 | SIMEX needs per-pair bootstrap SEs, which exist only for a subsample. The engine gives SE = NaN at d = 0. | Fisher SEs, calibrated on the subsample (inflation factor, AA–3Di correlation ρ). SE floor = 1/n_core. |
| M7 | The §11.4 paralog swap cannot be computed: the 1:1 sets drop multi-copy species, the new pairs are not computed, and the effect size is unknown. | Power = synthetic ×1.5 3Di inflation on 10% of lineages, 50 replicates. The swap becomes an exploratory check. |
| M8 | Selection "on AA only" uses TM-align D_AA, which is structural. "Clade-closed" is undefined. | Use sequence-only D_AA (Biopython; parasail is absent). The algorithm is given in the plan. |
| M9 | Resampling species duplicates tree tips, and SIMEX × 1,000 bootstrap replicates is unaffordable. | Pigeonhole weights with no refit. SIMEX inside 200 replicates (B = 10). |
| M10 | Foldseek's `lddt` averages aligned residues only, and LDDT is asymmetric. | Query = shorter protein; LDDT = Σ`lddtfull` / L_query. The in-house LDDT is validated to ≤ 0.01. |
| M11 | Columns used for the ML distances are unspecified, which moves k. | Core-masked columns as primary; all aligned columns as sensitivity. |
| M12 | Unconstrained WLS can give Σb_3Di ≤ 0, so the log is undefined. | Merge cells on AA only. Censor at 0.005 and count. |
| M13 | statsmodels MixedLM takes no weights; R and FastME are not installed. | Unweighted mixed model as primary. NJ (dendropy) replaces FastME. |

## SHOULD

1. Use one ratio direction internally (3Di/AA). Report k = 1/r.
2. Saturated pairs count as failures (weight 0) only inside the ≤10% HOG gate. Sensitivity: impute at the cap.
3. Restore the stats-M3 E1 rules:
   - D_AA,seq ∈ [0.05, 2];
   - ≥10 pairs per HOG, max/min ≥ 4;
   - η fixed (free as sensitivity);
   - a_h profiled.
4. §3 additivity: report it, do not gate on it.
5. The §11.5 misfit exclusion applies to E3–E5 only.
6. N0 has three states (accept, reject, inconclusive). §11.3 omits β > 1.25. Power gates acceptance only.
7. E3 CATH class = residue-weighted TED majority (≥ 60%). Unadjusted secondary test.
8. Split halves: alternating 20-residue blocks, because 3Di depends on neighbouring residues.
9. E6 residuals: use the E1 HOG curve. O1 = sister clade of the duplication level; O2 = next sister.
10. E3 counts as alignment-robust if k agrees within 25%.
11. Age draws: one lognormal draw per bootstrap replicate, with a monotone repair.
12. Seeds: `SeedSequence(20261030)` spawned per stage.
13. Store 20×20 count matrices. The confirmation script checks that tuning values are filled and the plan hash is recorded.
14. Engine settings:
    - `d_max = 20` (the default is 10);
    - Q.3Di.AF read in PAML ARND order;
    - validate against IQ-TREE.
15. Add an end-to-end synthetic test: an N0 world and a spike-in world.

No other contradictions found.
