# Post-hoc extra robustness checks: fold class, E6b and E1 under pLDDT control (2026-10-31)

**Everything in this note is POST HOC and exploratory.** None of it changes a pre-registered decision.

- Script: `scripts/posthoc/extra_checks.py`. Run with `CONFIRM_RUN=1`, confirm split, read-only. Production code in `src/strucevo/stats` was used: `prepare_curve_data`, `CurveEngine`/`two_pass`/`simex`/`fit_batch` and `run_E6`. The production branch tables `data/interim/branches/tdiaa*.parquet` were also used.
- Output: `results/posthoc/robustness_extra.json`. Log: `logs/posthoc/extra_checks.log`. Runtime is about 6 minutes.
- **Reproduction checks pass exactly:**
  - E3 k_all = 3.3218, k_internal = 4.6280, k_α = 4.677, k_β = 2.552;
  - E6b = 0.017024 (n = 6,027 duplication events);
  - E1 naive ĉ = 0.6695;
  - plan SIMEX (no D29 correction) = 0.7165, against 0.7166 in E1.json.

## 1. Fold-class contrast on internal edges only

- Method: k is a ratio of sums over the HOGs E3 keeps (4,891). Classes come from the E3 CATH rule.
- Uncertainty: HOG-only multinomial bootstrap, B = 1000, seed 20261030, 95% percentile CIs. This bootstrap has no panel-edge weights, so its CIs are narrower than those in E3.json.

| class | HOGs | mean HOG pLDDT | k all edges | k internal | k terminal | terminal share of 3Di length |
|---|---|---|---|---|---|---|
| α | 1,336 | 80.5 | 4.68 [4.52, 4.86] | **6.11 [5.94, 6.31]** | 4.18 | 0.74 |
| α/β | 2,261 | 84.3 | 3.16 | 4.31 | 2.84 | 0.78 |
| mixed | 375 | 80.3 | 3.18 | 4.60 | 2.79 | 0.79 |
| β | 696 | 81.8 | 2.55 [2.46, 2.65] | **3.85 [3.74, 3.97]** | 2.20 | 0.79 |
| **α/β ratio** | | Δ = −1.3 | 1.83 [1.74, 1.93] | **1.59 [1.52, 1.65]** | 1.90 | |

**α/β ratio within tertiles of HOG pLDDT.** Mean pLDDT is matched between classes within each tertile, to within 0.2.

| tertile | ratio, all edges | ratio, internal edges |
|---|---|---|
| T1 (< 78.7) | 2.11 | 1.72 |
| T2 | 1.72 | 1.52 |
| T3 (> 85.8) | 1.84 | 1.52 |

**Interpretation.**
- The α > β contrast survives on internal edges, which carry no tip correction. It shrinks somewhat, from a ratio of 1.83 to 1.59 (Δlog k from 0.61 to 0.46), but stays large and well resolved.
- Mean pLDDT is nearly the same in α and β HOGs: 80.5 vs 81.8. β HOGs have slightly *higher* pLDDT.
- Both classes have similar terminal shares of 3Di length (0.74 vs 0.79).
- The contrast also persists inside pLDDT-matched tertiles (internal ratio 1.5–1.7).
- So the tip / per-protein pLDDT artefact cannot explain most of the all-edge fold-class contrast. At most it inflates the all-edge ratio modestly.
- Caveats:
  - This is a descriptive contrast and is not adjusted for length or age stratum.
  - Internal k values are larger in both classes, so absolute k_α and k_β should be quoted with the edge type.

## 2. E6b (more-diverged-copy excess) adjusted for pLDDT

- Method: per-event excess, defined exactly as in production E6b (E1 reference curve at ĉ = 0.7226). It is regressed on:
  - Δ pLDDT = pLDDT(MDO) − pLDDT(LDO);
  - centred outgroup (O2) pLDDT;
  - in M3 only, centred mean copy pLDDT.
- Fit at the duplication-event level, n = 6,027.
- Uncertainty: HOG bootstrap, B = 1000, seed 20261030; 95% percentile CIs.

| model / stratum | n events | adjusted E6b (intercept or mean) | slope on Δ pLDDT (per pLDDT unit) |
|---|---|---|---|
| unadjusted (production E6b) | 6,027 | 0.0170 [0.0124, 0.0210] | — |
| M1: + Δ pLDDT | 6,027 | **0.0110 [0.0066, 0.0145]** | −0.0081 [−0.0090, −0.0072] |
| M2: + Δ pLDDT + outgroup pLDDT | 6,027 | **0.0110 [0.0066, 0.0145]** | −0.0081 [−0.0091, −0.0072]; outgroup −0.0003 [−0.0007, 0.0002] |
| M3: M2 + mean copy pLDDT | 6,027 | 0.0110 [0.0066, 0.0146] | −0.0081 |
| stratum \|Δ pLDDT\| < 2 | 2,780 | **0.0093 [0.0067, 0.0120]** | (mean Δ = −0.03) |
| stratum \|Δ pLDDT\| < 1 | 1,718 | 0.0090 [0.0061, 0.0121] | |
| stratum \|Δ pLDDT\| ≥ 2 | 3,247 | 0.0236 [0.0156, 0.0306] | |
| MDO has lower pLDDT | 3,351 | 0.0390 [0.0312, 0.0452] | |
| MDO has higher pLDDT | 2,676 | −0.0105 [−0.0144, −0.0062] | |

**Interpretation.**
- The MDO–LDO 3Di excess depends strongly on the pLDDT difference between the copies: about 0.008 of 3Di distance per pLDDT unit.
- The MDO has slightly lower pLDDT on average (mean Δ = −0.75; 56% of events are negative).
- Adjusting for Δ pLDDT lowers E6b by about a third (0.017 → 0.011). Restricting to copies with near-equal pLDDT gives 0.009.
- Both remain clearly positive, with CIs excluding 0. The direction of the positive control therefore survives pLDDT adjustment, but its magnitude is roughly halved.
- Caveat: lower pLDDT in the MDO could itself be a consequence of genuine structural divergence or degeneration (pLDDT as a mediator). If so, the adjusted value under-states the true excess. 0.009–0.011 is best read as a conservative lower range rather than a corrected estimate.

## 3. E1 curvature in high-confidence pairs

- Subset rule: pairs whose two proteins *both* have mean pLDDT ≥ 90 (or ≥ 85), with all production E1 filters re-applied to the subset.
- Uncertainty:
  - fixed-η estimators: production HOG × species pigeonhole bootstrap on the naive ĉ, B = 200, 95% CI;
  - SIMEX: no bootstrap. Its CI ≈ point ± 1.96 × the naive SE;
  - η_h estimator: HOG bootstrap, B = 200.
- The prot model was point-estimated only.

| | all (E1 data) | min pLDDT ≥ 85 | min pLDDT ≥ 90 |
|---|---|---|---|
| pairs / HOGs | 591,410 / 4,850 | 227,738 / 1,918 | 102,848 / 928 |
| mean min-pLDDT; median x | 82.5; 0.50 | 90.0; 0.55 | 93.0; 0.56 |
| **(a) fixed η = 0.0524 (pre-registered), naive** | 0.670 [0.621, 0.717] | 0.789 [0.749, 0.828] | **0.848 [0.812, 0.891]** |
| (a) fixed η = 0.0524, plan SIMEX (no D29 PB) | 0.717 | 0.864 | **0.943** (≈ ±0.04) |
| (a′) fixed η matched to subset floor (0.0524 / 0.0432 / 0.0398), naive / SIMEX | 0.670 / 0.717 | 0.731 / 0.796 | 0.764 / 0.842 |
| **(b) HOG baseline η_h** | 0.903 [0.881, 0.922] | 0.910 [0.880, 0.942] | **0.901 [0.865, 0.939]** |
| per-protein u_i + u_j (prot) | 0.728 | 0.741 | 0.726 |
| per-protein offset: mean u_i at c = 1 (median) | 0.045 (0.038) | 0.039 (0.034) | 0.035 (0.031) |
| u_i excess over η/2 = 0.0262 at c = 1 | 0.019 | 0.013 | **0.009** |
| u_i excess over the pLDDT-matched floor/2 at c = 1 | 0.019 | 0.017 | 0.015 |
| median η_h − η at c = 1 (HOG-baseline excess) | 0.032 | 0.020 | 0.013 |

**Interpretation.**
- The pre-specified prediction is met in direction. In high-pLDDT pairs the per-protein offset above the pre-registered floor shrinks (excess over η/2 falls from 0.019 to 0.009). The fixed-η ĉ moves towards 1 (naive 0.67 → 0.85; plan SIMEX 0.72 → 0.94). The HOG-baseline estimate stays at ≈ 0.90 in every subset.
- Two qualifications keep this short of a clean "the concavity is AF noise":
  1. Part of the rise in fixed-η ĉ comes from the pre-registered η (0.0524) being higher than the identical-sequence floor among high-pLDDT proteins (0.040). With a pLDDT-matched η, ĉ is 0.76 (naive) and 0.84 (SIMEX) at ≥ 90. That is still closer to 1 than the full data, but by about half as much.
  2. Measured against that matched floor, the per-protein excess falls only modestly (0.019 → 0.015 per protein). The per-protein model's ĉ also does not move (0.73 in every subset).
- So low-pLDDT proteins contribute to the fixed-η concavity, but a per-protein excess above the AF floor persists even at pLDDT ≥ 90.
- Overall, the pairwise data in high-confidence proteins remain compatible with c ≈ 0.75–1.0. This agrees with the earlier verification and does not settle linearity.
- Caveat: the high-pLDDT subset is a non-random set of HOGs (928 of 4,850), and its x range is similar to the full data.

## Limitations (all three checks)

- These analyses are post hoc, chosen after the artefact was identified, and are not multiplicity-controlled.
- The CIs are bootstrap percentile intervals with reduced B (200) for the curvature fits. The SIMEX points have no bootstrap.
- pLDDT is used as a covariate even though it may partly reflect real structural change. This applies especially to E6b.
