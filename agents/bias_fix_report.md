# Bias fix report: E1 (SIMEX curvature) and E4 (log ratio vs log age): statistician-engineer, 2026-10-30, 16:25–17:45

## Bottom line

**E1 (problem A).**
- **Cause.** The bias comes from the *size* of the AA measurement error. On real euk20 designs (n_core median 269; κ_AA ≈ 1.4) the naive ĉ is attenuated to 0.82, and the plan's quadratic SIMEX extrapolant then stops short.
  - Under a linear truth ĉ is 0.955–0.966.
  - With the oracle per-pair Σ_ij it is still 0.955.
  - With the AA error halved it is 1.002.
- **Ruled out.** Scalar κ/ρ, y-dependent weights, truncation, the x clip and the window proxy each move ĉ by ≤ 0.015.
- **Fix (D29).** Keep the plan's SIMEX on raw x, with per-pair Σ_ij via ρ(bin). Then subtract a **parametric-bootstrap estimate of the residual SIMEX bias**: 25 pseudo-worlds from the fitted model on the tree-smoothed AA design, each re-estimated with the identical estimator.
- **Production V3** (euk20 designs, tuned noise, edge-level error, 50 linear datasets):

  | Truth | Fixed ĉ | Plan ĉ (no correction) | Call rate |
  |---|---|---|---|
  | linear | **0.993** | 0.949 | false-concave **0/50** |
  | c = 0.8 | 0.788 | — | concave 100% |
  | c = 0.9 | 0.886 | — | concave 87% |

- **Limitation.** The correction is accurate near the calibrated noise level. Where the plan was already almost unbiased (lower noise; the stats-engine synthetic V3 world) it over-corrects by +0.01 to +0.025, which is the safe direction: false-concave stays at 0. At twice the calibrated variance it under-corrects (0.941).

**E4 (problem B).**
- **Cause.** Per-protein 3Di tip noise on terminal edges: it cannot be removed, and the branch-engine SEs cannot see it. In small young cells it produces censoring and Jensen bias. A ν error amplifies it.
- **Fix (D30).** An AA-only minimum-information merge rule.
- **Result** (50 null worlds, B = 100):

  | | Mean β̂ | Mean \|β̂\|/SE | False-positive rate (95%) |
  |---|---|---|---|
  | Plan rule | +0.058 | 1.60 | 20% |
  | D30 | **−0.0075** | **0.37** | **0%** |

- **True effects are recovered better:** −0.1 → −0.100 (plan cells: −0.050) and −0.2 → −0.190.
- **Cost:** the V6 permutation control in a world with a *true* effect fails, 11/20 (leak −0.013). PI decision needed (§4.1).

**Rejected along the way.**
- Tree-smoothed regressor x̂ (perfect with independent pair errors, but 0.900 with edge-level errors).
- Delta-corrected log ratio (unstable).

**Real-data structure finding.** The residual of aa_ml about the species-tree fit has only 13–21% of κ²SE². Real pair errors are mostly tree-additive, while every synthetic world used before today drew independent pair errors. A new edge-noise world was added to V3.

## 1. Diagnosis

Harnesses:
- `scripts/verify/{bias_common,e1_bias_diag,e4_bias_diag,v6_leak_diag,diag_summary}.py`
- Raw results (resumable): `data/interim/itest_euk20/results/tuning/bias_e{1,4}_diag.jsonl`
- Summarise with `python scripts/verify/diag_summary.py <jsonl> c_naive,c_simex,c_pb`

Designs:
- euk20 TEST root, split == 'tune': 96 real HOG designs (species sets, n_core, NNLS AA profiles, r_h).
- **No confirm data were read.**

### 1A. E1: linear truth; 32 (PB rows 16–24) datasets × 384 designs; SE of a mean ≈ 0.003

"base" is production-like:
- noise κ(bin), ρ(bin);
- Fisher SE columns at the **observed** distance;
- κ(bin of observed aa_ml) applied to the SE columns, as the driver does;
- scalar ρ; w = 1/(SE_y² + s²).

| Cause switched on/off | ĉ_naive | ĉ (SIMEX or fix) | Contribution |
|---|---|---|---|
| **base = plan estimator** | 0.823 | **0.966** | **−0.034 total** |
| integration agent's setting (SE at truth, scalar κ) | 0.837 | 0.958 | — |
| scalar-noise world + scalar estimator | 0.834 | 0.957 | κ/ρ scalar vs bin ≈ 0 |
| **oracle Σ_ij** (true per-pair κ, ρ, SE) | 0.840 | **0.955** | Σ misspecification is not the cause |
| per-pair ρ(bin) | 0.823 | 0.953 | −0.013 |
| weights from an x-only model / flat weights | 0.805 / 0.776 | 0.945 / 0.945 | −0.021: y-informative weights are not the cause |
| drop clipped SIMEX draws instead of clip (0.4% of draws) | 0.823 | 0.967 | +0.001 |
| aa_ml_seq window proxy shares the aa_ml error (worst case) | 0.829 | 0.977 | +0.012 |
| truncation (fail/saturation) | — | — | 0: no failures in the designs, no saturation at d ≤ 2 |
| η ± 0.01 (sensitivity; η is exact in the simulation) | 0.865 / 0.785 | 1.017 / 0.922 | **±0.05 per ±0.01** |
| **AA + 3Di SE halved / quartered** | 0.935 / 0.980 | **1.002 / 1.001** | **error magnitude = cause** |
| other extrapolants (cubic; quadratic in 1/c or log c) | — | 0.976–1.018 | config-dependent; not robust |
| **edge-noise world**, plan estimator (ρ scalar / pair) | 0.821 | 0.967 / 0.954 | plan bias is the same as with independent noise |
| rejected: tree-x̂ regressor + exact-covariance SIMEX, independent / **edge** world | 0.931 / 0.861 | 1.003 / **0.900** | ✗ |
| PB, design D* = observed x: pair / edge | — | 0.989 / 0.990 | under-corrects (c = 0.8: 0.789 / 0.785) |
| **PB, design D* = tree x̂ (adopted): pair / edge** | — | **1.010 / 0.998** | |
| … c = 0.8 (pair / edge); c = 0.9 (edge) | — | 0.804 / 0.790; 0.894 | |
| … edge world, aa_ml_seq shares the error / SE at truth (pair) | — | 0.995 / 0.995 | |
| … variance ÷ 2 (n_core × 2): pair / edge | 0.890 / 0.893 | plan 0.987 / 0.993 → PB 1.010 / 1.012 | over-corrects when the plan bias is small |
| … variance × 2 (n_core ÷ 2), edge | 0.725 | plan 0.872 → PB **0.941** | under-corrects at double noise |

**Mechanism.**
- The relative SD of x is 15–30%, largest at small x and on short cores. The SIMEX path c(λ) = 0.82, 0.76, 0.72, 0.69, 0.66 (λ = 0…2) is convex, so the quadratic extrapolation to λ = −1 stops 0.035–0.045 short.
- PB re-runs the identical estimator on pseudo-worlds whose truth is known (c = ĉ_SIMEX).
- The tree-smoothed design avoids the excess spread an observed-x design would add.
- The tree-x̂ *regressor* fails in edge-noise worlds because edge-level (shared) error passes through tree smoothing, while V_h assumes independent pair errors.

**Real-data structure check** (AA only, euk20 tune; used to choose the simulation model, never the rule):
- The residual variance of aa_ml about the species-tree WLS fit is **0.13–0.21 × κ²SE²(1 − h)** in every distance bin.
- So realised substitutions per edge, shared by every pair through that edge, dominate pair-level error.
- New edge-noise world: edge noise variance κ²·SE(b_e; L_h)², plus the non-additive pair remainder. Its residual ratio is 0.48, still less edge-dominated than real data.

### 1B. E4: null worlds; 24 × 400 HOGs on the real panel and tune species sets; tuned noise; point estimates

| Cause switched on/off | β̂ | SD | Censored | Reading |
|---|---|---|---|---|
| **base (plan cells)** | **+0.059** | 0.021 | 4.6% (> 2% flag) | |
| tip noise off | +0.020 | 0.016 | 1.9% | **tip noise ≈ +0.04** |
| measurement noise × 0.5 / × 0.25 | +0.062 / +0.063 | | 3.9% / 3.7% | alone ≈ 0 |
| measurement × 0.25 and tip off | +0.003 | 0.004 | 0.02% | remainder: measurement noise in small cells |
| merge threshold S_AA ≥ 0.2 | +0.031 | 0.019 | 0.9% | partial |
| covariate z off | +0.075 | 0.024 | | z is not the cause |
| drop censored cells | −0.007 | 0.018 | 0 | censoring + Jensen is the mechanism |
| ν + 0.005 / ν − 0.005 | **+0.141 / −0.016** | | 6.2% / 3.4% | extreme ν sensitivity (plan) |
| delta-corrected log ratio | unstable | | | rejected |
| edge-noise world, plan cells | +0.050 | 0.026 | 4.7% | |
| **D30, T = 3 / 4 / 6** | −0.009 / **−0.007** / −0.016 | 0.015 / 0.016 / 0.059 | 0.1% / 0.04% / 0 | T = 4 adopted |
| D30 (T = 4), edge-noise world | +0.012 | 0.017 | 0.03% | |
| D30 + edge-noise variance in V_t (edge / pair) | +0.015 / +0.014 | | | not adopted |
| D30, z off | −0.002 | 0.016 | | |
| D30, ν ± 0.005 | +0.031 / −0.042 | | | reduced but present |
| true effect −0.1: plan / D30 | −0.050 / **−0.093** | | | the plan cells also attenuate real effects |

**Mechanism.**
- Terminal 3Di edges absorb each protein's offset u_i − ν (SD 0.028). It is tree-additive, so it never appears in the SEs or in s².
- Young cells hold S_3Di ≈ 0.015–0.03 against an SD of about 0.03–0.04. That gives censoring (S_3Di ≤ 0 is set to 0.005) and E log S ≪ log E S.
- Young cells are therefore biased low, which gives β̂ > 0.
- D30 merges (plan direction) while r̄·S_AA < 4·√V_t, with V_t = Σ(b_se_3Di² + tip_sd²·terminal). The rule uses AA sums, 3Di SEs and the tip variance, never the cell's 3Di values: the D8 precedent.

## 2. Implementation

| File | Change |
|---|---|
| `src/strucevo/stats/pairwise.py` | `simex`: per-pair kx, ky, ρ arrays; `clip_mode` (diagnostic; the default is still the plan's clip). **New `simex_pb_correction`** (D29). `TreeX`/`xtree` (rejected regressor) is kept but unused by default. |
| `src/strucevo/stats/endpoints_curve.py` | `tree_x()` gives x̂ = X_h b̂_AA (PB design only). `fit_curve_full(..., pb_K, pb_Bs, x_design)`: c = c_SIMEX − bias; bootstrap reps are shifted by the fixed correction; SE = √(SE_boot² + MC²). `run_curve_endpoint(..., bf=, rho_pair=, pb_K=25, pb_Bs=20, pb_K_sens=10)`. New record keys: `c_simex_plan`, `pb_bias`, `pb_mc_se`, `pb_K`, `pb_design`, `rho_mean`. **`pb_K=0` = plan.** |
| `src/strucevo/stats/endpoints_branch.py` | `build_cells(rule, info_T, tip_var, r_bar, tip_all)` now also returns V_aa and V_t. **`run_E4(cell_rule='info', info_T=4.0, tip_var=…)` is the new default.** Sensitivity `cells_plan_rule` when `sensitivities=True`. |
| `src/strucevo/stats/driver.py` | `Tuning.tip_sd(aln)`: yaml `tip_sd.{aln}`, otherwise the floors.json floor_ident (q95 − q05)/3.29/√2. E1/E2 (primary and alternative alignments) pass `bf` (PB design) and per-pair `rho_pair`; C14 reuses `bf`; the α variants use the x design. Every `run_E4` passes `tip_var`. |
| `scripts/tuning/v3_tuned.py` (production V3) | New defaults: `--pb-K 25 --pb-design tree --rho-mode pair --se-at obs --noise-struct edge --e4-cell-rule info`; the plan versions remain as options. Reports the plan's uncorrected mean, the 95% FPR and \|β̂\|/SE. **At production tune size use `--pb-K 10`**: PB costs about 0.45 s per pseudo-fit per 20k pairs, and K = 25 would add 2–3 h to V3. |
| `scripts/verify/v3_curvature.py`, `v4_endtoend.py`, `v6_shuffled_ages.py` | D29/D30 by default (`--pb-K 0` / `--cell-rule aa` give the plan). V3 recomputes SEs at the observed distances. |
| `scripts/tuning/b_e4_timing.py`, `e5_power.py` | Pass `tip_var`. |
| `tests/test_stats_engine.py` | +2 tests (per-pair SIMEX arrays equal scalars; info-rule merge/keep). **18/18 pass.** |

The driver was checked end to end by the `--quick` dry run on the test root's tune split (`results/tuning/dryrun_biasfix/`): E1–E6 run, and the new keys are filled. The E2 PB correction there is −0.09; see §4.5.

## 3. Validation

| Check | Setting | Result | Target |
|---|---|---|---|
| **V3 production** (`v3_tuned`; test root `results/tuning/v3_tuned_pb_e1.json`) | 50 linear datasets × 384 euk20 designs, edge noise, B = 50, B_s = 20/5, K = 20 | ĉ **0.993** (SD 0.021, SE 0.042); plan uncorrected 0.949; **false-concave 0/50**; coverage of 1: 100% | ≤ 10%, \|bias\| < 0.01 ✅ |
| | c = 0.8 / 0.9 (15 each) | ĉ 0.788 / 0.886; concave **100% / 87%**; MDE \|c − 1\| = 0.10 | detects c = 0.8 ✅ |
| V3 synthetic (`results/verify/V3.json`; random panel, independent noise, κ 1.3, ρ 0.28, 40 × 300 HOGs) | linear | ĉ **1.025** (plan 1.009; SD 0.008); false-concave 0/40; 39 linear | FCR ✅; bias +0.025 ✗ (over-correction, safe side) |
| | c = 0.9 / 0.8 (15 each) | 0.924 / 0.822; concave 93% / 100% | ✅ |
| V4 counts world (`results/verify/V4_snr8.json`, D8 at 8) | N0 | E1 ĉ 1.003 [0.988, 1.017] "linear" (plan 0.994; PB bias −0.009); **E4 β 0.0086 [−0.005, 0.022]**, 1628 cells, 0 censored; E3/E6/spike as before | E1, E4 ✅; E5 margins ✗ (pre-existing power issue, D8 note) |
| **E4 null** (`v3_tuned`; `v3_tuned_fixed_e4.json`) | 50 null worlds × 400 HOGs, B = 100 | β̂ **−0.0075**, SE 0.044, mean \|β̂\|/SE **0.37**, FPR **0% / 0%** (95% / 98%) | \|β̂\| < 1 SE ✅; FPR ≤ 5% (bootstrap SE ≈ 2 × SD, conservative) |
| E4 null, plan cells (`v3_tuned_planrule_e4.json`) | same | +0.058, \|β̂\|/SE 1.60, FPR 20% / 8% | (before) |
| **E4 true effect** | −0.1 / −0.2, 50 worlds each | β̂ −0.100 / −0.190; detected 54% / 100% (95%), 28% / 84% (98%) | ✅ |
| E4 edge-noise null (point estimates, 24) | | +0.012 (SD 0.017) | ≈ 0.3 SE ✅ |
| **V6** (`results/verify/V6_edges.json`; D13, world with true effect −0.15) | D30 | true ages β −0.077, detected; permutations cover 0 in **11/20** (mean −0.013, SE 0.006) | ✗ (plan cells: 20/20, mean −0.0075, SE 0.010) |

## 4. PI decisions and caveats

1. **V6/D13 control under D30 (11/20).**
   - In a world with a *true* age effect, any merge rule that depends on edge properties leaks true ratio heterogeneity into permuted-age fits (`scripts/verify/v6_leak_diag.py`): plan −0.0075, D30 −0.013 (−0.0129 with tip variance on all edges, −0.012 at T = 3).
   - In a null world nothing leaks. At production SEs the plan rule would also fail this control.
   - **Options:**
     - (a) Keep D30 and judge the control by \|mean permuted β̂\| ≤ 0.2·\|β̂\| (C32 below).
     - (b) Keep the plan's cells and accept null bias +0.058, FPR 20%.
   - **I recommend (a).**
2. **E1 η sensitivity.** ±0.01 in η moves ĉ by ±0.05 under every estimator. η = 0.052 comes from identical-sequence pairs. Read "η free" as a primary robustness check.
3. **E4 ν sensitivity.** Under D30, ±0.005 in ν moves β̂ by +0.03 / −0.04. Recommend a new sensitivity with ν ± 0.005 (refit branches). The internal-edges-only sensitivity, which is free of tip noise, is the key check. Not implemented.
4. **D29 relies on the calibrated noise magnitude.**
   - Its error is ±0.01 near the calibrated level, +0.01 to +0.025 at lower noise and −0.06 at double the variance.
   - The plan's V3 guard (false-concave rate > 10% → inconclusive) stays in force, run with the D29 estimator on the real tune designs.
5. **E2** uses the same PB code, but it was **not separately validated** for lddt_D noise. In the euk20 dry run its correction (−0.09) is three times E1's.
6. The rejected code paths (tree-x̂ regressor `xtree`, `clip_mode`, `tip_all`) remain in place but are off by default.
7. PB draws from `stage_rng(seed_stage, *extra, 2)`, a sub-stream of the endpoint's own stage, so `config/seeds.yaml` is unchanged.
8. Nothing was committed. DEVIATIONS.md was not edited (rows below).

## 5. DEVIATIONS rows to add (pre-freeze; synthetic verification only; no confirm data; no endpoint value used for any choice)

| ID | Change | Reason | Effect |
|---|---|---|---|
| D29 | E1/E2: after the plan's SIMEX (per-pair Σ_ij: κ(bin) on the SE columns, ρ(bin) per pair), subtract a parametric-bootstrap estimate of its residual bias. K = 25 pseudo data sets are drawn from the fitted model (c = ĉ_SIMEX, a_h, η, s²) on the design D* = tree-smoothed AA distance X_h b̂_AA (branch-engine AA fit; observed aa_ml where there is no fit), with the same calibrated noise. Each is re-estimated by the identical estimator (two-pass + SIMEX, B_s = 20). ĉ = ĉ_SIMEX − (mean ĉ* − ĉ_SIMEX). The correction is held fixed in bootstrap replicates; SE = √(SE_boot² + MC²). The plan value is reported as `c_simex_plan`; sensitivities use K = 10. | V3 on euk20 tune designs (tuned κ(bin)/ρ(bin), real n_core): the plan's SIMEX gives ĉ = 0.95–0.97 under a linear truth (naive 0.82), so E1 would be called "concave" at production precision. The cause is the quadratic extrapolant under large AA error (oracle Σ 0.955; halved SE 1.002); κ/ρ form, weights, clipping, truncation and the window proxy are ≤ 0.015 each. With D29: 0.993 (production V3, edge-level error), 1.010 / 0.998 (independent / edge diagnostics); c = 0.8 → 0.79–0.80; false-concave 0/50. | estimand and regressor unchanged; removes the concave bias at the calibrated noise level (slight over-correction at lower noise); about 1–2 min per fit |
| D30 | E4 cells (§2.6): in addition to the AA rule, merge a cell (same direction rule) while r̄·S_AA < 4·√V_t, with V_t = Σ_edges(b_se_3Di² + tip_sd²·[terminal edge]). tip_sd = floor_ident tdi_ml (q95 − q05)/3.29/√2 (or tuning.yaml `tip_sd.{aln}`); r̄ = global Σb_3Di/Σb_AA (D8 precedent; no cell-level 3Di values). The plan rule is reported as sensitivity `cells_plan_rule`. | Null worlds (euk20 panel and species sets, tuned noise, 50 × 400 HOGs): the plan's cells give β̂ = +0.058, \|β̂\|/SE 1.6, 95% FPR 20%. Cause: per-protein 3Di tip noise on terminal edges, invisible to the SEs, which leads to censoring (4.6%) and Jensen bias in young cells. D30: β̂ −0.0075, \|β̂\|/SE 0.37, FPR 0%; true effects −0.1 / −0.2 are recovered as −0.100 / −0.190 (plan: −0.050 for −0.1). | about 40% as many cells; < 0.1% censored; ν sensitivity per ±0.005 drops from +0.14/−0.02 to +0.03/−0.04; the D13 control leaks more in worlds with a true effect (C32) |
| C31 | The V3 tuning simulation (`v3_tuned.py`) draws edge-level AA/3Di error (realised substitutions per HOG edge, variance κ²·SE(b_e)², plus the non-additive pair remainder) and reports Fisher SEs at the observed distances | On the tune split the aa_ml residual variance about the species-tree fit is 0.13–0.21 of κ²SE², so pair errors are mostly tree-additive; the ML engine reports SEs at the observed distance | test definition only |
| C32 *(only with PI option 4.1a)* | D13/V6 negative control in worlds with a true age effect passes if \|mean permuted β̂\| ≤ 0.2·\|β̂(true ages)\|, replacing "98% CI covers 0 in ≥ 18/20" | Merge rules that depend on edge properties leak true ratio heterogeneity into permuted fits (plan −0.0075, D30 −0.013, at β̂ = −0.077); smaller SEs make the coverage rule fail regardless of rule | test definition only |

## Progress log
- 16:25 start. 16:30–16:45 E1/E4 cause-switching diagnostics. 16:45 first E1 fix (tree-x̂ regressor) and D30 implemented; tests 18/18.
- 16:58 Quick dry run on euk20 tune: tree-x̂ and raw-x disagree in the opposite direction to the simulations. An AA-only check shows tree-additive (edge-level) error dominates the real data. 17:05 edge-noise world: the tree-x̂ regressor gives 0.900, so it is **rejected**.
- 17:12 PB correction implemented (tree design). 17:14–17:37 validations: v3_tuned E1 (PB, edge) and E4 (D30 + plan baseline); V3 synthetic; V4; V6; leak diagnostic; dry run.
