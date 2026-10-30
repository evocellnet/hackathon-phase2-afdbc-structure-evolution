# strucevo.stats: statistics engine report (statistical methods engineer, 2026-10-30)

**Bottom line.**
- `src/strucevo/stats/` implements ANALYSIS_PLAN §2.4–§2.6, §4 (E1–E6, including the E3 CATH contrast and E2 d\*) and §5 (Holm, §11.3 row ids).
- All 16 unit tests pass.
- These checks pass at the plan's settings: E1, E3, E4 and E6 recovery in V4, V3, and V2 unbiasedness.
- **Three verification checks fail as specified. Each has a proposed fix below, and each fix is tested.**
  1. V2 HC1 coverage is 0.85–0.88. Fix: **D2**, model-based SEs, which give 0.96.
  2. V4 E5 does not meet the N0 margins and the spike world is not detected. The cause is heavy-tailed lineage log-ratios. Fix: **D8**, a precision rule for lineages. With it, σ_L passes and the spike is detected. The β margin still fails, a power issue. See the E5 findings.
  3. V6 under the plan's node-age permutation covers 0 in 16/20. The permutation leaks age information through the tips. Fix: **D13**, an edge-level permutation, which gives 20/20.
- Runtime at production scale is minutes per endpoint (§4 below).
- Nothing under `results/endpoints/` was written. Synthetic outputs are in `results/verify/`.

## 1. API

```bash
source /opt/conda/etc/profile.d/conda.sh && conda activate <root>/conda_envs/main
python -m pytest tests/test_stats_engine.py -q        # 16 tests, ~5 s
python scripts/verify/v2_branch.py                     # ~1 min
python scripts/verify/v3_curvature.py [--kappa K --rho R --nu NU --tip_sd S]   # ~10.5 min
python scripts/verify/v4_endtoend.py [--min_t_snr 8] [--hogs N --obs counts|gauss]   # ~3.5 min (1000 HOGs, counts)
python scripts/verify/v6_shuffled_ages.py [--mode nodes|edges]                  # ~8 min
python scripts/verify/bench_stats.py                   # production-scale timing
```

| Module | Main entry points |
|---|---|
| `common` | `stage_rng(stage, *keys)` (SeedSequence(20261030), spawn key = config/seeds.yaml index); `load_pairs(aln, split)` and `load_branches`, which **raise `ConfirmSplitGuard` for split='confirm' unless `CONFIRM_RUN=1`** before reading any file; `guard_split(values)`, called by every engine on its input; `endpoint_record(**kw)` (the exact §4 key set comes first); `holm`; `normal_ci`/`pct_ci`; `write_json` |
| `trees` | `PanelTree.from_newick(path, ages={node:(med,lo,hi)})`; `prune(panel, S_h)` → `HogTree` (unary nodes suppressed, panel_child/parent, n_panel_edges, merged root D1); `design_matrix(ht, ia, ib)`; `nj_splits(D)`; `draw_ages(panel, rng)` (lognormal draw + top-down repair) |
| `branch` | `fit_branches(pairs, panel, metrics, nu, kappa, aln, hog_species, tau_misfit, se_method)` → `BranchFit`. It holds `.edges` (contract §1.3 columns plus `root_merged` and `share`), `.fits` (n_pairs_ok, n_edges, s2_hm, rss_rel, excluded, reason, fail_frac, misfit_ratio, var_sum, cov_sum_aa), `.pair_flags` (hidden paralogy) and `.hogs` (per-HOG X, A, scores). Also `wls_two_pass` (batched over metrics×alignments), `save_branch_tables` (writes `{aln}.parquet`, `{aln}_fit.parquet`, `{aln}_pairflags.parquet` under data/interim/branches) and `refit_with` |
| `pairwise` | `CurveEngine` (torch float64, GPU): profiled Q(c) for R replicate datasets at once; `two_pass`; `fit_batch` (grid → golden section, optional window with full-grid fallback); `simex` (λ grid, B_s, κx, κy, ρ, clip 1e-4, quadratic → λ=−1); `fit_eta_free`; `reference_fit` (independent numpy/scipy recomputation) |
| `bootstrap` | `multinomial_weights`, `Pigeonhole`, `pair_weights`, `edge_weights` (W_h × mean Λ over spanned panel edges), `replicate_groups` (new group ids for mixed models) |
| `reml` | `fast_mixed(y, x, z, groups)`: vectorised profiled REML for (1+x\|hog). Matches statsmodels to 1e-7 and is 15–40× faster |
| `endpoints_curve` | `run_curve_endpoint(tab, "E1"/"E2", eta, kx, ky, rho, excluded_hogs, B=200, B_s=50, B_s_boot=10, B_naive=1000, v3_false_concave, sensitivities, eta_raw, alignment_results)` → (record, curve{c, eta, a_h}, raw). Also `standard_sensitivities(...)`, `prepare_curve_data` (filter log → n_excluded), `curve_decision`, `d_star_boot` |
| `endpoints_branch` | `run_E3(bf, panel, tau_misfit, cath_class, hog_groups, B, alignment_k)`; `run_E4(bf, panel, tau_misfit, B, ages_override, ...)`; `build_cells` (§2.6 merge and censor rules); `fit_mixed` (statsmodels MixedLM, REML, unweighted, optimiser chain, random-intercept fallback flagged); `wls_sensitivity`; `pooled_slope` |
| `endpoints_lineage` | `run_E5(bf, panel, e4_rec, species_plddt, rho, nu_t, P, B, B_cv, spike_reps, guards={name: BranchFit}, e6_significant, min_t_snr)` → record with sigma_L, beta_contrast, power, floor_F, n_lineages, n_contrasts. Also `E5Data`, `parametric_floor`, `e5_core`, `spikein_power`, `e5_decide`, `e5_pvalues` |
| `endpoints_misc` | `run_E6(triplets, pairs, curve, B)`; `split_halves`; `apply_holm(dir)` (writes p_holm into E*.json, sets E4/E6 decisions, writes holm.json); `interpret`/`write_interpretation` (row ids only); `e4_robust_alpha_timescale` |
| `synth` | `random_panel`, `make_world` (N0, age effect, lineage spike), `observe_gauss` (ML-engine SE profile, κ, ρ, tip noise, optional power-law truth), `observe_counts` (simulate_counts with LG / Q.3Di.AF, d_max = 20), `make_triplets` (E6 N0 or duplication world) |

**E6 input schema** (the plan leaves it open). `triplets` has columns `hog, dcg_id, species, copy1, copy2, event_level, aa_pid_seq, aa_ml_seq, o1 (list of gene ids), o2 (list of gene ids)`. `pairs` has columns `gene_a, gene_b, status, aa_ml_hA, aa_ml_hB, tdi_ml_hA, tdi_ml_hB`. The upstream E6 pair builder must emit these.

## 2. Tests and verification

Each "Pass?" entry is judged against the plan's criterion.

| Check | Setting | Result | Pass? |
|---|---|---|---|
| Unit tests (`tests/test_stats_engine.py`) | confirm guard (load_pairs and engines), Holm, key set, seeds, pruning and merged root, exact WLS recovery including negative b, NJ = dendropy splits, tip correction, hidden-paralogy flag, GPU c vs numpy reference < 1e-6, batched/windowed = Brent, decision rules, pigeonhole weights, merge/censor rule, fast REML = statsmodels, split halves, interpretation rows | 16/16 | ✅ |
| **V2** branch engine (300 HOGs × 20 noise replicates, ρ = 0.3) | r(est, truth) | AA 0.994; 3Di 0.986 | AA ✅, 3Di ❌ (noise level, not bias) |
| | relative bias of Σb | AA 3e-4; 3Di 2e-4 | ✅ |
| | 95% coverage, **HC1 (plan)** | AA 0.884; 3Di 0.852 (κ = 1). At κ = 1.5 (SEs understated): 0.928 / 0.881 | ❌ |
| | 95% coverage, HC3 | 0.939 / 0.903 (κ = 1). 0.950 / 0.912 (κ = 1.5) | partial |
| | 95% coverage, model-based (XᵀWX)⁻¹ with w = 1/(SE²+s²) | **0.965 / 0.961** (κ = 1). 0.985 / 0.981 (κ = 1.5) | ✅ at calibrated κ → D2 |
| **V3** linear truth, 100 datasets × 200 HOGs (~18.5k pairs), κ = 1, ρ = 0.3, ν = 0.02, tip SD 0.01, B = 30, B_s = 20/5 | false-concave rate | **0 / 100** (all "linear"); mean ĉ_naive 0.972 → ĉ_SIMEX 1.004 (SD 0.005); 98% coverage of 1: 100% | ✅ |
| | power for "concave" | c = 0.97: 0/20; c = 0.94: 0/20; c = 0.90: 20/20. **Minimal detectable \|c−1\| ≈ 0.10** at this size | report |
| **V4** N0 world, counts (LG+G4 / Q.3Di.AF+G4), 150 species, 1000 HOGs, 118,870 pairs; plan rules | E1 98% CI covers 1 | ĉ = 0.995 [0.977, 1.012] (naive 0.936); "linear" | ✅ |
| | E3 (bonus) | k = 2.786 [2.701, 2.874], truth 2.755; σ_k 0.31 (truth ≈ 0.30) | ✅ |
| | E4 98% CI covers 0 | β = 0.0085 [−0.012, 0.029] | ✅ |
| | E6 (N0) CI covers 0 | E6a −0.003 [−0.010, 0.005]; E6b −0.001 [−0.011, 0.009] | ✅ |
| | E6 duplication world (bonus) | E6a 0.014, E6b 0.014, p = 9e-5 | ✅ detected |
| | E5 N0 margins (plan lineage rule) | β = 1.14 [0.90, 1.37]; σ²_obs 0.020 < F 0.058; σ_L upper 0.43 | ❌ |
| | Spike world detected (plan rule) | σ²_L lower99 −0.17 | ❌ |
| V4 with **D8** (`--min_t_snr 8`), same world | E5 σ_L 99% upper | **0.074** (< 0.10 ✅); F 0.008; σ²_obs 0.006 | σ ✅ |
| | β 98% CI | [0.89, 1.38] (HOG-only CI [0.99, 1.29]) | ❌ |
| | spike world | σ²_L 0.021, lower99 **0.0049 > 0** | ✅ detected |
| | spike-in power (procedure §4.5, 20 reps, B = 200) | f = 1.5: 0.35; f = 1: 0.00 | below 0.80 |
| V4 Gaussian, 5000 HOGs (~580k pairs), D8 at 8 | E1, E3, E4, E6 | all pass (ĉ 1.004 [0.998, 1.010]) | ✅ |
| | E5 | σ_L upper 0.057 ✅; β 0.955 [0.67, 1.24] (HOG-only [0.75, 1.16]) ❌; spike world detected ✅; power f = 1.5: 1.0 (10 reps); f = 1: 0.0 | β ❌ |
| V4 Gaussian, 5000 HOGs, D8 at 4 | E5 | σ_L upper 0.117 ❌; spike world not detected | ❌ |
| **V6** world with a true age effect (β_true(gen) = −0.15; fit −0.075 [−0.094, −0.055], detected), 20 permutations | **plan: internal-node ages permuted** | CI covers 0 in **16/20**; all 20 permuted β < 0 (mean −0.018) | ❌ (leak) |
| | **D13: edge midpoints τ permuted across all HOG edges** | **20/20** cover 0 (mean β −0.007) | ✅ |

Outputs:
- `results/verify/V2.json`, `V3.json` and `V4.json` (plan rules).
- `V4_snr4.json`, `V4_snr8.json`, `V4_gauss5000_snr4.json` and `V4_gauss5000_snr8.json`.
- `V6.json` (nodes) and `V6_edges.json`.
- `bench_stats.json`.
- `V4_endpoints/` holds E1/E3/E4/E5/E6.json and holm/interpretation from the last V4 run (gauss5000_snr8). These are synthetic and must not be pooled.
- Logs are in `logs/verify/`.

### E5 findings the PI should know about (they decide whether N0 can ever be accepted)
1. **The plan's lineage rule (≥ 50 HOGs and D_e ≥ 0.05) admits lineages whose 3Di sums are near zero.** These are young terminal lineages, for example 2–10 My cherries.
   - log R_e then has kurtosis 68, and lineages get censored at 0.005.
   - σ²_obs and F are each dominated by 2–3 lineages. In 3 of 10 F replicates σ² is ~0.5.
   - The result is σ²_obs < F, with SE(σ²_L) ≈ 0.05.
   - Fix D8 adds a selection rule that uses AA sums and SEs only: r̄·D_e/√ΣSE²_3Di ≥ 8. With it, kurtosis is 1.1 and the σ_L margin is met (0.057–0.074). The spike-in world is detected.
2. **β precision.** The 98% CI half-width is about 0.24 at 1000 HOGs and about 0.28 at 5000 HOGs with contrast-level resampling (0.20 HOG-only). That is too wide for ⊂ [0.8, 1.25] even under exact N0.
   - Precision depends on how much AA rates differ between sister lineages. The synthetic worlds used a lineage SD of 0.3.
   - **Measure it on the tuning set.** If it is similar, N0 acceptance is practically unreachable and E5 will be "inconclusive" whatever the truth. E5 *rejection* is still informative.
   - β̂ = 1.14 in the 1000-HOG counts world and 0.955 at 5000 HOGs. Upward bias from over-disattenuation is possible, because the pigeonhole variance of x_v is conservative. Check this with parametric N0 replicates on the tuning set.
3. **The V4 criterion "E5 accepts N0" cannot be met in an N0 world as written.** Acceptance also requires E6 to be significant, and E6 is null in that world (D18).

## 3. Runtime at production scale (A100 + 48 cores; `bench_stats.json`)

The benchmark used 5600 synthetic HOGs with 8–30 species (697k pairs; 582k E1 pairs). Rows marked "scaled" are extrapolated linearly in pairs to about 2.5M pairs.

| Step | Measured | Production estimate |
|---|---|---|
| Branch engine: 5600 HOGs × 3 metrics, incl. NJ misfit and NNLS | 19 s | ~1–2 min per alignment (×3 alignments) |
| E1/E2 point + SIMEX (B_s = 50) | 5.3 s | ~25 s (scaled) |
| E1/E2 bootstrap: 200 × (naive + SIMEX B_s = 10) + B_naive = 1000 | 1.17 s/replicate | **~4.5 min (measured) → ~18 min (scaled)** per metric × alignment |
| E3, B = 1000 | 1 s per 100 replicates | ~10 s |
| E4 point fit (statsmodels, 12.6k cells) | 14.5 s | ~30 s |
| E4 bootstrap with fast REML (0.93 s per replicate + cells) | — | **B_E4 = 1000 ≈ 16 min, so B_E4 = 1000 fits the 1 h rule.** statsmodels alone would take ~4 h |
| E5: floor P = 100 + B = 1000 + B_cv = 200 + spike-in 4 f × 50 × B = 200 | — | ~4 min |
| E6, B = 1000 | — | seconds |
| V4 end-to-end, 1000 HOGs, counts | 204 s | — |

## 4. Proposed DEVIATIONS and clarifications (for PI approval)

Items marked **D** change a plan rule. Items marked **C** fill a gap the plan leaves open; the implementation states its choice.

- **D1. Bifurcating HOG root.** Pairwise distances identify only the sum of the two root-adjacent edges.
  - These two edges share one WLS parameter. It is apportioned by panel time span (t_root − t_child), or 50/50 when ages are missing.
  - Such edges are flagged `root_merged` and have `panel_edge = NA`, so they are never E5 lineages.
  - E5 contrasts whose HOG node is a merged root are skipped.
- **D2. Per-edge SE.** HC1 under-covers (0.85–0.88) because of high-leverage tree designs. Propose `b_se` = model-based (XᵀWX)⁻¹ with w = 1/(SE_cal² + s²_hm), which gives 0.96.
  - It is conservative (0.98) if κ is under-calibrated.
  - `cov_aa` keeps the HC1 residual cross-product, scaled by √(hc1_AA·hc1_m), because a model-based cross-covariance would need ρ.
  - `se_method='hc1'|'hc3'|'model'` is selectable. **Primary inference is the pigeonhole bootstrap in every endpoint**, so per-edge SEs only feed E3 v_h, the E4 WLS sensitivity and nothing decisive.
- **D4. NJ.** NJ uses an in-house numpy implementation (identical splits to dendropy `nj_tree`, tested) for speed. AA pairs that failed are filled with the species-tree fitted values for the NJ input only.
- **D5. Weights in replicates.** In SIMEX and bootstrap replicates the pass-2 weights, including s², are held at their full-data values; bootstrap multiplies them by the pigeonhole pair weights.
  - Replicate fits use the grid restricted to ĉ ± 0.3, falling back to the full grid if the minimum is on the window edge, then a vectorised golden section. This equals Brent to 1e-6.
  - The point estimate uses the full grid plus scipy Brent.
- **D6. Fast REML in E4 bootstrap replicates.** The E4 point estimate and the sensitivity point estimates use statsmodels MixedLM; bootstrap replicates use `reml.fast_mixed` (Δβ ≈ 1e-7).
  - In N0 worlds statsmodels flags `re_slope_fallback`, because the slope variance sits at its boundary of 0. This is expected and is flagged as the plan requires.
- **D8 (important). E5 lineage precision rule.** Add r̄·D_e/√ΣSE²_3Di ≥ 8 (see §2). The threshold should be confirmed on the tuning set by spike-in power and the null rate.
- **C9. E5 resampling.**
  - Λ (lineage) weights cancel inside R_e, so they are used as weights in the across-lineage OLS and variance.
  - The β CI adds contrast-level multinomial weights. A HOG-only CI is also reported. **The PI should pick one.**
- **C10. E5 floor F.** The refit is exact-linear (b\* = A d\*) with fixed pass-2 weights.
  - r_h is the HOG ratio of tip-corrected sums, and x_e is log τ of the HOG edge.
  - Simulated 3Di distances include η = 2ν, so tip correction is consistent.
  - E4 fixed effects are held fixed in F, in the spike-in and in the guards.
- **C11. Spike-in implementation.** b_3Di is multiplied by f on simple edges of the chosen lineages. This is exactly the plan's "add (f − 1)·b̂ to crossing pairs and refit", because AX = I with fixed weights.
  - F, c_v and vx_v are held fixed, and E4 is not refitted per replicate.
  - 10% means 10% of the eligible E5 lineages.
- **C12. E4 cells.**
  - If merging leaves a single cell with 0 < S_AA < 0.05, the cell is kept.
  - Cells with S_AA ≤ 0 are dropped and counted (`cells_aa_nonpositive`).
  - Λ̄ edge weights also weight n_edges and the mean τ.
  - x is re-centred within each replicate; β does not depend on the centre.
- **D13. E4 negative control.** The plan's node-age permutation leaks age information: tips stay at 0, so terminal edges remain "young". In a world with a real age effect it fails (16/20).
  - Propose permuting edge midpoints τ across all HOG edges (`run_E4(tau_perm=...)`), which passes 20/20.
  - A small residual negative mean (−0.007) remains, likely a log-ratio/size effect, and should be watched on the tuning set.
- **C14. E1 sensitivity "failures imputed at the cap".** Implemented as y = 20 (d_max) with SE = the median SE of eligible pairs. The plan does not define the cap for this purpose. Imputing 20 inside a d ≤ 2 power fit is extreme, so **please specify**; a likely choice is the saturation threshold or the cap from §4.4b.
- **C15. E2 d\*.**
  - A and A′ are weighted medians with weights equal to the per-HOG information Σw x² (respectively Σw x^{2c′}), with w = 1/SE_y².
  - η_raw (the floor for 1 − LDDT_global) is an input.
  - d\* = NaN if no grid point exceeds 10%.
  - The CI uses 200 pigeonhole replicates.
- **C16. E3.**
  - v_h is computed by the delta method on HOG totals using the full sandwich (1ᵀV1), not the diagonal sum.
  - σ²_k uses HOGs with both sums > 0.
  - The k CI is computed on log k and back-transformed.
  - CATH classes are `alpha`/`beta` from the caller's map.
- **C17. E5 p-values** are computed on the σ² scale (σ_L = 0.10 ⇔ σ²_L = 0.01). The guards are internal-edges-only (built in) plus BranchFit variants supplied by the caller (hidden paralogy removed, single-copy HOGs).
- **D18. V4 E5 criterion.** In the N0 world V4 is judged on the E5 margins (β CI and σ_L bound) and on spike detection. E6 positivity is verified separately in a duplication world.
- **C19. Interpretation.**
  - The E1/E2 column uses the E1 class, unless E2 contradicts it; then the case is "discordant" and no curvature row matches.
  - All matching §11.3 rows are written, because the rows are not mutually exclusive, plus a `duplication` modifier.
  - E4 "≈ 0" means p_holm ≥ 0.05.
- **C20. E6.**
  - P1 is chosen by a seeded coin flip (stage `e6_p1_choice`).
  - O pairs that are missing are skipped.
  - MDO is chosen by mean `aa_ml_hA` to O1.
  - `replication_direction` reports the observed signs. PMID 40813247's direction must be filled in by the PI or a reviewer; I did not look it up.
- **C21. s² in the pairwise engine** uses unweighted means over eligible pairs, as written. The two-pass s² in the branch engine is per HOG × metric and uses the plan's formula as written. That formula is dimensionally odd (RSS_w/(n−p) is already a variance ratio) but harmless.

## 5. Not done / next steps
- No real data were touched. The E1/E3/E4 sensitivity *inputs* that need recomputed distances (α ×1.5 and ÷1.5, all columns, deep ages, other alignments) are accepted as variant tables, BranchFits or `ages_override`, but are not generated here.
- Rerun V3 with tuned κ, ρ, ν and tip SD on real HOG designs: `v3_curvature.py --kappa .. --rho ..`. Extend `--power_c` to c ∈ {0.92, 0.96} to refine the minimal detectable |c − 1|. Use B = 200 for the final run.
- V3's minimal detectable |β_E4| is not computed yet. The E4 bootstrap SE at 5000 HOGs is about 0.004 (CI half-width 0.009), which suggests |β| ≈ 0.01–0.015 is detectable.
- `scripts/endpoints/run_confirm.py` and `interpret.py` (CLI wrappers around the functions above) are not written.

## Files
- `src/strucevo/stats/{__init__,common,trees,branch,pairwise,bootstrap,reml,endpoints_curve,endpoints_branch,endpoints_lineage,endpoints_misc,synth}.py`
- `tests/test_stats_engine.py`
- `scripts/verify/{v2_branch,v3_curvature,v4_endtoend,v6_shuffled_ages,bench_stats}.py`
- `results/verify/*.json`, `logs/verify/*.log`
