# Review round 2: statistics (DESIGN.md v2)

Statistics agent · 2026-10-30. **Verdict:** right architecture, but the new tree-fit machinery adds four routes to spurious "saturation" or age trends: NNLS truncation, tip noise, branch-level power laws, and ratios of short branches. The fixes are estimator definitions only.

## 1. Round-1 items

| Item | Status | Remaining gap |
|---|---|---|
| S1 selection filter | Partial | Coverage rule still drops pairs (more at depth); handling in tree fits undefined (M8) |
| S2 errors-in-variables | Partial | Per-branch SEs and the SIMEX recipe undefined; x/y errors correlated (M1, M3) |
| S3 composition and dependence | Mostly | No within-HOG spread requirement; bootstrap not implementable as written (M6) |
| S4 decision rule | Mostly | No noise intercept; no "inconclusive" outcome (M3, M7) |
| S5 ratio estimator | Yes | E3 uses AA/3Di, E4/E5 use 3Di/AA |
| S6 RRT and dispersion | Partial | E5 has no margins and no disattenuation (M5) |
| S7 matching | Yes | Null simulation under a constant ratio not stated |
| S8 TimeTree circularity | Yes | — |
| S9 H1e | No | Post-2018-04-30 PDB entries; same-protein AF-vs-PDB floor |
| S10 power | No | Minimal detectable c and E4 slope from the S2 simulation |
| S11 pre-registration | Partial | Holm applied to endpoints with no p-value (M7) |

## 2. The new machinery

**NNLS and weights.** Fitch–Margoliash weights 1/d² understate the variance of ML distances near saturation, over-weighting deep pairs. Truncation at 0 biases short branches *upward*, more for the noisier metric (3Di/3D). So b_3Di/b_AA is inflated on short, young branches: a spurious negative E4 slope and spurious E1 concavity. Short-branch log-ratios are heavy-tailed. A minimum-length cut selects on noisy b̂_AA (regression to the mean) and correlates with age.

**Tip noise.** Per-protein noise (AF prediction) adds a constant to every distance involving that protein; least squares puts it on the terminal branch. Terminal branches are young, so this again manufactures an age trend.

**Branch-level curvature is not an estimand.** If b_3Di = a·b_AA^c held per branch with c ≠ 1, path sums would give aΣb^c ≠ a(Σb)^c, and adding a taxon (splitting a branch) would change c. Curvature belongs to the pairwise (path) map. Fitted patristic distances fail too: a concave map makes distances non-additive. Hence observed pairs (M3). For the same reason, the §3 non-additivity exclusion would remove exactly the saturating metrics: report it, do not gate on it.

## 3. Remaining changes

**M1 (must). Branch estimator.** Use unconstrained weighted least squares (WLS):

b̂ = (XᵀWX)⁻¹XᵀWd, with w_ij = 1/(SE²_ij + s²_m)

- SE²_ij from Fisher information (AA/3Di) or column bootstrap (3D); s²_m a misfit variance tuned on the 20% set.
- Per-branch SEs from the sandwich (XᵀWX)⁻¹XᵀWΣ̂WX(XᵀWX)⁻¹ (Σ̂ diagonal, plus AA–3Di within-pair covariance).
- Keep negative estimates for inference; NNLS for display only.
- Ratio estimands are ratios of *sums* over cells, never per-branch ratios.

**M2 (must). Tip noise.** Subtract ν_m = ½ × median identical-sequence pair distance (§4.4b) from each terminal branch; internal-branches-only sensitivity for E3–E5.

**M3 (must). Redefine E1/E2 on observed pairwise distances within HOGs:**

D_m,ij = η_m + a_h·D_AA,ij^c + e_ij

- η_m: pair noise floor, fixed from §4.4b (omitting it biases c below 1). a_h profiled per HOG; c common; weights 1/Var.
- Pairs with D_AA ∈ [0.05, 2]; HOGs with max/min D_AA ≥ 4.
- Pair-level SIMEX: bivariate noise λ·Σ_ij, λ ∈ {0.5, 1, 1.5, 2}, B = 50, quadratic extrapolant, frozen. §6.3 calibration includes tip noise and correlated errors.

**M4 (must). E4 specification.**

Y_hg = log(Σ_{e∈h,g} b̂_3Di / Σ_{e∈h,g} b̂_AA), Y_hg = μ + (β + β_h)·x_g + u_h + ε_hg

- g: 5 pre-fixed log-age bins by branch midpoint age √(t_parent·t_child); x_g centred log age; u_h, β_h random intercept and slope; weights 1/v_hg (delta method).
- Cells with Σb̂_AA < 0.05 merge with a neighbour (fixed rule, on x only). Fresh TimeTree age draw per bootstrap replicate.

**M5 (must). E5 specification.** Aggregate across HOGs for each species-tree branch e (lineage = child clade):

R_e = Σ_h b̂_3Di,he / Σ_h b̂_AA,he

- Residualise log R_e on the E4 fit and mean pLDDT; σ²_L = Var_e(log R_e) − mean(v_e), v_e = HOG-bootstrap variance.
- Sister contrasts on per-lineage AA and 3Di sums: x = log(ΣAA_A/ΣAA_B), y = log(Σ3Di_A/Σ3Di_B); disattenuated slope β̂ = (Ĉov(x,y) − c̄_err)/(V̂ar(x) − v̄_err,x).
- **N0 margins:** slope 90% CI ⊂ [0.8, 1.25] *and* σ_L upper 95% bound < 0.10 (log units).

**M6 (must). Bootstrap.** Species resampled with replacement give duplicate tips that cannot be fit. Use a pigeonhole bootstrap (Owen 2007) with multiplicity weights and no refit: HOG × species for pair estimands (weight m_i·m_j, no self-pairs), HOG × lineage for branch estimands. Mildly conservative; percentile CIs.

**M7 (must). Multiplicity.**
- E3 is estimation, outside Holm (only its CATH contrast is tested).
- One null and direction per endpoint: E1/E2 TOST plus one-sided c < 1; E4 two-sided; E5 equivalence; E6 one-sided.
- CI rules use Bonferroni-level CIs, 1 − 2α/m. E5's N0 claim is gated on E6 rejecting. E1/E2 otherwise "inconclusive".

**M8 (must). Missing and saturated pairs in the tree fits.**
- Choose each HOG's subset on x only: the largest clade-closed species set with all D_AA ≤ 2.
- Structural failures within it get weight 0 (sensitivity: impute at the floor). Exclude HOG×metric with >10% failures; report exclusion by root age.

**Should.** One ratio direction (log k everywhere); additivity reported, not gated; S7 null simulation; S9; S10.
