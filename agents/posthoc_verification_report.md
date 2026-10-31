# Independent verification of the post-hoc statistical review (2026-10-31)

**Scope.** I re-tested the main claims of `reviews/results_review_statistics.md` with my own code. Nothing in `scripts/review/` was reused.
- Data access and the production estimators were used read-only, through `strucevo.stats`: `EndpointDriver.load`, `fit_branches`, `prepare_curve_data` and `run_E4`.
- Branch fits were refitted in memory. Nothing in `data/interim`, `results/endpoints` or the plan was written.
- Everything here is **post hoc and exploratory**. The pre-registered decisions are unaffected.

**Files.**
- Scripts: `scripts/posthoc/{branch_checks.py, pairwise_checks.py, pairwise_checks2.py, sigma_k_internal.py, common.py}`.
- Outputs: `results/posthoc/*.json`.
- Logs: `logs/posthoc/`.

**Sanity checks.** The production numbers reproduce exactly:
- k_all = 3.3218;
- k_internal = 4.6280;
- E4 β = −0.15047, with 6,597 cells;
- naive fixed-η ĉ = 0.6696, with the same Q as the reviewer.

## Summary table

| # | Claim | Reviewer value | My value | Verdict | Notes |
|---|---|---|---|---|---|
| 1a | Non-identical orthologs carry a 3Di excess above the identical-sequence floor (η = 0.0524) | ≈0.014–0.016 per protein | Branch level: terminal-edge excess vs k_int is 0.0143, 0.0167, 0.0133 and 0.0083 per edge (age bins <30, 30–100, 100–300, 300–1000 My). Pairwise, median-binned: +0.009 per pair already at d_AA ≈ 0.003, rising to a plateau of +0.019–0.022 per pair (≈0.010–0.011 per protein) from d_AA ≈ 0.1 | **partly confirmed** | The excess is real and of that order. Pairwise, it is ~0.010 per protein (smaller than claimed); at branch level it is 0.013–0.017. It is **not a pure step**: about half appears at the first substitutions, and the rest builds up over d_AA ≈ 0–0.1. Identical sequences between species inside HOGs sit near the floor (median 0.0551 vs 0.0524). |
| 1b | The excess is a per-protein artefact (pLDDT, disorder, length) | not tested | Per-protein effects u_i from y = u_i + u_j + a_h x^c: Spearman with protein pLDDT is −0.18 overall and −0.12 within HOG. Length: +0.04 overall, −0.07 within HOG. u_i median by pLDDT quintile: 0.034, 0.032, 0.029, 0.026, 0.020. Identical-pair floor/2 by the same quintiles: 0.031, 0.030, 0.028, 0.024, 0.020 | **supported (new)** | The protein effects track pLDDT as the AF floor does. This is consistent with an AF-model-quality component. At ĉ = 0.73 the pLDDT-matched excess is only ≈0.001–0.003, so pairwise data alone cannot separate "per-protein excess" from "short-range concavity". |
| 1c | ν that equalises terminal and internal k | 0.040 | k_term = k_int (4.63) at ν ≈ 0.0397 (k_term 4.46 at ν = 0.039, 4.64 at ν = 0.040) | **confirmed** | |
| 2a | HOG-specific baselines η_h raise naive ĉ | 0.903 / 0.941 / 0.982 (windows 0.05–2 / 0.05–1 / 0.3–2) | 0.903 / 0.941 / 0.981. Fixed η gives 0.670 / 0.614 / 0.727. Q(c=1) − Q(ĉ) under η_h = 496 / 67 / 3. Same with the tree-smoothed x̂ (D29 x): 0.903 | **confirmed (reproduced)** | Median η_h at c = 1 is 0.085 (q10–q90 0.042–0.142), matching the reviewer. |
| 2b | The more principled model y = u_i + u_j + a_h x^c | (not done) | ĉ = **0.728** (0.05–2), 0.739 (0.05–1), 0.571 (0.3–2). On x̂: 0.749. Q falls from 304,761 (η_h) to 82,625. Q(c=1) − Q(ĉ) = 1,465 | **does not support ĉ ≈ 0.9–1** | Per-protein effects fit far better than η_h, and they put ĉ back at ≈0.73. The fit is unstable across windows (0.57 at 0.3–2), and it is the estimator most sensitive to x errors (row 2c). |
| 2c | Synthetic check on the real design (591k pairs): truth c = 0.72 or 1.0; per-protein offsets Gamma, mean 0.015; realistic tree-additive x errors (from production branch covariances), ± 17% iid remainder | "η_h recovers 0.720 in a noise-free world"; "fixed η is biased by offsets" | **Noise-free x:**<br>– fixed: 0.72→0.57, 1.0→0.78 with offsets; unbiased without;<br>– η_h: 0.72→0.72, 1.0→1.00;<br>– prot: 0.72→0.72, 1.0→1.00.<br>**Tree-additive x noise, offsets present:**<br>– fixed: 0.56 / 0.77;<br>– η_h: 0.71 / 0.975;<br>– prot: 0.69 / 0.94.<br>**Tree + iid noise, offsets present:**<br>– fixed: 0.54 / 0.745;<br>– η_h: 0.665 / 0.915;<br>– prot: 0.56 / 0.75.<br>**iid-only x noise:** all three collapse (prot 0.36 / 0.47) | **confirmed, with a caveat** | The fixed-η estimator is strongly biased towards concavity whenever per-protein offsets exist (Δĉ ≈ −0.15 to −0.22). η_h is nearly unbiased under both truths with realistic noise (−0.01 to −0.09). The c = 1 + offsets world reproduces the observed naive triplet on raw x (0.67 / 0.90 / 0.73) far better than any c = 0.72 world, which predicts η_h ≈ 0.67 against the observed 0.90. **Caveat:** on x̂, no single world matches all three estimators. The implied c is ≈0.86 (fixed), ≈0.94 (η_h) and ≈0.77 (prot), so the pairwise data bracket c at roughly 0.75–1.0 rather than pinning it at 0.9–1. |
| 3a | E4 internal-edge k is flat from 30 My | 4.47, 4.63, 4.94, 4.68 (30–100, 100–300, 300–1000, >1000 My); 3.13 at <30 My | 4.47, 4.63, 4.94, 4.68; 3.13 at <30 My. Terminal k: 0.726, 2.61, 3.57, 4.16 | **confirmed** | The young internal bin (<30 My) is also low (3.13), with an excess of only 0.0015 per edge. |
| 3b | E4 β depends on ν, with a sign change between 0.026 and ~0.039 | −0.35, −0.29, −0.15, +0.07, +0.08, −0.08 (ν = 0, 0.013, 0.026, 0.039, 0.040, 0.052) | ν = 0: −0.352. 0.013: −0.291. 0.021: −0.213. **0.026: −0.150**. 0.031: −0.080. 0.035: −0.014. 0.039: +0.063. 0.040: +0.082. 0.045: +0.168. 0.052: −0.050 (20% censored cells) | **confirmed** | The sign changes at ν ≈ 0.036. Above ν ≈ 0.045, censoring (12–20%) makes β unreliable. |
| 3c | E4 internal-edges-only | −0.047 [−0.128, 0.034] | −0.0468 (reproduced with run_E4; at B = 20 the 98% CI is [−0.125, 0.032]; E4.json gives [−0.128, 0.034]) | **confirmed** | Internal edges carry no tip correction, so this value does not depend on ν. |
| C35 | Pre-specified ν ± 0.005 sensitivity (missing from E4.json) | missing (point grid only) | ν = 0.0212: β = −0.213 [−0.264, −0.162]. ν = 0.0312: β = −0.080 [−0.133, −0.027], p = 4e-4 (B = 200, production run_E4) | **was missing; now computed** | **The sign holds within the pre-specified ±0.005 range.** But β moves by −0.063 / +0.070 per ∓0.005, about twice the +0.03 / −0.04 that C35 anticipated. β reaches zero at +0.010 (ν ≈ 0.036). |
| 4a | E3 internal-edge k | 4.63 [4.35, 4.93] | 4.628. My HOG-only bootstrap 98% interval is [4.56, 4.69] (narrower, because it has no panel-edge weights). The production E3.json value [4.35, 4.93] is the right one to cite | **confirmed** | |
| 4b | The σ_k estimator fails; the robust value is 0.3–0.5 | σ²_plan = −0.41; DL 0.32; precise HOGs 0.47; internal ≈0.2 | Plan-style −0.410 (v_h mean 0.82, median 0.0039, q99 0.36). DL τ = 0.323. Median-v moment 0.63. MAD-based 0.44. v_h < 0.02 (n = 4,214): 0.471, with k_h q10–q90 1.69–5.19. **Internal-only** (with proper internal-only v_h, my addition): DL 0.26; precise subsets 0.23–0.28 (n = 666–2,770) | **confirmed** | The internal-only σ_k is slightly higher than the reviewer's "≈0.2". The all-edge σ_k also contains variation in the tip-artefact share across HOGs, so **~0.25 (internal) is the better biological estimate**. |
| 4c | The E3 stratum trend is the tip artefact | internal 3.74, 4.71, 4.63, 4.65; terminal 1.54→3.53 | internal 3.74, 4.71, 4.63, 4.65; terminal 1.54, 1.82, 2.87, 3.53 | **confirmed** | |

## Interpretation

1. **E4 is not robust to the tip correction.** This is confirmed, and it is the strongest finding.
   - β crosses zero at ν ≈ 0.036, which lies inside the plausible ν range: the mean identical-pair floor/2 is 0.031, and the value that equalises terminal and internal k is 0.040.
   - The internal-only β is null.
   - The C35 sensitivity the plan pre-specified was missing. At ν ± 0.005, β stays significantly negative: −0.213 and −0.080, both with 98% CIs below 0.
   - So E4 passes its *pre-specified* tip check. It fails the designated key check (internal edges only) and any ν above ≈0.036.
2. **A per-protein excess exists, but it is smaller and less step-like than claimed.**
   - It is about 0.010 per protein from pairwise medians and 0.013–0.017 per protein from branch terminal edges.
   - It tracks low pLDDT in the same way as the identical-sequence floor.
   - Within-HOG identical-sequence pairs are close to the floor. About half the excess appears at the first substitution.
3. **E1 concavity: "fixed η is biased towards concavity" is confirmed, but "c ≈ 0.9–1" is only partly supported.**
   - With per-protein offsets, the fixed-η estimator is biased towards concavity, and the η_h estimator is nearly unbiased under both truths.
   - The more principled per-protein-effect model returns ĉ ≈ 0.73–0.75. In realistic synthetic worlds it is itself biased low (c = 1 → 0.75–0.94).
   - Taken together, the pairwise data are compatible with c ≈ 0.75–1.0. They do not establish linearity.
   - The branch-level evidence (internal k flat from 30 My to >1 Gy) favours near-linearity more directly than any pairwise fit.
   - **Recommended wording:** "pre-registered ĉ = 0.72 is not robust to the baseline specification; post hoc estimators give 0.73–0.98."
4. **E3:** k_internal = 4.63 [4.35, 4.93] and the failure of the plan's σ_k estimator are both confirmed. My internal-edge σ_k is ≈0.25 (DL 0.26).

## Appendix: E4 reproduction and C35 with bootstrap
Source: `results/posthoc/branch_checks.json` (keys `nu_grid`, `E4_prod_nu_reproduced` and `C35`).

**E4 reproduction at production ν (0.0262):**
- β = −0.15047, identical to E4.json;
- internal-only β = −0.0468.

**C35 (B = 200):**

| ν | β | 98% CI | p | SE |
|---|---|---|---|---|
| 0.0212 | −0.213 | [−0.264, −0.162] | <1e-6 | 0.022 |
| 0.0312 | −0.080 | [−0.133, −0.027] | 4.3e-4 | 0.023 |

**Limitations.**
- The synthetic worlds use iid y errors with independent per-protein offsets. They leave out AA/3Di error correlation (ρ ≈ 0.1–0.3).
- The pairwise bins condition on noisy observed d_AA.
- None of the naive pairwise estimators is SIMEX-corrected. Their absolute values are biased; only the comparisons between them carry information.
