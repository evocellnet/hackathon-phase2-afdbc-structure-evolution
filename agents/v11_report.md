# V11: independent re-derivation of E1–E6 point estimates (confirm run)

Verifier: an independent agent that did not write the project code. Date: 2026-10-31. Outcome: **PASS**. All 42 checks pass (`results/verify/V11.json`, `pass_all: true`).

**Method.** I wrote my own implementation from ANALYSIS_PLAN §2.5/§2.6/§4.1–4.6 and DEVIATIONS (D22 primary = tdiaa, D8, D31, D32, C24, C14/C37–C41). No project endpoint or driver code was called. From the project I used only:
- the data loaders `load_pairs` and `load_branches`;
- the seed helper `stage_rng`, used once, to reproduce the plan's seeded E6 P1 choice.

I also used these inputs: the confirm pair tables with κ(bin)·SE applied as in §2.3, the confirm branch tables, `hogs_endpoint_cols_confirm`, `e6_triplets.tsv` and `config/tuning.yaml`. I read the project code only to resolve gaps in the plan text: the a_h used by the E6 reference curve, how E5 picks the sister-contrast node and sorts kids by label, and the E5 z_e definition.

**Tolerances.**
- Deterministic quantities: 1e-6.
- MixedLM: 1e-4. In practice the difference was 2.3e-7.
- Monte-Carlo quantities (SIMEX, D31 PB bias, E5 c_v/vx_v): 3× the combined MC SE. I used my own seeds; I did not replicate the project's RNG stream.

Scripts: `scripts/verify/v11_common.py`, `v11_inputs.py`, `v11_curve.py`, `v11_curve_pb.py`, `v11_branch.py`, `v11_e5.py`, `v11_e6.py`, `v11_assemble.py`. Logs are in `logs/verify/V11_*.log`.

## Results

| endpoint | quantity | project | own | abs diff | tol | kind | pass | notes |
|---|---|---|---|---|---|---|---|---|
| E1 | n pairs in fit | 591410 | 591410 | 0.00e+00 | 0 | exact | PASS |  |
| E1 | n HOGs in fit | 4850 | 4850 | 0.00e+00 | 0 | exact | PASS |  |
| E1 | s2 (pass-1 excess var) | 0.00173713 | 0.00173713 | 3.14e-12 | 1e-06 | deterministic | PASS |  |
| E1 | c_naive | 0.669541 | 0.669541 | 1.02e-08 | 1e-06 | deterministic | PASS | Brent xatol; diff ~1e-8 |
| E1 | c_SIMEX (plan) | 0.716581 | 0.716323 | 2.58e-04 | 0.0015 | Monte Carlo (own seed) | PASS | independent SIMEX draws, B_s=50; tolerance 3 x approx MC SE |
| E1 | SIMEX extrapolation from project c(lambda) | 0.716581 | 0.716581 | 0.00e+00 | 1e-06 | deterministic | PASS |  |
| E1 | D31 PB bias | -0.00604073 | -0.0057879 | 2.53e-04 | 0.0011390115702158406 | Monte Carlo (own seed) | PASS | K=20 own vs K=25; z=0.67 |
| E1 | c (final = SIMEX - PB bias) | 0.722621 | 0.722111 | 5.11e-04 | 0.0018834402982535856 | Monte Carlo | PASS | own SIMEX + own PB |
| E1 | final = c_simex_plan - pb_bias (project parts) | 0.722621 | 0.722621 | 0.00e+00 | 1e-06 | arithmetic | PASS |  |
| E1 | decision from CI rule | concave | concave |  |  | rule | PASS | ci98 recomputed [0.6420, 0.8033]; p_conc 6.077e-16 |
| E2 | n pairs in fit | 594618 | 594618 | 0.00e+00 | 0 | exact | PASS |  |
| E2 | n HOGs in fit | 4891 | 4891 | 0.00e+00 | 0 | exact | PASS |  |
| E2 | s2 (pass-1 excess var) | 0.00558984 | 0.00558984 | 9.00e-14 | 1e-06 | deterministic | PASS |  |
| E2 | c_naive | 0.627349 | 0.627348 | 2.66e-08 | 1e-06 | deterministic | PASS | Brent xatol; diff ~1e-8 |
| E2 | c_SIMEX (plan) | 0.689242 | 0.689487 | 2.45e-04 | 0.0015 | Monte Carlo (own seed) | PASS | independent SIMEX draws, B_s=50; tolerance 3 x approx MC SE |
| E2 | SIMEX extrapolation from project c(lambda) | 0.689242 | 0.689242 | 0.00e+00 | 1e-06 | deterministic | PASS |  |
| E2 | D31 PB bias | -0.00958697 | -0.0101322 | 5.45e-04 | 0.0015465662304479095 | Monte Carlo (own seed) | PASS | K=20 own vs K=25; z=1.06 |
| E2 | c (final = SIMEX - PB bias) | 0.698829 | 0.699619 | 7.90e-04 | 0.00215449927016972 | Monte Carlo | PASS | own SIMEX + own PB |
| E2 | final = c_simex_plan - pb_bias (project parts) | 0.698829 | 0.698829 | 0.00e+00 | 1e-06 | arithmetic | PASS |  |
| E2 | decision from CI rule | concave | concave |  |  | rule | PASS | ci98 recomputed [0.6060, 0.7917]; p_conc 2.284e-14 |
| E3 | k (ratio of sums) | 3.32181 | 3.32181 | 0.00e+00 | 1e-06 | deterministic | PASS |  |
| E3 | k internal edges only | 4.62797 | 4.62797 | 0.00e+00 | 1e-06 | deterministic | PASS |  |
| E3 | sigma2_k | -0.410292 | -0.410292 | 1.30e-13 | 1e-06 | deterministic | PASS | v_h with HOG-total (co)variances var_sum/cov_sum_aa; diag-only b_se gives -2.122 |
| E3 | CATH delta log k_alpha - log k_beta | 0.605671 | 0.605671 | 2.22e-16 | 1e-06 | deterministic | PASS |  |
| E3 | n HOGs | 4891 | 4891 | 0.00e+00 | 0 | exact | PASS |  |
| E3 | n edges | 146546 | 146546 | 0.00e+00 | 0 | exact | PASS |  |
| E4 | beta (MixedLM REML) | -0.150469 | -0.150469 | 2.31e-07 | 0.0001 | deterministic (optimiser) | PASS | REML optimum; statsmodels optimiser tolerance ~1e-5 |
| E4 | mu | -1.38244 | -1.38244 | 5.93e-09 | 0.0001 | deterministic (optimiser) | PASS |  |
| E4 | gamma | -0.0743816 | -0.0743817 | 1.22e-07 | 0.0001 | deterministic (optimiser) | PASS |  |
| E4 | n cells | 6597 | 6597 | 0.00e+00 | 0 | exact | PASS |  |
| E4 | censored fraction | 0.0268304 | 0.0268304 | 0.00e+00 | 1e-06 | deterministic | PASS |  |
| E4 | beta, plan cells (sensitivity) | -0.258957 | -0.258957 | 5.40e-08 | 0.0001 | deterministic (optimiser) | PASS |  |
| E5 | n lineages | 162 | 162 | 0.00e+00 | 0 | exact | PASS |  |
| E5 | sigma2_obs | 0.111601 | 0.111601 | 3.88e-08 | 1e-06 | deterministic | PASS |  |
| E5 | sigma2_L = sigma2_obs - F (project F) | 0.0978404 | 0.0978405 | 3.88e-08 | 1e-06 | deterministic given F | PASS | floor F is a P=100 parametric bootstrap; not re-simulated |
| E5 | n contrasts (finite x_v) | 146 | 146 | 0.00e+00 | 0 | exact | PASS | 147 contrasts pass the >=50-HOG rule; 146 have finite x_v (project's count); 145 have finite x_v AND y_v and enter beta |
| E5 | beta contrast | 0.572332 | 0.582982 | 1.06e-02 | 0.02274424168812566 | Monte Carlo (c_v, vx_v from 200 own HOG-pigeonhole reps) | PASS | 5 own seeds sd 0.0076; naive (no disattenuation) 0.5112 |
| E6 | E6a | 0.0296833 | 0.0296833 | 1.97e-12 | 1e-06 | deterministic (seeded P1 draw) | PASS |  |
| E6 | E6b | 0.0170244 | 0.0170244 | 6.44e-13 | 1e-06 | deterministic | PASS |  |
| E6 | n groups (a) | 8969 | 8969 | 0.00e+00 | 0 | exact | PASS |  |
| E6 | n groups (b) | 6027 | 6027 | 0.00e+00 | 0 | exact | PASS |  |


The E4 Wald p-value was recomputed and is identical (1.28e-12). The Holm adjustment was recomputed and matches holm.json exactly.

## Input sanity
- **Split.** Only `split == confirm` was read: 5,600 HOGs, the same set as `split.tsv`. Every pair, branch and triplet row is confirm, and no tune HOG appears anywhere.
- **Primary pairs.** There are 866,959 primary pairs, which equals Σ C(|S_h|, 2) over the confirm HOGs. This holds per HOG (0 mismatches) in all three alignments, and there are no duplicate pair_ids. No pair was dropped silently.
- **E1 exclusion chain.** It reproduces exactly: 164,002 not 1:1, 78,220 outside aa_ml_seq, 11,168 in excluded HOGs, 1,322 with weight 0, leaving 591,410 pairs from 4,850 HOGs (E2: 594,618 pairs from 4,891 HOGs).
- **fail_frac.** Recomputed from the pairs, it matches the project's HOG columns exactly.
- **E6 triplets.** The confirm set has 29,117 triplets. Together with the 7,383 tune triplets this gives 36,500, matching `logs/filters/hogs_e6.tsv`.
- **E1 and E2 decisions.** The 98% CIs are [0.642, 0.803] and [0.606, 0.792]. The V3 false-concave rate is 0, and the 99% upper bound is below 1. Rule 3 therefore gives "concave", and the recomputed p_conc matches.

## Observations (none changes a result)
1. **E1/E2 `n_excluded` double count.** The `hog_lt_min_pairs_or_range` bucket also counts weight-0 pairs that are already in another bucket: 54 for E1 and 55 for E2. As a result the buckets plus n_units add up to more than the number of primary pairs. n_units itself is correct.
2. **D31 wording does not match the code.** DEVIATIONS says the PB uses "edge-level, tree-additive noise". The code uses a tree-additive design x̂ = X·b̂_AA with independent pair-level Σ_ij noise. My PB follows what the code does and reproduces pb_bias (z = 0.67 for E1, 1.06 for E2). The documentation should be corrected.
3. **Two HOGs with NaN misfit_ratio are kept in E3–E5.** HOG0022006 and HOG0102266 stay in because NaN is not greater than τ. This follows the rule as written; the report should mention it.
4. **E5 `n_contrasts` overstates the β count by one.** The reported 146 counts contrasts with a finite x_v. One of them has a non-finite y_v, so β is based on 145.
5. **E5 values I did not re-simulate:**
   - The floor F was not re-simulated; I checked σ²_L as σ²_obs − F using the project's F.
   - β_contrast is Monte Carlo through c_v and vx_v. Over 5 own seeds I get 0.583 ± 0.008, against the project's 0.572.
6. **E6a depends on the seeded P1 draw.** Averaged over both orientations, E6a is 0.0308 instead of 0.0297; the conclusion is the same.
7. **Branch tables were taken as given.** I used them as inputs and did not refit them from the pairs. I did check the tip correction and that τ is the edge midpoint.
