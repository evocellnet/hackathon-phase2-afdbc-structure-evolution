# ANALYSIS_PLAN.md: executable analysis specification

> **Status: FROZEN for the confirmation run once the tuning parameters are filled in §Tuning.** Until then only §Tuning may change. After the freeze, every change is logged in `DEVIATIONS.md` (date, reason, effect) and the affected result is labelled "deviation".
> (Written by the round-3 consistency reviewer; saved verbatim by the PI because the reviewer's harness blocked file writes. 2026-10-30.)

Sources, in order of precedence:
1. `reviews/round3_consistency.md` (R3-M\*, R3-S\*).
2. DESIGN.md v3 §11.
3. Non-superseded DESIGN §0–§10.

Where this plan is more specific than these sources, this plan governs.

## 0. Global conventions

| Item | Definition |
|---|---|
| Units | ML distances in expected substitutions/site; ages in Mya; LDDT and TM in [0, 1] |
| Metrics m | `aa_ml` (LG+Γ4); `tdi_ml` (Q.3Di.AF+Γ4); `lddt_D` = −ln((S − S0(L))/(1 − S0(L))), with S = global LDDT |
| Ratio direction | Internally r = Σb_3Di/Σb_AA. Report k = 1/r (Illergård's ratio) |
| Primary alignment | `tmalign` (Foldseek). Sensitivity alignments: `seqonly` (Biopython PairwiseAligner: global, BLOSUM62, gap open −11, extend −1, end gaps free) and `tdiaa` (Foldseek 3Di+AA) |
| Columns for ML distances | Core-masked: aligned columns with pLDDT ≥ 70 in both proteins |
| Family-wise error | α = 0.05 over F = {E1, E2, E4, E5, E6}; m = 5; α_e = 0.01. CI-based decisions use two-sided **98%** CIs and one-sided **99%** bounds. Holm is applied to the five endpoint p-values (§5) |
| CIs | Normal approximation θ̂ ± z·SE_boot; percentile CIs also reported. If the two disagree on a decision, flag `ci_disagree` |
| Bootstrap | Pigeonhole multiplicity weights, no tree refit (§2.4). B = 1000 unless stated |
| Seeds | `np.random.SeedSequence(20261030)`, spawned per stage name; the map is stored in `config/seeds.yaml` |
| Age bins | [0,30), [30,100), [100,300), [300,1000), [1000,∞) Mya, assigned on the arithmetic branch midpoint τ = (t_parent + t_child)/2. Bin 5 is also reported separately |
| Root-age strata (HOG sampling) | <100, 100–300, 300–1000, ≥1000 Mya |
| Licence | All tables containing ages live in `data/interim/timetree/` (gitignored). `results/` holds aggregates only |

## 1. Pipeline data contract

All tables are parquet. Required columns are non-null when `status == "ok"`. `scripts/verify/check_contract.py` validates names, dtypes, ranges and row counts.

### 1.1 `data/interim/hogs.parquet`: one row per candidate HOG

| Column | Type | Meaning |
|---|---|---|
| hog | str | Ancestra top-level HOG id |
| n_species_1to1 | int | Species with exactly one filtered copy |
| presample | list[int] | ≤ 60 species taxids (§1.6 step 2) |
| selected | list[int] | Final S_h (8–30 species), or empty |
| sel_clade | str | Panel node that defines S_h |
| root_age_mya | float | *Interim/timetree only.* Age of MRCA(S_h) |
| stratum | int | 1–4 |
| in_sample | bool | Among the ≤ 7,000 sampled HOGs |
| split | str | `tune` / `confirm` / `none` |
| cath_class | str | `alpha` / `beta` / `alphabeta` / `mixed` / `none` (§4.3) |
| mean_plddt | float | Over the selected proteins |
| misfit_ratio | float | RSS_rel(species tree) / RSS_rel(NJ tree); AA, tmalign |
| excl_misfit | bool | misfit_ratio > τ_misfit |
| fail_frac_{m}_{aln} | float | Failed-pair fraction for each metric × alignment |
| excl_{m}_{aln} | bool | fail_frac > 0.10 |

### 1.2 `data/interim/pairs/{aln}/*.parquet`: one row per pair × alignment (sharded, with `.done` sentinels)

| Column | Type | Meaning |
|---|---|---|
| pair_id | str | `{gene_a}__{gene_b}` |
| pair_set | str | `primary`, `e6_para`, `e6_orth`, `e6_ldo`, `floor_random`, `floor_ident`, `h1e_af`, `h1e_pdb`, `pclu` |
| hog, split | str | HOG id (`NA` for floor pairs); split as in §1.1 |
| gene_a, gene_b | str | **gene_a = the shorter protein** (ties broken lexicographically) = the Foldseek query |
| species_a, species_b | int | Taxids |
| relation, one_to_one, dcg_id, event_level, dup_branch_from | — | From `hog_pair_labels.py`; `relation` is `random` or `identical` for floor pairs |
| len_a, len_b, mean_plddt_a, mean_plddt_b | — | |
| n_aln, n_core, cov_short | int, int, float | Aligned columns; aligned columns with pLDDT ≥ 70 in both; n_aln/len_a |
| status | str | `ok`, `fail_align`, `fail_core` (n_core < 50), `fail_cov` (cov_short < 0.5) |
| aa_p, aa_ml, aa_ml_se, aa_sat | float/bool | Core columns. `d_max = 20`. SE floor = 1/n_core. `sat` = the engine's `saturated` flag |
| tdi_p, tdi_ml, tdi_ml_se, tdi_sat | float/bool | Same as above, for 3Di |
| aa_ml_seq, aa_pid_seq | float | AA ML distance and identity from the `seqonly` alignment. These are the selection variables and are present in every alignment table |
| lddt_global | float | Σ per-residue LDDT of aligned query residues / len_a |
| lddt_core | float | Mean per-residue LDDT over aligned core residues |
| lddt_D, lddt_D_se, lddt_sat | float/bool | `sat` if S ≤ q95 of S0 in the pair's length bin |
| tm_short | float | TM-score normalised by len_a (`qtmscore`) |
| rmsd_core | float | `kabsch_rmsd` over aligned core residues |
| aa_ml_hA, aa_ml_hB, tdi_ml_hA, tdi_ml_hB (+ `_se`) | float | Split-half distances; `e6_*` pairs only |
| boot_sub | bool | In the column-bootstrap subsample |

Raw counts go in `counts_*.npz`: 20×20 count matrices for AA and 3Di (core columns, plus halves for E6 pairs). ML distances are recomputed from these after α tuning, without re-aligning.

### 1.3 `data/interim/branches/{aln}.parquet`: one row per HOG × metric × HOG-tree edge

| Column | Meaning |
|---|---|
| hog, split, metric | |
| edge_id | Sorted taxids of the child clade |
| is_terminal | |
| panel_child, panel_parent | Panel MRCA of the child clade and of the parent clade |
| panel_edge | The panel edge id if the HOG edge equals exactly one panel edge; otherwise NA |
| n_panel_edges | |
| t_child, t_parent, tau | *Interim/timetree only* |
| b_raw | WLS estimate; may be negative |
| b | b_raw − ν_m on terminal edges, b_raw otherwise |
| b_se | Sandwich SE |
| cov_aa | Covariance with the AA estimate of the same edge |
| b_nnls | Display only |

A companion table `{aln}_fit.parquet` has one row per HOG × metric: n_pairs_ok, n_edges, s2_hm, rss_rel, excluded, reason.

### 1.4 `results/floors/floors.json`

- η_m: median distance over `floor_ident` pairs. ν_m = η_m/2; η_aa = 0.
- S0 table: median and q95 by length bin, for LDDT_global and TM.
- ρ_m(bin): correlation between AA and 3Di estimation errors.
- κ_m(bin): bootstrap SD / Fisher SE.
- SE model for lddt_D.

### 1.5 Panel tree and ages

- `data/interim/timetree/panel_tree.nwk`: collapsed per §2.2.
- `data/interim/timetree/node_ages.tsv`: node, median, lo95, hi95, source.
- `results/panel/species.tsv`: taxid, mean_plddt, n_proteins.

### 1.6 Construction order

Every step logs counts in → out to `logs/filters/*.tsv`.

1. **Gene-model filters.**
   - md5(UniProt sequence) = md5(AFDB v4 sequence).
   - Not a fragment: UniProt fragment flag, or AFDB cluFlag = 3.
   - Mean pLDDT ≥ 70.
   - Length within 0.5–2× the HOG median.
   - Form the 1:1 set: species with exactly one copy. Keep HOGs with ≥ 8 such species.
2. **Presample.** Choose ≤ 60 species by greedy panel-node coverage (seeded tie-break). Run the `seqonly` alignment on all presample pairs to get `aa_ml_seq`.
3. **Clade-closed selection.**
   - Visit panel internal nodes C ordered by: presample species in C (descending), then older first, then node id.
   - S = presample ∩ C. Accept the first C with |S| ≥ 8 and every pair `aa_ml_seq` ≤ 2. Saturated or failed pairs count as > 2.
   - If |S| > 30, reduce to 30 by greedy node coverage.
   - If no C qualifies, drop the HOG.
4. **Sample ≤ 7,000 HOGs.** Equal quotas per stratum; unused quota is redistributed proportionally (seeded). If primary + E6 pairs exceed 3.5M, scale all quotas down.
5. **Split 20/80** by seeded permutation within strata. Write `data/interim/split.tsv` and record its md5 in the ledger. Never regenerate it.
6. **Primary pairs**: all pairs within each S_h, for each of the three alignments.
7. **E6 pairs** (§4.6): drawn from the full filtered membership of the sampled HOGs.

## 2. Shared estimators

### 2.1 Distances and LDDT

- ML distances: `ml_distance(counts, model, alpha=α̂_m, ncat=4, d_max=20)`.
- Load Q.3Di.AF (`data/raw/models/Q.3Di.AF`, md5 8f897fff…) with `load_paml(path, alphabet=AA_ORDER)`. The file is in PAML ARND order, so encode 3Di strings with `AA_ORDER`.
- Global LDDT: query = the shorter protein; LDDT_global = Σ`lddtfull`/len_a. **Never use Foldseek's `lddt` field.**
- The in-house `strucevo.lddt` must match Foldseek within 0.01 (check V9). It is used for the `seqonly` and `tdiaa` alignments.

### 2.2 Species tree

- Collapse internodes shorter than 5 My, and internodes whose TimeTree CI overlaps the parent's, into polytomies.
- Use TimeTree API medians. Check that ages decrease from root to tips; log every fix.

### 2.3 SE calibration (once; tuning pairs plus floor sets)

- Take 20k primary pairs, stratified over 10 `aa_ml` bins.
- Run 100 column-bootstrap replicates per pair, resampling 20-residue blocks.
- Fit κ_m(bin) and ρ_m(bin). The SE used everywhere is κ × Fisher SE.
- SE of lddt_D: delta method applied to the bootstrap SD of S, modelled by length bin and S bin.

### 2.4 Pigeonhole bootstrap

Draw three sets of multinomial weights:
- W_h over HOGs;
- M_s over species;
- Λ_e over panel edges.

Apply them as follows:
- **Pair weight:** W_h·M_i·M_j.
- **Edge weight:** W_h × the mean Λ over the panel edges the HOG edge spans.
- **No refit.**
- **Mixed models:** replicate each HOG W_h times under new group ids.

### 2.5 Branch engine (E3–E5)

For each HOG h, metric m and alignment, work on T_h: the collapsed panel tree pruned to S_h, with unary nodes suppressed.

1. **Design matrix.** X is the path incidence matrix (pairs × edges). Use only `ok`, unsaturated pairs; failed pairs get weight 0. If n_ok < n_edges + 5, mark the fit "underdetermined" and exclude it.
2. **Weights in two passes.** First pass w = 1/SE². Then set s²_hm = max(0, RSS_w/(n_ok − n_edges)·mean(SE²) − mean(SE²)). Second pass w = 1/(SE² + s²_hm). Estimate b̂ = (XᵀWX)⁻¹XᵀWd. **No clipping** of negative values.
3. **Sandwich.** HC1: V = A diag(r²) Aᵀ · n/(n − p), with A = (XᵀWX)⁻¹XᵀW. The cross-metric covariance is A_AA diag(r_AA r_m) A_mᵀ; store its diagonal as `cov_aa`.
4. **Tip correction.** Subtract ν_m from terminal edges.
5. **Fit quality.** RSS_rel = Σw r²/Σw d².
6. **Misfit.** Fit the same AA WLS on an NJ topology (dendropy `nj_tree`; FastME is not installed). Exclude the HOG from E3–E5 if misfit_ratio > τ_misfit.
7. **Hidden paralogy.** Flag a pair when its AA residual / (1.4826·MAD) > 3. Used in sensitivity analyses only.

### 2.6 Cells (E4)

- Cell (h, g) = the edges of HOG h whose τ falls in bin g. Compute S_AA = Σb_AA and S_3Di = Σb_3Di.
- **Merge rule (based on AA only).** While S_AA < 0.05, merge the cell into the next-older bin; the oldest cell merges into the next-younger.
- If S_3Di ≤ 0 after merging, set it to 0.005 and mark the cell `censored`. If more than 2% of cells are censored, flag E4.

## 3. Tuning (20%) vs confirmation (80%) protocol

1. **Prerequisites.**
   - `split.tsv` exists.
   - Checks V1–V3, V5 and V10 pass.
   - Pair metrics are computed for all sampled HOGs.
   - No endpoint code has been run on `confirm`.
2. **Tune on `split == tune` only.**
   - α_AA and α_3Di: maximise the pooled likelihood with `profile_alpha` over tmalign primary core counts with `aa_ml_seq` ≤ 2; refine to ±0.05.
   - S0(L): ≈100k `floor_random` pairs, ≥ 2,000 per length bin.
   - τ_misfit = max(3, the 95th percentile of misfit_ratio under parametric simulation from the species-tree fit plus noise).
   - Measure η, ν, κ, ρ and the tmalign rate.
   - Set B_E4 to the largest of {1000, 500, 200} that keeps the E4 bootstrap ≤ 1 h.
3. **Also on tune:**
   - E5 spike-in power (§4.5);
   - the V3 curvature simulation;
   - the V4 end-to-end synthetic test;
   - dry runs of all endpoints.

   Dry-run outputs go to `results/tuning/` as exploratory and are never pooled with confirmation results.
4. **Freeze.**
   - Fill §Tuning and `config/tuning.yaml`.
   - Commit.
   - Record in `RESULTS_LEDGER.md` the commit hash and the md5 of `ANALYSIS_PLAN.md`, `tuning.yaml` and `split.tsv`.
5. **Confirm.**
   - `scripts/endpoints/run_confirm.py` refuses to run unless `tuning.yaml` has no nulls, the git tree is clean, the plan md5 matches the ledger, and `CONFIRM_RUN=1` is set.
   - `load_pairs(split="confirm")` raises unless `CONFIRM_RUN=1`.
   - Run once. Any rerun goes in `DEVIATIONS.md`.

## §Tuning (fill on the 20% set; then FROZEN)

| Parameter | Value | Source/date |
|---|---|---|
| α_AA (LG+Γ4) | `<TBD>` | |
| α_3Di (Q.3Di.AF+Γ4) | `<TBD>` | |
| S0(L) floor table (LDDT_global, TM; bins, median, q95) | `<TBD: floors.json#S0>` | |
| HOG-misfit threshold τ_misfit | `<TBD>` (≥ 3) | |
| η_tdi, η_lddtD (ν = η/2) | `<TBD>` | measured |
| κ_m(bin), ρ_m(bin) | `<TBD>` | measured |
| B_E4 | `<TBD>` | timing |
| HOGs per stratum; total primary + E6 pairs | `<TBD>` | |
| tmalign pairs/s | `<TBD>` | |
| E5 spike-in power at f = 1.5 (also 1.25, 2; null rate) | `<TBD>` | §4.5 |
| V3 false-concave rate; minimal detectable \|c−1\| and \|β_E4\| | `<TBD>` | §7 |

## 4. Endpoints

All endpoints use `split == confirm` and `aln == tmalign` unless stated otherwise.

Each endpoint writes `results/endpoints/E*.json` with the keys `{endpoint, estimand, estimate, se_boot, ci98, ci_pct98, one_sided_bound, p_value, p_holm, decision, n_hogs, n_units, n_excluded{reason:n}, sensitivity{name:{estimate, ci98, decision}}, alignment_robust, flags, plan_md5, tuning_md5, commit, input_md5}`.

### 4.1 E1: curvature of 3Di-ML vs AA-ML (pairwise engine)

**Inputs.** Primary pairs that are 1:1 orthologs. y = `tdi_ml`, x = `aa_ml`.

**Filters.**
- `aa_ml_seq` ∈ [0.05, 2].
- HOG not excluded for this metric (`excl_tdi_ml_tmalign` is false).
- Pairs that are not `ok` or are saturated get weight 0 and are counted.
- HOG retained if it has ≥ 10 eligible pairs and max(x)/min(x) ≥ 4.

**Model.** y = η_tdi + a_h·x^c + e, with η fixed.
- Weights: w = 1/(SE_y² + s²), where s² = max(0, mean r² − mean SE_y²) from a first pass with w = 1/SE_y².
- a_h(c) = Σw(y − η)x^c / Σw x^{2c}.
- Q(c) = Σw(y − η − a_h(c)x^c)².
- ĉ_naive = argmin Q: grid over 0.2–3 in steps of 0.01, then bounded Brent.

**SIMEX.**
- λ ∈ {0.5, 1, 1.5, 2}, B_s = 50.
- Add noise N(0, λΣ_ij), where Σ_ij has variances κ²SE_x² and κ²SE_y² and covariance ρκ²SE_xSE_y.
- Clip x at 1e-4 and count the clips.
- Fit a quadratic in λ (including λ = 0) and extrapolate to λ = −1.

**Uncertainty.** 200 pigeonhole replicates, each running SIMEX with B_s = 10. SE_boot = SD across replicates. The naive ĉ is also reported with B = 1000.

**Decision** (first matching rule):
1. V3 false-concave rate > 10% → *inconclusive*.
2. 98% CI ⊂ [0.9, 1.1] → *linear*.
3. 99% upper bound < 1 → *concave*.
4. 99% lower bound > 1 → *inconclusive (convex)*.
5. Otherwise → *inconclusive*.

**p-value for Holm.**
- p_lin = max(1 − Φ((ĉ − 0.9)/SE), 1 − Φ((1.1 − ĉ)/SE)).
- p_conc = Φ((ĉ − 1)/SE).
- p_E1 = the p-value of the declared class, or min(p_lin, p_conc) if inconclusive.

**Sensitivities.**
- η free.
- d_AA ≤ 1.
- Failures imputed at the cap.
- Hidden-paralogy pairs removed.
- Single-copy HOGs only.
- All aligned columns.
- α ×1.5 and ÷1.5.
- `seqonly` and `tdiaa` alignments. The result is alignment-robust if the class and the sign of (ĉ − 1) agree across all three.

**Checks.**
- The n entering the fit equals the filter log.
- An independent recomputation matches within 1e-6.
- The grid minimum is not at a boundary of the grid.

**Output:** `E1.json`.

### 4.2 E2: floor-corrected global LDDT

Same as E1 with y = `lddt_D`, η = η_lddtD, SE = `lddt_D_se`, ρ_lddt, and saturation = `lddt_sat`.

**d\*.**
- y_raw = 1 − lddt_global.
- Fit a linear model η + a_h·x on x ≤ 0.5 and a power model η + a'_h·x^{c'} on x ≤ 2. Both are profiled.
- A and A' = weighted medians of the per-HOG a_h and a'_h.
- d* = the smallest x on the grid 0.05:0.01:2 with |(η + A'x^{c'})/(η + Ax) − 1| > 0.10.
- 98% CI from 200 replicates. Stored as key `d_star`.

**Output:** `E2.json`.

### 4.3 E3: conservation ratio k (estimation only; outside Holm)

**Inputs.** Branch table, `aa_ml` and `tdi_ml`. HOGs excluded for misfit or failures are dropped.

**Estimands.**
- k = ΣΣb_AA / ΣΣb_3Di over all edges, with terminal edges tip-corrected.
- Per HOG, k_h. Between-HOG spread σ²_k = Var(log k_h) − mean v_h, where v_h comes from the delta method using b_se and cov_aa. Report σ_k = √max(0, σ²_k).

**CATH contrast.** Secondary; unadjusted, two-sided, α = 0.05.
- HOG class = the TED/CATH class that covers ≥ 60% of TED-domain residues; otherwise the HOG is `mixed` and excluded from the contrast.
- Δ = log k_α − log k_β, each as a ratio of sums.

**Uncertainty.** B = 1000.

**Also report.**
- k by length tertile, pLDDT tertile and stratum.
- k from internal edges only.
- Whether the 98% CI overlaps [3, 10].

**Alignment-robust** if k agrees within 25% across alignments.

**Output:** `E3.json`.

### 4.4 E4: structure/AA ratio vs age

**Variables.**
- y_hg = log(S_3Di/S_AA).
- x_hg = log(mean τ of the cell), centred.
- z_hg = log(S_AA/n_edges), centred.

**Model.** statsmodels `MixedLM`, REML, **unweighted**: y = μ + βx + γz + u_h + β_h x + ε, with random effects (1 + x | hog). If it does not converge, fall back to a random intercept only and flag it.

**Interpretation.** Use **β_k = −β**.

**Uncertainty.** B_E4 replicates. Each replicate:
- cluster-resamples HOGs under new ids;
- applies the edge weights inside the cell sums;
- draws one fresh set of ages: lognormal with σ = (ln hi − ln lo)/3.92, then a top-down repair child = min(child, parent − 0.1);
- recomputes the cells.

**Decision.** Wald p = 2(1 − Φ(|β̂|/SE)). Significant if p_holm < 0.05.

**Sensitivities.**
- WLS with weights 1/v_hg and HOG-clustered SEs.
- Pooled per-bin ratio-of-sums slope.
- Censored cells dropped.
- Internal edges only.
- Strassert 2021 / Betts 2018 deep ages.
- α ×1.5 and ÷1.5.
- Bin 5 excluded.
- `aa_ml_seq` ≤ 1.
- Other alignments.

"Robust to α and timescale" means the same sign and unadjusted p < 0.05 in both the α and the deep-age sensitivities.

**Negative control.** On the tuning set, 20 permutations of node ages: the 98% CI must cover 0 in at least 18.

**Output:** `E4.json`.

### 4.5 E5: lineage effects and N0

**Lineages.** Panel edges e with ≥ 50 HOGs contributing a simple edge mapped to e, and D_e = Σ_h b_AA ≥ 0.05.

**Lineage ratios.**
- R_e = Σ_h b_3Di / D_e.
- ℓ_e = log R_e − (μ̂ + β̂x_e + γ̂z_e), using the E4 fixed effects.
- ρ_e = the residual of an OLS of ℓ_e on the mean pLDDT of the descendant species.

**Excess variance.**
- σ̂²_obs = Var(ρ_e).
- Floor F: parametric bootstrap with P = 100. In each replicate, simulate d*_AA = X b̂_AA + ε and d*_3Di = X(r_h · e^{β̂(x_e − x̄_h)} · b̂_AA) + ε, with bivariate noise (variances κ²SE² + s², correlation ρ). Refit and recompute σ̂²_obs; F = the mean over replicates.
- σ²_L = σ̂²_obs − F.
- CI: B = 1000 with F fixed. Report the 99% one-sided bounds and σ_L = √max(0, ·).

**Sister contrasts.**
- Take internal nodes v whose children A and B each have ≥ 50 HOGs with species on both sides. For a polytomy, contrast consecutive children ordered by id.
- L_m,h(A) = the mean path length from the HOG node mapped to v down to the species of S_h in A.
- x_v = log(Σ_h L_AA(A)/Σ_h L_AA(B)); y_v = the same for 3Di. A = the child with the smaller id.
- β̂ = (Σx_v y_v − Σc_v)/(Σx_v² − Σvx_v), through the origin, where c_v and vx_v come from 200 HOG-pigeonhole replicates.
- CI: B = 1000 with c_v and vx_v held fixed.

**Spike-in power (tuning set).** 50 seeded replicates:
- Pick 10% of lineages.
- For every pair whose path crosses a simple edge mapped to a chosen lineage, add (f − 1)·b̂_3Di,he to y_tdi, with f = 1.5.
- Refit and run E5.
- Detection = the lower 99% bound of σ²_L > 0, or the 98% CI of β excludes 1. Power = the detection rate.
- Also report f = 1.25, f = 2, f = 1 (the null rate) and P(accept | f = 1.5).
- An MDO-reinsertion realism check (≤ 100k extra pairs) is exploratory.

**Guards.** Internal edges only; hidden-paralogy pairs removed; single-copy HOGs.

**Decision.**
- **N0 accepted** if all of:
  - the 98% CI of β ⊂ [0.8, 1.25];
  - the σ_L 99% upper bound < 0.10;
  - E6 is significant after Holm;
  - power ≥ 0.80.
- **N0 rejected** if either:
  - the 98% CI of β is disjoint from [0.8, 1.25]; or
  - the σ_L 99% lower bound > 0.10;

  and the rejection holds under every guard. If any guard fails, call it *artefact-suspect* instead.
- **Inconclusive** otherwise.

**p_E5.** Intersection–union test: max(p_{β<0.8}, p_{β>1.25}, p_σ), using one-sided Wald tests. If the decision is a rejection, use the smallest rejecting p instead.

**Output:** `E5.json`, with `sigma_L`, `beta_contrast`, `power`, `floor_F`, `n_lineages`, `n_contrasts`.

### 4.6 E6: duplication positive control

**Triplets.**
- P1 and P2: `in_paralog` genes in one species, in a sampled HOG.
- Paralog level P ≠ leaf (excludes recent in-paralogs).
- `aa_pid_seq` ≤ 0.95 and `aa_ml_seq` ≤ 2.
- O1 = species in the sister clade of P whose single gene is labelled ortholog to both copies.
- O2 = the same, but from the sister of parent(P).
- Require |O1| ≥ 1, and |O2| ≥ 1 for part (b).
- The copy called P1 is chosen at random (seeded).

**Halves.** Core columns, assigned by the paralog's residue index to alternating 20-residue blocks: A = odd blocks, B = even blocks.

**Reference curve.** f_h(x) = η_tdi + â_h·x^ĉ from E1. HOGs without an E1 fit are excluded. Residual r = y − f_h(x).

**(a) Paralog vs ortholog.**
- Use x = `aa_ml_hA` and y = `tdi_ml_hB`.
- Δ = r(P1, P2) − mean_{O1} r(P1, O).
- E6a = the mean over paralog groups (`dcg_id`) of the within-group mean Δ.

**(b) MDO vs LDO.**
- MDO = the copy with the larger mean `aa_ml_hA` to O1.
- Excess = mean_{O2}[r_B(MDO, O) − r_B(LDO, O)], where r_B uses x and y both from half B.
- E6b = the mean over paralog groups.

**Uncertainty.** HOG-only pigeonhole bootstrap, B = 1000. One-sided Wald tests give p_a and p_b.

**Decision.** p_E6 = min(1, 2·min(p_a, p_b)). The control is detected if p_holm < 0.05.

**Secondary.** Repeat with lddt_D on full columns (flagged as not split-half valid). First report `replication_direction` (the direction of the PMID 40813247 result).

**Checks.**
- No paralog level P is a leaf.
- O1 ∩ O2 = ∅.
- Halves are disjoint and together cover the core.
- In the V4 N0 world, the E6 CI covers 0.

**Output:** `E6.json`.

## 5. Multiplicity and interpretation

Outputs: `results/endpoints/holm.json` and `results/endpoints/interpretation.json`.

**Multiplicity.**
- Apply Holm at 0.05 to {p_E1, p_E2, p_E4, p_E5, p_E6}.
- E4 and E6 are decided on p_holm.
- E1, E2 and E5 are decided on their 98% CI rules. These are Bonferroni-level, so they are valid under Holm.
- E3 is not in the Holm set; its CATH contrast is tested unadjusted.

**Interpretation.**
- `scripts/endpoints/interpret.py` matches the endpoint decisions to a row of DESIGN §11.3 and writes the row id only, with no free text. The E4 column of that table is read on β_k.
- The joint interpretation requires power ≥ 0.80. Otherwise E5 acceptance becomes inconclusive, while the E1–E4 decisions and any E5 rejection still stand.
- If E6 > 0 and N0 is accepted, the "duplication" row applies.

## 6. Exploratory analyses (unadjusted 95% CIs; `results/exploratory/*.json`)

- **X1.** TM, RMSD, p-distances, `lddt_core` and raw 1 − S: curves and naive ĉ.
- **X2.** Root-to-tip and path distances vs age per metric; time-curvature fits; bin 5 shown separately.
- **X3.** Tajima-style relative-rate tests for AA and 3Di; effect sizes only.
- **X4 (H1e).** SIFTS ortholog pairs where both sides have X-ray structures (≤ 3 Å, coverage ≥ 70%, released after 2018-04-30).
  - Compare PDB–PDB vs AF–AF distances: paired differences, and their slope vs d_AA.
  - Floor: the same protein's AF vs PDB distance.
- **X5 (2-ii).** HOG ↔ AFDB-v3 cluster contingency table; logistic model split ~ TED architecture change + log length ratio + pLDDT + d_AA.
- **X6 (2-iii).** Actin and tubulin on a fixed IQ-TREE LG+G4 topology (keep nodes with UFBoot ≥ 95); branch engine; Asgard stem; 4 h; descriptive only.
- **X7.** P-clu vs ortholog curves.
- **X8.** Additivity (RSS_rel), reported only.
- **X9.** MDO-reinsertion check.

## 7. Verification (`make verify` → `results/verify/*.json`)

| ID | Check | Pass |
|---|---|---|
| V1 | ML engine, with AliSim LG+Γ and Q.3Di.AF+Γ | Bias < 3% for d ≤ 2; within 2% of IQ-TREE `.mldist`; identical sequences → 0; random → saturated |
| V2 | Branch engine on simulated additive distances plus noise | r > 0.99; unbiased; 95% sandwich coverage in [0.92, 0.97] |
| V3 | 100 simulated linear datasets with the tuned noise (κ, ρ, tip noise, real HOG designs) | False-concave rate ≤ 10%, otherwise E1/E2 inconclusive. Report the minimal detectable \|c − 1\| at 80% power |
| V4 | End-to-end synthetic: b_3Di = r_h·b_AA per HOG tree, `simulate_counts`, full pipeline | E1 CI covers 1; E4 CI covers 0; E5 accepts N0; the spike-in world is detected |
| V5 | Invariants | Pairs out = pairs in; self-pairs give TM = LDDT = 1; identical sequences give d_AA = 0; md5 identity; human–mouse age 87 ± 10 Mya; human–mouse 1:1 identity 80–90% |
| V6 | Negative controls | Random pairs reproduce S0; shuffled ages (§4.4) give no slope |
| V7 | Data contract | `check_contract.py` passes; no age columns under `results/` |
| V8 | Parsers | Foldseek column order; cigar-derived alignment reproduces `qaln`/`taln`; OrthoXML toy labels |
| V9 | LDDT | In-house vs `lddtfull`: \|Δ\| ≤ 0.01 on 1000 pairs |
| V10 | Model files | Q.3Di.AF md5 8f897fff…; letter-order check |
| V11 | Independent re-derivation | A non-author agent recomputes each E\*.json point estimate (tolerance 1e-6) |
| V12 | Ledger | Every reported number has a `RESULTS_LEDGER.md` line |

## 8. Execution traps

1. Never clip negative branch lengths. Never drop NaN silently. Exclude pairs only through the logged status/weight-0 rule.
2. Foldseek's `lddt` field is not global LDDT. The query must be the shorter protein.
3. Pass `d_max = 20` (the default is 10). The SE is NaN at d = 0: apply the 1/n_core floor.
4. Encode 3Di with `AA_ORDER` for Q.3Di.AF, not `THREEDI_ORDER`.
5. Do not resample species with replacement in tree fits.
6. `MixedLM` has no observation weights.
7. Ages stay in `data/interim/timetree/`.
8. The confirm split can be read only with `CONFIRM_RUN=1`.
9. Use β_k = −β for the §11.3 table.
10. Use `-e inf --tmscore-threshold 0` and exhaustive alignment of the pair list, and assert rows out = rows in.
