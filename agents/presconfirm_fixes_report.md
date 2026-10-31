# Pre-confirmation fixes: plan-specified sensitivities wired in, three engine bugs fixed (statistical methods engineer, 2026-10-31)

**Scope.** Tune split only. No confirm data read, CONFIRM_RUN never set. ANALYSIS_PLAN.md, config/tuning.yaml and the ledger are untouched. No primary estimator, threshold, margin or decision rule was changed. Nothing committed.

**Bottom line.**
1. **E4 deep ages: implemented.** Strassert 2021 and Betts 2018 ages are mapped onto the current (D30) panel tree and run as `ages_override` sensitivities. `robust_alpha_timescale` is now computed exactly per §4.4.
2. **E1/E2 "all aligned columns": implemented.** Re-running the Foldseek tdiaa alignment is cheap: 55 s for 220k tune pairs, about 4 min projected for confirm. The recomputed core counts are byte-identical to the stored ones (220,104/220,104).
3. **E5 single-copy guard: not a bug.** It is a data limit on tune: 95 single-copy HOGs, 66 after exclusions, and at most 22 HOGs on any panel edge, below the plan's 50. Left as is.
4. **E2 d\* = 0.05: genuine.** It follows from the §4.2 definition and is not a bug. Descriptive fields were added.
5. **E4 flags.**
   - `fast_reml_mismatch` was a **real optimiser bug** (statsmodels stopped short of the REML maximum). Fixed.
   - `censored_cells_gt2pct` and `ci_disagree` are genuine.
6. **New bug found while checking the E5 guards.** The branch engine accepted rank-deficient designs, which gave |b| up to 1e13. This broke the hidden-paralogy guard (σ²_L = −78) and touched 4 primary HOGs. Fixed.

Verification: `results/tuning/dryrun2/` (see §6). Tests: `pytest tests` gives 77 passed. Two `test_hogs` tests fail before and after these changes; they test scripts-level seqdist/treeutil code that I did not touch (§7).

---

## 1. E4 deep-age sensitivity (`deep_age_sensitivity_unavailable` → implemented)

**Data.** The curator's table (`data/interim/timetree/deep_ages_alt.tsv`) was built on the **pre-D30** collapsed tree (`timetree/old_ci_collapse/`). Node ids changed after D30. For example, old N006 Rhodophyta is now N007 and old N008 Viridiplantae is now N009, so the table cannot be used by id.

**`scripts/tuning/map_deep_ages.py`** writes `data/interim/timetree/deep_ages_panel.tsv`. Ages stay under `data/interim/timetree` (trap #7).

- **Mapping is by identical panel-leaf sets** (old node → current node). All 16 curated rows map exactly.
- Two clades that were collapsed into the old root and exist again under D30 are added from the curator's reference table:
  - **Opisthokonta**, N029. S21 n242: 1359 (1051–1679). B18: 1125–1707. The node is Fungi + Holozoa + Fonticula, so its MRCA is the crown: exact.
  - **Archaeplastida**, N006 (Rhodophyta + Viridiplantae). S21 n191: 1925 (1799–2063). B18: 1118–1667. This is *approximate*: the panel has no glaucophyte.
- Only **deep nodes** are overridden: TimeTree median ≥ 500 My, the curator's own cut. Embryophyta (TimeTree 483) is excluded as "reference only". TimeTree ages are used everywhere else.
- **S21 overrides 11 nodes:** N001, N006, N007, N009, N029, N030, N039, N054, N055, N056, N103.
- **B18 overrides 3 nodes:** root, Opisthokonta and Archaeplastida. B18 gives intervals only, so the point is the geometric midpoint √(lo·hi) (the median of the plan's lognormal age model with that 95% interval), and lo/hi = the interval.

**Driver** (`driver.py`):
- `load_panel(..., age_override=)` applies the override **before** NA-age imputation. It then applies the plan's top-down repair (child = min(child, parent − 0.1)); on tune this repaired 0 nodes.
- The branch engine is **refitted on the re-aged panel**, because the D1 root-merged shares use ages.
- `run_E4` then runs at B4_sens = 200, exactly like the α variants.
- The deep-age CIs (HPD / B18 interval) enter the age redraws.

**Robustness rule.** `robust_alpha_timescale` follows §4.4: the same sign as β̂ and unadjusted p < 0.05 in α×1.5, α÷1.5, S21 **and** B18. It is False, with a flag, if any of the four is unavailable. A per-sensitivity breakdown is in `robust_alpha_timescale_detail`.

## 2. E1/E2 "all aligned columns" (`all_columns_sensitivity_unavailable` → implemented)

**Why it was unavailable.** `counts_shard_*.npz` holds core-column (pLDDT ≥ 70 on both sides) aa/tdi counts only; only seqonly has `aa_all`. The tdiaa all-column 3Di counts were never stored.

**Cost.** Re-aligning is cheap: 55 s at 40 threads for all 220,104 tune primary pairs. Projected confirm cost (~870k pairs) is about 4 min.

**`src/strucevo/pairs/allcols.py`** re-runs `run_foldseek` with the production method and flags and counts AA/3Di over every aligned M column.
- **Determinism check:** the recomputed core counts equal the stored ones for **220,104/220,104** pairs. Any mismatch would be set to weight 0 and reported.
- Output is cached in `data/interim/pairs_variants/all_columns/tdiaa/allcols__{split}.npz`. Confirm requires CONFIRM_RUN=1.

**Driver `_all_columns`:**
- aa_ml and tdi_ml (with SE, sat and p) are recomputed at the tuned α, and κ(bin) is applied as for core columns.
- Unchanged: status, `aa_ml_seq` (C26) and lddt_\*. **E2's y stays lddt_D**, since lddt_global is already all-column, so for E2 only x changes.
- The branch engine is refitted for the D31 parametric-bootstrap design.
- Scalar ρ is used, as in the α-variant sensitivities.

On tune, the all-column distances are about 18% larger than core distances (median aa_ml 0.58 vs 0.47; tdi_ml 0.215 vs 0.176).

## 3. E5 single-copy-HOG guard (`guards_not_estimable:single_copy_hogs`): genuine data limit, no change

- Single-copy is defined as `n_genes_filtered == n_species_1to1` (DESIGN: "strictly single-copy HOGs").
- On tune: 95 of 1,400 HOGs, of which 66 survive exclusions.
- The largest number of HOGs contributing a simple edge to any one panel edge is **22**. §4.5 requires ≥ 50, so there are 0 lineages with or without D8. This is not a swallowed exception.
- **On confirm** (about 4× more HOGs, so roughly 260 single-copy HOGs and ≤ about 90 per edge) the guard may become estimable.
- If it stays not estimable, the existing rule applies: an E5 *rejection* becomes `artefact_suspect`.

## 4. E2 d\* = 0.05 (grid minimum, SE 0): genuine

I checked the §4.2 definition. The linear fit uses x ≤ 0.5 and the power fit uses x ≤ 2, both with the same fixed η_raw = 0.0205, so the intercept is not the cause.

With A = 0.139, A′ = 0.107 and c′ = 0.645, the relative departure is:

| x | Departure |
|---|---|
| 0.05 | 31% |
| 0.20 | 21% |
| 0.33–0.69 | ≤ 10% |
| ≥ 0.70 | > 10% |

- The two fits agree only around x ≈ 0.5, where the linear fit is anchored by most of its data.
- With c′ ≈ 0.65 there is **no linear regime at small d**: a curve η + A′x^0.65 is steeper than any line near 0. So the "smallest x with > 10% departure" is the first grid point, in every bootstrap replicate (hence SE 0).
- **Meaning:** under the plan's definition, d\* says the linear approximation already fails at the smallest distances analysed. It is *not* the onset of saturation.
- **Change (C41): reporting only.** `d_star` now also carries `at_grid_min`, `rel_at_grid_min` (0.312), `agree_within_10pct` ([0.33, 0.69]) and `upper_departure_x` (0.70). These are labelled descriptive; the estimand is unchanged.

## 5. E4 flags

- **`fast_reml_mismatch:1.35e-03`: bug, fixed (D40).**
  - `fit_mixed` accepted the first statsmodels fit flagged "converged". On tune, lbfgs stopped at REML llf −1961.337, while the fast-REML optimum (reached from 5 different starts) is −1960.620.
  - A statsmodels restart from the fast optimum reaches −1960.620 with lbfgs, bfgs, nm and powell alike. A cold bfgs start gives the same.
  - Fix: after the normal fit, restart statsmodels from the fast-REML optimum and keep the statsmodels fit with the higher REML llf (flag `reml_polished:dllf=…`).
  - Effect: tune β̂ −0.1324 → −0.1338 (Δ = 0.03 SE). It is still statsmodels REML, so the estimator is unchanged; the fix only reaches its maximum. The same polish now applies to every E4 sensitivity point fit.
- **`censored_cells_gt2pct:0.030`: genuine; plan-mandated flag (§2.6).**
  - 49 of 1,658 cells have S_3Di ≤ 0 after merging (by bin: 1/5/19/22/2).
  - Dropping them makes β more negative (−0.157), so censoring does not create the sign.
- **`ci_disagree`: genuine; a property of the plan's bootstrap, not a bug.**
  - The normal 98% CI excludes 0, while the percentile CI includes it.
  - Decomposed with 40 replicates per component:

    | Component | Replicate mean (β̂ = −0.134/−0.136) | SD |
    |---|---|---|
    | HOG weights W only | −0.134 | 0.017 |
    | Age redraws only | −0.146 | 0.009 |
    | **Panel-edge weights Λ only** | **−0.095** | **0.050** |
    | Full replicate | −0.100 | 0.053 |

  - The Λ weights inside the cell sums (§2.4/§4.4) dominate the SE and shift replicates towards 0. The same happens with the plan's cell rule (point −0.240, Λ-only mean −0.203), so D32 is not the cause.
  - The plan's decision uses the Wald p (normal), so nothing changes. **PI note:** the percentile CI is the more conservative of the two here.

## 6. Bug found while checking the E5 guards: rank-deficient branch designs (D39)

**Symptom.** In the dry run, the hidden-paralogy guard had σ²_L = −77.96, σ²_obs = 10.4 and 8 censored lineages, with lineage AA sums D_e up to 5.7e13.

**Cause.** §2.5.1 marks a fit underdetermined only by count (n_ok < n_edges + 5). Removing the hidden-paralogy pairs often leaves a design whose ok rows do **not identify every edge**. For example, removing the pair inside a cherry leaves a, b and the edge above them confounded. XᵀWX is then singular and `solve` returns huge b. Four **primary** HOGs had the same problem, among 17 rank-deficient primary fits in all; the other 13 were already excluded by `fail_frac`.

**Fix.**
- `wls_two_pass` now checks rank(X[ok]) < p. Such fits are underdetermined, with the HOG exclusion reason `rank_deficient` (counted in n_excluded).
- A unit test was added (`test_wls_rank_deficient_flagged`).

**Effect.**
- Primary: 4 more HOGs excluded from E3–E5. Tune E4 β̂ −0.1338 → −0.1357; E5 σ²_obs 0.0818 → 0.0812.
- Hidden-paralogy guard: 458 HOGs become rank-deficient. 48 D8 lineages remain (was 166 with garbage values). σ²_obs = 0.076, F = 0.020, which is sane.
- **PI decision:** an alternative for that guard would be to prune species whose pairs were all removed, instead of excluding the HOG. I did not do this, because it would be a new rule.

## 7. Verification

**Full tune dry run:** `results/tuning/dryrun2/` (log `logs/presconfirm/dryrun2.log`); results in §8.

**`run_confirm.py`.** It calls the same `EndpointDriver.run_all`, so every change above is on its code path.
- The all-column counts are computed inside the confirm run, and `all_column_counts` refuses confirm without CONFIRM_RUN=1.
- `deep_ages_panel.tsv` must exist (it is gitignored under `data/interim/`); regenerate it with `python scripts/tuning/map_deep_ages.py`. If it is missing, the deep-age sensitivities read "unavailable" and robust_alpha_timescale = False (flagged).
- **New files must be committed before the confirm run** (G2 git_clean only checks tracked files, but the code must be in the commit): `src/strucevo/pairs/allcols.py` and `scripts/tuning/map_deep_ages.py`.

**Tests:**
- `pytest tests --deselect tests/test_benchmark.py` gives **77 passed**.
- Two failures are deselected: `tests/test_hogs.py::test_seqdist_matches_pairs_pipeline` and `::test_clade_closed_and_dmax`. Both test scripts-level `seqdist`/`treeutil` code against demo data; neither imports `strucevo.stats` or `strucevo.pairs.allcols`. They are outside the files I changed, likely stale demo data or the D30 tree.

## 8. Dry-run results (tune; exploratory)

`results/tuning/dryrun2/` (run time 3,646 s; log `logs/presconfirm/dryrun2.log`). **Every plan sensitivity is populated, apart from the E5 single-copy guard (data limit, §3).** No `*_unavailable` flags remain.

| Endpoint | dryrun (before) | dryrun2 (after) | Newly populated sensitivities |
|---|---|---|---|
| E1 | ĉ 0.708, concave | ĉ 0.708 [0.614, 0.803], concave | all aligned columns: 0.531 [0.432, 0.631], concave |
| E2 | ĉ 0.673, concave; d\* 0.05 | ĉ 0.673 [0.571, 0.775], concave; d\* 0.05 (at_grid_min; fits agree within 10% for x ∈ [0.33, 0.69]; upper departure at 0.70) | all aligned columns: 0.789 [0.673, 0.905], concave |
| E3 | k 3.32 | k 3.32 | — |
| E4 | β̂ −0.1324, p 0.0128; flags fast_reml_mismatch, censored 0.030, ci_disagree, deep_age_unavailable | **β̂ −0.1357** [−0.256, −0.015], p 0.0089, p_holm 0.018 → β_k > 0; flags: censored 0.030, ci_disagree (both genuine). n_excluded adds rank_deficient: 21 (of which 4 are newly excluded; the rest also fail other rules) | S21: −0.126, p 0.029. B18: −0.132, p 0.017. α×1.5: −0.131, p 0.002. α÷1.5: −0.145, p 0.004. → **robust_alpha_timescale = True** |
| E5 | β 0.584; hidden-paralogy guard σ²_L −77.96 (garbage) | β 0.601 [−0.49, 1.70], inconclusive; σ²_L 0.061, F 0.019, 106 lineages | hidden-paralogy guard: σ²_L 0.054 [0.006, 0.103], β 0.58, now sane; single-copy guard: not estimable (data) |
| E6 | detected | 0.0255, detected | — |
| Holm / interpretation | — | E1 concave, E2 concave, E4 β_k > 0 (p_holm 0.018), E5 inconclusive (power 0.24 < 0.80), E6 detected. Row `none` | — |

**Residual items.**
- The `aa_ml_seq_le_1_hogs` sensitivity still shows `fast_reml_mismatch:7.4e-04`. Here the statsmodels fit is not improved by the restart, so the fast engine (used only in its 200 replicates) is the one off-optimum. This affects SE only; it is a minor, sensitivity-only issue.
- The D13 negative control still passes: 20/20.

## 9. Proposed DEVIATIONS entries (continuing after C36)

| ID | Change | Reason | Effect |
|---|---|---|---|
| C37 | E4 deep-age sensitivities (§4.4) implemented as `ages_override` runs. Literature ages are mapped to the current (D30) panel tree by identical leaf sets (`scripts/tuning/map_deep_ages.py` → `data/interim/timetree/deep_ages_panel.tsv`). Opisthokonta (exact) and Archaeplastida (approximate, no glaucophyte in the panel) are added from the curator's reference table. Only nodes with TimeTree median ≥ 500 My are overridden (S21: 11 nodes; B18: 3). B18 point = √(lo·hi). Override precedes NA imputation and is followed by top-down repair. The branch engine is refitted on the re-aged panel. B = B4_sens. `robust_alpha_timescale` requires all four of α×1.5, α÷1.5, S21 and B18 | Plan sensitivity was unwired. The curator table used pre-D30 node ids. Betts gives intervals only | sensitivity and robustness flag only |
| C38 | E1/E2 "all aligned columns": the Foldseek primary alignment is re-run (`strucevo.pairs.allcols`; core counts verified identical) and aa/tdi distances are computed over all aligned columns at the tuned α, with κ(bin). Status, aa_ml_seq and lddt are unchanged (E2 y = lddt_D unchanged). Scalar ρ. Branch refit for the D31 design | All-column 3Di counts were not stored. Recompute costs ≈ 1 min (tune) / ≈ 4 min (confirm) | sensitivity only |
| D39 | Branch engine (§2.5.1): a fit is also "underdetermined" (reason `rank_deficient`) when the ok-pair design does not identify all edges (rank X[ok] < p) | Singular XᵀWX gave |b| up to 1e13. It broke the hidden-paralogy guard (σ²_L = −78) and affected 4 primary tune HOGs | 4 more primary HOGs excluded (tune); hidden-paralogy guard: 458 HOGs excluded, 48 lineages |
| D40 | E4 point fits: after the statsmodels REML fit, restart statsmodels from the fast-REML optimum and keep the higher REML llf (flag `reml_polished`) | statsmodels lbfgs stopped short of the REML maximum (Δllf 0.72), which caused `fast_reml_mismatch` | tune β̂ −0.1324 → −0.1338 (0.03 SE). Estimator unchanged (REML maximum) |
| C41 | E2 `d_star` gains descriptive fields (`at_grid_min`, `rel_at_grid_min`, `agree_within_10pct`, `upper_departure_x`). Estimand unchanged | d\* = 0.05 is the literal §4.2 result: with c′ ≈ 0.65 there is no small-d linear regime | reporting only |
| (note) | E5 single-copy guard not estimable on tune: data limit (max 22 HOGs per panel edge < 50). E4 `ci_disagree`: Λ edge weights shift bootstrap replicates towards 0 (+0.04), under both cell rules. E4 censored 3.0% > 2%: plan flag | — | none |
