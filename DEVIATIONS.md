# DEVIATIONS from ANALYSIS_PLAN.md (after freeze) — date | change | reason | effect

## Pre-freeze amendments (2026-10-30 ~12:40; before any real-data endpoint was computed; justified by synthetic verification only)
Source: agents/stats_engine_report.md §4. PI decisions:
| ID | Change | Reason | Effect |
|---|---|---|---|
| D1 | Root-adjacent edges of a bifurcating HOG root are merged into one parameter, apportioned by time span; flagged `root_merged`, excluded from E5 lineages | not identifiable from pairwise distances | none on identifiable quantities |
| D2 | Per-edge `b_se` = model-based (XᵀWX)⁻¹ with w = 1/(SE²_cal + s²) instead of HC1 | HC1 under-covers (0.85–0.88) on high-leverage tree designs (V2); model-based gives 0.96 | affects only E3 v_h and the E4 WLS sensitivity; primary inference stays the pigeonhole bootstrap |
| D4 | In-house NJ (same splits as dendropy); failed AA pairs filled with species-tree fitted values for NJ input only | speed; NJ needs a complete matrix | misfit screen only |
| D5 | Pass-2 weights held at full-data values inside SIMEX/bootstrap replicates; fast windowed grid in replicates | speed; equivalent to 1e-6 | none |
| D6 | Fast in-house REML for E4 bootstrap replicates (matches statsmodels to 1e-7); statsmodels for point estimates | speed | none |
| D8 | E5 lineage eligibility adds a precision rule computed without the lineage's own 3Di sum: r̄·D_e/√ΣSE²_3Di ≥ 8 (threshold to be confirmed on tuning set by spike-in power and null rate) | the plan's rule admits near-zero-3Di young lineages (kurtosis 68) and breaks the σ_L floor; V4 spike missed without it | E5 detects the V4 spike; fewer lineages |
| C9 | β CI uses contrast-level multinomial weights (more conservative); HOG-only CI also reported | PI choice: conservative | wider CI |
| D13 | E4 negative control permutes edge midpoints τ across HOG edges instead of node ages | node-age permutation leaks age through tips (16/20 fail) | V6 passes 20/20 |
| C14 | "Failures imputed at cap" in E1/E2 sensitivity = the HOG's maximum eligible y (SE = median eligible SE) | the plan's cap (d_max = 20) is extreme inside a d ≤ 2 fit | sensitivity only |
| D18 | V4 N0 world is judged on the E5 margins plus spike detection; E6 positivity is verified in a separate duplication world | E6 is null by construction in an N0 world | test definition only |
| C19–C21 | Implementation choices for the interpretation table, E6 and s² as listed in the report; PMID 40813247 direction = MDO shows MORE structural divergence (LITERATURE.md §0.10) → `replication_direction` expected positive | gaps in plan | none |

**PI note (not a change):** synthetic worlds show the E5 β 98% CI half-width is ≈0.24–0.28, so N0 *acceptance* within [0.8, 1.25] is practically unreachable even under exact N0. The margins are NOT widened (no changing tests after the fact). E5 will be reported as inconclusive unless it rejects. Lineage AA-rate spread will be measured on the tuning set and reported.
