# Review round 1 — statistics (phylogenetic comparative methods, model selection)

Reviewer: statistics agent · 2026-10-30 · DESIGN.md v1

## 1. Verdict

**Major revision.** The questions and axis separation are sound. Five features can produce headline results artificially: (i) x-axis noise and alignment-success filtering bend curves toward "saturation"; (ii) the HOG mix changes along x; (iii) dispersion and relative-rate comparisons are not scale- or power-matched across metrics; (iv) MDO is selected and conditioned on the same noisy AA distance; (v) the metric × transform × model × subset space is not fixed (forking paths). Fixes need no new data, only changed estimands and decision rules.

## 2. Change proposals

**S1 (must) — Remove the structural selection hidden in the inclusion filter.** "≥ 50 aligned residues at pLDDT ≥ 70" depends on alignment success, which falls with divergence: the P-clu truncation §3 rejects, reintroduced. Filter pLDDT per protein only; keep failed alignments as censored (floor value or Tobit-type likelihood); report dropout vs AA distance. Make the primary 3D metric coverage-aware (TM by shorter *full* length; LDDT over all high-pLDDT query residues): aligned-only metrics measure a core that shrinks with divergence.

**S2 (must) — Errors-in-variables on the x-axis.** Var(d̂_AA) grows steeply with d, so regression dilution flattens the upper range, mimicking saturation; misalignment also correlates x and y errors (raises d_AA, lowers TM/LDDT). Required: (a) per-pair Var(d̂) from observed Fisher information (free in the GPU ML code); (b) pre-specified primary range d_AA ≤ 2.0, upper range secondary; (c) SIMEX or ODR with supplied error variances; (d) **simulation calibration**: truly linear response with the real x/y-noise, dropout and HOG mix, run through the full pipeline, report how often "saturating" wins. Treat TimeTree CIs as error SD on log-age.

**S3 (must) — Composition along x and the dependence structure.** Deep nodes are covered only by old, broad, usually slow HOGs, so pooled curves mix families along x (Simpson-type curvature). Primary analysis: within-HOG estimation (HOG random slope and curvature), restricted to HOGs spanning ≥ 3 depth bins. The HOG bootstrap handles within-HOG dependence but not cross-HOG dependence via shared species (fast lineages, species-specific AF quality): use crossed random effects (HOG, node/species) or a two-way pigeonhole bootstrap. The k = 3 pairs per node should use disjoint species; summarise by the mean of transformed distance, not the median of raw. **Should:** per-HOG branch lengths on the fixed species topology (ML for AA/3Di; non-negative least squares for TM/LDDT distances), analysing branch-wise ratios and root-to-node depths, which removes double-counted shared paths.

**S4 (must) — Replace AIC/BIC as the decision rule.** With dependent aggregated points the likelihood is not a true likelihood and BIC's *n* is undefined. Bounded raw metrics must curve, so "linear?" is meaningful only as (a) linearity of *floor-corrected* d = −ln((S−S₀)/(1−S₀)), S₀ from negative controls, and (b) location of d* on raw metrics. Decision rule: one curvature parameter (power exponent c) with two-way bootstrap CI; "linear" only if TOST equivalence gives c ∈ [0.9, 1.1]. Supporting: HOG-blocked out-of-range prediction split at a pre-fixed d_AA (not the median, which sits among shallow nodes); parametric-bootstrap LRT. Model variance explicitly (beta/quasi-binomial for S ∈ [0,1]; WLS with fitted variance function for distances). Keep the intercept as a tested prediction (time-dependency, AF floor; LIT §3 H-C3).

**S5 (must) — Estimator for the ratio (H1c).** Per-HOG TLS assumes equal error variances (false); ratios of noisy slopes are heavy-tailed. Primary: ratio of sums k = Σd_AA/Σd_3Di over the unsaturated range, two-way bootstrap CI; between-family variation via a hierarchical model on log k. Report the *local* ratio vs d_AA (drift = saturation diagnostic). Estimate Q.3Di and Γ α on a **held-out 20% of HOGs**, never on tested HOGs; report sensitivity of k to α.

**S6 (must) — Relative-rate and dispersion tests (H1d).** (a) Var(d_AO − d_BO) = Var(d_AO)+Var(d_BO)−2Cov: estimate by block (~10-residue) bootstrap over the triplet alignment for all metrics; one pre-specified triplet per (HOG, node), BH across HOGs. "Swapping ingroup labels" is just a sign flip; drop it. (b) **Rejection fractions are not comparable across metrics**: AA has more information, so structure looks "more clock-like" from low power alone. Compare standardised effect sizes, or regress Δstructure on ΔAA (slope 1 = slaved). (c) var/mean is scale-dependent and meaningful only for counts; use CV of root-to-tip lengths (S3) or dispersion on ML branch lengths. (d) Disattenuate lineage-effect correlations via split-half HOG reliability.

**S7 (must) — H2b/H2c matching on a noisy covariate.** Paralogs have larger true d_AA; at matched *observed* d_AA they still have larger *true* d_AA (regression to the mean), producing spurious structural "excess". In H2c, picking MDO by larger observed d_AO is a winner's curse, and correlated alignment errors bias in an unknown direction. Fix: **sample-split** — define MDO and the matching covariate from half the alignment columns (or multi-outgroup HOG distances), evaluate excess on independent information, apply S2's correction, match within HOG, require dose–response with asymmetry, calibrate the null by simulation under constant structure/AA ratio. Report global and per-domain metrics separately.

**S8 (should) — TimeTree circularity.** TimeTree ages are largely sequence-clock estimates, so AA-vs-time linearity is partly built in. Frame H1b as "structure vs consensus molecular time"; its non-redundant test is constancy of k across node ages.

**S9 (should) — H1e.** Prefer PDB entries released after the AF2 training cutoff (2018-04-30); calibrate with same-protein AF-vs-PDB (conformational noise floor); report an estimate with CI, not pass/fail.

**S10 (should) — Power.** For H1b effective N = distinct nodes (~150–200, mostly shallow); > 1 Ga nodes, where saturation lives, are few with wide CIs. Depth-balance the panel (≥ 10 nodes per log-age bin); use S2's simulation to report the minimal detectable curvature. H1e: ~150 pairs give 80% power for ΔTM = 0.02 at SD 0.08. H2c: report asymmetry-filtered N before testing (LDO strict set was 0.2%).

**S11 (must) — Pre-registration and multiplicity.** Seeded 20/80 HOG split: 20% for tuning (α, Q.3Di, thresholds, code), 80% for confirmatory endpoints only, analysis frozen before unblinding; Holm across the confirmatory family. Exploratory: supergroup splits, TM-align sensitivity, H1f, H2a, H2d, H2e, P-clu.

## 3. Recommended primary endpoints (confirmatory, 80% HOG set, Holm across E1–E6)

- **E1 (H1a):** Curvature exponent c of the HOG-level response of 3Di-ML distance vs AA-ML (LG+Γ4) distance, d_AA ≤ 2.0, corrected for errors in variables. Linear if the 90% CI ⊂ [0.9, 1.1].
- **E2 (H1a, 3D):** Same as E1 for floor-corrected LDDT distance (all high-pLDDT residues, S₀ from random pairs). Also d* for raw 1−LDDT with CI.
- **E3 (H1c):** k = Σd_AA/Σd_3Di (95% CI) plus the between-HOG SD of log k. Test the secondary claim α vs β (CATH class) as a fixed effect.
- **E4 (H1b/H1d):** Slope of log k against log node age; 0 = proportional clock.
- **E5 (H1d):** Slope of ΔStructure (3Di-ML) on ΔAA in triplet relative-rate contrasts; 1 = structure slaved to sequence; residual variance ratio vs AA.
- **E6 (H2b/H2c):** Paralog − ortholog structural residual (3Di-ML and LDDT) at sample-split-matched d_AA; MDO − LDO structural excess with dose–response.

Everything else (TM, RMSD, p-distances, aligned-only metrics, raw-metric time fits, other hypotheses) is exploratory with unadjusted CIs.
