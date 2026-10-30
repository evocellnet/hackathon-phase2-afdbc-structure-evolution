# DESIGN.md — Is protein structure a molecular clock? AFDB clusters revisited with orthology

Version: **v3** (after review rounds 1–2; reviews in `reviews/round{1,2}_*.md`). **§11 (v3 amendments) is binding and supersedes earlier text where they conflict.** · 2026-10-30 · PI: agentic (Claude)
Inputs: PROMPT.md, LITERATURE.md, agents/data_scout_report.md, agents/ancestra_report.md
Change log v1→v2 at the end (§10) — every change cites the review item that motivated it.

## 0. Questions and estimands

The brief fixes two goals. We turn them into **three headline estimands for Goal 1** and a **compact Goal-2 portfolio** whose first item doubles as the positive control for Goal 1.

### Goal 1 — Illergård, Ardell & Elofsson 2009 at AFDB scale, plus a real time axis

Illergård et al. showed that core structural change grows linearly with *ML amino-acid distance* (substitutions/site), not with time. So "linear with time" is a new question, not a replication.

- **H1a — response shape.** How does structural divergence grow with amino-acid (AA) divergence?
  - Estimand: the curvature exponent c of structure-vs-AA distance within ortholog families. A dose (AA distance) to response (structure) curve.
  - Measured for 3Di-ML distance and floor-corrected 3D distances (LDDT, TM).
  - Also d*: the AA distance at which the raw metric departs from linearity by >10%.
- **H1c — alphabet-matched conservation ratio.**
  - Estimand: k = Σ d_AA / Σ d_3Di over branches.
  - Both distances use ML corrections over 20 states (LG+Γ for AA, Q.3Di.AF+Γ for 3Di).
  - Report k with its 95% CI and its between-family spread (SD of log k).
  - Effects of fold class (CATH α, β, α/β), length, pLDDT and cluster age.
  - Comparison: Illergård's 3–10× and β > α.
- **H1b/H1d — time and the clock.**
  - (i) Does structure accumulate linearly with time relative to AA, i.e. is k constant across node ages? Estimand: slope of log k vs log node age; 0 means a proportional clock.
  - (ii) Is structural rate variation among lineages *more than* what AA rate variation predicts? Estimand: the lineage-effect variance of log(b_struct / b_AA) beyond the measured noise floor.
- **Neutral-structural-clock null (N0)** (external PI #3; molevo #4). Neutrality is *not* the same as clock-likeness. N0 states that structural change is a passive by-product of AA change:
  1. Each family's structural branch lengths are proportional to its AA branch lengths, with one family-specific ratio.
  2. That ratio is constant across lineages and node ages.
  3. There is no extra overdispersion beyond noise.

  Departures in the direction of extra structural *conservation* at depth (k rising with age, saturation) are interpreted as constraint. Lineage-specific excess is interpreted as selection or functional change. **Positive control:** N0 should be violated more for paralogs and for the more-diverged copy after duplication than for orthologs (Goal 2-i). If it is not, we cannot claim the test has power. We never claim "neutral" from a good fit alone.
- **H1e — AlphaFold artefact.**
  - Bound the bias of predicted structures using orthologs that have experimental PDB structures, compared PDB–PDB vs AF–AF at matched AA distance.
  - Add the AF noise floor (§4.4).
  - Qualifies all conclusions: they are about *predicted* structure evolution.
- **H1f — phylogenetic utility (if time).** Recovery of the species tree from structural vs AA distances, by node depth.

### Goal 2 — AFDB clusters in the light of evolution (cut to three items, external PI #8)

- **2-i Duplication against the ortholog null.** This merges the old H2b and H2c.
  - (a) Paralog vs ortholog structural residual at matched AA divergence, using pairs that share one protein (molevo #6).
  - (b) Least-diverged-ortholog (LDO) test with structure. For a duplicated pair (P1, P2), choose the more-diverged copy (MDO) with outgroup set O1. Measure its structural excess beyond AA expectation with an independent outgroup set O2 (molevo #7). Use column split-halves to avoid regression to the mean (stats S7).
  - First reproduce the direction of the PMID 40813247 result.
- **2-ii Cluster audit plus domain architecture.**
  - How Ancestra HOGs map onto AFDB v3 structural clusters: splits and merges.
  - Whether HOG splits across clusters are explained by changes in TED domain architecture, length change or low pLDDT.
  - The rate of domain-architecture change in orthologs vs time.
  - Cluster LCA age vs HOG root age.
- **2-iii One case study, time-boxed to 4 h: cytoskeletal folds.**
  - Families: actin/MreB/ParM/crenactin/Arp and tubulin/FtsZ/CetZ, across eukaryotes plus the prokaryote panel (Asgard archaea as the closest outgroup).
  - Measure the structure/AA branch ratio along a family tree, flagging branches with structural acceleration, e.g. at the eukaryote origin.
  - Kinases only if time remains.

## 1. Data

| Item | Source | Use |
|---|---|---|
| AFDB clusters v3 (paper release) | steineggerlab S3 `afdb-cluster` | cluster membership, LCA, dark flag; file 8 member→rep alignments (structure-selected comparison set, P-clu) |
| AFDB v4 structures | FoldComp `afdb_uniprot_v4`, fetched by HTTP range | coordinates, model sequence, pLDDT (v4 matches clusters v3) |
| Eukaryote panel | UniProt reference proteomes 2026_03, canonical; **~150–180 species**; gate: ≥80% of accessions have a v4 model with an identical sequence | Ancestra, Goal 1, Goal 2-i/ii |
| Prokaryote panel | ~100 reference proteomes, incl. Asgard | Goal 2-iii only (separate Ancestra run if time) |
| Divergence times | TimeTree: resolved dated topology from the TimeTree-of-Life 2015 newick, pruned to the panel; pairwise median ages and CIs from the TimeTree 6 API; sensitivity with deep-node ages from Strassert et al. 2021 and Betts et al. 2018 (molevo #5) | time axis, eukaryotes only. **Never committed (licence)** |
| 3Di substitution model | Q.3Di.AF (Edmond doi:10.17617/3.1MJJBH, md5 8f897ff…, CC0) | 3Di ML distances |
| AA model | LG (PAML format) | AA ML distances |
| Domains | TED (AFDB v4 choppings plus CATH labels), Pfam 38 | domain architecture, fold class |
| Experimental structures | PDB via SIFTS UniProt↔PDB mapping | H1e |
| Taxonomy | NCBI taxdump 2026-10-30 | lineages, Ancestra guide tree |

**D1 (revised).** Ancestra runs on the complete 2026_03 canonical proteomes (best orthology). The sequence–structure identity is enforced at the *protein* level: a protein enters any structural analysis only if md5(UniProt 2026 sequence) equals md5(AFDB v4 model sequence). Mismatch rates are reported per species. Proteins over 2,700 aa have no AFDB model, a documented bias.

## 2. Orthology

- **One Ancestra run** (commit ffb0bbe, `accurate`) on the eukaryote panel, guide tree from NCBI taxonomy. Launched on node 1; stage checkpoints make it resumable after the node move.
- Pair labels follow the deepest-common-group rule (`hog_pair_labels.py`).
- **D2:** Ancestra's duplication *ages* are not used (36% agreement with Ensembl). Paralog analyses match on AA distance and use explicit outgroups.
- **Gene-model filters** (molevo #6), applied to all analysed proteins:
  - not a fragment;
  - length within 0.5–2× the HOG median;
  - mean pLDDT ≥ 70 (per protein; pairs are *not* filtered on alignment outcome, stats S1).
- Sensitivity run: restrict to strictly single-copy HOGs.

## 3. Units of analysis and sampling (avoiding structural selection bias)

- **Ortholog family set.** One set per HOG: the 1:1 orthologs, one protein per species, with species that have more than one copy dropped for that HOG.
  - Keep HOGs with at least 8 species.
  - Cap at 40 species per HOG, chosen deterministically (seed 20261030) to maximise coverage of species-tree nodes.
  - Pairs are chosen by orthology only, never by structural similarity.
- **Per-HOG branch lengths on the fixed species tree** (external PI #1; molevo #2). For every metric m we compute all pairwise distances within the HOG set. We then fit branch lengths on the pruned, resolved TimeTree topology by non-negative least squares (Fitch–Margoliash weights), which gives b_m(branch, HOG). This single engine serves response shape, ratio, time and lineage analyses. It also handles shared branches, the main source of non-independence between pairs.
- **Additivity check.** Each metric's tree-fit residual (relative SS) is reported. A metric whose distances are not tree-additive (median relative residual >2× AA's) is reported but not used for clock claims.
- **Tuning vs confirmation split** (stats S11): HOGs are randomly split 20/80 (seed 20261030).
  - Everything tunable is set on the 20%: Γ α for AA and 3Di, metric floors S₀, filters, the curvature fitting procedure.
  - The 80% is then analysed once with frozen code (`ANALYSIS_PLAN.md`, hash recorded in `RESULTS_LEDGER.md`).
- **P-clu (comparison only).** AFDB-cluster member→representative alignments, structure-selected, used to show the selection effect. Not used to test H1.

## 4. Divergence measures

### 4.1 Alignments (molevo #1; algorithms C4)
- **Primary: Foldseek `tmalign`**, a coordinate-only structural alignment that favours neither the AA nor the 3Di letters, as in Illergård's structural-alignment approach. Both AA and 3Di distances for H1a/H1c, and all 3D metrics, come from it.
- **Sensitivity 1:** sequence-only alignment (parasail/MMseqs2 global or semi-global) for AA and 3Di distances.
- **Sensitivity 2:** Foldseek 3Di+AA alignment.
- The 3Di/AA ratio is known to move with alignment (0.56 vs 0.69 in the reviewer's benchmark), so k is reported for all three and the primary is pre-specified.
- All Foldseek pair runs use `-e inf --tmscore-threshold 0`, and exhaustive alignment of the supplied pair list (`tsv2db` → `tmalign`/`structurealign` → `convertalis`). Without this, pairs are silently dropped (algorithms C1).
- **Row-count invariant:** pairs out = pairs in, asserted.

### 4.2 Metrics per pair (post hoc, from alignment, coordinates and pLDDT)
- **Core mask:** residues with pLDDT ≥ 70 in both proteins.
- **AA:** p-distance; ML distance under LG+Γ4 (α tuned on the 20% set).
- **3Di:** p-distance; ML distance under Q.3Di.AF+Γ4 (separate α).
- **3D:**
  - global LDDT (over the full query, unaligned residues scoring 0) and core LDDT;
  - TM-score normalised by the shorter protein;
  - core RMSD (Kabsch over aligned core residues).
- **Corrected 3D distances** (molevo #2): D = −ln((S − S₀(L)) / (1 − S₀(L))), where S₀(L) is the length-dependent floor estimated from ~100k random cross-HOG pairs passed through the same pipeline (algorithms C5). Raw 1−S is kept as a secondary metric.
- **Saturation:** ML optimum at t_max = 20, or S ≤ S₀ + 2 SE. Such pairs are flagged, kept as censored, and never dropped silently.
- **Coverage rule:** ≥50 core aligned residues *and* coverage ≥ 0.5 of the shorter protein for ML distances. Failures are counted and reported per AA-distance bin (stats S1), with sensitivity analyses that include them as censored.

### 4.3 GPU ML distance engine
`src/strucevo/distances.py`:
- builds a 20×20 count matrix per pair;
- runs a batched ML fit on GPU with eigendecomposition and discrete Γ;
- reports SE from the Fisher information.

It is validated against AliSim and the IQ-TREE `.mldist` output (tests in `tests/`).

### 4.4 Noise floor (external PI #4)
- **(a) Estimation noise.** Per-pair column bootstrap (100 replicates on a subsample) gives an SE per metric.
- **(b) AF prediction noise.** Pairs of *identical* sequences with different accessions (e.g. strains) in AFDB v4: the 3D/3Di distance at d_AA = 0.
- **(c) Random floor S₀.**
- The clock comparison between metrics is made after subtracting (a) and (b).

## 5. Analyses and decision rules

### Primary endpoints
Confirmatory, on the 80% HOG set, with Holm correction across E1–E6 (stats §3).

| ID | Hypothesis | Estimand | Decision rule |
|---|---|---|---|
| E1 | H1a | curvature exponent c of HOG-level 3Di-ML vs AA-ML branch lengths, d_AA ≤ 2 | "linear" if 90% CI ⊂ [0.9, 1.1]; "concave (structure saturates)" if CI < 0.9. Corrected for errors in variables (SIMEX using per-branch SEs) |
| E2 | H1a | same as E1 for floor-corrected global LDDT distance | same rule; also d* for raw 1−LDDT, with CI |
| E3 | H1c | k = Σ b_AA / Σ b_3Di; SD of log k across HOGs; CATH class effect | report with 95% CI. Compare to 3–10× |
| E4 | H1b | slope of log(b_3Di / b_AA) vs log node age (TimeTree), mixed model with HOG random effect | 0 means proportional clock; CI excluding 0 means structure/AA ratio changes with depth |
| E5 | H1d / N0 | lineage-effect SD of log(b_struct / b_AA) minus the noise floor; slope of Δstructure on ΔAA in sister-lineage contrasts | slope 1 and excess SD ≈ 0 mean structure is slaved to sequence (N0) |
| E6 | 2-i | paralog − ortholog structural residual (3Di-ML, LDDT-corrected) at sample-split-matched d_AA; MDO − LDO structural excess with independent outgroups | > 0 means the positive control is detected |

Everything else is exploratory, with unadjusted CIs, and labelled as such: TM, RMSD, p-distances, aligned-only metrics, raw-metric time fits, time-curvature fits, H1e, H1f, 2-ii and 2-iii.

### Uncertainty
- Two-way bootstrap that resamples HOGs and species, 1,000 replicates (stats S3).
- Bootstrap replicate seeds are derived from 20261030.

### Linearity in time (exploratory figure plus E4)
- Branch-length sums from root to tip, or path distances between species pairs, against TimeTree age, per metric.
- AA itself is expected to be non-linear in time (rate variation), so the **curvature of structure relative to AA** is what is interpreted (external PI #2).
- Node ages: point estimates plus 200 draws from TimeTree CIs.
- Nodes ≥1 Ga are reported separately (external PI #7).

### Errors-in-variables check (stats S2)
- Simulate a truly linear structure response with our measured noise.
- Report how often the fitting procedure calls it "saturating".
- Report the false-curvature rate alongside E1/E2.

### Relative-rate tests (exploratory)
- Tajima/Wu–Li-style triplet tests, applied to AA and 3Di only (count-based models).
- Compared by effect size, not rejection counts (molevo #3; stats S6).

### H1e
- SIFTS ortholog pairs with experimental structures (X-ray ≤ 3 Å, coverage ≥ 70%) on both sides.
- Compare PDB–PDB vs AF–AF distances for the same pairs: paired differences, and slope vs d_AA.

### Goal 2-ii
- Contingency tables of HOG ↔ cluster mapping.
- Logistic model: HOG split across clusters ~ TED architecture change + length ratio + pLDDT + d_AA.

### Goal 2-iii
- Family extraction by Pfam/TED from both panels.
- AA and FoldMason 3Di MSAs; trees with IQ-TREE (LG / Q.3Di.AF).
- Branch-wise b_3Di / b_AA; accelerated branches are flagged against the HOG-level distribution of the ratio. Exploratory.

## 6. Controls and verification (written before results; executable)

`tests/` (pytest) and `scripts/verify/*.py`, all run by `make verify`. Each has a pass threshold.

1. **ML engine:**
   - AliSim LG+Γ and Q.3Di.AF+Γ: bias <3% for d ≤ 2;
   - agreement with IQ-TREE `.mldist` within 2%;
   - identical sequences → 0;
   - random sequences → flagged saturated.
2. **Tree-fit engine:** simulated additive distances on a known tree → branch lengths recovered (r > 0.99); added noise → unbiased.
3. **Curvature procedure:** synthetic linear, concave and convex data with realistic x-noise → classification accuracy reported. False-curvature rate on linear data must be ≤10% after the errors-in-variables correction, or E1/E2 are reported as inconclusive.
4. **Pipeline invariants** (asserted in code):
   - pairs out = pairs in;
   - self-pairs give TM = 1 and LDDT = 1;
   - identical-sequence pairs have d_AA = 0;
   - md5(sequence) of model equals that of UniProt for analysed proteins;
   - TimeTree ages in Mya, with human–mouse ≈ 87 ± 10;
   - human–mouse 1:1 ortholog median AA identity ≈ 85%.
5. **Negative controls:**
   - random non-homologous pairs give the floors;
   - shuffled node ages give E4 slope ≈ 0.
6. **Selection-effect demonstration:** P-clu vs ortholog curves.
7. **Independent review:** every E1–E6 result and every report figure is re-derived or audited by an agent that did not produce it, before it enters the report.

## 7. Process safeguards (AI-failure-mode review)

- **`ANALYSIS_PLAN.md`** is frozen before the 80% confirmation run; its git hash goes into `RESULTS_LEDGER.md`. The ledger has one line per reported number: value, script, commit, input md5.
- **No silent fallbacks.** Scripts fail loudly on missing data. Every filter logs counts in → out to `logs/filters/*.tsv`.
- **Licence.** TimeTree-derived tables live only in `data/interim/timetree/`, which is gitignored. A pre-commit check rejects files containing a `timetree` column header.
- **`.gitignore` fixed.** Large or raw data are ignored by directory, not by extension. Small result tables (`results/**.tsv|parquet` under 20 MB) are committed.
- **Node move protocol (`NODE_MOVE.md`).** Before 20:30:
  - stop long jobs at a checkpoint;
  - commit;
  - write a resume checklist (which jobs to restart, with which commands, and how to check completion sentinels).

  `scripts/elapsed.py` handles multiple START entries.
- **Sharded pair computations** with `.done` sentinels (row counts, md5). Reruns skip completed shards.

## 8. Compute plan

| Step | When | Est. wall |
|---|---|---|
| Euk panel, coverage gate, Ancestra input | node 1, now | ≤1 h |
| **Ancestra (euk ~150–180 proteomes, accurate)** | node 1, launch by ~13:00 | 3–6 h, checkpointed |
| FoldComp fetch plus Foldseek DB (createdb `--input-format 5`) | node 1, in parallel | 1–2 h |
| Distance engine, tree-fit engine and verification tests | node 1 | 2–3 h |
| P-clu baseline; noise floors; S₀ random pairs; SIFTS/PDB set | node 1 | 1–2 h |
| Pair metrics (tmalign ~265 pairs/s; ~3M pairs) | node 2 | ~3–4 h |
| Tree fits, E1–E6 (20% tuning, then 80%) | node 2 | 3–4 h |
| Goal 2-ii, 2-iii, H1e, H1f | node 2 | 8–10 h |
| Independent result review, report | node 2 | 8–10 h |

Seeds: global 20261030. Bootstrap and replicate seed list: `config/seeds.yaml` = [20261030, 1, 2, 3, 4, 5].

## 9. Risks and fallbacks

- **Ancestra exceeds the node-1 window.** Resume from its stage checkpoint on node 2. If it fails, fall back to the QfO-78 eukaryotes (already shown to run in less than 1 h).
- **TimeTree topology misses panel species.** Use the closest congener in TimeTree 2015 (documented mapping); otherwise drop the species from time analyses only.
- **Large AF artefact.** Conclusions are restricted to predicted structures; H1e gives the bound.
- **Underpowered positive control (E6 null).** No N0 claims (stated in advance).

## 10. Change log v1 → v2

| Change | Source |
|---|---|
| Per-HOG branch lengths on a fixed species tree replace pair medians and triplets | external PI #1; molevo #2; stats S3 |
| Pre-specified E1–E6 with decision rules; AIC/BIC dropped | stats S4, S11; external PI #2 |
| N0 neutral-clock null with paralogs and MDO as positive control | external PI #3; molevo #4 |
| Noise floor (estimation and AF prediction) | external PI #4 |
| TM-align primary alignment; sequence-only and 3Di+AA as sensitivity | molevo #1; algorithms C4; AI-FM |
| Floor-corrected 3D distances with S₀(L) from random pairs | molevo #2; algorithms C5 |
| Removed the "≥50 aligned residues" pair filter as a hidden selection; pLDDT filtered per protein | stats S1 |
| Errors-in-variables simulation and SIMEX | stats S2 |
| Two-way (HOG × species) bootstrap; 20/80 split; Holm correction | stats S3, S11 |
| Clock tests by effect size; label-swap null dropped; triplet tests restricted to count models | molevo #3; stats S6 |
| Resolved TimeTree topology; age CIs; alternative deep timescales; ≥1 Ga reported separately | molevo #5; external PI #7 |
| Pairs sharing one protein, gene-model filters, single-copy sensitivity, Peterson 2009 credited | molevo #6 |
| LDO with independent outgroups and split-half columns | molevo #7; stats S7 |
| Eukaryote panel ~150–180 with a coverage gate; Ancestra on 2026 proteomes; D1 enforced per protein | external PI #6; algorithms C2/C7; AI-FM |
| `-e inf`, row-count invariants, createdb from FoldComp v4 | algorithms C1–C3 |
| Goal 2 cut to 3 items | external PI #8 |
| ANALYSIS_PLAN, RESULTS_LEDGER, NODE_MOVE, .gitignore fix, licence guard | AI-FM review |

## 11. v3 amendments (binding) — after review round 2

Round-2 reviews: `reviews/round2_external_PI.md` (R1–R6), `reviews/round2_statistics.md` (M1–M8), `reviews/round2_molevo.md` (1–6).

### 11.1 Two estimation engines, each with a defined job (R1, M1, M3, molevo-1)
- **Pairwise engine** (E1, E2, Goal 2-i).
  - Input: observed within-HOG pairwise distances.
  - Curvature model: D_s = η + a_h · D_AA^c, with HOG-specific a_h as a random effect, d_AA ≤ 2, and one global c.
  - Correction: pair-level SIMEX using per-pair column-bootstrap SEs. It handles correlated x/y errors because AA and 3Di come from the same alignment.
  - Saturation shows up in total divergence, so curvature is never estimated on branch increments.
- **Branch engine** (E3, E4, E5).
  - Fit: plain, unconstrained weighted least squares on the fixed, collapsed species topology, with weights 1/(SE²_pair + s²_HOG).
  - Branch SEs: sandwich estimator.
  - Ratios are always **ratios of sums** over cells (HOG × lineage, HOG × age bin), never per-branch ratios. No minimum-branch-length cut.
  - **Tip-noise correction:** subtract the per-protein AF noise floor (from §4.4b and the identical-sequence pairs) from terminal branches. Terminal-branch estimates are also reported separately (M2).

### 11.2 Endpoints (replaces §5 table)

| ID | Estimand | Null / decision (pre-specified) |
|---|---|---|
| E1 | global c for 3Di-ML vs AA-ML (pairwise engine) | **Linear** = TOST equivalence, c ∈ [0.9, 1.1]. **Concave** = one-sided c < 1. Otherwise "inconclusive" |
| E2 | same as E1 for floor-corrected global-LDDT distance; plus d* for raw 1−LDDT | as E1 |
| E3 | k = Σb_AA / Σb_3Di (all HOGs), SD of log k across HOGs | estimation only (outside Holm); only its CATH-class contrast is tested (α vs β, two-sided) |
| E4 | mixed model on y = log(Σb_3Di / Σb_AA) per HOG × age bin, with log age as fixed effect, HOG random intercept and slope, adjusted for mean b_AA of the bin (R4) | slope = 0, two-sided |
| E5 | per-lineage R_e = Σ_h b_3Di / Σ_h b_AA, residualised on the E4 fit and pLDDT. σ²_L = Var(log R_e) − mean bootstrap variance. Disattenuated sister-contrast slope β | **N0 accepted** only if the β 90% CI ⊂ [0.8, 1.25] **and** the σ_L upper 95% bound < 0.10, **and** E6 rejects, **and** spike-in power ≥ 80% (11.4) |
| E6 | (a) paralog − ortholog residual of 3Di-ML at split-half-matched d_AA; (b) MDO − LDO 3Di excess with independent outgroups | one-sided > 0 |

**Multiplicity**
- Holm is applied over the tests of E1, E2, E4, E5 and E6.
- CI-based rules use Bonferroni-level CIs (1 − 2α/m).
- Everything else is exploratory.

### 11.3 Outcome → interpretation (pre-registered; molevo-2)
This table is used only if spike-in power ≥ 80%; otherwise the outcome is "inconclusive".

| E1/E2 | E4 | E5 | Interpretation |
|---|---|---|---|
| linear | ≈ 0 | N0 margins met | structure follows a rescaled AA clock: "proportional", **not** "neutral" |
| concave | > 0 | N0 met | saturation or under-correction; structure is clock-like only up to d* |
| linear | > 0, robust to α and timescale | any | epoch effect: deep branches are structurally conservative |
| any | < 0 | any | AA under-correction or early bursts; recheck with d_AA ≤ 1 |
| any | any | β < 0.8 or σ_L excess, robust to guards | structure-specific lineage rates; **N0 rejected** |
| any | any | excess not robust to guards | AF-quality or orthology artefact |

- "Rising k means constraint" is **removed**: constraint slows both alphabets.
- E6 > 0 with E5 meeting N0 is interpreted as: duplication, not ortholog lineage, changes the structure change per AA change.

### 11.4 Power: spike-in positive control for E5 (molevo-3)
- On the 20% tuning HOGs, swap a known ortholog for an out-paralog in a random 10% of lineages, giving a known structural perturbation.
- E5's ability to detect it (power) is estimated over 50 replicates.
- Power below 80% means E5 is reported as inconclusive.
- The errors-in-variables simulation (stats S2/S10) reports the minimal detectable departure of c from 1 and of the E4 slope from 0.

### 11.5 Species tree vs gene trees (molevo-4)
- **Short internodes:** collapse internodes shorter than 5 My or with contested placement (TimeTree CI overlap with the parent node) into polytomies before fitting.
- **HOGs that fit their own gene tree better:**
  - Per HOG, compare the AA least-squares fit on the species tree with an unconstrained FastME AA tree.
  - Exclude HOGs where the species-tree fit is >3× worse in relative SS; the threshold is set on the tuning set.
  - Report the excluded fraction.
- **Hidden paralogy:** pairs with d_AA more than 3 MAD above the HOG's expectation at their node age are flagged. A sensitivity analysis removes them.
- **E5 noise floor:** comes from a parametric bootstrap (simulate under the fitted tree plus measured noise).

### 11.6 Pair budget and HOG/pair exclusion rules (R2, R5, M8)
- **Size limits:** at most 30 species per HOG, chosen on AA information only (largest clade-closed set with all D_AA ≤ 2, maximal node coverage, seed 20261030). Up to ~7,000 HOGs, depth-stratified. Primary plus E6 pairs ≤ 3.5M.
- **Failed pairs:** structural failures within a chosen subset get weight 0 (sensitivity: imputed at the floor).
- **Exclusions:** a HOG × metric is excluded if >10% of its pairs fail. Exclusion rates are reported by root age.
- **No selective drops:** pairs are never dropped individually on a structural outcome.
- **Alignment robustness:** a result counts as "alignment-robust" only if the sign and decision of E1–E6 agree between TM-align (primary) and both sensitivity alignments.

### 11.7 Goal 2-i rules (molevo-5)
- **Paralog triplet:** two in-HOG paralogs P1, P2 in the same species, plus the closest outgroup species ortholog(s) from a sister clade (outgroup set O1).
- **Paralog vs ortholog pairs:** the ortholog pair shares P1 with the paralog pair (P1–O vs P1–P2).
- **MDO choice:** pick the more-diverged copy (MDO) with O1 on alignment columns from half A. Measure the excess with outgroup set O2, a different sister-clade species, on columns from half B.
- **Excluded:** pairs above 95% AA identity, and recent in-paralogs.
- **Primary metric:** 3Di-ML. LDDT-corrected distance is secondary.
- **Check:** reproduce the direction of the PMID 40813247 result first.

### 11.8 Case study (molevo-6)
- **Tree:** one fixed reference topology for the actin and tubulin superfamilies (AA ML tree, IQ-TREE LG+G4, UFBoot ≥ 95 nodes only).
- **Estimates:** b_3Di and b_AA come from the same branch engine on that topology.
- **Focal question:** the eukaryote stem relative to the Asgard homologs.
- **Descriptive only:** MreB/FtsZ/ParM, without tests.
- **Time box:** 4 h.

### 11.9 Ancestra contingency (R6)
- Use the `fast` preset if the projected `accurate` wall time is over 5 h.
- If stage S2 has not finished by 19:30, kill it cleanly. It is relaunched first on node 2, because it cannot resume mid-stage; NODE_MOVE.md records this.

### 11.10 H1e additions (stats S9)
- Restrict to PDB entries released after 2018-04-30, which AF2 did not train on, so the comparison with AF predictions is not circular.
- Add the same-protein AF-vs-PDB distance as the floor.
