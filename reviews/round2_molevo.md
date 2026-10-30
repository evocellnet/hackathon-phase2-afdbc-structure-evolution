# Round 2 review: molecular evolution (DESIGN v2)

**Verdict: accept with the must-changes in §3.**

## 1. Round-1 items

| Item | Adequate in v2? |
|---|---|
| M1 independent AA axis | Partly. TM-align for everything is defensible, but alignment moves k (0.69 vs 0.56). See S-c |
| M2 floors, additivity | Yes |
| M3 dispersion and relative-rate tests | Yes |
| M4 neutrality logic | Partly. The E4 reading is wrong and the control does not match E5 |
| M5 topology and time | Yes. Age monotonicity is unchecked (S-b) |
| M6 paralog design | Yes, but the triplet rule is missing |
| M7 LDO | Partly. No gene-conversion exclusion; O is not required to lie outside the duplication |
| S1, S3, S4 | Yes |
| S2 | Partly (S-a) |

## 2. Scrutiny

### 2.1 N0 and the positive control

- **Rising k ≠ constraint.** Under passive response, constraint lowers both rates proportionally and leaves k constant. A rising k comes instead from:
  - 3Di under-correction (sites depend on neighbours; the fold bounds the state space);
  - NNLS shrinking deep branches of a less additive metric;
  - errors in deep TimeTree ages.
- **Saturation is cumulative.** Curvature fitted on branch increments (E1/E2 as written) mixes saturation with tree shape.
- **The control misses E5.** E6 is a pair residual, E5 a variance of branch ratios. MDO acceleration may also be relaxed constraint, which is N0-consistent, so an E6 null is ambiguous.
- **No margins.** "≈ 0" and "slope 1" are undefined.

**Outcome → interpretation** (needs spike-in power ≥ 80%, else inconclusive):

| E1/E2 (pairwise c) | E4 (k vs age) | E5 (slope; excess SD) | Interpretation |
|---|---|---|---|
| Linear | ≈ 0 | ≈ 1; ≈ 0 | Rescaled AA clock. Say "proportional", not "neutral" |
| Concave | > 0 | ≈ 1; ≈ 0 | Saturation or under-correction; clock only up to d* |
| Linear | > 0 (robust to α and timescale) | Any | Epoch effect: deep branches structurally conservative |
| Any | < 0 | Any | AA under-correction or early bursts; check d_AA ≤ 1 |
| Any | Any | < 1 or excess > 0; robust to guards | Structure-specific lineage rates; N0 rejected |
| Any | Any | Excess > 0, not robust | AF-quality or orthology artefact |

E6 > 0 with E5 ≈ N0: duplication, not ortholog lineages, changes structure per AA change.

### 2.2 Gene trees vs the species tree

- **Why gene trees differ:** ILS on short internodes (Neoaves, placental base, eukaryote supergroups); hidden paralogy via differential loss (Natsidis 2021); rare protist HGT.
- **Effect:** misfit differs between metrics of unequal additivity, inflating E5's excess SD.

### 2.3 Goal 2-i without duplication ages

Implementable from `hog_pair_labels.py` and the OrthoXML nesting.
- **(a) Triplet (A1, B1, B2).** A1–B1 is a 1:1 ortholog. A1–B2's deepest group is a paralogGroup at least one node deeper than the A/B split. Timing errors near the split bias toward the null; the margin limits this.
- **(b) LDO.** For paralogGroup L (P1, P2), outgroups are species outside L labelled ortholog to *both* copies. O₁ and O₂ come from disjoint sub-clades. Require ≥2 species in L with both copies and ≥2 species per outgroup set. The contrast is age-matched by construction.
- O₁/O₂ removes only outgroup noise. P1/P2 noise needs column split-halves: valid for 3Di-ML, not cleanly for LDDT.

### 2.4 Cytoskeleton case study: ill-posed as written

- MreB/ParM/actin and FtsZ/tubulin are ancient paralogs.
- Tubulin–FtsZ identity is <20%, beyond E1's range.
- Eukaryotic crown branches are nearly invariant: ratios have near-zero denominators.
- AF conformations confound global metrics (actin cleft, tubulin curvature).
- Its reference distribution uses a different estimator (NNLS vs IQ-TREE).

## 3. Remaining changes

**Must**
1. **E1/E2 on cumulative distances.** Fit c on within-HOG pairwise or patristic distances (HOG random slope, two-way bootstrap). Keep branches for E3–E5; E4 age = branch midpoint.
2. **N0.** Drop "rising k = constraint" and adopt the §2.1 table. Freeze TOST margins on the 20% set (e.g. |E4| ≤ 0.05; E5 slope 1 ± 0.1; excess SD ≤ 0.1 log units).
3. **Spike-in control.** In tuning HOGs, swap one species' ortholog for its out-paralog or MDO. Report the detection rate and the minimal detectable lineage effect. Below 80% power, make no N0 statement.
4. **Discordance guards.**
   - Collapse internodes <10 Myr and contested nodes to polytomies; exclude them from E4/E5.
   - Exclude HOGs whose species-tree NNLS SS exceeds the BIONJ AA-tree SS beyond the parametric-bootstrap 95th percentile; report the count.
   - Show E5 without outlier terminal AA branches (hidden paralogy).
   - Take E5's noise floor from a parametric bootstrap with metric-specific residual variance.
5. **2-i rules.** Adopt the §2.3 definitions. Exclude P1–P2 pairs >95% identical (gene conversion). Primary E6(b) metric: 3Di-ML on split halves; LDDT secondary.
6. **2-iii reframe.**
   - One fixed AA topology for both alphabets, with the reference from the same estimator (pooled crown branches).
   - Stem defined by Asgard actins and tubulins; MreB/FtsZ/ParM descriptive only.
   - Core or per-domain metrics.
   - Verify the Asgard proteomes and their v4 coverage first.

**Should**
- **S-a.** Add species mean pLDDT (and MSA depth, if cheap) as an E5 covariate.
- **S-b.** Check that API ages are monotone on the 2015 topology; report the fixes.
- **S-c.** Call E1/E3 alignment-robust only if TM-align and sequence-only agree in direction and CI.
