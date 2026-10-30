# Round 2 review: external PI (DESIGN.md v2)

**Verdict: accept after targeted revisions.** v2 addresses most round-1 musts. Three problems remain:
- E1/E2 sit on branch lengths, where saturation is invisible.
- The pair budget is unbounded.
- NNLS zero branches break the log-ratio endpoints.

## 1. Round-1 must items

| Item | v2 | Gap |
|---|---|---|
| PI1 fixed-tree engine | yes | NNLS on pairwise ML, fine |
| PI2 estimands | yes | see R1 |
| PI3 N0 + positive control | partial | E6 does not show E5 power |
| PI4 noise floor | partial | identical-seq AF pairs share an MSA, so floor ≈0 |
| PI5/MO-M1/AI3 AA axis independent | partial | TM-align primary is a defensible deviation; needs R5 rule |
| PI6 panel/Ancestra | partial | no depth balance, `fast` switch or polyploid rule |
| PI7 time axis | partial | ≥5-study flag missing |
| PI8 Goal 2 → 3 | yes | |
| MO-M2/M3/M4/M5/M6 | yes | |
| MO-M7 LDO | partial | no >95% identity exclusion; outgroup rule undefined (D2) |
| ST-S1 hidden selection | partial | coverage rule drops pairs in primary |
| ST-S2 errors-in-variables | yes | per-branch SE source unstated |
| ST-S3 composition | partial | pooling of c across HOGs unspecified |
| ST-S4/S5/S6/S7/S11 | yes | Holm includes non-test E3 |
| AC-C1–C5 | yes | |
| AC-C8 node-1 order | partial | euk20 pilot and `fast` rule absent |
| AI1/AI2/AI5/AI6/AI8 | yes | fix LESSONS #2; add DEVIATIONS.md |
| AI4 parser contract tests | **no** | |
| AI7 resume | partial | Ancestra S2 is not resumable within the stage, though v2 says it is |

## 2. Whole-design assessment

**The engine fits some endpoints, not all.**
- It fits E3 (ratio of tree lengths), E4 and E5.
- It does **not** fit E1/E2. Saturation is a property of cumulative path distance. An additive fit turns it into misfit and shrunken deep branches, so branch-wise b_3Di vs b_AA measures local proportionality, not response shape.

**E4 largely restates E1.** If 3Di is concave in d_AA, k rises with age automatically. The new question is whether age predicts structure at fixed AA divergence.

**Zero branches.** NNLS gives exact zeros on short branches, so log(b_s/b_AA) is undefined, and dropping them selects on noise.

**Pair budget.**
- Likely >10k eligible HOGs × 780 pairs gives ≥8M pairs.
- At the measured `tmalign` rate of 265 pairs/s that is about 9 h, plus the E6 pairs.
- This is not budgeted.

**Over-engineered.**
- Tree-refitting two-way bootstrap × 200 age draws × SIMEX.
- H1f.
- A prokaryote panel used only for 2-iii.

**AI-execution risk.** There are many bespoke statistical parts: NNLS with missing cells, SIMEX on ratios, crossed random effects, TOST. Component tests will not catch integration errors. An end-to-end simulation (AliSim tree → pairs → engine → E1–E5 with known answers) is needed.

## 3. Remaining changes

### MUST

- **R1. E1/E2 on within-HOG pairwise distances.**
  - Fit log D = log a_h + c_h log d_AA over pairs with d_AA ≤ 2.
  - Random intercept and slope by HOG; report mean c.
  - Two-way bootstrap plus SIMEX.
  - Include only HOGs spanning ≥3 depth bins.
  - The branch engine serves E3–E5.
- **R2. Pair budget.**
  - Cap at 30 species per HOG and about 7k depth-stratified HOGs (seed 20261030).
  - Primary plus E6 pairs ≤3.5M.
  - Re-time `tmalign` on the 20% set before freezing.
- **R3. Zero branches and SEs.**
  - E4/E5 use either ratio-of-sums within age-bin strata, or branches with b_AA ≥ a threshold fixed on the 20% set (with sensitivity).
  - Branch SEs come from the NNLS sandwich with Fisher pair variances, validated in test 2.
- **R4. E4 confirmatory estimand.** Slope on log age with log b_AA as a covariate. The unconditional slope is descriptive.
- **R5. Pair failures and robustness.**
  - Never drop single pairs. A HOG enters E1–E5 only if ≥90% of its pairs pass, and exclusions are reported against depth.
  - Call a result "alignment-robust" only if sequence-only alignment gives the same E1 class and k within 25%.
- **R6. Ancestra.**
  - Use `fast` if the euk20 projection exceeds 5 h.
  - If S2 is unfinished by 19:30, kill it and relaunch first on node 2.
  - Pilot the metric pipeline on euk20 meanwhile.

### SHOULD

- **S1.** E5 power by simulation: inject ×1.5 and ×2 structural acceleration on 5% of branches. N0 statements are conditional on detecting it.
- **S2.** Noise floor from AF v4 vs v6 (if re-predicted) or multi-chain PDB entries. Treat identical-sequence pairs as a lower bound.
- **S3.** Apply Holm only to tests (E1/E2 TOST, E4, E5, E6), not E3. Set the E5 equivalence margin on the 20% set.
- **S4. Simplify.**
  - Bootstrap fitted branch tables without refitting.
  - One age draw per replicate.
  - Drop H1f.
  - 2-iii without a prokaryote panel (members from AFDB clusters or Pfam).
- **S5. Remaining fixes.**
  - Parser tests (AI4).
  - LDO: exclude pairs with >95% identity; define the outgroup as single-copy species outside the smallest clade containing both paralogs.
  - ≥5 birds and ≥5 Hymenoptera.
  - ≥5-study flag.
