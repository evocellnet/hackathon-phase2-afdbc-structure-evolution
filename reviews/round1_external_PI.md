# Review round 1: external PI (scientific strategy)

## 1. Overall verdict

**Major revision.** The data plumbing (v3 clusters, FoldComp v4, TED, a validated Ancestra, TimeTree) and D1 are sound. The weakness is in the logic. Goal 1 has six hypotheses that are not independent, and the neutral-clock question has no estimand or decision rule.

The problem is composition. If structural divergence S = f(AA distance) and AA distance = g(time), then H1b mostly follows from H1a and the known g. **The new information is in the residual**: does structure carry rate variation by lineage or family beyond what sequence explains, and is f the same across lineages? This decomposition is the single most important result available: is structure an independent clock, or only a compressed sequence clock? The design does not put it first.

Goal 2 has five items. H2c and H1f largely repeat published work. The ~300-proteome `accurate` Ancestra run is too large for the time before the node move.

## 2. Change proposals

**1. [MUST] One estimation engine: branch lengths per HOG on the fixed species topology.** For each HOG with ≥ 8 single-copy species, prune the species tree to those species and estimate branch lengths on it:
- AA: LG+Γ4 on a FAMSA MSA.
- 3Di: Q.3Di.AF+Γ4 on a FoldMason MSA (IQ-TREE `-te`).
- TM and LDDT: OLS/FastME on the fixed topology from −ln(TM) and corrected-LDDT matrices.

Each branch is then an independent contrast with a TimeTree duration. This replaces the pseudo-replicated (HOG, node) medians and the low-power triplet tests. Keep pairwise metrics only for Illergård-style curves and P-clu. Cost: a few CPU-hours for about 5k HOGs.

**2. [MUST] Three estimands with decision rules written in advance.**
- **E1, response shape f:** per-HOG c in S = a·d_AA^c, plus d*. Call it linear if the median c has a 95% CI inside [0.9, 1.1] and saturating models improve out-of-range prediction by ≤ 10%.
- **E2, rate ratio k:** tree length of 3Di / tree length of AA on a shared FoldMason MSA, by CATH class.
- **E3, structure-specific rate variation:** a mixed model log b_struct ~ log b_AA + log duration + (1|HOG) + (1|lineage). Report the leftover variance components, and test whether the log b_AA slope equals 1.

Report H1b as Δc, the curvature of structure against time relative to AA against time. With millions of points AIC/BIC always reject linearity, so decide on effect sizes with HOG-bootstrap CIs.

**3. [MUST] Make "neutral structural clock" a falsifiable null.** The null (passive response, after Echave): structural branch length ∝ AA branch length, with a ratio that varies by family but is the same across lineages. Goal 2 gives the **positive control**: paralogs and the more-diverged ortholog (MDO) should depart from the ortholog null. If they do not, the test has no power and no neutrality claim can be made.

**4. [MUST] Measure the noise floor before comparing clock-likeness.** 3Di, TM and LDDT carry noise from the prediction and from coordinates. A noisier metric looks "less clock-like" even when its true rate is constant.
- Estimate the noise from identical sequences with independent structures: several PDB entries or chains, or AF v4 vs v6 models if they were re-predicted.
- Compare the dispersion of each metric against a simulated strict clock with that noise added.

**5. [MUST] Keep the AA axis independent of the structural alignment.** For E1, compute d_AA from sequence-only MSAs and S from structural alignment. Using one Foldseek 3Di+AA alignment for both creates correlated error (LITERATURE §1a); a 5% TM-align check is not enough.

**6. [MUST] Shrink and re-target the panel; launch Ancestra on node 1.**
- Run eukaryotes only: about 150–180 species, `fast` unless a timing test shows `accurate` finishes in ≤ 6 h.
- Choose species for roughly log-uniform MRCA ages from 10 Mya to 1.8 Ga, not vertebrate-heavy.
- Require ≥ 90% of each proteome to have an AFDB v4 model with an identical sequence. Missing genes cause hidden paralogy to be called orthology.
- Avoid polyploid proteomes above 60k proteins.

Why: 200 eukaryotes plus 100 prokaryotes is about 4.4M proteins, the scale of the authors' 10.6 h run. Prokaryotes add no time information, and cross-domain orthology is dominated by HGT. Handle them through Pfam/TED/Foldseek sets for the case study.

**7. [MUST] Treat the time axis honestly.** TimeTree ages come largely from molecular clocks, so AA against time is partly circular.
- Propagate the age CIs with errors-in-variables fits or age resampling.
- Run primary fits only on nodes supported by ≥ 5 studies.
- Report nodes ≥ 1 Ga separately: their ±20% uncertainty can create or hide saturation.

**8. [MUST] Reduce Goal 2 to three items.**
- **(i) Merge H2b and H2c:** "Does duplication break the response?", measured as residuals from the E1/E3 ortholog null. This is also the positive control for proposal 3.
  - Define the MDO from AA branch lengths.
  - Add a symmetric control: define the MDO by structure and test AA.
  - H2c alone repeats PMID 40813247.
- **(ii) Merge H2a and H2d into a cluster audit:** HOG/cluster splits explained by TED architecture, length and pLDDT. Add one sharp question: is the cluster LCA older or younger than the HOG root, and how many of the 3.9% species-specific clusters are old HOGs?
- **(iii) One case study, time-boxed to 4 h:** the cytoskeleton (tubulin/FtsZ/CetZ with Asgard), where Ancestra is already tested. Kinases only if time remains.

**9. [SHOULD] Test the lineage effects against Daumiller et al. 2026.** Include ≥ 5 birds and ≥ 5 Hymenoptera, with their sister clades. Test whether their reported structural acceleration survives conditioning on AA branch length.

**10. [SHOULD] Re-scope H1e.**
- A few hundred close pairs rich in PDB templates only bound the bias for proteins with templates. They cannot show "no bias" in general.
- Add a check that E1 and E2 do not change across pLDDT strata, and use v6 MSA depth if it is cheap.
- Could: run ESMFold on about 5k proteins as an MSA-free comparison for f.

**11. [SHOULD] Check the k convention.** Normalise LG and Q.3Di.AF to one expected substitution per site. Report k under JC-20 as a sensitivity check. Do not estimate our own 3Di GTR unless Q.3Di.AF cannot be obtained.

**12. [SHOULD] Report representativeness.** Single-copy HOGs are biased toward slow-evolving housekeeping genes. Compare their GO slim profile with that of all HOGs.

**13. [COULD] Demote H1f.** Report it as a by-product of the proposal-1 trees; FoldTree has already benchmarked it.

**14. [COULD] Run P-clu (file 8) now,** during the Ancestra run, as the "selection creates saturation" figure.

## 3. Priority order

1. **Node 1, now:** set the panel and launch Ancestra (proposal 6). In parallel:
   - TimeTree node ages;
   - the Q.3Di.AF matrix;
   - AliSim validation of the branch-length engine;
   - the noise floor (proposal 4);
   - P-clu.
2. Fetch FoldComp structures for the HOGs; build FAMSA and FoldMason MSAs.
3. **Headline:** E1–E3 and the neutral-null decomposition (proposals 1–3, 5, 7). This must not slip.
4. Goal 2(i): duplication against the null, which is the positive control.
5. Goal 2(ii): the cluster audit.
6. The Daumiller test (proposal 9) and the H1e bounds (proposal 10).
7. The cytoskeleton case study, time-boxed.
8. H1f and ESMFold, only if time remains.
