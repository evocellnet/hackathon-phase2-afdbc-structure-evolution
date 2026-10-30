# Round 1 review: molecular evolution

Reviewer: molecular-evolution agent (Opus 5.5), 2026-10-30.

## 1. Verdict

**Sound overall, but revise before running.** H1a/H1b separation and avoiding Ancestra duplication timing are right. Five problems:
- the AA and structure axes share one alignment that is optimised for structure;
- the dispersion and relative-rate tests are not valid for continuous, bounded, non-additive metrics;
- "neutral" is conflated with "clock-like";
- the time axis uses an unresolved topology and ignores the uncertainty of node ages;
- the ortholog/paralog and LDO contrasts are exposed to known artefacts.

## 2. Change proposals

**M1 (must). Measure the AA axis independently of structure.** 3Di in the Foldseek score biases 3Di distance down (inflating k, H1c) and makes AA distance partly structural (circular for H1a); LITERATURE §1a flags it.
- Primary AA and 3Di distances: compute both on a sequence-only pairwise alignment (FAMSA, or the DIAMOND alignment), over the same columns.
- 3D metrics: use a geometry-only alignment (TM-align/US-align).
- Foldseek 3Di+AA: keep only as a sensitivity analysis.

**M2 (must). Floor-correct the bounded metrics and test additivity.**
- 1−TM and −ln(TM) do not correct for saturation. Use d = −ln((S−S₀)/(1−S₀)), with S₀ estimated from the random-pair negative control. For TM, estimate S₀ per length bin, because the random-pair floor depends on length.
- Pre-state: the clock claim is linearity in corrected distance.
- Clock tests assume additivity: per HOG and metric, fit branch lengths on the fixed species topology by least squares and report the residual percent SD. Exclude metrics that fail from the clock tests.

**M3 (must). Make the dispersion and relative-rate tests valid.**
- **Dispersion index.** R = var/mean is defined for counts and depends on units, so comparing R for TM or LDDT with R for AA is meaningless.
  - For AA, use R on expected substitution counts (Cutler 2000; Bedford et al. 2008).
  - For all metrics, use CV² of lineage rates minus sampling variance (per-residue bootstrap).
- **Relative-rate test.** Swapping ingroup labels only flips the sign of d(A,O)−d(B,O), and that symmetry is already the null hypothesis, so the planned permutation null is not valid. Use the distance tests of Wu & Li or Tajima, or the tests of Takezaki et al. 1995, with bootstrap variances.
- **Power.** Metrics with fewer effective characters have less power, so "fewer violations" does not mean "more clock-like". Compare rate-ratio effect sizes with CIs, not rejection fractions.

**M4 (must). Correct the logic used to interpret neutrality.** Clock-like does not mean neutral: strong purifying constraint lowers the rate and can still be perfectly clock-like.
- **Null model.** Use passive response (Echave): structural change is a by-product of AA substitutions, scaled by family. It predicts:
  - (a) a constant branch-wise structure/AA ratio within a HOG;
  - (b) structural lineage effects that correlate about 1 with AA lineage effects, including generation-time and population-size (Ne) shifts;
  - (c) no structural overdispersion beyond what AA overdispersion implies.
- Constraint = between-family ratio differences; selection/function change = episodic within-family branch outliers co-occurring with architecture/GO change (Pascual-García 2010, 2019).
- Write (a)–(c) into H1d, replacing "dispersion comparable to AA".

**M5 (must). Use a resolved, dated topology and propagate time uncertainty.**
- **Topology.** Deep NCBI nodes are polytomies, which pool MRCAs that differ by hundreds of Myr, and the eukaryotic root is contested. Use the TimeTree 5 panel subtree, checked against Burki et al. 2020, as both the Ancestra guide tree and the time axis.
- **Node ages.** Draw node ages from their CIs in each bootstrap replicate.
- **Alternative timescale.** The H1b verdict hinges on a few high-leverage nodes older than 1 Gyr. Re-run H1b with the Strassert et al. 2021 or Betts et al. 2018 timescale.

**M6 (must). Ortholog vs paralog (H2b).**
- **Matched comparison.** Out-paralog vs ortholog pairs from the same species pair and HOG, sharing one protein (A1–B1 vs A1–B2). Within-species paralogs differ from cross-species orthologs in annotation source and gene-model quality, not only in function (Altenhoff 2012; Weisman 2022).
- **Drop the in- vs out-paralog axis.** The split depends on duplication timing, which Ancestra gets right only about 36% of the time.
- **Gene-model filters.** Filter on length ratio, fragment flag, UniProt evidence (PE) level and pLDDT difference; a missing exon lowers TM.
- **Sensitivity for hidden paralogy.** Restrict to single-copy HOGs, because pseudo-orthologs inflate deep ortholog divergence (Natsidis 2021).
- **Credit.** Cite Peterson et al. 2009 as the direct precedent: orthologs are more structurally similar than paralogs at 30–70% identity.

**M7 (must). LDO/MDO bias (H2c).** If the MDO is chosen by AA distance to an outgroup O, and structural excess is then measured against the same O on a shared alignment, alignment errors inflate both and create a spurious excess.
- Choose the MDO with outgroup O₁ and measure the excess against an independent O₂ (cross-fitting), or use gene-tree branch lengths.
- Require O to lie outside the duplication on a reconciled gene tree (GeneRax subset).
- Exclude pairs above 95% identity (gene conversion).
- Positive control: reproduce Warwick Vesztrocy et al. 2025 matched-identity LDDT result.

**S1 (should). Make H1c truly alphabet-matched.**
- 20 = 20 states does not mean equal information. 3Di frequencies are skewed, so use the effective number of states K_eff = 1/Σπ² for both alphabets in the corrections.
- 3Di letters depend on spatial neighbours, so sites are not independent.
- 3Di has a noise floor. Estimate it from identical-sequence AF pairs and multiple PDB chains of the same protein, and fit k with an intercept.
- Estimate k at low divergence, and from tree-length ratios on a fixed topology (Q.3Di.AF vs LG+Γ).
- To bridge to Illergård, also compute their DSSP 3-state and accessibility 2-state descriptors against AA alphabets reduced to the same number of states.

**S2 (should). Estimate lineage effects branch-wise, not as per-species residuals.** Use the M2 branch fits and correlate AA and structural branch rates. Add a per-species mean-pLDDT covariate: AF confidence varies by lineage and would otherwise look like a structural lineage effect.

**S3 (should). Report aligned fraction and high-pLDDT coverage against time.** Local alignment and the pLDDT mask drop divergent regions, which mimics saturation. Add a variant that counts gaps as divergence (TM normalised by full length).

**S4 (should). Use errors-in-variables (Deming/SIMEX) for the H2b residual curves.** Restrict comparisons to AA-distance bins where orthologs and paralogs overlap; regression dilution plus curvature can otherwise create spurious residuals.

**C1 (could).** Structural expression–rate anticorrelation (Drummond 2005).

**C2 (could).** For vertebrates, compute structural change per synonymous substitution (dS from Ensembl CDS), a "structural ω" that removes generation-time effects (Ohta 1993).

## 3. Missing citations (PMIDs verified via PubMed)

Peterson 2009 Protein Sci (PMID 19472362, 10.1002/pro.143; H2b precedent) · Bedford, Wapinski & Hartl 2008 Genetics (18505862, 10.1534/genetics.108.089185; overdispersion vs Ne) · Cutler 2000 Genetics (10757779; R(t)) · Takezaki, Rzhetsky & Nei 1995 MBE (7476128; distance clock tests) · Wu & Li 1985 PNAS (3856856) and Tajima 1993 Genetics (8244016; relative-rate tests) · Natsidis 2021 iScience (33659875, 10.1016/j.isci.2021.102110; orthology errors) · Weisman, Murray & Eddy 2022 Curr Biol (35588743, 10.1016/j.cub.2022.04.085; annotation heterogeneity) · Studer & Robinson-Rechavi 2009 Trends Genet (19368988) and Rogozin 2014 GBE (24610837; ortholog conjecture) · Ohta 1993 PNAS (8248159; generation time) · Burki 2020 TREE (31606140; eukaryote tree/root) · Strassert 2021 Nat Commun (33767194) and Betts 2018 Nat Ecol Evol (30127539; deep timescales).
