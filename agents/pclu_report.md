# P-clu report: structure vs sequence inside AFDB clusters (v3 file 8)

Analyst run, 2026-10-30. Seed 20261030. Nothing committed.
Scripts: `scripts/pclu/build_pairs.py` (about 70 s) then `scripts/pclu/analyse_pclu.py` (about 9 min; the curve fitting dominates).
Tables: `results/pclu/*.tsv|json`. Figures: `figures/pclu/*.png` (300 dpi). Intermediate: `data/processed/pclu/pairs.parquet` (not for git).

## 0. Data and verified file layout

* File 8 has no header. Columns (verified against file 1 and file 2): `query` = cluster **representative**, `target` = cluster **member**, fident, alnlen, mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits, lddt, alntmscore.
  * `orientation_query_is_rep_frac` = 1.0 (every query equals the repId of its target in file 1). All 2,302,908 queries are in file 2. There are no duplicate pairs and no unmatched targets.
* 30,045,247 rows = 2,302,908 self-pairs (rep to rep) + **27,742,339 member-to-rep pairs** (all analysed; 2,302,908 clusters, nMem >= 2).
* Joins: file 1 gives the member taxId. The rep taxId is the file-1 row where member = rep. File 2 gives repLen, repPlddt, isDark, nMem and the LCA taxId. NCBI taxonomy comes from my own `nodes.dmp`/`merged.dmp`/`rankedlineage.dmp` parser (merged taxIds resolved; every taxId resolved). The "superkingdom" is the node with rank "domain" (taxIds 2, 2157, 2759, 10239).
* Sanity checks (`sanity_checks.json`):
  * All self-pairs have fident = 1.000 and TM = 1.000. LDDT >= 0.999 in 99.93% (minimum 0.871).
  * 1,723 self-pairs (0.07%) have an alignment length different from repLen, so some self-alignments are partial. TM is still exactly 1 for them (see the TM caveat below).
  * fident <= 1 - mismatch/alnlen holds for every pair.
  * Only 5,817 non-self pairs have fident = 1.

## 1. Joint distribution (figs `pclu_hexbin_fident_vs_tm_lddt.png`, `binned_curve_overall.tsv`)

* fident quantiles (1/5/25/50/75/95/99%): 0.076/0.102/0.153/0.209/0.295/0.482/0.687.
* 23.5% of pairs have fident < 0.15 and 4.4% have fident >= 0.5. Only 0.9% have fident >= 0.7 and 0.19% have fident >= 0.9.
* TM quantiles: 5% 0.22, median 0.63, 95% 0.92. 31% of pairs have TM < 0.5 and 10% have TM < 0.3.
* Pooled over all pairs, the rank correlations are weak to moderate. Spearman of fident with TM is 0.59, and of fident with LDDT is 0.82. Spearman of TM with LDDT is 0.72.
* **Shape of the median curve** (0.02 bins; IQR is wide throughout):
  * Median TM rises steeply from 0.39 at fident 0.10 to 0.66 at 0.20 to 0.22 and 0.79 at 0.35. **It then plateaus at 0.83 (IQR about 0.62 to 0.94) from fident 0.4 to 0.8**, and only rises above 0.85 at fident > 0.9. It is a hockey stick, not a line.
  * Median LDDT behaves differently. It keeps rising concavely: 0.52 at 0.10, 0.68 at 0.25, 0.83 at 0.4, 0.88 at 0.7 and 0.94 at 0.99. It has a much narrower IQR.
  * The plateau of TM at high identity is not a structural-clock effect. It is driven by the pLDDT and length of the representative (section 4), and by pairs with TM < 0.6 even at fident >= 0.8 (21% of those pairs).
* **Upturn at fident < 0.07 is an artefact.** There the median TM goes from 0.28 up to 0.45, and the median alignment is short (query coverage 0.06 to 0.35, median E-value 3 to 900). These are junk or very short alignments for which TM is high because it is computed on the aligned region only.
  * Evidence for the normalisation: 1,704 self-pairs are aligned over < 100% of the rep but have TM = 1.000. The file does not document the TM normalisation. This is inferred, not verified.
  * Do not interpret anything below fident of about 0.10. All fits start at fident >= 0.10.
* Discrete-fraction artefacts: the LDDT hexbin has vertical stripes at fident = 1/3, 1/2 and 2/3 (very short alignments), and there is a bump in counts at exactly fident = 0.50 (120,668 pairs in the 0.50 bin vs 86,434 and 85,256 in its neighbours). Both come from short alignments, not from a clustering threshold.

## 2. Curve-family fits (`curve_fits_*.tsv`; fig `pclu_curve_fits.png`)

Responses are 1-TM and 1-LDDT.
x axes: p = 1-fident; Poisson d = -ln(fident); Kimura d = -ln(1-p-0.2p^2), which is **undefined for fident < 0.146 (22% of pairs)**, so Kimura fits use fident >= 0.16 only.
Families:
* linear (with intercept), linear through the origin,
* saturating 1-exp (2 and 3 parameters),
* power,
* Chothia-Lesk exponential,
* hill sigmoid,
* **broken stick** (plateau + linear rise: y0 + s*max(0, x-xc)). I added this after seeing the plateau.

Fit data: (a) the 0.02-fident bin medians weighted by n (AIC with weights normalised to mean 1, so the effective N is the number of bins); (b) 1,000,000 random rows (OLS AIC; the dAIC values are huge and only ordinal). Windows: fident >= 0.10, >= 0.16, 0.10 to 0.50, >= 0.30.

Main findings on bin medians, fident >= 0.10:
* **No single smooth family fits.** The best by AIC is the broken stick for 1-TM on both p and Poisson d. For 1-TM versus Poisson d (dAIC vs the best):
  * linear: +115
  * Chothia-Lesk exponential: +18 (the best smooth family)
  * power and hill: +83 to +85
  * Plateau 0.18 followed by a rise of 0.39 per Poisson unit beyond d = 1.14 (fident about 0.32).
  * RMSE of the broken stick is 0.011, against 0.040 for linear.
* **Saturating (concave) models never win.** The 1-exp families collapse towards a line (A about 4000, b about 5e-5) and have dAIC > 90 for 1-TM. Within fident 0.10 to 1.0, structural divergence does not saturate on Poisson d. It accelerates with d.
  * The "saturation" seen at high identity is a floor in 1-TM, not a ceiling. (The true ceiling, TM about 0.17 for random pairs, is not reached because membership is selected.)
* For 1-LDDT versus Poisson d, the broken stick (0.14 + 0.31*(d-0.83)+) is best. Power is next (dAIC +9), then Chothia-Lesk exponential (+10) and hill (+11). Linear is +55.
* On **Kimura distance (fident >= 0.16)**, the plateau-plus-line fits 1-TM best (y0 = 0.17, slope 0.11 per Kimura unit beyond d = 1.04), and plain linear is within 2 AIC units of it (0.080 + 0.100 d, RMSE 0.015). 1-LDDT is best with an offset saturating curve (RMSE 0.008).
  * So Kimura correction makes the mid-range look almost linear. But Kimura diverges at fident = 0.146, which stretches the low-identity points. This is a property of the correction, not independent evidence for a clock.
* Random 1M-row fits (rows, not bins) agree on the ordering for the fident >= 0.10 window. The Chothia-Lesk exponential is the best of the smooth families for p and Poisson d, and linear is best for Kimura/1-TM. The row-level residual is large (RMSE about 0.17 for 1-TM), because the scatter around the median at fixed fident is very wide.
* Sensitivity (`curve_fits_cluster_balanced_max5.tsv`, `..._nominal_pass.tsv`, `..._repPlddt_ge90.tsv`; fig `pclu_sensitivity_curves.png`): the same conclusions hold in the cluster-balanced set (at most 5 random pairs per cluster, 6.5M pairs) and in the subset passing the nominal clustering criteria (qcov >= 0.9 and E < 0.01, 12.1M pairs). The broken stick wins for 1-TM in both.
  * In the rep pLDDT >= 90 subset (6.1M pairs) the plateau disappears: median TM rises from 0.40 to 0.99 in a smooth convex curve, and a power law in Poisson d (exponent 2.1) beats the broken stick.
  * Cluster-balanced medians are about 0.05 lower in TM at fident 0.4 to 0.8 (plateau 0.79 vs 0.83), because big clusters have better-ordered proteins.

## 3. Selection truncation (`truncation_summary.json`, `truncation_by_fident_bin.tsv`; fig `pclu_selection_truncation.png`)

* **File 8 alignments do not honour the nominal clustering filters** (90% bidirectional coverage, E < 0.01). They are re-computed alignments of member to representative, and many cluster members were joined via the cascade or set-cover rather than by direct hits to the representative.
  * Query (rep) coverage is below 0.9 for 52% of pairs and below 0.5 for 12.7% (rep length from file 2; the member length is not in the files, so target coverage cannot be computed).
  * E-value > 0.01 for 21.4%, > 1 for 9.8%, > 10 for 6.8%. Only **43.6%** of pairs satisfy both qcov >= 0.9 and E < 0.01.
  * Only 7.1% of pairs have E between 1e-3 and 1e-2, i.e. there is no pile-up at the E threshold.
* Consequence for the low-identity tail:
  * The tail is not truncated at a hard TM or identity floor. About 10% of pairs have TM < 0.3, and the fraction of pairs with E > 0.01 reaches 83% at fident 0.06 to 0.08. The tail is a mixture of genuine remote homologs, partial alignments and noise.
  * Membership is still a truncation of the pair population. All pairs belong to the same cluster, so no pairs with no detectable similarity are present (e.g. random pairs, TM about 0.17).
* **Truncation at high identity (by design):** members are AFDB50 representatives, so pairs with fident >= 0.5 are only 4.4% of the data (the 0.40 to 0.60 counts decay smoothly, apart from the artefact at exactly 0.50).
  * There are 627k pairs with fident in [0.45, 0.50); the density has no cliff at 0.5 because the Foldseek fident differs from the MMseqs2 identity used by AFDB50.
* The median TM, LDDT and the 5th TM percentile vary with fident bin in the same way across the filters. The IQR of TM never narrows below about 0.3 at fident >= 0.4.

## 4. Stratification (`stratum_summary.tsv`, `stratum_curve_fits.tsv`, `binned_curve_by_*.tsv`; figs `pclu_stratified_median_curves*.png`)

Medians of TM at fident >= 0.6 / at 0.2 to 0.3 (the second removes most of the composition effect on fident):
* **Pair superkingdom** (n = 23.2M same, 3.9M different, 0.66M unclassified or virus):
  * Different-superkingdom pairs have a *higher* TM at fixed fident: 0.96 vs 0.83 at fident >= 0.6, and 0.75 vs 0.70 at 0.2 to 0.3.
  * This is a selection effect. Cross-superkingdom clusters keep only conserved compact domains (sharing a structure across superkingdoms with high identity needs strong constraint).
  * Within same-superkingdom pairs, Euk-Euk is the lowest. Its TM plateau at fident >= 0.6 is 0.76, vs 0.89 for Bac-Bac and 0.91 for Arc-Arc. Many are multidomain or disordered proteins.
* **Cluster LCA rank** (age class; the LCA is the lowest standard-rank ancestor of the cluster's member taxa): plateau TM at fident >= 0.6 rises with age:
  * 0.75 (species/genus), 0.81 (family to class), 0.81 (phylum/kingdom), 0.84 (one superkingdom), 0.90 (cross-superkingdom).
  * 64% of pairs are in clusters spanning superkingdoms (LCA root/cellular). Median fident is lower there (0.19 vs 0.42 for species/genus clusters).
  * Younger clusters thus look *less* structurally conserved at equal identity, again because they are less constrained or lower in pLDDT, not because divergence accelerates with age.
  * LDDT curves are almost identical between age classes.
* **Dark vs annotated** (2.26M vs 25.5M pairs): dark clusters have lower TM at fixed fident (0.78 vs 0.85 at >= 0.6), and a higher fraction of E > 0.01 (34% vs 20%). LDDT is nearly identical.
* **Representative pLDDT** is the largest effect. TM at fident >= 0.6 is 0.58 for pLDDT < 60, 0.68 for 60 to 70, 0.79 for 70 to 80, 0.91 for 80 to 90 and 0.97 for >= 90. The high-identity plateau is mostly a pLDDT effect.
* **Length** (rep length): TM at fident >= 0.6 falls from 0.90 (< 100 aa) to 0.67 (>= 800 aa). The longest proteins also show the lowest TM at fident 0.2 to 0.3 (0.55 vs 0.77). LDDT is nearly length independent.
* **Filters and cluster size**: higher query coverage and lower E-value bins have higher TM at fixed fident. The E-value > 1 bins are a distinct low-TM population (median TM 0.29). Larger clusters have higher TM (0.94 vs 0.79 at fident >= 0.6 for >= 1000 vs < 5 members).
* The per-stratum AIC results (`stratum_curve_fits.tsv`: linear, saturating, power, Chothia-Lesk exponential on Poisson d, fident >= 0.10) show convex, not saturating, shapes in nearly every stratum. Linear is never the best family for 1-TM.

## 5. Caveats

1. The **TM normalisation** in `alntmscore` is undocumented here (inferred as aligned-region based; see section 1). It is probably not the usual TM normalised by the shorter chain. Do not compare absolute TM values with the orthology set unless the same Foldseek `alntmscore` is recomputed there. Use LDDT, which is more stable, as the primary comparison metric.
2. `fident` is the identity inside a Foldseek (3Di + amino acid) alignment, not from a sequence-only alignment (shared alignment error with the structural metrics; Illergard-style confound).
3. The pairs are a **star topology** (members vs one representative) and are not independent. Big clusters dominate (the top 1,000+ member clusters have 4.6M pairs). Cluster-balanced sensitivity runs are provided, but no hierarchical uncertainty or bootstrap over clusters was computed. The reported dAIC values are descriptive.
4. Only 4.4% of pairs have fident >= 0.5. High-identity behaviour is poorly sampled and mostly reflects pLDDT and length.
5. **No time axis here.** The LCA rank is a coarse proxy for age, and the pair shared-lineage rank (`pair_mrca`) uses rankedlineage names (many ranks are missing for bacteria and unclassified taxa, so "domain" means "share only the domain label").
6. 0.6M pairs (2.4%) involve unclassified or virus taxa and are shown separately.
7. This set is structure-selected. Its curves must be contrasted with orthology-selected pairs, not interpreted alone as evolutionary linearity.
8. Strata overlap (pLDDT, length, superkingdom, dark are correlated). The stratified tables are marginal, not mutually adjusted.
9. The 1M-row AIC differences are very large because N is huge. Use them only for ranking.
