<!--
SECOND HALF OF THE TECHNICAL REPORT (Abstract, §4–7, Data and code availability).
Figure files to copy into "technical report/figures/":
  Main (from figures/main/):  figA_shape.png, figB_selection_bias.png, figC_time.png, figD_ratio_k.png, figE_terminal_internal.png
  Goal 2 / bias checks (chosen by the writer; copy and rename as follows):
    figures/h1e/x4_paired_scatter_tdiaa.png      -> figures/figF_af_vs_pdb.png
    figures/goal2/hog_cluster_mapping.png        -> figures/figG_hog_cluster_mapping.png
    figures/goal2/hog_split_model.png            -> figures/figH_hog_split_model.png
    figures/casestudy/method_comparison.png      -> figures/figI_cytoskeleton_methods.png
Numbering continues from Figures 1–3 of the first half: Fig A–E = Figures 4–8, Fig F–I = Figures 9–12.
Every number is taken from results/endpoints/*.json (confirmatory), results/posthoc/, reviews/results_review_statistics.md,
agents/*_report.md or figures/main/key_numbers.json. Items marked [CHECK] need a human look before release.
-->

## Abstract

Protein structure is said to be "more conserved than sequence" and to change "linearly", but these claims rest on small sets of crystal structures, with amino-acid substitutions used as the clock. We revisited them with AlphaFold models for 180 eukaryotic proteomes, ortholog families inferred from sequence alone (Ancestra), and dated species trees. Amino-acid (AA) change was compared with change in Foldseek's 20-letter structural alphabet (3Di) under matched maximum-likelihood models. Six tests were pre-registered, tuned on 20% of 7,000 families and confirmed once on the other 80%.

The pre-registered tests found four things:
- Structure diverges more slowly than linearly with sequence (curvature exponent 0.72, 98% CI 0.64–0.80).
- Sequence changes 3.3 times faster than 3Di (α-helical families 4.7×, β-rich 2.6×), at the low end of the classic "3–10×".
- The sequence-to-structure ratio rises with branch age.
- Apparent lineage-specific structural rates were classed as an artefact.

A post-hoc review, checked by a second agent, traced most of these results to one cause: each AlphaFold model carries protein-specific noise (about 0.01–0.017 3Di substitutions per site, larger at low confidence) that the pre-registered noise floor missed. Allowing for it gives curvature estimates of 0.73–0.98. On internal tree branches the ratio stays at about 4.6 from 30 million to over 1,000 million years, and the lineage effects disappear. Our best reading is that 3Di structure behaves like a rescaled sequence clock; linearity and neutrality are not established. AlphaFold models also make orthologs look 30–40% more similar than crystal structures do.

For Goal 2:
- AFDB structural clusters split and merge ortholog groups, mainly because of length differences.
- The faster-evolving copy after duplication shows extra structural change, replicating the direction of a published result.
- Cytoskeleton and kinase case studies found no robust structural bursts.

---

## 4 Results: Goal 1, is structure a molecular clock?

**Main point.** The pre-registered analysis says that structure (as 3Di) saturates relative to sequence, that it changes about 3.3 times more slowly than sequence, and that the ratio grows with age. All of these decisions were computed correctly and re-derived independently (42/42 checks, V11). A critical review after confirmation found that most of them come from one source: per-model AlphaFold noise that the fixed noise floor did not remove. Once this is allowed for, the picture is closer to "structure follows sequence at a roughly constant ratio of about 4.6". In each subsection we first report the pre-registered result unchanged, and then the post-hoc interpretation, labelled as such.

**What ran.** The single confirmation run (2026-10-31, 1.8 h) passed all ten integrity guards, including checksums of the plan, tuning file and split, and a clean code tree. 4,891 of the 5,600 confirmation families entered the branch endpoints. The other 709 were excluded:
- 510 failed the misfit screen;
- 183 had too many failed pairs;
- 195 had trees that could not be solved.

These counts overlap. The Holm-adjusted p-values for E1, E2, E4, E5 and E6 were all below 10⁻⁴. The pre-registered interpretation row is **R6, "artefact suspect"**. It is used only for joint reading, because the E5 power condition (≥ 80%) was not met (power 0.24).

| Endpoint | Pre-registered result (confirmation set) | Decision |
|---|---|---|
| E1 3Di vs AA shape | ĉ = 0.72 [98% CI 0.64, 0.80] | concave |
| E2 LDDT distance vs AA shape | ĉ = 0.70 [0.61, 0.79] | concave |
| E3 ratio k = ΣAA/Σ3Di | 3.32 [3.09, 3.57] | estimate |
| E4 ratio vs branch age | β = −0.150 [−0.200, −0.101] | AA:3Di ratio rises with age |
| E5 lineage effects | slope 0.57 [0.10, 1.05]; excess lineage SD 0.31 | artefact suspect |
| E6 duplication control | E6a 0.030 [0.026, 0.033]; E6b 0.017 [0.012, 0.022] | detected |

### 4.1 Shape: concave by the pre-registered test, but the test is biased by per-model noise

**Pre-registered result (E1, E2).** The test fits y = η + a_h·x^c, where x is the AA distance and y is the structural distance. On that test, 3Di distance grows more slowly than AA distance.
- **Estimate.** ĉ = 0.723, 98% CI [0.642, 0.803], from 591,410 pairs in 4,850 families (Figure 4a).
  - The naive fit gave 0.670.
  - SIMEX raised it to 0.717.
  - The bias correction (D31) added 0.006.
- **Alignment robustness.** The decision is "concave" under all three alignments: sequence-only 0.80 [0.70, 0.90] and TM-align 0.80 [0.72, 0.89] (Figure 4b).
- **Planned sensitivities.** All of them stay concave:
  - α × 1.5: 0.75;
  - α ÷ 1.5: 0.72;
  - hidden-paralogy pairs removed: 0.72;
  - single-copy families: 0.62;
  - AA distance ≤ 1: 0.66;
  - all aligned columns rather than the confident core: 0.54.
- **E2 (floor-corrected LDDT distance).** ĉ = 0.70 [0.61, 0.79]; sequence-only 0.72, TM-align 0.32. E2's secondary quantity d* is the AA distance at which raw 1 − LDDT departs from a straight line by 10%. It fell at the lowest grid value (0.05), so there is no initial linear range.
- **Power.** In simulations on the real tuning designs, a truly linear world was never called concave (0/50; mean ĉ 1.03). So the test is not prone to false "concave" calls under the noise model it was tuned for. Its power is low, however: a true c = 0.90 was detected only 55% of the time.

**Post-hoc finding: the fixed floor η makes the estimator lean concave.** The model assumes that two AlphaFold models of orthologs differ by one global amount when no evolution has happened. That amount is η = 0.052 3Di substitutions per site, the median for pairs of identical sequences. Non-identical orthologs carry more than that.
- **Size of the excess.** Pairwise, the excess is about 0.010 per protein. On terminal branches of the tree it is 0.013–0.017.
- **Two components.** About half of it appears with the very first substitutions; the rest builds up by AA distance ≈ 0.1.
- **Link to model quality.** It tracks model confidence. Per-protein effects fall from 0.034 to 0.020 across pLDDT quintiles, mirroring the identical-sequence floor (0.031 → 0.020).

An extra offset at the low end of the curve with a fixed intercept is exactly what bends a power-law fit downwards. We checked this in two ways.
- **Real data.** Letting each family have its own baseline (η_h) raises the naive ĉ:
  - 0.67 → 0.90 for AA distance 0.05–2;
  - 0.94 for 0.05–1;
  - 0.98 for 0.3–2.
- **Synthetic data.** On the real design with realistic per-protein offsets, the fixed-η estimator turned a true c = 1 into 0.75–0.78. The family-baseline estimator stayed close to the truth (0.92–0.98).

The second agent also found counter-evidence. A model with one effect per protein, which fits the data far better, puts ĉ back at 0.73–0.75. That model is itself biased low in synthetic worlds (a true c = 1 is returned as 0.75–0.94). No single synthetic world reproduces all three estimators. The pairwise data therefore bracket c at roughly **0.75–1.0**; they do not pin it.

**What we conclude.** *Shown:* the pre-registered estimate ĉ = 0.72 is not robust to how the zero-divergence baseline is specified. Post-hoc estimators give 0.73–0.98. *Not shown:* that structure saturates. Two observations argue against saturation:
- Saturation should get worse with distance. Instead, the curve is *more* concave at short distances (0.66 for AA distance ≤ 1, against 0.72 overall).
- The branch-level evidence in §4.3 points to near-proportionality.

The step at small distances could be AlphaFold noise for non-identical sequences, or a real structural component that saturates quickly. The data cannot tell these apart.

![Figure 4](figures/figA_shape.png)

**Figure 4. In the pre-registered fit, structural divergence rises more slowly than sequence divergence, but the curve is close to a straight line plus an offset.** Data are 1:1 ortholog pairs from all 7,000 families (tuning and confirmation), with the E1/E2 filters applied.
- **(a)** 3Di ML distance against AA ML distance. Points are medians in 0.1-wide bins, and the band is the interquartile range (IQR). The black curve uses the confirmed exponent ĉ = 0.72. The grey dashed line is the linear reference (c = 1). The orange dotted line is the identical-sequence floor η = 0.052.
- **(b)** The same data under the three alignments, each with its own fitted curve.
- **(c, d)** As (a, b) for the floor-corrected LDDT distance (ĉ = 0.70).

Under TM-align the LDDT medians bend upward beyond AA distance ≈ 1.3, so its low ĉ (0.32) is not a robust shape. The curves are visual overlays fitted against the raw pair distance, not refits. Most of the visible concavity sits below AA distance 0.25, where per-model noise dominates (§4.1). Sample: 741,062 pairs in 6,067 families (3Di) and 744,561 pairs in 6,114 families (LDDT).

### 4.2 Conservation ratio: sequence changes about 3.3× (pre-registered) to 4.6× (internal branches, post hoc) faster than 3Di

**Pre-registered result (E3).** Summed over all branches of all family trees, AA branch length divided by 3Di branch length is **k = 3.32 [98% CI 3.09, 3.57]**. It is similar under the other alignments: sequence-only 3.06, TM-align 3.08. The CI overlaps Illergård et al.'s "3–10×" (Figure 7). Splitting the families shows four patterns:
- **Fold class.** α-helical families have k = 4.68 and β-rich families k = 2.55, using CATH classes from TED (n = 1,336 and 696 families). The difference in log k is 0.61 [95% CI 0.53, 0.68], unadjusted. β-rich proteins therefore change their 3Di states faster per substitution. Illergård et al. (2009) found the same direction for secondary-structure states, which is the only fold-class trend they reported as significant.
- **Length.** k is a little higher for short proteins (tertiles 3.72, 3.22, 3.16).
- **Confidence.** k is higher for the most confidently predicted families (pLDDT tertiles 3.16, 3.14, 3.57).
- **Age.** k rises with the age of the family's sampled clade, from 2.26 (< 100 Mya) to 3.74 (≥ 1,000 Mya).

The plan's estimate of the spread of k between families, σ_k, came out at 0 (variance −0.41).

**Post-hoc interpretation.** Two parts of E3 are artefacts, and one number is more defensible than the headline.
1. **The rise of k with age is the tip artefact.** The extra per-model noise from §4.1 sits on the *terminal* branches, which lead to single proteins. Terminal branches make up most of the young families. On internal branches alone k does not change across age strata (3.74, 4.71, 4.63, 4.65), whereas terminal-branch k rises from 1.54 to 3.53.
2. **σ_k = 0 is an estimator failure, not homogeneity.** The plan subtracted the *mean* sampling variance across families. That mean is dominated by a few imprecise families: mean 0.82 against a median of 0.0039. Robust estimators give the following, all on the log scale [CHECK: confirm σ_k is on the log scale in ANALYSIS_PLAN §4.3]:
   - a between-family SD of about 0.32 (DerSimonian–Laird);
   - about 0.47 among precisely estimated families, whose per-family k spans 1.7–5.2 from the 10th to the 90th percentile;
   - about 0.25 on internal branches.

   Families therefore differ about two- to threefold. This is smaller than, but in the spirit of, Illergård's roughly fivefold variation in "structural response".
3. **Internal-branch k = 4.63 [4.35, 4.93]** is our best central value. It is free of the tip correction. It was reported as a planned sensitivity, but choosing it as the headline is a post-hoc decision.

**Comparison with Illergård et al. (2009).** Both of our values, 3.3 and 4.6, lie inside "3–10×". The comparison is like-for-like in one respect: 3Di and AA both have 20 states and are scored with matched ML models (LG and Q.3Di.AF, each with Γ rate variation). Illergård obtained 3–10× only after grouping amino acids to match the number of states in each structural descriptor (Illergård 2010). It differs in three others:
- 3Di describes local tertiary contacts, not secondary structure, accessibility or RMSD;
- our structures are predictions;
- our families are eukaryotic orthologs, not SCOP domain pairs.

On AlphaFold's effect (§4.5), AF models compress the *offset* between orthologs but not their slope against AA distance; the slope is if anything 1.34× steeper. A ratio estimated from experimental structures could therefore differ in either direction. We did not measure it.

![Figure 7](figures/figD_ratio_k.png)

**Figure 7. On average sequence changes 3–5 times faster than 3Di structure, and faster still relative to structure in α-helical than in β-rich families.** k is the AA branch length divided by the 3Di branch length on the same family tree, after the tip correction; k > 1 means 3Di changes more slowly.
- **(a)** Per-family k for 5,915 families: median 3.24, IQR 2.38–4.34. The pooled ratio of sums is 3.32.
- **(b–d)** Pooled k on all data with 95% family-bootstrap CIs. Open diamonds are the confirmatory E3 values.
  - **(b)** By CATH class: α 4.67, α/β 3.16, mixed 3.21, β 2.55.
  - **(c)** By family age stratum: 2.19, 2.66, 3.27 and 3.74.
  - **(d)** All branches 3.32 against internal branches 4.63.

The trend with age in (c) disappears on internal branches (§4.3).

### 4.3 Time: the pre-registered age trend is real in the data, but it is carried by terminal branches

**Pre-registered result (E4).** The log ratio of 3Di to AA change falls with the log of branch age: **β = −0.150 [98% CI −0.200, −0.101]**, Holm p = 2.6 × 10⁻¹². In other words, the AA:3Di ratio rises with age (β_k > 0). The plan's required robustness checks all keep the sign and significance:
- α × 1.5: −0.138;
- α ÷ 1.5: −0.167;
- deep-node ages from Strassert et al. (2021): −0.142;
- deep-node ages from Betts et al. (2018): −0.150;
- sequence-only alignment: −0.110;
- TM-align: −0.151.

So the result is formally "robust to α and timescale". Calibration simulations showed that the estimator recovers only about 46% of a true slope. The underlying effect, if real, would therefore be larger than β̂.

The plan named one sensitivity as *the* key robustness check (C35): **internal branches only**. On internal branches β = −0.047 [−0.128, 0.034] (p = 0.18). A weighted least-squares version with family-clustered errors also gave a null result: −0.034 [−0.111, 0.044].

**Post-hoc interpretation: the trend is mostly a tip-correction artefact.** Terminal branches get a fixed subtraction ν = η/2 = 0.026 for prediction noise. If the true per-model noise is larger, the young terminal branches, which are short in AA, are left with too much 3Di length. Their ratio looks low, and β is pushed negative. Four lines of evidence support this reading.

1. **Edge-type decomposition** (confirmation set; Figure 8). k on internal branches is flat from 30 My onward. Terminal k climbs from below 1, meaning more 3Di than AA change, which is implausible for real evolution.

   | Branch midpoint age (My) | < 30 | 30–100 | 100–300 | 300–1,000 | > 1,000 |
   |---|---|---|---|---|---|
   | k, internal branches | 3.13 | 4.47 | 4.63 | 4.94 | 4.68 |
   | k, terminal branches | 0.73 | 2.61 | 3.57 | 4.16 | – |

2. **The excess on terminal branches is roughly constant per branch** (0.014, 0.017, 0.013, 0.008), as an additive per-tip artefact predicts.

3. **β depends strongly on ν.** The pre-specified ν ± 0.005 sensitivity (C35) was missing from the confirmation output. An independent agent computed it afterwards (deviation P1):
   - ν − 0.005: β = −0.213 [−0.264, −0.162];
   - ν + 0.005: β = −0.080 [−0.133, −0.027].

   The sign holds inside the pre-specified range. However, β moves about twice as much per step as anticipated, and it reaches zero at ν ≈ 0.036. That value lies inside the plausible range: the mean (not median) identical-pair floor gives ν = 0.031, and the ν that makes terminal and internal k equal is 0.040.

   | ν (tip correction) | 0 | 0.013 | 0.021 | **0.026 (plan)** | 0.031 | 0.035 | 0.040 | 0.045 |
   |---|---|---|---|---|---|---|---|---|
   | E4 β | −0.352 | −0.291 | −0.213 | **−0.150** | −0.080 | −0.014 | +0.082 | +0.168 |

   Above ν ≈ 0.045, 12–20% of cells are censored and β is unreliable. The internal-only β does not depend on ν, because internal branches get no correction.

4. **Simulated worlds.** We passed synthetic 3Di data through the production engine on the real confirmation design.
   - A world with curvature c = 0.72 reproduces the all-branch β (−0.149). It predicts an internal-only β of −0.138, which lies outside the observed CI.
   - A linear world with an extra 0.014 per tip reproduces both the all-branch β (−0.134) and the flat internal k, with internal β = +0.021 and terminal k 0.71 in the youngest bin.

   Only the tip-excess world matches everything observed.

The node-based view (Figure 5) shows the same thing from another angle. Over nodes aged 10–1,000 Mya, floor-corrected 3Di distance grows with roughly the square root of time (log–log slope 0.46), while AA distance grows almost linearly (0.85). At nodes younger than 30 My, the corrected 3Di distance even exceeds the AA distance. That is the regime in which per-model noise is larger than η.

The planned negative control for E4 permutes branch ages within families (D13). It was validated in simulation (V6) [CHECK: its result on the confirmation set is not in E4.json; report it or state that it was not run].

**What we conclude.** *Shown:* the pre-registered age trend in the AA:3Di ratio exists in the data and survives every pre-specified check, including ν ± 0.005. *Believed:* it is mostly an artefact of residual AlphaFold noise on terminal branches. On internal branches, the ratio of sequence to 3Di change is constant at about 4.5–4.9 from 30 My to more than 1 Gy. Any residual age trend is ≤ 0.05 in β (internal SE ≈ 0.035).

![Figure 5](figures/figC_time.png)

**Figure 5. Against time, sequence distance grows almost linearly but floor-corrected structural distance grows more slowly, mostly because of young, noisy comparisons.**
- **(a)** Distance between 1:1 orthologs against the age of their last common ancestor (TimeTree median, log axes), for 117 panel nodes with ≥ 5 families.
  - Each point is the median over families of the per-family median pair distance at that node.
  - Series: AA (blue), 3Di minus η (orange), LDDT distance minus η (green).
  - Lines are bin medians. The dashed line has slope 1, which a constant-rate clock would follow.
  - Nodes ≥ 1 Ga (shaded) were not used for the slopes.
- **(b)** AA/3Di ratio per node. Bin medians are 0.31 (< 10 My), 0.85, 1.99, 3.22 and 3.96 (300–1,000 My), and 3.93 at ≥ 1 Ga. The dashed line is k = 3.32.

Log–log slopes over 10–1,000 Mya: AA 0.85, 3Di 0.46, LDDT 0.50. The low ratio at young nodes reflects per-model noise that a single floor η does not remove (compare Figure 8). Per-node ages are not shown in tables, because of TimeTree licensing.

![Figure 8](figures/figE_terminal_internal.png)

**Figure 8. The age trend in the sequence-to-structure ratio is carried by terminal branches; on internal branches the ratio is flat from about 30 My onward.**
- **(a)** Pooled k against branch midpoint age (TimeTree, log axis) with 95% family-bootstrap CIs, for terminal branches (leading to one protein), internal branches and all branches.
  - Terminal k is 0.46 below 10 My and rises to about 4.2 by 300–1,000 My.
  - Internal k is about 4.5–5.0 from 30 My onward; the two youngest internal bins are noisy (2.6, 3.1).
- **(b)** Share of terminal branches per age bin: 88% below 10 My.

Data are all 6,112 families that passed the branch exclusions. The confirmation-only values are in the table in §4.3.

### 4.4 Lineages and the neutral clock: no evidence for lineage-specific structural rates on internal branches, at low power

**Pre-registered result (E5).** E5 asks whether some lineages change structure faster than their sequence change predicts. It uses 162 lineages and 146 sister-lineage contrasts. Both of its statistics depart from the neutral null (N0):
- The sister-contrast slope of 3Di change on AA change is 0.57 [98% CI 0.10, 1.05]; N0 needs 0.8–1.25.
- The excess lineage SD is σ_L = 0.31 (variance 98% CI [0.046, 0.150]); N0 needs < 0.10.

The pre-registered guard is to rerun the test without the parts most exposed to artefacts. On internal branches only, the excess disappears: slope 1.06 [0.79, 1.34], excess variance [−0.022, 0.117]. Under the plan, an excess that does not survive the guards is classed **"artefact suspect"**, not "lineage-specific rates".

**Power.** The plan required ≥ 80% power before any joint interpretation. A 1.5-fold 3Di acceleration spiked into 10% of lineages was detected in only 24% of tuning-set replicates, and the null rejection rate was 18%. As stated before the results, *accepting* N0 within [0.8, 1.25] was practically out of reach: even under exact N0 the CI half-width is about ±0.25. The margins were not widened.

**What can be said.** *Shown:* there is no evidence of structure-specific lineage rates on internal branches. The apparent excess sits on species-level terminal branches, where the per-model artefact of §4.1 lives. *Not shown:* that N0 holds. This is absence of evidence from a test that would miss a 1.5-fold effect three times out of four. Neutrality in the population-genetic sense was not tested at all: a constant ratio is equally expected if selection constrains both alphabets proportionally.

### 4.5 Two bias checks: AlphaFold compresses ortholog divergence, and in-cluster pairs describe a different range (exploratory)

**AlphaFold artefact (X4).** We compared AlphaFold models with experimental structures on the same residues. Ortholog pairs with X-ray structures on both sides, released after AlphaFold2's training cutoff, look **30–40% more similar as AF models than as crystal structures** (Figure 9). Data: 197 pairs in 160 families, 87 of them human–mouse; median AA distance 0.29.
- **Size of the compression.** Mean AF divergence relative to mean PDB divergence, with 95% CIs:
  - 3Di: 0.59 [0.52, 0.66];
  - 1 − LDDT: 0.70 [0.62, 0.76];
  - 1 − TM: 0.69 [0.60, 0.78];
  - core RMSD: 0.62 [0.51, 0.75].

  84–94% of pairs have smaller AF divergence. The conclusions are the same with TM-align, with chains released before the cutoff, and on tuning-split families only.
- **Shape of the bias.** It is mostly a **constant offset**, not a compression that grows with divergence. AF divergence rises with AA distance at least as fast as PDB divergence: the slope ratio for 3Di is 1.34 [1.15, 1.57]. So AlphaFold does not flatten the structure-versus-sequence curve; if anything, it removes an additive noise-like component that is present in experimental pairs.
- **Same-protein floor.** An AF model differs from its own crystal structure by a median 3Di distance of 0.241. That is *more* than two AF models of different orthologs differ (median 0.167), and close to half the PDB–PDB ortholog distance (0.344). This fits the idea that the offset reflects experimental and conformational variation plus AlphaFold "idealisation", rather than lost evolutionary signal.
- **Training cutoff.** Chains released before the cutoff show a slightly larger AF–PDB gap (+0.050 [0.015, 0.087] in 3Di on 50 matched pairs). This is consistent with mild memorisation, but it is confounded by entry age and construct.
- **Limits.** The sample is small and shallow: only about 10% of pairs have AA distance > 2. No PDB-versus-PDB floor for the same protein was computed. 118 of the 197 pairs belong to confirmation families (logged as C36; no tuning value used them).

**Selection effect of in-cluster pairs (X7, P-clu).** Pairs taken from inside AFDB clusters lie in a different, structurally pre-filtered range of identity from orthologs (Figure 6). The AFDB clusters release 27.7 M member-to-representative alignments.
- **Identity ranges.** Their median sequence identity is 0.21, and only 4.4% have identity ≥ 0.5. Our orthologs span 0.47–0.99 (5–95%).
- **Shape inside clusters.** Median TM-score rises steeply up to about 35% identity and then plateaus at about 0.83. The plateau is mainly a prediction-quality ceiling: at identity ≥ 0.6, median TM is 0.58 for representatives with pLDDT < 60 and 0.97 for pLDDT ≥ 90.
- **What looks conserved at low identity.** At low identity, cluster pairs still have median TM 0.66 at 21% identity. They look structurally conserved *because* membership required structural similarity.
- **Where the sets overlap.** At identity ≥ 0.4, orthologs are slightly *more* similar than cluster pairs (TM 0.86 vs 0.83 at 51%). This is at least partly because our orthologs are scored on confident core residues only.
- **Filters not honoured.** Only 43.6% of the released alignments meet the nominal clustering filters (≥ 90% coverage, E < 0.01). So the in-cluster sample is a mixture, not a clean truncation.

The honest summary: the overlap between the two pair sets is small, so Figure 6 shows *that* in-cluster curves describe another population, not *how much* selection distorts a given curve. The TM-scores are also normalised differently in the two sets.

![Figure 6](figures/figB_selection_bias.png)

**Figure 6. Pairs from inside AFDB structural clusters lie in a different identity range from orthologs, and their low-identity end looks conserved because it was selected for structural similarity.** Median structural similarity, with IQR bands, against sequence identity of aligned residues (0.02 bins). Blue: our 1:1 orthologs (873,628 pairs). Orange: AFDB cluster member-to-representative pairs (27.7 M; green dashed = eukaryote–eukaryote only).
- **(a)** TM-score: a hockey stick, rising to identity ≈ 0.35 and then plateauing at ≈ 0.83.
- **(b)** LDDT.
- **(c)** Identity distributions: cluster pairs are concentrated at 0.10–0.35, orthologs at 0.47–0.99.

Grey zone: cluster alignments below 10% identity are short or spurious. TM normalisation differs (orthologs: shorter protein; clusters: alignment length), and ortholog values use core residues with pLDDT ≥ 70.

![Figure 9](figures/figF_af_vs_pdb.png)

**Figure 9. AlphaFold models make orthologs look more similar than their experimental structures do, by a roughly constant offset.** Each point is one ortholog pair with X-ray structures on both sides. The x-axis is divergence between the two crystal structures, and the y-axis is divergence between the two AF models on exactly the same residues (3Di+AA alignment). Panels show 1 − TM-score, 1 − global LDDT, 1 − core LDDT, 3Di ML distance and core RMSD. Blue: chains released after AlphaFold2's training cutoff (n = 197, the primary set). Orange: earlier chains (n = 485). Points below the dashed identity line mean AF is more similar. Post-cutoff AF/PDB ratios: 3Di 0.59, LDDT 0.70, TM 0.69, RMSD 0.62.

---

## 5 Results: Goal 2, AFDB clusters in the light of evolution

**Main point.** Duplication leaves a measurable structural signature. The faster-evolving copy changes more in structure, which replicates the direction of Warwick Vesztrocy et al. (2025), though the size of the effect is comparable to AlphaFold noise. AFDB structural clusters are not orthology units: mostly because of length differences, they split and merge ortholog groups. In two case studies, the cytoskeletal folds and protein kinases, we found no robust "structural burst" on any lineage.

### 5.1 Duplication (E6): detected, but only the more-diverged-copy effect is robust to alignment

E6 was the positive control for Goal 1. It has two parts.

**E6a: paralogs vs orthologs.** At matched AA distance, paralogs diverge more in 3Di than orthologs.
- Excess: **0.030 3Di subst/site [98% CI 0.026, 0.033]**, from 8,969 paralog groups in 3,525 families. Pre-registered decision: "detected" (Holm p < 10⁻¹⁵).
- **Not alignment-robust.** Under the sequence-only alignment the excess is −0.011 [−0.024, 0.002]; under TM-align it is 0.027 [0.003, 0.051].
- **Post hoc, a second reason for caution.** E6a compares paralog pairs with the fixed-η reference curve of E1, which §4.1 shows is misspecified. We therefore do not claim E6a.

**E6b: the more-diverged copy (MDO).** After a duplication, the copy with the longer AA branch has more 3Di divergence than its sister copy, the least-diverged ortholog (LDO).
- Excess: **0.017 [0.012, 0.022]** (6,027 groups).
- The copies were chosen with one set of outgroups and one half of the alignment columns, and the excess was measured with a different outgroup set on the other half (§3.9). This prevents regression to the mean.
- It is positive under all three alignments (0.017, 0.032, 0.032) and with LDDT distance (0.035 [0.029, 0.040]).
- This **replicates the direction** of Warwick Vesztrocy et al. (2025, PMID 40813247). Using PANTHER gene trees and AFDB models, they found lower LDDT for the MDO than for the LDO in asymmetric duplications (median 0.723 vs 0.766). Our replication uses independent orthology (Ancestra), a different structural measure (3Di ML distance) and a split-half design.

**Caveat (post hoc).** The E6b effect, 0.017–0.035, is the same size as the per-model AlphaFold excess (§4.1). An MDO with lower pLDDT or a shallower input alignment would produce it without any real structural change. We did not adjust for the pLDDT difference between copies, so the direction can be claimed but its biological size cannot.

### 5.2 AFDB clusters vs orthology (X5): clusters split and merge ortholog groups, mainly by length

**Mapping.** Data: 71,961 ortholog groups (HOGs) with at least two members in non-singleton AFDB v3 clusters (Figure 10).
- **Splitting.** **70%** of these HOGs span more than one cluster. On average only 36% of a HOG's members are in its largest cluster. Splitting grows with HOG size, from 42% of two-member HOGs to 100% above 50 members.
- **Merging.** **48.5%** of clusters with ≥ 2 panel members contain more than one HOG.
- **Both at once.** 61% of HOGs, holding 94% of the members, are both split and merged. Only 14% map one-to-one.

**Why HOGs split.** A logistic model on 279,196 within-HOG pairs (48,130 HOGs; McFadden R² 0.40, AUC 0.90) finds that splitting is driven by:
- length difference: log-odds 4.06 per SD of |log length ratio| (R² 0.31 alone);
- AA divergence: 1.17 per SD (R² 0.20 alone).

A change in CATH domain architecture adds nothing once length is in the model (0.02, CI includes 0; Figure 11). This is what the 90% bidirectional-coverage rule of Foldseek clustering (Barrio-Hernandez et al. 2023) predicts: a domain gain or loss separates proteins through the length change it causes.

**Cluster ages are not orthology ages.** The AFDB clusters come with a taxonomic last common ancestor (LCA). For only **26.9%** of clusters does it coincide with the root of the dominant HOG, and 50.7% are at least consistent once censoring is counted. Two opposite effects are at work:
- **Splitting** makes the cluster LCA too young: restricted to our panel, 70% of cluster LCAs are younger than the HOG root.
- **Merging** and taxa outside the panel make it too old.

**Domain architecture over time.** Even nearly identical orthologs (AA distance < 0.02) differ in TED/CATH architecture 21–25% of the time. This is annotation noise. Architecture change is flat with time (0.24 below 100 Mya, 0.22 above 1,000 Mya) and rises only modestly with AA distance (0.22 to 0.36). Cluster membership, by contrast, separates sharply with time: 14% of ortholog pairs younger than 100 My are in different clusters, against 56% at 700–1,000 My. At matched AA distance, paralogs change architecture more often than orthologs (log-odds +0.23 [0.20, 0.26]).

All X5 analyses are exploratory, with unadjusted 95% CIs, and the pair-level analyses use only tuning-split families.

![Figure 10](figures/figG_hog_cluster_mapping.png)

**Figure 10. AFDB structural clusters cut across ortholog groups in both directions.**
- **(a)** Number of AFDB v3 clusters per HOG, for 71,961 HOGs with ≥ 2 eligible members (log axes). 70% span more than one cluster.
- **(b)** Splitting grows with HOG size. Blue: fraction of HOGs that span > 1 cluster. Orange: mean fraction of members in the largest cluster.
- **(c)** Number of HOGs per cluster, for 354,716 clusters; the last bar is ≥ 30. 48.5% of multi-member clusters contain more than one HOG.

![Figure 11](figures/figH_hog_split_model.png)

**Figure 11. Ortholog groups split across AFDB clusters mainly because their members differ in length, not because their domain architecture differs.**
- **(a)** Standardised log-odds (with 95% HOG-clustered CIs) that two members of the same HOG fall in different clusters, from a logistic model on 279,196 pairs in 48,130 split HOGs. Length ratio dominates (4.06 per SD), AA distance follows (1.17), and CATH architecture difference is near 0.
- **(b)** Fraction of pairs in different clusters against |log length ratio|, for pairs with the same (blue) or different (orange) CATH architecture. A different architecture adds only 5–9 percentage points within each length bin.

### 5.3 Case studies: no robust structural burst in the cytoskeleton or kinases; per-branch ratios depend on method

The case studies asked whether particular lineages carry unusually much structural change per sequence change. For each branch the measure is ρ: its 3Di/AA ratio relative to the whole tree, so that ρ = 1 is typical. All case-study results are exploratory and descriptive.

**Cytoskeleton (actin and tubulin folds; X6).** The result is negative. Data: 400 proteins per superfamily, from eukaryotes and 100 prokaryotes.
- **Actin.** The path from prokaryotic actins to the eukaryote + Asgard ancestor gives ρ = 1.40 with MAFFT alignments [0.88, 2.38] and 1.57 with FAMSA. With the structure-guided FoldMason alignment it gives 0.64, and with the pairwise branch engine 0.71.
- **Tubulin.** The stem gives 1.14 (MAFFT) against 0.56 (FoldMason).
- **Alignment error.** Sequence alignments lean high and structure-guided ones lean low, which is the signature of alignment error at depth.
- **Better models and topologies.** A site-heterogeneous AA model (LG+C20) brings the values to about 1. On constrained topologies, which AU tests do not reject (actin p = 0.061, tubulin p = 0.32), the eukaryote-only stem has zero length.
- **The one robust signal is not where expected.** It sits on the branch *below* crenactin (ρ = 2.40 [1.52, 5.71], MAFFT), not on the eukaryotic stem. Its 3Di changes concentrate in the hydrophobic plug, D-loop and C-terminal region. This is descriptive only.
- **Flagged tips.** Within eukaryotes, robust high flags are almost all single terminal branches (e.g. *Trichomonas* β- and γ-tubulin), and are enriched for low-pLDDT models. The pairwise engine supports fewer than half of them.

**Method dependence (Figure 12).** Per-branch AA lengths agree well between the IQ-TREE joint-likelihood fit and the pairwise branch engine (Spearman 0.74–0.80). Per-branch ρ barely agrees (Spearman 0.13 for actin and 0.29 for tubulin; log-scale Pearson correlations of 0.08 and 0.04 in the figure). Single-branch structural "accelerations" are therefore not reliable. This is independent support for the ratio-of-sums pooling used in Goal 1.

**Kinases (protein kinase-like fold; X8).** There is at most a weak tendency. Data: 500 kinase domains.
- **Tyrosine-kinase origin.** The two branches at the origin of holozoan tyrosine kinases have ρ = 1.68 (MAFFT), 1.47 (FAMSA) and 2.00 (FoldMason). That is the 88th–91st percentile of comparable branches, consistent in direction across alignments but never above the 95th percentile. It is a tendency, not a burst, and these branches are short (0.17 AA subst/site).
- **Prokaryotic vs eukaryotic kinases.** Internal branches of prokaryotic kinases have ratios 1.02–1.22× those of eukaryotic kinases, with every CI including 1.
- **Flagged branches.** Only 1 of 12 flagged branches is flagged with all three alignments.
- **Saturation.** Deep kinase-to-atypical-kinase pairs are AA-saturated (97% have p-distance > 0.75).

![Figure 12](figures/figI_cytoskeleton_methods.png)

**Figure 12. Branch-length estimates for single branches agree between methods for sequence but not for the structure-to-sequence ratio.** Internal branches of the actin (top) and tubulin (bottom) trees.
- **Columns 1–3:** IQ-TREE fixed-topology estimates (MAFFT alignment; x) against the pairwise branch engine (3Di+AA alignments; y), for AA length, 3Di length and their relative ratio ρ (log axes; r_log is the correlation on the log scale).
- **Column 4:** ρ for named clades under five methods (IQ-TREE with MAFFT, FAMSA or FoldMason; branch engine with pairwise or MSA distances), with bootstrap intervals where available.

AA lengths agree (r_log 0.74–0.80). ρ does not (0.04–0.08), so only path-level or strongly multiply-supported deviations are interpretable.

---

## 6 Discussion

**Main point.** At the scale of thousands of eukaryotic ortholog families, 3Di structure looks like a slowed-down copy of the sequence clock. On internal branches of the tree its rate relative to sequence is roughly constant, at about one-fifth. The pre-registered tests gave a more dramatic answer: saturation, an age trend and lineage effects. Those turned out to share one technical cause, residual per-model noise in AlphaFold predictions. The main lessons concern how to measure structural evolution with predicted structures, not a new law of structural evolution.

### 6.1 Answers to the three questions

**(1) Is structural divergence linear in sequence divergence?**
- *Shown (pre-registered):* the estimator we committed to finds sublinear growth (ĉ = 0.72 for 3Di, 0.70 for LDDT). It is robust to alignment and to every planned sensitivity.
- *Shown (post hoc, independently checked):* that estimator is biased towards concavity whenever per-protein offsets exist, and such offsets exist in our data. Alternative estimators give 0.73–0.98.
- *Believed:* roughly linear above a small step between non-identical orthologs. The step is partly AlphaFold noise that tracks pLDDT; its remainder is unexplained.

This agrees with Illergård et al. (2009) and Koehl & Levitt (2002), who found linear core change against ML distance. It disagrees with the saturation that Chothia & Lesk's exponential implies once it is rewritten on a substitution axis. Our range, however, stops at 2 AA substitutions per site, where saturation would start to matter.

**(2) Is structure 3–10× more conserved than sequence?**
- *Shown:* yes, within that range. k = 3.3 over all branches (pre-registered) and 4.6 on internal branches (planned sensitivity; chosen as headline post hoc). β-rich families change 3Di about 1.8× faster per substitution than α-helical families, as Illergård reported for secondary structure.
- *Believed:* families vary about two- to threefold around about 4.6.

The new element is the like-for-like comparison: two 20-state alphabets under matched ML models. Illergård could approximate this only by grouping amino acids. Our number is still specific to 3Di, which captures local tertiary contacts. A global measure such as TM-score or RMSD would give a different ratio, and so would experimental structures.

**(3) Is there a structural molecular clock, and is it neutral?**
- *Shown:* after the AlphaFold artefact is removed, there is no robust trend of the ratio with age, and no lineage-specific structural rate on internal branches. Both are absence-of-evidence results: the internal-branch age slope has SE ≈ 0.035, and E5's power was 0.24.
- *Not shown:* neutrality. A constant structure-to-sequence ratio fits the view that structural change is mostly a passive by-product of substitutions (Echave & Wilke 2017). It is equally expected if purifying selection acts on both proportionally.
- *Believed:* on internal branches, 3Di distance is a rescaled AA clock and adds little independent timing information.

This sits between the two published views. Pascual-García et al. (2019) found larger clock violations in structure, linked to function change. Moi et al. (2025) found structure-informed trees the most clock-like. Our design, with orthologs only and the ratio rather than absolute rates, removes much of the functional divergence that Pascual-García's superfamilies contained. The duplication signal (E6b) points in their direction: where function is expected to change, after duplication, extra structural change appears. At the depth of terminal branches, our data match the observation that 3Di signal is noisier than AA signal at short times (Fullmer et al. 2025, preprint). In our case the cause is identifiable: prediction noise.

### 6.2 What this means for structural phylogenetics and dating with AlphaFold structures

Five practical consequences follow.

1. **Terminal branches of 3Di trees are inflated.** Each AF model adds about 0.01–0.017 3Di subst/site beyond what an identical sequence would, and more at low pLDDT. On young splits this is comparable to, or larger than, the real signal. Below about 30 My, 3Di branch lengths are mostly noise. A single global floor taken from identical sequences under-corrects. Floors should be specific to each protein or family, or modelled as a tip-specific error term. FoldTree's filtering of low-pLDDT families (Moi et al. 2025) points the same way.
2. **Deep 3Di branches carry information proportional to AA branches.** At about 4.6 AA substitutions per 3Di substitution, a 3Di tree length of 0.5 corresponds to about 2.3 AA substitutions per site. That is where AA distances become unreliable. This supports using 3Di for remote homology and deep phylogeny (van Kempen et al. 2024; Puente-Lelievre et al. 2023). It does not support treating 3Di as an *independent* clock: its rate is tied to the sequence rate.
3. **AlphaFold models understate divergence between orthologs by 30–40%,** but do not flatten the curve. Absolute structural distances between AF models are therefore not comparable to distances between crystal structures. Slopes and ratios are less affected, although the experimental comparison set was small (197 pairs).
4. **AFDB clusters should not be used as orthology groups or for dating.** Membership is driven by the length-coverage rule. Cluster LCAs coincide with orthology-based ages for only about a quarter of clusters. And curves fitted inside clusters describe a structurally pre-selected population. This qualifies, rather than contradicts, Barrio-Hernandez et al. (2023), whose clusters were designed for structural search and annotation, not for evolutionary inference.
5. **Single-branch structural "bursts" should not be trusted without multiple methods.** In both case studies, per-branch ratios depended on the alignment and on the estimation method. Even strong flags were enriched for low-confidence models.

### 6.3 How the agentic process performed

**Main point.** The safeguards worked as designed for what they were designed to catch. Pre-registration protected the decisions. A pilot run and simulated data caught three flaws before confirmation. Independent re-derivation confirmed the arithmetic. But the most consequential flaw, a fixed noise floor taken from one kind of control, passed the design, three review rounds and the calibration simulations. It was found only by a critical statistical review after confirmation.

**What went well.**
- **Pre-registration and the deviation log.** Every change to the plan before confirmation (28 numbered items) was logged with its reason. The confirmation script refused to run unless the checksums of the plan, tuning file and split matched. The six decisions are therefore reported unchanged even where we now doubt their interpretation, and the reader can see exactly what was decided before and after the data.
- **The pilot on real data.** A full run on 20 proteomes found that TM-align misregisters close homologs (D22). No reviewer had predicted this.
- **Simulation before confirmation.** Synthetic worlds exposed two estimator biases, a false "concave" call (D31) and a spurious age slope (D32), and both were fixed before confirmation. Margins were not widened when simulations showed that E5 could not accept its null.
- **Independent checking.** One agent re-derived all 42 confirmatory numbers from the plan text with its own code (V11). Another tried to refute the post-hoc review and partly succeeded: it showed that "c ≈ 0.9–1" was overstated, and the report uses its range of 0.75–1.0.
- **Parallel specialists with written reports.** Twenty-two reports, read by the PI as summaries, kept the PI's context small. A node move after 10.9 h of compute needed no rebuilding.

**What went badly.**
- **The noise floor was a design flaw, and the safeguards did not catch it.** The design fixed η from identical-sequence pairs and treated it as known. The three review rounds questioned AlphaFold artefacts but not whether *non-identical* sequences carry extra model noise. The calibration simulations also used a fixed floor, so they could not reveal the problem: they verified the estimator under the noise model we had assumed. Pre-registration does not protect against a wrong model of the noise; it only guarantees that the wrong model is applied honestly.
- **The signal was visible before confirmation and was not acted on.** In the tuning dry run, the internal-branch-only E4 slope was already near zero (−0.014). The exploratory AF-versus-PDB analysis had already shown that AF noise behaves as an offset. Lesson: compare every key sensitivity with the primary result on the tuning data before freezing.
- **A planned sensitivity was missing from the output.** The ν ± 0.005 check (C35) was specified but not wired into the confirmation output. The review found the gap and it was computed afterwards (P1). V11 checked the numbers that existed, not whether all the specified numbers existed.
- **Some design choices made tests weak by construction.** E5 could not accept its null (power 0.24), and E1 detected c = 0.90 only 55% of the time. Both were known before confirmation and reported honestly. Ideally they would have changed the design rather than only its caveats.
- **Exploratory work competed with the critical path,** and one exploratory analysis touched confirmation families (C36, 118 pairs). Neither affected a decision.

---

## 7 Limitations

**Main point.** All conclusions concern AlphaFold models of well-predicted, mostly ordered eukaryotic proteins, measured with one structural alphabet over 0.05–2 AA substitutions per site. The two key post-hoc reinterpretations, the per-model noise and the internal-branch ratio, were checked by a second agent but were not pre-registered.

**Structures and measures**
- All structures are predictions. AF models compress ortholog divergence by 30–40% compared with crystal structures, and this was measured on only 197 pairs, mostly human–mouse, with a thin deep end.
- The per-model excess noise depends on pLDDT and is only partly removed by the floor η = 0.052. No pre-registered estimator modelled it, so the post-hoc estimates of curvature (0.73–0.98) and of internal-branch k (4.6) remain interpretations.
- 3Di describes local tertiary contacts. Ratios for global measures (TM-score, RMSD) or for whole-domain rearrangements may differ. The LDDT analysis (E2) agrees in direction but shares the noise-floor issue.
- Distances use only the confident core (both residues pLDDT ≥ 70). Loops, disorder and indels are excluded, and they are where much structural change happens.

**Sampling**
- The pLDDT ≥ 70 filter removed 38% of genes. The results describe ordered proteins.
- Only 224 families were in the youngest stratum (< 100 Mya), and terminal branches dominate young ages, so the shallow end of the time axis is thin and noisy.
- Sequence divergence was capped at 2 AA substitutions per site per family. Saturation beyond that range was not tested.
- Proteins longer than 2,700 residues have no AFDB model. Several model organisms (rat, zebrafish, chicken, *Xenopus*) are missing because their proteomes no longer match AFDB v4. Polyploid species are included, handled through the one-copy rule.
- The 4,891 analysed families were those that passed the misfit and fit-quality screens. Families with hidden paralogy or strong gene-tree discordance are under-represented.

**Orthology and time**
- Ancestra is unverified AI-generated software. Its ortholog calls checked well against Ensembl (148/149 human–mouse 1:1 recovered), but its duplication timing did not (36% agreement), so it was never used. 4% of pairs within one-to-one sets are labelled out-paralogs.
- TimeTree resolves deep eukaryotic splits crudely. 52% of species pairs fall at the root age, 6 internal nodes have imputed ages, and polytomies replace splits shorter than 5 My. The alternative deep-age calibrations (Strassert 2021, Betts 2018) changed E4 little, but they cover only 3–11 nodes.

**Statistics**
- The E1/E2 curvature test had low power (c = 0.90 detected 55% of the time). E4 recovers only about 46% of a true slope. E5 had power 0.24 and could not accept its null.
- E3's pre-registered between-family SD estimator failed (σ_k = 0). The robust values (0.25–0.47 on the log scale) are post hoc.
- E6 is not alignment-robust for paralogs vs orthologs (E6a). The more-diverged-copy effect (E6b) was not adjusted for differences in model quality between copies.
- The E4 age-permutation negative control is not in the confirmation output [CHECK].
- Minor issues found by V11 that do not change any result:
  - two families with an undefined misfit ratio were kept;
  - E5's contrast count is overstated by one (145 enter the slope);
  - the description of D31 does not match the code (the code simulates tree-additive design with pair-level noise).
- A small orientation asymmetry in the sequence-only alignment (C42) adds noise to the family-selection distance. The split was frozen and was not corrected.

**Goal 2**
- All Goal 2 analyses except E6 are exploratory, with unadjusted CIs, and were not independently re-derived.
- The cluster audit uses AFDB v3 clusters and TED/CATH annotations. Architecture labels differ for 21–25% of near-identical orthologs.
- The case studies rest on 400–500 sequences per family, one Asgard genome, no Asgard tubulin, and AA-saturated deep branches. Per-branch ratios depend on the method.

**Process**
- Everything was done by AI agents within 48 h of compute, and no human has verified the methods or results.

---

## Data and code availability

**Main point.** Code, configuration, pre-registration records, confirmatory results and agent reports are in the git repository. Raw data are not, but every dataset is listed with source, version, checksum and licence in DATA.md and can be re-downloaded with the provided scripts.

**Repository layout**

| Path | Contents |
|---|---|
| `PROMPT.md`, `RESEARCH_TREE.md`, `LESSONS.md` | Brief, research record with milestone decisions, process rules |
| `LITERATURE.md`, `DESIGN.md`, `reviews/` | Literature review, design (v3) and all design and results reviews |
| `ANALYSIS_PLAN.md`, `DEVIATIONS.md`, `RESULTS_LEDGER.md` | Frozen plan, every later change, and a ledger of reported numbers (script, commit, input checksum; plan, tuning and split checksums) |
| `config/` | Seeds (root 20261030) and `tuning.yaml` (all tuned values with their sources) |
| `src/strucevo/` | Distance library (ML distances on GPU, LDDT), pair pipeline, statistics engines |
| `scripts/` | Pipeline stages: `panel/`, `ancestra/`, `hogs/`, `pairs/`, `floors/`, `tuning/`, `endpoints/` (`run_confirm.py`), `verify/` (V11), `review/` and `posthoc/`, `exploratory/` (figures), `h1e/`, `pclu/`, `goal2/`, `casestudy/`, `kinase/` |
| `results/` | `endpoints/` (confirmatory JSON: E1–E6, Holm, interpretation, run guards), `verify/`, `posthoc/`, `tuning/`, and exploratory results by analysis |
| `figures/` | All figures, with the binned data behind the main figures (`figures/main/*.tsv`) |
| `agents/` | The written report of every subagent |
| `tests/` | Automated tests (distance library and pipeline invariants) |

**Data.** `DATA.md` lists each dataset: AFDB clusters v3/v6, AFDB v4 FoldComp archives, UniProt 2026_03 reference proteomes, NCBI Taxonomy, TED, Pfam 38.2, PDB/SIFTS, TimeTree and the deep-age literature tables. `data/raw/` and `data/interim/` are not committed; they are large and regenerable.

**TimeTree.** TimeTree data may not be redistributed. Divergence times and per-pair or per-node ages are therefore never committed. Committed tables carry age-bin labels only. To rerun, download the TimeTree 2015 tree and query the pairwise service yourself (`scripts/panel/`).

**Rerunning.**
1. Create the environment from `envs/` (conda; Foldseek with GPU support).
2. Download the data with `scripts/download_afdb_clusters.sh`, `scripts/download_domain_annotations.sh` and the panel scripts.
3. Run Ancestra (commit ffb0bbe), then the HOG build, pair metrics, floors and tuning stages in order. Each long stage writes a completion marker and can be resumed (see `NODE_MOVE.md`).
4. The confirmation run is `CONFIRM_RUN=1 python scripts/endpoints/run_confirm.py`. It refuses to run if the plan, tuning or split checksums differ from the ledger or the code tree is not clean, and it refuses to rerun without a logged deviation.
5. The main figures are produced by `CONFIRM_RUN=1 python scripts/exploratory/figures_main.py`.

[CHECK: list exact per-stage commands and approximate runtimes in README.md before release; Ancestra took 5.0 h wall time and 161 GB of memory; the confirmation run took 1.8 h.]

---

## References to add

These are not yet in the first half's reference list. Identifiers are from LITERATURE.md unless marked.

- Echave J, Wilke CO (2017). Biophysical models of protein evolution: understanding the patterns of evolutionary sequence divergence. *Annu Rev Biophys* 46:85–103. PMID 28301766. [CHECK: title, volume and pages; only the PMID is in LITERATURE.md]
- Fullmer MS, Puente-Lelievre C, Matzke NJ (2025). bioRxiv. doi:10.1101/2025.06.30.662300. [CHECK title; already listed in the first half with the same flag]
- Puente-Lelievre C, et al. (2023). Tertiary-interaction characters enable fast, model-based structural phylogenetics beyond the twilight zone. bioRxiv. doi:10.1101/2023.12.12.571181.
- Betts HC, et al. (2018). Integrated genomic and fossil evidence illuminates life's early evolution and eukaryote origin. *Nat Ecol Evol* 2:1556–1562. [CHECK: identifier from agents/deep_ages_report.md, not in LITERATURE.md]
- DerSimonian R, Laird N (1986). Meta-analysis in clinical trials. *Control Clin Trials* 7:177–188. PMID 3802833. [CHECK: not in LITERATURE.md]
- Le SQ, Gascuel O (2008). An improved general amino acid replacement matrix. *Mol Biol Evol* 25:1307–1320. PMID 18367465. [CHECK: not in LITERATURE.md; the LG model]
- Holm S (1979). A simple sequentially rejective multiple test procedure. *Scand J Stat* 6:65–70. [CHECK]
- Cook JR, Stefanski LA (1994). Simulation-extrapolation estimation in parametric measurement error models. *J Am Stat Assoc* 89:1314–1328. [CHECK; SIMEX]
- Shimodaira H (2002). An approximately unbiased test of phylogenetic tree selection. *Syst Biol* 51:492–508. PMID 12079646. [CHECK; AU test]
- Katoh K, Standley DM (2013). MAFFT multiple sequence alignment software version 7. *Mol Biol Evol* 30:772–780. PMID 23329690. [CHECK]
- Deorowicz S, Debudaj-Grabysz A, Gudyś A (2016). FAMSA: fast and accurate multiple sequence alignment of huge protein families. *Sci Rep* 6:33964. [CHECK; the case study used FAMSA 2]
- Minh BQ, et al. (2020). IQ-TREE 2. *Mol Biol Evol* 37:1530–1534. PMID 32011700. [CHECK; the case studies used IQ-TREE 3, whose citation should be confirmed]
- Zhang C, Shine M, Pyle AM, Zhang Y (2022). US-align. *Nat Methods* 19:1109–1115. [CHECK; used for case-study superpositions]
