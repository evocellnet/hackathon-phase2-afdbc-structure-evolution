# LITERATURE.md: protein structure evolution, AFDB clusters and orthology

Compiled 2026-10-30 by the literature agent for the phase2-afdbc-evolution project. Sources were PubMed (NCBI E-utilities and the PubMed connector), Europe PMC full text, bioRxiv, the publishers' sites and the live data portals. Every reference below was checked against PubMed or the publisher record (PMID and/or DOI given).

**How sure each claim is:** **[FT]** means read in the full text (or, for Illergård 2009, in the first author's PhD thesis summary). **[Abs]** means read in the abstract only. **[Live]** means checked directly on the data portal or API on 2026-10-30. **[Ours]** marks a derivation or inference made by this review, not by the cited authors.

---

## 0. Key takeaways for the design (read this first)

1. **Illergård et al. 2009 measured "time" as substitutions per site, not years.** Their "linear" result means: structural change in the protein core grows linearly with ML amino-acid distance (JTT+Γ) [FT-thesis]. It does not test linearity against geological time. Testing a structural *clock* against time is therefore a new question, not a replication.
2. **The "3–10×" number depends on the alphabet.** With raw discrete structural states, structure looked 15–30× more conserved. Only after grouping amino acids into as many states as the structural alphabet does the figure fall to 3–10× [FT-thesis]. Foldseek's 3Di alphabet has **20 states, the same number as amino acids**. This allows a direct like-for-like rate comparison that Illergård could only approximate.
3. **The AFDB cluster paper did no rate or divergence analysis.** Its only "evolutionary analysis" assigned each cluster a lowest common ancestor (LCA) from member presence: 23% cellular organisms, 16.1% Bacteria, 13.5% Eukaryota, 0.5% Archaea, 3.91% species-specific [FT]. It did not look at sequence versus structure divergence, time, orthology, gene trees or domain-architecture change. Those are open.
4. **AFDB clusters are a biased sample of structural divergence.** Cluster members are AFDB50 representatives, so pairs are mostly under 50% identity. Membership also required 90% bidirectional coverage and E < 0.01 [FT]. Any pair drawn from inside a cluster has passed a structural-similarity filter, which truncates divergence and can **create artificial saturation**. Linearity tests must use pairs defined independently of structure (orthologs from Ancestra, OMA or sequence), with clusters used only for stratification.
5. **The released cluster files already contain member-to-representative alignments** with `fident`, `lddt` and `alntmscore` (v3 file 8, 948 MB) [Live]. There is also a **v6 cluster release** on AFDB v6 (UniProt 2025_03; 241M structures) [Live]. IDs differ between v3 (the paper) and v6.
6. **The evidence on structural clocks is mixed.** Pascual-García et al. 2019 found more and larger clock violations in structure than in sequence, linked to function change [Abs]. Moi et al. 2025 (FoldTree) found trees built from statistically corrected 3Di+AA distances had the *lowest* root-to-tip variance [FT]. Fullmer et al. 2025 (bioRxiv) found 3Di *less* clock-like than amino acids in histones [preprint, via web summary]. Whether structure is "clock-like" depends on the metric and needs a direct test.
7. **Saturation has to be modelled explicitly.** Chothia & Lesk's RMSD = 0.40·e^(1.87H) [Abs], once rewritten with Poisson distance, levels off near 2.6 Å [Ours]. FoldTree's correction sends the distance to infinity as Fident approaches ~7%, and as LDDT or TM approach 0 [FT]. Linear and saturating models should be compared formally (AIC/BIC), separately for each metric.
8. **AlphaFold2 adds its own bias.** AF2 is largely insensitive to point mutations [Abs: Buel & Walters 2022; Pak et al. 2023] and draws on overlapping MSAs and templates for homologs. Predicted structures of close orthologs may therefore be artificially similar, which would inflate "structure is more conserved". A PDB-only control subset and an MSA-depth covariate (AFDB v6 now ships MSAs) are needed.
9. **Usable divergence times exist mainly for eukaryotes.** The TimeTree 5 API works, e.g. human–mouse 83.5 Mya [Live]. For prokaryotes it is sparse and highly uncertain: E. coli–B. subtilis 1,909 Mya with a CI of 659–3,158 from 2 studies, and no result for human–E. coli [Live]. Prokaryotes need GTDB tree branch lengths (substitutions, not time) or relative timing.
10. **The least-diverged-ortholog (LDO) paper (PMID 40813247) already uses AFDB structures and Foldseek LDDT**, and finds excess structural divergence in the faster-evolving paralog [FT]. A structure-first LDO test using Ancestra HOGs is a natural extension, but it is not structurally novel in principle.

---

## 1. The two anchor papers

### 1a. Illergård, Ardell & Elofsson 2009

**Citation:** Illergård K, Ardell DH, Elofsson A. *Structure is three to ten times more conserved than sequence—a study of structural response in protein cores.* Proteins 77(3):499–508 (2009). DOI 10.1002/prot.22458, PMID 19507241. The title in the task brief is correct. The article is closed access. Method details below come from Illergård's 2010 Stockholm University PhD thesis (*On the effects of structure and function on protein evolution*, DiVA urn:nbn:se:su:diva-35872, Ch. 8–9), which summarises this paper as Paper I.

**Data and methods [FT-thesis]:**
- **Dataset.** Pairs of domains from **SCOP 1.74**, counted as homologous if in the same SCOP family or superfamily. After filtering for structural quality, the set was homology-reduced to **one protein pair per family or superfamily**. The thesis does not give the exact number of pairs, and we did not verify it.
- **Alignment.** Structural alignment with **STRUCTAL** (Subbiah, Laurents & Levitt 1993). Using structural rather than sequence alignments allowed distant pairs to be included.
- **Divergence axis.** ML evolutionary distance (substitutions per site) under **JTT+Γ** with among-site rate variation, computed with **TREE-PUZZLE 5.2**. They compared this against fraction identity.
- **Structural measures**, restricted to the protein **core**: (i) RMSD; (ii) discrete per-residue descriptors, namely 3-state secondary structure (STRIDE: helix/sheet/coil), 2-state relative accessibility (Naccess 2.1.1, 1.4 Å probe, buried/exposed) and residue contacts. For the descriptors, the number of state changes was counted along the alignment, like a "structural p-distance". The precise core definition is not given in the thesis and was not verified.
- **"Structural response"** is defined as structural change per unit of sequence change. It is the slope of structural divergence against evolutionary distance.

**Findings:**
- Plotting against evolutionary distance instead of fraction identity changes the RMSD relation **"from exponential to linear"** [FT-thesis]. The abstract says that on average protein cores evolve linearly with substitutions per site for all measures [Abs].
- With raw discrete descriptors, structure is **15–30× more conserved** than sequence. After grouping amino acids into the same number of states as each structural alphabet, it is **3–10× more conserved** [FT-thesis]. The headline number is therefore alphabet-dependent.
- **Large variation:** domain families differ about **fivefold** in structural response [Abs].
- Of the structural and functional categories tested, only one trend was significant: **β-rich domains change secondary-structure state faster than α-rich ones**. The authors attribute this to shorter strands and to changes at the ends of secondary-structure elements [Abs].
- The thesis says the cause of the large variance "remains an open question".

**Limitations relevant to us:**
1. There is no geological time. "Linear" refers to the sequence-distance axis.
2. Only pairwise comparisons, with no phylogenetic correction. A single pair per family limits within-family power.
3. Only cores of crystallised SCOP domains. Loops, indels, multi-domain context, orientation between domains and disorder were all excluded.
4. Evolutionary distance was estimated *from structural alignments*. Sequence distance and structural distance therefore share alignment error, which can induce correlation.
5. Data from 2009: SCOP 1.74 and a few thousand families at most. It is biased toward crystallisable, well-studied proteins.

**Implications and testable hypotheses:**
- **H-I1 (replication).** In AFDB-derived ortholog and homolog pairs, core structural divergence is linear in ML amino-acid distance up to some d*, and saturates beyond it. Estimate d* for each metric.
- **H-I2 (alphabet-matched rate ratio).** Compare 3Di and amino acids directly, both 20 states, using the same alignment and the same distance correction. The prediction is a structure/sequence rate ratio somewhere in 3–10×, possibly different for 3Di.
- **H-I3.** About fivefold variation in structural response between families. The β > α trend should replicate.
- **Pitfall.** Do not estimate the sequence distance and the structural distance from one shared alignment without a control. Repeat the analysis with a sequence-only MSA, such as Ancestra/FAMSA, for the sequence axis.

### 1b. Barrio-Hernandez et al. 2023: AFDB clusters

**Citation:** Barrio-Hernandez I, Yeo J, Jänes J, Mirdita M, Gilchrist CLM, Wein T, Varadi M, Velankar S, Beltrao P, Steinegger M. *Clustering predicted structures at the scale of the known protein universe.* Nature 622:637–645 (2023). DOI 10.1038/s41586-023-06510-w, PMID 37704730, PMC10584675. Open access under CC-BY 4.0. Read in full text [FT].

**Pipeline [FT]:**
- **Input:** 214,684,311 structures from the AFDB "UniProt v3" set.
- **Step 1, sequence clustering.** MMseqs2 at 50% identity and 90% bidirectional coverage gives 52,327,413 clusters (AFDB50). The member with the highest pLDDT becomes the representative.
- **Step 2, structural clustering.** Foldseek cluster (Linclust adapted to 3Di; k-mer m = 300, k = 10; then cascaded prefilter, structural Smith–Waterman and set-cover clustering) groups the AFDB50 representatives into 18,661,407 clusters. Criteria: E < 0.01 and 90% coverage of both sequences, with **no sequence-identity threshold**. Run time was 129 h on 64 cores.
- **Step 3.** UniProt "fragment" entries are removed, leaving **2,302,908 non-singleton clusters**. These average 13.05 members and mean pLDDT 71.64. The 13,012,338 singletons average pLDDT 58.95. The current download README says 2.27M after a 2023 fix to an erroneous merge in the fragment step [Live].

**Quality checks:**
- Median within-cluster LDDT is 0.77 and median TM-score 0.71 (member versus representative).
- 66.4% of Pfam-annotated clusters are 100% Pfam-consistent. As LDDT rises from 0.5 to 0.9, perfect consistency rises from 29.2% to 93.8% for Pfam and from 40.3% to 81.3% for EC level 4.
- ECOD benchmark: pairwise consistency within clusters is 99.6% (A-group), 98.6% (X), 97.4% (H), 96.8% (T) and 72.8% (F). With E relaxed to 10, these fall to 69.7/55.7/53.3/51.9/36.6%. The strict E-value is what keeps clusters mostly homologous.

**Dark clusters:** 711,705 clusters (30.9%) have no PDB, Pfam or TIGRFAM signal, but they hold only **4.06% of AFDB proteins**. A high-confidence subset (33,842 clusters with pLDDT > 90) was examined for pockets (AutoSite) and function (DeepFRI). 32.7% of the pockets covered more than 40% of the protein, a sign of non-compact, probably wrong predictions. Membrane and transporter functions were enriched.

**Evolutionary and taxonomic analysis actually performed:**
- **Method.** Cluster LCA with `mmseqs lca` over member taxIDs, shown as Pavian Sankey plots.
- **Results.** 23% of clusters map to cellular organisms (532,478 clusters span the tree of life), 16.1% to Bacteria, 13.5% to Eukaryota, 0.5% to Archaea, and **3.91% are species-specific**. Species-specific clusters are smaller, 56% dark, about 40 aa shorter and have similar pLDDT (69.35 vs 71.73). The authors describe them as possible low-quality predictions or de novo genes.
- **Human-specific clusters.** Of 13, 9 have pLDDT < 70; the rest are artefacts or fragments. **No evidence for new human-specific structural clusters.**
- **Cross-kingdom examples.** BPI vs bacterial proteins (TM 0.81) and the E. coli YceB match; AIM2-like proteins in Lokiarchaeota and Clostridium; histone-like proteins; gasdermin-like domains; a Frag1/DUF998 link; anthrax_toxA-like domains in algae.
- **Domain families.** Built from local Foldseek hits among cluster representatives (boundary clustering plus walktrap): 7,599 families match Pfam, 2,032 match DUFs and 10,722 are unannotated, among families with ≥5 members.

**What they did *not* do** (our assessment from the full text):
- No relationship between sequence divergence and structural divergence.
- No rates and no calibration to time.
- No gene trees, orthology/paralogy, duplication or loss.
- No test of whether LCA ages are artefacts of HGT, contamination or homology-detection failure.
- No domain-architecture evolution beyond grouping domain families.
- No 3Di-based phylogenetics.

Their discussion states, without testing it, that "protein structure is conserved for longer periods of evolutionary time than protein sequences".

**Released data [Live, https://afdb-cluster.steineggerlab.workers.dev, CC-BY 4.0; files on `steineggerlab.s3.amazonaws.com/afdb-cluster/`]:**

| v3 (paper) file | Content | Size |
|---|---|---|
| 1-AFDBClusters-entryId_repId_taxId.tsv.gz | members of non-singleton clusters, with taxID | 261.8 MB |
| 2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz | one row per cluster: dark flag, size, lengths, pLDDT, LCA | 39.1 MB |
| 3-sapId_sapGO_repId_cluFlag_LCAtaxId.tsv.gz | human proteins with GO terms and cluster | 2.1 MB |
| 4-domain-clustering.zip | domain-family analysis | 4.7 MB |
| 5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz | all 214M entries mapped (cluFlag: 1 AFDB50, 2 AFDB cluster, 3 fragment-removed, 4 singleton) | 1.6 GB |
| 6-all-vs-all-similarity-…-Dec16_2024.tsv.gz | cluster-representative all-vs-all E-values | 4.5 GB |
| 7-AFDB50-repId_memId.tsv.gz | AFDB50 sequence-cluster membership | 1.2 GB |
| **8-cluster_alignment_tmscore_lddt.tar.gz** | **member→representative alignments: fident, alnlen, mismatch, gapopen, coordinates, evalue, bits, lddt, alntmscore**; per-cluster mean LDDT/TM | 948.3 MB |
| commands.gz | commands used | 2.8 kB |

The v6 release (based on AFDB v6, updated 2026-03 to 2026-08) has files 1, 2, 5, 6 and 7, at 1.5 GB, 47.6 MB, 1.9 GB, 4.1 GB and 1.7 GB. No file 8 is listed for v6. The interactive site is https://cluster.foldseek.com (HTTP 200 on 2026-10-30).

**Implications:**
- File 8 gives within-cluster (`fident`, `lddt`, `alntmscore`) triplets for about 30M member–representative pairs at no compute cost. It can be used for quick exploratory curves.
- These pairs come in a **star topology** (members versus representative), are **structurally pre-filtered** and are **not phylogenetically independent**. Use them for exploration and stratification, not as the main test.
- The cluster LCA gives a crude age for each cluster, which can serve as a covariate.

---

## 2. Classic relationships between structural and sequence divergence

- **Chothia & Lesk 1986.** *The relation between the divergence of sequence and structure in proteins*, EMBO J 5:823–826. DOI 10.1002/j.1460-2075.1986.tb04288.x, PMID 3709526 [Abs]. Core RMSD Δ(Å) = **0.40·e^(1.87H)**, where H is the fraction of mutated residues. At about 20% identity, the common core can be under half of each molecule.
  - *[Ours]* With Poisson H = 1 − e^(−d), this becomes Δ(d) = 0.40·exp(1.87(1 − e^(−d))). The initial slope is about 0.75 Å per substitution per site, and the curve **saturates at about 2.6 Å**. "Exponential in identity" is thus roughly linear for small d and saturating for large d. It is consistent with Illergård's linear range but warns that the relation will flatten.
- **Wood & Pearson 1999.** *Evolution of protein sequences and structures*, J Mol Biol 291:977–995. DOI 10.1006/jmbi.1999.2972, PMID 10452901 [Abs]. In 36 families, DALI structural similarity is linear in Smith–Waterman sequence similarity (18/36 families with r ≥ 0.878). Higher-order terms improved the fit by less than 7% in 27/36. **"Structural mutation sensitivity" varies about fourfold** between families and is unrelated to structural class, mutation rate or size.
- **Koehl & Levitt 2002.** *Sequence variations within protein families are linearly related to structural variations*, J Mol Biol 323:551–562. DOI 10.1016/S0022-2836(02)00971-3, PMID 12381308 [Abs]. Twelve families. cRMS is linear in a sequence-space distance built from designed sequences. Related: Koehl & Levitt 2002 PNAS 99:1280 (PMID 11805293) on the sequence space allowed by a fold.
- **Panchenko & Madej 2005.** BMC Evol Biol 5:10. DOI 10.1186/1471-2148-5-10, PMID 15691378 [Abs]. Loop structural similarity is also linearly correlated with sequence similarity in most families. Loops are not random coils, and loop-based distances can give better trees. This argues against core-only analysis.
- **Grishin 2001.** *Fold change in evolution of protein structures*, J Struct Biol 134:167–185. DOI 10.1006/jsbi.2001.4335, PMID 11551177 [Abs]. Homologs can change fold through indels, circular permutation and β-sheet rearrangement. "Structure is conserved" has exceptions, so a nonlinear tail is expected.
- **Pascual-García et al. 2010.** *Quantifying the evolutionary divergence of protein structures: the role of function change and function conservation*, Proteins 78:181–196. DOI 10.1002/prot.22616, PMID 19830831 [Abs]. Globins, aldolases, P-loop and NADP-binding superfamilies. They introduce a **contact divergence** score that tracks sequence divergence better than earlier measures. Main findings:
  1. Small structural and sequence divergences are proportional, consistent with a clock.
  2. Function conservation strongly limits structural divergence.
  3. Structure diverges faster per unit of sequence for pairs with *different* function; they suggest positive selection.
  4. At low identity, structural diversity "explodes".
  5. Large indels almost always come with function change.

  We did not verify the exact formula; we believe it is a log-transformed, null-corrected contact overlap. Related: Pascual-García et al. 2009 PLoS Comput Biol 5:e1000331 (PMID 19325884), on discrete versus continuous structure space.
- **Sikosek & Chan 2014.** *Biophysics of protein evolution and evolutionary protein biophysics*, J R Soc Interface 11:20140419. DOI 10.1098/rsif.2014.0419, PMID 25165599 [Abs]. A review of marginal stability, robustness, epistasis and fold switching. Background on why structural response to substitutions is not uniform.
- **The Echave line of work, a mutational and neutral view of structural divergence:**
  - Echave & Fernández 2010, Proteins 78:173 (PMID 19731380): the linearly forced elastic network model (LFENM) reproduces structural divergence as the physical response to random mutations, so the conserved core matches mutationally robust regions.
  - Marcos & Echave 2020, Curr Res Struct Biol 2:156 (PMID 34235475): the site-by-site profile of structural divergence is *shaped by mutation and scaled by selection*.
  - Echave & Carpentier 2025, Proteins 94:1030 (PMID 41387199) and 2026, Mol Biol Evol 43(7) msag162 (PMID 42501309): in 34 enzyme families, mutational, stability and activity constraints each contribute, and their balance varies widely between families.
  - Echave & Wilke 2017, Annu Rev Biophys 46:85 (PMID 28301766): a review of biophysical models of sequence-rate variation.
- **Metric references:**
  - TM-score: Zhang & Skolnick 2004, Proteins 57:702 (PMID 15476259).
  - TM-align: Zhang & Skolnick 2005, NAR 33:2302 (PMID 15849316).
  - TM-score 0.5 as a same-fold threshold, P = 5.5×10⁻⁷: Xu & Zhang 2010, Bioinformatics 26:889 (PMID 20164152).
  - lDDT: Mariani et al. 2013, Bioinformatics 29:2722 (PMID 23986568).
  - Twilight zone: Rost 1999, Protein Eng 12:85 (PMID 10195279).

**Implications:**
1. Both linear and saturating forms are documented. Fit both, family by family, and report d*.
2. Expect 4–5-fold variation in slope between families (Wood & Pearson; Illergård).
3. Function change speeds up structural divergence (Pascual-García 2010). This gives a direct prediction for paralogs and for LDO/MDO pairs.
4. Loops and indels carry signal (Panchenko), so analyse both core and whole-domain.
5. RMSD needs core-only comparisons. TM-score and LDDT have hard bounds (TM for unrelated pairs is roughly 0.2–0.3; both are ≤ 1), so they saturate by construction and need a transform before any linearity test.

---

## 3. Molecular clock, neutral theory, and saturation for sequence and structure

**Sequence clock and neutral theory:**
- Zuckerkandl & Pauling 1965 ("Evolutionary divergence and convergence in proteins", in *Evolving Genes and Proteins*, Academic Press, pp. 97–166; book chapter, no PMID).
- Kimura 1968, *Evolutionary rate at the molecular level*, Nature 217:624 (DOI 10.1038/217624a0, PMID 5637732): under neutrality the substitution rate equals the mutation rate, which gives a clock.
- Ohta 1973, Nature 246:96 (PMID 4585855): nearly neutral theory, in which rates depend on population size and generation time, so the clock varies by lineage.
- Reviews: Bromham & Penny 2003, Nat Rev Genet 4:216 (PMID 12610526); Kumar 2005, Nat Rev Genet 6:654 (PMID 16136655).
- **Time-dependency of rates:** Ho et al. 2005, MBE 22:1561 (PMID 15814826); Ho et al. 2015, Mol Ecol 24:6007 (PMID 26769402). Short-timescale rates are inflated by segregating slightly deleterious polymorphisms.
- **Variation in rate between proteins:** mainly expression and misfolding avoidance (Drummond et al. 2005, PNAS 102:14338, PMID 16176987). Lineage-specific rate shifts covary across interacting genes (evolutionary rate covariation, ERC: Clark et al. 2012, Genome Res 22:714, PMID 22287101).

**Structural clock, and the only direct test we found:** Pascual-García A, Arenas M, Bastolla U. *The molecular clock in the evolution of protein structures.* Syst Biol 68:987–1002 (2019). DOI 10.1093/sysbio/syz022, PMID 31111152 [Abs]. Code at https://github.com/ugobas/Molecular_clock. Main points:
- A new "Significant Clock Violations" test applied to three superfamilies.
- Clock violations are **generally consistent between sequence and structure, but larger and more significant in structure**.
- Function change (GO, InterPro) is associated with large structural violations.
- About one third of significant violations appear only in structure.
- For close pairs, violations are significant in sequence but not in structure. They attribute this to time-dependent rates from slightly deleterious substitutions that do not affect structure.

**Saturation corrections:**
- **Amino acids.** p-distance → Poisson d = −ln(1−p) → JC-type for 20 states, d = −(19/20)·ln(1 − 20p/19) → Kimura's empirical protein distance d = −ln(1 − p − 0.2p²) (Kimura 1983, *The Neutral Theory of Molecular Evolution*, CUP; book) → gamma-corrected distances (Grishin 1995, J Mol Evol 41:675, PMID 18345592) → full ML distances (JTT/LG/WAG + Γ), as used by Illergård.
- **Structural similarity S ∈ [S₀, 1].** The generic correction is d = −ln((S − S₀)/(1 − S₀)) [Ours, standard form]. FoldTree applies a correction of this type (Supplementary Eqs 3–12). Its Fident distance goes to infinity as identity approaches about 7%, and the corrected LDDT and TM distances go to infinity as similarity approaches 0 [FT]. Pascual-García's contact divergence is of the same type [Abs; formula not verified].
- **3Di.** These are 20-state sequences, so the same Poisson/JC/ML machinery applies. ML needs a 3Di Q-matrix (§4).

**Implications and hypotheses:**
- **H-C1.** After correcting for saturation, the slope of structural divergence against time (TimeTree, eukaryotic orthologs) is constant within families. Test by comparing linear, saturating and piecewise fits, and with relative-rate (triplet) tests that use an outgroup species.
- **H-C2 (neutral prediction).** If structural change is mostly a passive response to neutral substitutions (Echave), then structural rate ∝ sequence rate across families, so the ratio is constant. Excess structural divergence given sequence divergence marks function change or selection (Pascual-García 2010, 2019).
- **H-C3.** At short times the sequence clock looks faster (time-dependency), but structure does not (Pascual-García 2019). This predicts a positive intercept for sequence distance against time and not for structure.
- **Pitfall.** Using sequence distance as the "time" axis and then claiming a structural clock is circular. Keep the axes separate: (A) structure vs sequence, (B) sequence vs time, (C) structure vs time.

---

## 4. Structure-based phylogenetics

- **Foldseek and 3Di.** van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M. *Fast and accurate protein structure search with Foldseek.* Nat Biotechnol 42:243–246 (2024; online 2023). DOI 10.1038/s41587-023-01773-0, PMID 37156916.
  - 3Di is a 20-letter alphabet describing tertiary interactions per residue.
  - Search is 4–5 orders of magnitude faster than earlier tools, with 86%, 88% and 133% of the sensitivity of Dali, TM-align and CE respectively [Abs].
  - Foldseek reports `fident`, `lddt`, `alntmscore` and `prob` per alignment.
- **FoldTree.** Moi D, Bernard C, Steinegger M, Nevers Y, Langleib M, Dessimoz C. *Structural phylogenetics unravels the evolutionary diversification of communication systems in gram-positive bacteria and their viruses.* Nat Struct Mol Biol 32:2492–2502 (2025). DOI 10.1038/s41594-025-01649-8, PMID 41073779, PMC12700828 [FT]. Bioinformatics preprint 2023: 10.1101/2023.09.19.558401.
  - **Method.** Foldseek all-vs-all, then distances from **Fident of 3Di+AA alignments, statistically corrected**, then QuickTree NJ, then MAD rooting.
  - **Benchmarks.** OMA HOG families (n = 4,592, relatively close) and CATH families (n = 488, more divergent). Scored by **taxonomic congruence score (TCS)**, ASTRAL quartet support and **root-to-tip variance** as a measure of clock adherence.
  - **Results:**
    1. FoldTree gets the highest share of top-TCS trees, especially in divergent families.
    2. Trees from **LDDT- or TM-based distances underperform**, which the authors attribute to conformational variation confounding rigid or global metrics.
    3. Using only 3Di for both alignment and distance did not help, so the gain comes mostly from better homologous-residue alignment.
    4. FoldTree trees have the **lowest root-to-tip variance**, i.e. are most clock-like, both before and after correction.
    5. Filtering out families with mean pLDDT < 60 or high pLDDT variance improves performance.
  - The authors themselves say that clock adherence on its own is a weak quality indicator.
- **3DiPhy.** Puente-Lelievre C, et al. (2023, bioRxiv). *Tertiary-interaction characters enable fast, model-based structural phylogenetics beyond the twilight zone.* DOI 10.1101/2023.12.12.571181. Turns the Foldseek 3Di scoring matrix into an ML substitution model. We did not confirm a peer-reviewed version; the BEAST2 blog (Dec 2025) still cites the preprint.
  - Review: Puente-Lelievre C, Malik AJ, Douglas J. *Protein structural phylogenetics.* Genome Biol Evol 17:evaf139 (2025). DOI 10.1093/gbe/evaf139, PMID 40839422 [FT]. It explicitly says "Do protein structures obey the molecular clock?" is unresolved: Pascual-García 2019 found more clock violations in structure, whereas Moi et al. found 3Di trees more ultrametric. "The question … depends on both the dataset and how 3D structures are interpreted."
  - Fullmer MS, Puente-Lelievre C, Matzke NJ (2025, bioRxiv 10.1101/2025.06.30.662300). Adding 3Di to amino-acid data improves resolution, but less so for short, α-helical proteins such as histones. In BEAST GammaSpike tests the 3Di data *accepted* punctuated (non-clock) rate changes while amino acids did not reject them, the opposite of what was expected [web summary of preprint; not read in full].
  - BEAST2 now has a 3Di package ("FoldBeast") with the Foldseek, Garg–Hochberg AF and LLM matrices [Live, beast2.org blog].
- **Q.3Di substitution matrices.** Garg SG, Hochberg GKA. *A general substitution matrix for structural phylogenetics.* Mol Biol Evol 42:msaf124 (2025). DOI 10.1093/molbev/msaf124, PMID 40476610, PMC12198762 [FT].
  - **Q.3Di.AF** was estimated with QMaker from 1,660 Foldseek clusters of AFDB-SwissProt (50–2,500 members each). **Q.3Di.LLM** was estimated from 6,653 Pfam families translated to 3Di by ProtT5.
  - The two matrices are practically identical. On held-out Pfam 3Di MSAs, 3DiPhy is preferred in only 278 alignments against 6,267 that prefer Q.3Di.*.
  - EF-Tu/EF-G 3Di trees root the tree of life between Archaea and Bacteria. In the amino-acid tree the branch between the paralogs is BL = 3.284, which they consider "probably saturated".
  - **Key caveat:** when sequences are modelled in two different conformations (KaiB, Rpc10), the 3Di trees **split by conformation, with long spurious branches**. Check conformational homogeneity before building trees.
  - Note: the Q.3Di name exists (Q.3Di.AF, Q.3Di.LLM). We found no model called "Q3Di" other than these.
- **FoldMason.** Gilchrist CLM, Mirdita M, Steinegger M. *Multiple protein structure alignment at scale with FoldMason.* Science 391:485–488 (2026). DOI 10.1126/science.ads6733, PMID 41610233 [Abs]. Progressive multiple structural alignment (MSTA) built on Foldseek and TM-align. Scales to hundreds of thousands of structures, is about 100× faster than earlier MSTA tools, and gives per-column confidence. It supported phylogenetics beyond the twilight zone in a flavivirus glycoprotein case study.
- **Other approaches.** Malik AJ, Poole AM, Allison JR. *Structural phylogenetics with confidence.* MBE 37:2711 (2020). DOI 10.1093/molbev/msaa100, PMID 32302382: structural trees are sensitive to length differences, and they propose bootstrap-like support from MD or Monte Carlo shape fluctuations. Also Structome (Malik et al. 2023, Bioinform Adv, PMID 38046099). LoL-align (Reifenrath et al. 2026, Nat Methods, DOI 10.1038/s41592-026-03249-w), a probabilistic aligner now inside Foldseek.

**Reported linearity and saturation, in short:**
- Structural distances "are less likely to saturate over time" (Moi et al., Discussion).
- Uncorrected distances "may saturate at a maximum finite distance" (Moi et al., Supp.).
- Amino-acid trees are "taller" in substitutions per site than 3Di trees (BEAST2 blog).
- Garg & Hochberg argue that 3Di overcomes saturation for paralogs that diverged near LUCA.
- **No study we found fits structural distance against calibrated time across many families.** This is a clear gap (§10).

**Implications:**
1. Use Fident(3Di+AA), corrected 3Di p-distance, ML 3Di distance (Q.3Di.AF) and LDDT/TM as **separate metrics**, and expect different saturation behaviour.
2. Use ML distances, IQ-TREE with Q.3Di.AF against LG/JTT+Γ, on FoldMason MSAs to get tree lengths for 3Di and AA on the same topology. The ratio of tree lengths is an alphabet-matched version of Illergård's rate ratio.
3. Filter by pLDDT, as FoldTree does (mean ≥ 60–70).
4. Flag conformational heterogeneity.

---

## 5. Orthology, paralogy, the ortholog conjecture, and LDO

- **Ortholog conjecture:**
  - Nehrt NL, Clark WT, Radivojac P, Hahn MW. PLoS Comput Biol 7:e1002073 (2011). DOI 10.1371/journal.pcbi.1002073, PMID 21695233 [Abs]. In human and mouse, **paralogs (especially within-species) were better predictors of function than orthologs**, against the conjecture. They argue that cellular context matters more than sequence.
  - Altenhoff AM, Studer RA, Robinson-Rechavi M, Dessimoz C. PLoS Comput Biol 8:e1002514 (2012). DOI 10.1371/journal.pcbi.1002514, PMID 22615551 [Abs]. After controlling for GO biases (authorship, term frequency, background, propagation), orthologs are **weakly but significantly** more functionally similar than paralogs, with only a weak decline as sequence divergence grows.
  - Chen X, Zhang J. PLoS Comput Biol 8:e1002784 (2012). PMID 23209392. GO cannot test the conjecture, but RNA-seq supports it.
  - Kryuchkova-Mostacci N, Robinson-Rechavi M. PLoS Comput Biol 12:e1005274 (2016). PMID 28030541. Tissue specificity diverges slowly between orthologs and quickly between paralogs.
  - Stamboulian M, Guerrero RF, Hahn MW, Radivojac P. Bioinformatics 36:i219 (2020). DOI 10.1093/bioinformatics/btaa468, PMID 32657391 [Abs]. Uses 40k+ experimentally annotated proteins in human/mouse and S. cerevisiae/S. pombe. **Strong evidence against the conjecture for function prediction**; questions about the evolution of function "remain difficult". Using all homologs maximises accuracy.
  - Reviews: Koonin 2005, Annu Rev Genet 39:309 (PMID 16285863); Gabaldón & Koonin 2013, Nat Rev Genet 14:360 (PMID 23552219).
  - Duplication fate: Lynch & Conery 2000, Science 290:1151 (PMID 11073452). Asymmetric divergence: Conant & Wagner 2003, Genome Res 13:2052 (PMID 12952876).
- **LDO conjecture, PMID 40813247 (title in the brief verified):** Warwick Vesztrocy A, Glover N, Thomas PD, Dessimoz C, Julca I. *Unveiling the functional fate of duplicated genes through expression profiling and structural analysis.* Genome Res 35:2273–2284 (2025). DOI 10.1101/gr.280166.124, PMC12487826 [FT].
  - **Claim:** after duplication, the copy with the shorter branch (the **least diverged ortholog**, LDO) keeps the ancestral function. The longer-branch copy (MDO) specialises.
  - **Data:** PANTHER v18, 143 species, 15,693 gene trees. Expected branch lengths come from species-tree branch lengths times a family-specific factor, which is a rate × time decomposition. Duplications are classed by significant rate shifts.
  - About 82% are normal–normal and about 13% normal–long. In total **86% symmetric and 14% asymmetric**. About 70% of duplications are on terminal branches; 68% are in Viridiplantae.
  - **Structure:** AFDB v4 models for 1,829,120 of 1,968,858 proteins; Foldseek 8.ef4e960 alignment LDDT; 926,814 paralog pairs; 857,663 paralog–outgroup comparisons.
  - **Design:** compare MDO with outgroup og1 and LDO with outgroup og2, with og2 chosen so that the sequence identities differ by at most 20 points; robustness checked from 5% to 50%. The strict matched set was only 1,319 of 548,818 (0.2%).
  - **Result:** asymmetric duplications have median LDDT 0.723 (MDO–og1) vs 0.766 (LDO–og2); symmetric ones 0.851 vs 0.856. Excess structural divergence of the MDO grows with branch-length asymmetry.
  - **Expression** (16 animals, 20 plants): the MDO has lower co-expression with orthologs and higher tissue specificity (τ).
  - **Caveat:** the authors note that "structural predictions often rely on sequence information", so structural similarity is confounded with sequence similarity. This is why they match on sequence.
- **Structural divergence of paralogs in the AF era:**
  - Dennler O, Ryan CJ. NAR Genom Bioinform 7:lqaf051 (2025). DOI 10.1093/nargab/lqaf051, PMID 40290317 [Abs]. In yeast and human paralogs, AF structural similarity and PLM similarity are **not redundant with sequence identity** for predicting shared function (protein–protein interactions, synthetic lethality).
  - Nakamura K, Siltberg-Liberles J. J Mol Evol 94:149 (2026). PMID 41606170. Calmodulin superfamily (448 proteins): structure and disorder track sequence phylogeny except in clades with recent duplications, and structural divergence is a "symptom" of functional divergence [Abs].
  - Zhang S, Zhang T, Fu Y. Genome Res (2023). DOI 10.1101/gr.277771.123 [Abs]. Human, E. coli and M. jannaschii proteomes: twilight-zone structural homolog pairs cover 8–32% of genes.
  - Beyond LDO, we found **no large-scale, time-calibrated, structure-based test of the ortholog conjecture.**

**Implications and hypotheses:**
- **H-O1.** At matched amino-acid distance, paralog pairs show greater structural divergence than ortholog pairs. Test with Ancestra HOGs and pair-level matching or regression with sequence distance as covariate.
- **H-O2 (structural LDO).** Within HOGs with in-paralogs, the longer-branch copy shows larger LDDT/TM/3Di divergence to the outgroup orthologs than the LDO, beyond what sequence explains.
- **H-O3.** Ortholog pairs give clean time points: their divergence time is the speciation time from TimeTree. Paralog divergence times are bracketed by the HOG duplication level.
- **Pitfalls.**
  1. Orthology errors put paralogs into "ortholog" sets.
  2. One-to-one ortholog filtering selects slowly evolving families.
  3. AF2 similarity is inflated for close homologs, and the LDO paper's sequence matching is the minimum control.
  4. Branch-length asymmetry is noisy for short branches.

---

## 6. Ancestra (evocellnet/ancestra-mvp-ai) [Live, README; not installed]

- MIT-licensed. Created 2026-10-28 as an AI-generated output of the same ETH hackathon, and labelled as such with the human disclaimer. **None of its results have been verified by a human.**
- **What it does:** reconstructs **hierarchical orthologous groups (HOGs)** level by level up a rooted guide species tree:
  1. Reference-free family step: DIAMOND clustering, connected components, then Leiden splitting of giant components.
  2. Similarity computed only within families.
  3. Bottom-up merge at each node by mutual best matches, absorption of in-paralogs and repair rules. Ambiguous many-to-many cases are resolved with small local trees (FAMSA + FastTree in `accurate`, NJ in `fast`).
  4. An optional synteny layer uses gene order.
- An optional **Foldseek 3Di "family rescue"** exists but is **off by default**; it lowered OrthoBench F by 5.0 points.
- **Outputs:** OrthoXML 0.5 (`hogs.orthoxml`, which loads in pyham), `pairs.tsv` (gene_a, gene_b, level) and a provenance `manifest.json`.
- **Inputs:** one FASTA per species, file names matching the leaves of a Newick guide tree. A helper builds the tree from the NCBI taxdump.
- **Reported performance:** Open OrthoBench test F of 76.2 (`accurate`) and 74.8 (`fast`), against OrthoFinder 3.2 at 69.7 and SonicParanoid2 at 70.4. The confidence intervals overlap, so it is not shown to be better than those two. A QfO 2022 run took about 30 GB RAM and a few hours for 78 proteomes and 989k proteins.
- **Implications.**
  1. HOG levels give duplication timing (the taxonomic level) and ortholog/paralog labels for every pair.
  2. Ancestra is unbenchmarked on deep prokaryotic and eukaryotic sampling at AFDB scale. Run it on a **species subset of UniProt reference proteomes**, and cross-check a sample against OMA (Altenhoff et al. 2024, NAR 52:D513, PMID 37962356) or PANTHER.
  3. Keep the structure rescue off, so that orthology is defined independently of structure. Otherwise the analysis becomes circular.

---

## 7. Candidate case studies (brief)

**Domain architecture:**
- **TED.** Lau AM, Bordin N, Kandathil SM, Sillitoe I, Waman VP, Wells J, Orengo CA, Jones DT. *Exploring structural diversity across the protein universe with The Encyclopedia of Domains.* Science 386:eadq4946 (2024). DOI 10.1126/science.adq4946, PMID 39480926 [Abs]. Nearly 365M domains in AFDB v4, more than 100M beyond sequence-based methods, covering more than 1M taxa. 77% of non-redundant domains resemble known CATH superfamilies. More than 10,000 new superfamily–superfamily interactions.
  - Data on Zenodo 10.5281/zenodo.13908086 (CC-BY 4.0) [Live]. Includes `ted_365m.domain_summary.cath.globularity.taxid.tsv.gz` (19.9 GB) and `ted_214m_per_chain_segmentation.tsv.gz` (9.6 GB), among others.
- Classic domain-architecture evolution: Björklund AK, Ekman D, Light S, Frey-Skött J, Elofsson A. J Mol Biol 353:911 (2005), PMID 16198373. Indels at the N and C termini dominate, domain exchange is rare, and a "domain distance" agrees with sequence similarity. Moore AD et al. Trends Biochem Sci 33:444 (2008), PMID 18656364.
- CATH: Sillitoe et al. 2021, NAR 49:D266 (PMID 33237325). Pfam: Mistry et al. 2021, NAR 49:D412 (PMID 33125078).
- **Implication.** AFDB clusters need 90% coverage, so **changes in domain architecture fall between clusters**. Use TED segmentation per chain on HOG members to score gain, loss and rearrangement of domains, and test whether they coincide with accelerated structural divergence (Pascual-García 2010: large indels go with function change).

**Kinases:**
- Hanks-type (ePK) kinases are monophyletic and present in all three domains (Stancik IA et al. J Mol Biol 430:27, 2018, PMID 29138003). ePK-like kinases (ELKs) are as common as histidine kinases in ocean metagenomes: 27,677 ePKs and 18,699 ELKs in 20 families, and all but one of 10 key conserved residues is lost in some family (Kannan N, Taylor SS, Zhai Y, Venter JC, Manning G. PLoS Biol 5:e17, 2007, PMID 17355172).
- Other sources:
  - "Universal core" of the protein kinase-like fold and a combined sequence+structure tree (Scheeff & Bourne 2005, PLoS Comput Biol 1:e49, PMID 16244704).
  - Bacterial and archaeal ePK-like kinases (Leonard, Aravind & Koonin 1998, Genome Res 8:1038, PMID 9799791).
  - Human kinome (Manning et al. 2002, Science 298:1912, PMID 12471243).
  - Choanoflagellate tyrosine kinases (Manning et al. 2008, PNAS 105:9674, PMID 18621719).
  - Pseudokinases across the tree of life (Kwon et al. 2019, Sci Signal, PMID 31015289).
  - KinOrtho (Huang et al. 2021, BMC Bioinformatics, PMID 34537014).
  - Evolution of the kinase switch (Taylor, Meharena & Kornev 2019, IUBMB Life 71:672, PMID 31059206).
- Pkinase (PF00069) is the most frequent Pfam among the AFDB-cluster domain families (Barrio-Hernandez 2023, Fig. 4g).
- **Hypothesis.** Structural divergence of the kinase core is accelerated on the branches where tyrosine kinases and ePK regulatory elements (activation loop, αC–β4) arose, compared with the deep ELK split.

**Cytoskeleton:**
- FtsZ is a tubulin homolog (Löwe & Amos 1998, Nature 391:203, PMID 9428770). MreB is an actin homolog (van den Ent, Amos & Löwe 2001, Nature 413:39, PMID 11544518).
- Prokaryotic cytoskeletons are compositionally plastic, and most eukaryotic complexity predates the last eukaryotic common ancestor (LECA) (Wickstead & Gull 2011, J Cell Biol 194:513, PMID 21859859). Kinesin evolution: Wickstead, Gull & Richards 2010 (PMID 20423470).
- Asgard archaea: Zaremba-Niedzwiedzka et al. 2017, Nature 541:353 (PMID 28077874); Asgard profilins regulate actin (Akıl & Robinson 2018, Nature 562:439, PMID 30283132); Odinarchaeota tubulin (Akıl et al. 2022, Sci Adv 8:eabm2225, PMID 35333570); Heimdallarchaeial ancestry of eukaryotes (Eme et al. 2023, Nature 618:992, PMID 37316666).
- Independently, Daumiller, Li & Bryant 2026 (§9) report that the most structurally conserved eukaryotic "Strict Core" is dominated by cytoskeletal proteins.
- **Hypothesis.** Actin and tubulin are very slow in both sequence and structure within eukaryotes, but the prokaryote-to-eukaryote transition (MreB/FtsZ → actin/tubulin) is a burst of structural change. Concretely: a large structural-distance jump on the stem branch relative to sequence change.

---

## 8. Divergence-time resources (needed for any "linear with time" claim)

- **TimeTree 5.** Kumar S, Suleski M, Craig JM, Kasprowicz AE, Sanderford M, Li M, Stecher G, Hedges SB. MBE 39:msac174 (2022). DOI 10.1093/molbev/msac174, PMID 35932227 [FT].
  - 137,306 species from 4,075 studies. Article licence CC-BY 4.0. Data downloadable "for research and teaching".
  - **API (tested [Live]):** `https://api.timetree.org/pairwise/{taxidA}/{taxidB}[/csv|/json]` and `/mrca/id/{id+id+…}/json`. Inputs are NCBI taxIDs.
  - Human–mouse: 83.5 Mya (CI 80.1–86.9; 75 studies). Human–S. cerevisiae: 1,218 Mya (CI 1,113–1,323; adjusted 1,164; 17 studies). E. coli–B. subtilis: 1,909 Mya (CI 659–3,158; 2 studies). Human–E. coli: **empty**. Strain-level taxIDs (E. coli K-12 MG1655, 511145) are **missing**, so map to species.
  - A global Newick from 2015 is directly accessible (`timetree.org/public/data/TimetreeOfLife2015.nwk`, 50,646 leaves). TTOL5 subtrees can be downloaded through the GUI.
- **NCBI Taxonomy.** Schoch CL et al. Database 2020:baaa062. PMID 32761142. Gives topology and ranks, no times. The Ancestra guide tree builder uses it. Public domain.
- **GTDB.** Parks DH et al. NAR 50:D785 (2022), PMID 34520557. GTDB-Tk v2 (PMID 36218463). **Latest release R11-RS232 (15 Apr 2026): 901,341 genomes in 199,923 species clusters; bac120 and ar53 trees** [Live]. Branch lengths are in substitutions, not time. Licence: CC BY-SA 4.0 according to the GTDB site (not re-verified).
- **Other options:** species-tree branch lengths in substitutions, as the LDO paper did with the PANTHER species tree; relative time from HOG levels (Ancestra/OMA), i.e. an ordinal age; cluster LCA (Barrio-Hernandez); Hug et al. 2016 tree of life (Nat Microbiol 1:16048, PMID 27572647).
- **Pitfalls.**
  1. Prokaryote divergence times are poorly constrained, and HGT breaks gene/species concordance.
  2. Deep eukaryotic nodes carry ±10–20% CIs.
  3. Pairwise times are not independent: many pairs share the same node. Analyse at the level of nodes, or use phylogenetic mixed models.
  4. Using sequence-based branch lengths as "time" when testing a sequence-versus-structure relationship is circular.
- **Recommended design** [Ours]:
  - Primary time test on **eukaryotic one-to-one orthologs** with TimeTree node ages (roughly 10–1,600 Mya).
  - Prokaryotic and deep analyses only against sequence distance (Illergård-style), or ordinal HOG level.
  - For ortholog pairs, divergence time is the TimeTree MRCA age of the two species.

---

## 9. Recent (2023–2026) work using AF/AFDB at scale

- **Daumiller D, Li Q, Bryant P.** *The structural history of Eukarya.* Sci Adv 12:eaeh2524 (2026, published 2026-10-28). DOI 10.1126/sciadv.aeh2524, PMID 42950045 [Abs]; preprint 10.64898/2026.02.04.703859 [web summary].
  - 1,542 eukaryotic proteomes from AFDB, about 300 trillion Foldseek alignments, scored by a "structural bit score" (bits × geometric mean of LDDT and TM).
  - Proteome distance trees built by NJ.
  - **Findings:** a "Strict Core" dominated by cytoskeleton; structural acceleration in birds and ants that is unrelated to proteome size (ρ = 0.04).
  - **No time calibration and no linearity test**, according to the preprint summary. This is the closest competitor; our time-calibrated, ortholog-level analysis would complement it.
- **Lau et al. 2024 (TED)** — see §7.
- Durairaj J, … Schwede T, Pereira J. *Uncovering new families and folds in the natural protein universe.* Nature 622:646 (2023). DOI 10.1038/s41586-023-06622-3, PMID 37704037 [Abs]. Sequence-similarity network of AFDB90v4 (uniprot3d.org), the β-flower fold, the TumE–TumA toxin–antitoxin system.
- Bordin N et al. *AlphaFold2 reveals commonalities and novelties in protein structure space for 21 model organisms.* Commun Biol 6:160 (2023). DOI 10.1038/s42003-023-04488-9, PMID 36755055.
- Zhang, Zhang & Fu 2023, Genome Res (§5). Warwick Vesztrocy et al. 2025 (§5). Garg & Hochberg 2025 (§4). Moi et al. 2025 (§4).
- Bastida A, Muñoz AM, Egea-Cortines M (2026, bioRxiv 10.64898/2026.04.29.721599) [web summary]. Secondary-structure distance (S2D) in 96 MADS-box proteins: the relation between S2D and sequence distance is **highly dispersed, not linear**, and depends on the domain (the K domain is variable).
- Echave & Carpentier 2025 and 2026 (§2).
- **AFDB resources:** Varadi et al. 2022, NAR 50:D439 (PMID 34791371); Varadi et al. 2024, NAR 52:D368 (PMID 37933859); Bertoni et al. 2026, NAR (DOI 10.1093/nar/gkaf1226, PMID 41273079).
  - **AFDB v6:** UniProt 2025_03, **241,070,489 structures**, adds isoforms and **per-entry A3M MSAs plus `msa_depths.csv`** [Live, CHANGELOG]. Since v4: 65.7M added and 39.3M removed.
  - The AFDB is ageing relative to UniProt: 631 of 20,504 human entries conflict with UniProt 2025_03 (Tsitsa et al. 2025, NSMB 32:2374, PMID 41372615).
- **AlphaFold caveats:**
  - Jumper et al. 2021, Nature 596:583 (PMID 34265844).
  - AF2 does not capture the structural effects of missense mutations (Buel & Walters 2022, NSMB 29:1, PMID 35046575; Pak et al. 2023, PLoS One 18:e0282689, PMID 36928239).
  - AF2 fails on fold switching (Chakravarty & Porter 2022, Protein Sci 31:e4353, PMID 35634782).
  - Predictions are hypotheses, not experimental structures (Terwilliger et al. 2024, Nat Methods 21:110, PMID 38036854).
- **Gene-age pitfalls:** lineage-specific genes are often artefacts of failing to detect homology (Weisman, Murray & Eddy 2020, PLoS Biol 18:e3000862, PMID 33137085); phylostratigraphic bias (Moyers & Zhang 2015, MBE 32:258, PMID 25312911). GenEra adds Foldseek to gene dating (Barrera-Redondo et al. 2023, Genome Biol 24:54, PMID 36964572). These bear directly on the 3.91% species-specific AFDB clusters.

**What we could not find:** a structure-based molecular clock calibrated against TimeTree across thousands of families, or an AFDB-wide comparison of the relative rates of 3Di and amino-acid evolution. Our searches covered PubMed, bioRxiv via web search, and Google-like web search. It is possible we missed work that is not indexed under these terms.

---

## 10. Open gaps we can address

1. **Is there a structural clock against real time?** Fit structural divergence (corrected TM, LDDT, 3Di-ML, Fident) against TimeTree ages for eukaryotic one-to-one orthologs across thousands of HOGs. Compare linear, saturating and piecewise models. Estimate per-family rates and their dispersion (overdispersion of the clock) for sequence and structure side by side.
2. **An alphabet-matched rate ratio.** Re-estimate Illergård's "3–10×" with 20-state 3Di against 20-state amino acids under matched ML models (Q.3Di.AF vs LG/JTT+Γ) on shared topologies, at AFDB scale and stratified by fold class (β > α?), length, pLDDT and cluster LCA.
3. **Saturation and dynamic range.** Locate the sequence distance d* at which each structural metric departs from linearity. Test whether 3Di extends the usable range for phylogenetic signal beyond amino acids (FoldTree and Garg–Hochberg claim this, but did not quantify it against time).
4. **Neutral-response vs selection decomposition.** Treat residual structural divergence given sequence divergence as an index of selection or function change. Test whether it is higher for paralogs than for orthologs, for MDO than for LDO, for function-changing pairs (GO/EC discordance) and for pairs that changed domain architecture.
5. **Lineage-specific structural acceleration.** Relative-rate tests on ortholog triplets to find branches with structural speed-up beyond sequence change. Candidate examples: kinases at the origin of tyrosine kinases or in specific bacterial phyla, and the prokaryote–eukaryote transition in the cytoskeleton. Compare with Daumiller et al.'s acceleration in birds and ants.
6. **Evolutionary audit of AFDB clusters.** Check how cluster LCAs and species-specific clusters hold up when orthology (Ancestra) and homology-detection failure are taken into account. Quantify how much cluster splitting is driven by changes in domain architecture.
7. **Calibrating for AF2 artefacts.** Compare PDB–PDB pairs with AF–AF pairs of the same proteins, and include MSA depth from AFDB v6 as a covariate. Estimate how much "structure is more conserved" is a prediction artefact.

---

## 11. Datasets and resources

| Name | URL | Version (as of 2026-10-30) | Size | Licence |
|---|---|---|---|---|
| AFDB clusters (Barrio-Hernandez 2023) | https://afdb-cluster.steineggerlab.workers.dev ; https://cluster.foldseek.com | v3 (paper, AFDB "UniProt v3"/v4 IDs); v6 (2026-03, AFDB v6) | v3: 0.26–4.5 GB per file (~9.5 GB total); v6: ~9.3 GB | CC-BY 4.0 |
| AFDB cluster alignments (file 8) | …/afdb-cluster/v3/8-cluster_alignment_tmscore_lddt.tar.gz | v3 (added 2025-09-12) | 948 MB | CC-BY 4.0 |
| AlphaFold DB | https://alphafold.ebi.ac.uk ; ftp.ebi.ac.uk/pub/databases/alphafold ; gs://public-datasets-deepmind-alphafold | v6 (UniProt 2025_03), 241,070,489 models, with MSAs | sequences.fasta 110 GB; full models tens of TB (use Foldseek DBs / proteome tars) | CC-BY 4.0 |
| Foldseek / Foldseek cluster | https://github.com/steineggerlab/foldseek | Paper used 4.645b789; LDO used 8.ef4e960; check latest | — | GPLv3 |
| FoldMason | https://github.com/steineggerlab/foldmason | Science 2026 | — | GPLv3 (verify) |
| FoldTree | https://github.com/DessimozLab/fold_tree | NSMB 2025 | — | verify on repo |
| Q.3Di.AF / Q.3Di.LLM matrices | supplementary data of Garg & Hochberg 2025 (MBE msaf124) | 2025 | small | per journal/supplement (verify) |
| TED domains | https://ted.cathdb.info ; Zenodo 10.5281/zenodo.13908086 | AFDB v4, 2024-10-31 | ~70 GB total; key files 2.9–19.9 GB | CC-BY 4.0 |
| CATH | https://www.cathdb.info | v4.3/4.4 (verify latest) | — | CC-BY 4.0 (verify) |
| Pfam / InterPro | https://www.ebi.ac.uk/interpro | Pfam 37.x (paper used 34.0/35.0) | — | CC0 |
| TimeTree 5 | https://timetree.org ; API https://api.timetree.org | TToL5 (2022) | API; 2015 global nwk 50,646 tips | article CC-BY 4.0; data "for research and teaching" (verify terms) |
| NCBI Taxonomy | https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | daily | ~60 MB | public domain |
| GTDB | https://gtdb.ecogenomic.org | R11-RS232 (2026-04-15) | trees < 50 MB | CC BY-SA 4.0 (verify) |
| Ancestra | https://github.com/evocellnet/ancestra-mvp-ai | freeze-v1 (2026-10-28) | — | MIT |
| OMA (orthology cross-check) | https://omabrowser.org | All.May2026 referenced by Ancestra README | — | verify (CC-BY 4.0 expected) |
| PANTHER (LDO data) | https://pantherdb.org | v18 used by LDO paper | — | verify |
| UniProt reference proteomes | https://www.uniprot.org | 2025_03 matches AFDB v6 | — | CC-BY 4.0 |
| SCOP / ECOD (benchmarks) | https://scop.mrc-lmb.cam.ac.uk ; http://prodata.swmed.edu/ecod | ECOD F99 20230309 used by the cluster paper | — | verify |
| Molecular_clock (clock-violation test) | https://github.com/ugobas/Molecular_clock | Syst Biol 2019 | — | verify |

---

## References (verified PMID/DOI; grouped by section)

*Anchors.* Illergård, Ardell & Elofsson 2009 Proteins 77:499 (PMID 19507241, 10.1002/prot.22458). Illergård 2010 PhD thesis, Stockholm University (urn:nbn:se:su:diva-35872). Barrio-Hernandez et al. 2023 Nature 622:637 (PMID 37704730, 10.1038/s41586-023-06510-w).

*Sequence–structure.* Chothia & Lesk 1986 (PMID 3709526). Wood & Pearson 1999 (PMID 10452901). Koehl & Levitt 2002 JMB (PMID 12381308) and PNAS (PMID 11805293). Panchenko & Madej 2005 (PMID 15691378). Grishin 2001 (PMID 11551177). Pascual-García et al. 2009 (PMID 19325884) and 2010 (PMID 19830831). Sikosek & Chan 2014 (PMID 25165599). Echave & Fernández 2010 (PMID 19731380). Marcos & Echave 2015 (PMID 25922797) and 2020 (PMID 34235475). Echave & Wilke 2017 (PMID 28301766). Echave & Carpentier 2025 (PMID 41387199) and 2026 (PMID 42501309). Zhang & Skolnick 2004 (PMID 15476259) and 2005 (PMID 15849316). Xu & Zhang 2010 (PMID 20164152). Mariani et al. 2013 (PMID 23986568). Rost 1999 (PMID 10195279).

*Clock.* Kimura 1968 (PMID 5637732). Ohta 1973 (PMID 4585855). Bromham & Penny 2003 (PMID 12610526). Kumar 2005 (PMID 16136655). Ho et al. 2005 (PMID 15814826) and 2015 (PMID 26769402). Drummond et al. 2005 (PMID 16176987). Clark et al. 2012 (PMID 22287101). Grishin 1995 (PMID 18345592). Pascual-García, Arenas & Bastolla 2019 (PMID 31111152).

*Structural phylogenetics.* van Kempen et al. 2024 (PMID 37156916). Moi et al. 2025 (PMID 41073779). Garg & Hochberg 2025 (PMID 40476610). Gilchrist, Mirdita & Steinegger 2026 (PMID 41610233). Puente-Lelievre, Malik & Douglas 2025 (PMID 40839422). Puente-Lelievre et al. 2023 bioRxiv 10.1101/2023.12.12.571181. Fullmer et al. 2025 bioRxiv 10.1101/2025.06.30.662300. Malik, Poole & Allison 2020 (PMID 32302382). Malik et al. 2023 (PMID 38046099). Reifenrath et al. 2026 (10.1038/s41592-026-03249-w). Steinegger & Söding 2018 (PMID 29959318).

*Orthology.* Nehrt et al. 2011 (PMID 21695233). Altenhoff et al. 2012 (PMID 22615551). Chen & Zhang 2012 (PMID 23209392). Kryuchkova-Mostacci & Robinson-Rechavi 2016 (PMID 28030541). Stamboulian et al. 2020 (PMID 32657391). Warwick Vesztrocy et al. 2025 (PMID 40813247). Dennler & Ryan 2025 (PMID 40290317). Nakamura & Siltberg-Liberles 2026 (PMID 41606170). Koonin 2005 (PMID 16285863). Gabaldón & Koonin 2013 (PMID 23552219). Lynch & Conery 2000 (PMID 11073452). Conant & Wagner 2003 (PMID 12952876). Altenhoff et al. 2024 OMA (PMID 37962356). Thomas et al. 2022 PANTHER (PMID 34717010).

*Case studies.* Lau et al. 2024 TED (PMID 39480926). Björklund et al. 2005 (PMID 16198373). Moore et al. 2008 (PMID 18656364). Sillitoe et al. 2021 (PMID 33237325). Mistry et al. 2021 (PMID 33125078). Kannan et al. 2007 (PMID 17355172). Scheeff & Bourne 2005 (PMID 16244704). Leonard et al. 1998 (PMID 9799791). Manning et al. 2002 (PMID 12471243) and 2008 (PMID 18621719). Stancik et al. 2018 (PMID 29138003). Kwon et al. 2019 (PMID 31015289). Huang et al. 2021 (PMID 34537014). Taylor et al. 2019 (PMID 31059206). Löwe & Amos 1998 (PMID 9428770). van den Ent et al. 2001 (PMID 11544518). Wickstead & Gull 2011 (PMID 21859859). Wickstead et al. 2010 (PMID 20423470). Zaremba-Niedzwiedzka et al. 2017 (PMID 28077874). Akıl & Robinson 2018 (PMID 30283132). Akıl et al. 2022 (PMID 35333570). Eme et al. 2023 (PMID 37316666).

*Time and taxonomy.* Kumar et al. 2022 TimeTree 5 (PMID 35932227). Schoch et al. 2020 (PMID 32761142). Parks et al. 2022 (PMID 34520557). Chaumeil et al. 2022 (PMID 36218463). Hug et al. 2016 (PMID 27572647).

*AF-scale and AF caveats.* Daumiller, Li & Bryant 2026 (PMID 42950045). Durairaj et al. 2023 (PMID 37704037). Bordin et al. 2023 (PMID 36755055). Zhang, Zhang & Fu 2023 Genome Res (10.1101/gr.277771.123). Bastida et al. 2026 bioRxiv (10.64898/2026.04.29.721599). Varadi et al. 2022 (PMID 34791371) and 2024 (PMID 37933859). Bertoni et al. 2026 (PMID 41273079). Tsitsa et al. 2025 (PMID 41372615). Jumper et al. 2021 (PMID 34265844). Buel & Walters 2022 (PMID 35046575). Pak et al. 2023 (PMID 36928239). Chakravarty & Porter 2022 (PMID 35634782). Terwilliger et al. 2024 (PMID 38036854). Weisman et al. 2020 (PMID 33137085). Moyers & Zhang 2015 (PMID 25312911). Barrera-Redondo et al. 2023 (PMID 36964572).
