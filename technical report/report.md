(Human-written text starts)
This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects (LINK). All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication.

[Project-specific text] The goal of this specific project was to investigate protein 3d structure evolution. It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents, and finished within 48 hours on a computational cluster with 48 cpus, 500 RAM, and a GPU.

Currently none of the results have been verified by a human. We plan to verify the agent's methodology and benchmarking, and we will update this section accordingly.

The rest of this repository is entirely AI generated. It includes code for implementation of the research and a AI-written report in the style of a journal article, detailing the final method, the path followed to reach it, and benchmarking against other methods.
(Human-written text ends)

---


# Is protein structure a molecular clock? Revisiting structural divergence with predicted structures, orthology and divergence times

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


## 1 Introduction

**Main point.** Many papers say that protein structure is "more conserved than sequence" and that it "evolves linearly". These statements come from small sets of experimental structures, and they measure evolution in amino-acid substitutions, not in years. The AlphaFold Database now holds predicted structures for almost every known protein, and ortholog families can now be built across hundreds of proteomes. This makes it possible to ask, at scale, whether structural change behaves like a molecular clock. We ask that question, and we also ask what the AlphaFold Database clusters tell us about evolution.

### 1.1 Background: clocks, sequences and structures

A **molecular clock** is the idea that molecular change accumulates at a roughly constant rate, so that the amount of change between two species is proportional to the time since they split (Zuckerkandl & Pauling 1965). Neutral theory explains why this can happen. If most substitutions do not affect fitness, the substitution rate equals the mutation rate (Kimura 1968). Real clocks are imperfect: rates vary between lineages (Ohta 1973; Bromham & Penny 2003).

Protein structure changes more slowly than sequence. Chothia & Lesk (1986) showed that structural deviation of the protein core rises as sequence identity falls. Wood & Pearson (1999) and Koehl & Levitt (2002) found near-linear relations between structural and sequence divergence within families. The reference point for this project is **Illergård, Ardell & Elofsson (2009)**. They compared pairs of SCOP domains and found that core structural change grows linearly with the maximum-likelihood amino-acid distance. They also found that structure is "three to ten times more conserved than sequence", and that families differ about fivefold in their "structural response".

Three points about that study shape our design. All three come from the literature review in LITERATURE.md.

1. **The axis was substitutions, not time.** Illergård et al. plotted structural change against the number of amino-acid substitutions per site. Their "linear" result therefore says that structure changes in proportion to sequence change. It does not say that structure changes in proportion to geological time. Testing a structural clock against dated species splits is a new question, not a replication.
2. **The "3–10×" figure depends on the alphabet.** With raw structural states, structure looked 15–30× more conserved. The 3–10× figure appeared only after the amino acids were grouped into as few states as the structural descriptors had (Illergård 2010, PhD thesis). A fair rate comparison needs two alphabets of the same size, scored with the same kind of model.
3. **The evidence on structural clocks is mixed.** Pascual-García et al. (2019) found more and larger clock violations in structure than in sequence, linked to changes in function. Moi et al. (2025) found that trees built from structure-informed distances were the most clock-like of the methods they compared. Fullmer et al. (2025, preprint) found the opposite pattern in histones. A recent review calls the question unresolved (Puente-Lelievre, Malik & Douglas 2025).

### 1.2 New resources and the gaps they leave

**AlphaFold Database (AFDB) clusters.** Barrio-Hernandez et al. (2023) clustered about 214 million predicted structures from the AFDB. They first grouped sequences at 50% identity, then grouped the resulting representatives by structural similarity with Foldseek. This gave 2.3 million non-singleton structural clusters. Their evolutionary analysis was limited to assigning each cluster the last common ancestor (LCA) of its member species. They did not relate sequence divergence to structural divergence, did not use divergence times, and did not separate orthologs from paralogs. Their discussion says that structure "is conserved for longer periods of evolutionary time than protein sequences", but this was not tested.

**A matched alphabet.** Foldseek describes each residue's local 3D environment with one of 20 letters, called **3Di** (van Kempen et al. 2024). Because 3Di has 20 states, like the amino-acid alphabet, the two can be compared like for like. Garg & Hochberg (2025) estimated a substitution model for 3Di (Q.3Di.AF). With it, a 3Di distance can be computed with the same maximum-likelihood machinery as an amino-acid distance. This gives an alphabet-matched version of Illergård's ratio that was not available in 2009.

**What remained unknown.** Our literature search found no study that:

- fits structural divergence against calibrated divergence times across thousands of families;
- compares 3Di and amino-acid rates under matched models at AFDB scale;
- or tests whether the structure-to-sequence ratio stays constant across lineages and across depths of the tree.

The closest work, Daumiller, Li & Bryant (2026), compared 1,542 eukaryotic proteomes by structural similarity, but without time calibration or a linearity test.

### 1.3 Three pitfalls

**Selection bias inside clusters.** A pair of proteins taken from inside an AFDB cluster has already passed a structural-similarity filter. The filter required alignment coverage of 90% of both proteins and an E-value below 0.01 (Barrio-Hernandez et al. 2023). Very diverged pairs are therefore missing from such a sample, and this alone can make structural divergence look as if it levels off. A fair test must choose protein pairs without looking at their structures, for example as **orthologs** defined from sequences.

**Prediction artefacts.** AlphaFold2 is largely insensitive to point mutations (Buel & Walters 2022; Pak et al. 2023). It also draws on overlapping sequence alignments and templates for related proteins (Jumper et al. 2021). Predictions for two close orthologs may therefore be more similar than the real structures are. Any conclusion is first a conclusion about *predicted* structures, and the size of this bias has to be measured.

**Shared alignment error.** If the sequence distance and the structural distance come from the same alignment, an alignment error affects both and can create a correlation between them (noted for Illergård et al. 2009 in LITERATURE.md §1a). We therefore computed every distance under three different alignments (§3.6).

### 1.4 The two goals

**Goal 1: Is protein structure a molecular clock?** We revisit Illergård et al. (2009) at AFDB scale, using eukaryotic ortholog families chosen without regard to structure. We ask four things:

- **Shape.** How does structural divergence grow with amino-acid divergence: linearly, or does it level off (saturate)?
- **Ratio.** How much slower is 3Di change than amino-acid change, when both are measured with matched 20-state models?
- **Time.** Is that ratio the same for young and old branches of the tree of life?
- **Lineages.** Do some lineages change structure faster than their sequence change predicts?

Together these test a **neutral structural clock** null hypothesis: that structural change is a passive by-product of sequence change, with one constant ratio per family.

**Goal 2: AFDB clusters in the light of evolution.** This goal is more open. We chose three items:

- **Duplication.** Do paralogs diverge more in structure than orthologs at the same sequence distance? Does the faster-evolving copy after a duplication diverge more in structure? This follows the least-diverged-ortholog idea of Warwick Vesztrocy et al. (2025), building on the ortholog conjecture debate (Nehrt et al. 2011; Altenhoff et al. 2012). It also serves as the positive control for Goal 1.
- **Cluster audit.** How do orthology groups map onto AFDB structural clusters, and how often does domain architecture change?
- **Case study.** The cytoskeletal actin and tubulin folds.

---

## 2 How the project was run

**Main point.** An AI "PI" agent ran the project with AI subagents, under rules meant to stop an autonomous system from fooling itself:

- independent design review before any data analysis;
- a pre-registered analysis plan, with tuning and confirmation on separate data;
- a public log of every change to the plan;
- verification on simulated data with known answers before trusting any estimator.

**Agents.** A PI agent (Claude Opus 5.5) planned the work, made decisions and kept the research record. It delegated tasks to subagents (Claude Opus 5.5 or Sonnet 5.5 [CHECK: model mix per agent]). Each subagent returned a written report. 17 such reports were in `agents/` when this draft was written. The tasks were literature, data scouting, Ancestra set-up, the species panel, the distance library, the pair pipeline, the statistics engine, noise floors, HOG construction, integration, bias fixes, deep node ages, the AFDB-cluster baseline, the PDB control (two reports), the cluster audit and the case study. The PI read the summaries rather than raw logs, and recorded milestones, decisions and status in a research tree (RESEARCH_TREE.md).

**Timeline and compute.**
- The project started on 2026-10-30 at 09:22 on one cluster node (48 CPUs, 500 GB RAM, one A100 GPU).
- The first node had a 12-hour limit. Work was checkpointed and moved to a second node after 10.9 h of compute, within a total cap of 48 h of compute.
- Long jobs had completion markers and written resume rules (NODE_MOVE.md).

**Design by independent critique.** The PI wrote a design (DESIGN.md) from the literature and data-scouting reports. Independent reviewer agents then critiqued it in three rounds:

- **Round 1:** five reviewers (an external PI, molecular evolution, statistics, algorithms and compute, and AI failure modes).
- **Round 2:** three reviewers (external PI, statistics, molecular evolution).
- **Round 3:** one consistency reviewer. This reviewer turned the design into an executable analysis plan (ANALYSIS_PLAN.md).

Each reviewer either proposed concrete changes within scope or said none were needed. The version 1 → 2 change log alone lists 17 changes, each tied to a review item. The reviews replaced pair-based summaries with branch lengths on a species tree, added pre-specified decision rules, added a positive control, and cut Goal 2 down to three items.

**Pre-registration.** Before any result was computed, the 7,000 sampled gene families were split at random into a **tuning set (20%, 1,400 families)** and a **confirmation set (80%, 5,600 families)**. The split file's checksum is recorded in RESULTS_LEDGER.md and is never regenerated.
- Every adjustable setting is chosen on the tuning set only. This includes rate-variation parameters, noise floors, error calibration and exclusion thresholds.
- The settings are then frozen, and the confirmation set is analysed once.
- The confirmation script refuses to run unless the plan's checksum matches the ledger, the code tree is clean and an explicit flag is set.

**Deviation log.** Every change to the frozen plan goes into DEVIATIONS.md with its reason and effect. When this draft was written, the log had 26 entries (28 numbered items). All of them were made before any confirmatory endpoint was computed. Most were justified by synthetic tests or by pipeline checks on a 20-species test run. Three were substantive:

- **D22:** the primary alignment changed (§3.6).
- **D31 and D32:** two estimator biases were found in simulation and corrected (§3.10).

Where a fix replaced an original criterion, the original result is reported too.

**Verification.** Each estimator was tested on simulated data with a known answer before it was used on real data (§3.10). An agent that had not produced a result is assigned to re-derive each confirmatory number (check V11).

**Process lessons.** LESSONS.md keeps short process rules with evidence, such as "run the whole pipeline on a small real dataset before freezing" (this caught the D22 problem, which no reviewer had predicted). It is pruned at milestones and reviewed by a fresh agent.

---

## 3 Data and methods

**Main point.** We built ortholog families from 180 eukaryotic proteomes using sequences only. We sampled 7,000 families and computed sequence, 3Di and 3D distances for 1.4 million protein pairs from AlphaFold structures. We then estimated two things: how structure grows with sequence divergence, and how the structure-to-sequence ratio varies with time and lineage. The whole pipeline is summarised in Figure 1.

![Figure 1](figures/fig1_workflow.png)

**Figure 1. Workflow.** The species panel and predicted structures (1–2) feed orthology inference (3). Ortholog families are filtered and sampled without looking at structure (4), then split once into tuning and confirmation sets (5). All pairs within each family (6) are aligned three ways (7) and scored for amino-acid, 3Di and 3D divergence (8). Noise floors (9) come from random and identical-sequence pairs. Two estimation engines (10a, 10b) feed six pre-registered endpoints (11). Counts are from the project's filter logs.

### 3.1 Terms used throughout

| Term | Meaning |
|---|---|
| **Ortholog / paralog** | Two genes are **orthologs** if they descend from one ancestral gene separated by a speciation, e.g. the same gene in human and mouse. They are **paralogs** if they were separated by a gene duplication. An **in-paralog** duplicated after the speciation in question; an **out-paralog** before it. |
| **HOG** | A hierarchical orthologous group: all the genes descending from one ancestral gene at a given node of the species tree. Here, "family" means a top-level HOG. |
| **AA** | Amino acid (protein sequence). |
| **pLDDT** | AlphaFold's per-residue confidence score (0–100). Residues with pLDDT ≥ 70 are usually modelled reliably. |
| **3Di** | Foldseek's 20-letter structural alphabet. Each residue gets one letter describing its local 3D interaction with a neighbouring residue (van Kempen et al. 2024). |
| **TM-score** | A global structural similarity score between 0 and 1 (Zhang & Skolnick 2004). Unrelated proteins typically score about 0.2–0.3. |
| **LDDT** | Local distance difference test: the fraction of inter-residue distances that are preserved between two structures, between 0 and 1 (Mariani et al. 2013). |
| **RMSD** | Root-mean-square deviation (Å) between aligned atoms after optimal superposition. |
| **Substitution model** | A table of how often each letter is replaced by each other letter over evolutionary time. We use LG for amino acids and Q.3Di.AF for 3Di (Garg & Hochberg 2025). [CHECK: add LG reference; not in LITERATURE.md] |
| **Γ (gamma) rate variation** | Different sites in a protein evolve at different speeds. A gamma distribution with shape α models this: small α means strong variation. |
| **ML distance** | The maximum-likelihood estimate of the expected number of substitutions per site separating two sequences, given a substitution model. Unlike the raw fraction of differing sites (the p-distance), it corrects for multiple hits at the same site. |
| **Saturation** | When so many changes have occurred that further change no longer increases the observed difference. The distance then becomes uninformative. |
| **Branch length** | The amount of change along one branch of a tree, here in substitutions per site. |
| **Divergence time** | The age of the split between two species, in millions of years ago (Mya), here from TimeTree (Kumar et al. 2022). |

### 3.2 Species panel

**We chose 180 eukaryotes that span the eukaryotic tree and whose proteins nearly all have an AlphaFold model with an identical sequence** (Figure 2a).

- **Source:** UniProt release 2026_03 reference proteomes, one canonical protein per gene; 3,460,494 proteins in total.
- **Coverage gate:** a species was eligible only if at least 80% of its proteins had an AFDB version 4 model with exactly the same sequence (checked by sequence checksum).
- **Completeness gates:** BUSCO genome completeness thresholds by group (relaxed for protists, whose BUSCO scores are systematically low).
- **Dating:** each species had to be datable in TimeTree.
- **Selection order:**
  1. the eukaryotes of the Quest for Orthologs 2022 benchmark set that passed the gates;
  2. 33 sister-species pairs, giving very young splits;
  3. one high-quality species per phylum;
  4. a greedy fill aimed at an even spread of split ages on a log scale.
- **Composition:** 66 animals, 40 fungi, 20 land plants and 54 species from other groups, covering all major eukaryotic supergroups.

Some model organisms could not be included because their current proteomes no longer match the AFDB v4 sequences, among them rat, zebrafish, chicken and *Xenopus*. A separate panel of 100 prokaryotes (60 bacteria, 40 archaea) was built for the cytoskeleton case study only. Its only Asgard archaeon is *Prometheoarchaeum syntrophicum* MK-D1.

**Divergence times.** Times came from TimeTree (Kumar et al. 2022): the TimeTree-of-Life 2015 tree supplied the topology, and the TimeTree pairwise service (version 20260820) supplied the ages.
- 15,755 of the 16,110 species pairs were dated.
- The undated pairs mostly involve three metamonads and *Fonticula*.
- After pruning to the panel, the tree has 166 internal nodes. Internodes shorter than 5 million years were collapsed into **polytomies** (nodes with more than two children), leaving 123 internal nodes.
- 117 internal nodes are dated; 6 small nodes are undated.
- Ages were made monotone from root to tips (34 fixes, all deep nodes).
- The root is dated at 1,455 Mya.
- Check: human–mouse = 83.5 Mya (95% interval 80.1–86.9), which passes the plan's check of 87 ± 10 Mya.

A first version also collapsed nodes whose age intervals overlapped their parent's. This merged well-supported clades such as Opisthokonta into the root, so it was dropped before any family was built (deviation D30).

Pairs are dominated by deep splits (Figure 2c): 52% of dated pairs fall at the root age. Analyses therefore work with branches of the tree, not with pairs. TimeTree data cannot be redistributed, so per-pair ages are kept out of committed files. Only binned counts are shown.

![Figure 2](figures/fig2_panel.png)

**Figure 2. Species panel.** (a) Number of species per group (180 eukaryotes). (b) Ages of the 117 dated internal nodes of the collapsed panel tree. (c) Ages of all dated species pairs (log scale on both axes): most pairs share a few deep splits. (d) Gene families per age stratum, where age is that of the clade each family was sampled from. Grey bars are eligible families; coloured bars are sampled families, split into confirmation (80%) and tuning (20%). The youngest stratum is small (224 families).

### 3.3 Predicted structures

**Every analysed protein has an AlphaFold v4 model whose sequence is identical to the 2026 UniProt sequence.**

- Of the 3.46 M panel proteins, 3,290,300 (95.1%) have an AFDB v4 model, and 94.9% of all proteins match exactly. Models were fetched from the compressed FoldComp archive of AFDB v4 with zero failures.
- Proteins longer than 2,700 residues have no AFDB model (16,160 proteins), so very long proteins such as titins and dyneins are absent.
- For each structure we kept the coordinates, the per-residue pLDDT and the 3Di string computed by Foldseek (version 10.941cd33).
- AFDB v4 is the release that the AFDB cluster paper used (cluster release v3).

### 3.4 Orthology with Ancestra

**Orthologs and paralogs were inferred from sequences alone**, so that orthology is independent of the structures we later compare.

- **Tool:** Ancestra (commit ffb0bbe, `accurate` preset; github.com/evocellnet/ancestra-mvp-ai). Ancestra is a graph-based method that builds HOGs level by level up a rooted guide tree, here the NCBI taxonomy (Schoch et al. 2020).
  - It finds families with DIAMOND searches.
  - It merges genes at each node by mutual best matches.
  - It resolves ambiguous cases with small local gene trees.
  - Its optional structure-based "rescue" step was left off.
- **Provenance:** Ancestra is itself an AI-generated output of the same hackathon and has not been verified by a human. We therefore checked it.
  - Against Ensembl, it recovered 148 of 149 human–mouse one-to-one orthologs, and none of its 187 sampled human–mouse ortholog calls was an Ensembl paralog.
  - Its *timing* of duplications agreed with Ensembl for only about 36% of human paralog pairs. We therefore never use Ancestra's duplication ages (decision D2). Paralog analyses instead match on sequence distance and use explicit outgroups.
- **Run:** one run on all 180 proteomes took 5.0 h wall time (84 CPU-hours, 161 GB peak memory) and produced 109,598 top-level HOGs.
- **Pair labels:** every pair of genes within a HOG is labelled ortholog or paralog by the deepest group that contains both.

### 3.5 Family sampling and the tuning/confirmation split

**We sampled families by sequence information only, keeping one gene per species and avoiding families whose sequences are too diverged to measure.**

**Gene filters** (counts are genes in top-level HOGs):

| Step | Genes remaining | Families remaining |
|---|---|---|
| Input | 2,924,444 | 109,598 |
| Model sequence identical to UniProt | 2,761,752 | |
| Not a fragment | 2,624,408 | |
| Mean pLDDT ≥ 70 | 1,632,375 | 57,348 |
| Length within 0.5–2× the family median | 1,521,329 | |
| ≥ 8 species with exactly one copy | | **12,336** |

The pLDDT filter removes 38% of genes and is the largest cut. The analysis is therefore about well-predicted, mostly ordered proteins. This is a stated limitation.

**One-to-one sets.** Within each family we keep only species with exactly one gene, so that each species contributes a single ortholog. This does not guarantee orthology: 4% of pairs in these sets are labelled out-paralogs by Ancestra, a sign of hidden paralogy or of differential gene loss. They are kept and flagged, and removed in a sensitivity analysis.

**Choosing species within a family.** For each family we searched the species tree for the largest clade in which:
- at least 8 species are present;
- every pair of genes is at most 2 amino-acid substitutions per site apart.

The distance used here comes from a sequence-only alignment with fixed settings (α = 1, never re-tuned), so selection cannot depend on structure or on tuned parameters. The 2-substitution limit keeps every pair in the range where ML distances are reliable. A family keeps at most 30 species, chosen to cover as many tree nodes as possible. 8,211 families qualified.

**Sampling.** Families were assigned to strata by the age of their selected clade: < 100, 100–300, 300–1,000 and ≥ 1,000 Mya. We sampled 7,000 families with equal quotas per stratum, and unused quota moved to the other strata.
- The youngest stratum had only 224 eligible families, so the shallow end of the time axis is thin (Figure 2d).
- Final sample: 224, 756, 5,064 and 956 families.

**Split.** Within strata, families were split at random (seed 20261030) into **1,400 tuning** and **5,600 confirmation** families.

**Pairs.**
- **Primary set:** all pairs within each family, 1,087,063 pairs.
- **Paralog-control set (Goal 2, endpoint E6):** 36,500 duplication triplets in 12,808 paralog groups from 5,001 families, each made of two paralogs in one species plus orthologs in outgroup species. These give 312,674 further pairs.
- **Cap:** each triplet has at most 5 outgroup species per side, and each paralog group at most 5 triplets (C25). Without the cap these pairs alone would have exceeded the budget of 3.5 M pairs.

### 3.6 Divergence measures and the three alignments

**Each pair is aligned three ways, and from each alignment we compute amino-acid, 3Di and 3D divergence over the reliably predicted core.**

**Alignments.**

| Alignment | What it is | Role |
|---|---|---|
| **3Di+AA** | Foldseek's alignment, which scores 3Di and amino-acid letters together | Primary |
| **Sequence-only** | Global alignment with BLOSUM62 (Biopython), which ignores structure | Sensitivity |
| **TM-align** | Foldseek's TM-align mode, which superimposes coordinates and ignores letters | Sensitivity |

A conclusion counts as **alignment-robust** only if all three alignments agree on its sign and decision.

**Why the primary changed (deviation D22).** The design first named TM-align as primary, because it uses coordinates only and so favours neither alphabet. A full pipeline run on a 20-species test set, made before any endpoint was computed, showed three problems.

1. TM-align misaligned residues between close relatives. For 5,295 human–mouse orthologs, the median core amino-acid p-distance was 0.091 under TM-align vs 0.064 under the other two alignments, and median LDDT was 0.83 vs 0.93.
2. 8 of 104 pairs of *identical* sequences received a non-zero amino-acid distance (up to 3.2).
3. TM-align left more pairs with too little coverage: 10.7% of ortholog pairs and 24.7% of paralog pairs, vs 3.9% and 9.8% for 3Di+AA.

An error of this kind adds distance at short range and would bias the shape estimate. The 3Di+AA alignment uses both alphabets, so its bias on the 3Di-to-AA ratio should be more balanced than that of the sequence-only alignment, which favours amino acids. TM-align stays as a sensitivity alignment.

**Core.** Distances use only aligned columns where both residues have pLDDT ≥ 70. A pair is "ok" if it has at least 50 such columns and its alignment covers at least half of the shorter protein. Failed pairs are counted and given zero weight. They are never silently dropped.

**Measures per pair.**
- **AA distance:** ML distance under LG+Γ4.
- **3Di distance:** ML distance under Q.3Di.AF+Γ4. Both are computed on a GPU from 20×20 counts of aligned letter pairs, with a standard error (SE) from the likelihood curvature. The Γ shapes were tuned on the tuning set from the sequence-only alignment: α_AA = 1.64, α_3Di = 1.60 (C29). On the 3Di+AA alignment, the 3Di shape ran to the edge of the search grid, an artefact of column placement that is optimised for 3Di.
- **Global LDDT** (normalised by the shorter protein), **TM-score** and **core RMSD.** LDDT and TM are similarities bounded between a chance level and 1, so they saturate by construction. LDDT is therefore converted to a distance with a floor correction: D = −ln[(S − S₀)/(1 − S₀)], where S₀ is the chance similarity of unrelated proteins of that length (§3.7). This form follows the corrections used by FoldTree (Moi et al. 2025).
- **SE calibration.** On the tuning set, a block bootstrap over alignment columns was used to calibrate the SEs. Bootstrap SDs were 1.0–1.6 times the model SEs. AA and 3Di errors were correlated at 0.08–0.33 for the same pair, because both come from one alignment.

### 3.7 Noise floors

**Two reference sets tell us what "no evolution" and "no relationship" look like.**

**Chance level (random pairs).** 100,000 pairs of proteins were drawn so that the two proteins:
- come from different species;
- sit in different AFDB clusters;
- share no Pfam family or clan (Mistry et al. 2021).

Their similarity gives the chance level S₀ in six length bins. Under TM-align, random pairs have median TM-scores of 0.22–0.33, as expected for unrelated proteins.

**Zero divergence (identical sequences).** 17,768 pairs of different AFDB entries with *identical* sequences (almost all from different species) show how much two AlphaFold predictions of the same sequence differ.
- They are not identical: median TM-score 0.950, and only 8.8% of pairs reach TM ≥ 0.999.
- Their median 3Di distance is 0.052 substitutions per site. It rises with protein length, from 0.037 to 0.066.
- Their median 1 − LDDT is 0.020.

This prediction noise η is a floor under every structural distance. It is included as an intercept in the shape model and subtracted from the terminal branches (the branches leading to single species) in the branch model.

### 3.8 Two estimation engines

**We use two complementary estimators: one fits the shape of the structure-versus-sequence curve from pairwise distances, the other places each family's distances onto the dated species tree.**

**Pairwise engine (shape; endpoints E1, E2, E6).**
- **Data.** Within each family, each pair gives one point: x = AA distance, y = 3Di distance (or floor-corrected LDDT distance).
- **Model.** We fit y = η + a_h · x^c.
  - η is the prediction-noise floor.
  - a_h is a family-specific scale.
  - c is a single **curvature exponent** shared by all families.
  - c = 1 means structure grows in proportion to sequence (linear). c < 1 means structural change slows relative to sequence change (concave, i.e. saturating).
- **Why not branch increments.** Saturation shows up in total divergence, not in increments, so the shape is estimated on whole pairwise distances.
- **Errors-in-variables problem.** Both x and y are measured with error. Error in x flattens a fitted curve and makes a truly linear relation look concave.
- **Correction.** We correct for this with SIMEX (simulation–extrapolation) [CHECK: add SIMEX reference; not in LITERATURE.md]:
  1. Add extra, known amounts of noise to the data.
  2. Watch how the estimate drifts as noise increases.
  3. Extrapolate back to zero noise.
  4. Remove the bias that remains after SIMEX with a parametric bootstrap (D31, §3.10).
- **Uncertainty.** A bootstrap resamples families and species.

**Branch engine (ratio, time, lineages; endpoints E3, E4, E5).**
- **Idea.** The pairwise distances within a family are fitted as sums of branch lengths on the fixed, dated species tree, pruned to that family's species. Think of a tree drawn with unknown branch lengths. Each pairwise distance is the length of the path between two tips. With many pairs, the individual branch lengths can be solved for by weighted least squares.
- **Output.** Doing this separately for AA and for 3Di gives, for every branch of every family, an AA length and a 3Di length on the same branch.
- **Uses.**
  - Ratios of these lengths on young and old branches show whether the 3Di-to-AA ratio changes with time.
  - Ratios on particular lineages show whether some lineages change structure faster than sequence.
  - Fitting on the shared species tree also handles the main non-independence among pairs: many pairs share the same deep branches.
- **Safeguards.**
  - Ratios are always ratios of sums over many branches, never ratios of single short, noisy branches.
  - Families whose distances fit their own sequence-based tree more than three times better than the species tree are excluded from these endpoints. This is a guard against hidden paralogy and gene-tree discordance.

[Figure 3 to come: schematic of the two engines, showing a pairwise curve with c = 1 vs c < 1, and a family tree with AA and 3Di branch lengths.]

### 3.9 Endpoints and the pre-registered interpretation

**Six endpoints were fixed before analysis.** Each has a decision rule. Five of them (E1, E2, E4, E5, E6) are tested together with a Holm correction at an overall 5% error rate [CHECK: add Holm reference]. CI-based rules use 98% intervals. E3 is reported as an estimate.

| Endpoint | Question in plain words | Decision rule |
|---|---|---|
| **E1** | Is 3Di divergence linear in AA divergence (c ≈ 1), or does it saturate (c < 1)? | "Linear" if the 98% CI of c lies within [0.9, 1.1]; "concave" if the upper bound is below 1; otherwise inconclusive |
| **E2** | The same question for floor-corrected LDDT distance; also the AA distance d* at which raw 1 − LDDT departs from a straight line by more than 10% | As E1 |
| **E3** | How many times faster does sequence change than 3Di? k = total AA branch length / total 3Di branch length, its spread across families, and the α-helical vs β-sheet contrast | Estimate with 95% CI; compared with Illergård's 3–10× |
| **E4** | Does the 3Di/AA ratio change with branch age? Slope of log ratio vs log age in a mixed model over family × age-bin cells | Slope = 0 tested two-sided |
| **E5** | Do lineages differ in their 3Di/AA ratio beyond noise? Is the sister-lineage slope of 3Di vs AA change equal to 1? | "Neutral clock (N0) accepted" only if the slope's 98% CI lies within [0.8, 1.25], the excess lineage SD is below 0.10, E6 detects the positive control, and power is ≥ 80% |
| **E6** | Positive control: do paralogs diverge more in 3Di than orthologs at matched AA distance, and does the more-diverged copy after duplication diverge more in 3Di? | One-sided > 0 |

**E6 design.** E6 adapts the least-diverged-ortholog design of Warwick Vesztrocy et al. (2025) to avoid regression to the mean.
- **More-diverged copy.** The more-diverged copy (MDO) is chosen using one set of outgroup species and one half of the alignment columns (alternating 20-residue blocks).
- **Excess.** Its structural excess is then measured with a different outgroup set on the other half of the columns.
- **Direction check.** The direction of the published result should replicate.

**Interpretation table, in plain words.** The table was written before any result and is used only if the E5 spike-in power is at least 80%.

| If E1/E2 say | and E4 says | and E5 says | we conclude |
|---|---|---|---|
| linear | no age trend | N0 margins met | Structure follows a rescaled sequence clock. We call this "proportional", not "neutral": a good fit does not prove neutrality. |
| concave | ratio of AA to 3Di rises with age | N0 met | Structure saturates, or AA distances are under-corrected. Structure is clock-like only up to d*. |
| linear | AA/3Di rises with age, robust to α and to alternative deep ages | any | An epoch effect: deep branches are structurally more conservative. |
| any | AA/3Di falls with age | any | AA under-correction or early structural bursts; recheck with AA distance ≤ 1. |
| any | any | slope < 0.8 or excess lineage variance, robust to guards | Lineage-specific structural rates. N0 is rejected. |
| any | any | excess not robust to guards | An AlphaFold-quality or orthology artefact. |

The design explicitly removed the reading "a rising ratio means constraint", because constraint slows both alphabets. If E6 is positive and N0 holds, the reading is that duplication, not ortholog lineage, changes how much structure changes per sequence change.

**Known limitation, stated before the results.** Simulations showed that the slope's 98% CI in E5 is about ±0.24–0.28 wide even when N0 is exactly true. *Accepting* N0 within [0.8, 1.25] is therefore practically unreachable. The margins were not widened. E5 can still reject N0, but otherwise it will be reported as inconclusive.

### 3.10 Controls and verification

**Every estimator was checked on simulated data with a known answer, and the main sources of bias have explicit controls.**

**Known-answer tests.**
- **Distance library (33 automated tests).** On data simulated with AliSim, ML distances were unbiased within 1.6% for distances up to 2 [CHECK: add AliSim/IQ-TREE reference]. They matched IQ-TREE's own pairwise distances within 0.004%.
- **Pair pipeline.**
  - Every input pair yields exactly one output row.
  - A protein compared with itself scores TM = LDDT = 1.
  - Identical sequences get AA distance 0 (except under TM-align, see D22).
  - The in-house LDDT matches Foldseek's within 0.01.
- **Branch engine (V2).** The engine recovers known branch lengths without bias. Model-based SEs gave 96% coverage of 95% intervals, whereas the planned robust SEs gave only 85–88% (D2).
- **Curvature (V3) and the D31 fix.**
  - In simulations with the real tuning-set designs and calibrated noise, the planned SIMEX estimator returned c ≈ 0.95 under a truly linear relation. This would have produced a false "concave" call.
  - With the D31 bias correction, c = 0.993 and 0 of 50 linear datasets were called concave.
  - True curvature was recovered (c = 0.8 → 0.788, detected in all datasets).
  - The plan's rule still applies: if more than 10% of linear simulations are called concave on the real designs, E1/E2 are declared inconclusive.
- **Age slope (E4) and the D32 fix.** In simulated worlds with no age effect, the planned cell construction gave a spurious slope (+0.058, false-positive rate 20%). This came from noise on young terminal branches. A minimum-information rule for merging cells (D32) reduced this to −0.0075 with 0% false positives and recovered true effects (−0.1 → −0.100).
- **End-to-end (V4).** Synthetic families are simulated under N0 and passed through the whole pipeline. The endpoints should recover the truth, and a lineage "spike" should be detected.

**Controls on real data.**
- **Positive control:** E6 (duplication).
- **Negative control for E4:** branch midpoint ages are permuted among a family's branches, so no age signal should remain (D13).
- **Power check for E5:** a known 1.5-fold 3Di acceleration is spiked into 10% of lineages in the tuning families. E5's ability to detect it is measured over 50 replicates.
- **Sensitivity analyses:** the other two alignments; α × 1.5 and ÷ 1.5; AA distance ≤ 1; removal of hidden-paralogy pairs; single-copy families only; internal branches only; and alternative deep-node ages from Strassert et al. (2021) [CHECK: identifier from agents/deep_ages_report.md, not in LITERATURE.md].

**Bias checks on structures and samples (exploratory).**
- **AlphaFold artefact.** Ortholog pairs with experimental X-ray structures on both sides, released after AlphaFold2's training cutoff, compare AF–AF divergence with PDB–PDB divergence on the same residues.
- **Selection bias.** Pairs from inside AFDB clusters (about 27.7 M member–representative alignments released with the clusters) are compared with our ortholog pairs to show the effect of structural pre-selection.

**Goal 2 analyses (exploratory).** These have unadjusted CIs and are labelled exploratory.
- **Cluster audit:** maps every HOG onto AFDB v3 clusters and models why HOGs split across clusters. Predictors are length ratio, AA distance, pLDDT, and change in TED domain architecture (Lau et al. 2024).
- **Cytoskeleton case study:** estimates AA and 3Di branch lengths for actin- and tubulin-fold proteins from eukaryotes and the prokaryote panel. It uses fixed reference trees and FoldMason structural alignments (Gilchrist, Mirdita & Steinegger 2026).

**Reproducibility.**
- Seeds derive from 20261030.
- Every filter logs counts in and out.
- Every reported number gets a line in RESULTS_LEDGER.md (script, commit, input checksum).
- Licence-restricted TimeTree data never enter committed files.

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

## References

Identifiers are as verified in LITERATURE.md, except where marked [CHECK].

- Altenhoff AM, Studer RA, Robinson-Rechavi M, Dessimoz C (2012). Resolving the ortholog conjecture: orthologs tend to be weakly, but significantly, more similar in function than paralogs. *PLoS Comput Biol* 8:e1002514. doi:10.1371/journal.pcbi.1002514. PMID 22615551.
- Ancestra (2026). evocellnet/ancestra-mvp-ai, commit ffb0bbe (AI-generated hackathon software; not peer reviewed). https://github.com/evocellnet/ancestra-mvp-ai
- Barrio-Hernandez I, Yeo J, Jänes J, Mirdita M, Gilchrist CLM, Wein T, Varadi M, Velankar S, Beltrao P, Steinegger M (2023). Clustering predicted structures at the scale of the known protein universe. *Nature* 622:637–645. doi:10.1038/s41586-023-06510-w. PMID 37704730.
- Bromham L, Penny D (2003). The modern molecular clock. *Nat Rev Genet* 4:216–224. PMID 12610526. [CHECK title/pages]
- Buel GR, Walters KJ (2022). Can AlphaFold2 predict the impact of missense mutations on structure? *Nat Struct Mol Biol* 29:1–2. PMID 35046575.
- Chothia C, Lesk AM (1986). The relation between the divergence of sequence and structure in proteins. *EMBO J* 5:823–826. doi:10.1002/j.1460-2075.1986.tb04288.x. PMID 3709526.
- Daumiller D, Li Q, Bryant P (2026). The structural history of Eukarya. *Sci Adv* 12:eaeh2524. doi:10.1126/sciadv.aeh2524. PMID 42950045.
- Fullmer MS, Puente-Lelievre C, Matzke NJ (2025). bioRxiv preprint. doi:10.1101/2025.06.30.662300. [CHECK title]
- Garg SG, Hochberg GKA (2025). A general substitution matrix for structural phylogenetics. *Mol Biol Evol* 42:msaf124. doi:10.1093/molbev/msaf124. PMID 40476610.
- Gilchrist CLM, Mirdita M, Steinegger M (2026). Multiple protein structure alignment at scale with FoldMason. *Science* 391:485–488. doi:10.1126/science.ads6733. PMID 41610233.
- Illergård K (2010). *On the effects of structure and function on protein evolution.* PhD thesis, Stockholm University. urn:nbn:se:su:diva-35872.
- Illergård K, Ardell DH, Elofsson A (2009). Structure is three to ten times more conserved than sequence—a study of structural response in protein cores. *Proteins* 77:499–508. doi:10.1002/prot.22458. PMID 19507241.
- Jumper J, et al. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature* 596:583–589. PMID 34265844.
- Kimura M (1968). Evolutionary rate at the molecular level. *Nature* 217:624–626. doi:10.1038/217624a0. PMID 5637732.
- Koehl P, Levitt M (2002). Sequence variations within protein families are linearly related to structural variations. *J Mol Biol* 323:551–562. doi:10.1016/S0022-2836(02)00971-3. PMID 12381308.
- Kumar S, Suleski M, Craig JM, Kasprowicz AE, Sanderford M, Li M, Stecher G, Hedges SB (2022). TimeTree 5: an expanded resource for species divergence times. *Mol Biol Evol* 39:msac174. doi:10.1093/molbev/msac174. PMID 35932227.
- Lau AM, Bordin N, Kandathil SM, Sillitoe I, Waman VP, Wells J, Orengo CA, Jones DT (2024). Exploring structural diversity across the protein universe with The Encyclopedia of Domains. *Science* 386:eadq4946. doi:10.1126/science.adq4946. PMID 39480926.
- Mariani V, Biasini M, Barbato A, Schwede T (2013). lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests. *Bioinformatics* 29:2722–2728. PMID 23986568.
- Mistry J, et al. (2021). Pfam: the protein families database in 2021. *Nucleic Acids Res* 49:D412–D419. PMID 33125078.
- Moi D, Bernard C, Steinegger M, Nevers Y, Langleib M, Dessimoz C (2025). Structural phylogenetics unravels the evolutionary diversification of communication systems in gram-positive bacteria and their viruses. *Nat Struct Mol Biol* 32:2492–2502. doi:10.1038/s41594-025-01649-8. PMID 41073779.
- Nehrt NL, Clark WT, Radivojac P, Hahn MW (2011). Testing the ortholog conjecture with comparative functional genomic data from mammals. *PLoS Comput Biol* 7:e1002073. doi:10.1371/journal.pcbi.1002073. PMID 21695233.
- Ohta T (1973). Slightly deleterious mutant substitutions in evolution. *Nature* 246:96–98. PMID 4585855.
- Pak MA, et al. (2023). Using AlphaFold to predict the impact of single mutations on protein stability and function. *PLoS One* 18:e0282689. PMID 36928239.
- Pascual-García A, Arenas M, Bastolla U (2019). The molecular clock in the evolution of protein structures. *Syst Biol* 68:987–1002. doi:10.1093/sysbio/syz022. PMID 31111152.
- Puente-Lelievre C, Malik AJ, Douglas J (2025). Protein structural phylogenetics. *Genome Biol Evol* 17:evaf139. doi:10.1093/gbe/evaf139. PMID 40839422.
- Schoch CL, et al. (2020). NCBI Taxonomy: a comprehensive update on curation, resources and tools. *Database* 2020:baaa062. PMID 32761142.
- Strassert JFH, Irisarri I, Williams TA, Burki F (2021). A molecular timescale for eukaryote evolution with implications for the origin of red algal-derived plastids. *Nat Commun* 12:1879. doi:10.1038/s41467-021-22044-z. PMID 33767194. [CHECK: verified in agents/deep_ages_report.md, not in LITERATURE.md]
- van Kempen M, Kim SS, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Söding J, Steinegger M (2024). Fast and accurate protein structure search with Foldseek. *Nat Biotechnol* 42:243–246. doi:10.1038/s41587-023-01773-0. PMID 37156916.
- Warwick Vesztrocy A, Glover N, Thomas PD, Dessimoz C, Julca I (2025). Unveiling the functional fate of duplicated genes through expression profiling and structural analysis. *Genome Res* 35:2273–2284. doi:10.1101/gr.280166.124. PMID 40813247.
- Wood TC, Pearson WR (1999). Evolution of protein sequences and structures. *J Mol Biol* 291:977–995. doi:10.1006/jmbi.1999.2972. PMID 10452901.
- Zhang Y, Skolnick J (2004). Scoring function for automated assessment of protein structure template quality. *Proteins* 57:702–710. PMID 15476259.
- Zuckerkandl E, Pauling L (1965). Evolutionary divergence and convergence in proteins. In *Evolving Genes and Proteins*, Academic Press, pp. 97–166.

[References still to add with verified identifiers: LG model; gamma rate heterogeneity; SIMEX; Holm correction; AliSim/IQ-TREE; Biopython; BUSCO; Betts et al. 2018 (doi from agents/deep_ages_report.md).]


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
