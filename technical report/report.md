(Human-written text starts)
This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects (LINK). All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication.

[Project-specific text] The goal of this specific project was to investigate protein 3d structure evolution. It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents, and finished within 48 hours on a computational cluster with 48 cpus, 500 RAM, and a GPU.

Currently none of the results have been verified by a human. We plan to verify the agent's methodology and benchmarking, and we will update this section accordingly.

The rest of this repository is entirely AI generated. It includes code for implementation of the research and a AI-written report in the style of a journal article, detailing the final method, the path followed to reach it, and benchmarking against other methods.
(Human-written text ends)

---

> **DRAFT, first half (Introduction, project conduct, Data and methods).** Written before the confirmatory results exist. Results, Discussion and Limitations are placeholders. Numbers marked [CHECK] need confirmation against the project files before release.

# Is protein structure a molecular clock? Revisiting structural divergence with predicted structures, orthology and divergence times

## Abstract

[to be written after results]

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

[PLACEHOLDER: to be written after the confirmation run.]
- 4.1 Response shape (E1, E2, d*); alignment robustness. [Figure 4 to come: 3Di and LDDT distance vs AA distance with fitted curves.]
- 4.2 Conservation ratio k (E3) vs Illergård's 3–10×; fold-class contrast.
- 4.3 Ratio vs node age (E4); negative control; deep-age sensitivity. [Figure 5 to come.]
- 4.4 Lineage effects and the N0 test (E5); spike-in power.
- 4.5 Prediction artefact (AF vs PDB, exploratory) and the selection effect of in-cluster pairs (exploratory). [Figure 6 to come.]

## 5 Results: Goal 2, AFDB clusters in the light of evolution

[PLACEHOLDER]
- 5.1 Duplication positive control (E6) and the structural least-diverged-ortholog test.
- 5.2 Cluster audit: HOGs vs AFDB clusters, domain architecture vs time (exploratory). [Figure 7 to come.]
- 5.3 Cytoskeleton case study (exploratory). [Figure 8 to come.]

## 6 Discussion

[PLACEHOLDER]

## 7 Limitations

[PLACEHOLDER. Items already known:
- the pLDDT ≥ 70 filter removes 38% of genes, so results concern ordered proteins;
- the youngest age stratum is small (224 families);
- TimeTree resolves deep eukaryotic splits crudely;
- there are no AFDB models for proteins over 2,700 residues;
- polyploid species are in the panel;
- Ancestra is unverified by humans, and its duplication timing is unreliable;
- E5 acceptance is practically unreachable;
- all structures are predictions.]

## Data and code availability

[PLACEHOLDER: repository layout, DATA.md, licences; TimeTree data not redistributed.]

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
