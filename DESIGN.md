# DESIGN.md — Is protein structure a molecular clock? AFDB clusters revisited with orthology

Version: v1 (draft for review round 1) · 2026-10-30 · PI: agentic (Claude)
Inputs: PROMPT.md, LITERATURE.md, agents/data_scout_report.md, agents/ancestra_report.md

## 0. Questions (fixed by the brief) and how we turn them into testable hypotheses

**Goal 1 — revisit Illergård, Ardell & Elofsson 2009 (structure 3–10× more conserved; linear structural response) at AFDB scale.**
Illergård et al. showed that core structural change grows *linearly with ML amino-acid distance* (not with time). We therefore test two separable claims:

- **H1a (response linearity):** structural divergence is linear in amino-acid (AA) ML distance over the full range of eukaryotic ortholog divergence, for each structural metric. Alternative: saturating (exponential-to-plateau) or power-law.
- **H1b (clock linearity in time):** structural divergence of orthologs is linear in species divergence time (TimeTree). Alternative: saturating in time.
- **H1c (alphabet-matched conservation ratio):** with Foldseek's 20-state 3Di alphabet vs 20 amino acids under analogous ML corrections, structure (3Di) evolves k× slower than sequence; estimate k with CI, and its variation across fold classes / families (Illergård: 3–10×, ~5-fold between-family variation, β faster than α).
- **H1d (clock-likeness / neutrality):** the variance of structural divergence among lineages at fixed time is (i) overdispersed relative to a Poisson clock, and (ii) more or less overdispersed than sequence. Relative-rate tests on species triplets; per-family index of dispersion. A *neutral* structural clock predicts: rate ∝ substitution rate (constant structure/AA ratio across lineages), no saturation before AA saturation, dispersion comparable to AA. Departures = evidence for structural constraint/selection.
- **H1e (AlphaFold artefact):** AF2-predicted ortholog pairs are more structurally similar than experimental (PDB) pairs of the same orthologs at matched sequence divergence. Quantify the bias to qualify all conclusions.
- **H1f (phylogenetic utility, secondary):** trees built from structural distances (3Di-ML, corrected TM/LDDT; FoldTree-style) recover the species tree as well as / better than AA at deep nodes and worse at shallow nodes.

**Goal 2 — AFDB clusters in the light of evolution (open-ended; pick tractable, verifiable sub-questions).**

- **H2a (cluster audit):** AFDB structural clusters vs Ancestra hierarchical orthologous groups (HOGs): how often a HOG splits across clusters and a cluster merges HOGs; whether splits are explained by domain-architecture change (TED), length change, low pLDDT/disorder. Cluster LCA age vs HOG root age.
- **H2b (structural ortholog conjecture):** at matched AA divergence, paralog pairs are more structurally divergent than ortholog pairs (residual structural divergence > 0). Control for the fact that paralog time is unknown by matching on AA distance, not time.
- **H2c (structural least-diverged-ortholog, LDO):** after a duplication, the paralog with the longer AA branch (more diverged from an outgroup ortholog; "MDO") also carries excess structural divergence and more domain-architecture change than predicted by its AA divergence alone. Extends PMID 40813247 using Ancestra HOGs and 3Di/LDDT/TM.
- **H2d (domain architecture):** rate of domain-architecture change (TED/CATH domain strings) vs time for orthologs, and its contribution to global structural divergence (global vs per-domain TM).
- **H2e (case studies, exploratory):** (i) protein kinase fold (ePK, PF00069/PF07714; ELKs in prokaryotes) and (ii) cytoskeletal folds (actin/MreB/ParM/crenactin; tubulin/FtsZ/CetZ) — lineage-specific structural acceleration relative to AA (branches where structure changes faster than predicted by sequence), origins of eukaryote-specific structural features (Asgard archaea as outgroup).

## 1. Data

| Item | Source | Use |
|---|---|---|
| AFDB clusters v3 (paper release) | steineggerlab S3 `afdb-cluster` | cluster membership, LCA, dark flag; file 8 member→rep alignments (fident, LDDT, TM) for a *structure-selected* comparison set |
| AFDB v4 structures | FoldComp `afdb_uniprot_v4` by HTTP range | coordinates, sequence, pLDDT for panel proteins (v4 matches clusters v3) |
| Species panel | UniProt reference proteomes (2026_03), ~200 eukaryotes + ~100 prokaryotes | orthology (Ancestra), all pairwise analyses |
| Divergence times | TimeTree 6 API (pairwise median age + CI) | time axis (eukaryotes only; prokaryote times too uncertain → used only in relative analyses) |
| Domains | TED (AFDB v4 choppings + CATH labels), Pfam 38 | domain architecture |
| Experimental structures | PDB via SIFTS UniProt↔PDB mapping | AF artefact control (H1e) |
| Taxonomy | NCBI taxdump 2026-10-30 | guide tree for Ancestra, lineage labels |

Key decision (D1): **use the sequence embedded in the AFDB v4 model as the protein sequence everywhere** (Ancestra input, AA distances), so sequence and structure are identical entities; proteins without a v4 model are dropped (biases against >2,700 aa proteins; documented).

## 2. Orthology

- Run Ancestra (commit ffb0bbe, `accurate` preset) **once** on the whole panel (~300 proteomes) with the NCBI-taxonomy guide tree. Expected ≲ 4–8 h wall, ≲ 100 GB RAM (extrapolated from 20 proteomes = 17 min and 1,000 proteomes = 64 CPU-h).
- Pair labels from `hog_pair_labels.py` (deepest-common-group rule): ortholog / in-paralog / out-paralog; flag 1:1 orthologs within the HOG at the species-pair's MRCA.
- Known weakness: duplication timing (~36% agreement with Ensembl). Decision (D2): **use duplication timing only in coarse bins; base paralog analyses on AA-distance matching, never on inferred duplication age.** Ortholog calls (precision ~0.99 vs Ensembl) are used as the backbone.
- Validation: Ensembl Compara comparison on human/mouse (done by scout) + QfO-style species-tree discordance check on our panel (fraction of 1:1 ortholog gene trees concordant with species tree is not computed; instead we check that ortholog AA distance increases monotonically with TimeTree age — a sanity check, not a test).

## 3. Pair sampling (avoid structural selection bias)

- **Primary set (P-orth):** eukaryotic 1:1 orthologs, both proteins with mean pLDDT ≥ 70 and ≥ 50 aligned residues at pLDDT ≥ 70. Pairs are selected by orthology only — **never by structural similarity** (LITERATURE §0.4: in-cluster pairs are truncated by the clustering filter and can create artificial saturation).
- Non-independence: many species pairs share an MRCA node. **Unit of analysis = (HOG, MRCA node)**: for each HOG and each internal node of the species tree, sample up to k = 3 species pairs spanning that node (seeded), and summarise by the median. Fits use node-level points with HOG random effects; uncertainty by **HOG-level bootstrap** (1,000 replicates, seeds fixed).
- **Comparison set (P-clu):** member→representative alignments from AFDB clusters file 8 (structure-selected). Used to demonstrate the selection effect, not to test H1.
- **Paralog set (P-para):** in-paralogs and out-paralogs from the same HOGs, same pLDDT filters.
- Prokaryote and cross-domain pairs: for H1a/H1c (response vs AA distance) and Goal 2, not for time.

## 4. Divergence measures (per pair)

All computed on one alignment per pair to avoid alignment-method artefacts, with sensitivity analysis on a second alignment:

- Primary alignment: Foldseek (3Di+AA, `--alignment-type 2`, exhaustive pairwise on the given pairs). Secondary: TM-align mode (`--alignment-type 1`) / US-align on a 5% subsample.
- Restrict to residues with pLDDT ≥ 70 in both proteins (core-like, mirroring Illergård's core focus); also report all-residue values.
- **AA:** p-distance; ML distance under LG+Γ4 (α fixed per-run from a global fit) — implemented as batched PyTorch on GPU; validated against IQ-TREE pairwise distances and AliSim simulations with known branch lengths.
- **3Di:** p-distance; ML distance under a 3Di substitution model (Q.3Di.AF if obtainable from Puente-Lelievre/Garg et al., else a GTR-20 matrix estimated from close ortholog pairs) — same code path as AA.
- **3D:** TM-score (normalised by shorter protein and by aligned length), LDDT (Foldseek), core RMSD over aligned high-pLDDT residues; transformed distances: 1−TM, −ln(TM), 1−LDDT, −ln(LDDT) (FoldTree-style saturation-correctable forms).
- Per pair also: lengths, length ratio, mean pLDDT, MSA depth not available for v4 (D3: use pLDDT as the confidence covariate; MSA depth only if AFDB v6 MSAs can be fetched cheaply — optional).

## 5. Analyses and statistical tests

**H1a/H1b (linearity):** for each metric y and each axis x ∈ {AA-ML distance, TimeTree age}: fit (i) linear through origin y = bx, (ii) linear with intercept, (iii) saturating y = a(1−e^{−bx}), (iv) power y = ax^c. Fit on (HOG,node) medians with mixed effects (HOG random slope) or, if unstable, on node-level medians with HOG bootstrap. Compare by AIC/BIC and by out-of-range prediction (fit on x < median, predict x > median). Report the x at which saturating and linear models diverge by > 10% (d*). Separate fits by supergroup to check generality.
**H1c (ratio):** per HOG, slope of 3Di-ML vs AA-ML (total least squares through origin); distribution across HOGs; stratify by CATH class (α, β, α/β), length, mean pLDDT, disorder fraction, cluster LCA age; compare to Illergård's 3–10×.
**H1d (clock):** (i) relative-rate (Tajima-style / distance-based) tests on triplets (two ingroup species + outgroup) per HOG for AA and for each structural metric; fraction rejecting at FDR 5%, compared to the fraction expected under a permutation null (swap ingroup labels). (ii) Index of dispersion R = var/mean of root-to-tip-like distances per HOG for AA vs 3Di vs TM (FoldTree's root-to-tip variance analogue). (iii) Lineage effects: per-species mean residual (from HOG×node model) for AA and for structure; correlation across species (if structure is slaved to sequence, lineage effects correlate ~1).
**Neutral-clock interpretation:** stated as predictions above (H1d); we will not claim neutrality from fit quality alone.
**H1e (AF artefact):** ortholog pairs where both proteins have experimental structures (SIFTS, X-ray ≤ 3 Å, coverage ≥ 70%): compare TM/LDDT/3Di distances for PDB–PDB vs AF–AF of the same pairs at matched AA distance; paired Wilcoxon + slope comparison. Expected N: hundreds of pairs (human/mouse/yeast/Arabidopsis/E. coli-heavy).
**H1f (trees):** for ~500 HOGs with ≥ 20 species of 1:1 orthologs: NJ/FastME trees from AA-ML, 3Di-ML, −ln(TM), −ln(LDDT), and FoldTree combined; normalised RF and quartet distance to the NCBI/TimeTree species tree, stratified by node depth.
**H2a:** contingency of HOG ↔ cluster mapping; logistic model of "HOG split across clusters" ~ domain-architecture change + length ratio + pLDDT + AA distance.
**H2b:** residual structural divergence from the ortholog response curve (H1a fit), paralogs vs orthologs, matched on AA-ML distance (coarsened exact matching in bins) and pLDDT; HOG-stratified permutation test.
**H2c:** duplication triplets (paralog A, paralog B, nearest outgroup ortholog O). MDO = the one with larger AA-ML distance to O. Test whether structural distance to O is larger for MDO beyond AA expectation; sign test across triplets and HOG-bootstrap; domain-architecture change rate MDO vs LDO.
**H2d:** domain strings from TED (CATH superfamily labels in N→C order); architecture-change indicator per ortholog pair vs time (logistic), vs AA distance; global vs per-domain TM.
**H2e:** case-study families extracted by Pfam/TED from the panel incl. prokaryotes; FoldTree-style + AA trees; branch-wise structure/AA ratio; flag accelerated branches (outliers of ratio, bootstrap support). Exploratory; reported as such.

## 6. Controls and verification (written before results)

- **Toy data with known answers:** AliSim (IQ-TREE 3) simulated AA alignments along known branch lengths → our ML distance code must recover d (bias < 2% for d < 2). Synthetic linear vs saturating data → model selection must recover the generating model.
- **Positive control:** same protein vs itself (distance 0); AF models of identical sequences from different strains (≈0).
- **Negative control:** random non-homologous pairs from the panel (TM ~0.3, LDDT low) → defines the saturation floor for each metric; permuted time labels → no linear trend.
- **Structural-selection control:** P-clu vs P-orth curves (expect artificial saturation in P-clu).
- **Sanity:** AA distance increases monotonically with TimeTree age across nodes; human–mouse ortholog median AA identity ≈ 85%.
- Every key result reviewed by an agent that did not produce it before it enters the report.

## 7. Compute plan and milestones

| Step | Est. wall | Resumable |
|---|---|---|
| Panel + proteomes + FoldComp fetch (~3–4 M structures) | 2–3 h | yes (per species) |
| Ancestra on ~300 proteomes | 4–8 h | yes (stage checkpoints) |
| Pair sampling + Foldseek alignments (~5–10 M pairs; GPU/CPU) | 2–4 h | yes (chunks) |
| ML distances (GPU) | < 1 h | yes (chunks) |
| Statistics, figures | rolling | — |

Node 1 (until 20:52): design review, panel, Ancestra launched and checkpointed; AFDB-clusters-only analyses (P-clu, cluster LCA stats) and verification toy tests run in parallel. Node 2 (48 h): pair metrics, H1–H2 analyses, review, report.

Seeds: global seed 20261030; bootstrap/replicate seeds [1, 2, 3, 4, 5]; recorded in `config/seeds.yaml`.

## 8. Risks

- Ancestra runtime/memory blow-up on 300 proteomes → fall back to QfO-78 + selected close-pair species (~120).
- TimeTree gaps → drop species without ages from H1b only.
- AF2 artefact large → conclusions restricted to "predicted structure evolution"; H1e quantifies.
- Duplication timing unreliable → design avoids it (D2).
- Agent failure modes: silent data leakage between P-clu and P-orth, metrics computed on mismatched sequences (D1 prevents), unchecked parsers → toy tests and independent review.
