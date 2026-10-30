# Case study X6 (Goal 2-iii): actin- and tubulin-fold superfamilies — preparation report

Author: evolutionary-bioinformatics agent (Opus 5.5). Date: 2026-10-30, 12:52–15:00 CET, node 1, ≤ 8 threads
(Ancestra was using 40). Scripts: `scripts/casestudy/01–11`. Seed 20261030. Nothing committed; no `git add`.
Status: **exploratory/descriptive** (DESIGN §11.8). Nothing here has been independently re-derived yet (§6 item 7).

## TL;DR

- **Membership:** 3,437 actin-superfamily and 2,930 tubulin-superfamily proteins from 3.56M md5-matched panel proteins
  (hmmsearch `--cut_ga`, 26 InterPro Pfam HMMs with Hsp70/PilM/misato competitors). After pLDDT ≥ 70 + length
  filters, 3,033 and 2,623 remain. Subfamilies were assigned by DIAMOND best hit to 1,147 curated Swiss-Prot
  references.
- **Asgard (MK-D1 only):** 6 actins, no tubulin. Also present: crenactin, Ta0583 actins, salactin/volactin, BtubA/B
  (*Prosthecobacter*), 3 artubulins, CetZ. ParM is absent.
- **Sample:** 400 per superfamily (all 136/141 prokaryotic+Asgard members kept; eukaryotes stratified over 15
  supergroups).
- **Alignments and trees:** MAFFT L-INS-i, FoldMason, and 3Di projected on the MAFFT columns; IQ-TREE3 LG+G4 with
  UFBoot 1000. 220/397 (actin) and 231/397 (tubulin) internal nodes have UFBoot ≥ 95; the collapsed reference
  topologies are written.
- **Checks:** FtsZ │ CetZ+tubulin, FtsA, MreB, and crenactin-sister-to-Euk+Asgard all behave as expected. Asgard actins
  nest *inside* the eukaryotic Arp radiation. Artubulins and BtubA/B nest inside eukaryotic tubulins. Conventional
  actin and several Arps are not coherent.
- **First look (IQ-TREE fixed topology):** the actin Euk+Asgard stem has a *low* relative 3Di/AA ratio (ρ = 0.33,
  11th percentile), so there is no structural burst. For the tubulin stem ρ = 1.21 (MAFFT) vs 0.60 (FoldMason), which
  is alignment-dependent and AA-saturated.

## 1. Membership (step 1)

**Sequence universe.** All panel proteins (180 eukaryotes + 100 prokaryotes) whose AFDB v4 model sequence has the
same md5 as the UniProt 2026_03 sequence: 3,558,525 proteins (`01_build_seqdb.py`; FASTA from
`data/interim/structures/v4_fasta/`, i.e. the model sequences, identical by construction).

**HMMs.** `data/raw/pfam/` has no `Pfam-A.hmm`, so 26 HMMs were fetched from the InterPro API
(`/api/entry/pfam/<acc>?annotation=hmm`, 2026-10-30) into `data/interim/casestudy/pfam/` (concatenated
`casestudy.hmm`, md5 2bfa493429b507c3963d250d7c34822c). hmmsearch (HMMER 3, `--cut_ga`, 8 chunks × 1 CPU).

DATA.md-ready line:

| Item | Source | Version | Location | md5 | Licence | Notes |
|---|---|---|---|---|---|---|
| Pfam HMMs for X6 (26 families) | InterPro API `entry/pfam/<acc>?annotation=hmm` | current InterPro Pfam (HMM DATE fields 2025-05 … 2026-01; local metadata says Pfam 38.2) | data/interim/casestudy/pfam/ | casestudy.hmm 2bfa4934… | CC0 | PF00022.27 Actin, PF06723.20 MreB_Mbl, PF14450.13 FtsA, PF06406.18 StbA_N, PF21523.4 ParM_N, PF17989.8 ALP_N, PF21522.4 MreB-like_C, PF22128.3 Alp7A_like_C, PF25216.1 Volactin, PF25229.2 Salactin, PF17003.13 Actin_micro, PF02491.26 SHS2_FTSA; PF00091.32 Tubulin, PF03953.23 Tubulin_C, PF12327.14 FtsZ_C, PF21011.3 CetZ_C, PF13809.13 Tubulin_2, PF14881.13 Tubulin_3, PF21493.4/PF22334.2/PF22340.2/PF22453.2 TubZ C-termini, PF28909.1 FtsZ_C_2; competitors PF00012.27 HSP70, PF11104.15 PilM_2, PF10644.16 Misat_Tub_SegII. Fetched 2026-10-30 by X6 agent |
| Swiss-Prot reference set for subfamily labels | UniProt REST, reviewed entries of families "actin", "FtsA/MreB", "tubulin", "FtsZ", "CetZ" + P11904 (ParM) + Q8GCC5/Q8GCC1 (BtubA/B, unreviewed) | UniProt 2026-10-30 | data/interim/casestudy/refs/ | refs.fasta 39e9cef6… | CC-BY 4.0 | 1,147 sequences, 32 labels |

**Membership rules** (`03_annotate.py`):
- *Actin superfamily*: ≥ 1 GA hit to one of the 11 cytoskeletal actin-fold families above **and** best such
  sequence score > best HSP70/PilM_2 score (the MreB_Mbl HMM hits ~2,200 Hsp70s above its GA; these are removed).
- *Tubulin superfamily*: GA hit to the GTPase domain PF00091 **and** tubulin score > misato score. Proteins with only
  Tubulin_3 (+ Misat_Tub_SegII; misato/DML1, 135 hits) or only Tubulin_2 (9 large bacterial ~1,100-aa proteins) are
  excluded and listed in `data/interim/casestudy/hmm_hits_not_members.tsv`.
- *Quality*: md5 match (by construction), mean pLDDT ≥ 70, length 0.5–2× the subfamily median (`pass`).
- *Subfamily*: DIAMOND blastp `--ultra-sensitive` best hit (max bitscore, e ≤ 1e-5) against the curated reference set,
  restricted to references of the same superfamily. Labels come from UniProt subfamily annotation ("ARPn subfamily",
  "TubZ subfamily", …), else from rules on the recommended name (centractin→Arp1, ACTL6A/B→Arp4, ACTL7/ACTRT/ACTL9/10→
  "actin-like_metazoan_orphan", etc.; code in `02_build_refs.py`). Overrides by taxonomy: Asgard best hits to any
  eukaryotic label → `Asgard_actin`; Nitrososphaerota tubulins → `artubulin`; bacterial best hit to Ta0583 →
  `ALP_bacterial`; eukaryotic tubulins with best-hit identity < 40 % → `tubulin_divergent` (no ζ/η tubulin is in
  Swiss-Prot); eukaryotic MreB/FtsA → `euk_MreB`/`euk_FtsA` (flagged contaminants: e.g. *Glossina*, *Darwinula*,
  *Nematostella* — bacterial-type MreB with Hsp70-like Pfam architecture). No-hit proteins are labelled by their best
  Pfam family (Salactin/Volactin in *Haloferax*/*Halobacterium*, an FtsA, an ALP).
- Second-best label and its bitscore are kept (`second_label`, `second_bits`) to show ambiguity.

Full table: `results/casestudy/membership_all.tsv` (6,367 proteins, one row each, with Pfam architecture, best
reference, % identity, pLDDT, filters). Summary: `results/casestudy/membership_summary.tsv`.

Counts per (superfamily, domain, subfamily); `sampled` = in the downsampled tree set:

| superfamily | domain | subfamily | hits | pass (pLDDT≥70, length) | species | sampled | sampled_species | supergroups |
|---|---|---|---|---|---|---|---|---|
| actin | Archaea | MreB | 2 | 2 | 2 | 2 | 2 | 1 |
| actin | Archaea | Ta0583_archaeal_actin | 3 | 3 | 3 | 3 | 3 | 1 |
| actin | Archaea | crenactin | 1 | 1 | 1 | 1 | 1 | 1 |
| actin | Archaea | salactin | 2 | 2 | 2 | 2 | 2 | 1 |
| actin | Archaea | volactin | 1 | 1 | 1 | 1 | 1 | 1 |
| actin | Asgard | Asgard_actin | 6 | 6 | 1 | 6 | 1 | 1 |
| actin | Bacteria | ALP_bacterial | 2 | 2 | 2 | 2 | 2 | 2 |
| actin | Bacteria | FtsA | 53 | 53 | 49 | 53 | 49 | 42 |
| actin | Bacteria | MreB | 66 | 66 | 48 | 66 | 48 | 45 |
| actin | Eukaryota | Arp1 | 186 | 179 | 138 | 21 | 21 | 9 |
| actin | Eukaryota | Arp10 | 129 | 107 | 121 | 16 | 16 | 7 |
| actin | Eukaryota | Arp2 | 226 | 214 | 167 | 29 | 29 | 14 |
| actin | Eukaryota | Arp3 | 258 | 226 | 170 | 32 | 32 | 14 |
| actin | Eukaryota | Arp4 | 271 | 243 | 167 | 35 | 33 | 15 |
| actin | Eukaryota | Arp5 | 116 | 106 | 89 | 15 | 15 | 8 |
| actin | Eukaryota | Arp6 | 182 | 167 | 158 | 31 | 31 | 14 |
| actin | Eukaryota | Arp7 | 1 | 1 | 1 | 1 | 1 | 1 |
| actin | Eukaryota | Arp7_plant | 43 | 41 | 28 | 7 | 7 | 3 |
| actin | Eukaryota | Arp8 | 128 | 114 | 91 | 17 | 16 | 10 |
| actin | Eukaryota | Arp9 | 24 | 17 | 23 | 3 | 3 | 1 |
| actin | Eukaryota | actin | 1515 | 1282 | 180 | 42 | 36 | 15 |
| actin | Eukaryota | actin-like_metazoan_orphan | 145 | 137 | 29 | 6 | 6 | 4 |
| actin | Eukaryota | actin-like_other | 51 | 50 | 43 | 6 | 6 | 5 |
| actin | Eukaryota | actin_POTE | 13 | 5 | 5 | 2 | 2 | 1 |
| actin | Eukaryota | actin_divergent | 4 | 1 | 4 | 1 | 1 | 1 |
| actin | Eukaryota | euk_FtsA | 2 | 2 | 2 | 0 | 0 | 0 |
| actin | Eukaryota | euk_MreB | 7 | 5 | 6 | 0 | 0 | 0 |
| tubulin | Archaea | CetZ | 15 | 15 | 8 | 15 | 8 | 2 |
| tubulin | Archaea | FtsZ | 64 | 64 | 30 | 64 | 30 | 9 |
| tubulin | Archaea | artubulin | 3 | 3 | 3 | 3 | 3 | 1 |
| tubulin | Asgard | Asgard_FtsZ | 1 | 1 | 1 | 1 | 1 | 1 |
| tubulin | Bacteria | BtubA | 1 | 1 | 1 | 1 | 1 | 1 |
| tubulin | Bacteria | BtubB | 1 | 1 | 1 | 1 | 1 | 1 |
| tubulin | Bacteria | FtsZ | 59 | 56 | 55 | 56 | 52 | 43 |
| tubulin | Eukaryota | euk_FtsZ | 166 | 137 | 45 | 24 | 21 | 8 |
| tubulin | Eukaryota | tubulin_alpha | 1055 | 961 | 178 | 62 | 55 | 15 |
| tubulin | Eukaryota | tubulin_beta | 1074 | 941 | 180 | 62 | 55 | 15 |
| tubulin | Eukaryota | tubulin_delta | 44 | 43 | 39 | 10 | 10 | 5 |
| tubulin | Eukaryota | tubulin_divergent | 142 | 111 | 60 | 16 | 16 | 12 |
| tubulin | Eukaryota | tubulin_epsilon | 79 | 73 | 69 | 32 | 32 | 11 |
| tubulin | Eukaryota | tubulin_gamma | 226 | 216 | 173 | 53 | 52 | 15 |

Notes:
- **Asgard:** the only Asgard proteome (*Prometheoarchaeum syntrophicum* MK-D1) contributes 6 actin-superfamily
  proteins (Lokiactin-type; best hits eukaryotic actin at 21–61 % identity, one 557 aa) and **no tubulin** (only one
  FtsZ). The tubulin "Asgard stem" question therefore cannot be asked with this panel (Odinarchaeota tubulin is not in
  the panel); the tubulin focal branch is the tubulin-vs-CetZ/FtsZ stem.
- **ParM:** no ParM proper (plasmid-borne) in the reference proteomes; 2 bacterial ALPs (ALP_N/ParM_N/Alp7A_like_C
  architecture) are kept as `ALP_bacterial`.
- **Crenactin:** *Pyrobaculum calidifontis* A3MWN5 (Swiss-Prot reference itself). *Thermoproteus tenax* has no member
  passing the rules. Thermoplasmatota carry Ta0583-type archaeal actins (3); *Haloferax*/*Halobacterium* carry
  salactin/volactin (3); two archaeal MreBs (*Archaeoglobus*, *Methanobrevibacter*).
- **BtubA/BtubB:** both present in *Prosthecobacter dejongeii* (A0A7W7YIU5, A0A7W7YJ10).
- **Artubulins:** 3 Nitrososphaerota proteins (*Nitrosotenuis*, *Nitrosarchaeum*, *Nitrosotalea*) best-hit γ-tubulin
  at 25–27 % identity — the artubulins of Yutin & Koonin 2012.
- **Organellar FtsZ** in 45 eukaryotes (plants, algae, stramenopiles, a few others), labelled `euk_FtsZ`.
- pLDDT ≥ 70 plus the length filter remove 11–12 % of eukaryotic members (actin 3,301 → 2,897; tubulin 2,786 → 2,482) and ≈ 0 % of prokaryotic
  ones (3 bacterial FtsZ, incl. the *Prosthecobacter* FtsZ at mean pLDDT 65).


## 2. Downsampling (step 2; `04_sample.py`)

≤ 400 per superfamily. All `pass` prokaryotic and Asgard members kept (actin 136, tubulin 141). Eukaryotes fill the
rest by stratified round-robin over (subfamily × supergroup) cells, one random member per cell per round
(preferring species not yet drawn for that subfamily; seed 20261030), with caps on lineage-restricted labels
(POTE 3, metazoan orphan actin-likes 6, other actin-likes 6, divergent actin 2, organellar FtsZ 24, divergent tubulins
16). `euk_MreB`/`euk_FtsA` are not sampled. Result: `results/casestudy/sample_members.tsv`;
FASTA `data/interim/casestudy/{actin,tubulin}_sample.fasta`.

Eukaryotic sampled members per supergroup:

| group | actin | tubulin |
|---|---|---|
| Alveolata | 23 | 21 |
| Amoebozoa | 24 | 20 |
| Charophytes | 8 | 9 |
| Chlorophyta | 19 | 24 |
| Cryptista | 10 | 5 |
| Discoba | 19 | 23 |
| Fungi | 30 | 19 |
| Haptista | 8 | 9 |
| Holozoa_other | 19 | 13 |
| Land_plants | 22 | 22 |
| Metamonada | 12 | 15 |
| Metazoa | 29 | 31 |
| Rhizaria | 9 | 5 |
| Rhodophyta | 11 | 15 |
| Stramenopiles | 21 | 28 |


## 3. Alignments and reference topology (step 3)

- **Structures:** the 800 AFDB v4 models extracted from the panel FoldComp DBs (`05_structures.py`); Foldseek
  `createdb` gives per-protein 3Di strings; AA from the PDBs equals the sampled sequence for 800/800 (asserted).
- **AA MSA:** MAFFT L-INS-i (`--localpair --maxiterate 1000`, 4 threads).
- **Structure MSA:** FoldMason `easy-msa` (AA + 3Di outputs on shared columns).
- **3Di on the AA columns:** each sequence's 3Di letters are projected residue-by-residue onto the MAFFT columns
  (`06_alignments.py`).
- **Column filter:** occupancy ≥ 50 %, identical for AA and 3Di.

| superfamily | source | n_seq | cols_full | cols_kept | mean_gap_frac_kept | median_seq_len |
|---|---|---|---|---|---|---|
| actin | linsi | 400 | 4673 | 370 | 0.11 | 404 |
| actin | foldmason | 400 | 3218 | 311 | 0.346 | 404 |
| tubulin | linsi | 400 | 3054 | 427 | 0.134 | 444 |
| tubulin | foldmason | 400 | 2056 | 435 | 0.168 | 444 |

The MAFFT alignments are very long before filtering (many lineage-specific insertions); 370/427 columns are kept. FoldMason keeps a similar number but the actin FoldMason core is gappier (35 % gaps in kept columns).


- **AA ML tree:** IQ-TREE 3, LG+G4, UFBoot 1000, seed 20261030, on the MAFFT columns
  (`data/interim/casestudy/iqtree/{sf}_linsi_aa_ml.*`).
- **Reference topology (§11.8):** internal branches with UFBoot < 95 collapsed:
  `results/casestudy/{sf}_reference_topology_ufb95.nwk`; the full binary tree with UFBoot labels is
  `results/casestudy/{sf}_aa_ml_ufboot.nwk`.

| | actin | tubulin |
|---|---|---|
| columns | 370 | 427 |
| logL (LG+G4) | −203,366.4 | −151,669.1 |
| Γ α (AA) | 2.74 | 1.67 |
| internal branches with UFBoot ≥ 95 | 220 / 397 | 231 / 397 |
| tree length AA / 3Di (fixed topology) | 245.3 / 99.3 | 149.9 / 74.5 |
| IQ-TREE wall time (4 threads) | 1 h 26 min (UFBoot ran to the 1,000-iteration cap) | 45 min |

Figures: `figures/casestudy/{actin,tubulin}_ml_tree.{pdf,png}` (midpoint root, for display only; strips for domain
and subfamily; dots mark UFBoot ≥ 95).

## 4. Checks (expected groupings)

`results/casestudy/clade_checks_{sf}_linsi.tsv`: for each label set, the smallest side of any bipartition that
contains all its members, the number of intruders and the UFBoot of that split. Backbone views:
`python scripts/casestudy/11_backbone.py actin FtsA` / `… tubulin FtsZ,euk_FtsZ,Asgard_FtsZ`.

Per-subfamily coherence (largest clade made only of that label; `results/casestudy/subfamily_coherence.tsv`)
and key groupings (`clade_checks_*.tsv`):

**Tubulin: passes the expected checks.**
- FtsZ (bacterial, archaeal, organellar, and the single Asgard FtsZ) forms one split (145 tips, UFBoot 99).
- CetZ is monophyletic (UFBoot 96) and sister to the tubulin clade (CetZ+tubulin split UFBoot 99), as published.
- α- and γ-tubulins are each monophyletic (UFBoot 100); β 57/62 in one clade (UFBoot 78); δ, ε and the
  "divergent" bin are less coherent, as expected for label-by-best-hit.
- Eukaryotic tubulins are **not** a clade on their own. The clade of eukaryotic tubulins + 3 artubulins +
  BtubA/BtubB has UFBoot 100 (240 tips). Inside it, the artubulins sit next to the γ/α/β group (UFBoot 81) and BtubA/B
  sit among α/β (BtubB sister to α at UFBoot 42). BtubA/B nesting fits horizontal transfer from eukaryotes;
  artubulins nested inside rather than at the base disagrees with Yutin & Koonin 2012 and is a probable LBA artefact
  (artubulin stem b_AA = 2.13, the longest supported internal branch).
- Organellar FtsZ sits among bacterial FtsZ, but not in one cyanobacterial clade (largest pure euk_FtsZ clade
  9/24): plastid and mitochondrial FtsZ are mixed in the sample, and the internal support is low.

**Actin: partly passes.**
- FtsA is monophyletic (53/53, UFBoot 100). Bacterial MreB is one clade apart from one sequence (67/68, UFBoot 86).
  The two archaeal MreBs fall in an archaeal/ALP group with the Ta0583-type actins.
- Eukaryotic actin + Arps + the 6 Asgard actins form one clade (270 tips, UFBoot 100). Crenactin (*Pyrobaculum*) is
  its sister (UFBoot 100). Then come the bacterial ALPs, the Ta0583 archaeal actins and the archaeal MreBs, then MreB,
  FtsA and the salactin/volactin pair. This matches the published order (actin-proximal crenactin; MreB/FtsA/ParM
  outside).
- **But the Asgard actins are nested within the eukaryotic radiation, not sister to it.** Two are sister to an Arp2
  pair, one is sister to Arp3 (UFBoot 98), two sit near conventional actin/plant Arp7, and one (1.44 subs/site) sits at
  the base of most Arps. So eukaryotic actin+Arps alone is not a split, and "the eukaryote stem relative to Asgard"
  can only be measured on the Euk+Asgard stem.
- Arp2, Arp3, Arp6, Arp10, plant Arp7 and Arp9 are coherent (≥ 0.91). Arp1 is 0.81, Arp4 0.60, Arp5 0.47 and Arp8
  0.29; conventional actin is only 0.19, its members interleaved with Arp1/centractin and metazoan orphan
  actin-likes inside a poorly supported crown. Single long-branch protist sequences move around: e.g. one Arp4 at the
  base of all Arps (1.57 subs/site).
- ParM proper is absent from the panel.

**Verdict:** the deep backbone matches expectations in both superfamilies: MreB/FtsA/ALPs outside, crenactin as
sister to Euk(+Asgard), and the FtsZ | CetZ+tubulin split. The within-eukaryote actin/Arp backbone and the placement
of the Asgard actins and artubulins are not reliable at this sampling and model (LG+G4, 370 columns). The
reference topology collapsed at UFBoot < 95 removes most of the problematic nodes, but not the Asgard-actin
placements, which are supported at 94–100.

## 5. First look: branch-wise b_3Di / b_AA on the fixed topology

Method (`09_fixed_topology.sh`, `07_branch_ratio.py`): IQ-TREE `-te <AA ML tree>` with LG+G4 on the AA columns and
`Q.3Di.AF+G4` (`data/raw/models/Q.3Di.AF`) on the 3Di letters of the **same** columns; branches matched by
bipartition. Because the 3Di and AA scales differ, I report ρ_e = (b_3Di,e / b_AA,e) / (Σb_3Di / Σb_AA). Only
internal branches with UFBoot ≥ 95 are interpreted; branches with UFBoot < 95 stay in the binary tree as nuisance
parameters (IQ-TREE `-te` resolves polytomies at random, so the collapsed tree cannot be fixed directly; the WLS
branch engine can use the collapsed tree natively, see §6). Percentiles are among supported internal branches with
b_AA ≥ 0.02. Sensitivity: the same topology with FoldMason (structure-guided) columns.
Tables: `results/casestudy/branches_{sf}_{linsi,foldmason}_binary.tsv`, summaries
`branch_summary_{sf}_{src}_binary.json`. Figure: `figures/casestudy/branch_b3di_vs_baa.{pdf,png}`.

Focal branches (MAFFT columns; ρ percentiles among 209 actin / 193 tubulin supported internal branches):

| superfamily | focal branch | UFBoot | b_AA | b_3Di | ρ | ρ pct | b_AA pct |
|---|---|---|---|---|---|---|---|
| actin | Euk+Asgard stem (vs crenactin and the rest) | 100 | 0.386 | 0.052 | **0.33** | 11 | 84 |
| tubulin | tubulin clade stem (euk + artubulin + BtubA/B vs CetZ/FtsZ) | 100 | 1.485 | 0.896 | **1.21** | 84 | 99 |
| tubulin | FtsZ │ CetZ+tubulin split | 99 | 0.399 | 0.174 | 0.88 | 65 | 91 |
| tubulin | artubulin stem | 100 | 2.131 | 0.786 | 0.74 | 50 | 99.5 |

The eukaryote-only stem does not exist as a split in either tree (Asgard actins, and artubulins + BtubA/B, are nested
inside).

Background (supported internal branches, median ρ [IQR]):

| | within-eukaryote | within-prokaryote/mixed |
|---|---|---|
| actin | 0.77 [0.50–1.09] (n = 131) | 0.76 [0.51–1.19] (n = 74) |
| tubulin | 0.75 [0.44–1.02] (n = 104) | 0.72 [0.49–1.01] (n = 88) |

**Sensitivity, FoldMason columns (tubulin):** the tubulin-clade stem becomes b_AA 3.20, b_3Di 0.68, ρ 0.60 (35th
percentile). The structure-guided alignment makes the deep AA branch twice as long; the 3Di branch changes little.
For actin, FoldMason placed 67 prokaryotic sequences (50 FtsA, 15 MreB, 2 salactin) in their own columns. Those rows
are all gaps after the 50 % occupancy filter, so IQ-TREE refused the alignment. FtsA's swapped subdomain IC is the
likely cause; this needs a per-domain or `--refine-iters`/guide-tree FoldMason run.

**Reading (first look, descriptive):**
- The hypothesis of a structural "burst" on the prokaryote→eukaryote stem (LITERATURE §7, more 3Di change per AA
  change) is **not supported for actin**. The Euk+Asgard stem has *less* 3Di change per AA substitution than 89 % of
  supported branches (ρ = 0.33): the actin fold was already in place in the Asgard–eukaryote ancestor.
- **Tubulin is ambiguous.** ρ = 1.21 (84th percentile) with MAFFT columns but 0.60 (35th) with FoldMason columns.
  The stem's AA length is at the saturation edge (99th percentile b_AA), so its ratio depends on the alignment.
- Within-eukaryote and within-prokaryote branches have the same median ρ in both superfamilies.

## 6. Plan for the branch-engine estimate (later)

The engine in `src/strucevo/stats/branch.py` (`fit_branches`) fits pairwise distances on a fixed tree by
two-pass WLS with HC1 sandwich SEs; the design matrix is built from a `PanelTree`. For X6:

1. **Tree:** `results/casestudy/{sf}_reference_topology_ufb95.nwk` (UFBoot < 95 collapsed; the WLS design handles
   polytomies natively). Leaves are proteins, not species: build a `PanelTree` with pseudo-taxids (row index of
   `sample_members.tsv`), all ages NaN (no time axis for X6), and treat each superfamily as one "HOG"
   (`hog_species` = all 400 leaves).
2. **Pairwise distances** on the shared MAFFT columns (`aln/{sf}_linsi_{aa,3di}.fasta`): `aa_ml` = ML distance under
   LG+Γ and `tdi_ml` = Q.3Di.AF+Γ (`strucevo.distances.ml_distance`, with pair counts from `pair_counts`, Γ α fixed
   to the IQ-TREE estimates on the fixed tree: AA actin 2.74, tubulin 1.67, 3Di actin 1.97, tubulin 1.49), with Fisher SEs and saturation
   flags; optional `lddt_D` from the structures already on disk (`data/interim/casestudy/pdb/`). ~79,800 pairs per
   superfamily; design 79,800 × ≤ 797 — seconds.
3. **Saturation is the main problem.** `results/casestudy/pairwise_pdist_by_class.tsv`: eukaryote–prokaryote pairs
   have median AA p-distance 0.83 (actin) / 0.84 (tubulin), and 99 % of them exceed 0.75. Their 3Di p-distances are
   lower (median 0.62–0.64; < 1 % above 0.75). The pairwise AA ML distances that cross the focal stem will be mostly
   saturated or very imprecise, so the WLS AA fit of the stem will lean on few pairs. Mitigations to pre-register:
   (a) report the stem from both engines (pairwise WLS and the joint-likelihood IQ-TREE `-te` used here);
   (b) restrict the stem estimate to pairs whose AA SE/d < 0.5; (c) add a 3Di-only and an LDDT-based version, which
   saturate later.
4. **Outputs:** b_AA, b_3Di, their covariance and ρ per branch with sandwich CIs; flag "accelerated" branches as
   those whose ρ lies above the 95th percentile of the HOG-level ρ distribution from the main analysis
   (DESIGN §5, Goal 2-iii), and the eukaryote stem against within-eukaryote and within-prokaryote branches.
   Descriptive only (§11.8).
5. **Sensitivity:** FoldMason columns (already prepared), the downsampling seed, and AU/constraint trees for the
   alternative placements of Asgard actins and artubulins/BtubAB (§4 here).

## 7. Caveats

- **One Asgard genome, no Asgard tubulin.** MK-D1 gives the actin side only. Adding *Lokiarchaeum ossiferum* or
  Heimdall-/Odinarchaeota would require non-AFDB structures (predict with ColabFold/ESMFold) — a protocol deviation.
- **Long-branch attraction at the stems.** The focal branches are the longest supported internal branches in both
  trees (≥ 98th percentile of b_AA). Their placement (artubulin, BtubA/B, Asgard actins) carries LBA risk under
  LG+G4; site-heterogeneous models (LG+C20/C60+F+G, PMSF) would be the right check but cost hours at 400 taxa.
- **Saturation of b_AA on deep branches**, and alignment dependence (MAFFT vs FoldMason) of the deep AA branch
  lengths: the ratio at the stem moves by a factor of ~2 between alignments (§5). Any claim about the stem must be
  stable across both.
- **Scale:** ρ is relative to the tree-wide ratio, which is dominated by terminal branches (71–77 % of AA tree length).
- **Contaminant-like eukaryotic MreB/FtsA** (7 + 2) were excluded from sampling, not from the membership table.
- **Subfamily labels** are best-hit labels; eukaryotic "tubulin_divergent" and "actin-like_other/orphan" are
  heterogeneous bins, and Arp7/Arp9 (fungal SWI/SNF Arps) and plant Arp7 have few references.
- **Model fit:** Q.3Di.AF was trained on AFDB-derived 3Di alignments; its use on MAFFT (sequence) columns is a
  projection, not a structure alignment. The FoldMason run is the structure-consistent counterpart.
- UFBoot ≥ 95 is used as support; the trees were not checked with SH-aLRT (not requested).

## 8. Files

- Scripts: `scripts/casestudy/01_build_seqdb.py`, `02_build_refs.py`, `03_annotate.py`, `04_sample.py`,
  `05_structures.py`, `06_alignments.py`, `07_branch_ratio.py <sf> <src>`, `08_figures.py`,
  `09_fixed_topology.sh`, `10_pairwise_saturation.py`, `11_backbone.py`, `12_purity.py`.
- Small results (`results/casestudy/`): `membership_all.tsv`, `membership_summary.tsv`, `sample_members.tsv`,
  `alignment_stats.tsv`, `{sf}_aa_ml_ufboot.nwk`, `{sf}_reference_topology_ufb95.nwk`, `clade_checks_{sf}_*.tsv`,
  `subfamily_coherence.tsv`, `backbone_{sf}.txt`, `branches_{sf}_{src}_binary.tsv`,
  `branch_summary_{sf}_{src}_binary.json`, `pairwise_pdist_by_class.tsv`.
- Figures (`figures/casestudy/`): `actin_ml_tree`, `tubulin_ml_tree`, `branch_b3di_vs_baa` (.pdf/.png).
- Big (`data/interim/casestudy/`): `panel_all_v4match.fasta` (3.56M seqs), `hmm_chunks/`, `hmm_domtbl_all.tsv.gz`,
  `pfam/`, `refs/`, `members_all.fasta`, `members_vs_refs.tsv`, `{sf}_sample.fasta`, `{sf}_linsi.fasta`, `pdb/`,
  `fsdb/`, `{sf}_3di.fasta`, `foldmason/`, `aln/`, `iqtree/`.
- **DATA.md:** not edited. The two rows in §1 are ready to paste.
