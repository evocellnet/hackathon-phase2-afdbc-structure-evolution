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

---

## Node 2 results (2026-10-30, 21:54–23:10 CET, ≤ 12 threads; exploratory/descriptive, DESIGN §11.8)

Scripts `scripts/casestudy/13–22`; tables `results/casestudy/`; figures `figures/casestudy/`. Nothing committed.
Ledger: 15 new lines in RESULTS_LEDGER.md (X6, reviewed-by = pending).

### TL;DR
- **No robust structural acceleration on either eukaryotic stem.** The 3Di/AA ratio ρ (relative to the tree) on the
  stems depends on the alignment and the method:
  - Sequence MSAs (MAFFT, FAMSA) give moderately high values: actin path to the Euk+Asgard ancestor ρ = 1.40 / 1.57,
    tubulin path to the tubulin-clade ancestor ρ = 1.14 / 1.31, at the 87th–93rd percentile, with bootstrap CIs
    covering 1 or close to it.
  - Structure-guided FoldMason, the pairwise branch engine (tdiaa) and a site-heterogeneous AA model (LG+C20) give
    ≤ 1: 0.56–1.15 across methods.
  - On constrained topologies (Asgard sister to eukaryotes; eukaryotic tubulins excluding artubulins/BtubAB, both not
    rejected by AU) the eukaryote-only stem has b_AA ≈ 0, and the path ρ is unchanged (actin 0.96–1.12, tubulin
    0.55–1.33).
- **Where the actin stem signal sits:** the high-ratio branch is the one *below crenactin* (crenactin+Euk+Asgard
  vs other prokaryotic actins; ρ = 2.40 [1.52, 5.71], robust in MAFFT), not the Euk+Asgard stem (ρ = 0.33). Its
  structure-shift columns concentrate in the hydrophobic plug, the D-loop and the C-terminal region. The two
  adjacent short branches trade length (FAMSA puts b_AA ≈ 0.02 on it), so only the path value is interpretable.
- **Within eukaryotes:** within-eukaryote and prokaryote/mixed internal branches have the same median ρ (actin
  0.77 vs 0.76, tubulin 0.75 vs 0.72). Robust high flags are almost all *terminal* branches of single proteins:
  - Tubulin: *Trichomonas* β and γ (ρ ≈ 9.4), *Paramecium* β, *Hibiscus* α.
  - Actin: *Capsella* Arp6, *Dictyostelium purpureum* Arp10.
  - Two clades: the Arp10 clade and a δ/divergent tubulin clade.

  None of these is supported by the pairwise engine at its own 97.5 % level in more than half the cases.
  Low-pLDDT models are over-represented among them (actin terminal ρ vs pLDDT Spearman −0.15, p = 0.004).
  *Giardia*, apicomplexan, fungal and plant actins/tubulins do **not** stand out at supergroup level.
- **Method agreement is the main caveat:**
  - Per-branch b_AA agrees between IQ-TREE and the pairwise engine (Spearman 0.74–0.76), but per-branch ρ does not
    (0.13 actin, 0.29 tubulin).
  - Between alignments ρ agrees moderately (0.61–0.84).
  - Only large, multiply-supported deviations should be read.

### 1. Methods added on node 2
- **Alignments (15):**
  - FAMSA 2 (default; 379/438 columns kept at ≥ 50 % occupancy).
  - FoldMason guided by the AA ML tree (`fmg`, actin; 358 columns, 0 isolated rows). The node-1 default FoldMason
    actin MSA isolated 67–69 FtsA/MreB rows.
  - For tubulin the guided run isolated rows instead, so the default FoldMason MSA is used.
  - 3Di is projected onto the AA columns, as in 06.
- **IQ-TREE on the fixed topology (16):**
  - `-te` on the binary AA ML tree; LG+G4 for AA and Q.3Di.AF+G4 for 3Di on the same columns.
  - Column bootstrap: B = 100 (ML tree) or B = 50 (constrained trees). The same resampled columns are used for AA
    and 3Di, with α fixed at the full-data value.
  - In the tubulin AA replicates IQ-TREE re-attaches sequences that became identical as cherries (RF 2–4 in about
    1/3 of replicates). Those near-zero splits are NaN in that replicate.
  - Saturation and LBA check: LG+C20+G4 AA lengths on the same topology. These are on average 1.47× (actin) and
    1.43× (tubulin) the LG lengths, but 1.98× and 1.88× on the stems.
- **Branch engine (14):**
  - Store built from the 800 AFDB models.
  - All 79,800 within-superfamily pairs, Foldseek `structurealign` 3Di+AA (the project's tdiaa).
  - Counts from `strucevo.pairs.posthoc`; aa_ml/tdi_ml at the tuned α (1.64/1.60), SE × κ (tdiaa), lddt_D with the
    tdiaa S0/SE model.
  - `fit_branches` on the UFBoot ≥ 95 collapsed topology (polytomies native), se_method = model (D2),
    ν = η/2 = 0.026 (3Di).
  - A second variant uses pairwise ML distances on the MAFFT columns (`msa`).
  - Fits: rss_rel ≈ 0.01, misfit ratio 1.13 (actin) and 1.68 (tubulin); 2.5k hidden-paralogy pair flags per
    superfamily.
  - Wall time 2–3 min per superfamily.
- **Ratio and flags (17):**
  - ρ_e = (b_3Di/b_AA)/(Σb_3Di/Σb_AA).
  - Comparable branches: supported, same kind (internal or terminal), 0.02 ≤ b_AA ≤ 2. b_AA > 2 is labelled
    "saturated" and excluded.
  - Flag: log ρ outside the 2.5–97.5 % range of the 50 comparable branches nearest in log b_AA.
  - "Robust": flagged in MAFFT, with bootstrap p_beyond ≥ 0.8, and flagged in the same direction in ≥ 2 of 3
    alignments.
  - Stem paths: the sum over consecutive stem branches, because adjacent short deep branches trade length between
    methods.
- **Constraint/LBA (13, 22):** constrained ML searches (IQ-TREE `-g`, LG+G4):
  - actin: ((Euk actin+Arps),(Asgard actins)), with prokaryotes free;
  - tubulin: eukaryotic tubulins form a clade that excludes the artubulins and BtubA/B.

  AU test (RELL 10,000):
  - actin ΔlnL 112.8, p_AU 0.061;
  - tubulin ΔlnL 31.7, p_AU 0.32.

  Neither constrained topology is rejected.
- **Structure mapping (18):**
  - Clade vs sister (the sister side is widened to ≥ 5 sequences).
  - Representatives: the highest-pLDDT member of each side's modal subfamily, superposed with US-align, giving
    per-residue CA deviations.
  - Per MSA column: total-variation distance of 3Di (TV_3Di) and AA distributions between the two sides.
    "Shift" = TV_3Di ≥ 0.5.
  - Reference numbering via `mafft --add --keeplength`: rabbit skeletal α-actin P68135 (mature numbering) and pig
    α-tubulin P02550.
  - Approximate literature loop boundaries plus UniProt features.
  - Control: 8 random unflagged eukaryotic tubulin tips (seed 20261030).

### 2. Focal question (i): eukaryotic stems vs prokaryotic homologues

| | MAFFT (IQ) | FAMSA (IQ) | FoldMason (IQ) | engine tdiaa | engine MSA | MAFFT, LG+C20 AA |
|---|---|---|---|---|---|---|
| actin, Euk+Asgard stem alone (b_AA 0.39) | 0.33 [0.00, 0.93] | 0.24 | 0.18 | 0.45 | 0.47 | 0.25 |
| actin, crenactin+Euk+Asgard stem (b_AA 0.41) | **2.40 [1.52, 5.71], robust high** | 30 (b_AA 0.02) | 1.01 | 0.85 | 1.78 | 2.21 |
| actin path, both branches (b_AA 0.80) | 1.40 [0.88, 2.38], 87th pct | 1.57 [1.06, 3.23] | 0.64 [0.33, 1.93] | 0.71 | 1.25 | 1.15 |
| actin path on constrained tree (Euk stem b_AA 0.0008) | 0.96 [0.79, 1.35] | 1.12 [0.82, 1.47] | 0.98 [0.78, 1.25] | – | – | – |
| tubulin, tubulin-clade stem (b_AA 1.49) | 1.21 [0.93, 1.75] | 1.20 | 0.60 | 0.52 | 0.93 | 0.92 |
| tubulin path CetZ+tub → tubulin clade (b_AA 1.9) | 1.14 [0.87, 1.66], 81st pct | 1.31 [0.92, 1.81] | 0.56 [0.42, 0.85] | 0.52 | 0.85 | 0.91 |
| tubulin path on constrained tree (Euk stem b_AA < 0.001) | 1.20 [0.93, 1.68] | 1.33 [0.85, 1.80] | 0.55 [0.41, 0.84] | – | – | – |

Reading:
- No method gives a stem ρ outside the comparable-length window except the MAFFT branch below crenactin. Any
  structural change concentrated at the prokaryote→eukaryote transition is therefore **not detectable robustly**
  with this sampling.
- The sequence MSAs lean high and the structure-guided and pairwise-structural methods lean low. This is the
  signature of alignment error at deep divergence:
  - structure-guided columns put more AA substitutions on deep branches (FoldMason path b_AA 4.2 vs 1.9 for MAFFT);
  - sequence MSAs misalign 3Di-variable regions.
- The engine's lddt_D/AA ratio is high on these stems (2.7). However, lddt_D for deep pairs sits near the S0 floor,
  so this is not taken as evidence.
- **AA saturation:** the tubulin path is at b_AA 1.9–4.2 and the artubulin stem (b_AA 2.13) is excluded as
  saturated. LG+C20 lengthens the stems more than the rest of the tree (1.9–2.0× vs 1.45×). Under C20 the
  sequence-MSA ρ falls to about 1 (actin path 1.15, tubulin 0.91).
- **LBA:** forcing the textbook placements does not reveal a eukaryote-specific stem. The constrained Euk-only stem
  has zero length in both superfamilies, i.e. the data place all stem change below the Asgard/artubulin divergence.

### 3. Focal question (ii): lineage-specific changes within eukaryotes

`results/casestudy/lineage_rho.tsv` (per supergroup and genus) and `flagged_branches.tsv` (all robust flags):
- Supergroup medians of ρ over comparable within-eukaryote branches range 0.82–1.32 (actin) and 0.82–1.75
  (tubulin). Land plants and Rhizaria have the highest tubulin medians (n = 11 and 5). No supergroup has more than
  2 robust high flags.
- Genera:
  - *Giardia*: actin 1.12, tubulin 0.89 (n = 3 each).
  - *Plasmodium*: actin 1.02, tubulin 2.09 (n = 6, 0 flags).
  - *Cryptosporidium*: 0.87.
  - *Trichomonas*: tubulin 1.68, with 2 robust high flags.
  - *Leishmania*: tubulin 0.46.

  The classic divergent protist actins (*Giardia*, *Plasmodium*) and apicomplexan tubulins are long in AA, but
  **not** structurally accelerated relative to AA.
- Robust within-eukaryote high flags:
  - Actin (7): mostly terminal branches. These are *Capsella* Arp6 (b_AA 0.64, ρ 2.0), *D. purpureum* Arp10,
    *Thalassiosira* Arp10, *Monodelphis* Arp5, a cavefish Arp2 and *Leishmania* Arp1 (b_AA 0.02). The one
    internal branch is the Arp10 clade (18 tips, ρ 3.9).
  - Tubulin (10): terminal *Trichomonas vaginalis* β (A2EEB7) and γ (A2G7V3), *Paramecium* β (b_AA 0.03),
    *Hibiscus* α, *Marchantia*, two macaque tubulins, *Gonapodya* and *Planoprotostelium*, plus one δ/divergent
    clade (17 tips).
  - Outside eukaryotes, the Brugia "euk_FtsZ" tip (ρ 4.9) is presumably the Wolbachia FtsZ (endosymbiont-derived
    sequence; not checked).
- **Caveats for terminal flags:** terminal 3Di lengths include AF prediction noise. Many flagged tips have mean
  pLDDT 70–80, below the typical 85–92, and actin terminal ρ correlates negatively with pLDDT. The pairwise
  engine, which subtracts the tip floor ν, puts only 3/6 (actin) and 4/10 (tubulin) robust terminal flags beyond its
  own 97.5 %. For the two *Trichomonas* tubulins the engine gives a *negative* terminal 3Di length. Treat these as
  candidates, not findings.

### 4. Structural mapping of selected branches

Tables: `structmap_summary.tsv`, `structmap_{actin,tubulin}_known_regions.tsv` (per functional region) and
`structmap_*_regions.tsv` (runs of shift columns, with UniProt features).
- **Actin, branch below crenactin / Euk+Asgard stem:**
  - Representatives: *T. brucei* actin vs *Sphaerobacter* MreB (TM 0.74, RMSD 3.1 Å, 13 % identity); 47 % of
    columns are shift columns.
  - The shift share is higher than the core average in the **hydrophobic plug** (res. 262–274: 6/7 columns, TV_3Di
    0.77 vs core 0.47, CA deviation 3.0 Å), the **D-loop** (38–52: 0.60; CA 6.1 Å vs core 2.1 Å) and the
    **C-terminal region** (349–375: 0.67; CA 6.0 Å).
  - It is lower in the nucleotide P-loops (0.12–0.33) and the 137–144 hinge (0.12).
  - Consistent with the filament-specific elements (H-plug, D-loop) distinguishing actin from MreB-type filaments,
    while the nucleotide pocket is conserved.
- **Actin Arp10 clade (ρ 3.9):** *Anopheles* Arp10 vs *T. brucei* actin. Shifts are in the sensor loop (70–78,
  9/9 columns, CA 5.6 Å), the SD2 base (60–69, 0.78) and the H-plug (0.70).
- ***Tritrichomonas* actin (robust low, ρ 0.38):** few shifts (20 %), TM 0.93 vs *Entamoeba*. The shifts sit in
  the D-loop (TV 0.63 vs core 0.25), i.e. AA change outpaces 3Di change elsewhere.
- **Tubulin clade stem:**
  - *T. cruzi* α vs *Pyrococcus* CetZ (TM 0.64, RMSD 3.3 Å); 54 % shift columns.
  - Highest in the **M-loop** (272–288: 0.88, CA 6.6 Å), the **S9–S10 loop** (1.00, CA 6.9 Å) and the **T7 loop**
    (CA 4.3 Å). These are the lateral/longitudinal protofilament contacts.
- ***Trichomonas* β (519) and γ (600), δ/divergent clade (766):** the **M-loop** is the most shifted region in all
  three (frac_shift 0.92 / 1.00 / 0.88; CA 2.1 / 4.2 / 3.1 Å). The T7 loop (β, δ), H8 (γ, δ) and the H1–S2 loop
  (δ) are also shifted.
- **Control:** in 8 unflagged tubulin tips the M-loop is also above the core share (0.37 vs 0.21). The M-loop is
  intrinsically variable, but much more so on the flagged branches.
- **Limits:**
  - One representative per side, an AF model each.
  - Loop boundaries are approximate literature ranges.
  - No test was done; these are descriptions.

### 5. Figures
- `figures/casestudy/{actin,tubulin}_tree_rho.{pdf,png}`: ML tree, x = AA length, branches coloured by log2 ρ.
  Grey = not comparable; dashed = saturated; stars = robust flags.
- `figures/casestudy/rho_distributions.{pdf,png}`: ρ vs b_AA with the comparable-length window and the named clades
  with bootstrap CIs. The display is clipped to [0.01, 100]: b_3Di at IQ-TREE's 1e-6 floor gives ρ ≈ 1e-5.
- `figures/casestudy/method_comparison.{pdf,png}`: IQ-TREE vs engine for b_AA, b_3Di and ρ, and named clades across
  methods.

### 6. Files
- Scripts: `13_constraint_trees.py`, `14_pairs_engine.py`, `15_alignments_node2.py`, `16_iqtree_fixed.py`,
  `17_ratio_analysis.py <sf> <ml|cons>`, `18_structural_map.py <sf> <branch_id…>`, `19_figures_node2.py`,
  `20_lineage_summary.py`, `21_run_after_cons.sh`, `22_au_test.sh`.
- Results:
  - `branch_ratios_{sf}_{ml,cons}.tsv`: per bipartition, all methods, windows, flags, bootstrap CIs, clade content.
  - `branch_ratio_summary_{sf}_{ml,cons}.json`, `stem_paths_{sf}_{ml,cons}.tsv`.
  - `engine_branches_{sf}_{tdiaa,msa}.tsv`, `engine_summary_{sf}.json`.
  - `lineage_rho.tsv`, `flagged_branches.tsv`.
  - `structmap_*` (focal, plus a `_control_` set).
  - `alignment_stats.tsv` (now with famsa/fmg).
- Big files (`data/interim/casestudy/`): `famsa/`, `foldmason/guide/`, `aln/*_{famsa,fmg}_*`, `iqfix/` (full and
  boot treefiles), `constraint/` (constrained trees, AU), `pairs/` (store, counts, distances), `structmap/`,
  `refs/P68135*`, `refs/P02550*`.

### 7. Caveats and what remains
- Per-branch ρ depends on the method. IQ-TREE joint-likelihood ρ and pairwise-engine ρ are nearly uncorrelated
  across internal branches, which reflects the large 3Di SEs on short internal edges. Only path-level or strongly
  flagged values are interpretable.
- One Asgard genome and no Asgard tubulin; ParM is absent. The tubulin "eukaryote stem" is in practice the tubulin
  clade stem.
- Not done:
  - PMSF/C60 tree *searches*; only C20 branch lengths on fixed trees were run.
  - Re-sampling of the 400-sequence subsets.
  - Per-region tests.
  - An independent re-derivation (V11).
- Operational notes:
  - The early-launched `21_run_after_cons.sh` started on intermediate constrained treefiles. Those outputs were
    deleted and the constrained-tree runs were redone after the searches had finished. The constrained trees used
    are final: actin lnL −203,478.5, tubulin −151,701.2.
  - The MAFFT boot-tree topology check is in §1.
