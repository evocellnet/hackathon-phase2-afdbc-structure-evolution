# Case study X7 (Goal 2-iii): protein kinase-like (PKL) fold — structural vs sequence change per branch

Author: evolutionary-bioinformatics agent (Opus 5.5). Date: 2026-10-30, 22:24 → 00:40 CET, node 2, ≤ 8 threads.
Scripts: `scripts/kinase/01–12`. Seed 20261030. Nothing committed; no `git add`.
Status: **exploratory/descriptive** (DESIGN §11.8). Not independently re-derived. Methods copied from X6
(`agents/casestudy_report.md`) with the changes listed in §2.

## TL;DR

- **Membership:** 92,042 panel proteins (180 euk + 100 prok, md5-matched) carry ≥ 1 GA hit to one of the 48 Pfam
  clan CL0016 (PKL) HMMs; 81,943 have Pkinase/PK_Tyr_Ser-Thr as their best PKL domain (ePK). The kinase domain is
  well predicted: domain pLDDT ≥ 70 removes < 1 % of ePK domains; the HMM-coverage ≥ 0.70 filter removes 17 %
  (fragments/pseudokinase pieces).
- **Tree set:** 500 kinase domains (domain only, envelope coordinates): 263 eukaryotic ePKs (≈ 33 per group:
  AGC, CAMK, CK1, CMGC, STE, TK, TKL, Other, + 10 "TK-like" non-holozoan), 85 prokaryotic/Asgard ePK-like, 152
  ELK/atypical (RIO, ABC1/UbiB, APH, choline kinase, Kdo, Bud32).
- **Tree:** IQ-TREE3 LG+G4, UFBoot 1000, 244 MAFFT L-INS-i columns; 234/497 internal nodes have UFBoot ≥ 95.
  Holozoan TKs are one clade (32 TK, UFBoot 96; 33 + one *Guillardia* TK-like, UFBoot 100). CK1 groups with the
  ELKs (probable long-branch artefact).
- **(i) TK origin:** the two branches at the base of the TK clade together have ρ = 1.68 (MAFFT), 1.47 (FAMSA),
  2.00 (FoldMason). That is the 88th–91st percentile of single supported within-eukaryote ePK branches and about
  the 95th percentile of pairs of such branches. The direction is the same with all three alignments, but the
  branches are short (b_AA ≈ 0.17) and the value sits at the edge of the background. **This is a weak signal, not
  a burst.** The other group stems are not exceptional: AGC 0.42, CAMK 0.90, CK1 1.18 (FoldMason 0.19), CMGC
  0.58–1.67, TKL 0.62, Other 0.21.
- **(ii) Prokaryotic vs eukaryotic ePKs:** internal branches of prokaryotic ePKs have pooled ρ 1.02–1.22× that of
  eukaryotic ePK internal branches. All branch-bootstrap 95 % CIs include 1 (e.g. MAFFT 1.11 [0.84, 1.41]).
  Terminal branches give 0.94–1.07×. Prokaryotic ELK internal branches have the *lowest* ρ (pooled 0.71–0.76).
  **We see no evidence that prokaryotic ePK-like or ELK kinases change more in structure per AA change.**
- **Caveats dominate:** ePK–ELK pairs are AA-saturated (97 % of pairs have p-distance > 0.75). Per-branch ρ agrees
  only moderately between alignments (Spearman 0.74 FAMSA, 0.59 FoldMason). Of the 12 "ρ-high" flagged branches,
  only 1 is flagged with all three alignments.

## 1. Membership (`01_hmmsearch.py`, `03_annotate.py`)

- **Sequence universe:** the X6 set of 3,558,525 panel proteins whose AFDB v4 model sequence equals the UniProt
  2026_03 sequence. The same 8 FASTA chunks were reused.
- **HMMs:** all 48 families of Pfam clan CL0016 (from the Pfam 38.2 clan table in `data/raw/pfam`), fetched from the
  InterPro API (`entry/pfam/<acc>?annotation=hmm`, 2026-10-30, one transient HTTP 500 retried) into
  `data/interim/kinase/pfam/` (`pkl.hmm` md5 e9f29042dadb33fb8f0b6010b01ed8ae).
- **Search:** `hmmsearch --cut_ga`, 8 chunks × 1 CPU, about 5 min.
- **Kinase domain:** the single highest-scoring CL0016 domain hit (clan competition); its family is the protein's
  family. Boundaries are the HMM envelope. Multi-kinase proteins (JAK, RSK, …) contribute only their best domain.
- **Quality filters:** HMM coverage ≥ 0.70; **domain** mean pLDDT ≥ 70 (per-residue pLDDT from
  `panel_*_v4_plddt.npz`, averaged over the envelope); domain length ≤ 2× the family median.
- **ePK group labels:** DIAMOND blastp `--ultra-sensitive` best hit (bitscore) of the domain sequence against 496
  human kinase domains (483 reviewed human proteins of the "protein kinase superfamily"; UniProt REST 2026-10-30;
  one reference per "Protein kinase" DOMAIN feature).
  - The group label is parsed from the UniProt SIMILARITY family string: AGC/CAMK/CK1/CMGC/STE/TKL Ser/Thr family;
    "Tyr protein kinase" → TK.
  - Other ePK families (generic "Ser/Thr protein kinase family", NEK, BUD32, PAN3) → Other. RIOK → RIO; ADCK → ABC1.
    Alpha-type kinases are excluded.
  - This is an operational proxy for the Manning groups, not Manning's own table.
  - Script: `02_refs.py`; files `data/interim/kinase/refs/`.
- **Two relabels**, decided after inspection:
  - **Eukaryotic best-hit-TK outside Holozoa → `TKlike_nonholozoan`.** These have median best-hit identity 30 %,
    versus 74 % for metazoan TKs: 676 domains in plants, amoebozoans, fungi, stramenopiles, etc. True PTKs are
    holozoan-specific.
  - **Eukaryotic "Kdo" best hits are all TP53RK/Bud32-like** (best human hit TP53RK, 51 domains). They are reported
    as Kdo/Bud32.
- **Label ambiguity** (fraction with a bitscore margin < 10 to the second-best group, sampled ePKs): TKL 0.59,
  STE 0.39, Other 0.19, AGC 0.09, CAMK 0.09, CK1/CMGC/TK 0. **TKL and STE labels are therefore soft.**
- **Prokaryotic ePKs** are kept as one class. Their best human hits are Other/CAMK (mostly STK36/PNCK/NEK-like).

**Counts:** `results/kinase/membership_summary.tsv`, with the full per-protein table in `membership_all.tsv.gz`.
Domains passing the filters:

| class | domains passing | | class | domains passing |
|---|---|---|---|---|
| TKL | 21,982 | | RIO (euk / arc / bac) | 459 / 52 / 4 |
| CAMK | 10,268 | | ABC1 (euk / bac / arc) | 1,239 / 59 / 20 |
| Other | 7,970 | | APH (euk / bac / arc) | 452 / 121 / 12 |
| STE | 7,950 | | Choline kinase (euk) | 528 |
| CMGC | 7,668 | | Kdo (bac) / Bud32-like (euk) | 14 / 51 |
| TK (Holozoa) | 4,409 | | other CL0016 families | recorded, not sampled |
| AGC | 4,358 | | | |
| CK1 | 1,936 | | | |
| TK-like non-holozoan | 676 | | | |
| prokaryotic ePK | 722 (bac 694, arc 23, Asgard 5) | | | |

## 2. Downsampling, structures, alignments (`04`–`06`)

**Sample.** 500 domains, seed 20261030, stratified round-robin over cells with species de-duplication
(`results/kinase/sample_members.tsv`):
- Prokaryotic ePK: ≤ 80 over (domain, phylum), giving 80 from 38 phyla, plus all 5 Asgard (MK-D1) ePKs.
- Each ELK class × domain: ≤ 12 over supergroups, giving 152 in total. No prokaryotic choline kinase passed the
  filters.
- TK-like non-holozoan: 10.
- Eukaryotic ePK groups: 32–33 each over supergroups.
- TK draws only from Metazoa (17) and Holozoa_other (16: 9 *Capsaspora*, 5 *Monosiga*, 2 *Sphaeroforma*).

**Domain pLDDT of the sample:** median 88, min 71.

**Structures.** AFDB v4 models taken from the panel FoldComp DBs and **cropped to the domain envelope**. Foldseek 3Di
was computed on the cropped domains, so the projected 3Di and FoldMason use the same structural input. AA from the
PDB equals the sampled sequence for 500/500 (asserted).

**Alignments.** Columns with occupancy ≥ 50 % are kept, identically for AA and 3Di; no row falls below 10 %
occupancy.

| source | columns full → kept | mean gap fraction (kept) |
|---|---|---|
| MAFFT L-INS-i (primary) | 2,088 → 244 | 0.155 |
| FAMSA 2 | 1,663 → 252 | 0.097 |
| FoldMason easy-msa (default; no guide tree needed, unlike X6 actin) | 1,346 → 242 | 0.228 |

## 3. Tree and fixed-topology branch lengths (`07_iqtree.sh`, `08_branch_ratio.py`)

**ML tree:** IQ-TREE3, LG+G4, UFBoot 1000, 8 threads, 1 h 36 min (ran to the 1,000-iteration cap). Fit: logL
−192,559; α 1.36; AA tree length 339 subs/site over 500 tips. 234 of 497 internal nodes have UFBoot ≥ 95.

**Fixed topology:** branch lengths re-estimated on the same tree for each column set, with LG+G4 for AA and
Q.3Di.AF+G4 for 3Di (same columns). Tree lengths AA / 3Di: MAFFT 339 / 197; FAMSA 429 / 272; FoldMason 487 / 209.
FoldMason lengthens the AA tree by 44 %; the 3Di tree changes little.

**Ratio:** ρ_e = (b_3Di/b_AA)_e / (Σb_3Di/Σb_AA). Only supported branches (terminal, or UFBoot ≥ 95) with b_AA ≥ 0.02
are interpreted.

**Rooting (for classifying clades only):** on the MRCA clade of all ELK domains. That clade also contains all 33 CK1
domains and 4 prokaryotic ePKs (2 ePK "Other"). The CK1+ELK grouping disagrees with the accepted ePK monophyly
(CK1 is a bona fide ePK group) and is most likely long-branch attraction/saturation under LG+G4 on 244 columns.
Rooting therefore puts CK1 on the "outgroup" side. Its stem is still a pure, supported clade (33/33, UFBoot 100).

**Topology checks:**
- **TK:** holozoan TK MRCA = 34 tips (33 TK + 1 *Guillardia theta* TK-like), UFBoot 100. The pure TK clade is 32
  tips, UFBoot 96. Its sister is STE/TKL (UFBoot 88).
- **Other groups** are not monophyletic as wholes, but their largest pure clades are 27/33 (AGC, UFBoot 99),
  25/33 (CMGC), 14 (STE), 11 (CAMK), 6 (Other) and 5 (TKL). These labels come from best hits, and TKL/STE/Other are
  ambiguous (§1).
- **Prokaryotic ePKs:** the largest pure clade is 64/85 (UFBoot 91), with no eukaryotic group nested inside.
- **Asgard (MK-D1):** the ePK-like domains sit with bacterial ePKs. RIO and APH domains sit inside their own families.
- **RIO, ABC1 and choline kinase** are each one clade (RIO 31/31 but UFBoot 20; ABC1 36/36, UFBoot 96; ChoK 12/12,
  UFBoot 100). APH and Kdo/Bud32 are not monophyletic.

**Tables:** `results/kinase/branches_{linsi,famsa,foldmason}.tsv`, `branch_summary_*.json`, `class_rho_*.tsv`,
`class_rho_baa_binned_*.tsv`, `group_stems_*.tsv`, `stems_by_source.tsv`, `class_rho_by_source.tsv`.

## 4. Question (i): TK origin and group stems

The background is the supported internal within-eukaryote ePK branches with b_AA ≥ 0.02. On MAFFT columns this is
n = 127, median ρ 0.73 [IQR 0.39–1.00], q95 2.07.

| stem (UFBoot) | ρ MAFFT | ρ FAMSA | ρ FoldMason | pct vs background (M / F / FM) | b_AA (MAFFT) |
|---|---|---|---|---|---|
| TK holozoan MRCA (100) | 1.51 | 0.98 | 2.22 | 87 / 78 / 90 | 0.101 |
| TK pure clade (96) | 1.93 | 8.87* | 1.41 | 94 / 100 / 76 | 0.065 |
| **TK origin segment (2 branches pooled)** | **1.68** | **1.47** | **2.00** | 91 / 90 / 88 (pairs: 95 / 95 / 96) | 0.166 |
| CK1 (100) | 1.18 | 0.91 | 0.19 | 79 / 73 / 5 | 0.523 |
| CMGC MRCA, 33/36 (85, unsupported) | 1.67 | 1.67 | 1.57 | 91 / 92 / 80 | 0.177 |
| AGC pure 27 (99) | 0.42 | 0.45 | 0.33 | 27 / 31 / 11 | 0.149 |
| CAMK pure 11 (92) | 0.90 | 0.65 | 1.36 | 66 / 54 / 75 | 0.051 |
| TKL pure 5 (100) | 0.62 | 0.22 | 0.28 | 43 / 11 / 8 | 0.234 |
| Other pure 6 (94) | 0.21 | 0.12 | 0.41 | 11 / 9 / 19 | 0.154 |

\* FAMSA gives this branch b_AA = 0.011, so the ratio is unstable. That is why the two-branch segment
(`results/kinase/tk_segment_by_source.tsv`, `12_focal_extra.py`) is the preferred summary.

**Reading.** The branches leading to holozoan TKs carry somewhat more 3Di change per AA change than the typical
within-eukaryote ePK branch. The direction holds with all three alignments, at about the 88th–91st percentile of
single branches. It is not an outlier: it never exceeds q95 of single branches, and it sits at about the 95th
percentile of random pairs of branches. The AA length of the segment is modest (0.17 subs/site), so the estimate
is noisy; no column bootstrap was run in the time box. No other group stem is consistently high. CK1 even flips
from 1.18 to 0.19 with FoldMason, and its stem is the long, LBA-prone branch noted above. **Verdict:** at most a
weak, alignment-consistent tendency for TK origin to be structurally "faster" than background. There is no
evidence of a structural burst at the origin of a eukaryote-specific group.

## 5. Question (ii): prokaryotic vs eukaryotic ePK-like / ELK

Pooled ρ (Σb_3Di/Σb_AA, relative). Branch bootstrap B = 2000; the CI is for the ratio prokaryotic/eukaryotic
(`results/kinase/prok_vs_euk_by_source.tsv`).

| branches | MAFFT | FAMSA | FoldMason |
|---|---|---|---|
| internal prok ePK / euk ePK (n 37 / 127) | 0.86 / 0.77 → 1.11 [0.84, 1.41] | 0.83 / 0.68 → 1.22 [0.98, 1.48] | 1.14 / 0.96 → 1.19 [0.97, 1.44] |
| same, b_AA 0.05–0.5 only | 1.02 [0.78, 1.35] | 1.20 [0.93, 1.50] | 1.15 [0.88, 1.44] |
| terminal prok ePK / euk ePK (85 / 269) | 1.07 [0.98, 1.17] | 1.04 [0.95, 1.13] | 0.99 [0.91, 1.08] |
| terminal, b_AA 0.05–0.5 | 0.97 [0.84, 1.12] | 1.01 [0.85, 1.18] | 0.94 [0.82, 1.09] |
| internal prok ELK (pooled ρ, n 28) | 0.71 | 0.74 | 0.76 |
| internal euk ELK (n 20–21) | 1.13 | 1.17 | 0.86 |

**Reading.** Internal branches of prokaryotic ePKs are 2–22 % "more structural" per AA change than eukaryotic ePK
internal branches with all three alignments. Every CI includes 1, and the excess shrinks when b_AA is matched.
Terminal branches show nothing. Prokaryotic ELKs (APH, Kdo, RIO, ABC1) are, if anything, *less* structurally
divergent per AA change. **We do not support** "prokaryotic ELK/ePK-like proteins diverge more in structure per AA
change". Within the sensitivity of this design (about ±30 %), the ratio looks the same across the PKL fold. Note
that branches are not independent (shared tree), so the bootstrap CIs are approximate.

## 6. Flagged branches

12 supported internal branches have ρ > q95 (2.14) on MAFFT columns (`branches_linsi.tsv`, flag `rho_high`):
- 6 within eukaryotic ePK: shallow pairs/triplets of Other, CMGC, CAMK and CK1;
- 2 within prokaryotic ePK;
- 3 within or around eukaryotic ELK;
- the ABC1 clade stem (b_AA 0.54, ρ 2.58; with FoldMason its b_AA hits the 10.0 upper bound, i.e. saturated).

Agreement across alignments is poor. FAMSA shares 5/12 flags and FoldMason 2/12. Only one branch is flagged with
all three: a CMGC pair, ρ 3.1 / 2.6 / 3.3 (`flagged_rho_high_all_sources.tsv`). Most single-branch flags are
therefore alignment noise on short branches.

Saturation/LBA-risk branches (b_AA > 1): the Kdo pure-clade stem (1.09). With FAMSA and FoldMason, the RIO1, ABC1
and choline kinase stems all exceed 1–10 subs/site.

## 7. Caveats

- **Saturation.** Median AA p-distance is 0.76 between ePK groups and 0.81 between ePKs and ELKs; 60 % and 97 % of
  those pairs exceed 0.75 (`pairwise_pdist_by_class_linsi.tsv`). 3Di saturates later (median 0.52 / 0.70; ≤ 18 %
  above 0.75). Deep AA branch lengths, and hence deep ρ, are therefore biased: ρ goes up when b_AA is underestimated.
  This argues for caution on any deep stem, including the TK and ELK stems.
- **LBA.** CK1 clusters with the ELKs. Group monophyly is partial. LG+G4 on 244 columns is too simple for
  ePK-vs-ELK depths; site-heterogeneous models (C20/C60, PMSF) were not run.
- **Alignment dependence.** Per-branch ρ Spearman between alignments is 0.74 (FAMSA) and 0.59 (FoldMason) on
  supported internal branches, and 0.80 / 0.61 on terminal branches. Conclusions above use only directions that hold
  with all three alignments.
- **Labels.** Groups come from best hits to human references. TKL/STE/Other are soft, and non-human-like protist
  kinases get forced labels. The TK clade includes a *Guillardia* TK-like domain.
- **Sampling.** The TK group is half non-metazoan holozoans by design (supergroup stratification), so the "TK stem"
  is the holozoan PTK stem. The prokaryotic set is small (85) and bacterial-dominated. The Asgard data come from one
  genome.
- **What is not done:** no column bootstrap of ρ (X6's `16_iqtree_fixed.py` pattern would take about 1 h at 8
  threads); no pairwise WLS engine; no constraint/AU tests of CK1 placement; domain cropping makes 3Di at the domain
  termini differ slightly from full-chain 3Di; regulatory elements (activation loop, αC–β4) were not analysed
  separately.

## 8. Files

- **Scripts** (`scripts/kinase/`): `01_hmmsearch.py`, `02_refs.py`, `03_annotate.py`, `04_sample.py`,
  `05_structures.py`, `06_alignments.py run|filter`, `07_iqtree.sh ml|fixed`, `08_branch_ratio.py <src>`,
  `09_figures.py`, `10_pairwise_pdist.py <src>`, `11_compare_sources.py`, `12_focal_extra.py`.
  Logs are in `logs/kinase/`.
- **Results** (`results/kinase/`): `membership_all.tsv.gz`, `membership_summary.tsv`, `sample_members.tsv`,
  `alignment_stats.tsv`, `branches_*.tsv`, `branch_summary_*.json`, `class_rho_*.tsv`, `group_stems_*.tsv`,
  `stems_by_source.tsv`, `class_rho_by_source.tsv`, `tk_segment_by_source.tsv`, `prok_vs_euk_by_source.tsv`,
  `flagged_rho_high_all_sources.tsv`, `sensitivity_summary.json`, `pairwise_pdist_by_class_*.tsv`.
- **Figures** (`figures/kinase/`):
  - `pkl_tree_rho` — tree, branches coloured by log2 ρ, strips for category and group;
  - `rho_by_group` — ρ by branch class, and by ePK group with the stems marked;
  - `rho_alignment_sensitivity`.
- **Big** (`data/interim/kinase/`): `pfam/`, `hmm/`, `hmm_domtbl_all.tsv.gz`, `refs/`, `pkl_domains.fasta`,
  `pkl_vs_refs.tsv`, `pkl_sample.fasta` (md5 13706d36…), `pdb/` (cropped), `fsdb/`, `pkl_3di.fasta`, `aln/`,
  `foldmason/`, `iqtree/`.
- **DATA.md lines** (not edited):
  - Pfam CL0016 HMMs (48; InterPro API, Pfam 38.2, CC0, `pkl.hmm` md5 e9f29042…);
  - human kinase reference set (UniProt REST 2026-10-30, reviewed human "protein kinase superfamily", CC-BY 4.0,
    `human_pks_uniprot.tsv` md5 cbffc07e…, `refs_kd.fasta` md5 8f91ad9f…).

