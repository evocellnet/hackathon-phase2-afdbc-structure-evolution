# Species panel, proteomes, AFDB v4 structures and TimeTree ages: report

Author: data-engineer agent (Opus 5.5). Date: 2026-10-30, about 10:30–12:20 CET, node 1 (project wall clock 1.1 → 2.8 h).
Scripts: `scripts/panel/` (`common.py`, `select_species.py`, `download_proteomes.sh`, `map_proteomes_to_afdb.py`, `fetch_foldcomp.py`, `extract_structures.py`, `timetree_pairwise.py`, `sanity_checks.py`, `afdb_v4_tables.py`, `ted_md5_table.py`).
Environment: `conda_envs/main`. Seed: 20261030.
Logs: `logs/panel_*.log`.

## TL;DR

| item | eukaryotes | prokaryotes |
|---|---|---|
| species | 180 | 100 (60 bacteria in 50 phyla, 40 archaea in 10 phyla) |
| canonical proteins (UniProt 2026_03) | 3,460,494 | 277,544 |
| with an AFDB v4 model | 3,290,300 (95.1%) | 275,566 (99.3%) |
| v4 model sequence identical to the 2026 sequence | 3,283,644 (94.9% of all proteins) | 274,881 (99.0%) |
| v4 model exists but sequence changed since 2021 | 6,656 (0.20% of models) | 685 (0.25%) |
| no v4 model | 170,194, of which 16,160 are longer than 2,700 aa | 1,978, of which 175 are longer than 2,700 aa |
| structures fetched / failed | 3,290,300 / 0 | 275,566 / 0 |
| per-species exact-match coverage (min / median) | 0.831 / 0.955 | 0.810 / 0.996 |
| mean pLDDT (protein mean) | 70.4 (55% of models have mean ≥ 70) | 85.0 |

- **Ancestra inputs:** `data/interim/ancestra_input/` (eukaryotes) and `data/interim/ancestra_input_prok/` (prokaryotes). Each holds one `<UPID>_<taxid>.fasta` per species, with the 2026 canonical sequences and the accession as header. This follows PI decision v2.
- **Structures:**
  - FoldComp DBs: `data/interim/structures/panel_v4` and `panel_prok_v4` (each with `.index`, `.lookup` and `.dbtype`).
  - Foldseek DBs: `data/interim/foldseek/panel_db` and `panel_prok_db`, built directly from the FoldComp DBs.
- **Per-protein table:** `data/processed/panel/protein_table.parquet` (and `_prok`). Columns: accession, species, `len_2026`, `md5_2026`, `has_v4`, `len_v4`, `md5_v4`, `seq_match`, `mean_plddt`, and the fractions of residues with pLDDT ≥ 70 and ≥ 90.
- **TimeTree:** 16,110 eukaryotic species pairs were queried. 15,755 were dated and **355 got no answer**: 109 involve *Fonticula*, and 203 involve the three metamonads (68 for *T. vaginalis*, 68 for *Tritrichomonas*, 67 for *Giardia*).

## 1. Species selection (`select_species.py`)

### Source and gates

The source is the 37,190 UniProt reference proteomes of release 2026_03. Metadata come from the proteomes REST API; note that `query=reference:true` is required, because `proteome_type:1` returns nothing. One proteome is kept per NCBI species.

Gates:

- **AFDB v4 exact coverage ≥ 0.80.** This is the fraction of canonical accessions whose AFDB v4 model has an identical sequence. It follows PI decision v2.
  - To apply this gate before selecting, I downloaded all 1,516 plausible eukaryotic candidates (7 GB, 4 min) and mapped them exactly. These were proteomes passing BUSCO with a taxId-based coverage proxy ≥ 0.5, plus every proteome of a requested sister species.
  - The taxId proxy is unreliable in both directions. *E. histolytica* has a proxy of 0.001 but 0.93 exact coverage; rat has a proxy of 0.65 but 0.555 exact.
  - **Distribution over the 1,516 candidates:**

    | exact coverage ≥ | 0.5 | 0.7 | 0.8 | 0.9 | 0.95 |
    |---|---|---|---|---|---|
    | candidates passing | 1,368 | 1,340 | 1,321 | 1,216 | 871 |

    At 0.80 every supergroup is still covered. 0.90 would have removed, among others, *Chlamydomonas* (0.83), chimpanzee (0.85), *Anopheles gambiae* (0.85) and several *Plasmodium* species.
- **BUSCO complete:**
  - ≥ 80% for Metazoa, Fungi and Viridiplantae.
  - ≥ 50% for protist supergroups, whose BUSCO eukaryota scores are systematically low.
  - ≥ 15% for Metamonada (*Giardia* 23.9%).
  - QfO-2022 proteomes are exempt.
  - One documented override: *C. briggsae* UP000008549, where UniProt reports C:0.9%, a metadata artefact (exact coverage 0.96).
  - 17 panel species have BUSCO below 80%. They are listed in `species_panel.tsv`; they are mostly protists, plus *Ciona* (78.7%) and *Monosiga* (78.8%).
- **TimeTree:** each species must be dated by the pairwise API against its closest pool member, or against a second, more distant one. 1,158 of 1,290 species passed.

### Algorithm

1. **Mandatory species:** the QfO-2022 eukaryotes that pass the gates (38 of 48), plus 33 requested sister pairs.
2. **Breadth:** the best-quality species for each NCBI phylum, using class level within Chordata, Arthropoda, Ascomycota, Basidiomycota, Streptophyta, Chlorophyta and Apicomplexa. This added 60 species.
3. **Greedy fill (34 species), aimed at log-uniform node ages:**
   - Adding species c creates a new divergence node of age min over selected s of age(c, s).
   - At each step we pick the candidate whose new node falls into the least-filled log10 age bin (1 Mya to 2 Gya, bin width 0.3).
   - Constraints: supergroup quotas, at most 3 species per genus, and a small quality term.
   - Ages used: TimeTree pairwise ages for congeneric pairs, and median ages of NCBI MRCA nodes otherwise.

The result is reproducible: rerunning the script gives a byte-identical panel.

### Composition (180 eukaryotes)

| group | species |
|---|---|
| Metazoa | 66 |
| Fungi | 40 |
| Land plants | 20 |
| Charophyte (*Klebsormidium*) | 1 |
| Chlorophyta | 7 |
| Rhodophyta | 4 |
| Alveolata | 10 |
| Stramenopiles | 8 |
| Rhizaria (*Plasmodiophora*) | 1 |
| Discoba | 8 |
| Metamonada (*Giardia*, *Trichomonas*, *Tritrichomonas*) | 3 |
| Amoebozoa | 6 |
| Other opisthokonts (choanoflagellates, *Capsaspora*, *Sphaeroforma*, *Fonticula*) | 4 |
| Haptista (*Chrysochromulina*) | 1 |
| Cryptista (*Guillardia*) | 1 |

- **Metazoa breadth:** Porifera is absent (it failed the gates). Present: Placozoa (*Trichoplax*), Cnidaria, Rotifera, Tardigrada, Brachiopoda, Mollusca, Annelida, Platyhelminthes (*Macrostomum*), Nematoda, Echinodermata, Cephalochordata, Tunicata, and the vertebrate and arthropod classes.
- **Fungi breadth:** Chytridiomycota, Neocallimastigomycota, Blastocladiomycota, Mucoromycota, Mortierellomycota, Glomeromycota, Zoopagomycota, Entomophthoromycota, and the classes of Ascomycota and Basidiomycota.

**Sister pairs (examples):**

| pair | age (Mya) |
|---|---|
| *Gossypium darwinii* / *G. barbadense* | 0.4 |
| *Oryza sativa* / *O. glaberrima* | 0.9 |
| *Anopheles gambiae* / *A. arabiensis* | 1.3 |
| *Leishmania infantum* / *L. donovani* | 1.4 |
| *Neurospora crassa* / *N. tetrasperma* | 1.9 |
| *D. simulans* / *D. sechellia* | 2.1 |
| human / chimpanzee | 6.5 |
| *D. melanogaster* / *D. simulans* | 6.2 |

Further pairs: human/chimpanzee/bonobo/gorilla, two *Macaca*, mouse/*Mus caroli*, *Arabidopsis thaliana*/*A. lyrata*, *Brassica rapa*/*B. oleracea*, *C. elegans*/*briggsae*/*nigoni*, *Brugia*, *Plasmodium falciparum*/*P. reichenowi*, *Theileria*, *Candida albicans*/*C. dubliniensis*, *Cryptococcus*, *Aspergillus*, *Phytophthora*, *Thalassiosira*, *Entamoeba*, *Dictyostelium*, *Naegleria* and *Trypanosoma*.

**Not available:**

| species or pair | reason |
|---|---|
| rat (QfO) | exact coverage 0.555, because its 2022+ proteome has new accessions |
| *S. paradoxus*, *Mus spretus* | no reference proteome |
| *S. mikatae*, *Solanum pennellii*, *Hammondia*, *Neospora* | v4 coverage near 0 (proteomes newer than 2021) |
| *Cryptosporidium hominis* | v4 coverage 0.58 |
| *Schistosoma* | BUSCO below 80% |
| zebrafish, chicken, *Xenopus*, cow, maize, *Branchiostoma floridae* (all QfO) | exact coverage 0.08–0.72 |
| dog (QfO) | proteome no longer in the release |
| *Galdieria* pair, *Cryptosporidium* pair | TimeTree treats the two species as the same taxon |

### Prokaryotes (100), a separate set

- **Bacteria (60):** the 22 QfO bacteria in the release, plus case-study species (*Caulobacter*, *Myxococcus*, *Nostoc*, *Gemmata*, ***Prosthecobacter dejongeii*** (BtubA/B), *Thermus*, *Rickettsia*, *S. aureus*), plus one species per further phylum. The phylum round-robin covers 50 phyla, including many single-species candidate phyla and one CPR (Patescibacteria) genome.
- **Archaea (40):**
  - QfO archaea.
  - Case studies: *Haloferax volcanii* (CetZ), *Pyrobaculum calidifontis* and *Thermoproteus tenax* (crenactin), *Sulfolobus*, *Thermoplasma*, *Nanoarchaeum*, *Aeropyrum*, *Archaeoglobus*, *Pyrococcus*, *Methanobrevibacter*.
  - Thermoproteota 8, Nitrososphaerota 7, Methanobacteriota 8, Thermoplasmatota 7, Halobacteriota 2, and 5 DPANN.
  - **Asgard: only *Prometheoarchaeum syntrophicum* MK-D1.** The only other Asgard reference proteome, *Lokiarchaeum ossiferum* (2022), has no AFDB v4 models.
- BUSCO is relaxed to 20% for Asgard, DPANN and CPR. Asgard's coverage gate is relaxed to 0.5, but MK-D1 reaches ≥ 0.8 anyway.

## 2. Mapping to AFDB v4

- **Lookup:** the FoldComp `afdb_uniprot_v4.lookup` has 214,683,829 entries. All are `AF-<acc>-F1-model_v4` (no F2+ fragments).
- **Exact sequence check without fetching:** the TED per-chain table lists md5(sequence) for every v4 model. I verified this on a decoded entry. Over 3.57M proteins, the TED md5 match and the md5 of the FoldComp-decoded sequence **disagree 0 times**.
- **Missing models:**
  - Almost all proteins longer than 2,700 aa lack a model: 16,160 of 16,172 eukaryotic ones. These are titin- and dynein-type giants, so the panel is biased against very long proteins.
  - The other 154k missing proteins (4.5%) are accessions created after 2021_04.
- **Changed models:** 6,656 v4 models have a sequence that differs from the 2026 sequence; 78% of these have the same length (point revisions).

## 3. Structures

- **Fetch:** `fetch_foldcomp.py` uses HTTP range requests on S3 with 8 processes × 96 threads. It ran at about 4,400 entries/s: 3.29M entries (20.4 GB) in 12.3 min with 0 failures. Prokaryotes: 0.28M entries in 1.8 min.
- **Resumability:** the fetch writes per-shard data plus done-lists, then merges them into one FoldComp DB sorted by key. The shards were deleted after the merge was verified.
- **Extraction** (`extract_structures.py`): 3.29M entries decoded with 0 errors in 6.7 min on 32 processes. Outputs:
  - `*_proteins.parquet`: v4 sequence, length, md5, mean pLDDT, and fractions of residues with pLDDT ≥ 70 and ≥ 90.
  - `*_plddt.npz`: per-residue pLDDT in float16. FoldComp stores B-factors lossily, at about 0.1 precision.
  - `v4_fasta/<leaf>.fasta`: v4 model sequences, if a strictly v4-consistent sequence set is needed.
- **Foldseek DB:** `foldseek createdb <FoldComp DB> out` works directly (format auto-detected; `--input-format 5` is the explicit option). 3,290,300 entries in 10.4 min, peak 24 GB; 0 ignored.

## 4. TimeTree ages (`timetree_pairwise.py`)

Run with API version 20260820 on 2026-10-30. The ages file is `data/interim/timetree/pairwise_ages.tsv` and is **not committed**.

- **Age definition:** `age_mya` = `adjusted_age` (the TimeTree estimated time) when it is above 0; otherwise the median of studies (173 pairs). The median of studies, its CI and the number of studies are kept as separate columns.
- **Range:** 0.38 Mya (*Gossypium*) to 1,516 Mya, the TimeTree eukaryote root.

The age distribution is **not uniform when counted per pair, but roughly uniform per distinct node**:

| bin (Mya) | 1–2 | 2–4 | 4–8 | 8–16 | 16–32 | 32–63 | 63–126 | 126–251 | 251–501 | 501–1000 | 1000–1995 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| pairs | 9 | 8 | 14 | 23 | 44 | 32 | 120 | 113 | 458 | 2,491 | 12,443 |
| distinct ages (nodes) | 8 | 8 | 11 | 7 | 9 | 10 | 9 | 10 | 32 | 23 | 11 |
| nearest-neighbour age per species | 15 | 16 | 19 | 8 | 8 | 11 | 8 | 11 | 35 | 29 | 20 |

- There are 138 distinct ages and 117 distinct NCBI MRCA nodes. The KS distance from log-uniform is 0.28 for distinct ages and 0.82 for pairs.
- **Pairs are dominated by deep splits.** 79% of pairs are ≥ 1 Gya, and 8,172 pairs share the single root age of 1,515.8 Mya. TimeTree resolves deep eukaryote splits crudely: many supergroup pairs share its root age of 1,516 Mya, and another large set of deep pairs shares 1,164 Mya (2,760 pairs).
- Only 40 pairs are under 10 Mya and 168 are 10–100 Mya.
- Analyses should therefore weight by node, or subsample pairs per node. Raw pair-level regression would be swamped by root pairs.
- **No answer for 355 pairs:**
  - *Fonticula alba* vs all opisthokonts (109).
  - *Giardia*, *Trichomonas* and *Tritrichomonas* vs most others (203). TimeTree has no metamonad data for these.
  - A few individual pairs (Ustilaginomycotina vs Pucciniomycotina; *Xylona*, *Chloropicon*, *Pycnococcus*).
  - These need fixed literature ages (e.g. the eukaryote root, about 1.5–1.8 Gya) or exclusion. They are listed in `results/panel/timetree_age_summary.json`.

## 5. Sanity checks (`results/panel/sanity_checks_{panel_v4,panel_prok_v4}.json`)

- **Random 20 proteins per set (seed 20261030):** the live UniProt REST sequence equals the local 2026 FASTA in 20/20, and equals the FoldComp model sequence in 20/20, for both sets. The first run had 3 failures, all transient HTTP 500 errors from UniProt; I added retries.
- **Mismatch rate** (model differs from the 2026 sequence) among fetched models: 0.20% for eukaryotes and 0.25% for prokaryotes. Two independent checks (TED md5 and FoldComp sequence) agree on all 3.57M proteins.
- **Completeness:** every protein with a v4 lookup entry was fetched and decoded (0 missing).
- **Per-species proteins:**

  | set | min | median | max |
  |---|---|---|---|
  | eukaryotes | 2,787 | 15,661 | 82,543 |
  | prokaryotes | 483 | 2,248 | 8,253 |

  The eukaryotic maximum is a polyploid plant. Every species is at or above 0.80 exact-match coverage.
- **Lowest exact-match coverage:** *Chlamydomonas reinhardtii* 0.831, *P. reichenowi* 0.841, chimpanzee 0.848, *A. gambiae* 0.850.
- **Lowest mean pLDDT:** *Brachionus* 57, *Thalassiosira oceanica* 60, *Glossina* 62.

## 6. Problems and caveats for the PI

1. **Polyploid and hybrid genomes are in the panel:**
   - allotetraploid cottons: *Gossypium barbadense*, *G. darwinii*, *G. mustelinum*, *G. tomentosum* (AD genomes)
   - goldfish *Carassius auratus*, and tetraploid *Sinocyclocheilus* cavefish (2 or 3 species)
   - *Glycine max* / *G. soja*
   - mesopolyploid *Brassica*

   Their homeologues add many in-paralogs and "orthologs" between subgenomes. The 0.38 Mya *Gossypium* pair is an inter-species pair, but each genome carries A and D homeologues that diverged about 5–10 Mya. I did not exclude polyploids because Ancestra had already started. **Recommendation:** flag these species in downstream analyses and exclude them from shallow-pair rate estimates, or rerun Ancestra without them.
2. The greedy fill favoured vertebrate and plant genera for the 1–250 Mya bins (primates, bats, cichlids, cyprinids, *Gossypium*, *Glycine*, *Oryza*). Shallow pairs outside vertebrates and plants are *Drosophila*, *Anopheles*, *Caenorhabditis*, *Brugia*, *Plasmodium*, *Theileria*, *Leishmania*, *Trypanosoma*, *Naegleria*, *Neurospora*, *Candida*, *Cryptococcus*, *Aspergillus*, *Bipolaris*, *Phytophthora*, *Thalassiosira*, *Entamoeba* and *Dictyostelium*. This allows between-lineage rate comparisons, but each such clade has only 1–3 pairs.
3. Ten QfO eukaryotes are absent, mostly because of low v4 coverage (rat, zebrafish, chicken, *Xenopus*, cow, maize, *B. floridae*), plus 3 proteomes no longer in the release. The panel is therefore not the QfO benchmark set.
4. The TimeTree root age (1,516 Mya) is applied to 52% of dated pairs. Within-root structure is lost. Metamonads and *Fonticula* are undated against most species.
5. AFDB v4 has no models for proteins longer than 2,700 aa (16k proteins). Analyses of very long proteins (dyneins, titins, polyketide synthases) cannot use these structures.
6. `protein_table.parquet` is 190 MB. That is too large to commit; it is reproducible from the scripts. Committed files: `species_panel*.tsv`, `results/panel/*.tsv|json` (aggregate only), and `scripts/panel/*`. TimeTree-derived per-pair data stay in `data/interim/` and are not committed.

## 7. Reproduce

```bash
source /opt/conda/etc/profile.d/conda.sh && conda activate conda_envs/main
python scripts/panel/afdb_v4_tables.py; python scripts/panel/ted_md5_table.py
curl -s -o data/raw/proteomes/meta/refprot_meta.tsv "https://rest.uniprot.org/proteomes/stream?query=reference:true&fields=upid,organism,organism_id,protein_count,busco,lineage,cpd&format=tsv"
python scripts/panel/select_species.py --set euk-candidates
bash scripts/panel/download_proteomes.sh data/interim/panel/euk_candidates.tsv
python scripts/panel/map_proteomes_to_afdb.py
python scripts/panel/select_species.py --set euk --n-euk 180 --min-exact 0.80
# prokaryotes: iterate select(prok) -> download -> map until no new downloads (converged after 5 rounds)
python scripts/panel/map_proteomes_to_afdb.py --fasta-only --ancestra-panel data/processed/panel/species_panel.tsv --ancestra-out data/interim/ancestra_input
python scripts/panel/fetch_foldcomp.py --panel data/processed/panel/species_panel.tsv --name panel_v4
python scripts/panel/extract_structures.py --name panel_v4 --panel data/processed/panel/species_panel.tsv
foldseek createdb data/interim/structures/panel_v4 data/interim/foldseek/panel_db --threads 24
python scripts/panel/timetree_pairwise.py
python scripts/panel/sanity_checks.py panel_v4
```

The metadata snapshot is dated (2026-10-30). The UniProt `current_release` will move on, so pin to the 2026_03 FTP paths or keep `data/raw/proteomes/`.
