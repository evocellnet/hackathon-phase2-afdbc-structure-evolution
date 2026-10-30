# Data scout report: AFDB clusters, structures, sequences, taxonomy, divergence times and domains

Author: data-scientist agent. Date: 2026-10-30, about 09:25–09:55 CET on node 1.
Scripts: `scripts/download_afdb_clusters.sh`, `scripts/download_domain_annotations.sh`, `scripts/afdb_cluster_stats.py`, `scripts/afdb_v3_v6_overlap.py`, `scripts/timetree_coverage.py`.
Python environment for these scripts: `.venv_datascout/` (polars 2.1, pyarrow, foldcomp). It is separate from the shared conda environment.
Machine-readable results are in `results/data_scout/*.json` and `*_taxspan.parquet`. Provenance is tracked in `DATA.md`.

---

## 0. Executive summary

* **The AFDB Clusters bucket has two releases.**
  * **v3** is the paper release, built on AFDB v4 with 214.7M entries. It has 2,302,908 non-singleton clusters.
  * **v6** is a full re-clustering of AFDB v6 (UniProt 2025_03, 241M entries). It has 2,597,930 clusters. It was added in March 2026 and updated in August 2026.
* **All core tables are downloaded** to `data/raw/afdb_clusters/{v3,v6}/` (19.2 GB compressed in total):
  * membership tables (entry → representative → taxId)
  * cluster statistics (pLDDT, LCA, dark status)
  * all-vs-all E-values between representatives
  * AFDB50 mapping
  * **v3 only: 30.0M member-vs-representative Foldseek alignments with fident, lDDT and TM-score.** This file (file 8) answers the "sequence identity vs structural similarity" question of Elofsson 2009 directly, at scale and without computing any alignments.
* **Use v3 as the primary dataset.** It matches the paper and has complete taxonomy (only 152 entries lack a taxId). It is also the only release with the alignment file.
  * v6 has a serious defect: **82.4M of its 241M entries have taxId = 0**, and 322k clusters have an LCA of 0. Also, 21,936 entries are listed twice.
  * Only 213k of the 2.30M v3 representatives are still representatives in v6. The two releases are different clusterings and should not be mixed.
* **Structures: the recommended route is FoldComp `afdb_uniprot_v4` with HTTP range requests.**
  * The data file is 1.0 TB, but it never needs to be downloaded whole. We download only the `.index` (5.6 GB) and `.lookup` (8.4 GB) files, then fetch individual entries by byte range from S3.
  * Tested: about 170 structures/s from one Python process with 128 threads, all 2,000 test fetches successful, decompressed to PDB.
  * Estimated cost: 1M structures ≈ 1.6 h single-process, or ≈ 20–30 min with 4–8 processes, using ≈ 5 GB of disk (FoldComp format).
  * For all-vs-all Foldseek work on AFDB50, download the prebuilt Foldseek DB `afdb50` (114 GB, v6).
  * Per-entry download from the AFDB website is not usable at scale: the server allows a burst of about 200 requests, then returns HTTP 429 for about 5 minutes. Only v6 files are served. The Google Cloud bucket now denies anonymous reads of objects.
* **Taxonomy is downloaded.** NCBI taxdump and new_taxdump were downloaded on 2026-10-30, and the md5 checks passed.
* **Divergence times:**
  * The TimeTree of Life 2015 newick is downloaded. It has 50,641 species, but these are 99% eukaryotes (only 152 bacteria and 25 archaea).
  * TimeTree 6 has a working, undocumented REST API: `/api/pairwise/{taxid}/{taxid}` returns a CSV of the median divergence time with its confidence interval. It runs at about 17 queries/s with 8 parallel requests. It returns nothing for very deep splits (e.g. human vs E. coli, Bacteria vs Eukaryota), so the deepest nodes must be fixed from the literature.
  * Licence: TimeTree data are for personal research and teaching only and must not be redistributed.
* **Domains:**
  * Pfam 38.2 `Pfam-A.regions.tsv.gz` (4.7 GB, built on UniProtKB 2025_03) is downloaded, and its md5 matches the official checksum.
  * The TED v1 consensus domain summary (19.9 GB) and per-chain segmentation (9.6 GB) from Zenodo are downloaded, and their md5s match Zenodo.
  * InterPro 110.0 `protein2ipr.dat.gz` (12 GB) is not downloaded; it is optional.

---

## 1. Contents of the AFDB Clusters bucket

Index page: https://afdb-cluster.steineggerlab.workers.dev/ redirects to https://steineggerlab.s3.amazonaws.com/afdb-cluster/.
The S3 listing (`?prefix=afdb-cluster/`) shows 14 objects and is not truncated. A snapshot is saved in `data/raw/afdb_clusters/s3_listing.xml`.
All objects were re-uploaded on 2026-10-07.
Licence: all files are under CC-BY-4.0. Cite Barrio-Hernandez et al., Nature 2023.

The index README lists these updates:
* 2023-03-19: file 5 added.
* 2023-05-25: fix to the fragment merge.
* 2024-12: files 6 and 7, with representative IDs substituted.
* 2025-09-12: file 8 added.
* 2026-03-04: v6 release.
* 2026-03-10 to 2026-08-10: several cluFlag, LCA and repLen fixes to v6.

### 1.1 v3 files (paper release; AFDB v4 = UniProt 2021_04)

| # | file | size | rows | columns (no header) | content |
|---|---|---|---|---|---|
| 1 | 1-AFDBClusters-entryId_repId_taxId.tsv.gz | 274.5 MB | 30,045,247 | memberId, repId, taxId | Members of the 2.30M non-singleton, fragment-free Foldseek clusters. The members are AFDB50 representatives. |
| 2 | 2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz | 41.0 MB | 2,302,908 | repId, isDark (0/1), nMem (members without fragments), repLen, avgLen, repPlddt, avgPlddt, LCAtaxId | One row per cluster. |
| 3 | 3-sapId_sapGO_repId_cluFlag_LCAtaxId.tsv.gz | 2.2 MB | – | human acc, GO terms (`;`-separated), repId, cluFlag, LCA taxId | Every human protein, with its cluster and GO terms. |
| 4 | 4-domain-clustering.zip | 4.9 MB | – | `Structure_groups_info.csv` (group_ID, connected_group, counts, top-2 Pfam and their frequencies); `Structures_from_groups.csv` (structure_ID = acc_domainIdx, start, end, group_ID, uniprot, pfam_all, pfam_overlap_rel); `Headers.xlsx` | The paper's domain-level clustering. Group "1;1" is the 153k-structure kinase / zinc-finger supergroup. |
| 5 | 5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz | 1.68 GB | 214,684,311 | repId, entryId, cluFlag, taxId | All AFDB v4 entries. cluFlag: 1 = AFDB50 member (162.4M), 2 = member of an AFDB cluster (30.0M), 3 = fragment (9.3M), 4 = singleton (13.0M). |
| 6 | 6-all-vs-all-similarity-queryId_targetId_eValue-Dec16_2024.tsv.gz | 4.82 GB | 600,195,809 | queryId, targetId, E-value | Foldseek all-vs-all search between cluster representatives. **E-value only: no TM-score or identity.** |
| 7 | 7-AFDB50-repId_memId.tsv.gz | 1.27 GB | 214,684,311 | repId, memId (third column empty) | MMseqs2 AFDB50 clustering of all entries. |
| 8 | 8-cluster_alignment_tmscore_lddt.tar.gz | 994 MB (3.9 GB uncompressed) | 30,045,247 | query(rep) target(member) fident alnlen mismatch gapopen qstart qend tstart tend evalue bits lddt alntmscore. IDs look like `AF-<acc>-F1-model_v3.cif`. Also two per-cluster files (repId, sum, nMem, mean) for lDDT and TM-score. | **A Foldseek alignment of every cluster member against its representative.** |
| – | commands.gz | 2.9 KB | – | shell commands | Pipeline: fragment removal → singleton removal → LCA (mmseqs lca) → dark labelling (Pfam via MMseqs, PDB via Foldseek, Pfam/TIGRFAM via FTP) → file assembly. |

Quick look at file 8. We excluded the 2.30M self-alignments, leaving 27.7M member–representative pairs.

| fident bin | n pairs | median TM-score | median lDDT |
|---|---|---|---|
| < 0.2 | 12.8M | 0.505 | 0.520 |
| 0.2–0.3 | 8.25M | 0.706 | 0.683 |
| 0.3–0.4 | 3.84M | 0.799 | 0.777 |
| 0.4–0.5 | 1.66M | 0.827 | 0.830 |
| 0.5–0.7 | 0.96M | 0.826 | 0.863 |
| 0.7–0.9 | 0.20M | 0.830 | 0.899 |
| ≥ 0.9 | 0.05M | 0.865 | 0.937 |

Caveats for file 8:
* The members are AFDB50 representatives, so the pairs are **depleted of pairs above 50% sequence identity** by design.
* The fident value comes from a Foldseek alignment (3Di + amino acids), not from a sequence-only alignment.
* The Foldseek clustering used 90% coverage and E < 0.01, so pairs below TM ≈ 0.3–0.5 are truncated (selection bias).
* To cover the high-identity end, align AFDB50 members (cluFlag 1) ourselves.
* To cover the low-similarity end, use the all-vs-all hits in file 6 and realign them.

### 1.2 v6 files (AFDB v6 = UniProt 2025_03)

| # | file | size | rows | columns | notes |
|---|---|---|---|---|---|
| 1 | 1-AFDBClusters-repId_entryId_cluFlag_taxId.tsv.gz | 1.65 GB | 199,373,567 | repId, entryId, cluFlag (1 = Foldseek-clustered, i.e. AFDB50 rep: 43.1M; 2 = AFDB50 member: 156.3M), taxId | **Column order differs from v3 file 1, and the flag meanings differ too.** |
| 2 | 2-repId_…_LCAtaxId.tsv.gz | 49.9 MB | 2,597,930 | same as v3 | 1,080,862 dark clusters. 322k clusters have LCAtaxId = 0. |
| 5 | 5-allmembers-…tsv.gz | 2.08 GB | 241,070,489 (241,030,641 unique) | repId, entryId, cluFlag (1/2), taxId | **82,392,903 rows have taxId = 0.** About 59.4M of these can be filled from v3; the rest need AFDB metadata (§4). |
| 6 | 6-all-vs-all-…tsv.gz | 4.39 GB | 552,702,072 | query, target, E-value | |
| 7 | 7-AFDB50-repId_memId.tsv.gz | 1.79 GB | 241,070,489 | repId, memId as `AF-<acc>-F1-model_v6` | |

v6 has **no** human file, no domain file and **no alignment (TM-score/lDDT) file**.

Overlap between v3 and v6 (`results/data_scout/afdb_v3_v6_overlap.json`):
* 175.4M of the 214.7M v3 entries (81.7%) are still in v6.
* 152.6M of the 188.3M members of v3 non-singleton clusters are still in v6.
* 2.03M of the 2.30M v3 representatives still exist as AFDB v6 entries, but only 213k are still cluster representatives in v6.

This matches the AFDB changelog: v6 added 65.7M entries, removed 39.3M, changed 277k and left 175.1M unchanged.

### 1.3 Checksums

The local md5 values are in `data/raw/afdb_clusters/md5sums.txt` and are copied into `DATA.md`.
For files 3, 4 and `commands.gz` (single-part S3 ETags), the ETag equals the local md5, which confirms those downloads.
The other files were uploaded in multiple parts, so their ETags are not md5s. They decompressed without errors.

---

## 2. Quick statistics of the cluster tables

Script: `scripts/afdb_cluster_stats.py`.
Taxonomy comes from NCBI new_taxdump, with merged taxIds resolved. "Domain" here means the three cellular domains plus viruses.
"Other" covers taxa outside these four groups, such as unclassified or metagenome taxa, and LCA = root ("cellular organisms").
The full numbers are in `results/data_scout/afdb_cluster_stats_{v3,v6}.json`.

| statistic | v3 | v6 |
|---|---|---|
| AFDB entries (file 5) | 214,684,311 | 241,070,489 rows (241.03M unique) |
| entries with taxId missing (0) | 152 | 82,392,903 |
| distinct taxIds | 1,006,700 | 1,033,731 |
| non-singleton clusters (file 2) | 2,302,908 | 2,597,930 |
| dark clusters | 711,705 (30.9%) | 1,080,862 (41.6%) |
| sum of nMem (cluster members that are AFDB50 representatives) | 30,045,247 | 42,139,046 |
| all entries mapped to these clusters (incl. AFDB50 members, plus fragments in v3) | 188,252,820 | 199,373,567 |
| nMem median / 90th / 99th / 99.9th percentile / max | 3 / 19 / 153 / 932 / 38,901 | 3 / 19 / 194 / 1,404 / 90,299 |
| clusters with nMem 2–4 / 5–9 / 10–19 / 20–49 / 50–99 / 100–999 / ≥1000 | 1.43M / 432k / 220k / 138k / 46k / 36k / 2.1k | 1.65M / 466k / 235k / 147k / 53k / 48k / 4.0k |
| avgPlddt median | 74.0 | 71.4 |
| clusters with avgPlddt ≥ 70 / ≥ 90 | 1.37M / 170k | 1.38M / 166k |
| LCA = cellular organisms (root) | 529,373 | 279,782 |
| LCA within Bacteria / Eukaryota / Archaea | 896k / 832k / 38k | 718k / 1.25M / 30k |
| clusters whose members (incl. AFDB50 members) span ≥ 2 cellular domains | 605,881 | 279,807* |
| …all 3 domains (Bacteria + Archaea + Eukaryota) | 103,776 | 33,622* |
| ≥ 2 NCBI kingdoms (bacteria now have kingdoms, e.g. Pseudomonadati and Bacillati, so this is inflated) | 1,302,998 | 808,570* |
| ≥ 5 phyla | 486,389 | 248,509* |
| clusters with ≥ 1 member from the 16 model organisms below | 152,380 | 133,565 |
| clusters with ≥ 3 different model organisms | 25,542 | 23,345 |
| clusters containing both human and E. coli K-12 | 332 | 513 |
| clusters containing both human and S. cerevisiae | 1,093 | 1,388 |

\* The v6 span numbers are underestimates, because a third of v6 entries have taxId 0.

Other v3 subsets:
* 381,584 clusters span at least 2 domains and have avgPlddt ≥ 70.
* 231,545 clusters have nMem ≥ 10, avgPlddt ≥ 70 and at least 5 phyla.
* The clusters spanning all 3 domains have median nMem 22 and median avgPlddt 83.

Model organisms in v3 (human 9606, mouse 10090, rat, zebrafish, fly, worm, Arabidopsis, S. cerevisiae 559292, S. pombe, Dictyostelium, E. coli K-12 83333, M. jannaschii, C. albicans, maize, rice, soybean). The table gives the number of clusters containing the organism, then the number of its proteins inside clusters out of all its AFDB proteins:

| organism | clusters | proteins in clusters / all AFDB proteins |
|---|---|---|
| human | 33,987 | 145,625 / 186,020 |
| mouse | 25,014 | 59,096 / 79,301 |
| Arabidopsis | 23,714 | 113,752 / 132,903 |
| S. cerevisiae | 4,567 | 5,480 / 6,720 |
| E. coli K-12 | 3,349 | 4,272 / 4,584 |
| M. jannaschii | 1,466 | 1,612 / 1,773 |

Only 15,031 clusters contain both human and mouse proteins. Human–Arabidopsis: 2,767; human–fly: 3,732.
In other words, these Foldseek clusters (90% bidirectional coverage) are **not orthogroups**. They split multi-domain proteins and lineage-specific architectures.
For orthology we must use Ancestra or another orthology tool, not cluster membership.

---

## 3. Structures: how to get coordinates for 100k–2M proteins

| option | what | size | version | measured speed | verdict |
|---|---|---|---|---|---|
| (a) AFDB per-entry files `https://alphafold.ebi.ac.uk/files/AF-<acc>-F1-model_v6.{pdb,cif,bcif}` | full-atom PDB/mmCIF, plus PAE, confidence JSON and (new in v6) the MSA in a3m format | ~50–100 KB per entry | **only v6 is served.** v1–v5 return 404. Requests without a User-Agent header get 403. | A burst of about 200 requests at up to 30/s succeeds, then **HTTP 429 for about 5 minutes** (tested twice). That leaves a sustained rate below 1 request/s, so 1M entries would take more than 10 days. | Use only for a few thousand entries, or for PAE/MSA of a small selected set. |
| (b) FoldComp DBs, https://foldcomp.steineggerlab.workers.dev/ (S3 `steineggerlab/foldcomp/`) | Compressed full-atom models (≈ 4.7 KB per entry); the sequence is recoverable from them | `afdb_uniprot_v4` 1.0 TB (+ .index 5.6 GB + .lookup 8.4 GB); `afdb_rep_v4` 9.5 GB (2.3M AFDB-cluster reps); `afdb_rep_dark_v4` 2.2 GB; `afdb_swissprot_v4` 2.8 GB; `highquality_clust30` 112 GB; `esmatlas` 1.7 TB | **v4.** This matches the v3 clusters, and v4 also contains the 39M entries removed in v6. | **HTTP Range requests work** (S3 returns 206). Using `afdb_swissprot_v4.index`, 2,000 random entries took 13.9 s with 32 threads (72/s) and 6.0 s with 128 threads (167/s). All succeeded and decompressed with `foldcomp.decompress`. Bulk S3 download runs at roughly 2–40 MB/s per stream and about 140 MB/s with 6 streams, varying with load. | **Recommended for coordinates.** Download the index and lookup (14 GB, about 5 min), then range-fetch the selected entries. 1M entries ≈ 1.6 h with one process, ≈ 20–30 min with 4–8 processes. Store the results as a local FoldComp DB, which `foldseek createdb` and `foldcomp` read directly. `afdb_rep_v4` (9.5 GB) gives all 2.3M representatives in one download. |
| (c) Foldseek prebuilt DBs (`foldseek databases`; files at https://foldseek.steineggerlab.workers.dev/) | Foldseek DB: amino acids, 3Di, Cα coordinates, headers, `_mapping` (entry → taxId) and `_taxonomy` (verified by unpacking `afdb_proteome`) | `afdb` (Alphafold/UniProt) **457.7 GB**; `afdb50` (Alphafold/UniProt50) **114.1 GB**; `afdb50clusearch` 346.9 GB (cluster-search variant); `afdb_proteome` 1.7 GB (downloaded: 590,435 structures from 48 proteomes); `afdb_swissprot` 1.5 GB; `teddb` 320 GB and `teddb_afdb50` 59.5 GB (TED domains, AFDB v4); `pdb100` 2.2 GB; `cath50` 0.97 GB | AFDB DBs are **v6 2025_03** | One stream reached about 40 MB/s, so afdb50 would take about 50 min single-stream or about 15–20 min with parallel ranges. | **Recommended for any Foldseek search or all-vs-all** (no need to run `createdb` on millions of PDBs). Use `afdb50` (114 GB) if we search beyond the cluster representatives; this needs PI approval because it is over 50 GB. No full-atom coordinates (Cα only), which is enough for TM-score and lDDT in Foldseek. Note: v6 IDs, while the v3 clusters use v4. |
| (d) Google Cloud `gs://public-datasets-deepmind-alphafold-v4` | per-entry cif/json files, proteome tar files (`proteomes/proteome-tax_id-*-0_v4.tar`), `accession_ids.csv`, `sequences.fasta` (99 GB) | ~23 TB | v4 (plus v3) | **Listing works anonymously, but object reads return 403 ("Anonymous caller does not have storage.objects.get")**, even for objects that are listed. | Would need an authenticated Google Cloud account (likely requester-pays). Not needed. |
| (e) EBI FTP https://ftp.ebi.ac.uk/pub/databases/alphafold/ | `latest/` and `v6/` proteome tar files (model organisms and global-health set); `swissprot_pdb_v6.tar` 28.6 GB; `v4/` archives | 0.2–5 GB per proteome | v6 (plus older versions) | About 5.6 MB/s single stream (500 MB in 90 s) | Fine for a few proteomes. Slow for bulk downloads. |

**Recommendation for structures:**
1. **Download now** (14 GB, under the 50 GB limit): the FoldComp `afdb_uniprot_v4.index` and `.lookup`, plus `afdb_rep_v4` (9.5 GB).
2. Write `scripts/fetch_foldcomp_subset.py`, which looks up the selected accessions and range-fetches them in parallel into a local FoldComp DB.
3. Run Foldseek (`createdb` on the FoldComp subset) for pairwise TM-score, lDDT and 3Di alignments, and alignment-free 3Di comparisons, within the chosen sets (orthogroups or clusters).
4. Request the 114 GB `afdb50` Foldseek DB only if we need searches against the whole AFDB50.

---

## 4. Sequences

| source | size | content | notes |
|---|---|---|---|
| AFDB `sequences.fasta`, https://ftp.ebi.ac.uk/pub/databases/alphafold/sequences.fasta | **110 GB** (2026-02-12) | One line per sequence. Header: `>AFDB:AF-<acc>-F1 <name> UA= UI= OS= OX=<taxid> GN=` | v6. Takes about 5.5 h single-stream from EBI, or under 1 h with about 16 parallel ranges. |
| AFDB release metadata, `release/metadata/gdm_2025_09.csv` | **95 GB** | modelEntityId, uniprotAccession, providerId, globalMetricValue (mean pLDDT), sequenceStart, sequenceEnd, **taxId**, organism, **sequence**, fractionPlddtLow, fractionPlddtConfident, fractionPlddtVeryHigh | **The best single table for v6**: sequence, taxId and pLDDT for each entry. It also fixes the v6 taxId = 0 problem. We could stream and parse it without storing it, keeping only the columns we need (a few GB). Other files: `atbc_2026_02.csv` (5 GB, new provider "ATBC" models of metagenome-derived proteins) and `bfvd/ntdx/vr3d` (small). |
| `accession_ids.csv` (EBI FTP) | 8.1 GB | acc, first residue, last residue, AFDB ID, latest version | |
| `diffs.ndjson.gz` | 294 MB | ADDED / REMOVED / CHANGED between v4 and v6 | Useful to see which v3 members changed. |
| UniProtKB 2026_03: `uniprot_trembl.fasta.gz` 38 GB, `uniprot_sprot.fasta.gz` 89 MB | | | **Newer than AFDB.** Some v3 accessions are obsolete: UniProt REST returned nothing for 4 of 4 tested v6 accessions that lack a taxId. |
| FoldComp entries (option 3b) | – | amino-acid sequence embedded in each entry | Free if we fetch structures anyway. |
| Foldseek `afdb50` DB | 114 GB | amino-acid sequences plus 3Di of the AFDB50 representatives | If downloaded for search, it also gives sequences. |

**Recommendation for sequences:** take sequences of the selected proteins from the FoldComp entries that we fetch anyway.
For the full-scale v3 sequence work, download `sequences.fasta` (110 GB; needs PI approval) with parallel range requests. Alternatively, stream `gdm_2025_09.csv` once and keep only accession, taxId, pLDDT, length and sequence for cluster members.

---

## 5. Taxonomy and divergence times

### NCBI Taxonomy (downloaded)

* Files: `data/raw/ncbi_taxonomy/taxdump.tar.gz` (80 MB) and `new_taxdump.tar.gz` (163 MB), from the FTP on 2026-10-30. The md5 checks passed and both are extracted.
* There are 3,024,644 nodes.
* Note: NCBI now uses the rank "domain" (and "realm" for viruses) instead of "superkingdom", and bacteria now have "kingdom" ranks.
* `rankedlineage.dmp` columns: taxId, name, species, genus, family, order, class, phylum, kingdom, domain/realm.
* Licence: public domain.

### TimeTree

TimeTree 6 was released in August 2026; the API reports version 20260820.

* **Bulk downloads**:
  * `https://timetree.org/public/data/TimetreeOfLife2015.nwk` (Hedges et al. 2015), downloaded to `data/raw/timetree/`. It has 50,641 tips.
    * 49,928 tips map to NCBI scientific names.
    * Of these, 49,751 are eukaryotes (31,952 Metazoa, 16,267 Viridiplantae, 926 Fungi), 152 are bacteria and 25 are archaea.
  * TimeTree 6 clade trees (Mammals, Tetrapods, Amphibians, Reptiles, Gymnosperms) from `/downloads`, also downloaded. They contain vertebrates and plants only.
  * No TimeTree 5 or 6 tree of all life is offered for download (e.g. `/downloads/TimeTree5.nwk` returns an HTML page).
* **Coverage of AFDB v3** (`results/data_scout/timetree_coverage.json`):
  * 40,826 of 966,981 AFDB species are in TTOL2015: 40,675 eukaryotes, 129 bacteria and 22 archaea.
  * These species hold 30.4M of the 214.7M proteins.
  * 592,761 clusters have at least 2 TTOL species, and 126,949 have at least 10.
  * In short, TTOL2015 is adequate for eukaryotic clusters and useless for prokaryotes at species level.
* **REST API** (undocumented but stable; it accepts NCBI taxIds):
  * `GET https://timetree.org/api/pairwise/{taxidA}/{taxidB}` returns CSV with these columns: `taxon_a_id, taxon_b_id, scientific_name_a, scientific_name_b, all_total (n studies), precomputed_age, precomputed_ci_low, precomputed_ci_high, adjusted_age`.
    * Examples: human–mouse 83.5 Ma (CI 80.1–86.9, 75 studies); human–Arabidopsis 1,455 Ma; E. coli–B. subtilis 1,909 Ma (adjusted 2,789).
    * Strain taxIds are mapped to species automatically (1310605 → 470).
    * **Empty results** for human–E. coli, Bacteria–Eukaryota, Archaea–Bacteria, human–Archaea and E. coli K-12–yeast.
  * `GET https://timetree.org/api/taxon/{taxid}` returns JSON with the per-study crown-age estimates (`study_times`). Examples: Mammalia; Eukaryota (1004–1781 Ma); Bacteria (≈ 3300–3435 Ma); Archaea (4147 Ma). Some nodes, such as Bacillota and "cellular organisms", have none.
  * Speed: about 0.33 s per call, or about 17 calls/s with 8 parallel requests. 100k pairs would take about 1.6 h.
* **Recommended strategy:** assign each pair of proteins the age of the NCBI node that is the MRCA of their two species.
  1. Collect the set of MRCA nodes needed (thousands, not millions).
  2. Query `/api/pairwise` once per node, using a representative species pair from either side of that node.
  3. Fix the deepest splits from the literature: LUCA ≈ 4.2 Ga (Moody et al. 2024 Nat Ecol Evol); Archaea–Eukaryota (Asgard) ≈ 2.7–3.5 Ga; LECA ≈ 1.6–1.9 Ga (Strassert et al. 2021).
  4. For prokaryotes, consider other timetrees with wider coverage: GTDB-based dated trees (Betts et al. 2018; Moody et al. 2024 supplementary), or `datelife` (R).
  * Use the CI or the number of studies as an uncertainty weight.
* **Licence**: "Data … provided openly for personal research and teaching use … Redistribution of TimeTree data and its transformations are not permitted."
  * Do not commit the TimeTree newicks or tables derived from the API to git. Commit only the script that regenerates them.
  * Cite Kumar et al. 2022 MBE (TimeTree 5) and Hedges et al. 2015.

---

## 6. Domain annotations

| dataset | URL | files and sizes | version | licence | status |
|---|---|---|---|---|---|
| **Pfam-A.regions** | https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.regions.tsv.gz | 4.7 GB. Columns: pfamseq_acc, seq_version, crc64, md5, pfamA_acc, seq_start, seq_end, ali_start, ali_end. Plus `Pfam-A.clans.tsv.gz` (541 KB) and `pfamA.txt.gz` (13 MB). | Pfam 38.2 (2026-01; UniProtKB 2025_03, the same as AFDB v6), 30,134 families | CC0 | **Downloaded** to `data/raw/pfam/`; md5 matches the official checksums |
| **TED** domain summary | https://zenodo.org/records/13908086 (newest of 4 versions, 2024-10-31) | `ted_365m.domain_summary.cath.globularity.taxid.tsv.gz` **19.9 GB**. Columns: `AF-<acc>-F1-model_v4_TEDxx`, md5, consensus level (high/medium), chopping (e.g. `11-41_290-389`), nres, n segments, domain pLDDT, packing and SSE metrics, proteome ID, CATH label (e.g. 3.40.640.10), CATH level (H/T/…), assignment method (foldseek/foldclass), norm_rg, globularity score, -, species, lineage. | TED v1, AFDB v4 | CC-BY-4.0 | **Downloaded** to `data/raw/ted/`; md5 matches Zenodo |
| TED per-chain segmentation | same record | `ted_214m_per_chain_segmentation.tsv.gz` 9.6 GB (AF ID, md5, nres, number of high/medium/low-consensus domains, choppings by consensus level, proteome) | TED v1 | CC-BY-4.0 | **Downloaded**; md5 matches Zenodo |
| TED other files | same record | `ted_324m_seq_clustering.cathlabels.tsv.gz` 3.3 GB (domain → sequence-cluster representative and CATH label); three per-method choppings at 7–8 GB each (chainsaw, merizo, unidoc); redundant-set summaries at about 1.7 GB | | CC-BY-4.0 | Not needed now |
| TED Foldseek DBs | foldseek.steineggerlab.workers.dev | `teddb_afdb50.tar.gz` 59.5 GB; `teddb.tar.gz` 320 GB | AFDB v4 | CC-BY-4.0 | Only if we need domain-level structural search |
| InterPro | https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/ | `protein2ipr.dat.gz` **12 GB** (acc, IPR, name, member DB signature, start, end); a parquet version `parquet/uniprotkb_matches_*.snappy.parquet` (26 × ~0.5 GB ≈ 14 GB); `interpro2go` 2.9 MB; `ParentChildTreeFile.txt`; `match_complete.xml.gz` 53 GB | InterPro 110.0 (2026-09-03) | CC0 | Not downloaded. Optional; Pfam plus TED is enough for domain architecture. |

Note that the v3 cluster file 4 already contains the paper's domain-level clustering, with Pfam labels.

---

## 7. Recommended acquisition plan and time estimates

| step | data | size | time | status |
|---|---|---|---|---|
| 1 | AFDB Clusters v3 + v6 (all files) | 19.2 GB | about 3 min | **done** |
| 2 | NCBI taxdump + new_taxdump | 0.24 GB | under 1 min | **done** |
| 3 | TimeTree 2015 + TT6 clade newicks | 4 MB | – | **done** |
| 4 | Foldseek `afdb_proteome` (model-organism 3Di/Cα DB) | 1.7 GB | about 1 min | **done** |
| 5 | Pfam 38.2 regions + clans | 4.7 GB | about 15 min at 5 MB/s | **done** (md5 verified) |
| 6 | TED domain summary + per-chain segmentation | 29.5 GB | about 5 min | **done** (md5 verified against Zenodo) |
| 7 | FoldComp `afdb_uniprot_v4.index` + `.lookup` (+ `afdb_rep_v4` 9.5 GB) | 14 (+9.5) GB | about 5–10 min | **recommend now** |
| 8 | Range-fetch selected structures (FoldComp v4) | about 5 KB per protein | about 170/s per process; 1M proteins in about 30 min with 6 processes | after protein selection |
| 9 | TimeTree node ages via the API for the needed MRCA nodes | KB | about 10 min per 10k nodes | after species selection |
| 10 (needs approval, over 50 GB) | Foldseek `afdb50` DB | 114 GB | about 20–50 min | only if we search the whole AFDB50 |
| 11 (needs approval, over 50 GB) | AFDB `gdm_2025_09.csv` (stream-parse, keep columns) or `sequences.fasta` | 95 / 110 GB | about 1–5 h from EBI | only if v6 or full-scale sequences are needed |

Suggested analysis design based on what is available:
1. **Linearity of structure versus sequence divergence (Elofsson 2009 revisited).** v3 file 8 provides 27.7M (fident, TM-score, lDDT) pairs, each with the taxIds of the representative and the member. Combine these with MRCA ages from TimeTree.
   * Supplement the high-identity range by realigning AFDB50 members (cluFlag 1) against their representatives with Foldseek on FoldComp subsets.
   * Restrict to pLDDT ≥ 70 (avgPlddt in file 2; per-entry pLDDT from FoldComp or the gdm metadata).
2. **Orthology.** Use Ancestra on proteomes; AFDB clusters are not orthogroups (§2). Then fetch the orthologs' structures from FoldComp v4 by accession.
3. **Domain architecture.** Use TED choppings (AFDB v4, which matches v3) and Pfam 38.2 regions.
4. **Keep v6 as a robustness check** for re-clustering effects only after filling its missing taxIds.

## 8. Licences

| source | licence |
|---|---|
| AFDB / AFDB Clusters / FoldComp / Foldseek DBs | CC-BY-4.0 |
| UniProt | CC-BY-4.0 |
| TED | CC-BY-4.0 |
| Pfam / InterPro | CC0 |
| NCBI Taxonomy | public domain |
| TimeTree | personal research and teaching use; no redistribution; cite Kumar et al. 2022 and Hedges et al. 2015 |

## 9. Gotchas found

* The AFDB website needs a User-Agent header (without one it returns 403). It rate-limits to about 200 requests per ~5 minutes (HTTP 429), and serves only v6 files.
* The GCS AlphaFold bucket allows anonymous listing but not anonymous reads.
* v3 and v6 cluster files use different column orders and different cluFlag meanings in file 1. v6 file 5 uses flags 1/2, while v3 uses 1–4 with different meanings.
* v6 has missing taxIds and duplicate entries (see §1.2).
* In polars ≥ 2.0, unnamed columns are `column_0…` (0-based); this caused one bug, which has been fixed.
* EBI FTP is slow per stream (about 5 MB/s); use parallel range requests for files over 10 GB.
* Pairs in file 8 are within-cluster only, so they are biased towards high TM-score, and towards identity below 50% because of the AFDB50 preclustering.
