# DATA.md — datasets used in this project

Each row lists a dataset we downloaded (or that is in progress or recommended), where it came from, its version, where it is stored locally, its checksum and its license.
Raw data lives under `data/raw/` and is **never committed to git**: the files are large, and TimeTree data may not be redistributed.
Full details and column descriptions are in `agents/data_scout_report.md`.
Download scripts: `scripts/download_afdb_clusters.sh` and `scripts/download_domain_annotations.sh`.

Checksum notes. "md5" is the local `md5sum` of the downloaded file. On S3, a single-part ETag equals the md5, and where that applies we confirmed the two match (files 3, 4 and `commands.gz`). Multipart ETags (format `xxx-N`) are not md5s, so for those files only the local md5 is recorded.

| dataset | source URL | version/date | local path | md5 | license | notes |
|---|---|---|---|---|---|---|
| AFDB Clusters v3 — file 1: AFDB-cluster members (entryId, repId, taxId) | https://steineggerlab.s3.amazonaws.com/afdb-cluster/v3/1-AFDBClusters-entryId_repId_taxId.tsv.gz | v3 (paper release; S3 upload 2026-10-07) | data/raw/afdb_clusters/v3/ | 62dcd980bb31b62a90a23cc09f728a12 | CC-BY-4.0 | 30.0M members (AFDB50 reps) of 2.30M clusters |
| AFDB Clusters v3 — file 2: cluster summary (repId, isDark, nMem, repLen, avgLen, repPlddt, avgPlddt, LCAtaxId) | .../v3/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz | v3 | data/raw/afdb_clusters/v3/ | a0243674279d84c02e3f1ffe65649b2b | CC-BY-4.0 | 2,302,908 clusters, of which 711,705 are dark |
| AFDB Clusters v3 — file 3: human proteins (sapId, GO, repId, cluFlag, LCA) | .../v3/3-sapId_sapGO_repId_cluFlag_LCAtaxId.tsv.gz | v3 | data/raw/afdb_clusters/v3/ | d43d9e14d714830646cfdb48730f7bff (= ETag) | CC-BY-4.0 | |
| AFDB Clusters v3 — file 4: domain clustering | .../v3/4-domain-clustering.zip | v3 | data/raw/afdb_clusters/v3/ (also extracted) | 29b914c58ac15d8be32b4d761e56519f (= ETag) | CC-BY-4.0 | Domain groups with Pfam annotation |
| AFDB Clusters v3 — file 5: all 214M AFDB v4 entries (repId, entryId, cluFlag, taxId) | .../v3/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz | v3 (AFDB v4) | data/raw/afdb_clusters/v3/ | 09bc2517fa3793ffd2884bff1dd3069f | CC-BY-4.0 | cluFlag 1 = AFDB50 member, 2 = AFDB-cluster member, 3 = fragment, 4 = singleton |
| AFDB Clusters v3 — file 6: all-vs-all similarity between representatives (E-value only) | .../v3/6-all-vs-all-similarity-queryId_targetId_eValue-Dec16_2024.tsv.gz | Dec 16 2024 | data/raw/afdb_clusters/v3/ | 06ecb926be335e9eaa27b68798249cfb | CC-BY-4.0 | 600,195,809 rows |
| AFDB Clusters v3 — file 7: AFDB50 (MMseqs2) repId → memId | .../v3/7-AFDB50-repId_memId.tsv.gz | Dec 13 2024 | data/raw/afdb_clusters/v3/ | 89d6be22a45f6df57d56012768a5de85 | CC-BY-4.0 | 214,684,311 rows |
| AFDB Clusters v3 — file 8: member-vs-representative alignments (fident, alnlen, E-value, bits, lDDT, TM-score) plus per-cluster averages | .../v3/8-cluster_alignment_tmscore_lddt.tar.gz | Sep 12 2025 | data/raw/afdb_clusters/v3/8-cluster_alignment/ | 798f35a111fb2adb0c4a39c849a9b302 | CC-BY-4.0 | 30,045,247 rows. **Key data for the sequence-vs-structure analysis** |
| AFDB Clusters v3 — commands.gz | .../v3/commands.gz | v3 | data/raw/afdb_clusters/v3/ | 1b46e24a461b1dd8270d7f3eabc5363e | CC-BY-4.0 | The authors' processing commands |
| AFDB Clusters v6 — file 1 (repId, entryId, cluFlag, taxId) | .../v6/1-AFDBClusters-repId_entryId_cluFlag_taxId.tsv.gz | v6 (AFDB v6 / UniProt 2025_03; updated 2026-08-10) | data/raw/afdb_clusters/v6/ | 93e89ac938b68cd37db6590307b1452d | CC-BY-4.0 | 199.4M rows; cluFlag 1 = Foldseek-clustered, 2 = AFDB50 member |
| AFDB Clusters v6 — file 2: cluster summary | .../v6/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz | v6 | data/raw/afdb_clusters/v6/ | fc42913f44d7ce67da2d0d43ea9ca95e | CC-BY-4.0 | 2,597,930 clusters; 322k have LCAtaxId = 0 |
| AFDB Clusters v6 — file 5: all 241M entries | .../v6/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz | v6 | data/raw/afdb_clusters/v6/ | f263bed7394777ab7d0033b94737431c | CC-BY-4.0 | **82.4M entries have taxId = 0** and 21,936 entries are listed twice |
| AFDB Clusters v6 — file 6: all-vs-all E-values | .../v6/6-all-vs-all-similarity-queryId_targetId_eValue.tsv.gz | v6 | data/raw/afdb_clusters/v6/ | c32cb51bfa29f3cabf11536d073a069a | CC-BY-4.0 | 552,702,072 rows |
| AFDB Clusters v6 — file 7: AFDB50 | .../v6/7-AFDB50-repId_memId.tsv.gz | v6 | data/raw/afdb_clusters/v6/ | 96a76ac9848ef717ef415f5844453d3b | CC-BY-4.0 | IDs given as AF-…-F1-model_v6 |
| S3 bucket listing snapshot | https://steineggerlab.s3.amazonaws.com/?prefix=afdb-cluster/ | 2026-10-30 | data/raw/afdb_clusters/s3_listing.xml | – | – | Provenance (ETags) |
| NCBI Taxonomy taxdump | https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | 2026-10-30 (Last-Modified 07:29 GMT) | data/raw/ncbi_taxonomy/taxdump/ | see taxdump.tar.gz.md5 (verified OK) | Public domain (NCBI) | 3,024,644 nodes. Ranks now use "domain" and "realm" in place of "superkingdom" |
| NCBI Taxonomy new_taxdump (rankedlineage, taxidlineage, fullnamelineage) | https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz | 2026-10-30 | data/raw/ncbi_taxonomy/new_taxdump/ | see new_taxdump.tar.gz.md5 (verified OK) | Public domain | Used for domain, kingdom and phylum lookups |
| TimeTree of Life 2015 species timetree (Hedges et al. 2015 MBE) | https://timetree.org/public/data/TimetreeOfLife2015.nwk | 2015 | data/raw/timetree/ | 2b746660008c8f185533b24a4b008b77 | TimeTree terms: personal research and teaching use only; **no redistribution** | 50,641 tips (49,751 eukaryotes, 152 bacteria, 25 archaea) |
| TimeTree 6 clade trees (Mammals, Tetrapods, Amphibians, Reptiles, Gymnosperms) | https://timetree.org/downloads/biodiversity/TT6_*.nwk | TimeTree 6 (Aug 2026) | data/raw/timetree/ | see report | TimeTree terms (no redistribution) | Vertebrate and plant clades only |
| TimeTree pairwise API (not downloaded) | https://timetree.org/api/pairwise/{taxidA}/{taxidB} ; https://timetree.org/api/taxon/{taxid} | TimeTree 6 (version 20260820) | – | – | TimeTree terms | Accepts NCBI taxIds and returns CSV. About 17 queries/s with 8 parallel requests. Returns nothing for very deep splits (e.g. human vs E. coli) |
| Foldseek DB afdb_proteome (AA + 3Di + Cα + headers + taxonomy for 48 model-organism and global-health proteomes) | https://foldseek.steineggerlab.workers.dev/afdb_proteome.tar.gz | v6 2025_03 | data/raw/foldseek_db/afdb_proteome/ | see data/raw/foldseek_db/md5sums.txt | CC-BY-4.0 | 590,435 structures |
| Pfam-A.regions (per-UniProt Pfam domain coordinates) + clans | https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/ | Pfam 38.2 (2026-01; UniProtKB 2025_03) | data/raw/pfam/ | regions f22fb1efefd1e52d6d9c64a8f10003ba; clans a70eeffcd8105ff094ea5834456d8bab; pfamA.txt 6c62b48afa7228e93958aa0532a271ab (all match the official md5_checksums) | CC0 | 4.7 GiB; download complete |
| TED consensus domain summary (CATH labels, pLDDT, taxonomy) | https://zenodo.org/records/13908086/files/ted_365m.domain_summary.cath.globularity.taxid.tsv.gz | TED v1 (Zenodo 2024-10-31; based on AFDB v4) | data/raw/ted/ | Zenodo md5 51468b8ec6c74388085ebe0f66f973c5 | CC-BY-4.0 | 19.9 GB; downloaded, local md5 matches Zenodo |
| TED per-chain segmentation | https://zenodo.org/records/13908086/files/ted_214m_per_chain_segmentation.tsv.gz | TED v1 | data/raw/ted/ | Zenodo md5 99b6afad9d670aa51813397d58f8e7c6 | CC-BY-4.0 | 9.6 GB; downloaded, local md5 matches Zenodo |

## Derived data (reproducible from the raw files with scripts)
| item | script | path |
|---|---|---|
| Per-cluster taxonomic span (number of domains, kingdoms and phyla, model-organism counts) joined with the file-2 statistics | scripts/afdb_cluster_stats.py v3/v6 | results/data_scout/afdb_clusters_{v3,v6}_taxspan.parquet ; afdb_cluster_stats_{v3,v6}.json |
| Overlap between v3 and v6 | scripts/afdb_v3_v6_overlap.py | results/data_scout/afdb_v3_v6_overlap.json |
| Mapping of TimeTree 2015 tips to NCBI taxIds, with coverage | scripts/timetree_coverage.py | data/interim/timetree/ttol2015_tips_taxid.tsv ; results/data_scout/timetree_coverage.json |
