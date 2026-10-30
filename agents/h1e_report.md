# H1e experimental-structure control set (data engineering report)

Date 2026-10-30. Scope: DESIGN.md §5 "H1e" and §11.10. Nothing was committed to git.

## 1. Raw data (data/raw/pdb_meta/, md5s in `md5sums.txt`, rows added to DATA.md)
| file | source | note |
|---|---|---|
| `pdb_chain_uniprot.tsv.gz` | SIFTS flatfile (2026-10-25; PDB 43.26, UniProt 2026.03) | 1,049,232 segment rows, 246,567 entries. Columns: PDB, CHAIN (author id), SP_PRIMARY, RES_BEG/END, PDB_BEG/END, SP_BEG/END |
| `pdb_entry_meta.tsv.gz` | RCSB GraphQL (bulk, 6 threads, 43 s) via `scripts/h1e/fetch_pdb_meta.py` | method, resolution, **initial release date**, deposition date |
| `entries.idx`, `resolu.idx` | wwPDB | cross-check only. `entries.idx` has the *accession (deposition)* date, not the release date, so it cannot be used for the AF2 cutoff; this is why RCSB GraphQL was chosen. Resolutions match `resolu.idx` for all X-ray entries (0 mismatches of 197,710) |
| `pdb_chain_modeled.tsv.gz` | RCSB GraphQL (`fetch_pdb_modeled.py`) | modeled-residue count per chain (669,083 chains) |
| `uniprot_3d_<taxid>.tsv`, `uniprot_3d_genes_<taxid>.tsv` | UniProt REST | length, taxon, gene symbol for accessions with `structure_3d:true` in the 11 taxa |

Licenses: PDB CC0, SIFTS CC-BY 4.0, UniProt CC-BY 4.0.

**Important caveat found and fixed.** SIFTS segments span SEQRES, including residues with no coordinates. Coverage computed from SP_BEG/SP_END alone is inflated (example: Q92917, 7dvq chain y, maps 476 residues but only 13 have coordinates). Coverage is therefore `min(SIFTS-mapped residues, modeled_residue_count) / UniProt length`. The uncorrected value is kept as `coverage_sifts_construct`. Validated below.

## 2. Table: data/processed/h1e/pdb_chains.tsv (built by `scripts/h1e/build_pdb_chains.py`)
- Filters: X-ray, or EM with resolution <= 3.0 A (flag in `method`); resolution <= 3.0 A; canonical UniProt accessions only (isoform accessions dropped).
- One chain per accession, chosen by highest observed coverage, then best resolution, then newest release. Only chains with coverage >= 0.70 are eligible where the UniProt length is known.
- Columns: uniprot, taxid, organism_id, reviewed, pdb_id, chain (author id), method, resolution, release_date, `post_af2_cutoff` (release > 2018-04-30), sp_beg, sp_end, n_sifts, modeled_residue_count, n_mapped, n_segments, uniprot_length, coverage, coverage_sifts_construct, qualifies.
- UniProt length is known only for the 11 taxa (UniProt REST) and the 48 afdb_proteome proteomes (Foldseek index, length = index col 3 - 2). Rows outside these keep coverage and `qualifies` empty (40,834 of 53,342 rows; ranked by n_mapped). 12,508 rows have known length and qualify (8,828 X-ray, 3,680 EM; 50% are post-cutoff).
- **Second table `pdb_chains_post_af2.tsv`**: best chain among entries released after 2018-04-30 only (DESIGN §11.10 requires this, and the best overall chain is often an older entry). 7,317 qualifying accessions.
- Not filtered: engineered mutants, chimeras, ligand-bound or complex chains (SIFTS does not flag them). Identity against UniProt is reported at extraction time (section 3).

## 3. Script: scripts/h1e/fetch_pdb_chains.py
`--acc A B ... | --acc-file list --table pdb_chains.tsv --out DIR --workers N`. Downloads mmCIF from files.rcsb.org and the per-residue SIFTS XML (ftp.ebi.ac.uk/.../sifts/xml/). It keeps the chain's polymer residues that SIFTS maps to the accession, renumbers them to UniProt numbering (residue-level via label_seq_id, so insertion codes and engineered indels are handled; the segment-level table has 2% of segments where RES and SP lengths differ, so the XML route was chosen over offsets), keeps the first altloc, drops H, water, ligands and unmapped tags, renames MSE to MET, and writes chain `A`. Output is `<acc>_<pdb>_<chain>.pdb` (`.cif` only if a UniProt number > 9999; none occurred) plus `fetch_summary.tsv` (identity against the UniProt FASTA at the same positions, extracted coverage).
- **Test, 20 accessions** (P69905, P68871 plus 18 random human/mouse/yeast/E. coli; `data/interim/h1e/chains_test/`): all 20 ok; 17/20 identity 1.000, others 0.993 (P69905, 1bab), 0.987 (P80075), 0.888 (P18531, engineered). `foldseek createdb` read all 20 without errors.
- P00533 (EGFR) has **no qualifying chain** (best X-ray kinase construct 1m17 covers 26% of 1210 aa). It was tested directly: 312 residues, identity 1.000.
- **Validation, 300 random qualifying chains** (human, mouse, yeast, E. coli, A. thaliana; `chains_val300/`): 300/300 extracted; identity mean 0.996, median 1.000, 91.7% >= 0.99, 97.7% >= 0.95, minimum 0.77 (engineered/chimeric constructs). Observed extracted coverage: median 0.91, 99% >= 0.70 (3/300 between 0.60 and 0.70); mean bias of the table coverage vs exact extracted coverage -0.002. Run time 100 s for 300 with 8 workers.

## 4. Counts of accessions with a qualifying chain (coverage >= 0.70, res <= 3.0 A)
`taxon_counts.tsv`; n_any = X-ray + EM, "post" = best chain among entries released after 2018-04-30.
| taxon | n_any | n_X-ray | post-AF2 any | post-AF2 X-ray | fraction of AFDB proteome entries (any / post) |
|---|---|---|---|---|---|
| human 9606 | 3502 | 2317 | 2524 | 1241 | 8.3% / 6.0% |
| mouse 10090 | 881 | 638 | 475 | 215 | 3.0% / 1.6% |
| rat 10116 | 298 | 200 | 143 | 43 | 1.3% / 0.6% |
| yeast 559292 | 1034 | 563 | 649 | 133 | 17.0% / 10.7% |
| A. thaliana 3702 | 616 | 411 | 385 | 177 | 2.1% / 1.3% |
| D. melanogaster 7227 | 210 | 130 | 136 | 54 | 1.4% / 0.9% |
| C. elegans 6239 | 205 | 103 | 132 | 29 | 1.0% / 0.6% |
| zebrafish 7955 | 132 | 59 | 90 | 17 | 0.5% / 0.3% |
| E. coli K-12 83333 | 1426 | 1206 | 616 | 357 | 32.5% / 14.1% |
| chicken 9031 | 173 | 112 | 96 | 35 | n/a (not in afdb_proteome; 0.4% of the 43.7k refprot count) |
| bovine 9913 | 299 | 188 | 177 | 57 | n/a (not in afdb_proteome) |
Denominators are AFDB proteome entries (Swiss-Prot + TrEMBL, so human is 42k, not 20k); yeast is S288C only.

## 5. Forecast of ortholog pairs with PDB structures on both sides
Gene-symbol-matched (upper-case) pairs of qualifying accessions, vertebrates (`pair_forecast.tsv`). These are lower bounds on true orthologs, since symbol matching misses renamed genes, and they include paralog-collapsed symbols.
| pair | both qualifying (any) | both from post-AF2 entries only |
|---|---|---|
| human-mouse | 513 | 304 |
| human-rat | 216 | 100 |
| human-bovine | 213 | 131 |
| human-chicken | 83 | 54 |
| human-zebrafish | 75 | 60 |
| mouse-rat | 97 | 68 |
| mouse-bovine | 97 | 71 |
| mouse-chicken | 52 | 42 |
| mouse-zebrafish | 67 | 58 |
| rat-bovine / rat-chicken / rat-zebrafish | 40 / 32 / 34 | 24 / 25 / 32 |
| bovine-chicken / bovine-zebrafish / chicken-zebrafish | 11 / 5 / 35 | 3 / 4 / 34 |

Rough order of magnitude for deep pairs (a guess, not computed, since no orthology table exists yet): human vs yeast, A. thaliana, D. melanogaster or C. elegans have of the order of 100-300 pairs (any), 50-150 post-cutoff; E. coli vs eukaryotes below about 100. Structured proteins are strongly enriched in both species (human-mouse observed 513 pairs against about 170 expected under independence, roughly 3x), so the true counts are above a naive product. Practical conclusion: after the post-2018 restriction, expect about 300-500 human-mouse-class pairs, on the order of 1000-2000 pairs in total across the panel for closely related taxa, and about 100-300 for each deep pair. This is enough for a paired AF-vs-PDB comparison, but the d_AA range will be dominated by short distances. A fair slope-vs-d_AA test needs the deep pairs (yeast, plant, fly, worm), where the counts are small.

## 6. Caveats and next steps
- UniProt length for non-panel accessions is missing; fill it from the FoldComp v4 models when the real ortholog list exists, then re-run `build_pdb_chains.py` (about 1 min).
- Coverage uses the RCSB modeled-residue count, an upper bound of observed mapped residues (it counts tags too). The validation shows 1% of chains fall slightly under 70%.
- Ranking is by coverage first, so EM and X-ray compete equally; filter `method == "X-ray"` for the strict control.
- The AF2 training cutoff is applied on the PDB initial release date. A post-2018 structure can still be a close homolog of a pre-2018 one; §11.10 allows that.
- Helper scripts: `fetch_pdb_meta.py`, `fetch_pdb_modeled.py`, `fetch_uniprot_lengths.py [genes]`, `build_pdb_chains.py`, `forecast_pairs.py` in `scripts/h1e/`. Logs in `logs/h1e_*.log`.
