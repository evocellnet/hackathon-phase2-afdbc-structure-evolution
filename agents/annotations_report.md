# Per-protein annotation of the panel (AFDB v3 clusters, TED, Pfam): report

Date: 2026-10-30. Scripts: `scripts/annot/{clu_extract,afdb_cluster_annot,ted_extract,pfam_extract,build_protein_annot}.py`. Logs: `logs/annot_*.log`.
Output: `data/processed/annotations/protein_annot.parquet` (3,738,038 rows = all eukaryote + prokaryote panel proteins, 62 columns, one row per `acc`; not committed, 327 MB).
Summary (committable): `results/annotations/summary.json` (blocks: all, euk, prok, per superkingdom, and `*_seq_match_only`). Intermediates: `data/interim/annotations/`.
`clu_extract.py` was run once with the scratchpad as output directory (it writes `f5_panel.parquet`); `afdb_cluster_annot.py` takes that directory as argument. The other scripts take no arguments except `ted_extract.py domains|chains`.

## 1. AFDB v3 cluster columns (`afdb_*`)

Column semantics were checked on the data, not assumed (the README is not shipped; the scout report was the starting point).

- **File 5 (all members) is the authority.** `repId` is the representative of the entry's *Foldseek (structural) cluster*. For cluFlag 1 (AFDB50 non-representative member) it is the structural-cluster rep of the member's AFDB50 representative. The AFDB50 map (file 7) is not consistent with file 5 for about 8% of cluFlag-1 entries (rep IDs substituted in Dec 2024), so file 7 is not used.
- **cluFlag (file 5, v3):** 1 = AFDB50 member (non-representative), 2 = AFDB50 representative that is a member of a non-singleton Foldseek cluster (always present in file 1), 3 = fragment, 4 = singleton (`repId` = itself, never in file 2). This matches the scout report. File 3 (human only) uses a different, partly different flag scheme (not used).
- **Non-singleton cluster (`afdb_in_cluster`)** = `repId` is a row of file 2 (2,302,908 clusters, all nMem >= 2). Verified: every flag-2 entry, 1,940,203 of 2,129,291 flag-1 entries (the other 189,088 belong to an AFDB50 group whose representative is a Foldseek singleton) and 55,522 of 119,019 fragments are in such a cluster; flag 4 never.
- `afdb_cluster_nmem` = file 2 nMem (cluster members = AFDB50 representatives, fragments excluded; equals the file-1 row count for every panel cluster). `afdb_cluster_nentries` = all AFDB entries assigned to the rep in file 5 (adds AFDB50 members and fragments). `afdb_dark` = file-2 isDark. Also `afdb_rep_len/avg_len/rep_plddt/avg_plddt`.
- **LCA:** `afdb_lca_taxid` (file 2) with scientific name, NCBI rank (taxdump 2026-10-30) and `afdb_lca_group` (Eukaryota / Bacteria / Archaea / cellular_organisms = LCA at the root of cellular life / Viruses / root_other). NCBI names the 131567 rank "cellular root"; it is the LCA of 52% of clustered eukaryote panel proteins' clusters... see summary (the LCA of clusters spanning several domains of life).
- Proteins without an AFDB v4 model (172,172) are absent from file 5: all `afdb_*` null, `afdb_in_cluster` = false. Use `afdb_cluflag is not null` to restrict to entries (the summary reports both denominators).

## 2. TED columns (`ted_*`)

- Domain summary columns were verified on the data (21 columns; col 3 consensus high/medium, col 4 chopping, col 5 nres, col 7 domain pLDDT, col 14 CATH label, col 15 level). Only `high` and `medium` consensus domains exist in this file; low-consensus counts per chain are in `ted_n_low_consensus_chain` (not used for architectures). TED is computed on the AFDB v4 model, so use it only for `seq_match` proteins (99.79% of modelled proteins; `seq_match_ted` == `seq_match`).
- Checks: n_high + n_medium of the per-chain file equals `ted_n_domains` for all 3,565,908 modelled proteins; chain nres equals `len_v4` for all.
- Fields: `ted_n_domains` (null if no model, 0 if none), `ted_n_high/medium`, `ted_n_cath` (labelled domains), `ted_n_cath_H` (labelled to full C.A.T.H), `ted_boundaries` ("31-207;703-793", segments of discontinuous domains joined by `_`), `ted_arch_cath` (labels N->C by domain start, unlabelled domains omitted, e.g. `3.80.20.20|3.30.200.20|1.10.510.10`), `ted_arch_all` (same with `NA` placeholders), `ted_consensus`, `ted_dom_residues`, `ted_frac_in_domain` (= domain residues / `len_v4`, 0 if modelled but no domain), `ted_mean_dom_plddt` (residue-weighted), `ted_min_dom_plddt`.
- **CATH label depth:** labels are given at the deepest level available per domain: 4 levels (C.A.T.H) for 'H' assignments, 3 levels (C.A.T) for 'T' assignments; `-` = unassigned. 69%/80% (euk/prok) of modelled proteins have >= 1 labelled domain.
- **Class rule (ANALYSIS_PLAN 4.3):** class digit 1 / 2 / 3 -> `alpha` / `beta` / `alphabeta` if it covers >= 60% of the TED-domain residues of the protein (denominator = all TED domains, so unlabelled domains count against), else `mixed`; `none` = modelled with no TED domain; `no_model` = no AFDB v4 model. Class 4 (few secondary structures) and rare other digits count in the denominator and can only produce `mixed`.
  - **Caveat for E3:** under this literal rule `mixed` contains many proteins with no CATH-labelled domain at all (euk 20.3% of modelled proteins; prok 16.7%), plus genuinely multi-class ones (euk 5.5%). `ted_class_labeled` is the sensitivity version (denominator = labelled residues; `unlabeled` if none). Also given: `ted_frac_class{1..4}`, `ted_frac_unlabeled`, `ted_top_class_digit`, `ted_top_class_frac`, `ted_class_digits_present`.
  - At HOG level the plan's >= 60% rule should be applied to the HOG's pooled TED residues; protein-level columns here allow that aggregation (`ted_frac_class*` x `ted_dom_residues`).

## 3. Pfam columns (`pfam_*`)

- Pfam 38.2 regions (UniProtKB 2025_03), matched by accession. `pfam_seq_match_2026` = the Pfam sequence md5 equals md5 of the 2026 sequence (99.8% of proteins with hits), so envelope coordinates apply to the 2026 sequence; `pfam_has_hit`, `pfam_n_regions`, `pfam_n_families`, `pfam_arch` (accessions N->C by envelope start, all Pfam-A types, nested/overlapping regions kept), `pfam_arch_domain_repeat` (only types Domain and Repeat), `pfam_boundaries` (envelope), `pfam_frac_covered` (union of envelopes / `len_2026`, only for matching md5).
- Limitation: Pfam-A.regions lists only hits, so a protein without hit cannot be distinguished from "not in Pfam sequence database"; since 99.8% of hit-bearing panel proteins have the identical sequence, the panel is well covered, but absence of a hit is not proof of a protein without domain.

## 4. Coverage (details in summary.json)

| | eukaryotes | prokaryotes |
|---|---|---|
| proteins / with v4 model | 3,460,494 / 95.1% | 277,544 / 99.3% |
| modelled with >= 1 TED domain | 89.2% | 96.4% |
| modelled with >= 1 CATH-labelled domain | 68.9% | 79.7% |
| mean TED domains per modelled protein / mean fraction in domain | 1.64 / 0.56 | 1.76 / 0.79 |
| class (of modelled): alpha / beta / alphabeta / mixed / none | 19.6 / 9.8 / 24.0 / 35.8 / 10.8 % | 13.9 / 8.4 / 44.0 / 30.1 / 3.6 % |
| of AFDB entries: in non-singleton cluster | 79.8% | 90.9% |
| cluFlag 1 / 2 / 3 / 4 (counts) | 1,937,805 / 815,195 / 118,356 / 418,986 | 191,486 / 65,751 / 663 / 17,666 |
| clusters: median nMem / dark | 56 / 5.2% | 317 / 4.8% |
| proteins with a Pfam hit | 61.9% | 67.2% |
| TED domain or Pfam hit (modelled) | 92.3% | 97.3% |

Spot checks: human EGFR P00533 -> TED `3.80.20.20|3.30.200.20|1.10.510.10`, Pfam `PF01030|PF00757|PF01030|PF14843|PF21314|PF07714`, cluster LCA Bilateria; E. coli HslV P0A7B8 -> `3.60.20.10` (Ntn hydrolase), class alphabeta, Pfam PF00227, cluster LCA cellular root.

## 5. Notes and decisions

- TED and cluster annotation exist only for proteins with an AFDB v4 model; the structural analyses should use `seq_match` proteins (md5 equality), as per D1.
- The TED consensus levels of the domain file are `high`/`medium` only; no filtering by pLDDT was applied (domain pLDDT is provided).
- Per-domain `cath` labels at T level (3 levels) are kept as labels of the form `C.A.T`; architectures therefore mix 3- and 4-level labels.
