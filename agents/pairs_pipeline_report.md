# Pairs pipeline: report (pipeline engineer, 2026-10-30)

**Bottom line.** `src/strucevo/pairs/` and `scripts/pairs/run_pairs.py` produce the ANALYSIS_PLAN §1.2 per-pair table for `seqonly`, `tmalign` and `tdiaa`. The output is sharded and resumable, with `.done` sentinels and `counts_*.npz`.

- **Verification.**
  - V5 (pair part), V8 and V9 pass on the euk20 demo (56,000 pairs).
  - The 13 pytest tests in `tests/test_pairs.py` pass.
  - ML distances recomputed from counts reproduce the table bit-for-bit.
- **Demo run.** All three alignments of 56k pairs took **4.5 min wall** on the shared node (load 30–70).
- **Projection for 3.5M pairs:** about **4.6 h** on 48 cores. TM-align alone accounts for about 3.0 h of that.

**Two findings the PI should see before freezing §Tuning:**
1. **TM-align alignments are noisy at the residue level for close homologs.**
   - On 5,295 human–mouse 1:1 orthologs, the median core AA p-distance is **0.091 under tmalign vs 0.064 under seqonly/tdiaa**.
   - Median global LDDT is **0.83 vs 0.93**.
   - For 25% of these pairs, the tmalign p-distance exceeds the seqonly one by more than 0.02.
   - 8 of 104 pairs of *identical* sequences (different accessions) get a core d_AA > 0 under tmalign, up to 3.2. They get 0 under the other two alignments.
   - `--tmalign-fast 0` changes nothing (medians identical) and is 10× slower.
2. **TM-align leaves distant residue pairs unaligned**, as standalone `TMalign` does. Coverage is therefore lower: the fail_cov rate is 10.7% for orthologs and **24.7% for paralogs** under tmalign, vs 3.9%/9.8% under tdiaa and 1.0%/1.7% under seqonly.
   - This interacts with the >10%-failed-pairs HOG exclusion rule (§11.6) and with E6.

## Usage

```bash
source /opt/conda/etc/profile.d/conda.sh && conda activate conda_envs/main
# 0) structures: FoldComp range-request fetch of the needed AFDB v4 entries (TSV acc, fc_key, entry)
python scripts/pairs/fetch_structures.py accessions.tsv data/interim/structures/<name>     # 70,760 entries in 25 s
# 1) pair table (builds the store once from the FoldComp DB, then runs seqonly -> tmalign -> tdiaa per shard)
python scripts/pairs/run_pairs.py --pairs pairs.tsv --foldcomp data/interim/structures/<name> \
    --store data/interim/structures/<name>_store --out data/interim/pairs --n-shards 200 \
    [--alns seqonly tmalign tdiaa] [--alpha-aa 1.0 --alpha-tdi 1.0] [--shards 3 7]
# 2) after §Tuning: recompute ML distances from counts (no realignment). Run seqonly first;
#    --s0 fills lddt_D/lddt_sat from floors.json
python scripts/pairs/recompute_ml.py --out data/interim/pairs --aln seqonly --alpha-aa A --alpha-tdi B [--s0 results/floors/floors.json]
python scripts/pairs/recompute_ml.py --out data/interim/pairs --aln tmalign --alpha-aa A --alpha-tdi B [--s0 ...]
# 3) verification
python scripts/verify/v8_v9_foldseek.py --store S --pairs pairs.tsv --n 1000 --json results/verify/V8_V9.json
python scripts/verify/v5_pairs_invariants.py --pairs-in pairs.tsv --out data/interim/pairs --store S \
    --fasta-glob '<proteome fastas>' --json results/verify/V5_pairs.json
python -m pytest tests/test_pairs.py -q
```

**Pair-list input.** A TSV with a header.
- Required columns: `gene_a`, `gene_b` (UniProt accessions).
- Optional metadata, passed through: `pair_set`, `hog`, `split`, `species_a`, `species_b`, `relation`, `one_to_one`, `dcg_id`, `event_level`, `dup_branch_from`.
- The pipeline reorients each pair so that gene_a is the shorter protein (ties broken lexicographically), swapping the species columns with it.
- Duplicate `pair_id`s raise an error.
- `scripts/pairs/make_demo_pairs.py` shows how to build the list from `hog_pair_labels.py` output.

## Implementation (what each piece does)

| Module | Content |
|---|---|
| `strucevo/pairs/store.py` | `foldseek createdb <FoldComp> --input-format 5`; then python `foldcomp.get_data` decodes CA atoms (atom 1 of each residue, using per-residue heavy-atom counts) and pLDDT (b-factors). **Validated against PDB-text parsing on 300 entries: 0 mismatches.** The Foldseek AA string must equal the FoldComp residues; the AA from `db` and the 3Di from `db_ss` are stored as ragged `.npy` arrays (memory-mapped, shared by forked workers). `meta.parquet` holds key, name, acc, length, offset, mean_plddt and md5(seq). 70,760 entries in 21 s. |
| `strucevo/pairs/align.py` | Fake prefilter via `tsv2db` (keys sorted by query). `tmalign -a 1 --tmscore-threshold 0 [--tmalign-fast 1]`; `structurealign -a 1 -e inf --tmscore-threshold 0 --lddt-threshold 0 --alignment-type 2`; `convertalis` with `query,target,qstart,qend,tstart,tend,qlen,tlen,alnlen,fident,cigar,lddt,alntmscore,qtmscore,ttmscore,rmsd,evalue,bits,lddtfull`. `seqonly`: Biopython `PairwiseAligner` (global, BLOSUM62, `open_gap_score=-11`, `extend_gap_score=-1`, `end_gap_score=0`; first optimal path; non-standard letters → X). |
| `strucevo/lddt.py` | In-house CA lDDT, Foldseek convention: query is the reference; R0 = 15 Å; neighbours are all query residues, and unaligned neighbours count as not preserved; thresholds 0.5/1/2/4 Å. Returns per-residue scores for the aligned residues. |
| `strucevo/pairs/posthoc.py` | Per pair: cigar → (qi, ti); core mask (pLDDT ≥ 70 in both); n_aln, n_core, cov_short; 20×20 AA and 3Di counts over core columns (both alphabets encoded with `AA_ORDER`, trap #4); for `e6_*` pairs, halves A/B = `(qi // 20) % 2`; for seqonly, AA counts over all aligned columns. lddt_global = Σ per-residue / len_a. lddt_core = mean per-residue over aligned core columns. rmsd_core = `strucevo.distances.kabsch_rmsd`. 48 forked workers, one BLAS thread each. |
| `strucevo/pairs/ml.py` | `ml_from_counts`: LG+Γ4 and Q.3Di.AF+Γ4 (md5 checked; loaded with `AA_ORDER`), d_max = 20, SE = max(SE, 1/n_core) with NaN → 1/n_core, α a parameter (placeholder 1.0). GPU, about 1.5 s per 56k pairs for all 6 count sets. `lddt_D(S, len_a, s0_table)` gives the floor-corrected LDDT distance plus lddt_sat. |
| `strucevo/pairs/pipeline.py` | Shard = md5(hog), or md5(pair_id) when hog = NA, modulo n_shards. Per shard: seqonly first (it supplies `aa_ml_seq`/`aa_pid_seq` to all tables), then tmalign and tdiaa. Writes `{aln}/shard_k.parquet`, `{aln}/counts_shard_k.npz` and `{aln}/shard_k.done` (JSON: rows in/out, parquet and npz md5, status counts, per-stage timings, α). A rerun skips any shard whose sentinel md5 matches (verified: 48/48 skipped in 4 s; a tampered sentinel triggers recomputation). |

### Row integrity
- Every input pair yields one row in every table.
- A pair with no structure in the store, or no Foldseek output row, gets `status = fail_align`.
- Foldseek rows for pairs that were not requested, or duplicated rows, raise an error.
- `run_pairs.py` asserts rows out = rows in, per shard and overall.
- Demo: 56,000 → 56,000 for each alignment, with 0 fail_align.

### Status rule
Statuses are assigned in priority order: `fail_align` > `fail_core` (n_core < 50) > `fail_cov` (cov_short < 0.5) > `ok`. ML values are still computed for failed pairs whenever counts exist; filter on `status`.

### Schema check
The columns, in §1.2 order, are:

`pair_id, pair_set, hog, split, gene_a, gene_b, species_a, species_b, relation, one_to_one, dcg_id, event_level, dup_branch_from, len_a, len_b, mean_plddt_a, mean_plddt_b, n_aln, n_core, cov_short, status, aa_p, aa_ml, aa_ml_se, aa_sat, tdi_p, tdi_ml, tdi_ml_se, tdi_sat, aa_ml_seq, aa_pid_seq, lddt_global, lddt_core, lddt_D, lddt_D_se, lddt_sat, tm_short, rmsd_core, aa_ml_hA(_se), aa_ml_hB(_se), tdi_ml_hA(_se), tdi_ml_hB(_se), boot_sub`

Choices made where §1.2 was silent. Please confirm or overrule:

| Item | Choice |
|---|---|
| `aa_ml_seq`, `aa_pid_seq` | From **all** aligned columns of the seqonly alignment, no pLDDT mask (selection variables that need no structure). pid = identities / aligned columns with two valid letters. |
| `tm_short` for seqonly | NaN: no TM-score is computed for a fixed sequence alignment. |
| `lddt_global` for tmalign/tdiaa | From Foldseek `lddtfull`, as the task specifies. §2.1 says to use in-house for tdiaa; the two agree to ≤ 5e-4 (V9), so the choice is immaterial. |
| `lddt_core` | All-neighbour per-residue lDDT averaged over aligned core columns. This is not the reviewer's variant (c), which also restricts the neighbours. |
| `lddt_D`, `lddt_D_se`, `lddt_sat` | Left NaN/null until floors exist; `recompute_ml.py --s0` fills lddt_D and lddt_sat. lddt_D is NaN only where it is undefined (S ≤ S0 median). Pairs with S ≤ q95 are flagged `lddt_sat` but keep their value (censored, never dropped). |
| `boot_sub` | False placeholder until §2.3 picks the subsample. |
| Seqonly gap scores | Biopython semantics: a gap of length k scores −11 − (k−1). The BLAST "11/1" convention would be `open_gap_score = −12`. **Please confirm which convention is intended.** |
| `--tmscore-threshold 0` | Means TM > 0, so a TM = 0 pair would be dropped. It would show up as fail_align; there were none in the demo. |
| Demo self-pairs | Use `pair_set = v5_self`, `relation = self`, which is not a contract value. Demo only. |

## Verification (demo, `data/interim/pairs_demo/verify/`)

| Check | Result |
|---|---|
| V5 rows | 56,000 in → 56,000 out per alignment; pair_id unique; gene_a shorter in 100% of pairs. **Pass.** |
| V5 self-pairs (1,000) | lddt_global = 1.000, tm_short = 1.000 (tmalign, tdiaa), aa_ml = tdi_ml = aa_ml_seq = 0, rmsd_core ≈ 1e-14. **Pass.** 35 self-pairs are `fail_core` (low-pLDDT proteins). |
| V5 identical sequences (104 non-self pairs, plus the self-pairs) | aa_ml_seq = 0 in all tables. **Pass.** Core aa_ml is 0 for seqonly and tdiaa but **> 0 for 8/104 under tmalign** (TM-align misregistration; reported as `identical_warn`). |
| V5 human–mouse 1:1 identity (5,295 pairs) | median aa_pid_seq **0.884** (target 0.80–0.90). **Pass.** Median aa_ml_seq is 0.137. Median core aa_ml is 0.069 (seqonly), 0.104 (tmalign), 0.069 (tdiaa). |
| V5 md5 identity | md5(UniProt 2026 FASTA) = md5(AFDB v4 sequence) for 70,304 of 70,760 proteins; **mismatch rate 0.64%**. The D1 filter (§1.6.1) must drop these; the pipeline only reports them. |
| V8 parsers (1,000 pairs × 2 Foldseek modes) | Column order and names; qlen/tlen = store lengths; qseq/tseq = store sequences; start/end bounds; lddtfull count = number of M columns; cigar-reconstructed qaln/taln identical to Foldseek output for **1000/1000** pairs. Convention confirmed: `I` consumes the query, `D` the target. **Pass.** |
| V9 lDDT (1,000 pairs × 2) | Global \|Δ\| median 3e-5, max **5.3e-4** (tmalign) and 5.0e-4 (tdiaa). The aligned-residue mean vs Foldseek `lddt` differs by at most 9.5e-4. The worst per-residue difference is 0.05–0.08 (a few residues; lddtfull has 3 decimals). **Pass** (criterion 0.01). |
| pytest `tests/test_pairs.py` | 13 pass, about 20 s. The tests cover: cigar parsing and string reconstruction; lDDT vs a brute-force reference and under rigid motion; seqonly free end gaps and exact gap scoring; post-hoc self-pair metrics and the half-block partition (hA + hB = core); Q.3Di.AF md5 and order; SE floor and empty counts; `lddt_D`; pair orientation, tie-break and duplicate detection; end-to-end on the demo store (schema order, fail_align for a missing structure, self-pairs, halves only for e6, resume and skip, and recomputation of a tampered shard). |
| recompute_ml | Recomputation at α = 1 from `counts_*.npz` reproduces every ML column exactly (max \|Δ\| = 0, identical NaN pattern). |

## Demo sanity numbers (euk20; 50k orthologs, 5k paralogs, 1k self-pairs)

| Alignment | Relation | ok | fail_core | fail_cov | Median cov_short | Median TM | Median LDDT_global | Median rmsd_core | Median aa_ml (ok) | Median tdi_ml (ok) | Σtdi_p/Σaa_p |
|---|---|---|---|---|---|---|---|---|---|---|---|
| seqonly | ortholog | 0.88 | 0.108 | 0.010 | 0.95 | — | 0.64 | 4.28 | 1.86 | 0.55 | 0.65 |
| tmalign | ortholog | 0.82 | 0.073 | 0.107 | 0.84 | 0.75 | 0.63 | 2.00 | 2.30 | 0.55 | 0.60 |
| tdiaa | ortholog | 0.88 | 0.082 | 0.039 | 0.92 | 0.69 | 0.66 | 3.02 | 1.98 | 0.42 | 0.54 |
| seqonly | paralog | 0.83 | 0.156 | 0.017 | 0.91 | — | 0.54 | 4.96 | 2.09 | 0.69 | 0.65 |
| tmalign | paralog | 0.67 | 0.086 | 0.247 | 0.70 | 0.58 | 0.47 | 2.16 | 2.65 | 0.69 | 0.62 |
| tdiaa | paralog | 0.79 | 0.114 | 0.098 | 0.82 | 0.49 | 0.55 | 3.29 | 2.15 | 0.48 | 0.53 |

- euk20 orthologs are mostly deep: plants vs animals vs fungi, with median pairwise identity 0.34. So most pairs lie beyond d_AA = 2. That is expected for a random draw and is not the production pair distribution.
- Saturation is rare: the `aa_sat` fraction among ok pairs is ≤ 3.1%.
- The 3Di/AA p-ratio again depends on the alignment (0.54–0.65; the reviewer's C4 found 0.56/0.69).
- Full tables are in `data/interim/pairs_demo/sanity.txt`.

## Throughput (56,000 pairs, 16 shards, 48 threads, shared node at load 30–70; upper bounds)

| Stage | Time | Pairs/s |
|---|---|---|
| FoldComp fetch (70,760 structures, range requests, 8 × 64 threads) | 25 s | 2,800 structures/s |
| Store build (createdb 8.8 s + decode 10 s) | 21 s | 3,400 structures/s |
| seqonly: Biopython align + in-house lDDT + counts + Kabsch | 28 s | 1,980 |
| seqonly, shard total including ML, I/O and npz | 40 s | **1,400** |
| tmalign: Foldseek `tmalign` | 144 s | 389 |
| tmalign: convertalis (with lddtfull) | 18.6 s | 3,000 |
| tmalign: post-hoc + ML | 3.4 s | 16,000 |
| tmalign, shard total | 174 s | **322** |
| tdiaa: structurealign | 17.8 s | 3,150 |
| tdiaa: convertalis | 20 s | 2,800 |
| tdiaa, shard total | 49 s | **1,140** |
| GPU ML (6 count sets per shard) | about 1.5 s per alignment | about 40k pairs/s including halves and overhead |

- **Fixed in this run:** the first run spent 3.5× longer in seqonly. The cause was 48 forked workers each running a multi-threaded BLAS. Workers now run one thread each (threadpoolctl, plus env vars in `run_pairs.py`); Foldseek still gets `--threads`.
- `--tmalign-fast 0`: 5,295 pairs took 269 s vs 26 s with fast mode, with identical medians. Keep fast mode.

**Projection for 3.5M pairs** (same length mix; HOG members have length within 0.5–2× of the median, so lengths are similar or shorter):

| Step | Time |
|---|---|
| Fetch, about 0.3M proteins | about 2 min |
| Store build | about 2 min |
| seqonly | 42 min |
| tmalign | **3.0 h** |
| tdiaa | 51 min |
| **Total** | **about 4.6 h** wall on 48 cores |

- Disk: about 0.4 KB per pair per alignment (parquet + npz), so about 4.5 GB in total. The store needs about 25 bytes per residue (arrays plus Foldseek DB): 857 MB for 34M residues in the demo.
- If time is short:
  - Run tdiaa on a stratified subsample (it is a sensitivity analysis), saving about 50 min.
  - Run shards of the three alignments on two nodes in parallel. Shards are independent; `--shards` selects a subset.
  - TM-align is the long pole and cannot be sped up with the GPU (review C11).

## Files

- **Code:**
  - `src/strucevo/lddt.py`
  - `src/strucevo/pairs/{__init__,store,align,posthoc,ml,pipeline}.py`
  - `scripts/pairs/{run_pairs,recompute_ml,fetch_structures,make_demo_pairs}.py`
  - `scripts/verify/{v5_pairs_invariants,v8_v9_foldseek}.py`
  - `tests/test_pairs.py`
- **Demo data** (`data/interim/pairs_demo/`, regenerable, not for git):
  - `pairs.tsv`, `accessions.tsv`
  - `fc` (FoldComp, 512 MB)
  - `store/` (857 MB)
  - `pairs/{seqonly,tmalign,tdiaa}/` (73 MB)
  - `verify/{V5_pairs,V8_V9}.json`
  - `sanity.txt`, `run.log`, `run_v0_oversubscribed.log`
