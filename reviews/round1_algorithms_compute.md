# Review round 1: algorithms and compute (reviewer: structural-bioinformatics / HPC)

Benchmarks: 2026-10-30, 10:33–10:45, shared node at load 40–85 (timings are upper bounds). Scripts are in `data/interim/review_bench/`: `extract.py`, `posthoc.py`, `lddt_var.py`, `ml/mldist.py` and `ml/test_ml.py`.

## 1. Verdict

**Feasible; accept with changes.** Pairwise metrics are much cheaper than §7 says (Foldseek realigns about 37k homolog pairs/s; GPU ML distances take seconds per million pairs). The long pole is Ancestra plus the structure fetch. Four gaps can bias H1a/H1c:

- **(i)** silent pair loss from the default E-value filter;
- **(ii)** local, query-referenced metrics (lDDT averages over aligned residues only);
- **(iii)** the 3Di/AA ratio depends on the alignment method (measured);
- **(iv)** no pLDDT masking inside Foldseek (cheap post hoc; validated).

Q.3Di.AF is available and verified (CC0).

## 2. Change proposals

**C1 (must). Pair-list alignment recipe and integrity checks.** Use the recipe in §3.

- Always pass `-e inf`. Without it, only 673 of 10,000 random pairs survived the default `-e 10` filter.
- Assert that the number of rows out equals the number of pairs in.
- Drop `--exact-tmscore`: it changed nothing in convertalis output (r = 1.000). `--mask-bfactor-threshold` only masks seeding.

**C2 (must). Build the structure DB from FoldComp v4, never from the prebuilt Foldseek AFDB DBs.**

- The prebuilt DBs are v6 without pLDDT (breaks D1). `createdb --input-format 5` reads a `createsubdb` FoldComp subset directly, keeping `AF-…-model_v4` names (2,000 entries in 2.5 s).
- FoldComp keeps per-residue pLDDT (0.01 resolution). Extract it once, with CA coordinates, into a ragged parquet/npz.

**C3 (must). Define the masked and global metric variants exactly, and compute them post hoc from one parsed alignment.**

- convertalis emits `qstart,tstart,cigar` plus scalars (verified: `I` consumes the query).
- **Foldseek lDDT, reproduced exactly** (median |Δ| = 0.0000, r = 1.000): CA-only, query as reference, R0 = 15 Å, neighbours from *all* query residues (unaligned count as not preserved), averaged *over aligned query residues only*. So it is asymmetric and inflated for partial alignments.
- **Variants to report:**
  - **(a) lDDT_aln:** Foldseek's value.
  - **(b) lDDT_global:** averaged over all query residues, with unaligned residues scoring 0. Compute it in both directions and take the mean. Use this as the primary variant for "global divergence".
  - **(c) lDDT_core:** neighbours and averaging restricted to residues with pLDDT ≥ 70 in both proteins.
- **RMSD.** Foldseek `rmsd` = Kabsch over all aligned CA (my Kabsch matched, median |Δ| = 0.000). Core RMSD = Kabsch over aligned pairs with pLDDT ≥ 70 in both; a CPU 3×3 SVD suffices, GPU unnecessary.
- AA/3Di masked p-distances and 20×20 count matrices come from the same index pairs (3Di strings from `db_ss`).
- **Cost:** 1 ms/pair on one core (L = 277) → about 4–8 CPU-h for 5M pairs at L ≈ 550, i.e. 10–15 min on 48 processes.

**C4 (must). Make H1c (the 3Di/AA rate ratio) robust to the alignment method.**

- On 50k human–zebrafish/Arabidopsis homolog pairs, the ratio (1 − 3Di identity)/(1 − AA identity) was:

  | Alignment | Ratio |
  |---|---|
  | Foldseek 3Di+AA | 0.56 |
  | Foldseek `tmalign` (coordinates only) | 0.69 |

- The 3Di+AA objective rewards 3Di matches, inflating 3Di conservation. Use **coordinate-only (TM-align) alignments as primary for H1c**, 3Di+AA as sensitivity. Take TM-score from TM-align too: it is a median 0.087 higher than the 3Di+AA TM-score (r = 0.95), with alignments a median 130 residues longer.

**C5 (must). Add a coverage and floor rule before any ML distance.**

- Local alignments of *random* pairs give fident 0.11–0.21 over 30–190 residues: finite "distances" without homology, i.e. artificial saturation.
- A priori: aligned core residues ≥ 50 *and* coverage ≥ 0.5 of the shorter protein; flag saturation when the ML optimum hits t_max = 20; take each metric's floor from about 100k random cross-HOG pairs run through the same pipeline.

**C6 (should). Use Q.3Di.AF as the 3Di model, with a per-alphabet Γ.**

- **Source:** Edmond doi:10.17617/3.1MJJBH, `Q.3Di.AF` (md5 8f897fffd14f247afbee44dad2fce2ce), `Q.3Di.LLM` (md5 3a6c99bd…); CC0; PAML format, ARNDCQEGHILKMFPSTWYV order over 3Di letters. Downloaded to `data/interim/review_bench/`; move to `data/raw/` and DATA.md.
- Not built into IQ-TREE 3.1.4, but `-m Q.3Di.AF+G4` with the file works. Foldseek `3di.out` is log-odds, not a rate matrix; no own GTR-20 unless Q.3Di.AF fails.
- **Γ:** one α per alphabet from IQ-TREE on FoldMason MSAs of the H1f HOGs, fixed topology. k depends on α, so report k for α ∈ {0.5, 1, 2, ∞}.
- **Cross-check k** with tree-length ratios (Q.3Di.AF+G4 vs LG+G4, same MSA and topology) on the H1f HOGs.

**C7 (should). Pairwise ML on GPU: yes, as prototyped.**

- Eigendecompose Q once; log(π_i·mean_k P_ij(t·r_k)) on a 512-point log grid; log-likelihood = counts (B×400) × table; three per-pair refinements.
- On 20k simulated pairs (300 sites, t = 0.02–3) the median relative bias was < 0.5% (LG and Q.3Di.AF). On an AliSim 4-taxon tree (2,000 sites) it recovered patristic distances within 2–6%. 1M pairs take 4.1 s on the A100.
- Do not call IQ-TREE per pair. Fix: the 1e-4 grid floor must allow 0.

**C8 (must). Change the order of node-1 work: put the structures and Ancestra on the critical path now, and pilot on euk20.**

- **Panel coverage rule.** D1 requires v4 models. Proteomes from 2026_03 contain accessions added after 2021_04 that have no v4 model. Before fixing the panel, compute the v4 coverage of each proteome against the local `afdb_uniprot_v4.lookup`, and drop species below 80% coverage.
- **Decouple Ancestra from the fetch.** Run Ancestra on the 2026_03 canonical FASTAs right away, with no wait for the fetch. Enforce D1 at the pair level: require the md5 of the UniProt sequence to equal the md5 of the FoldComp sequence, and report the mismatch rate.
- **Ancestra runtime.** Scaling from euk20 (357k proteins, 6.1 CPU-h) with the authors' exponent of 1.35 gives about 130 CPU-h (`accurate`) for about 3.4M proteins. That is about 3–5 h wall at 48 threads, matching the design's 4–8 h. However, S2 is not resumable *within* the stage.
- **Ancestra launch rule.** Launch by about 13:00. Cap the panel at about 3.5M proteins. Use `fast` if the projected wall time is over 6 h.
- **Pilot.** The existing `results/ancestra/runs/euk20_accurate` (12M labelled pairs) is a ready end-to-end pilot for the whole metric pipeline (C1–C7 and the toy tests) on node 1.

**C9 (should). Checkpointing.**

- About 200 shards by hash(HOG), each `tsv2db → structurealign → convertalis → post-hoc → parquet`, with a `.done` sentinel (row counts, md5); reruns skip done shards. Do not keep `qaln/taln/lddtfull` text (2.3 GB per 852k pairs).

**C10 (could). Run TM-align for all primary pairs; keep US-align for validation only.**

- Foldseek `tmalign` ran 50k long pairs in 190 s (about 265/s under load), so about 3M P-orth pairs take about 3 h on node 2. Validate against `USalign` on 1k pairs.

**C11 (could). No GPU for Foldseek.** The GPU only speeds up the prefilter (5k × 56k search: 50 s GPU vs 49 s CPU), which given-pair alignment does not use.

## 3. Recommended pipeline (measured)

```bash
# DB (once): FoldComp v4 subset -> Foldseek DB (+ pLDDT/CA ragged arrays via python foldcomp)
foldseek createsubdb keys.txt data/raw/foldcomp/<local_v4_subset> fc_sub --subdb-mode 0
foldseek createdb fc_sub structdb --input-format 5 --threads 48          # 2,000 entries: 2.5 s (8 thr)
# per shard: pair list = DB keys "qkey\ttkey\t0\t0", sorted by qkey
foldseek tsv2db pairs.tsv pref --output-dbtype 7                          # 1M pairs: 0.15 s
foldseek structurealign structdb structdb pref aln -a 1 -e inf --alignment-type 2 --threads 48
foldseek tmalign        structdb structdb pref tma -a 1 --threads 48      # primary for C4/C10
foldseek convertalis structdb structdb aln aln.tsv --format-output \
  query,target,qstart,qend,tstart,tend,qlen,tlen,alnlen,fident,cigar,lddt,alntmscore,qtmscore,ttmscore,rmsd,evalue,bits
# 3Di strings: read structdb_ss directly (or convertalis on *_ss DBs with symlinked _h for qaln/taln)
python posthoc.py  # cigar -> index pairs; masks; lDDT a/b/c; core Kabsch; AA & 3Di 20x20 counts
python mldist.py   # GPU, LG+G4(alpha_AA), Q.3Di.AF+G4(alpha_3Di)
```

| Step | Measured | Projection for 5–10M pairs |
|---|---|---|
| structurealign, random pairs (mean length 277) | 1M in 3.4 s | < 1 min |
| structurealign, homolog pairs (mean length 621) | 852k in 22.7 s | 2–5 min |
| convertalis with all string fields (852k pairs) | 83 s, 2.3 GB | Use the cigar-only output |
| Foldseek tmalign (50k long pairs) | 190 s | 3M pairs in about 3 h |
| Post-hoc lDDT/RMSD (NumPy, 1 core) | 1 ms per pair (L = 277) | 10–15 min on 48 processes |
| GPU ML distance | 1M in 4.1 s | < 2 min for both alphabets |
| FoldComp fetch (data scout) | 170 structures/s per process | 3.5M in about 1–1.5 h with 8 processes |
| Ancestra on about 3.4M proteins | scaled from euk20 | 3–5 h wall (`accurate`) |

**Corrected §7.** On node 1: the panel with the v4-coverage gate, the FoldComp fetch and DB build, the Ancestra launch, and the euk20 pilot of the metric pipeline. On node 2: everything else. The pairwise metrics need only about 4–5 h in total, including TM-align, so H1f trees and the H2e case studies have room.
