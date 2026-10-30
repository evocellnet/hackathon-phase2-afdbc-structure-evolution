# Panel tree and HOG construction: report (pipeline engineer, 2026-10-30, project clock 3.5 → 3.9 h)

**Bottom line.** The panel tree is built and dated, `scripts/hogs/build_hogs.py` implements ANALYSIS_PLAN §1.6 end to end, and the euk20 test passes. The script can run on euk180 as soon as `data/interim/ancestra_euk180/ancestra/hogs.orthoxml` appears.

- On the euk20 run: 11/11 pytest pass, lazy and exhaustive selection give byte-identical outputs, and the pair lists run through `run_pairs.py` with matching `aa_ml_seq`.
- **One PI decision is needed before the first real run: E6 caps** (§4). Without them, E6 probably exceeds the 3.5M pair budget on its own, and the plan's budget rule would then shrink the HOG sample.

## Commands

```bash
source /opt/conda/etc/profile.d/conda.sh && conda activate conda_envs/main
# done once (already run): per-protein fragment/cluFlag table, panel tree
python scripts/hogs/gene_flags.py          # -> data/interim/hogs/gene_flags.parquet (~5 min)
python scripts/hogs/panel_tree.py          # -> data/interim/timetree/{panel_tree.nwk,node_ages.tsv,...} (10 s)
# euk180 (after Ancestra finishes; 48 threads)
X=data/interim/ancestra_euk180/ancestra/hogs.orthoxml
python scripts/hogs/build_hogs.py --orthoxml $X --threads 48 --projection-only   # optional dry run: counts only, no split / pair lists
python scripts/hogs/build_hogs.py --orthoxml $X --threads 48 [--e6-max-o 5 --e6-max-per-dcg 5]
python -m pytest tests/test_hogs.py -q
```

- The selection distances are cached in `data/interim/hogs/seqdist_cache_alpha1.parquet`, so a real run after a projection run does not recompute them.
- `split.tsv` is written by the first non-projection run. **Fix the E6 flags before that run.** Any later change to the sample makes the script stop with "exists and differs — never regenerate".

## 1. Panel tree (`scripts/hogs/panel_tree.py`; outputs in `data/interim/timetree/`)

**Topology** comes from the TimeTree-of-Life 2015 newick, induced on one tip per panel species.

- **Direct match:** 91 of 180 species are tips themselves.
- **Proxies:** for the others, the deepest NCBI taxon X that contains 2015 tips is used. The mapping table is `panel_tree_mapping.tsv`.

  | Case | Rule | Species |
  |---|---|---|
  | `proxy_sole` | No other panel species in X, so any tip of X places the species exactly | 43 |
  | `proxy_group` | A "neutral" tip of X (one that sits in no sub-taxon holding another panel species) | 16 |
  | `shared_tip` | Congeners on one tip | 14 |
  | `ncbi_attach` | Polytomy child of the panel MRCA of the deepest NCBI taxon with another panel species. Used when no neutral tip exists or when the species is undated against all members of X | 16 |

  Examples of `ncbi_attach`: *Fonticula* → Opisthokonta, *Tritrichomonas* → Parabasalia, *Brachionus* → Spiralia, *Wallemia* → Basidiomycota.
- **Re-resolution:** within every proxy group, the 2015 placement is arbitrary. Each such group whose members form a clade and are all dated pairwise is therefore re-resolved by UPGMA on TimeTree 6 pairwise times (ties ≤ 0.01 My become polytomies). This was done for 19 groups, all small (2–6 species), e.g. Saccharomycotina, Spiralia and Sinocyclocheilus. 14 groups were skipped (undated pairs or not monophyletic); they keep the 2015/NCBI structure. Everything is logged in `panel_tree_fixes.tsv`.

**Ages.** For each internal node I took all spanning species pairs from the cached TimeTree 6 API ages (`pairwise_ages.tsv`).

- `median` = TimeTree `median_age` (the API median of study estimates) of the representative pair, which is the spanning pair at the lower median. lo95/hi95 = that pair's CI.
- Single-study pairs carry CI = (0, 0), which I treat as no CI. `adjusted_age`, n_span and min/max are kept as extra columns.
- Monotonicity was fixed top-down by setting a child to its parent's age: **34 fixes**, all deep (e.g. Ecdysozoa 853 > Protostomia 608; Streptophyta 1131 > Viridiplantae 974). TimeTree's deep study medians are noisy, e.g. *Entamoeba*–*Gossypium* is 928 Mya from a single study.
- **Collapse:** 166 → 100 internal nodes. 38 nodes were collapsed for internodes < 5 My (including all zero-length nodes left by the monotone fixes) and 28 for CI overlap with the parent.
  - The CI-overlap rule collapses much of the deep structure. Opisthokonta, Amorphea and several fungal nodes join the root polytomy at 1455 Mya.

**Undated nodes.** The 355 undated pairs are placed by NCBI taxonomy, and **6 nodes have age NA**. Each holds ≤ 4 species, so none can be the MRCA of an 8-species S_h.

| Node | Species |
|---|---|
| N003 | *Naegleria* + Parabasalia |
| N004 | *Naegleria* (2 species) |
| N005 | *T. vaginalis* + *Tritrichomonas* |
| N031 | Pucciniomycotina + Ustilaginomycotina |
| N085 | *Entamoeba* (2 species) |
| N094 | Ochrophyta subset |

- I found no literature age that is "clearly documented" enough to fill these, so none was added.
- *Fonticula* sits in the root polytomy, by NCBI placement under Opisthokonta, which itself collapsed into the root. Its pairs therefore get the root age as an upper bound, flagged `attach=ncbi_polytomy` in `node_ages.tsv`. Holomycota is not used because NCBI has no such node; PI may prefer a Holomycota node with an NA age.
- Several congeneric pairs are also undated in TimeTree 6 (*C. latens*, *T. rangeli*, *N. fowleri*, *E. dispar*, *G. tomentosum*). They are attached to their genus MRCA.

**Checks.**

- 180 leaves equal the panel.
- Ages are strictly monotone, and every internode is ≥ 5 My (test).
- **Human–mouse = N052 (Boreoeutheria) 89.5 Mya [86.0, 92.9]**, which passes 87 ± 10. The direct human–mouse API median is 83.5.

## 2. HOG construction (`scripts/hogs/build_hogs.py`)

The script follows §1.6 exactly. Implementation choices:

- **OrthoXML** is read with Ancestra's own reader (imported from `tools/ancestra-mvp-ai`). Top-level HOGs only. Genes outside any group are not candidates.
- **Gene filters, in this order:**
  1. md5 match (`protein_table.seq_match`; genes not in the panel table fail here).
  2. Fragment: UniProt `(Fragment)` in the 2026_03 FASTA description, or AFDB v3 file-5 cluFlag = 3.
     - Both sources are available. `scripts/hogs/gene_flags.py` gives 196,393 UniProt fragments and 118,356 cluFlag-3 entries among the 3.46M panel proteins. 170,152 panel accessions are absent from file 5, because they have no v4 model; they already fail the md5 filter.
  3. Mean pLDDT ≥ 70.
  4. Length within 0.5–2× the median of the HOG members that passed 1–3.
  5. 1:1 species = species with exactly one *filtered* copy; HOGs need ≥ 8.
- **Presample:** greedy panel-node coverage. The gain is the number of uncovered panel-tree ancestors; ties are broken by `SeedSequence(20261030, spawn_key=(stage, crc32(hog)))`, with stages from `config/seeds.yaml`.
- **aa_ml_seq:** `scripts/hogs/seqdist.py` calls the pairs pipeline's own seqonly functions (`seqonly_pairs` → all-column AA counts → `ml_from_counts`) on sequences only, so no structures are needed. It equals the pipeline's `aa_ml_seq` to within 2e-15 on 300 demo pairs (test) and on 704 pairs run through `run_pairs.py`.
- **Lazy selection (default):** pairs are computed on demand, deepest species pairs first. A clade is rejected at its first failing pair, and a failing pair is reused for every larger clade.
  - On euk20 this gives **identical hogs, split and pair lists (byte-identical) to `--exhaustive`** (all presample pairs, the plan's literal recipe), with 101k vs 128k pairs computed.
  - For euk180 the savings should be large. With 60-species presamples, exhaustive is ≤ 1,770 pairs/HOG × ~20k HOGs ≈ 35–45M pairs (~3 h on 48 cores). I estimate lazy at 5–10M (~30–45 min). Throughput is ~80 pairs/s/core.
- **Sampling:**
  - Equal quotas per stratum. Unused quota is redistributed to the other strata in proportion to their spare HOGs (largest-remainder rounding, seeded).
  - Each stratum takes a prefix of a seeded permutation. A budget rescale is therefore nested (it takes a shorter prefix).
  - The split is a within-stratum seeded permutation with round(0.2·n) tune. An existing `split.tsv` is checked and never rewritten; a test with a tampered file shows the script aborts.
- **Labels:** an on-demand deepest-common-group labeller (`HogIndex`), proven identical to `hog_pair_labels.label_pairs` on all pairs of 150 euk20 HOGs (> 10k pairs, test). A full O(n²) label table is never built for euk180.
- **E6 (§4.6):** uses the filtered membership of the sampled HOGs and the OrthoXML guide tree, because that tree defines P.
  - P1/P2 is an `in_paralog` pair with P not a leaf, `aa_ml_seq` ≤ 2 and `aa_pid_seq` ≤ 0.95.
  - O1 / O2 are the species in the sister clade(s) of P and of parent(P) with exactly one filtered gene in the HOG, that gene being an ortholog of both copies. |O1| ≥ 1 is required.
  - P1 is a seeded coin flip.
  - Pair sets: `e6_para` (P1,P2); `e6_orth` (each copy × O1); `e6_ldo` (each copy × O2). Pairs are de-duplicated, and pairs already in `primary` are dropped. `e6_triplets.tsv` holds the triplets.
- **`hogs.parquet`:** the §1.1 columns plus `selected_genes`, `presample_genes`, `sel_mrca`, `n_genes(_filtered)`, `n_primary_pairs` and `n_e6_pairs`.
  - `cath_class`, `misfit_ratio`, `excl_*` and `fail_frac_*` are placeholders (null) for the later steps.
  - `root_age_mya` is the panel median age of MRCA(S_h). The file lives in gitignored `data/interim/`, as the contract specifies.

## 3. Tests and dry run on euk20 (`data/interim/hogs_test_euk20/`)

17 of the 20 euk20 species are in the panel. Zebrafish, chicken and *Xenopus* are not, so their genes fail at the md5 step.

**`tests/test_hogs.py`: 11 passed (75 s).** It covers:

- panel tree: leaves, monotone ages and ≥ 5 My internodes, human–mouse age, newick round-trip
- seqdist vs pipeline
- HogIndex vs label_pairs
- contract columns and ranges
- clade-closure, with all S_h pairs recomputed ≤ 2
- primary count = ΣC(n,2) and unique pairs
- E6 invariants: no leaf P; P1/P2 in the same species; O1 ∩ O2 = ∅; relations; no overlap with primary
- split fractions and split.tsv consistency

**Pipeline compatibility.** `run_pairs.py --alns seqonly tdiaa` on 896 listed pairs (40 HOGs, demo store) gave 896 → 896 rows, all ok. The pipeline's `aa_ml_seq` equals the selection value (|Δ| ≤ 9e-16), and E6 split halves are filled.

**Filter counts** (`logs/hogs_genes.tsv`; the production run writes `logs/filters/hogs_*.tsv`):

| Step | Genes in → out | HOGs in → out |
|---|---|---|
| Top-level HOG members | 263,013 | 19,968 |
| md5 match (62,617 non-panel genes) | → 188,064 | → 19,411 |
| Not fragment | → 174,723 | → 19,183 |
| Mean pLDDT ≥ 70 | → 125,036 | → 14,619 |
| Length 0.5–2× median | → 120,516 | → 14,619 |
| ≥ 8 1:1 species | 23,564 1:1 genes | → **2,109** |
| Clade-closed selection | | → **479** |

- **Strata:** all 479 are in stratum 4. In euk20 the only panel clade with ≥ 8 of the 17 species is the root (1455 Mya); Opisthokonta has 7.
- **Pairs:** mean |S_h| is 11.6 and the split is 96 tune / 383 confirm. Primary = 30,961 pairs. **4.2% of primary pairs are labelled `out_paralog`**: 1:1 species sets do not guarantee orthology (differential loss).
- **E6:** 1,449 candidate paralog pairs → 292 fail `aa_ml_seq`, 39 fail pid → 1,118 triplets (666 with O2) in 180 paralog groups and 133 HOGs. That gives 6,724 pairs (1,118 para, 3,070 orth, 2,536 ldo), i.e. **14 E6 pairs per sampled HOG**.
- **Runtime:** 75 s wall on 8 threads; 20 s with the cache.

## 4. Projection for euk180, and decisions for the PI

**Primary pairs.** The panel clades with ≥ 8 leaves, by stratum, are:

| Stratum | Clades | C(min(n,30), 2) |
|---|---|---|
| < 100 Mya | 4 (Catarrhini 10, Schizophora 8, Boreoeutheria 14, malvids 9) | 28–91 |
| 100–300 | 7 | 28–153 |
| 300–1000 | 17 (Gnathostomata 27, Embryophyta 20, Fungi 36–40, Opisthokonta 67–69, …) | 36–435 |
| ≥ 1000 | root only | 435 |

- Stratum 1 will be under-filled, because clade order (largest first) sends any HOG that qualifies deeper to the deeper clade. Its unused quota moves to strata 2–3.
- Stratum 4 will hold only ultra-conserved HOGs.
- Expected primary ≈ 7,000 × ~150–250 ≈ **1.0–1.8M pairs**.

**E6 is the risk.** On euk20 each triplet had a median of 2 O1 and 1 O2 species, but up to 12–13. One paralog group had **191 triplets**: every pair of n in-paralogs, in each species. With 9× more species, O sets and retaining species both grow roughly 9×, so uncapped E6 is about **0.5–1k pairs per sampled HOG, i.e. several million in total**.

- Under the plan's literal rule ("scale all quotas down"), E6 would then shrink the HOG sample, possibly below ~3–4k HOGs.
- **Proposed deviation (needs sign-off and a DEVIATIONS entry):** `--e6-max-o 5 --e6-max-per-dcg 5`, i.e. at most 5 O1 and 5 O2 species per triplet and at most 5 triplets per paralog group, chosen by seeded subsets. That gives ≤ 21 pairs per triplet. On euk20 a cap of 3/3 cut E6 from 6,724 to 2,724 pairs with every paralog group kept (180/180).
- E6a and E6b average over paralog groups, so the estimand is unchanged; precision per group drops slightly.

**Run the `--projection-only` mode first.** It prints the exact eligible counts per stratum and the primary and E6 pair totals, both uncapped and capped (rerun it with the flags; the cache makes this cheap). Then run the real build once.

**Other items for the PI:**

1. **`aa_ml_seq` uses α_LG = 1.0** (placeholder), the same as the pipeline. Once §Tuning sets α, `recompute_ml.py` will change the stored `aa_ml_seq` and some selected pairs could exceed 2. Recommendation: define the selection variable at α = 1 (record as a C-item), and keep the cache file name, which encodes α.
2. **C23 says seqonly gap open = −12, but `src/strucevo/pairs/align.py` still has −11.** `seqdist` imports the pipeline aligner, so it follows any change automatically. Changing it **after** the HOG build would, however, make the stored selection distances inconsistent with the pipeline (and would fail `test_seqdist_matches_pairs_pipeline` against the old demo). **Change it before the euk180 build**, then delete the cache.
3. **Primary pairs include `out_paralog` relations** (4% on euk20). They are kept, with `relation` filled, as the plan does not exclude them. They will inflate misfit. An option would be to require all S_h pairs to be orthologs, which would be a plan change.
4. **Panel-tree deep structure is coarse** because of the CI-overlap rule (100 internal nodes; large root polytomy). This reduces the candidate clades for selection at > 1 Gya, but only the root lies there anyway.
5. Polyploid species (*Gossypium* ×4, *Sinocyclocheilus* ×3, *Glycine*, *Carassius*) add in-paralogs. Their homeologue pairs are mostly leaf-level or genus-level and are excluded by the "P ≠ leaf" rule only when P is a leaf. *Gossypium*-level P passes. Consider flagging these in E6.

## Files

- **Scripts:**
  - `scripts/hogs/{panel_tree.py, build_hogs.py, seqdist.py, treeutil.py, taxo.py, gene_flags.py}`
  - `tests/test_hogs.py`
- **Data (gitignored):**
  - `data/interim/timetree/{panel_tree.nwk, panel_tree_uncollapsed.nwk, node_ages.tsv, panel_tree_mapping.tsv, panel_tree_fixes.tsv, panel_tree_summary.json}`
  - `data/interim/hogs/gene_flags.parquet`
  - `data/interim/hogs_test_euk20/` (test build: `hogs.parquet`, `pairlists/`, `split.tsv` for the **test only**, `logs/`)
  - `data/interim/hogs_test_euk20_exh/` (exhaustive comparison)
- **Logs:**
  - `logs/filters/panel_tree.tsv`
  - `logs/hogs_gene_flags.log`
- Nothing was committed. The production `data/interim/split.tsv` has not been written.
