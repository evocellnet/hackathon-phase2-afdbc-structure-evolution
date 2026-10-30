# X4 / H1e analysis: do AlphaFold models distort ortholog structural divergence? (EXPLORATORY)

Analyst, node 2, 2026-10-30. ANALYSIS_PLAN §6 X4, DESIGN H1e + §11.10. Unadjusted 95% CIs (percentile, cluster bootstrap over top-level HOGs, B = 2000, seed SeedSequence(20261030, spawn_key=(1010,))). Nothing committed.

## Bottom line
- **AF models make orthologs look 30–40% more similar than their experimental structures do.** Same residues, same core mask, and post-cutoff chains (tdiaa, n = 197 pairs / 160 HOGs). The ratio of mean d_AF to mean d_PDB is:
  - d_3Di: 0.59 [0.52, 0.66]
  - 1 − LDDT: 0.70 [0.62, 0.76]
  - 1 − TM: 0.69 [0.60, 0.78]
  - core RMSD: 0.62 [0.51, 0.75]
- 84–94% of pairs have Δ < 0. The result is the same under tmalign, for pre-cutoff and any-date chains, and on tune+none HOGs only.
- **The bias does not grow with sequence divergence. It is mostly a constant offset that shrinks slightly at large d_AA.** The slope of Δ against d_AA is ≥ 0 for every metric:
  - d_3Di: +0.052 [+0.026, +0.073]
  - 1 − LDDT: +0.012 [+0.005, +0.018]
  - 1 − TM: +0.004 [−0.005, +0.013]
- On the same pairs, d_AF rises with d_AA at least as fast as d_PDB: the slope ratio for 3Di is 1.34 [1.15, 1.57]. AF does not flatten the structure-vs-sequence curve. It removes an additive, noise-like component that is present in experimental pairs.
  - On tune+none HOGs only (n = 79), the slopes are ~0 (3Di +0.014 [−0.034, +0.047]). On any-date chains (n = 795), they are positive again (+0.047 [+0.036, +0.057]).
  - So "no growth in the bias" is robust. "The bias shrinks with d_AA" is suggestive only.
- **The same-protein floor is large.** Matched AF model vs its own PDB chain, post-cutoff, median:
  - d_3Di 0.241 [0.232, 0.256]
  - 1 − LDDT 0.043
  - 1 − TM 0.018
  - core RMSD 0.83 Å
  - For comparison, the AF identical-sequence floor η_3Di is 0.052. Two AF models of *different orthologs* are closer in 3Di (median 0.167) than one AF model is to its own crystal structure. The PDB–PDB excess (median 0.344) is of the order of two experimental-vs-ideal deviations. This is consistent with the offset reflecting experimental/conformational noise plus AF idealisation, not lost evolutionary signal.
- **Post- vs pre-AF2-cutoff.** The conclusions are the same, but the bias is somewhat smaller with post-cutoff structures:
  - On the same 50 ortholog pairs with both kinds of chain, ΔΔ(post − pre) for d_3Di is +0.050 [+0.015, +0.087], for 1 − TM +0.022 [+0.003, +0.049], and for 1 − LDDT +0.013 [+0.002, +0.025].
  - The post floor is slightly higher than the pre floor (1 − TM 0.018 vs 0.015; 3Di 0.241 vs 0.235).
  - This is consistent with mild memorisation of training structures. It is confounded by entry age, resolution and construct.

**Implication for the project.** The structural distances in all endpoints come from AF models. Absolute structural divergences (and η, S0) are therefore about 0.6–0.7× what experimental structures give. The slope against d_AA, and so the curvature and ratio endpoints, is not compressed (slope ratio ≥ 1). The "AF homogenises homologs" artefact is real as an offset. It does not create apparent saturation at depth in this sample, but d_AA > 2 has only ~10% of pairs.

## 1. Data and pair sets (`scripts/h1e/prep_x4.py`, log `logs/h1e_x4_prep*.log`)
- **Ortholog pairs.** Ancestra euk180 1:1 orthologs from different species, labelled with HogIndex as in node-1 `run_h1e.py`.
- **Chain rule.** X-ray, ≤ 3.0 Å, observed coverage = min(SIFTS, RCSB modeled) / UniProt length ≥ 0.70, ranked by coverage, then resolution, then newest. Three chain choices:
  - `post`: released after 2018-04-30, the §11.10 primary set. It is identical to the node-1 rule and gives 2,389 accessions.
  - `pre`: released on or before 2018-04-30 (4,043 accessions).
  - `any`: any date (5,251 accessions).
- **Candidates.** 4,629 pairs of qualifying proteins: 820 1:1 orthologs, 1,695 orthologs that are not 1:1, and 2,114 out-paralogs.
- **Verification of the node-1 set.** The post set reproduces the node-1 h1e set: 201 pairs vs 204, the difference being ≥ 30 matched residues and AF availability.
- **Fetch fix.** The 5 node-1 fetch failures were a race on the shared download cache (two threads fetching the same entry, e.g. 8oi5). Fetching is now grouped by entry; 1,277/1,277 chains were fetched.
- **Pairs with both sides OK (tdiaa):**

| set | pairs | HOGs |
|---|---|---|
| post | 197 | 160 |
| pre | 485 | 358 |
| any | 795 | 510 |

- **Composition of the post set.** It is dominated by human–mouse (87) and yeast pairs (237561–559292: 21; human–yeast: 11).
- **d_AA (AA ML distance) quantiles** for post, p10/p50/p90/max: 0.014 / 0.29 / 2.2 / 4.0. The deep end is thin.

## 2. Residue matching (the key fix vs node 1)
- **Matched residue set.** For each protein × chain, M = UniProt positions where the PDB residue has N/CA/C atoms, the residue exists in the AFDB v4 model, and the PDB residue equals the AF amino acid.
- **Two structures per protein:**
  - `X_` = PDB coordinates on M.
  - `A_` = AF coordinates on M.
  - Both carry the AF pLDDT as their B-factor and enter the store as kind "af". The core mask (pLDDT ≥ 70 in both) is therefore identical, and so are the sequence, the chain breaks (which matters for 3Di), and the TM/LDDT normalisation lengths.
- **Tests.**
  - On 5 chains: identical residue lists and pLDDT.
  - Floor rows: aa_ml ≈ 0 (max 0.03) and n_aln ≈ len.
  - Median aa_ml gap between AF and PDB pairs: 0.000.
  - 1,590 chain tags were matched. 12 had fewer than 30 residues and 8 were excluded for > 5% AA mismatches. The median fraction of AF residues kept is 0.93.
- **Why matching was needed.** Without matching (node-1 style, at tuned α), the sign of ΔTM flips (median +0.009 vs −0.011 matched) and ΔLDDT is halved (−0.017 vs −0.032).
- **Pipeline settings.** `run_pairs.py` was run on all three alignments at the tuned α (AA 1.64, 3Di 1.60), with 12 threads. Output is in `data/interim/h1e_x4/pairs/`.

## 3. Outputs
- **Tables:**
  - `results/h1e/x4_paired_summary.tsv`: every alignment × chain choice × subset (all / tune+none) × metric. Columns: n, medians, mean and median Δ with CI, ratio, fraction with Δ < 0, Wilcoxon p (unclustered), and the slopes of Δ, d_AF and d_PDB against d_AA, plus the slope ratio.
  - `results/h1e/x4_floor.tsv`
  - `results/h1e/x4_post_vs_pre_same_pairs.tsv`
  - `results/exploratory/x4.json` (also copied to `results/h1e/x4_results.json`)
- **Figures (300 dpi):**
  - `figures/h1e/x4_paired_scatter_{tdiaa,tmalign}.png`: d_AF vs d_PDB with the identity line.
  - `figures/h1e/x4_delta_vs_daa_{tdiaa,tmalign}.png`: Δ vs d_AA with OLS slopes and CIs.
  - `figures/h1e/x4_floor_tdiaa.png`
- **Ledger.** 15 lines appended to RESULTS_LEDGER.md, reviewed-by "pending".
- Per-pair rows are kept in `data/interim/h1e_x4/`.

## 4. Caveats and decisions for the PI
1. **Confirm split.** 118 of the 197 post pairs (459 of 795 any-date pairs) are in confirm-split HOGs.
   - I read the separate H1e pair outputs, never `data/interim/pairs`, and CONFIRM_RUN was not set. The node-1 summary had already aggregated these pairs.
   - X4 feeds no tuning value. The tune+none-only rows are reported for every quantity, and the offset conclusions hold there.
   - Strictly, though, this is a look at AF–AF d_3Di against d_AA on 118 confirm HOGs. Please log it as a C-item, or tell me to restrict X4 to tune+none.
2. **Seed.** Stage 1010 is outside `config/seeds.yaml`, as with the other H1e seeds (integration report item 5).
3. **Paired samples.** The post vs pre comparison uses different PDB entries for the same protein. Entry age, resolution, ligand state and construct are confounded with AF2 training exposure.
4. **No PDB–PDB floor.** There is no same-protein PDB–PDB floor (two independent crystal structures). It would separate "AF idealisation" from "experimental noise" and is a cheap next step with the `pre`/`post` chain pairs of the same protein.
5. **Thin deep end.** The deep d_AA range is sparse (d_AA > 2: ~10% of pairs). The slope CIs are dominated by human–mouse and yeast–yeast pairs.
6. **Possible artefact.** Engineered constructs and complexes are not filtered. The AF–PDB floor tail (TM > 0.3, RMSD > 10 Å) probably reflects domain motions and alternative conformations.
