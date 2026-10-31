# Main figures A–E: figure and analysis engineer report (2026-10-31)

**Scope.** These figures are exploratory. They use all data (tune + confirm, `CONFIRM_RUN=1`) and were made after the single confirmation run. Nothing was committed.

**Files.**
- Script: `scripts/exploratory/figures_main.py`. Run with `CONFIRM_RUN=1 python scripts/exploratory/figures_main.py [--only A,B,C,D,E] [--B 1000]`. It takes about 30 s in total.
- Outputs in `figures/main/`:
  - `figA_shape.png` and `figA_shape_binned.tsv`
  - `figB_selection_bias.png` and `figB_selection_bias.tsv`
  - `figC_time.png` and `figC_time_binned.tsv`
  - `figD_ratio_k.png` and `figD_ratio_k.tsv`
  - `figE_terminal_internal.png` and `figE_terminal_internal.tsv`
  - `key_numbers.json`
- All PNGs are 300 dpi.

**Licence.**
- The TSVs hold bin-level aggregates only. In Fig C, a bin's mean age is given only when the bin spans at least 3 nodes.
- The node-level table behind the Fig C scatter (one TimeTree age per node) was written to the session scratchpad, not to the repository.

**Checks.**
- I recomputed the pooled k from `data/interim/branches/tdiaa*.parquet` on the confirm split, using the E3 exclusion rules (fit reason, `excluded`, or misfit > τ = 3). It gives 3.32181, which equals E3.json exactly.
- For the Fig A filters I reused `prepare_curve_data`, so they are the E1/E2 filters.

**Style.** Charts follow the dataviz reference palette in light mode, with at most 3 categorical hues per panel (blue, orange, aqua; violet for the single series in Fig C b). Text is in ink colours.

---

## Fig A: shape of structural divergence against sequence divergence (Goal 1)

**Caption draft.**
- For 1:1 ortholog pairs (all 7,000 sampled HOGs, E1/E2 filters), the panels plot structural distance against the maximum-likelihood (ML) amino-acid (AA) distance between the two proteins, in expected substitutions per site.
  - E1/E2 filters: selection distance `aa_ml_seq` between 0.05 and 2; failed and saturated pairs removed; HOGs with ≥ 10 pairs and a distance range ≥ 4×.
- **(a)** 3Di distance: the ML distance between the two proteins' Foldseek 3Di structural-alphabet strings. **(c)** LDDT distance D = −ln((S − S0)/(1 − S0)), where S is the global LDDT similarity and S0 is the LDDT of unrelated proteins of the same length.
- Points are medians in 0.1-wide AA-distance bins, and the band is the interquartile range (IQR) across pairs.
- The black curve is the confirmed fit y = η + A·x^ĉ:
  - ĉ is the confirmatory curvature: E1 ĉ = 0.72 for 3Di, E2 ĉ = 0.70 for LDDT.
  - A is the pair-weighted median of the per-HOG amplitudes a_h, re-estimated with ĉ held fixed.
- The grey dashed line is the linear reference (c = 1, refitted the same way).
- The orange dotted line is the prediction-noise floor η. This is the distance between two AlphaFold models of the *same* sequence: 0.052 for 3Di and 0.025 for D.
- An exponent below 1 means structure diverges quickly at first and then more slowly relative to sequence. The relationship is concave.
- **(b, d)** The same analysis under two other ways of aligning each pair: a sequence-only alignment and the TM-align structural superposition. Points are bin medians and lines are each alignment's own fitted curve.
- Primary-alignment sample: 741,062 pairs in 6,067 HOGs (3Di) and 744,561 pairs in 6,114 HOGs (LDDT).

**Key numbers.**

| Panel | Alignment | ĉ | A (fixed ĉ) | Linear A (c = 1) | η |
|---|---|---|---|---|---|
| 3Di | 3Di+AA | 0.723 | 0.222 | 0.247 | 0.0524 |
| 3Di | seqonly | 0.801 | 0.249 | | |
| 3Di | TM-align | 0.804 | 0.245 | | |
| LDDT | 3Di+AA | 0.699 | 0.212 | 0.235 | 0.0253 |
| LDDT | seqonly | 0.719 | | | |
| LDDT | TM-align | 0.319 | | | |

**Caveats to state.**
- ĉ was estimated on the confirm split with the tree-smoothed AA distance (D29) and a SIMEX bias correction. The curve shown here is drawn against the raw pair AA distance over all data, so it is a visual overlay, not a refit.
- With TM-align, LDDT D medians bend *upward* beyond AA distance of about 1.3, so the low ĉ = 0.32 reflects the young-pair end. The TM-align LDDT shape is not robust; only the "concave" decision is.

## Fig B: what structural pre-selection does to the curve (X7)

**Caption draft.**
- The panels plot median structural similarity, with IQR bands, against sequence identity of the aligned residues (0.02 bins) for two sets of pairs:
  - our 1:1 orthologs (blue; identity over core columns with pLDDT ≥ 70; 873,628 pairs);
  - AlphaFold DB cluster member-vs-representative pairs ("P-clu", orange; 27.7 M pairs; green dashed = eukaryote–eukaryote pairs only).
- **(a)** TM-score. **(b)** LDDT.
- **(c)** Distribution of identity in each set.
- Orthologs span 0.37–1.0 identity (5–95%: 0.47–0.99). Cluster pairs are concentrated at 0.10–0.35, because AFDB clusters are built from proteins that already passed a structural-similarity threshold.
- In the low-identity range covered only by the cluster pairs, median TM is still 0.66 at 21% identity and 0.78 at 31%. The curve is a hockey stick:
  - a steep rise up to about 35% identity;
  - then a plateau at TM ≈ 0.83.
- The low-identity end looks structurally conserved *because* membership required structural similarity.
- Where the two sets overlap (identity ≥ 0.4), orthologs are slightly *more* similar than cluster pairs:
  - TM 0.86 vs 0.83 at 51% identity, and 0.89 vs 0.83 at 71%;
  - LDDT 0.90 vs 0.85, and 0.96 vs 0.89.
- Grey zone: cluster alignments below 10% identity are short or spurious (`pclu_report` §1).

**Caveats to state.**
- The TM-scores are normalised differently: orthologs by the shorter protein, P-clu by alignment length.
- Ortholog identity and LDDT use core residues only (pLDDT ≥ 70), while P-clu uses all aligned residues. The ortholog-above-P-clu offset at high identity is therefore at least partly a confidence and disorder effect, not a selection effect.
- The two curves barely overlap. The main message is that a curve fitted on cluster pairs describes a different, structurally pre-filtered range of identity from the one where our orthologs lie.

## Fig C: divergence against time (X2)

**Caption draft.**
- **(a)** Path distance between 1:1 orthologs against the age of their last common ancestor (TimeTree median, log axis).
  - Each point is one panel node: the median over HOGs of each HOG's median distance across pairs whose common ancestor is that node.
  - Point area grows with the number of HOGs. Nodes with ≥ 5 HOGs are shown: 117 nodes, of which 6 are ≥ 1 Ga.
  - Metrics: AA ML distance (blue), 3Di ML distance minus floor η (orange), LDDT D minus η (green).
  - Lines are medians of the node values in log-age bins (0–10, 10–30, 30–100, 100–300, 300–1000 Mya).
  - Nodes ≥ 1 Ga (shaded, open circles) are shown separately and were not used for the slopes. Nodes whose floor-corrected value is ≤ 0 cannot be drawn on a log axis.
  - The dashed grey line has slope 1, which a constant-rate clock would follow.
- **(b)** Ratio of AA to floor-corrected 3Di path distance per node (ratio of HOG-summed per-HOG means). The dashed line is the confirmatory branch-based E3 estimate, k = 3.32.

**Key numbers.**
- Log–log slopes over nodes 10–1000 Mya: AA 0.85, 3Di 0.46, LDDT D 0.50.
- AA/3Di ratio, median of nodes per bin: 0.31 (< 10 My), 0.85 (10–30), 1.99 (30–100), 3.22 (100–300), 3.96 (300–1000), 3.93 (≥ 1 Ga, 6 nodes).

**Interpretation and caveats.**
- Structural distance grows roughly as the square root of time, and sequence distance close to linearly. The ratio therefore rises with age and plateaus at about 4 beyond 300 My.
- For the youngest nodes (< 30 My), the floor-corrected 3Di distance *exceeds* the AA distance. This is the regime where per-protein prediction noise is larger than the median floor η. The rise in the ratio at young ages therefore mixes biology with incompletely removed noise; Fig E shows this directly.

## Fig D: the conservation ratio k (E3)

**Caption draft.**
- k is the AA branch length divided by the 3Di branch length, both on the same HOG tree (tip-corrected). k > 1 means structure, as encoded by 3Di, changes more slowly than sequence.
- **(a)** Per-HOG k_h for 5,915 HOGs with both sums positive; 197 HOGs have a non-positive sum.
  - Median 3.24, IQR 2.38–4.34, 10–90% 1.70–6.45.
  - The pooled ratio of sums is k = 3.32.
- **(b–d)** Pooled k over all data, with 95% HOG-bootstrap CIs (B = 1000), and the confirmatory E3 values (open diamonds; confirm split only).
  - **(b)** By CATH class: α 4.67 [4.51, 4.82], α/β 3.16, mixed 3.21, β 2.55 [2.46, 2.63], no CATH 3.54. α is about 1.8× β.
  - **(c)** By HOG root age: < 100 Mya 2.19, 100–300 2.66, 300–1000 3.27, ≥ 1000 3.74.
  - **(d)** All edges 3.32 [3.28, 3.36] vs internal edges only 4.63 [4.58, 4.68].
- 6,112 HOGs (4,891 confirm + 1,221 tune) passed the E3 exclusions.

## Fig E: terminal-edge dependence (E4/E5 caveat)

**Caption draft.**
- **(a)** Pooled k against branch midpoint age (TimeTree), with 95% HOG-bootstrap CIs. Branches are split into:
  - terminal edges, which lead to an extant protein;
  - internal edges;
  - all edges (grey dotted).
- **(b)** Share of terminal edges in each age bin.
- Terminal edges have a very low k at young ages: 0.46 below 10 My, and 1.60 at 10–30 My. Their 3Di length still contains per-protein prediction noise after the tip correction. Their k rises steadily to about 4.2 by 300–1000 My.
- Internal edges are much flatter. Apart from the noisy < 30 My bins (2.6 and 3.1), they sit at about 4.5–5.0 from 30 My onward.
- Young age bins consist mostly of terminal edges (88% below 10 My), so the pooled "k increases with age" trend is largely carried by terminal edges.
- This matches the confirmatory sensitivities:
  - E4: β = −0.150 with all edges, but −0.047 (p = 0.18) on internal edges only;
  - E5: "artefact_suspect" with all edges, "no_excess" on internal edges.

**Statement for the report.** The confirmed age trend in the AA:3Di ratio (E4, β_k > 0) cannot be separated from residual AlphaFold/3Di prediction noise on terminal branches. On internal branches, k is approximately constant at about 4.6–5.0 from 30 My to beyond 1 Ga.
