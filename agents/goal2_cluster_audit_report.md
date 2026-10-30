# Goal 2-ii (X5): AFDB v3 structural clusters in the light of orthology. Report

Analyst agent, node 2, 2026-10-30 (about 21:53 to 22:40). **Exploratory analysis.** CIs are 95% and unadjusted. Every number is in RESULTS_LEDGER.md (reviewed-by = pending).

**Scripts** (`scripts/goal2/`, seed 20261030, ≤12 threads):

| Script | What it does |
|---|---|
| `g2_01_hog_extract.py` | ancestra env: OrthoXML → members + thinned labelled pairs |
| `g2_02_mapping.py` | items 1 and 3 |
| `g2_03_pairs_dist.py` | pair subsample + aa_ml_seq |
| `g2_04_tune_pairs.py` | tune-split pairs + architecture |
| `g2_05_arch_vs_time.py` | item 4 |
| `g2_06_split_model.py` | item 2 |
| `g2_07_fig_mapping.py` | figures for items 1 and 3 |

**Outputs:**
- Tables: `results/goal2/` (no ages, only age-bin labels).
- Figures: `figures/goal2/` (300 dpi).
- Interim data with ages: `data/interim/goal2/`, never committed.
- Logs: `logs/goal2_*.log`.

## Definitions
- **Eligible member.** The protein has an AFDB v4 model with md5-identical sequence (seq_match), is not a fragment (cluFlag ≠ 3), and is in a non-singleton AFDB v3 Foldseek cluster. The cluster is `afdb_cluster_rep`.
- **HOG membership.** Of the 2.92 M genes in the 109,598 top-level HOGs:
  - 2.39 M are eligible;
  - 299 k are structural singletons;
  - 75 k are fragments;
  - 163 k have no model or a sequence mismatch.
- **Sensitivity "with singletons".** Each singleton or AFDB50 group counts as its own unit.

## 1. HOG ↔ cluster mapping
71,961 HOGs have at least 2 eligible members. **70.3%** of them span more than one cluster (median 2 clusters, 90th percentile 20, 99th percentile 74).
- The member-weighted fraction of members in a HOG's largest cluster is **0.36**.
- Splitting depends strongly on HOG size: 42% of 2-member HOGs are split, and 100% of HOGs with more than 50 members.

Panel members fall into 354,716 clusters. The median is 1 HOG per cluster (90th percentile 3). **48.5%** of clusters with at least 2 members contain more than one HOG.

Contingency, as a fraction of HOGs:

| Category | Fraction of HOGs | Fraction of members |
|---|---|---|
| 1 HOG ↔ 1 cluster | 0.142 | 0.013 |
| HOG split only | 0.096 | 0.029 |
| Cluster merges HOGs only | 0.154 | 0.014 |
| Both split and merged | **0.607** | **0.943** |

With singletons as their own units: 68.7% of HOGs are split and 20.0% are 1:1.

**Reading.** AFDB clusters are not orthology units: they cut across HOGs in both directions. Splitting dominates in large HOGs.

Figure: `hog_cluster_mapping.png`. Tables: `contingency_*`, `split_by_hog_size.tsv`.

## 2. Why HOGs split (pair-level logistic regression, HOG-clustered SEs)
**Data.** Within-HOG pairs, at most 8 per HOG, both eligible, from HOGs that span more than one cluster: **279,196 pairs in 48,130 HOGs**. 75% of these pairs are in different clusters.

**Covariates:**
- CATH-string difference;
- an indicator for no CATH label on either protein;
- |Δ number of TED domains|;
- |log length ratio| (len_v4);
- mean pLDDT;
- AA divergence (see below);
- relation: ortholog, out-paralog or in-paralog.

**AA divergence** is `aa_ml_seq`: the pairs pipeline's own seqonly global alignment with LG+Γ ML, α = 1, the selection metric, via `scripts/hogs/seqdist.py`. 349 k pairs were computed fresh in 127 s, and 15 k were taken from the cache. The value is capped at 5, with a saturation flag (3.6% of pairs).

Standardized log-odds per SD, or 0→1 for binary terms:

| Term | Log-odds [95% CI] |
|---|---|
| \|log length ratio\| | **4.06 [4.01, 4.12]** |
| d_AA | 1.17 [1.15, 1.19] |
| pLDDT | −0.14 [−0.15, −0.12] |
| out-paralog vs ortholog | +0.17 [0.13, 0.22] |
| in-paralog vs ortholog | −0.26 [−0.30, −0.22] |
| \|Δ number of domains\| | 0.03 [0.01, 0.05] |
| CATH architecture differs | 0.018 [−0.010, 0.046] |

Model fit: McFadden R² 0.40, AUC 0.90.

Each block on its own (McFadden R²):

| Block | R² alone |
|---|---|
| Length | 0.305 |
| AA divergence | 0.200 |
| Architecture | 0.041 |
| pLDDT | 0.008 |
| Relation | 0.004 |

Dropping the length block costs 57,045 units of deviance; dropping the architecture block costs 40.

**Length ratio dominates, then AA divergence.** Architecture change matters only through the length change it causes:
- Descriptively, a different architecture adds 5–9 percentage points within each length bin (`split_desc_arch_len.tsv`).
- This is expected from Foldseek-cluster's 90% bidirectional-coverage criterion.

The results are the same for orthologs only and for pairs where both proteins are CATH-labelled (`split_model_summary.json`).

Figure: `hog_split_model.png`.

## 3. Cluster LCA vs HOG root (asymmetry handled explicitly)
**Method.** The cluster LCA is AFDB's own `mmseqs lca`, computed over all AFDB members, as in the paper.
- A eukaryotic LCA is projected onto the panel: the guide-tree MRCA of the panel species that descend from that taxon.
- A non-eukaryotic LCA (cellular organisms, Bacteria, …) is labelled **pre-eukaryotic**.
- A HOG root can be at most the Eukaryota root, because the panel is eukaryote-only. A pre-eukaryotic LCA paired with a HOG rooted at Eukaryota is therefore **"consistent (censored)"**: it cannot be refuted.
- Both nodes lie on the lineage of the shared member, so they are always nested. Only 0.04% are incomparable, which comes from taxonomy-version differences.

**Results:**

| Unit | Same node | Censored-consistent | Cluster LCA younger | Cluster LCA older, within Eukaryota | Pre-euk. LCA but HOG younger |
|---|---|---|---|---|---|
| Clusters vs dominant HOG (354,716) | 26.9% | 23.8% | 24.2% | 14.1% | 10.8% |
| HOGs vs largest cluster (81,786) | 18.7% | 12.7% | 7.1% | 32.0% | 29.4% |

**About half of clusters (50.7%) have an LCA consistent with the orthology-based age**, and only 26.9% are exactly consistent. Among clusters with a eukaryotic LCA, 41% are on the same node.

**Decomposition** (clusters; `lca_decomposition_clusters.tsv`):
- **Effect of non-panel taxa.** The full-AFDB LCA is older than the panel-restricted LCA for 79.8% of clusters: 45% within Eukaryota and 35% pre-eukaryotic. They are equal for 19%.
- **Effect of the partition.** The panel-restricted LCA is younger than the dominant HOG's root for 70%, equal for 22% and older for 8%.

So two effects work against each other:
- **Cluster splitting makes the LCA too young.** A cluster covers only a subclade of its HOG.
- **Taxon sampling and cross-HOG merging make it too old.** These merges come from HOG root-splits, deep paralogs and cross-kingdom homologs.

The paper's LCA ages are therefore not orthology ages.

Figure: `cluster_lca_vs_hog_root.png`.

## 4. Domain-architecture change vs time and AA distance
**Data.** Tune-split primary 1:1 ortholog pairs: 178,593 pairs in 1,400 HOGs, read from the seqonly shards with a `split == 'tune'` predicate (asserted). 86% of pairs have both proteins CATH-labelled. Ages are the panel-tree MRCA medians.

**Annotation noise floor.** Even nearly identical orthologs (d_AA < 0.02, n = 9,818) differ in architecture:

| Metric | Fraction that differs |
|---|---|
| CATH string | 0.245 |
| Topology-level (C.A.T) string | 0.215 |
| Number of domains | 0.176 |
| AFDB cluster | 0.071 |

TED chopping and label variation therefore dominate the signal. I added the C.A.T truncation to remove differences in H- vs T-level label depth.

**Against time, C.A.T architecture change is flat:**

| Age | Fraction that differs |
|---|---|
| <100 Mya | 0.239 [0.217, 0.262] |
| 300–500 Mya | 0.280 |
| >1000 Mya | 0.221 [0.170, 0.272] |

In a logit model, the age effect per 100 My is −0.043 [−0.073, −0.013] after adjusting for d_AA.

**Against AA distance, the rise is modest:** from 0.22 at d < 0.25 to 0.36 at d of 1.5–2. The logit coefficient for d_AA is +0.52 [0.33, 0.71].

**Cluster split, by contrast, rises steeply with time:**

| Age | Fraction in different clusters |
|---|---|
| <100 Mya | 0.142 [0.128, 0.156] |
| 700–1000 Mya | 0.559 [0.519, 0.601] |

**Paralogs vs orthologs at matched d_AA** (quadratic in d_AA, C.A.T architecture). Paralogs change architecture more often:
- all-HOG pair sample (212,502 pairs): log-odds +0.23 [0.20, 0.26];
- tune E6 in-paralogs vs E6 orthologs: +0.15 [0.06, 0.24].

Figure: `arch_change_vs_time_and_daa.png`.

## Caveats
- **Pair sampling.** Labelled pairs were thinned to ≤300 per HOG, and then ≤8 per HOG. HOGs with more than 600 genes (435 HOGs, 516 k genes) were not used for the pair analyses (items 2 and 4b), only for the mapping.
- **d_AA.** It is a global-alignment distance. For pairs that differ in architecture, it reflects only the aligned part.
- **Model specification.** The model is pair-level, so splitting is modelled conditionally on the HOG. Length enters linearly in logit; the near-separation above |log ratio| 0.4 is real (99.9% of those pairs are in different clusters).
- **LCA source.** The AFDB LCA comes from file 2. The panel-restricted LCA uses all eligible panel members of the cluster, including AFDB50 members, which explains the 0.7% where the full LCA is younger than the panel-restricted one.
- **No confirm-split rows were read.**
