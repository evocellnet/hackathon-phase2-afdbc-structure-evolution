# Ancestra: installation, characterisation and scaling report

Author: methods-engineer agent (Opus 5.5), 2026-10-30, about 09:24-11:00 CET, node of SLURM job 66670604.
Everything here was run on the project node, which other agents were also using. The 1-minute load was 5-42 during the runs, so wall times are indicative only.

## TL;DR

- **Ancestra works out of the box.** The pinned conda environment built in 80 s. The unit tests gave 147 passed and 6 skipped. The toy example reproduces the README numbers exactly.
- **What it is.** Ancestra is a graph-based, bottom-up method that builds hierarchical orthologous groups (HOGs) along a rooted guide species tree, in the lineage of SYNERGY, Hieranoid and FastOMA. It does not reconcile full gene trees: it builds small local FAMSA+FastTree trees only for ambiguous clusters and does not save them. It uses sequences only. Structure enters only through an optional 3Di "family rescue", which is off by default and lowered precision in the authors' benchmark. It never uses a GPU.
- **Output.** OrthoXML with nested `orthologGroup` (speciation, with a TaxRange) and `paralogGroup` (duplication, with the level directly below the duplication branch). From this, **every pair of proteins in a HOG can be labelled ortholog or paralog, with the species-tree node of the speciation or duplication.** I wrote and validated a script for this (`scripts/ancestra/hog_pair_labels.py`). No branch lengths or gene trees are output.
- **Accuracy we observed.**
  - Ortholog calls are reliable. Human–mouse vs Ensembl Compara: 148/149 Ensembl one-to-one orthologs recovered. Of the 187 human–mouse orthologs Ancestra called for 200 sampled genes, none is an Ensembl paralog.
  - **Duplication timing is noisy.** Only about 36% of human paralog pairs have a duplication interval compatible with Ensembl's. Errors go both ways (35% younger, 29% older). Clear errors include MAPK1/MAPK3, ACTB/ACTG1 and HOXA1/HOXB1.
  - Duplications older than the guide-tree root, and deep "merge splits", both show up as *separate root HOGs*, and the two cannot be told apart.
- **Cost.** The cost sits in within-family all-vs-all alignment (DIAMOND, roughly quadratic in family size; about 2e-4 CPU-s per pair). Measured runs:

  | Input | Wall time | CPU time | Peak RAM |
  |---|---|---|---|
  | One family: 14,358 tubulins, 1,000 species (default) | 18.5 min | 7.6 CPU-h | 17 GB |
  | Same family, single-family table, `fast` preset | 13.4 min | 5.3 CPU-h | 28 GB |
  | 20 full eukaryotic proteomes (357k proteins) | 17 min | 6.1 CPU-h | 8.6 GB |

  The default workflow splits a family above 7,500 genes and then re-joins it, which aligns it twice. Running with a precomputed single-family table avoids this and saved 28% of the time (see §3).
- **Recommendation for the study.** Run Ancestra **once on whole reference proteomes** of a few hundred species. The authors ran 1,000 proteomes / 4.75 M proteins in 10.6 h, 64 CPU-h and 221 GB, which fits this 500 GB node. Then join the HOGs to AFDB cluster membership through UniProt accessions. Running one Ancestra job per AFDB cluster is not viable for 2.3 M clusters (§5).

---

## 1. What Ancestra is

Source: <https://github.com/evocellnet/ancestra-mvp-ai>.

- **Commit:** `ffb0bbef28bc31d05755264bbe0094e4c0d4a5fa` (2026-10-27). This is the only branch (`main`); there are no tags in the public repository.
- **Code:** `ancestra/` is a byte-identical copy of the frozen method (`freeze-v1`) of an AI-agent hackathon project (Beltrao lab). It is not peer reviewed and its results are not human-verified. The README and `docs/DISCLOSURES.md` list protocol deviations.
- **Report:** a technical report PDF is in `technical_report/`. No external paper or preprint is referenced; it cites SYNERGY, Hieranoid, OrthoFinder, OMA/FastOMA, QfO and OrthoBench as prior work.

### Method

The method is graph-based, progressive and runs per level of the species tree. These are the pipeline stages (`ancestra/pipeline.py`):

| Stage | What it does | Tools |
|---|---|---|
| `input` | Proteomes go in as one FASTA per species (file stem = species = leaf of the guide tree). The rooted guide tree has named internal nodes, which become the HOG levels; leaves without a proteome are pruned. | – |
| `families` (S0) | Reference-free homology families:<br>• DIAMOND cluster (30% id / 50% cover)<br>• all genes vs representatives (`--sensitive`, e ≤ 1e-5, ≥30% coverage of the shorter sequence)<br>• connected components<br>• Leiden split of components larger than 7,500 genes<br>Can be replaced by `--families TSV`, e.g. AFDB clusters. | DIAMOND 2.2.8, igraph, leidenalg |
| `scores` (S2) | Within-family all-vs-all (DIAMOND block-vs-block, `-k 0 --sensitive`). MMseqs2 re-aligns the pairs DIAMOND missed, in families of 1,000 or fewer. Families above 5,000 keep only the top 2 hits per other species. A "cross-family rescue" re-joins S0 families linked by a strong hit and re-scores them. | DIAMOND, MMseqs2 18 |
| `merge` (S3) | Bottom-up over the species tree, for all families at once. At each node, the ancestral genes of the child clades are merged:<br>• R1: mutual-best matches (tolerance 0.95)<br>• R2: in-paralog absorption (InParanoid rule)<br>• R3: fission repair<br>Ambiguous many-to-many clusters are resolved by a local tree labelled by species overlap: FAMSA + FastTree `-lg` in `accurate`, neighbour joining in `fast`. Optional gene-order (synteny) evidence. Numba engine. | FAMSA 2.4.1, FastTree 2.2.0 |
| `output` | `hogs.orthoxml` (OrthoXML 0.5, XSD-validated, loads in pyham), `pairs.tsv` (ortholog pairs only, QfO flattening), `orthobench_ogs.txt`, `manifest.json` (versions, md5s, per-stage time and RSS). Every stage is checkpointed and skipped on rerun. | lxml |

Answers to the questions in the brief:

- **Tree-based?** Only locally. There is no full gene-tree inference and no reconciliation.
- **Reconciliation with a species tree?** Yes, implicitly. The species tree drives the merge order, and each merge event is labelled speciation or duplication relative to it.
- **Structure-aware?** Essentially no.
  - `--family-rescue 3di --structure-3di FASTA` uses Foldseek 3Di strings to attach small orphan families. It is off by default; on OrthoBench test it gave F -5.0 [-13.1, 0.0].
  - A `foldseek3di` scoring backend exists only as a stub that raises `NotImplementedError`.
- **GPU?** None. `grep -ri "gpu\|cuda" ancestra/` returns nothing. The README mentions a `v2-gpu` branch (CUDA Smith-Waterman, 1.18x speed-up) in a private research repository; it is not in the public repository.

### Inputs

| Input | Format |
|---|---|
| Protein sequences | FASTA, one file per species, `.gz` allowed. `sp\|ACC\|NAME` headers are reduced to ACC. IDs must be unique across species. One protein per gene is expected. |
| Guide tree | Rooted Newick; multifurcations are merged k-way. Can be built from NCBI taxonomy with `ancestra.speciestree`. |
| Optional | `--families` gene→family TSV; `--gene-order` parquet tables; 3Di FASTA. |

Structures are not an input. I added a wrapper mode that extracts amino-acid and 3Di sequences from PDB/mmCIF with Foldseek.

### Outputs relevant to us

- **Pair relations.** Orthologs and paralogs, derived from the HOG nesting.
- **Speciation node of each ortholog pair.** This always equals the MRCA of the two species (checked: 100%).
- **Duplication node of each paralog pair.** The paralogGroup level P is the node directly below the duplication branch. Pairs in different root HOGs are not related in the file.
- **Not output:** gene trees, branch lengths, support values, or ancestral sequences.
- **Pairwise similarities.** Pairwise alignment statistics (bits, %id, coverage) for all scored within-family pairs are in `stages/scores/scores.npz`. This can serve as a ready-made sequence-divergence matrix to sit beside structural distances. Integer gene IDs map to accessions through `stages/input/genes.tsv`.

### Scaling (authors' numbers)

| Run | Proteins | Wall time | CPU time | Peak RAM | Notes |
|---|---|---|---|---|---|
| QfO 78 species | 989k | – | 18.5 CPU-h (`accurate`) | – | – |
| 1,000 UniProt proteomes, `fast` | 4.75 M | 10.6 h | 64 CPU-h | 221 GB | Exponent 1.35 in CPU-h vs proteins |

S0 and S2 dominate at that scale.

## 2. Installation (reproducible, about 2 min)

```bash
cd /scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution
mkdir -p tools && git clone https://github.com/evocellnet/ancestra-mvp-ai tools/ancestra-mvp-ai
git -C tools/ancestra-mvp-ai checkout ffb0bbef28bc31d05755264bbe0094e4c0d4a5fa
mamba env create -p conda_envs/ancestra -f tools/ancestra-mvp-ai/environment.yml -y   # 80 s, pinned versions
conda_envs/ancestra/bin/pip install -e tools/ancestra-mvp-ai --no-deps
export PATH=$PWD/conda_envs/ancestra/bin:$PATH          # diamond, mmseqs, famsa, FastTree, foldseek must be on PATH
cd tools/ancestra-mvp-ai && python -m pytest -q -m "not slow" -p no:cacheprovider    # 147 passed, 6 skipped, 57 s
bash example/run_example.sh 4                                                          # 22 s, 0.43 GB
```

Pinned tool versions: Python 3.11.18, DIAMOND 2.2.8, MMseqs2 18.8cc5c, FAMSA 2.4.1, FastTree 2.2.0, Foldseek 10.941cd33, pyham 1.1.11. NCBI taxdump comes from `data/raw/ncbi_taxonomy/taxdump.tar.gz` (md5 a4ea5ad4…), extracted to `data/interim/ancestra/taxdump/`.

**Tests.** 147 passed, 6 skipped. The README says 148/5; the extra skip is `test_geneloss` because `miniprot` is not in `environment.yml`, which is harmless. The other skips need the benchmark data (QfO webservice, OrthoBench, S4 3Di model).

**Toy example.** Simulated data, 6 species, 1,828 proteins, exact truth. Ortholog-pair precision/recall/F1 was 0.9976/0.9122/0.9530 (`accurate`) and 0.9979/0.9148/0.9545 (`fast`), identical to the README. The OrthoXML verifier passed XSD, gene-reference integrity, TaxRange = LCA, pair count and pyham.

## 3. Runtime measurements on this node (48 threads unless stated)

`scripts/ancestra/summarize_runs.py` produces `results/ancestra/runtime_table.tsv`. Wall and CPU time are from `/usr/bin/time -v` on the whole `ancestra run`; peak RSS is the largest single process.

| Run | Preset | Threads | Species | Proteins | S0 families (max size) | Wall | CPU-h | Peak RSS | S0 s | S2 s | Merge s (of which local trees) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Toy example | accurate | 4 | 6 | 1,828 | – | 22 s (both presets + checks) | – | 0.43 GB | 5.9 | 2.4 | 1.2 |
| Tubulin PF00091, 50 species | accurate | 8 | 50 | 828 | 3 (822) | 16 s | 0.02 | 0.4 GB | 4.5 | 9.9 | 0.6 |
| Tubulin, 200 species | accurate | 8 | 200 | 3,320 | 1 (3,320) | 2.7 min | 0.28 | 2.9 GB | 8.1 | 148 | 4.8 |
| Tubulin, 1,000 species | accurate | 48 | 1,000 | 14,358 | 4 (6,250) → rescue rejoined | 18.5 min | 7.6 | 17.3 GB | 24 | 1,062 (DIAMOND 627 + rescue 300) | 20 (11) |
| Tubulin, 1,000 species, `--families` one family | fast | 48 | 1,000 | 14,358 | 1 (14,358) | 13.4 min | 5.3 | 28.1 GB | 0 | 784 (DIAMOND 670; cap dropped 68 M of 92 M pairs) | 10 (0.2) |
| 20 eukaryotic reference proteomes (UniProt 2026_03) | accurate | 48 | 20 | 357,201 | 67,203 (7,432) | 16.8 min | 6.1 | 8.6 GB | 216 | 334 | 443 (419) |
| Structures plumbing (6 PDB tubulin entries, 19 chains) | accurate | 4 | 4 | 19 | 3 | 6 s | – | 0.2 GB | 4.2 | 0.6 | 0.6 |

Pair labelling (`hog_pair_labels.py`) is single-threaded:

- euk20: 12 M labelled pairs in 63 s, 0.4 GB.
- tubulin 1,000 species with `--families-from-run`: 47 M pairs in 311 s.

### Datasets

- **Proteomes.** UniProt reference proteomes, current release 2026_03, downloaded by `scripts/ancestra/download_qfo_species_proteomes.sh` for the 78 leaves of Ancestra's QfO-2022 tree.
  - 72 were found. Six UP IDs are no longer in the current release; Ancestra simply prunes those leaves.
  - The 20 eukaryotes used span human, mouse, chicken, *Xenopus*, zebrafish, fly, mosquito, worm, *Nematostella*, *Monosiga*, three yeasts, *Batrachochytrium*, *Arabidopsis*, rice, *Physcomitrium*, *Chlamydomonas*, *Dictyostelium* and *Plasmodium*.
- **Family.** All non-fragment UniProtKB entries with Pfam PF00091 (tubulin) in eukaryotic reference proteomes, release 2026_03: 44,936 proteins in 3,114 species. Nested random species subsets of 50, 200 and 1,000 were drawn (seed 11) with `scripts/ancestra/make_family_subsets.py`.
- **Guide trees.** Induced NCBI taxonomy trees.

### Bottleneck analysis

1. **Single family.** About 95% of the time is S2: DIAMOND `--sensitive -k 0` all-vs-all inside the family, roughly quadratic in size. Measured CPU per aligned pair is about 1.9e-4 CPU-s (single-family run: 19,185 CPU-s for 103 M pairs). Memory also grows with the number of hits: a single 14k family peaked at 28 GB, against 17 GB when S0 split it in four. From 3,320 to 14,358 members, CPU time rose 27x for 19x more pairs. S0 and the merge are negligible (24 s and 20 s for 14k genes and 1,000 species).
2. **Double alignment of big families (inefficiency).** S0 splits components larger than 7,500 with Leiden. The S2 cross-family rescue then sees strong hits between the parts, joins them and re-scores them all. That is a second all-vs-all of 7,401 genes (300 s here). It also means the rescue undoes the size cap. For per-family runs, pass `--families` (one family) to skip both.
3. **Whole proteomes.** Time splits roughly evenly between S0 (21%), S2 (33%) and merge (44%). In `accurate`, 95% of the merge time is FAMSA+FastTree local trees; the authors measured that `fast` (neighbour joining) saves about 30–40% of CPU-h. Peak memory is in S2 (8.6 GB here; 221 GB at 4.75 M proteins).
4. **GPU.** None of the stages uses the GPU. GPU help would require code changes, for example a GPU MMseqs2 prefilter or Foldseek-GPU for a structural score backend. The authors' own GPU branch gave only 1.18x and is not public. **Do not plan on the A100 for Ancestra.**

## 4. Labelling protein pairs as ortholog / paralog, with event node

Script: `scripts/ancestra/hog_pair_labels.py OUT/hogs.orthoxml pairs.tsv.gz [--families-from-run OUT] [--between-species-only]`. It uses Ancestra's own OrthoXML reader and the QfO/pyham "deepest common group" (DCG) rule.

| Column | Meaning |
|---|---|
| `relation` | `ortholog` (DCG is an orthologGroup, different species)<br>`in_paralog` (DCG is a paralogGroup, same species)<br>`out_paralog` (DCG is a paralogGroup, different species, so the duplication predates their speciation)<br>`inconsistent` (same species under an orthologGroup; QfO drops these)<br>`root_split` (optional: same S0/S2 family but different root HOG) |
| `event_level` | For orthologs: the speciation node (TaxRange). For paralogs: the node P directly below the branch carrying the duplication. |
| `dup_branch_from` | parent(P) in the guide tree, so the duplication lies on the branch parent(P) → P. |
| `enclosing_og_level` | Level of the enclosing orthologGroup (the duplication is younger than that speciation). |
| `species_lca`, `lca_ok` | MRCA of the two species; whether TaxRange == MRCA. |
| `one_to_one` | 1:1 orthologs: each side of the DCG holds exactly one gene of the pair's species. |
| `dcg_id`, `hog` | IDs of the DCG and the root HOG. Pairs from the same duplication event share a `dcg_id`. |

### Validation

**Toy simulation with known truth** (`results/ancestra/toy/validation.txt`, script `validate_labels_toy.py`).

- My ortholog and paralog sets equal `ancestra.orthoxml.hogs_to_pairs` exactly, and the ortholog `event_level` equals the species MRCA in 100% of pairs.
- Against the simulator truth:

  | Label | Precision | Recall |
  |---|---|---|
  | Ortholog | 0.998 | 0.912 |
  | Paralog | 0.962 | 0.869 |

- The duplication level is correct for 95.9% of pairs that are paralogs in both. This is an easy simulation.

**Real data: human/mouse in euk20 vs Ensembl Compara** (gene-tree reconciliation; 200 random human genes via the REST API; `results/ancestra/ensembl_check/`, scripts `compare_ensembl_compara.py` and `ensembl_level_compat.py`).

- **Human–mouse orthologs.**
  - Ensembl one-to-one: 148/149 called ortholog by Ancestra, and 1 not in the same HOG.
  - Ensembl one-to-many: 9/14 ortholog, 3 out-paralog.
  - Ensembl many-to-many: 13/46 ortholog, 22 out-paralog, 11 separate HOGs.
  - In reverse, of the 187 human–mouse orthologs Ancestra called for the sampled genes, 170 (91%) are Ensembl orthologs, 0 are Ensembl paralogs, and 17 could not be mapped.
  - **So ortholog calls are precise; recall is lower in many-to-many families.**
- **Human paralog pairs listed by Ensembl** (2,249).
  - 855 are paralogs in the same Ancestra HOG.
  - 1,394 (62%) fall in *different root HOGs*. Most are Ensembl "other_paralog" pairs at Opisthokonta or Bilateria, i.e. superfamily-level paralogy that Ancestra leaves unconnected.
- **Duplication timing** for the 855 shared paralog pairs. Both sides were converted to intervals on the human lineage, because Ensembl and our 20-species tree sample species differently:

  | Verdict | Share |
  |---|---|
  | Compatible | 35.8% |
  | Ancestra younger | 35.4% |
  | Ancestra older | 28.8% |

  Exact node agreement before the interval conversion was 21%.
  - Part of the "older" set reflects Ensembl's sparse outgroups: it has one yeast and no cnidarian or choanoflagellate, so "Bilateria" in Ensembl absorbs older duplications.
  - Part of the "younger" set comes from Ancestra placing old duplications late, when a paralog pair looks species- or clade-specific.
  - Spot checks:

    | Pair | Ancestra | Expected | Verdict |
    |---|---|---|---|
    | HBA/HBB | Duplication below Euteleostomi | Same | Correct |
    | MAPK1/MAPK3 | Below Euarchontoglires | Vertebrate 2R | Too young |
    | ACTB/ACTG1 | Human-specific | Vertebrate | Too young; near-identical sequences |
    | HOXA1/HOXB1 | Below Opisthokonta | 2R | Too old; homeobox hidden paralogy |
    | TUBA1A vs TUBB | Separate root HOGs | Pre-LECA duplication | Correct |
    | CDK1 vs yeast CDC28/cdc2 | 1:1 orthologs at Opisthokonta | Same | Correct |

**Tubulin family, 1,000 species: biological sanity.**

- Root HOGs at Eukaryota separate α (4,496 members, 950 species), β (3,748), γ, δ and ε, as expected from the pre-LECA duplications.
- There is fragmentation:
  - 83% of α and 71% of β are in their largest HOG.
  - γ is split into two large HOGs at Eukaryota.
  - Pezizomycotina α forms its own root HOG.
  - 4% of proteins are root singletons.
- Consequence: between subfamilies, `root_split` mixes true ancient paralogs with fission artefacts.

### How to use this for the structural-divergence study

- **Orthologs:** use `relation == ortholog`, preferably `one_to_one == 1`. The divergence time proxy is the speciation node (`event_level`), mapped to TimeTree ages. This is the robust part.
- **Paralogs:**
  - **In-HOG paralog pairs** are reliable as paralogs.
  - **Duplication timing:** treat `event_level` / `dup_branch_from` as a noisy estimate. Use it only in coarse bins (e.g. lineage-specific vs. deep). Validate a subset against Ensembl Compara or PANTHER, or re-time duplications with a proper gene tree (IQ-TREE/GeneRax on the HOG) for the families we analyse in depth.
  - **Root-split pairs** within one AFDB cluster are "paralogs older than the root, or method fission". Label them `ancient_or_unresolved`, not paralog.
- **Least-diverged ortholog analysis.** Ancestra gives the orthology relation, but not which paralog diverged more. Compute that from sequence distances (`scores.npz`: %id and bits for every within-family pair) or from structural distances (Foldseek/TM-score) per paralogGroup.

## 5. Extrapolation to the planned study

Cost model for one-family runs, from the measured S2 cost:

```
CPU-s ≈ (1.9–2.5)e-4 × n²/2 + about 10 s overhead
RAM grows about as n²: 28 GB at n = 14k
```

| Scenario | Estimate (one Ancestra run per family) |
|---|---|
| (a) 1 family, 1,000 members, about 1,000 species | about 3 CPU-min; well under 1 min wall at 48 threads; < 2 GB |
| (a) 1 family, 10,000 members, about 1,000 species | about 2.6–3.5 CPU-h (up to 1.4x more if the Leiden split + rescue double pass triggers; avoid with `--families`); about 6–8 min wall; about 14 GB |
| (a) 1 family, 45,000 members (all eukaryotic tubulins, 3,114 species) | about 55–70 CPU-h; 1.5–2 h wall; S2 memory about 250–300 GB if RAM scales with n² (28 GB at 14k), which is at the edge of the node. Not run. Subsample species or split by subfamily. |
| (b) All AFDB clusters v3: 2.3 M clusters, 30 M members, Σn² = 2.4e10 | Alignment alone: about 1,300–1,700 CPU-h. Per-run overhead × 2.3 M runs is prohibitive (about 6,000 CPU-h). The 26 clusters with more than 10k members contribute a third of Σn², and the largest (39k) would need about 200 GB. **Not feasible as per-cluster runs.** |
| (b) Thousands of families, batched (`--families` = AFDB cluster IDs, proteomes = union of members split by species) | This is Ancestra's normal multi-family mode, and the cost becomes that of a proteome-scale run. The authors' 1,000-proteome run took 64 CPU-h and 221 GB; ours with 20 proteomes took 6 CPU-h. |

The single-family estimates use the measured 1.9e-4 CPU-s per pair from the single-family run in §3, with 2.5e-4 as an upper bound for runs that hit the double pass.

### Recommended design

1. Pick about 200–1,000 UniProt reference proteomes that are spread across the tree of life and have AFDB models. Use canonical proteomes only, because Ancestra expects one protein per gene; KW-1185 family exports include TrEMBL redundancy that inflates in-paralogs.
2. Run Ancestra once on the whole proteomes. Expected cost at 48 threads:

   | Proteomes | Wall time | RAM | Preset |
   |---|---|---|---|
   | 250 | about 2 h | about 25 GB | accurate |
   | 1,000 | about 6–10 h | about 220 GB | fast |

3. Label pairs with `hog_pair_labels.py`.
4. Intersect with AFDB cluster membership (`5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz`) to get ortholog and paralog pairs inside and across structural clusters.
5. Compute structural distances only for the pairs needed (Foldseek on GPU or CPU).

This keeps orthology inference in the regime the method was benchmarked in (QfO, OrthoBench).

An alternative that stays per-cluster: give Ancestra precomputed AFDB-cluster families with `--families` and batch about 1,000–10,000 clusters per run. This works, but the species tree becomes huge (AFDB covers about 1 M taxa) and accuracy has not been benchmarked in that regime.

## 6. Bugs, limitations, workarounds

1. **Leiden split + rescue double-scoring.** For families larger than 7,500, S0 splits the family and S2 re-joins it and re-scores it (observed: +300 s, 28% of the run). Workaround: `--families` with precomputed families (measured: 13.4 min instead of 18.5 min, with 226 instead of 216 root HOGs), or `--config` with `{"scoring": {"family_rescue": false}}`; the latter is untested here.
2. **Cap rule for families above 5,000.** The cap (top 2 hits per other species) drops 47% of the scored pairs at 14k. All pairs are still aligned first, so the cap saves no time; it only thins the graph.
3. **Ancient duplications vs fissions are indistinguishable.** Both become separate root HOGs. Paralogs older than the guide-tree root are never labelled. Use `--families-from-run` (`root_split`) and treat these pairs as unresolved.
4. **Duplication-level accuracy on real data is modest** (about 36% interval-compatible with Ensembl). Errors are typical for near-identical paralogs (gene conversion, e.g. ACTB/ACTG1) and short or repetitive domains (homeobox).
5. **No gene trees or branch lengths.** Local trees are internal and discarded. For rate-based analyses, build trees per HOG separately (FAMSA + IQ-TREE/FastTree, already in the environment).
6. **Output depends slightly on Newick child order** (authors' disclosure 11). Keep the guide tree fixed and record it.
7. **Needs a rooted guide tree.**
   - NCBI multifurcations become k-way merges.
   - A wrong tree costs accuracy; the authors' random-tree control lost 16 F points.
   - When a species' taxid is an ancestor of another sampled taxid (e.g. a species plus its strain), `ncbi_tree` still works because the leaf hangs off the internal node. This was not exercised here.
8. **AFDB file server.** It returned 403/429 after a 16-way parallel burst of `alphafold.ebi.ac.uk/files/...` requests; the API endpoint stayed reachable. Use bulk sources (Foldseek AFDB databases, Google Cloud) or the data agent's downloads. The structure-input mode of my wrapper was therefore only tested on 6 PDB entries (plumbing).
9. **Six QfO UniProt proteome IDs are not in release 2026_03** (UP000000561, UP000001055, UP000001542, UP000001555, UP000002254, UP000002521). They are pruned automatically.
10. **Minor.**
    - `miniprot` is missing from the environment (one test skipped).
    - The README test count is off by one.
    - `pairs.tsv` lists orthologs only; use my script for paralogs.
    - The full all-pairs label table with `root_split` is O(Σ family size²). It is fine for one family and huge for proteomes, so it is off by default in `run_ancestra.py`.

**Fallbacks** (not needed, since Ancestra works): OrthoFinder 3 (gene trees plus duplication nodes with support) or species-overlap labelling of FastTree/IQ-TREE gene trees per AFDB cluster against NCBI taxonomy. OrthoFinder would also provide resolved gene trees, which Ancestra lacks. For the duplication-timing analyses, a gene-tree reconciliation step on a subset of families is advisable anyway.

## 7. Files

| Path | Content |
|---|---|
| `tools/ancestra-mvp-ai/` | Clone at `ffb0bbe`; `example/out/` holds the toy outputs |
| `conda_envs/ancestra/` | Environment |
| `scripts/ancestra/run_ancestra.py` | Wrapper: `--proteomes DIR` / `--fasta FILE` (taxid from `OX=` or `--taxmap`) / `--structures DIR --taxmap`; builds the NCBI guide tree, runs Ancestra under `/usr/bin/time`, writes `pair_labels.tsv.gz` and `timing.json`. Options: `--families`, `--preset`, `--root-split`, `--extra` |
| `scripts/ancestra/hog_pair_labels.py` | Pair labeller (§4) |
| `scripts/ancestra/validate_labels_toy.py` | Validation on toy truth |
| `scripts/ancestra/compare_ensembl_compara.py`, `scripts/ancestra/ensembl_level_compat.py` | External check against Ensembl Compara |
| `scripts/ancestra/summarize_runs.py` | Runtime table |
| `scripts/ancestra/download_qfo_species_proteomes.sh`, `scripts/ancestra/make_family_subsets.py` | Data preparation |
| `data/raw/ancestra_test/` | `proteomes/` (72 UniProt reference proteomes), `families/PF00091_euk_refprot.tsv.gz` (UniProt 2026_03), `ensembl/` (UniProt→Ensembl maps) |
| `data/interim/ancestra/` | `taxdump/`, `euk20/` (symlinks), `families/` (tubulin subsets), `structures_test/` (6 PDB entries) |
| `results/ancestra/runs/<run>/` | `ancestra/` (`hogs.orthoxml`, `pairs.tsv`, `manifest.json`, `stages/`), `pair_labels.tsv.gz`, `ancestra.log`, `guide_tree.nwk` |
| `results/ancestra/runtime_table.tsv` | Runtime table |
| `results/ancestra/toy/` | Toy validation |
| `results/ancestra/ensembl_check/` | Ensembl comparison |
| `logs/ancestra_*.log` | Logs |

Example invocations:

```bash
P=conda_envs/ancestra/bin/python
# one family / AFDB cluster (FASTA with OX= taxids), skip S0 by declaring one family:
grep ">" fam.fasta | awk '{print substr($1,2)"\tfam"}' > fam.tsv
$P scripts/ancestra/run_ancestra.py --fasta fam.fasta --families $PWD/fam.tsv --preset fast --threads 48 --root-split --out results/ancestra/runs/famX
# whole proteomes (UP..._taxid.fasta.gz files):
$P scripts/ancestra/run_ancestra.py --proteomes data/interim/ancestra/euk20 --threads 48 --out results/ancestra/runs/euk20_accurate
# structures (AFDB PDB/mmCIF) + taxmap "<AF id or accession>\t<taxid>":
$P scripts/ancestra/run_ancestra.py --structures DIR --taxmap map.tsv --out OUT
```

**Storage note.** `results/ancestra/runs/*/ancestra/stages/` holds about 8 GB of checkpoints (score tables, merge pickles). Do not commit these; keep only `hogs.orthoxml`, `manifest.json`, `pair_labels.tsv.gz` and `timing.json`, or delete the stages after use.
