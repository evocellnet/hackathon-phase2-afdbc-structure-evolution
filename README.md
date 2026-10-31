(Human-written text starts)
This is an exercise in assessing the bioinformatics capabilities of agentic AI. It is part of a hackathon run by the Beltrao lab at ETH Zurich, where we used Claude to run bioinformatics projects (LINK). All results in this project were generated entirely by AI with minimal human intervention. This is not intended to be viewed as a traditional publication.

[Project-specific text] The goal of this specific project was to investigate protein 3d structure evolution. It was executed by Claude Opus 5.5 with delegation to Opus 5.5/Sonnet 5.5 subagents, and finished within 48 hours on a computational cluster with 48 cpus, 500 RAM, and a GPU.

Currently none of the results have been verified by a human. We plan to verify the agent's methodology and benchmarking, and we will update this section accordingly.

The rest of this repository is entirely AI generated. It includes code for implementation of the research and a AI-written report in the style of a journal article, detailing the final method, the path followed to reach it, and benchmarking against other methods.
(Human-written text ends)

---

# Is protein structure a molecular clock? AFDB-cluster structural divergence, revisited with orthology and divergence times

**Read first:** [`technical report/report.md`](technical%20report/report.md), the full write-up with figures.

## What was done

We revisited two papers:

- **Illergård, Ardell & Elofsson 2009.** Their finding: protein structure is "3–10× more conserved than sequence" and changes linearly with sequence divergence.
- **Barrio-Hernandez et al. 2023.** The AlphaFold Database (AFDB) structural clusters.

For this we built a new pipeline:

- **Species panel:** 180 eukaryotic proteomes.
- **Structures:** about 3.3 M AFDB v4 models.
- **Orthology:** ortholog families inferred from sequence alone with Ancestra (109,598 hierarchical orthologous groups, HOGs).
- **Divergence times:** a dated species tree from TimeTree.
- **Pairs:** 1.4 M protein pairs from 7,000 sampled families, each scored for:
  - amino-acid divergence (ML, LG+Γ);
  - Foldseek 3Di structural-alphabet divergence (ML, Q.3Di.AF+Γ);
  - 3D divergence (LDDT, TM-score).
- **Pre-registered tests:** six endpoints, tuned on 20% of families and confirmed once on the other 80%.

## Main results (details, caveats and figures in the report)

| Question | Pre-registered result (confirmation set) | Post-hoc reading (independently verified) |
|---|---|---|
| Is structural divergence linear in sequence divergence? | Concave: curvature exponent ĉ = 0.72 (98% CI 0.64–0.80), robust to alignment | The estimator is biased toward concavity by protein-specific AlphaFold noise. Corrected estimates span 0.73–0.98. Mild saturation at most; linearity is not excluded |
| How much more conserved is 3Di than sequence? | k = 3.3 (3.09–3.57); α-helical 4.7 vs β-rich 2.6 | Internal branches: k = 4.6 (4.35–4.93) |
| Does the ratio change with time (clock)? | Rises with branch age; robust to timescale and α | Mostly a terminal-branch / tip-correction artefact. Internal branches are flat from 30 My to over 1 Gy, consistent with structure following a rescaled sequence clock |
| Lineage-specific structural rates (neutral-clock null)? | "Artefact-suspect": the excess disappears on internal branches; power only 0.24 | No robust evidence for lineage-specific structural rates; neutrality is not tested by this |
| Duplication (positive control) | Detected; the more-diverged copy shows extra structural change, replicating the direction of PMID 40813247 | The paralog-vs-ortholog part is not alignment-robust |
| AlphaFold artefact (exploratory) | AF models make orthologs 30–40% more structurally similar than experimental structures do (constant offset) | — |
| AFDB clusters vs orthology (exploratory) | 70% of ortholog groups with clustered members are split across clusters, mainly by length/coverage | — |
| Case studies (exploratory) | Cytoskeleton (actin/tubulin): no robust structural bursts. Kinases: weak excess at the tyrosine-kinase origin | — |

## Repository map

| Path | Content |
|---|---|
| `PROMPT.md` | the original brief (verbatim) |
| `RESEARCH_TREE.md` | research record: milestones, decisions, status (start here to resume) |
| `LITERATURE.md` | verified literature and dataset review |
| `DESIGN.md` | study design v3 (three independent review rounds; `reviews/round*_*.md`) |
| `ANALYSIS_PLAN.md` | frozen, executable analysis plan (pre-registration); freeze md5s in `RESULTS_LEDGER.md` |
| `DEVIATIONS.md` | every change to the plan, with reason and effect (pre-freeze, pre-confirm, post-hoc) |
| `RESULTS_LEDGER.md` | one line per reported number (script, commit, input md5, reviewer) |
| `DATA.md` | data sources, versions, checksums, licences |
| `LESSONS.md` | process rules for running agentic projects (with evidence) |
| `agents/` | reports written by each subagent (methods, checks, numbers) |
| `reviews/` | design reviews, results review, report review, retrospective |
| `src/strucevo/` | Python package: ML distances (GPU), pair pipeline, statistics engines |
| `scripts/` | all pipeline, tuning, endpoint, verification, exploratory and figure scripts |
| `tests/` | pytest suite (distance library, pair pipeline, statistics engine) |
| `config/` | seeds and the frozen tuning values |
| `results/endpoints/` | confirmatory endpoint outputs E1–E6 (+ Holm, interpretation) |
| `results/verify/`, `results/posthoc/` | verification (V2–V11) and post-hoc analyses |
| `results/*`, `figures/*` | exploratory analyses (P-clu, H1e, Goal 2, case studies) |
| `data/processed/` | small committed tables (species panel, PDB control set) |
| `envs/` | conda environment spec (`envs/README.md`) |
| `technical report/` | report and its figures |

Not committed: raw and intermediate data (`data/raw`, `data/interim`; regenerable, several hundred GB) and all **TimeTree-derived ages** (licence forbids redistribution; regenerate with `scripts/panel/timetree_pairwise.py`).

## Reproducing

1. **Environment.** Build with `envs/README.md` (conda + pip torch cu128). Ancestra is pinned to commit ffb0bbe and installed per `agents/ancestra_report.md`.
2. **Data.** Download with `scripts/download_afdb_clusters.sh` and `scripts/download_domain_annotations.sh`, plus the FoldComp index (see `DATA.md`).
3. **Panel and structures.** Build them with `scripts/panel/` (`agents/panel_report.md`).
4. **Orthology.** `conda_envs/ancestra/bin/python scripts/ancestra/run_ancestra.py --proteomes data/interim/ancestra_input --out data/interim/ancestra_euk180 --preset accurate --no-labels`
5. **Tree, families and pairs.** Run `scripts/hogs/panel_tree.py`, `scripts/hogs/build_hogs.py` (seeded; the split md5 is in the ledger), then `scripts/run_after_ancestra.sh`.
6. **Tuning (tune split).** `bash scripts/tuning/run_tuning.sh`
7. **Endpoints.**
   - Dry run: `python scripts/endpoints/run_dryrun.py`
   - Confirmation: `CONFIRM_RUN=1 python scripts/endpoints/run_confirm.py`. It refuses to run unless the plan, tuning and split md5s match the ledger.
8. **Verification.** `pytest tests`, `scripts/verify/`, `scripts/posthoc/`

Seeds are global 20261030 (`config/seeds.yaml`). Hardware used: 48 CPUs, 500 GB RAM, 1× A100 40 GB. Compute: about 24 h within a 48 h cap, across two nodes.
