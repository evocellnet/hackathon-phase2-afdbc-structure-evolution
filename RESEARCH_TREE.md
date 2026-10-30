# Research tree — AFDB-cluster structural divergence & the structural molecular clock

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### M2 — 2026-10-30 ~11:10 (t = 1.8 h)
The design passed three review rounds: 5 reviewers in round 1, 3 in round 2, and a round-3 consistency check. DESIGN.md v3 is final, and ANALYSIS_PLAN.md is the executable, frozen specification; only its tuning values remain open. The core design is two engines.
- Within-family pairwise distances test the response shape (does structure saturate relative to sequence?).
- Per-family branch lengths on the fixed species tree estimate the 3Di/AA conservation ratio, its trend with node age, and lineage effects (the neutral-structural-clock null N0).
- Paralogs and the more-diverged copy after duplication (MDO) serve as the positive control.

Statistical safeguards: tuning on 20% of families and confirmation once on 80%, Holm correction over five endpoints, pre-registered interpretation rows, a spike-in power check, and a simulation of false-curvature rates.

The distance library is built and verified (33 tests; ML distances unbiased within 1.6%; 1.3 M pairs/s on GPU). Running in parallel: species panel construction, the pair-metrics pipeline (demo on the 20-proteome Ancestra run), the statistical engines (synthetic checks V2–V4), the PDB control set, and the P-clu baseline. Next: launch Ancestra on the eukaryote panel as soon as it is final (target ≤ 13:00), then the node-move checkpoint by 20:30.

**Decision notes M2:**
- Primary alignment is coordinate-only TM-align, so neither the amino-acid nor the 3Di alphabet is favoured; the sequence-only and 3Di+AA alignments are sensitivity analyses.
- AIC/BIC were dropped for equivalence-type decision rules with pre-set margins.
- Goal 2 was cut to three items: duplication as positive control, cluster audit with domain architecture, and a cytoskeleton case study (time-boxed).
- The reviewer-written ANALYSIS_PLAN.md was accepted wholesale. Any later simplification will be logged in DEVIATIONS.md, not edited silently.

### M1 — 2026-10-30 ~10:45 (t = 1.4 h)
Scouting is complete. The literature review (LITERATURE.md) changed the framing in three ways. Illergård et al. 2009 measured "time" as substitutions per site, so testing linearity against geological time is a new question. Pairs taken from inside AFDB clusters are truncated by the structural clustering filter and can create artificial saturation, so the linearity test must use orthologs chosen without looking at structure. TimeTree ages are usable mainly for eukaryotes. The data are in place: AFDB clusters v3 and v6, the FoldComp v4 index plus all 2.3 M cluster representatives, TED, Pfam, NCBI taxonomy and TimeTree (API plus the 2015 tree). Ancestra installs and passes its tests. Its ortholog calls are reliable (148/149 Ensembl human–mouse 1:1 recovered), but its duplication timing is not (36% agreement).

DESIGN.md v1 is drafted and in review round 1 with 5 reviewers: external PI, molecular evolution, statistics, algorithms/compute and AI failure modes. Three design-independent tasks run in parallel: building the species panel (~300 proteomes plus structures plus TimeTree ages), a tested distance library (ML distances on GPU, AliSim checks), and the in-cluster (P-clu) baseline from AFDB-cluster alignments.

**Decision notes M1:**
- D1: use the AFDB v4 model sequence as the sequence everywhere.
- D2: never rely on Ancestra duplication ages; match paralogs on amino-acid distance instead.
- D3: run Ancestra once on whole proteomes, not per cluster (cost).
- Approved the 25 GB FoldComp downloads. Deferred afdb50 (114 GB) and the 95 GB metadata table as not needed.
- Panel building started before the design review closes, because species choice is robust to likely design changes.

### M0 — 2026-10-30 09:30 (t = 0.1 h)
Project launched on node 1 (12 h node budget, vacate by 20:52; total 48 h). Workspace was empty; network access confirmed to Foldseek cluster site, AFDB-cluster S3 bucket, EBI FTP and GitHub. No pre-installed bioinformatics tools; a conda env is being built.

Four scouting agents run in parallel: literature (LITERATURE.md), data acquisition (AFDB cluster tables, structures, taxonomy, divergence times, domains), Ancestra install/characterisation, compute environment. Next milestone: draft DESIGN.md from their reports, then up to three rounds of independent review.

**Decision note M0:** Scout before design. Divergence times (TimeTree) and Ancestra feasibility are the two critical unknowns for Goal 1 (linearity in time) and Goal 2 (ortholog/paralog). Design waits for both.

## Tree
- [x] 0. Setup: env (conda_envs/main; foldseek GPU works), record files, git
- [x] 1. Literature & datasets (LITERATURE.md)
- [x] 2. Data acquisition scouting (DATA.md; agents/data_scout_report.md)
- [x] 3. Ancestra feasibility (agents/ancestra_report.md): works; orthologs reliable, duplication timing not
- [x] 4. DESIGN.md v3 + 3 review rounds; ANALYSIS_PLAN.md frozen pending tuning
- [~] 4d. Pair-metrics pipeline (src/strucevo/pairs)
- [~] 4e. Statistical engines E1–E6 + synthetic verification V2–V4
- [~] 4f. PDB/SIFTS control set (H1e)
- [~] 4a. Species panel + proteomes + AFDB v4 structures + TimeTree ages
- [x] 4b. Distance library with toy/simulation tests (src/strucevo): 33/33 pass
- [x] 4c. P-clu baseline (agents/pclu_report.md): 27.7M member→rep pairs. Median TM 0.39 at 10% id → plateau ≈0.83 from ≥40% id; LDDT rises smoothly 0.52→0.94. Plateau tracks rep pLDDT/length (at id ≥0.6, TM 0.58 for pLDDT<60 vs 0.97 for ≥90), i.e. a prediction-quality ceiling, not saturation. Older (cross-superkingdom) clusters have higher TM at equal identity [exploratory; confounded with pLDDT → check in X7]. Artefact below 7% id (short alignments).
- [ ] 5. Goal 1 — Is structural divergence linear in time? (Elofsson 2009 revisited at AFDB scale)
- [ ] 6. Goal 2 — Evolution of AFDB clusters (orthologs/paralogs, domains, case studies)
- [ ] 7. Verification & independent review of results
- [ ] 8. Technical report
