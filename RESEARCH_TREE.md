# Research tree — AFDB-cluster structural divergence & the structural molecular clock

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

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
- [~] 4. DESIGN.md + review rounds (≤3) — round 1 in progress
- [~] 4a. Species panel + proteomes + AFDB v4 structures + TimeTree ages
- [~] 4b. Distance library with toy/simulation tests (src/strucevo)
- [~] 4c. P-clu baseline: in-cluster fident vs TM/LDDT (AFDB clusters file 8)
- [ ] 5. Goal 1 — Is structural divergence linear in time? (Elofsson 2009 revisited at AFDB scale)
- [ ] 6. Goal 2 — Evolution of AFDB clusters (orthologs/paralogs, domains, case studies)
- [ ] 7. Verification & independent review of results
- [ ] 8. Technical report
