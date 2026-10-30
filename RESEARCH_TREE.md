# Research tree — AFDB-cluster structural divergence & the structural molecular clock

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### M0 — 2026-10-30 09:30 (t = 0.1 h)
Project launched on node 1 (12 h node budget, vacate by 20:52; total 48 h). Workspace was empty; network access confirmed to Foldseek cluster site, AFDB-cluster S3 bucket, EBI FTP and GitHub. No pre-installed bioinformatics tools; a conda env is being built.

Four scouting agents run in parallel: literature (LITERATURE.md), data acquisition (AFDB cluster tables, structures, taxonomy, divergence times, domains), Ancestra install/characterisation, compute environment. Next milestone: draft DESIGN.md from their reports, then up to three rounds of independent review.

**Decision note M0:** Scout before design. Divergence times (TimeTree) and Ancestra feasibility are the two critical unknowns for Goal 1 (linearity in time) and Goal 2 (ortholog/paralog). Design waits for both.

## Tree
- [~] 0. Setup: env, record files, git
- [~] 1. Literature & datasets (LITERATURE.md)
- [~] 2. Data acquisition scouting (DATA.md)
- [~] 3. Ancestra feasibility
- [ ] 4. DESIGN.md + review rounds (≤3)
- [ ] 5. Goal 1 — Is structural divergence linear in time? (Elofsson 2009 revisited at AFDB scale)
- [ ] 6. Goal 2 — Evolution of AFDB clusters (orthologs/paralogs, domains, case studies)
- [ ] 7. Verification & independent review of results
- [ ] 8. Technical report
