# Research tree — AFDB-cluster structural divergence & the structural molecular clock

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (latest first)

### M3: 2026-10-30 ~14:00 (t = 4.6 h; node-1 stop target 20:30)
**Orthology.** Ancestra has been running on the 180-eukaryote panel since 11:18. The families stage finished at 12:52. The scores stage is at about 13/34 DIAMOND blocks, ETA ~15:25. The merge stage follows, so the whole run should end around 17:30–18:30. Go/no-go checks are at 15:00, 17:00 and 18:30. If the merge stage has not finished by 19:30, it restarts first on node 2.

**Ready to run on Ancestra's output.** The infrastructure is built and verified:
- species panel, structures and the dated, collapsed panel tree;
- HOG builder (tested on euk20);
- pair-metrics pipeline (V5/V8/V9 pass);
- statistical engines (V2–V4 and V6 pass after the logged pre-freeze fixes);
- noise floors;
- per-protein cluster and domain annotations;
- PDB control set.

Waiting for Ancestra: HOG building, then pair metrics (≤ 3.5M pairs per alignment, ~4–5 h, mostly node 2), then tuning on 20% of HOGs, freeze, and a single confirmation run.

**Corrections to M2:**
- The **primary alignment is Foldseek 3Di+AA**, not TM-align (D22). TM-align fails the identical-sequence invariant and mis-registers close homologs.
- The E5 neutral-clock acceptance is practically unreachable at the pre-set margins; E5 can still reject. The margins are kept unchanged.

**Decision notes M3:**
- All agent "decisions for you" were resolved immediately and logged in DEVIATIONS.md (D1–D22, C9–C28). Every one was made before any endpoint was computed.
- The pair budget is fixed at ≤ 3.5M pairs per alignment (primary + E6).
- Node 1 runs only light or GPU work until Ancestra finishes, to avoid CPU contention.
- LESSONS.md was pruned to 12 rules after a fresh-agent retrospective.

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
- [x] 4e. Statistical engines E1–E6 + synthetic verification (V2 cov 0.96 with D2; V3 false-concave 0/100, min detectable |c−1|≈0.10; V4 recovers truth; V6 20/20 with D13). [!] E5 N0 acceptance practically unreachable (β CI ±0.24); margins kept, E5 likely inconclusive unless it rejects
- [x] 4f. PDB/SIFTS control set (H1e): human 3502 qualifying chains (2524 post-AF2); ~513 human–mouse pairs with PDB on both sides (304 post-AF2); deep pairs ~100–300
- [x] 4a. Species panel (agents/panel_report.md): 180 eukaryotes (Metazoa 66, Fungi 40, plants 20, protists across all supergroups) + 100 prokaryotes (60 bacteria / 40 archaea, Asgard = MK-D1 only). Gate: ≥80% exact-sequence v4 coverage. 3.29M + 0.28M structures fetched, 0 errors; 94.9% exact-sequence match. TimeTree: 15,755/16,110 pairs dated, 138 distinct ages (0.4–1,516 Mya); 52% of pairs at the root age → analyse by node/edge, not pair (already in design). Missing: rat/zebrafish/chicken/Xenopus (low v4 coverage). Polyploids included (cotton×4, soybean, Brassica, goldfish), handled by the 1:1-copy rule; limitation noted.
- [~] 4g. Ancestra on 180 eukaryotes (accurate, launched 11:18). Families stage done 12:52 (94 min, 37 GB; single-threaded Python passes over 139M DIAMOND hits dominate)
- [x] 4h. Pair-metrics pipeline (agents/pairs_pipeline_report.md): 3 alignments, V5/V8/V9 pass, 56k demo pairs in 4.5 min. [!] TM-align mis-registers close homologs (human–mouse core p 0.091 vs 0.064) and fails V5 → D22: primary alignment = Foldseek 3Di+AA; tmalign/seqonly as brackets
- [x] 4i. Panel tree (180 leaves, 100 internal nodes after collapsing; human–mouse 89.5 Mya) + HOG builder tested on euk20 (projection euk180: 1.0–1.8M primary pairs)
- [x] 4j. Noise floors: 100k random pairs (TM 0.22–0.33 tmalign) + 17.8k identical-sequence pairs. [!] AF prediction noise: identical sequences median TM 0.95, η_3Di = 0.053 (rising with length), η(1−LDDT) = 0.02
- [x] 4k. Goal 2-ii per-protein annotation: 3.74M proteins; 80% (euk) in non-singleton AFDB clusters; TED domain in 89%; classes α 19.6% β 9.8% α/β 24.0%
- [~] 4l. Goal 2-iii cytoskeleton (agents/casestudy_report.md): 3,033 actin-fold + 2,623 tubulin-fold proteins; 400-seq trees (UFBoot ≥95 at ~55% of nodes). First look (exploratory): actin Euk+Asgard stem shows no structural burst (3Di/AA ratio relative ρ=0.33, 11th pct); tubulin stem alignment-dependent (ρ 1.21 MAFFT vs 0.60 FoldMason; AA branch near saturation). Caveats: Asgard = MK-D1 only (6 actins, 0 tubulins); Asgard actins nest in Arps and artubulins/BtubA/B in eukaryotic tubulins (likely LBA). TODO node 2: branch engine on fixed topology, constrained-topology sensitivity.
- [x] 4b. Distance library with toy/simulation tests (src/strucevo): 33/33 pass
- [x] 4c. P-clu baseline (agents/pclu_report.md): 27.7M member→rep pairs. Median TM 0.39 at 10% id → plateau ≈0.83 from ≥40% id; LDDT rises smoothly 0.52→0.94. Plateau tracks rep pLDDT/length (at id ≥0.6, TM 0.58 for pLDDT<60 vs 0.97 for ≥90), i.e. a prediction-quality ceiling, not saturation. Older (cross-superkingdom) clusters have higher TM at equal identity [exploratory; confounded with pLDDT → check in X7]. Artefact below 7% id (short alignments).
- [ ] 5. Goal 1 — Is structural divergence linear in time? (Elofsson 2009 revisited at AFDB scale)
- [ ] 6. Goal 2 — Evolution of AFDB clusters (orthologs/paralogs, domains, case studies)
- [ ] 7. Verification & independent review of results
- [ ] 8. Technical report
