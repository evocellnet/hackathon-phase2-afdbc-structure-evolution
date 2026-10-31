# Editing notes (editor pass, 2026-10-31)

**Length.** ~15.8k → ~12.2k words (prose words incl. references; `wc -w` 12.7k counts table pipes). Repetition between Methods, Results and Discussion removed; all figures and key numbers kept; each subsection still separates the pre-registered result from the post-hoc interpretation.

**[CHECK] markers (26 → 0).**
- Model mix: PI Claude Opus 5.5; Opus 5.5 subagents for most technical work; Sonnet 5.5 for environment build, AFDB-cluster baseline, annotations, noise floors, PDB preparation, AI-failure-mode review and retrospective (per PI instruction; agent report headers record Opus 5.5 where stated).
- σ_k scale: confirmed as SD of log k (ANALYSIS_PLAN §4.3: Var(log k_h)); post-hoc values checked in results/posthoc/sigma_k_internal.json and agents/posthoc_verification_report.md.
- E4 permuted-age negative control (D13): not in the confirmation output. What exists: tuning dry runs (results/tuning/dryrun*/E4.json) pass 20/20; synthetic V6 (results/verify/V6.json, V6_edges.json): plan permutation 16/20, D13 20/20 under plan cells, 11/20 under D32 cells (leak −0.013 vs −0.077), hence C34. Reported as such in §3.9 and §7, with a note that this control cannot detect the terminal-branch artefact.
- README per-stage commands: claim removed; runtimes stated in "Rerunning".
- All references verified (PubMed/bioRxiv/publisher); missing ones added: LG, Yang 1994 Γ, SIMEX (Cook & Stefanski), Holm, DerSimonian–Laird, AliSim, IQ-TREE 3, MAFFT, FAMSA, AU test, Biopython, BUSCO, Betts 2018 (PMID 30127539), Strassert 2021. US-align dropped (not cited in text).
- Fullmer et al. 2025: title verified; the claim was rephrased to what the abstract supports (3Di adds less information in short, α-helical proteins such as histones), replacing "less clock-like"/"noisier at short times".

**References.** "References to add" merged into one alphabetical list (47 entries); every in-text citation has an entry and vice versa (checked by script).

**§2 rewritten.** Timeline table (node 1 → node 2; ~22 h node time within the 48-h cap; one confirmation run 06:00–07:47), V11 42/42, post-hoc review plus independent verification, final DEVIATIONS.md counts (34 rows: 32 pre-confirmation rows = 34 numbered items, plus P1, P2).

**Number checks.** Spot-checked against results/endpoints/*.json, DEVIATIONS.md, RESEARCH_TREE.md, START.txt and CONFIRM_RUN.json. One addition: E5's excess persists when hidden-paralogy pairs are removed (E5.json), now stated.

**Figures.** Placeholder "Figure 3 (to come)" removed; figures renumbered 1–11 in order of appearance (figA→3, figD→4, figE→5, figC→6, figF→7, figB→8, figG→9, figH→10, figI→11); all paths exist and each figure is cited in the text.

**Terminology.** Internal labels removed (X4–X8, P-clu); HOG, η, ν, k, c, ρ, MDO/LDO, N0 and deviation-ID prefixes defined once; "3Di+AA alignment" used throughout.
