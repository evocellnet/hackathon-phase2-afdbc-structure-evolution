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

---

# Changelog: response to independent review (reviews/report_review.md, 2026-10-31)

Review saved verbatim to `reviews/report_review.md`. All Must, Should, Cut and Add items applied. Disclaimer block unchanged.

**Must**
1. Abstract: per-protein offset stated as an observation ("behaves like AlphaFold model noise; a fast-saturating real structural component cannot be excluded"); §4 main point, §6 main point and §7 softened the same way.
2. Curvature "0.73–0.98" → "bracket between about 0.75 and 1.0" (Abstract, §6.1(1), §7), as in posthoc_verification_report.
3. E5: Abstract "no lineage-specific excess is detected (a test with 24% power)"; §4.4 rewritten: σ_L departs from N0, slope 0.57 [0.10, 1.05] neither meets nor rejects the margin (p_reject(β<0.8) = 0.131, E5.json); single-copy guard uninformative (slope CI [−25.9, 26.6], excess variance 0.095 vs 0.098 full data); "artefact suspect" rests mainly on the internal-branch guard.
4. §1.4: definition of "clock" added; Abstract now "changes in proportion to sequence (about one 3Di substitution per 4.6 AA substitutions on internal branches)"; §6.1(3) wording aligned; §4.3 notes AA grows sublinearly with time (slope 0.85, key_numbers.json).
5. Discussion main point: "(post hoc)" added.
6. §5.2 corrected (checked in agents/goal2_cluster_audit_report.md): split model (279,196 pairs, 48,130 HOGs) and paralog–ortholog architecture contrast (212,502 pairs) used all HOGs; only the against-time analyses (annotation floor, architecture vs age, cluster separation vs age) used the 178,593 tuning-set one-to-one pairs.
7. Old Figure 6 (now 7) reduced to panel (a); panel (b) (AA/(3Di − η)) dropped; caption and text softened to "consistent with extra per-model noise in young comparisons (§4.3)".

**Should**
8. Fold-class caveat (all branches, not recomputed on internal branches) in Abstract bullet, §4.2, §6.1(2), §7.
9. E6b caveat added to §6.1; Abstract names Warwick Vesztrocy et al. (2025), "direction only".
10. Case-study statistics unified on Spearman, recomputed on the plotted subset (internal branches, positive lengths under both methods; results/casestudy/branch_ratios_*_ml.tsv): AA 0.77 (actin) / 0.83 (tubulin), ρ 0.16 / 0.32. The case-study report's 0.74–0.76 / 0.13 / 0.29 could not be reproduced exactly (subset not documented); other subsets tried give 0.68–0.83 and 0.14–0.38, stated in the text. Figure 12 regenerated (scripts/casestudy/19_figures_node2.py `method_fig`): 2×3 layout, plain labels, Spearman in titles, unreadable column 4 dropped.
11. Figure 3 caption: plotted data are tuning + confirmation (6,067 / 6,114 families) vs 4,850 in the confirmatory fit; TM-align 3Di medians bend upward beyond AA distance 1.5 (figA_shape_binned.tsv).
12. Figure 4 regenerated (figures_main.py figD, same seed; key_numbers.json byte-identical in content): "HOG" → "family", panel c "By age of sampled clade" (= sampling stratum), panel d "branches"; caption states all-data n (α 1,661, β 867) vs confirmation-only text n (1,336, 696).
13. Figure 8 (AF vs PDB) redrawn with new script scripts/h1e/fig_report_x4.py (same pair construction and mismatch filter as analyse_h1e.py; n = 197 post, 484–485 pre): main figure = 3Di and 1 − global LDDT, enlarged; all five measures moved to new Figure S1. "tdiaa" removed.
14. Independence sentence added to §2.
15. Illergård 2009 read via abstract + 2010 thesis (LITERATURE.md) stated in §1.1; Chothia & Lesk substitution-axis rewrite marked as our own derivation (LITERATURE.md [Ours]).
16. §5.2: possible contribution of Ancestra orthology errors (fusions, fragments) to length-driven splitting.
17. Internal IDs (D22, D31, D32, C35, C36, C42, P1, V11) removed from Results/Discussion/Limitations; kept in §2–§3 where defined. "Cell" defined in §3.8, "censored" in §4.3 (ANALYSIS_PLAN.md §2.6).
18. Deviations: "34 entries (one, C19–C21, covers three items): 32 before, 2 after". reviews/ now holds 13 files: 12 earlier reviews + this report review (stated as such).
19. E6a: "Holm-adjusted p reported as 0, i.e. below numerical precision" (E6.json p_holm = 0.0).

**Cut**
- §2 timeline table → three bullet lines plus total compute.
- Old Figure 6 → panel (a) only (scripts/exploratory/figures_main.py figC now single-panel).
- §3.9 "Reproducibility" merged into Data availability "Reproducing the analysis".
- §3.5 D22 detail → one sentence.

**Add**
- New Figure 5: tip-artefact schematic (`technical report/figures/make_figT_tip_artefact.py`, illustrative k = 4.6, δ = 0.014; no data). Figures renumbered: old 5–11 → 6–12; new Figure S1.
- Plain-language reading of k ≈ 4.6 in §4.2.
- Practical sentence in Abstract (inflated 3Di terminal branches; AFDB clusters not ortholog groups).
- §6.1(1): experiments that would settle the open question (same-protein PDB–PDB floor; AF models of identical sequences with varied MSA depth).

**Files changed:** report.md; figures figD_ratio_k.png, figC_time.png, figF_af_vs_pdb.png, figI_cytoskeleton_methods.png (regenerated), figT_tip_artefact.png, figS1_af_vs_pdb_all.png (new); scripts/exploratory/figures_main.py, scripts/casestudy/19_figures_node2.py (labels/layout only), scripts/h1e/fig_report_x4.py (new). Plotted data unchanged. Not committed.
