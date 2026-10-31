# Review of "technical report/report.md" (independent reviewer, 2026-10-31)

## Verdict
**Minor-to-moderate revision.** Honest, unusually transparent, mostly accurate; pre-registered decisions reported unchanged, post-hoc labelled in Results. 60+ numbers spot-checked against results/endpoints/*.json, results/posthoc/, posthoc_verification_report, results_review_statistics, DEVIATIONS, X4 and cluster-audit reports — nearly all match. Problems: (1) Abstract and Discussion state the "AlphaFold noise" mechanism as fact, which §4.1 says the data cannot settle; (2) E5 is overstated in two places; (3) "clock" is used loosely — only proportionality to sequence was tested, not time-constancy; (4) a few factual/caption mismatches. Accessibility good at section level; internal IDs/symbols still make the text heavy.

## Must
1. Abstract cause stated as fact. Current: "Each AlphaFold model carries protein-specific noise (about 0.01–0.017 …) that the pre-registered noise floor missed." Replace with: "Non-identical orthologs carry an extra per-protein 3Di offset (about 0.01–0.017 substitutions per site, larger at low model confidence) that the pre-registered noise floor did not include. It behaves like AlphaFold model noise, but a fast-saturating real structural component cannot be excluded."
2. Abstract curvature range "0.73–0.98" comes from naive (non-SIMEX) fits; one window gave 0.57; verifier concluded "bracket … roughly 0.75–1.0". Replace with "Post-hoc estimators that allow for it bracket the curvature between about 0.75 and 1.0." Same in §6.1(1) and §7.
3. E5 wording: Abstract "the lineage effects disappear" → "no lineage-specific excess is detected (a test with 24% power)". §4.4 "both statistics depart from N0": slope 98% CI [0.10,1.05] overlaps the margin (p_reject(β<0.8)=0.13) → "the excess lineage SD departs from N0 (σ_L = 0.31; N0 needs < 0.10); the slope, 0.57 [0.10, 1.05], is too imprecise either to meet or to reject the N0 margin." §4.4 "it also disappears in single-copy families": CI [−25.9, 26.6], point excess 0.095 = full-data → "in single-copy families the test is uninformative"; note "artefact suspect" rests mainly on the internal-branch guard.
4. What "clock" means. Add to §1.4: "Throughout, 'clock' means that the 3Di/AA ratio is constant across branch ages and lineages. Whether AA change itself is clock-like in time was not tested; a structural 'rescaled clock' inherits the sequence clock's irregularities." Abstract: "behaves like a rescaled sequence clock" → "changes in proportion to sequence (about one 3Di substitution per 4.6 AA substitutions on internal branches)". (AA distance grows sublinearly with time: Fig 6 slope 0.85.)
5. Discussion main point: add "(post hoc)" to "on internal branches its rate relative to sequence is roughly constant, at about one-fifth".
6. §5.2 "pair-level analyses on tuning-set families only" is wrong: split model (279,196 pairs) and paralog-architecture contrast (212,502 pairs) used all families; only architecture-vs-time is tuning-only. Correct and note.
7. Figure 6 caption: panel (b) plots AA/(3Di − η), not "AA/3Di"; "mostly because of young, noisy comparisons" not shown — redraw after subtracting estimated per-protein excess or soften to "consistent with extra per-model noise in young comparisons (§4.3)".

## Should
8. Fold-class contrast (α 4.7× vs β 2.6×) includes terminal branches; not recomputed on internal branches. Add caveat (part could be tip artefact if α models are more confident); caveat or drop in Abstract.
9. E6b in Discussion: after "The duplication result (E6b) points their way" add "although its size equals the per-model offset and was not adjusted for pLDDT differences between copies". Abstract: name Warwick Vesztrocy et al. (2025) and say "direction only".
10. Case-study statistics mixed: text "Spearman 0.74–0.80" vs report 0.74–0.76; text ρ agreement 0.13/0.29 vs Fig 11 caption log-r 0.04–0.08 — use one statistic. Fig 11: fix clipped titles and internal labels (b_AA, r_log, tdiaa); drop or redraw unreadable column 4.
11. Figure 3 caption: plotted data are tuning+confirmation (6,067/6,114 families) unlike confirmatory fit (4,850); TM-align 3Di medians also bend upward beyond 1.5.
12. Figure 4: panel n's (α 1,661, β 867) differ from text 1,336/696 — state all data; panel c "HOG root age" vs text "age of the sampled clade"; panel d "edges" vs text "branches".
13. Figure 7: five-panel strip too small, shows "tdiaa"; enlarge 3Di and 1−LDDT panels, move the rest to a supplement.
14. Independence of checks: add "Independence means separate code and context, not a different model or a human; shared blind spots (like the noise-floor assumption, missed by every agent) remain possible."
15. Illergård credit: state the closed-access 2009 paper was read only via abstract and the 2010 thesis summary (LITERATURE.md); mark the Chothia & Lesk "rewritten on a substitution axis" argument as the project's own derivation.
16. Cluster audit: note length-driven splitting may partly reflect errors in unverified Ancestra orthology (fusions, fragments).
17. Internal IDs: in Results/Discussion replace D22/D30–D32, C25/C29/C34–C36/C42, P1/P2, V11 with plain words; define "cell" and "censored" at first use.
18. Deviation count: "34 entries (one covers three items), 32 made before the confirmation run and 2 after it"; reviews/ holds 12 reviews, not 13 (recount after saving this review).
19. E6a p-value: E6.json reports p = 0 (below numerical precision) — write that instead of "Holm p < 10⁻¹⁵".

## Cut
- §2 timeline table → three lines plus total compute.
- Figure 6 largely redundant with Fig 5 and least defensible: move to supplement or keep only panel (a) with fixed wording.
- "Reproducibility" and "Rerunning" duplicate each other: keep one.
- D22 detail in §3.5 → one sentence.

## Add
- A tip-artefact schematic: terminal branches with a fixed 3Di offset, showing why this lowers k on young branches and pushes the age slope negative.
- Plain-language reading of k ≈ 4.6: for every 4–5 amino-acid substitutions per site, about one residue's local 3D-contact letter changes.
- One practical sentence in the Abstract: 3Di terminal branches from AlphaFold models are inflated, and AFDB clusters are not ortholog groups.
- What would settle the open question: a same-protein PDB–PDB floor, or AlphaFold models of identical sequences built with varied MSA depth.
