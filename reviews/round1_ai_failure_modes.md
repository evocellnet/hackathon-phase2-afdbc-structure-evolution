# Round 1 review: AI-agent failure modes (DESIGN.md v1)

## 1. Verdict
**Major revision, fixable within scope.** The science is well framed and aware of selection bias. The execution risks are: a UniProt/AFDB version mismatch that can silently gut the panel, no frozen primary endpoints (large forking-path space), a verification plan with no executable tests or thresholds, a .gitignore that will silently drop committed results, TimeTree licence exposure, and no concrete resume protocol for the node move.

## 2. Change proposals

1. **[MUST] Resolve the version mismatch (D1).** The panel uses UniProt reference proteomes 2026_03, but structures come from AFDB v4 (UniProt 2021_04). LESSONS #2 says "AFDB is now v6". Accessions and sequences differ, so "drop proteins without a v4 model" removes an unknown, species-dependent fraction. Ancestra then runs on incomplete proteomes, where missing genes look like losses. Build the panel from AFDB-v4 entries per taxId, or write `results/panel/coverage.tsv` (per species: proteome size, with v4 model, sequence-identical, dropped by reason). Assert kept fraction ≥ a preset threshold (e.g. 80%) or exclude the species with a logged reason. Correct LESSONS #2 to say which route serves which version.

2. **[MUST] Frozen analysis plan.** Create `config/analysis_plan.yaml`, git-tagged `prereg-v1` with its sha256 in RESEARCH_TREE, before any pair metrics exist. Fix one primary structural and one primary AA metric for H1a–c, all thresholds (pLDDT 70, ≥50 residues, k=3, bins, d* = 10%), the model set, the rule for "linear" (ΔBIC plus out-of-range error), the FDR level, and what outcome counts as "not a clock". Everything else (other metrics, strata, H2e) is auto-labelled exploratory. Post-freeze changes go in append-only `DEVIATIONS.md`. Currently 6+ metrics × 2 axes × 4 models × strata is a p-hacking garden.

3. **[MUST] Break alignment circularity.** §4 takes AA distance from the same structure-biased Foldseek 3Di+AA alignment. At deep divergence this biases AA identity and couples AA and 3Di mechanically, which can create or erase "linearity". Make a sequence-only alignment (FAMSA/MAFFT global) the primary AA axis; keep the structure-alignment AA distance as a sensitivity analysis.

4. **[MUST] Parser/column contract tests** (`src/schemas.py`, `tests/test_parsers.py`) on a 1,000-row real slice.
   - AFDB file 8: query = rep, target = member; IDs end `-model_v3.cif`. Strip the suffix and map to v4 accessions, with join rate ≥99% asserted.
   - Foldseek: always pass an explicit `--format-output`. `alntmscore` is not "shorter-protein" or "aligned-length" normalised. If both normalisations are wanted, compute them with US-align and cross-check on 50 pairs (|ΔTM| < 0.02). Test `--alignment-type` and the lDDT column on a known pair.
   - TimeTree: the plan picks one field (`precomputed_age` vs `adjusted_age`). Assert unit = Mya and 0 < t < 4600. Write unmatched species to `timetree_unmatched.tsv` instead of dropping them silently.
   - OrthoXML: parse only through Ancestra's reader/`hog_pair_labels.py`; no ad-hoc XML parsing.

5. **[MUST] TimeTree licence guard.** Redistribution of data *and transformations* is prohibited, and `.gitignore` does not cover `data/processed/` or results with ages. Keep age-bearing tables in gitignored `results/restricted/`. Commit only scripts and aggregate fit parameters. Add a pre-commit `scripts/check_commit.sh` that fails on tracked files matching `*timetree*` or with `age|mya` columns.

6. **[MUST] Results ledger and gitignore fix.** `.gitignore` ignores `*.parquet`, `*.gz` and `*.npz`, so "committed results" can vanish. Commit small tables as `.tsv`. Add `results/LEDGER.tsv` (id, hypothesis, script@commit, input md5s, seeds, output path, status, reviewer), and `scripts/ledger_check.py` asserting each output is tracked or listed as regenerable. Every number in DESIGN or the report cites a ledger id. Scout claims (e.g. "precision ~0.99", "17 queries/s") enter DESIGN only with a log citation. Commit at each milestone and before the node move.

7. **[MUST] Checkpoint/resume protocol.**
   - `STATE.md` is rewritten at each milestone with the running commands, expected outputs and resume steps.
   - Long jobs run via `scripts/run_stage.sh <stage>`, which is idempotent and writes `.done` with md5s.
   - No job starts after ~20:00 on node 1 unless it can finish or checkpoint by 20:30. Confirm Ancestra resumes *within* its DIAMOND stage, or split the run by clade into chunks under 1 h.
   - `scripts/resume_check.sh` verifies envs, GPU and ledger md5s, then reruns only missing stages.
   - Fix `elapsed.py`: total time includes the inter-node gap and node 2's budget is not reported. Write the second START line as the first action on node 2.
   - The panel fetch (2–3 h) followed by Ancestra (4–8 h) leaves almost no slack on node 1. Run a 60-species Ancestra dry run first and size the panel from the measured cost.

8. **[MUST] Make §6 executable.** Create `tests/` plus `scripts/verify.sh` (non-zero exit on failure, log committed) with these tests:
   - **AliSim recovery:** simulate under a different model (WAG+Γ with indels) and under LG, because LG-on-LG is circular. Report bias for the mismatched case. Compare against IQ-TREE pairwise distances (|Δ| < 0.02).
   - **3Di model fit:** fitted on held-out HOGs, not the analysed pairs.
   - **Model selection:** recovers the true model ≥90% over 200 synthetic replicates, including a mixed case with power at d*.
   - **Controls:** identity pair = 0; shuffled and non-homologous pairs define floors, and the plan states how floors enter the −ln transforms.
   - **Smoke test:** a 5-species, 200-HOG slice with human–mouse identity ≈ 85% and pair counts equal to `hogs_to_pairs`.

9. **[SHOULD] Independent-review protocol.** The reviewer gets code and raw outputs, not the producer's summary. It recomputes ≥20 random pairs by an independent route (US-align vs Foldseek; IQ-TREE vs GPU ML) and re-derives each headline number from the ledger. It files `reviews/result_<id>.md` as PASS/FAIL and may not edit results. Nothing is marked [!] in RESEARCH_TREE before PASS.

10. **[SHOULD] Species-tree/time mapping.** Ancestra's guide tree (NCBI taxonomy, with polytomies) differs from TimeTree's topology, so "(HOG, MRCA node)" must not be matched by node name. Query TimeTree per species pair and record `n_studies`. Flag pairs whose taxonomic LCA differs from the TimeTree LCA and exclude them from H1b with a count. Fill deep nodes missing from the API only from a cited file, never from agent memory.

11. **[SHOULD] Silent-default audit.** All thresholds, seeds and paths come from `config/*.yaml`. Scripts echo the resolved config and fail on missing keys. Shell scripts use `set -euo pipefail`. Every dropped row is counted in `drop_log.tsv` with a reason code, and the dropped fractions are reported.

12. **[SHOULD] H1d/H1e guards.**
    - H1e: pre-specify matching on chain state, ligand and resolution. Check the SIFTS mapping against sequence identity ≥95%.
    - H1d: swapping ingroup labels is not a null for relative-rate tests (it is symmetric). Calibrate the false-positive rate with AliSim on ultrametric trees.

13. **[COULD] Literature/overclaim check.** A fresh agent verifies cited PMIDs/DOIs and marks sources not read in full. The report says "predicted-structure evolution" and lists negative and inconclusive hypotheses.
