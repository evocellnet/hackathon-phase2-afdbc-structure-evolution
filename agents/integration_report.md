# Integration report: post-Ancestra pipeline, tuning, endpoints, exploratory (integration engineer, 2026-10-30, 15:20–16:30)

**Bottom line.**
- Every stage was built and run end to end on euk20 test roots: orchestration (a–d) → tuning (15 steps) → §Tuning table filler → tune dry run → guarded confirm run → interpretation → X1/X2/X7.
- Nothing touched the production `data/interim/{split.tsv,hogs.parquet,pairs/}`, `results/endpoints/` or `ANALYSIS_PLAN.md`. Nothing was committed.
- `scripts/run_after_ancestra.sh` is ready for euk180. It adopts the PI's hand-run build_hogs output and never starts build_hogs itself when `ALLOW_BUILD=0` is set.
- **Three findings need a PI decision before the freeze (§4):**
  1. α_3Di at the grid boundary.
  2. E1 biased towards "concave" under a linear truth.
  3. E4 biased under no age effect.

Built in parallel by five sub-agents (H1e; tuning α/SE; tuning misfit/power/V3/yaml; endpoints; exploratory) against one contract: every script takes `--root` and resolves paths in a mini-project layout.

## 1. Commands (production; env `conda activate conda_envs/main`)

```bash
# pair metrics + H1e (after the PI's build_hogs finished; adopts its outputs, never rebuilds)
ALLOW_BUILD=0 THREADS=48 scripts/run_after_ancestra.sh          # steps a b1 b2 c0 c1 c2 c3 c4 d, resumable
scripts/run_after_ancestra.sh status                              # progress + ETA per alignment/split
# tuning (split == tune only); V3 at plan B=200 takes ~3-3.5 h, B=100 ~1.7 h
PAIRS_THREADS=48 PAIRS_RUN_THREADS=48 V3_ARGS="--B 100" bash scripts/tuning/run_tuning.sh
python scripts/tuning/fill_tuning_table.py --plan ANALYSIS_PLAN.md --tuning config/tuning.yaml   # only the §Tuning table
python scripts/endpoints/run_dryrun.py                            # E1-E6 dry run on tune -> results/tuning/dryrun/
python scripts/exploratory/{x1_curves,x2_age,x7_pclu}.py --aln tdiaa   # (+ seqonly, tmalign) tune only
# freeze: commit; append the ledger rows printed by
python scripts/endpoints/run_confirm.py --print-freeze-rows
CONFIRM_RUN=1 python scripts/endpoints/run_confirm.py && python scripts/endpoints/interpret.py --dir results/endpoints
```

## 2. Orchestration: `scripts/run_after_ancestra.sh`

Every step has a sentinel in `$ROOT/logs/pipeline/<step>.done`, and each pair shard has its own `.done` (md5-checked).

| Step | What it does |
|---|---|
| a | Checks the OrthoXML. **No all-pairs label table is built.** build_hogs labels pairs on demand (HogIndex, proven identical to `label_pairs`), and H1e labels only its own candidates. |
| b1, b2 | build_hogs, projection run then real run, with C25 caps 5/5. The projection prints counts and an ETA for every pair step. If a complete build exists it is adopted; a running build_hogs is waited for. |
| c0 | `scripts/pipeline/prep_pair_lists.py` writes per-split lists and the accession list, then fetches FoldComp structures. |
| c1–c4 | Tune: seqonly+tdiaa, then tmalign. Confirm: the same. ~15k pairs per shard. Tune shards are 0…, confirm shards 1000… (new `run_pairs.py --shard-offset`), all in `data/interim/pairs/{aln}/`. Each step prints an ETA from the measured rate (tune measures confirm), then the measured rates. |
| d | `scripts/h1e/run_h1e.py` |

- **Order deviation (minor).** The PI asked for tdiaa → seqonly → tmalign. The tdiaa table joins `aa_ml_seq` from the seqonly shard, so seqonly (the cheapest stage) runs immediately before tdiaa within each shard. tmalign runs last.
- **Tested on fresh test root `data/interim/itest_euk20_orch`** (8 threads):

  | Step | Time |
  |---|---|
  | a | 0 s |
  | b1 | 21 s |
  | b2 | 21 s |
  | c0 | 17 s |
  | c1 | 30 s |
  | c2 | 68 s |
  | c3 | 74 s |
  | c4 | 295 s |
  | d | 146 s |

- Rerun behaviour: everything skips on rerun, and adoption of an existing build was tested.
- Measured rates per 8 threads (34,925 pairs): seqonly 650–700, tdiaa 800–1,050, tmalign 87–100 pairs/s.

**Projection at 48 threads.** Linear scaling from 8 threads gives a lower bound; the demo rates at 48 threads under load (1,400 / 1,140 / 322 pairs/s) give an upper bound. Rescale with the projection pair total; halve the throughput at 24 threads.

| Pairs | seqonly | tdiaa | tmalign | Total |
|---|---|---|---|---|
| 2.5M | 11–30 min | 8–37 min | 1.2–2.2 h | **≈ 1.6–3.2 h** |

Fetch and store build take a few minutes.

## 3. H1e (`scripts/h1e/run_h1e.py`, `src/strucevo/pairs/store_files.py`)

- **Chain rule.** The plan's rule is applied to the raw SIFTS/RCSB tables: X-ray, ≤ 3 Å, released after 2018-04-30, coverage ≥ 0.70. The precomputed `qualifies` column was incomplete.
- **Labels.** Only candidate pairs are labelled, with Ancestra 1:1 orthologs, via HogIndex.
- **Store.** One combined AF + PDB store. PDB residues get pseudo-pLDDT = 100.
- **Pair sets:** `h1e_af`, `h1e_pdb`, and a new `h1e_afpdb` floor set (contract addition). Output goes to `data/interim/pairs_h1e/`.
- **euk20 test:** 160 pairs (299 proteins, cap 300), all three alignments, 0 fail_align, 2.4 min including downloads.
- **euk180 estimate:** about 350 proteins and ≤ 830 pairs, 2–5 min.
- **Open issues:**
  - PDB–PDB pairs use all aligned columns, while AF–AF pairs use the pLDDT ≥ 70 core. A matched-core comparison is needed for X4.
  - Unobserved PDB residues make 3Di gaps, which likely inflate PDB d_3Di (median PDB−AF Δd_3Di = +0.15).

## 4. Tuning (`scripts/tuning/`, `run_tuning.sh` with per-step sentinels; fragments in `results/tuning/*.json` → `write_tuning_yaml.py` → `config/tuning.yaml`)

**Steps:**
- profile_alpha validation (AliSim + `simulate_counts`)
- tune_alpha
- floors_tuned: η and S0 at the tuned α; C28
- select_bootsub, plus a boot_sub run with the new `run_pairs.py --dump-columns` (per-column AA/3Di/pLDDT/lDDT; outputs byte-identical when off)
- se_calibration: 20-residue block column bootstrap → κ(bin), ρ(bin), lddt_D SE model
- apply_tuning: recompute in place at the tuned α. **`aa_ml_seq` stays at α = 1 (C26)**; `recompute_ml.py` was fixed accordingly. It also fills lddt_D/lddt_D_se and can write α×1.5 / α÷1.5 variants to `data/interim/pairs_variants/`.
- tau_misfit, b_e4_timing, e5_power, v3_tuned, sample_throughput

`fill_tuning_table.py` changes only the 11 §Tuning rows (diff-verified, idempotent).

**Full chain on the test root:** 481 s, 0 nulls. Main values (euk20, all stratum 4, so not representative):

| Quantity | Value |
|---|---|
| profile_alpha recovery | within 2.7% (α 0.5/1/2; LG and Q.3Di.AF); validated, no fix needed |
| α_AA | 1.51 |
| α_3Di | **1000 (grid bound)** |
| κ_AA | 1.03–1.45 |
| κ_3Di | 0.98–1.23 |
| ρ_AA,3Di | 0.08–0.35 |
| η_tdi | 0.052 |
| τ_misfit | 3 |
| B_E4 | 1000 |
| E5 power at f = 1.5 | 0.8–1.0 (synthetic world of tune size; real euk20 tune has 2 lineages) |

**Production projections:**

| Step | Time |
|---|---|
| tune_alpha | ~8 min |
| boot_sub alignment | ~2 min at 48 threads |
| SE bootstrap | < 1 min (GPU) |
| apply_tuning | ~5–10 min |
| τ / E4 timing / E5 | ~1–3 min each |
| **V3 at B = 200** | **3–3.5 h** (1.7 h at B = 100). Long pole; run it on node 2 or in parallel with the dry run. |

## 5. Endpoints (`src/strucevo/stats/driver.py`, `scripts/endpoints/{run_dryrun,run_confirm,interpret}.py`)

**Guards.** `run_confirm.py` checks ten guards before reading any data:
- tuning.yaml has no nulls;
- `git status --untracked-files=no` is clean, and tuning.yaml is tracked;
- the md5s of the plan, tuning.yaml and split.tsv match the ledger rows `freeze: …`;
- the split md5 in tuning.yaml matches;
- `CONFIRM_RUN=1`;
- this is not a rerun (`--rerun-logged-deviation ID` must name an ID in DEVIATIONS.md).

Tested: on the project root it refuses (exit 2) and lists every failure. `--test-root` is rejected on the project root.

**Sensitivities.** All §4 sensitivities run as variant inputs: other alignments, α×1.5 / α÷1.5 (variant tables), d_AA ≤ 1, hidden paralogy (C27), single-copy HOGs, C14 imputation, η free, and for E4: WLS, pooled slope, censored cells, internal edges, bin 5, aa_ml_seq ≤ 1. The E4 negative control (D13) runs on tune only. Writes E*.json (§4 key set), holm.json and interpretation.json.

**Test-root confirm (euk20, 383 HOGs).** Values are meaningless for euk20, but the plumbing works:

| Endpoint | Result |
|---|---|
| E1 | concave, ĉ = 0.86 |
| E2 | inconclusive |
| E3 | k = 4.5 |
| E4 | β = −0.23 |
| E5 | inconclusive |
| E6 | not detected |

- Interpretation row: none.
- Runtimes: 460 s at plan B; 57 s with `--quick`. Tune dry run: 485 s (56 s with `--quick`).
- **Production projection:** confirm ≈ 3.5 h (E1/E2 ≈ 70 min each, E4 ≈ 40 min); tune dry run ≈ 1–1.5 h.

**Unavailable sensitivities (flagged in the JSON):**
- E4 deep ages: no Strassert 2021 / Betts 2018 table is in the repo, so "robust to α and timescale" cannot be reached.
- All-aligned-columns: 3Di all-column counts are not stored.

**Mapping choices** (for a C-item): κ(bin) multiplies the SE columns before the engines; ρ is one count-weighted value per endpoint; HOG-level exclusion columns go to `data/interim/hogs_endpoint_cols_{split}.parquet` rather than rewriting hogs.parquet.

## 6. Exploratory (`scripts/exploratory/{x1_curves,x2_age,x7_pclu}.py`; tune only; `results/exploratory/X*_<aln>.json`, `figures/exploratory/`)

- X1: naive ĉ per metric, with 95% CIs from a HOG bootstrap.
- X2: path distance vs 2×MRCA age, plus node-to-tip. Root-to-tip is not identifiable on euk20, where every root is 1455 Mya.
- X7: ortholog curves built the same way as the P-clu binned curves. Orthologs sit about +0.09 TM above P-clu overall, but close to the P-clu rep-pLDDT ≥ 90 curve.
- Ages appear only as aggregates (bins with ≥ 3 nodes).
- Runtime per alignment: 10 s / 50 s / 12 s on euk20. Production tune: about 1 / 3–5 / 1 min.

## 7. Decisions for the PI (before freeze)

1. **α_3Di on tdiaa runs to the grid bound (1000).** seqonly and tmalign give 1.5–2. The likely cause is that tdiaa uses 3Di in the alignment score. As a result, α×1.5 is a no-op for 3Di. Options: take α_3Di from seqonly (needs a DEVIATIONS entry), or accept α = ∞.
2. **E1 SIMEX bias.** On real tune designs with calibrated κ(bin)/ρ(bin), a linear truth gives ĉ ≈ 0.94 (naive 0.82), with SE 0.018 at 1,440 HOGs. V3 will likely fail at production scale, which would make E1/E2 "inconclusive" by rule 1. The bias is probably heteroscedastic/binned noise that the scalar-ρ SIMEX does not remove. Synthetic V3 with scalar noise was unbiased (1.004). Needs investigation on the real tune split before the freeze.
3. **E4 null bias.** A null world gives β̂ = +0.05 to +0.08 (SE 0.02), significant in 5 of 6 worlds. This is likely a log-of-small-sums effect in young cells, and the D13 permutation control would not catch it.
4. The E5 β CI is uninformative at tune size (D8 note stands).
5. Exploratory and H1e seeds use spawn keys outside `config/seeds.yaml` (16, 1001–1007). Add the stages, or record them as a C-item.
6. Minor:
   - The engine counts saturated pairs as failures in fail_frac.
   - The 6 NA-age panel nodes are imputed (flag `age_imputed:6`).

## Files

| Kind | Paths |
|---|---|
| New scripts | `scripts/run_after_ancestra.sh`; `scripts/pipeline/{prep_pair_lists,eta}.py`; `scripts/h1e/run_h1e.py`; `scripts/tuning/*` (16 files); `scripts/endpoints/*`; `scripts/exploratory/*` |
| New library code and tests | `src/strucevo/pairs/store_files.py`; `src/strucevo/stats/driver.py`; `tests/test_tuning.py`; `config/tuning.yaml` (null template) |
| Modified | `scripts/pairs/{run_pairs,recompute_ml}.py`; `src/strucevo/pairs/{pipeline,posthoc,ml}.py`; `src/strucevo/stats/endpoints_lineage.py` (NaN-safe species pLDDT); `scripts/h1e/fetch_pdb_chains.py` |
| Test roots (gitignored) | `data/interim/itest_euk20/`, `data/interim/itest_euk20_orch/` |

pytest: `tests/test_pairs.py`, `tests/test_tuning.py` and `tests/test_stats_engine.py` give 33 passed.
