# strucevo distance library: report (methods developer, 2026-10-30)

**Bottom line.** `src/strucevo/` gives pairwise sequence distances (p, Poisson, Kimura, and batched maximum-likelihood (ML) distance on the GPU under any reversible 20-state model with +G/+I), TM/LDDT transforms, batched Kabsch RMSD and y-vs-x curve models. All 33 pytest tests pass. On AliSim data the ML distance is unbiased within 1.6% for every d ≤ 2. It matches IQ-TREE's own pairwise ML distances to within 0.004% (worst case over 6,400 pairs). On the A100 it processes **1.27 M pairs/s**.

Global seed: 20261030 (sub-seeds are derived as `seed + offset` in each test).

## Install / use

```bash
source /opt/conda/etc/profile.d/conda.sh && conda activate <root>/conda_envs/main
pip install -e . --no-deps          # already done in conda_envs/main; tests also work via tests/conftest.py sys.path
python -m pytest tests -q            # 33 tests, ~2 min (benchmark + AliSim are marked `slow`; -m "not slow" ≈ 20 s)
python tests/alisim_sim.py --jobs 32 # regenerates AliSim/IQ-TREE data -> data/interim/tests/alisim/ (~5 min, cached)
```

## API

`strucevo.distances`
| Function | Purpose |
|---|---|
| `pair_counts(seqs_a, seqs_b, alphabet=AA_ORDER, masks=None)` | Builds a (N,20,20) table of aligned states for each pair. It counts only columns where both sequences have a valid letter (gaps, X, B, Z and similar are dropped) and the optional per-column mask is True. CPU speed ≈ 45 k pairs/s per core (L = 400). |
| `plddt_mask(plddt_a, plddt_b, 70)` | Per-column mask that keeps columns where both residues have pLDDT ≥ 70. |
| `simple_distances(counts)` / `p_distance`, `poisson_distance`, `kimura_distance` | p, −ln(1−p), and Kimura −ln(1−p−0.2p²). Undefined cases return NaN and set a flag. |
| `load_model('LG'/'WAG')`, `load_paml(path, alphabet)`, `RateModel(S, pi, alphabet)` | Reversible models. `.Q()` is normalised to mean rate 1. Also `.reorder(alphabet)` and `.to_paml()`. |
| `ml_distance(counts, model, alpha=1.0, ncat=4, pinv=0.0, d_max=10, device='cuda')` | Batched ML distance. Returns `MLResult` with fields `d, se, saturated, sat_cap, sat_identity, lr_unrelated, at_zero, n_sites, identity, random_identity, loglik, converged`. |
| `profile_alpha(counts, model, alphas)` | Total profile log-likelihood over a grid of Gamma shape values, for choosing one global α (see caveats). |
| `estimate_gtr(counts, alphabet, alpha)` / `counting_gtr(...)` | Estimates a GTR-20 model, for example for 3Di when no published matrix is available. |
| `simulate_counts(model, d, length, alpha, pinv)` | Simulates count matrices, for tests and benchmarks. |
| `one_minus(x)`, `neg_log(x)` | TM-score / LDDT transforms: 1−x and −ln x (clipped at 1e-6). |
| `kabsch_rmsd(X, Y, mask=None, return_transform=False)` | Batched (B,L,3) Kabsch RMSD in torch (CPU or GPU) with masks; proper rotations only. |

`THREEDI_ORDER = "ACDEFGHIKLMNPQRSTVWY"` is Foldseek's 3Di letter order. `AA_ORDER = "ARNDCQEGHILKMFPSTWYV"` is the PAML order.

`strucevo.models`: `fit_model(name, x, y)` takes `name` ∈ {`lin0` (y = bx), `lin` (y = a + bx), `sat` (y = a(1−e^{−bx})), `power` (y = ax^c)}. Related functions are `fit_models`, `compare` (ΔAIC and Akaike weights), `best_model(criterion='aic'|'bic')`, `out_of_range_error(x, y, q)` (fit on x < q-quantile, predict above it; returns RMSE, MAE and bias) and `divergence_x(f_lin, f_sat, rel=0.1)` (the d* of DESIGN §5). Nonlinear fits use variable projection: a log-grid over the nonlinear parameter with a closed-form amplitude, then `least_squares` refinement. This makes them robust to x scale; this is tested over a range of b from 0.002 to 60.

### How the ML distance works
- **Inputs.** Q comes from S and π, normalised so the mean rate is 1. The discrete Gamma uses category means (Yang 1994), which is the IQ-TREE/AliSim default. With +I, the variable-site rates are scaled by 1/(1−pinv).
- **Objective.** The log-likelihood is L(d) = Σ_ij N_ij log(π_i P_ij(d)), with P(d) = pinv·I + Σ_k w_k V e^{Λ r_k d} W. The symmetric eigendecomposition of Q is computed once.
- **Search.**
  1. L is evaluated on a 400-point log grid (0 and 1e-5 to d_max) as a single (N×400)·(400×G) matmul.
  2. The grid maximum gives a bracket.
  3. A safeguarded Newton search runs inside the bracket, using analytic P′ and P″ and falling back to bisection.
  4. All work is in float64, in chunks of 200 k pairs.
- **SE.** The standard error is 1/√(observed information) at d̂. It is NaN at the boundaries (d = 0 or d = d_max).
- **Saturation flag.** A pair is flagged `saturated` when any one of these holds:
  - d̂ hits d_max (`sat_cap`);
  - identity ≤ the identity expected for a random pair with the same composition (`sat_identity`);
  - the likelihood-ratio statistic against "unrelated" (P_ij = π_j) is below 2.706 (5% level, boundary χ² mixture).

  The function returns capped values and flags rather than inf.

### Matrix provenance
`src/strucevo/data/lg.dat` and `wag.dat` come from the PAML GitHub repository (abacus-gene/paml, `dat/`), fetched 2026-10-30. The sha256 hashes are in `src/strucevo/data/PROVENANCE.txt`. The WAG file is numerically identical to easel's bundled `wag.dat`. **No 3Di matrix is bundled; none was invented.** To use one, load a published matrix with `load_paml(path, THREEDI_ORDER)` after checking the file's letter order, or estimate one with `estimate_gtr`.

## Test results

### 1. AliSim validation (`tests/alisim_sim.py`, `tests/test_alisim_ml.py`)
- **Simulation.** Two-taxon trees `(A:d,B:0)`, simulated under `LG+G4{1}` and `WAG+G4{1}`, with L ∈ {300, 1000}, d ∈ {0.05 … 3} and 200 replicates per cell (6,400 pairs).
- **IQ-TREE reference.** IQ-TREE 3.1.4 requires ≥ 3 sequences, so I stacked 20 replicate pairs per alignment. I ran `-m X+G4{1} -t BIONJ -n 0`, which writes `.mldist` under the fixed model, and used only the within-pair entries.
- **Criteria (fixed before running):** |bias| < 3% for d ≤ 2; median and 99th-percentile |rel. diff| to IQ-TREE < 2%; no saturation flags; SE within 15% of the empirical SD.
- **Result:** all criteria pass.

| Model | L | d | mean d̂ (ours) | mean d̂ (IQ-TREE) | bias % | MC-SE % | median \|diff\| vs IQ-TREE | max \|diff\| vs IQ-TREE | mean SE / empirical SD | mean p | mean Kimura |
|---|---|---|---|---|---|---|---|---|---|---|---|
| LG | 300 | 0.05 | 0.0495 | 0.0495 | -1.00 | 1.99 | 1.1e-04% | 3.5e-03% | 0.96 | 0.047 | 0.048 |
| LG | 300 | 0.1 | 0.1008 | 0.1008 | +0.75 | 1.38 | 3.9e-05% | 1.1e-03% | 1.04 | 0.090 | 0.096 |
| LG | 300 | 0.25 | 0.2523 | 0.2523 | +0.91 | 1.07 | 2.1e-05% | 3.8e-04% | 0.99 | 0.194 | 0.226 |
| LG | 300 | 0.5 | 0.5028 | 0.5028 | +0.56 | 0.85 | 5.4e-06% | 1.5e-04% | 1.06 | 0.317 | 0.412 |
| LG | 300 | 1 | 1.0008 | 1.0008 | +0.08 | 0.82 | 2.9e-06% | 9.0e-05% | 1.02 | 0.464 | 0.710 |
| LG | 300 | 1.5 | 1.4999 | 1.4999 | -0.01 | 0.84 | 3.6e-06% | 5.4e-05% | 1.00 | 0.553 | 0.957 |
| LG | 300 | 2 | 2.0270 | 2.0270 | +1.35 | 0.84 | 1.6e-06% | 5.9e-05% | 1.05 | 0.614 | 1.176 |
| LG | 300 | 3 | 3.0108 | 3.0108 | +0.36 | 0.88 | 1.9e-06% | 3.1e-05% | 1.05 | 0.689 | 1.545 |
| LG | 1000 | 0.05 | 0.0500 | 0.0500 | -0.09 | 1.07 | 1.1e-04% | 2.4e-03% | 0.99 | 0.047 | 0.049 |
| LG | 1000 | 0.1 | 0.1004 | 0.1004 | +0.39 | 0.82 | 3.9e-05% | 1.1e-03% | 0.96 | 0.090 | 0.096 |
| LG | 1000 | 0.25 | 0.2482 | 0.2482 | -0.71 | 0.59 | 2.8e-05% | 2.7e-04% | 0.97 | 0.193 | 0.224 |
| LG | 1000 | 0.5 | 0.5015 | 0.5015 | +0.30 | 0.49 | 5.2e-06% | 1.1e-05% | 1.00 | 0.317 | 0.412 |
| LG | 1000 | 1 | 1.0036 | 1.0036 | +0.36 | 0.47 | 2.9e-06% | 7.7e-05% | 0.99 | 0.464 | 0.709 |
| LG | 1000 | 1.5 | 1.5109 | 1.5109 | +0.73 | 0.47 | 5.6e-06% | 4.8e-05% | 0.99 | 0.553 | 0.954 |
| LG | 1000 | 2 | 1.9925 | 1.9925 | -0.38 | 0.47 | 1.5e-06% | 5.1e-05% | 1.00 | 0.612 | 1.162 |
| LG | 1000 | 3 | 2.9932 | 2.9932 | -0.23 | 0.50 | 2.2e-06% | 2.7e-05% | 0.99 | 0.690 | 1.540 |
| WAG | 300 | 0.05 | 0.0508 | 0.0508 | +1.56 | 1.81 | 1.3e-04% | 2.7e-03% | 1.06 | 0.048 | 0.050 |
| WAG | 300 | 0.1 | 0.0988 | 0.0988 | -1.15 | 1.42 | 4.2e-05% | 1.2e-03% | 1.00 | 0.089 | 0.095 |
| WAG | 300 | 0.25 | 0.2536 | 0.2536 | +1.45 | 1.04 | 2.7e-05% | 5.4e-04% | 1.01 | 0.197 | 0.229 |
| WAG | 300 | 0.5 | 0.5040 | 0.5040 | +0.81 | 0.84 | 4.9e-06% | 2.0e-04% | 1.06 | 0.323 | 0.422 |
| WAG | 300 | 1 | 1.0046 | 1.0046 | +0.46 | 0.82 | 3.4e-06% | 1.2e-04% | 1.04 | 0.473 | 0.733 |
| WAG | 300 | 1.5 | 1.5043 | 1.5043 | +0.29 | 0.88 | 2.7e-06% | 5.5e-05% | 0.98 | 0.559 | 0.977 |
| WAG | 300 | 2 | 2.0263 | 2.0263 | +1.31 | 0.84 | 2.0e-06% | 5.1e-05% | 1.06 | 0.624 | 1.219 |
| WAG | 300 | 3 | 3.0698 | 3.0698 | +2.33 | 1.02 | 1.6e-06% | 3.6e-05% | 0.95 | 0.702 | 1.627 |
| WAG | 1000 | 0.05 | 0.0502 | 0.0502 | +0.45 | 1.01 | 8.9e-05% | 2.2e-03% | 1.05 | 0.048 | 0.049 |
| WAG | 1000 | 0.1 | 0.0987 | 0.0987 | -1.29 | 0.79 | 1.2e-04% | 1.1e-03% | 0.99 | 0.089 | 0.095 |
| WAG | 1000 | 0.25 | 0.2494 | 0.2494 | -0.23 | 0.58 | 2.3e-05% | 1.8e-04% | 0.97 | 0.195 | 0.227 |
| WAG | 1000 | 0.5 | 0.4996 | 0.4996 | -0.09 | 0.48 | 5.3e-06% | 1.1e-05% | 1.00 | 0.321 | 0.418 |
| WAG | 1000 | 1 | 1.0139 | 1.0139 | +1.39 | 0.47 | 2.4e-06% | 1.2e-04% | 0.99 | 0.474 | 0.733 |
| WAG | 1000 | 1.5 | 1.4929 | 1.4929 | -0.48 | 0.44 | 2.7e-06% | 2.4e-05% | 1.05 | 0.560 | 0.976 |
| WAG | 1000 | 2 | 2.0041 | 2.0041 | +0.21 | 0.44 | 1.6e-06% | 5.0e-05% | 1.11 | 0.622 | 1.204 |
| WAG | 1000 | 3 | 3.0082 | 3.0082 | +0.27 | 0.52 | 9.9e-07% | 2.0e-05% | 1.00 | 0.701 | 1.610 |

MC-SE is the Monte-Carlo standard error of the mean, relative to d. Most biases are within 2 MC-SE of 0. The largest are WAG L = 1000, d = 1 at +1.39 ± 0.47% (3.0 MC-SE, the only one of 32 cells above 2.3 MC-SE) and WAG L = 300, d = 3 at +2.3 ± 1.0%. IQ-TREE gives identical numbers in both cells, so these are properties of the ML estimator or Monte-Carlo noise, not implementation errors. The last column is a useful side result for the PI: Kimura's correction recovers only about 50% of the true distance at d = 3 under LG+G4, and p-distance only about 23%. **Do not use p or Kimura as the time axis beyond d ≈ 0.3.**

### 2. Unit tests (`tests/test_distances.py`, 22 tests, ~10 s)
| Test | Result |
|---|---|
| Model properties: rows of Q sum to 0, mean rate 1, detailed balance, πP = π (LG, WAG); LG spot values; PAML round-trip | pass |
| Γ category rates for α = 1 equal the published values 0.1370, 0.4768, 1.0000, 2.3863; mean rate 1 with +I | pass |
| Gaps and X are excluded; per-column mask works; pLDDT mask handles NaN; p, Poisson and Kimura give the expected values; Kimura NaN case is flagged | pass |
| Identical sequences give d = 0, `at_zero`, not saturated; a pair with no comparable sites gives NaN and a flag | pass |
| **Random unrelated pairs** (sampled from LG π; L = 100/300/1000; 1000 pairs each) are flagged saturated in **100%** of cases. All hit d_max; 51% also have identity ≤ the random expectation. | pass (criterion ≥ 95%) |
| Related pairs at d = 3, L = 300 are flagged saturated in less than 1% of cases (observed 0%) | pass |
| GPU result matches an independent scipy reference (`expm` + bounded optimiser) under WAG+G4(α = 0.7)+I(0.1) to within 1e-4 relative | pass |
| SE ≈ empirical SD of d̂ (within 8%) at d = 0.1, 0.5 and 1.5 | pass |
| +I path is unbiased (< 3%) at α = 0.8, pinv = 0.15 | pass |
| **3Di alphabet**: random synthetic GTR over `ACDEFGHIKLMNPQRSTVWY`; sequences simulated in Python, then strings → `pair_counts` → ML. d̂ is unbiased within 3% for d = 0.1 to 2. | pass |
| **`estimate_gtr`** on 4000 simulated pairs (synthetic GTR, d ~ U(0.05, 0.6), L = 300, α = 1 known): correlation of log S with the truth is 0.995 for ML and 0.979 for the counting estimator. Held-out distances under the estimated model have bias −0.1%, −0.1%, +0.1% and +1.2% at d = 0.1, 0.5, 1 and 2. The counting estimator alone is biased by −1.4% to −4.1%. | pass (criteria r > 0.9, \|bias\| < 5%) |
| Kabsch: random rotation plus translation gives RMSD < 1e-8 and recovers the transform; a mirror image is not superimposable; isotropic noise σ = 0.5 Å (L = 500) gives mean RMSD equal to σ√((3L−6)/L) within 1%; results match a per-structure numpy SVD reference to 1e-9; masked corrupted atoms are ignored; CPU and GPU agree; a symmetric-stretch analytic case gives exactly √(2e²/6) | pass |

### 3. Curve models (`tests/test_models.py`, 6 tests)
- Noiseless parameter recovery for all four models (rel 1e-5), including extreme x scales.
- Exact AIC and BIC formulas: k + 1 parameters with σ² counted.
- Out-of-range error: on saturating data, `lin0` fitted on the lower half over-predicts the upper half (bias > 0.05), and `sat` has the lower RMSE.
- `divergence_x` matches the analytic root.
- **AIC model selection** (100 seeds, n = 100 points, x ~ U(0.02, 3), σ = 0.05; generators y = 0.3x and y = 0.9(1−e^{−0.8x})):

| Comparison | linear data, correct | saturating data, correct |
|---|---|---|
| AIC, lin0 vs sat (asserted ≥ 90%) | **95 / 100** | **100 / 100** |
| BIC, lin0 vs sat | 99 / 100 | 100 / 100 |
| AIC, 4 models (reported only) | lin0 76, lin 16, power 7, sat 1 | sat 100 |

In the 4-model setting, linear data choose a model of the linear family (lin0 / lin / power with c ≈ 1) 99% of the time. Picking a nested model with one extra parameter about 16% of the time is expected behaviour of AIC; use BIC or Akaike weights when deciding "linear vs saturating".

### 4. Throughput (`tests/test_benchmark.py`)
The benchmark used 1,000,000 realistic count matrices: LG+G4 at d ~ U(0.01, 4), L ~ U(50, 800), with float64 throughout.

| Device | Setting | Pairs/s |
|---|---|---|
| A100-PCIE-40GB | chunk 100 k | 1.13 M |
| A100-PCIE-40GB | chunk 200 k (default) | 1.24 M |
| A100-PCIE-40GB | chunk 400 k | **1.27 M** (0.79 s for 1 M pairs) |
| CPU, 48 threads (torch) | — | 8.3 k |

Of the 1 M pairs, 100% converged and 0.002% were flagged saturated. Building the count matrices from strings with `pair_counts` runs at about 45 k pairs/s per core and is therefore the bottleneck: parallelise it over processes, or build counts directly from encoded int arrays.

## Caveats
1. **α is fixed by the user.** Pairwise data identify the Gamma shape only weakly. `profile_alpha` exists but was not validated. DESIGN's "global α fit" should come from an IQ-TREE fit on a few multi-sequence ortholog alignments, not from pairs.
2. **No 3Di matrix is bundled.** Estimating a 3Di GTR with `estimate_gtr` requires α to be specified, and both it and the counting estimator were validated only on synthetic data from a reversible stationary model. 3Di "evolution" may violate both assumptions, and the counting estimator is biased by 1–4% in downstream distances. Prefer a published matrix if one can be obtained.
3. **The model's π is used, not each pair's composition.** Compositionally skewed pairs (common for 3Di, possible for low-complexity proteins) will be biased. `sat_identity` uses the pair's own composition.
4. **Saturation.** Under α = 1, random pairs hit the cap d_max = 10; IQ-TREE caps at 9. Estimates above about 4–5 have large SE, so use the `se` and `saturated` flags rather than raw large d. The flag combines three criteria; for strict analyses, filter on the components.
5. **The ML distance treats the alignment as known.** Gap columns are discarded, so alignment errors at high divergence bias d̂ downward. The DESIGN's sensitivity analysis with a second alignment remains necessary.
6. AliSim uses the same discrete-Γ model as the estimator, so the simulation test is a correctness check, not a robustness check against model misspecification.

## Files
- `src/strucevo/{__init__,distances,models}.py`, `src/strucevo/data/{lg.dat,wag.dat,PROVENANCE.txt}`, `pyproject.toml`
- `tests/{conftest,alisim_sim,test_alisim_ml,test_distances,test_models,test_benchmark}.py`
- Regenerable outputs (gitignored): `data/interim/tests/alisim/` (45 MB) and `data/interim/tests/results/*.json`
