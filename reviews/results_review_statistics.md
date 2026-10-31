# Statistical review of the confirmatory results (independent reviewer, 2026-10-31)

All diagnostics are read-only on the confirm split and are point estimates unless stated. Scripts are in `scripts/review/` and outputs in `logs/review/`.

## Headline

The pre-registered decisions are computed correctly and should stand as reported. One mechanism explains four results:
- E3's low k and its rise with stratum;
- E4;
- E5's excess;
- most of E1/E2's concavity.

That mechanism is **a per-protein 3Di excess of ≈ 0.014–0.016 3Di-ML units that the identical-sequence floor (ν = η/2 = 0.026) does not remove**. Two independent estimates agree:
- branch level: the ν that equalises terminal and internal k is 0.040;
- pairwise level: the median HOG-specific intercept is 0.085, against η = 0.052.

So "3Di saturates" (E1/E2), "deep branches are structurally conservative" (E4) and "k = 3.3" are not supported substantively.

## Q1. E4 is driven by terminal edges

**Edge-type decomposition** (k = ΣAA/Σ3Di, tip-corrected):

| age bin (My) | terminal k | internal k | ν share of terminal raw 3Di |
|---|---|---|---|
| <30 | **0.73** | 3.13 | 61% |
| 30–100 | 2.61 | 4.47 | 41% |
| 100–300 | 3.57 | 4.63 | 31% |
| 300–1000 | 4.16 | 4.94 | 24% |
| >1000 | – | 4.68 | – |

- Internal k is flat (4.5–4.9) from 30 My to more than 1 Gy.
- Young terminal edges show more 3Di change than AA change.
- The terminal excess per edge is roughly constant: 0.014, 0.017, 0.013, 0.008. That points to an additive per-tip artefact.

**E4 point β against ν** (production value −0.150):

| ν | 0 | 0.013 | **0.026** | 0.039 | 0.040 | 0.052 |
|---|---|---|---|---|---|---|
| β | −0.35 | −0.29 | **−0.15** | **+0.07** | **+0.08** | −0.08 (20% censored) |

**The sign flips inside the plausible range of ν.**
- The *mean* identical-pair tdi_ml is 0.0626, against a median of 0.0524. The mean is the relevant value for ratio-of-sums, and it gives ν = 0.031.
- Non-identical orthologs probably have noisier AF models than identical sequences, which share their MSA.
- The ν ± 0.005 sensitivity promised in C35 is **missing from E4.json**.
- Fine grid (`e4_nu_grid.json`), monotone in ν:

| ν | 0.021 | 0.026 | 0.029 | 0.031 (mean-based) | 0.034 | 0.036 | 0.040 |
|---|---|---|---|---|---|---|---|
| β | −0.213 | −0.150 | −0.116 | **−0.079** | −0.034 | +0.010 | +0.081 |

  - Under C35's ±0.005, β ranges from −0.21 to −0.10.
  - With the mean-based ν (0.031), β halves to −0.079.
  - β reaches zero at ν ≈ 0.036.
  - Internal-only β stays at −0.047 for every ν, as it should, since internal edges are not tip-corrected.

**By stratum.** β is −0.40 in stratum 2, −0.19 in stratum 3 and −0.12 in stratum 4. It is largest where terminal edges are shortest.

**Internal-only fit.** β = −0.047, 98% CI [−0.128, 0.034].
- With SE ≈ 0.035, a true −0.15 would have been detected at more than 4 SE.
- So this is a change in the point estimate, not a loss of power.
- The tune dry run already showed it: internal-only β was −0.014.

**Is E4 the time-domain image of E1?**
- Analytically: under D_3Di ∝ D_AA^c on a clock-like tree, the branch ratio scales as τ^(c−1). That gives a slope of −0.28; times the V3 attenuation of 0.46, the expected β̂ is ≈ −0.13.
- Simulation: synthetic 3Di = η + a_h·d_AA^c on the real confirm design, run through the production branch engine and cells:

| world | β all | β internal | internal k, bin 0→4 |
|---|---|---|---|
| c = 1 | +0.012 | +0.004 | flat |
| c = 0.72 | **−0.149** | **−0.138** | rises 2.4→5.3 |
| c = 1 + 0.014/tip | **−0.134** | **+0.021** | flat; terminal bin-0 k **0.71** |
| observed | −0.150 | −0.047 | flat (4.5–4.9 after 30 My) |

- Both artefact worlds reproduce the all-edge β.
- Only the tip-excess world reproduces both the flat internal k and the terminal k of 0.71–0.73.
- Pure concavity predicts an internal β of −0.138, which lies outside the observed 98% CI.
- **Verdict:** E4 is mostly tip artefact. Any residual is ≤ −0.05, which is consistent with c ≈ 0.9–1.0.

## Q2. E1/E2 concavity: a robust estimator output, not saturation

**AA under-correction.** Heterotachy or an LG+Γ misfit compresses large d_AA, so it biases towards *convexity*, not concavity. The α×1.5 and ÷1.5 sensitivities move ĉ by ≤ 0.03.

**3Di under-correction or saturation.** This would predict concavity that grows with d. The opposite is observed:
- d_AA ≤ 1 is *more* concave: 0.66, and 0.53 for E2.
- All aligned columns, which are noisier, give 0.54.

**AF homogenisation (X4).**
- A multiplicative compression is absorbed by a_h.
- An offset that shrinks with d_AA makes AF curves *more convex* than PDB curves (AF/PDB slope ratio 1.34).
- So X4 cannot create concavity (n = 197, thin beyond d_AA = 2).

**SIMEX clipping.** 218k clips out of 118M draws (0.18%).
- SIMEX moves ĉ up: naive 0.670 → 0.717.
- The window x ∈ [0.3, 2] has no clipping and gives naive ĉ = 0.73.
- Negligible.

**The η floor is the problem.** With one global η the power law wins: the Q profile has its minimum at η ≈ 0.052–0.06, and the plan's "η free" check found 0.053. But the floor varies by family: by length alone, η ranges 0.037–0.066. Allowing a **HOG-specific intercept** changes the naive fits:

| window | ĉ (η fixed) | ĉ (η_h free) | Q: η fixed → η_h free | Q(c=1) − Q(ĉ), η_h free |
|---|---|---|---|---|
| 0.05–2 | 0.670 | **0.903** | 346,897 → 304,764 | 496 |
| 0.05–1 | 0.615 | **0.941** | 244,672 → 216,321 | 69 |
| 0.3–2 | 0.725 | **0.982** | 174,856 → 156,004 | 3 |

**Checks on the η_h-free fit.**
- Median η_h at c = 1 is 0.085 (q10–q90 0.042–0.142).
- A noise-free c = 0.72 world is recovered exactly (0.720), so the free intercept does not force c towards 1.
- A noise-free *linear* world with a +0.028 offset gives ĉ = 0.80 under the production fixed-η fit; a +0.06 offset gives 0.62.
- Naive ĉ is biased low: in V3, linear truth gave naive 0.89 against SIMEX 1.03. So the η_h-free results are compatible with **linearity above d_AA ≈ 0.05**.

**What the elevated baseline (~0.03 per non-identical pair) could be.** Either:
1. AF noise for different sequences, above the identical-sequence floor (an artefact); or
2. a genuine structural component that saturates quickly below d_AA ≈ 0.05.

The data cannot distinguish these. E2 shows the same signature: it is more concave at d ≤ 1, and d* sits at the grid minimum. E2 was not refitted here.

## Q3. E5 and E6

**E5.** *artefact_suspect* is correct.
- With internal edges only, the excess disappears: σ²_L 98% CI [−0.022, 0.117], β = 1.06 [0.79, 1.34].
- Power is 0.24, so neither N0 nor its rejection can be claimed.
- What can be claimed: no evidence of structure-specific lineage rates on internal branches. The apparent excess sits on species-level terminal edges.

**E6a.**
- It flips under seqonly: −0.011 [−0.024, 0.002].
- It compares pairs at different d_AA against the misspecified fixed-η curve, so it is not claimable.

**E6b.** MDO − LDO is positive under all three alignments (0.017, 0.032, 0.032), and the lddt_D version gives 0.035.
- It can be claimed directionally, as consistent with PMID 40813247.
- But its size matches the per-protein AF excess. A more-diverged copy with lower pLDDT or a shallower MSA would produce it.
- An adjustment for the MDO–LDO pLDDT difference is required before a biological reading.
- "E6 detected" must be labelled **not alignment-robust**.

## Q4. E3

**σ_k = 0 is an estimator failure.** The plan subtracts the *mean* v_h, which is outlier-dominated: mean 0.82, median 0.0039, q99 0.36. As a result σ²_k is −0.41 on the confirm split but +0.375 on tune.

Better estimates:
- DerSimonian–Laird τ = 0.32.
- Precise HOGs (v_h < 0.02, n = 4,214): σ_k = 0.47, with k_h q10–q90 = 1.7–5.2.
- Internal edges only: ≈ 0.2. This is a lower bound, because the plan's v_h uses whole-HOG variances.

So between-family variation of about 2–3-fold is real, and it does not contradict Illergård.

**The stratum trend is the tip artefact.**
- Internal k is flat across strata: 3.74, 4.71, 4.63, 4.65.
- Terminal k rises from 1.54 to 3.53.

**k_internal = 4.63 [4.35, 4.93] is the defensible central value (post hoc).**
- The CATH α/β contrast (4.68 vs 2.55) should be recomputed on internal edges.
- X4's slope ratio suggests k from experimental structures could be higher.

## Q5. Defensible conclusions

**(a) Linear?**
- *Shown:*
  - The pre-registered estimator finds sublinear growth (ĉ 0.70–0.72), robust to alignment and every planned sensitivity.
  - Post hoc, the curvature largely disappears with family-specific baselines (0.90–0.98).
  - The structure/AA ratio on internal branches is constant from 30 My to more than 1 Gy.
- *Believed:* roughly linear above a small baseline jump between non-identical orthologs. The cause of that jump is unresolved.

**(b) 3–10× more conserved?**
- *Shown:* k = 3.3 [3.1, 3.6] as planned.
- *Believed:* about 4.6 (internal edges), varying 2–3-fold across families. Both values are inside Illergård's range.

**(c) Clock or neutral?**
- *Shown:* no robust depth trend and no lineage-specific structural rate on internal branches. Both are absence-of-evidence results (E5 power 0.24; internal E4 SE 0.035).
- Neutrality is not addressed.
- *Believed:* a rescaled AA clock on internal branches.

## Caveats that MUST appear in the report

1. The tip correction under-corrects by about 0.014 per protein. E3's k and stratum trend, E4 and E5 all depend on it. Show the edge-type and ν tables, and add the missing C35 ν ± 0.005 sensitivity.
2. E4 is not robust to the tip correction. Its sign flips at ν ≥ 0.039, and the internal-only fit is not significant.
3. E1/E2 concavity depends on the global identical-sequence η. With HOG-specific baselines ĉ is 0.90–0.98. Do not read it as saturation: d ≤ 1 being more concave argues against saturation.
4. σ_k = 0 is an estimator failure. Report σ_k ≈ 0.3–0.5 (all edges) and ≈ 0.2 (internal edges), labelled post hoc.
5. E6a is not alignment-robust. E6b is not adjusted for model-quality differences between copies.
6. Power limits: E5 power 0.24; c = 0.90 detected only 55% of the time; E4 attenuation ≈ 0.46.
7. AF distances are 0.6–0.7× PDB distances (X4: n = 197, mostly human–mouse).
8. Everything in this review is post hoc. The pre-registered row stays R6_artefact_suspect.
