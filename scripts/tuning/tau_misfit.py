"""§3.2 tau_misfit = max(3, q95 of misfit_ratio under parametric simulation from the species-tree fit plus noise).

1. Branch engine (alignment --aln, default tdiaa; metric aa_ml) on split == 'tune' primary pairs, calibrated SEs
   (kappa x Fisher SE), D2 model SEs. Observed misfit_ratio = RSS_rel(species tree) / RSS_rel(NJ) (branch.py §2.5.6).
2. For each fitted HOG (>= 4 species, not underdetermined), --reps replicates on the REAL design (same S_h, same
   pairs, same failed pairs): d* = X b_hat_AA + e, e ~ N(0, SE_cal^2 + s2_hm). Refit the two-pass WLS on the species
   tree, NJ on d* (failed pairs filled with fitted values, D4), WLS on the NJ topology -> simulated misfit_ratio.
3. tau_misfit = max(3, pooled q95). Also reported: per-HOG q95 distribution, observed distribution and the fraction
   of HOGs excluded at tau.
Writes ROOT/results/tuning/tau_misfit.json.
"""
import argparse
import time

import numpy as np

import tuning_common as tc
from strucevo.stats.branch import wls_two_pass
from strucevo.stats.trees import nj_splits

ap = tc.add_common_args(argparse.ArgumentParser())
ap.add_argument("--reps", type=int, default=50, help="parametric replicates per HOG")
a = ap.parse_args()
t0 = time.time()
root = a.root
tuned = tc.load_tuned(root, a.aln)
bf, panel, tab, info = tc.fit_tune_branches(root, a.aln, tuned, metrics=["aa_ml", "tdi_ml"])
ka = bf.metrics.index("aa_ml")
obs, sim_all, per_hog_q95, hogs_used = [], [], [], []
n_skip = {"underdetermined": 0, "lt4_species": 0}
t1 = time.time()
for i, hf in enumerate(bf.hogs):
    if hf.fit["under"][ka]:
        n_skip["underdetermined"] += 1; continue
    if len(hf.ht.species) < 4:
        n_skip["lt4_species"] += 1; continue
    rng = tc.tuning_rng("tau_misfit_sim", i)
    X, se, ok = hf.X, hf.se[ka], hf.ok[ka]
    b, s2 = hf.fit["b"][ka], hf.fit["s2"][ka]
    mu = X @ b
    sd = np.sqrt(np.where(ok, se * se + s2, 0.0))
    R = a.reps
    dsim = mu[None] + sd[None] * rng.standard_normal((R, len(mu)))
    okR = np.repeat(ok[None], R, 0)
    seR = np.repeat(se[None], R, 0)
    fs = wls_two_pass(X, dsim, seR, okR)
    nsp = len(hf.ht.species)
    ratios = np.full(R, np.nan)
    for r in range(R):
        if fs["under"][r]:
            continue
        fitted = X @ fs["b"][r]
        D = np.zeros((nsp, nsp))
        dfill = np.where(ok, dsim[r], fitted)
        D[hf.ia, hf.ib] = dfill; D[hf.ib, hf.ia] = dfill
        M = nj_splits(D)
        Xnj = (M[hf.ia] ^ M[hf.ib]).astype(float)
        fnj = wls_two_pass(Xnj, dsim[r:r + 1], se[None], ok[None])
        if np.isfinite(fnj["rss_rel"][0]) and fnj["rss_rel"][0] > 0:
            ratios[r] = fs["rss_rel"][r] / fnj["rss_rel"][0]
    sim_all.append(ratios)
    per_hog_q95.append(np.nanquantile(ratios, 0.95) if np.isfinite(ratios).any() else np.nan)
    obs.append(hf.misfit_ratio)
    hogs_used.append(hf.hog)
t_sim = time.time() - t1
sim = np.concatenate(sim_all) if sim_all else np.zeros(0)
sim = sim[np.isfinite(sim)]
obs = np.asarray(obs, float)
q = lambda v, p: float(np.quantile(v, p)) if len(v) else None
q95 = q(sim, 0.95)
tau = max(3.0, q95) if q95 is not None else None
obs_f = obs[np.isfinite(obs)]
phq = np.asarray(per_hog_q95, float)
rec = tc.header(root, "tau_misfit", a)
rec.update({
    "tau_misfit": tau, "rule": "max(3, q95 of simulated misfit_ratio)", "q95_sim": q95,
    "sim_quantiles": {p: q(sim, float(p)) for p in ("0.5", "0.9", "0.95", "0.99")},
    "n_sim": int(len(sim)), "reps_per_hog": a.reps, "n_hogs_sim": len(hogs_used), "n_skip": n_skip,
    "per_hog_q95_quantiles": {p: (float(np.nanquantile(phq, float(p))) if np.isfinite(phq).any() else None)
                              for p in ("0.05", "0.5", "0.95")},
    "observed": {"n": int(len(obs_f)), "quantiles": {p: q(obs_f, float(p)) for p in ("0.05", "0.25", "0.5", "0.75", "0.95", "0.99")},
                 "max": float(obs_f.max()) if len(obs_f) else None,
                 "frac_gt_3": float(np.mean(obs_f > 3)) if len(obs_f) else None,
                 "frac_excluded_at_tau": float(np.mean(obs_f > tau)) if (len(obs_f) and tau) else None,
                 "n_excluded_at_tau": int(np.sum(obs_f > tau)) if tau else None,
                 "frac_obs_gt_own_sim_q95": float(np.nanmean(obs > phq))},
    "noise": {"kappa": tuned["sources"]["kappa"], "s2_hm": "per-HOG two-pass s2 (branch engine)",
              "provisional": tuned["provisional"], "summary": tuned["summary"]},
    "fit_info": info, "runtime_s": {"total": time.time() - t0, "simulation": t_sim,
                                   "per_hog_sim_s": t_sim / max(1, len(hogs_used))},
})
tc.dump(rec, tc.out_path(root, "tau_misfit"))
print({k: rec[k] for k in ("tau_misfit", "q95_sim", "n_sim")}, rec["observed"], rec["runtime_s"])
