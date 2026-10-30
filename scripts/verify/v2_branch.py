"""V2 (ANALYSIS_PLAN §7): branch engine on simulated additive distances plus noise.
Pass: r > 0.99 (estimate vs truth), unbiased (|mean error| < 2 MC-SE and < 1% of mean length),
95% HC1 sandwich coverage in [0.92, 0.97].  Writes results/verify/V2.json."""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from strucevo.stats.common import stage_rng, write_json
from strucevo.stats.synth import random_panel, make_world, observe_gauss
from strucevo.stats.branch import fit_branches

ap = argparse.ArgumentParser()
ap.add_argument("--hogs", type=int, default=300)
ap.add_argument("--reps", type=int, default=20)
ap.add_argument("--rho", type=float, default=0.3)
ap.add_argument("--out", default=str(ROOT / "results" / "verify" / "V2.json"))
a = ap.parse_args()
t0 = time.time()
rng = stage_rng("v4_sim", 2)
panel = random_panel(100, rng=rng)
world = make_world(panel, a.hogs, rng=rng)
res = {"check": "V2", "n_hogs": a.hogs, "reps": a.reps, "rho": a.rho, "scenarios": {}}
for kap in (1.0, 1.5):
  for method in ("hc1", "hc3", "model"):
    rngs = stage_rng("v4_sim", 3)            # same noise across methods
    est, tru, se, metric_of = [], [], [], []
    for rep in range(a.reps):
        tab = observe_gauss(world, rng=rngs, rho=a.rho, nu_tdi=0.02, tip_sd=0.0, kappa=kap)
        bf = fit_branches(tab, panel, metrics=["aa_ml", "tdi_ml"], nu={"tdi_ml": 0.02}, misfit=False,
                          display_nnls=False, se_method=method)
        for hf in bf.hogs:
            ia, ib, Da, D3, be_aa, be_3, ht = world.pairs_true[hf.hog]
            for k, (m, be) in enumerate((("aa_ml", be_aa), ("tdi_ml", be_3))):
                if hf.fit["under"][k]:
                    continue
                truth = np.zeros(ht.n_params)
                np.add.at(truth, ht.param_of_edge, be)
                if m == "tdi_ml":   # observed 3Di distances include eta = 2 nu = nu on each terminal parameter
                    term = np.zeros(ht.n_params); np.add.at(term, ht.param_of_edge, ht.is_terminal * 0.02)
                    truth = truth + term
                b = hf.fit["b"][k]; s_ = np.sqrt(np.diag(hf.fit["V"][k]))
                est.append(b); tru.append(truth); se.append(s_); metric_of += [m] * len(b)
    est, tru, se = map(np.concatenate, (est, tru, se)); metric_of = np.array(metric_of)
    sc = res["scenarios"].setdefault(f"kappa={kap}", {})
    for m in ("aa_ml", "tdi_ml"):
        s = metric_of == m
        err = est[s] - tru[s]
        r = np.corrcoef(est[s], tru[s])[0, 1]
        bias = err.mean(); mcse = err.std() / np.sqrt(s.sum())
        cov = np.mean(np.abs(err) <= 1.959964 * se[s])
        relbias = err.sum() / tru[s].sum()
        sc.setdefault(method, {})[m] = {"n_params": int(s.sum()), "r": float(r), "mean_error": float(bias),
            "mcse": float(mcse), "rel_bias_of_sum": float(relbias), "coverage95": float(cov),
            "pass_r": bool(r > 0.99), "pass_unbiased": bool(abs(bias) < 2.5 * mcse or abs(relbias) < 0.01),
            "pass_coverage": bool(0.92 <= cov <= 0.97)}
prim = res["scenarios"]["kappa=1.0"]
res["pass_plan_hc1"] = all(v["pass_r"] and v["pass_unbiased"] and v["pass_coverage"] for v in prim["hc1"].values())
res["pass_unbiased_all"] = all(v["pass_unbiased"] for sc in res["scenarios"].values() for mm in sc.values() for v in mm.values())
res["runtime_s"] = time.time() - t0
write_json(res, a.out)
print(json.dumps(res, indent=1))
