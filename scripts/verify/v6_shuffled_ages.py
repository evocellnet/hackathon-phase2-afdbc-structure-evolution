"""V6 (ANALYSIS_PLAN §4.4 negative control, §7 V6): shuffled node ages give no E4 slope.

Synthetic world WITH a true age effect (y = log(3Di/AA) slope beta_true on centred log age) so that the control
is non-trivial: (i) true ages -> E4 must detect the slope; (ii) 20 permutations of internal-node ages (median, lo,
hi permuted together across internal nodes; draws without the monotone repair) -> the 98% CI must cover 0 in
>= 18/20. Writes results/verify/V6.json.
"""
import argparse, json, sys, time, warnings
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
from strucevo.stats.common import stage_rng, write_json
from strucevo.stats.synth import random_panel, make_world, observe_gauss
from strucevo.stats.branch import fit_branches
from strucevo.stats.endpoints_branch import run_E4

ap = argparse.ArgumentParser()
ap.add_argument("--hogs", type=int, default=600)
ap.add_argument("--beta_true", type=float, default=-0.15)
ap.add_argument("--perms", type=int, default=20)
ap.add_argument("--B", type=int, default=40, help="E4 bootstrap per run (plan: B_E4)")
ap.add_argument("--mode", choices=["nodes", "edges"], default="nodes", help="nodes = plan (internal node ages permuted); edges = proposed D13 (edge midpoints permuted across all HOG edges)")
ap.add_argument("--out", default=str(ROOT / "results" / "verify" / "V6.json"))
a = ap.parse_args()
t0 = time.time()
panel = random_panel(150, rng=stage_rng("v4_sim", 10))
w = make_world(panel, a.hogs, rng=stage_rng("v4_sim", 11), age_beta=a.beta_true)
tab = observe_gauss(w, rng=stage_rng("v4_sim", 12), rho=0.3, nu_tdi=0.02, tip_sd=0.01)
bf = fit_branches(tab, panel, metrics=["aa_ml", "tdi_ml"], nu={"tdi_ml": 0.02}, se_method="model", misfit=False)
r0 = run_E4(bf, panel, B=a.B, sensitivities=False)
res = {"check": "V6", "settings": vars(a),
       "true_ages": {"beta": r0["estimate"], "ci98": r0["ci98"], "p": r0["p_value"],
                     "detected": not (r0["ci98"][0] <= 0 <= r0["ci98"][1])}}
internal = np.flatnonzero(~panel.is_leaf)
perm_res = []
for k in range(a.perms):
    rng = stage_rng("age_perm", k)
    if a.mode == "nodes":
        perm = rng.permutation(internal)
        age, lo, hi = panel.age.copy(), panel.lo.copy(), panel.hi.copy()
        age[internal], lo[internal], hi[internal] = panel.age[perm], panel.lo[perm], panel.hi[perm]
        r = run_E4(bf, panel, B=a.B, sensitivities=False, ages_override=(age, lo, hi), repair_draws=False)
    else:
        nE = int((bf.edges["metric"] == "aa_ml").sum())
        r = run_E4(bf, panel, B=a.B, sensitivities=False, tau_perm=rng.permutation(nE))
    perm_res.append({"beta": r["estimate"], "ci98": r["ci98"], "covers_0": bool(r["ci98"][0] <= 0 <= r["ci98"][1])})
n_cov = sum(p["covers_0"] for p in perm_res)
res["permutations"] = perm_res
res["n_cover_0"] = n_cov
res["pass"] = bool(n_cov >= int(np.ceil(0.9 * a.perms)) and res["true_ages"]["detected"])
res["runtime_s"] = time.time() - t0
write_json(res, a.out)
print(json.dumps({k: v for k, v in res.items() if k != "permutations"}, indent=1))
