"""V3 (ANALYSIS_PLAN §7): E1/E2 curvature engine under linear truth with realistic noise.

N datasets (plan: 100) of linear truth y = eta + r_h x on synthetic HOG trees (placeholder for real HOG designs),
realistic ML-engine SEs, correlated x/y errors (kappa, rho), per-protein tip noise; full E1 estimator (two-pass,
SIMEX, pigeonhole bootstrap, decision rules). Reports the false-concave rate (pass: <= 10%; otherwise E1/E2 are
'inconclusive') and the minimal detectable |c-1| at 80% power from concave-truth worlds.
Placeholder noise kappa=1, rho=0.3 until §Tuning; rerun with --kappa/--rho/--nu/--tip_sd once measured.
"""
import argparse, json, sys, time, warnings
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
from strucevo.stats.common import stage_rng, write_json, Z98
from strucevo.stats.synth import random_panel, make_world, observe_gauss
from strucevo.stats.endpoints_curve import prepare_curve_data, fit_curve_full, curve_decision

ap = argparse.ArgumentParser()
ap.add_argument("--datasets", type=int, default=100)
ap.add_argument("--hogs", type=int, default=200)
ap.add_argument("--kappa", type=float, default=1.0)
ap.add_argument("--rho", type=float, default=0.3)
ap.add_argument("--nu", type=float, default=0.02)
ap.add_argument("--tip_sd", type=float, default=0.01)
ap.add_argument("--B", type=int, default=30, help="bootstrap reps (plan: 200)")
ap.add_argument("--Bs", type=int, default=20, help="SIMEX B_s (plan: 50)")
ap.add_argument("--Bs_boot", type=int, default=5, help="SIMEX B_s inside bootstrap (plan: 10)")
ap.add_argument("--power_c", type=float, nargs="*", default=[0.97, 0.94, 0.90])
ap.add_argument("--power_datasets", type=int, default=20)
ap.add_argument("--out", default=str(ROOT / "results" / "verify" / "V3.json"))
a = ap.parse_args()
t0 = time.time()
panel = random_panel(150, rng=stage_rng("v3_sim", 0))


def one(i, c_true=None, tag=0):
    r = stage_rng("v3_sim", 1 + tag, i)
    w = make_world(panel, a.hogs, rng=r)
    tab = observe_gauss(w, rng=r, kappa=a.kappa, rho=a.rho, nu_tdi=a.nu, tip_sd=a.tip_sd, c_true=c_true)
    data, info, idx, sp = prepare_curve_data(tab, "E1")
    res = fit_curve_full(data, 2 * a.nu, a.kappa, a.kappa, a.rho, B=a.B, B_s=a.Bs, B_s_boot=a.Bs_boot, B_naive=a.B,
                         seed_stage="v3_sim", sp=sp, seed_extra=(100 + tag, i))
    cls, p, ex = curve_decision(res["c"], res["se_boot"])
    return {"c_naive": res["c_naive"], "c_simex": res["c"], "se": res["se_boot"], "class": cls,
            "n_pairs": info["n_pairs_fit"], "n_clip": res.get("n_clip", 0)}


runs = [one(i) for i in range(a.datasets)]
cls = [r["class"] for r in runs]
cs = np.array([r["c_simex"] for r in runs]); cn = np.array([r["c_naive"] for r in runs]); se = np.array([r["se"] for r in runs])
res = {"check": "V3", "settings": vars(a), "n_datasets": a.datasets,
       "false_concave_rate": float(np.mean([c == "concave" for c in cls])),
       "class_counts": {k: cls.count(k) for k in set(cls)},
       "mean_c_naive": float(cn.mean()), "mean_c_simex": float(cs.mean()), "sd_c_simex": float(cs.std(ddof=1)),
       "mean_se_boot": float(se.mean()), "coverage98_of_1": float(np.mean(np.abs(cs - 1) <= Z98 * se)),
       "mean_n_pairs": float(np.mean([r["n_pairs"] for r in runs]))}
power = {}
for j, c in enumerate(a.power_c):
    rr = [one(i, c_true=c, tag=10 + j) for i in range(a.power_datasets)]
    power[str(c)] = {"concave_rate": float(np.mean([r["class"] == "concave" for r in rr])),
                     "mean_c_simex": float(np.mean([r["c_simex"] for r in rr])),
                     "mean_c_naive": float(np.mean([r["c_naive"] for r in rr]))}
res["power"] = power
det = [abs(1 - float(c)) for c, v in power.items() if v["concave_rate"] >= 0.80]
res["min_detectable_abs_c_minus_1"] = min(det) if det else None
res["pass"] = res["false_concave_rate"] <= 0.10
res["runtime_s"] = time.time() - t0
write_json(res, a.out)
print(json.dumps(res, indent=1))
