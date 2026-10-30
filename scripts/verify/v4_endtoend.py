"""V4 (ANALYSIS_PLAN §7): end-to-end synthetic world.

N0 world: b_3Di = r_h * b_AA on random HOG trees from a synthetic dated ~100-species tree; AA counts simulated
under LG+G4, 3Di counts under Q.3Di.AF+G4 (strucevo.simulate_counts), per-protein tip noise on 3Di; ML distances
(d_max=20); then branch engine -> E1, E3, E4, E5, E6 (N0 paralogs) -> Holm -> interpretation.
Spike-in world: same, with b_3Di x f on a random 10% of panel edges -> E5 must detect.
Duplication world (bonus): P2 accelerated -> E6 detected.

Pass: E1 98% CI covers 1; E4 98% CI covers 0; E5 N0 margins met (beta CI in [0.8,1.25] and sigma_L 99% upper
< 0.10); spike world detected; E6 (N0 world) CI covers 0.
Writes results/verify/V4.json and endpoint JSONs under results/verify/V4_endpoints/ (never results/endpoints/).
"""
import argparse, json, sys, time, warnings
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
from strucevo.stats.common import stage_rng, write_json
from strucevo.stats.synth import random_panel, make_world, observe_counts, observe_gauss, make_triplets
from strucevo.stats.branch import fit_branches
from strucevo.stats.endpoints_curve import run_curve_endpoint
from strucevo.stats.endpoints_branch import run_E3, run_E4
from strucevo.stats.endpoints_lineage import run_E5
from strucevo.stats.endpoints_misc import run_E6, apply_holm, write_interpretation

ap = argparse.ArgumentParser()
ap.add_argument("--hogs", type=int, default=1000)
ap.add_argument("--species", type=int, default=150)
ap.add_argument("--obs", choices=["counts", "gauss"], default="counts")
ap.add_argument("--B1", type=int, default=60, help="E1 bootstrap (plan: 200)")
ap.add_argument("--Bs", type=int, default=20, help="E1 SIMEX B_s (plan: 50)")
ap.add_argument("--B4", type=int, default=100, help="E4 bootstrap (plan: B_E4)")
ap.add_argument("--B", type=int, default=300, help="E3/E5/E6 bootstrap (plan: 1000)")
ap.add_argument("--P", type=int, default=50, help="E5 floor replicates (plan: 100)")
ap.add_argument("--spike_reps", type=int, default=20, help="spike-in replicates (plan: 50)")
ap.add_argument("--triplets", type=int, default=400)
ap.add_argument("--spike_f", type=float, default=1.5)
ap.add_argument("--min_t_snr", type=float, default=None, help="proposed deviation D8 (plan: None)")
ap.add_argument("--out", default=str(ROOT / "results" / "verify" / "V4.json"))
a = ap.parse_args()
t0 = time.time()
NU = 0.02
outdir = ROOT / "results" / "verify" / "V4_endpoints"
rng = stage_rng("v4_sim", 0)
panel = random_panel(a.species, rng=rng)
obs = observe_counts if a.obs == "counts" else observe_gauss
RHO = 0.0 if a.obs == "counts" else 0.3          # counts: independent AA/3Di draws; gauss default rho=0.3


def world_and_fit(spike_frac=0.0, seed_extra=0):
    r = stage_rng("v4_sim", 1)                                  # same world draws for N0 and spike
    w = make_world(panel, a.hogs, rng=r, spike_frac=spike_frac, spike_f=a.spike_f)
    tab = obs(w, rng=stage_rng("v4_sim", 2, seed_extra), nu_tdi=NU, tip_sd=0.01)
    bf = fit_branches(tab, panel, metrics=["aa_ml", "tdi_ml"], nu={"tdi_ml": NU}, se_method="model", misfit=True)
    return w, tab, bf


res = {"check": "V4", "settings": vars(a), "timing": {}}
w, tab, bf = world_and_fit()
res["timing"]["world_fit_s"] = time.time() - t0
true_k = sum(w.pairs_true[h][4].sum() for h in w.hogs) / sum(w.pairs_true[h][5].sum() for h in w.hogs)
res["n_pairs"] = int(len(tab["hog"])); res["frac_sat_tdi"] = float(np.mean(tab["tdi_sat"]))
# ---------------- E1
t = time.time()
rec1, curve, _ = run_curve_endpoint(tab, "E1", eta=2 * NU, kx=1.0, ky=1.0, rho=RHO, B=a.B1, B_s=a.Bs, B_s_boot=5,
                                     B_naive=a.B1, v3_false_concave=None)
write_json(rec1, outdir / "E1.json"); res["timing"]["E1_s"] = time.time() - t
# ---------------- E3
t = time.time()
rec3 = run_E3(bf, panel, B=a.B); write_json(rec3, outdir / "E3.json"); res["timing"]["E3_s"] = time.time() - t
# ---------------- E4
t = time.time()
rec4 = run_E4(bf, panel, B=a.B4, sensitivities=False); write_json(rec4, outdir / "E4.json"); res["timing"]["E4_s"] = time.time() - t
# ---------------- E6 (N0 paralogs)
t = time.time()
trip, ppairs = make_triplets(w, a.triplets, rng=stage_rng("v4_sim", 3), nu_tdi=NU)
rec6 = run_E6(trip, ppairs, curve, B=a.B); write_json(rec6, outdir / "E6.json"); res["timing"]["E6_s"] = time.time() - t
# ---------------- E5
t = time.time()
rec5, dat, core = run_E5(bf, panel, rec4, rho=RHO, nu_t=NU, P=a.P, B=a.B, B_cv=min(200, a.B), spike_reps=a.spike_reps,
                         spike_B=min(200, a.B), spike_fs=(1.5, 1.0), e6_significant=False, min_t_snr=a.min_t_snr)
write_json(rec5, outdir / "E5.json"); res["timing"]["E5_s"] = time.time() - t
holm = apply_holm(outdir); interp = write_interpretation(outdir)
# ---------------- spike world
t = time.time()
ws, tabs, bfs = world_and_fit(spike_frac=0.10, seed_extra=1)
rec4s = run_E4(bfs, panel, B=max(20, a.B4 // 4), sensitivities=False)
rec5s, _, cores = run_E5(bfs, panel, rec4s, rho=RHO, nu_t=NU, P=a.P, B=a.B, B_cv=min(200, a.B), spike_reps=0,
                         spike_fs=(), power_override=1.0, min_t_snr=a.min_t_snr)
res["timing"]["spike_world_s"] = time.time() - t
spike_detect = bool(cores["sigma2_L_lower99"] > 0 or not (cores["beta_ci98"][0] <= 1 <= cores["beta_ci98"][1]))
# ---------------- duplication world (bonus positive control for E6)
t = time.time()
tripd, ppd = make_triplets(w, a.triplets, rng=stage_rng("v4_sim", 4), m_aa=1.5, m_3=1.3, nu_tdi=NU)
rec6d = run_E6(tripd, ppd, curve, B=a.B)
res["timing"]["dup_world_s"] = time.time() - t

cov = lambda ci, v: bool(ci[0] <= v <= ci[1])
res["E1"] = {"c": rec1["estimate"], "ci98": rec1["ci98"], "c_naive": rec1["c_naive"], "decision": rec1["decision"],
             "covers_1": cov(rec1["ci98"], 1.0), "n_pairs": rec1["n_units"], "n_hogs": rec1["n_hogs"], "flags": rec1["flags"]}
res["E3"] = {"k": rec3["estimate"], "ci98": rec3["ci98"], "true_k": true_k, "covers_true": cov(rec3["ci98"], true_k),
             "sigma_k": rec3["sigma_k"], "true_sigma_k_approx": 0.3}
res["E4"] = {"beta": rec4["estimate"], "ci98": rec4["ci98"], "covers_0": cov(rec4["ci98"], 0.0), "p": rec4["p_value"],
             "flags": rec4["flags"], "n_cells": rec4["n_units"]}
res["E5"] = {"beta": rec5["estimate"], "beta_ci98": rec5["ci98"], "beta_ci98_hog_only": rec5["beta_ci98_hog_only"], "sigma2_L": rec5["sigma2_L"],
             "sigma_L_upper99": rec5["one_sided_bound"]["sigma_L_upper99"], "floor_F": rec5["floor_F"],
             "sigma2_obs": rec5["sigma2_obs"], "n_lineages": rec5["n_lineages"], "n_contrasts": rec5["n_contrasts"],
             "margins_met": rec5["decision_detail"]["margins_met"], "decision_with_gates": rec5["decision"],
             "power": rec5["power"]}
res["E5_spike_world"] = {"sigma2_L": cores["sigma2_L"], "sigma2_L_lower99": cores["sigma2_L_lower99"],
                         "beta_ci98": cores["beta_ci98"], "detected": spike_detect, "floor_F": rec5s["floor_F"],
                         "sigma2_obs": rec5s["sigma2_obs"]}
res["E6_N0"] = {"E6a": rec6["E6a"], "ci98_a": rec6["ci98"], "E6b": rec6["E6b"], "ci98_b": rec6["ci98_b"],
                "covers_0": cov(rec6["ci98"], 0.0) and cov(rec6["ci98_b"], 0.0), "n": rec6["n_units"], "p": rec6["p_value"]}
res["E6_dup_world"] = {"E6a": rec6d["E6a"], "E6b": rec6d["E6b"], "p": rec6d["p_value"], "detected_unadj": rec6d["p_value"] < 0.05}
res["holm"] = holm; res["interpretation"] = interp
res["pass"] = {"E1_covers_1": res["E1"]["covers_1"], "E4_covers_0": res["E4"]["covers_0"],
               "E5_N0_margins_met": res["E5"]["margins_met"], "spike_detected": spike_detect,
               "E6_N0_covers_0": res["E6_N0"]["covers_0"]}
res["pass_all"] = all(res["pass"].values())
res["runtime_s"] = time.time() - t0
write_json(res, a.out)
print(json.dumps(res, indent=1, default=float))
