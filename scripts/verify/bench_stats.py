"""Full-scale runtime benchmark of the statistics layer (synthetic Gaussian data at production size).
Times: branch engine (HOGs x 3 metrics), E1 point + SIMEX + per-bootstrap-replicate cost, E4 cell build + fast
REML per replicate and statsmodels point fit, E5 floor/bootstrap. Extrapolates to the plan's B values.
Writes results/verify/bench_stats.json."""
import argparse, json, sys, time, warnings
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
from strucevo.stats.common import stage_rng, write_json
from strucevo.stats.synth import random_panel, make_world, observe_gauss
from strucevo.stats.branch import fit_branches
from strucevo.stats.endpoints_curve import prepare_curve_data, fit_curve_full
from strucevo.stats.endpoints_branch import edge_view, build_cells, cell_design, fit_mixed, run_E3
from strucevo.stats.reml import fast_mixed
from strucevo.stats.endpoints_lineage import E5Data, parametric_floor, e5_core

ap = argparse.ArgumentParser()
ap.add_argument("--hogs", type=int, default=5600)
ap.add_argument("--out", default=str(ROOT / "results" / "verify" / "bench_stats.json"))
a = ap.parse_args()
T = {}
t = time.time()
panel = random_panel(170, rng=stage_rng("v4_sim", 20))
w = make_world(panel, a.hogs, rng=stage_rng("v4_sim", 21), min_sp=8, max_sp=30)
tab = observe_gauss(w, rng=stage_rng("v4_sim", 22))
tab["lddt_D"] = tab["tdi_ml"] * 0.7; tab["lddt_D_se"] = tab["tdi_ml_se"]; tab["lddt_sat"] = tab["tdi_sat"]
T["build_world_s"] = time.time() - t
T["n_pairs"] = int(len(tab["hog"]))
t = time.time()
bf = fit_branches(tab, panel, metrics=["aa_ml", "tdi_ml", "lddt_D"], nu={"tdi_ml": 0.02, "lddt_D": 0.01}, se_method="model")
T["branch_engine_3metrics_s"] = time.time() - t
# E1 (pairwise)
data, info, idx, sp = prepare_curve_data(tab, "E1")
T["E1_n_pairs"] = info["n_pairs_fit"]
t = time.time()
r = fit_curve_full(data, 0.04, 1.0, 1.0, 0.3, B=2, B_s=50, B_s_boot=10, B_naive=2, sp=sp)
T["E1_point_simex_plus_2boot_s"] = time.time() - t
t = time.time()
r = fit_curve_full(data, 0.04, 1.0, 1.0, 0.3, B=0, B_s=50, B_s_boot=10, B_naive=0, sp=sp)
T["E1_point_simex_s"] = time.time() - t
per_boot = (T["E1_point_simex_plus_2boot_s"] - T["E1_point_simex_s"]) / 2
T["E1_per_boot_rep_s"] = per_boot
T["E1_full_estimate_s(B=200,Bs=50,Bs_boot=10,+B_naive=1000)"] = T["E1_point_simex_s"] + 200 * per_boot + 800 * per_boot / 41
# E3
t = time.time(); run_E3(bf, panel, B=100); T["E3_B100_s"] = time.time() - t
# E4
ev = edge_view(bf, panel)
t = time.time(); c = build_cells(ev, panel.age); y, x, z, _, _ = cell_design(c); T["E4_cells_s"] = time.time() - t
T["E4_n_cells"] = int(len(y))
t = time.time(); f = fit_mixed(y, x, z, c["hog"]); T["E4_statsmodels_point_s"] = time.time() - t
t = time.time(); g = fast_mixed(y, x, z, c["hog"]); T["E4_fastreml_fit_s"] = time.time() - t
T["E4_fast_vs_statsmodels_beta_diff"] = abs(f["beta"] - g["beta"])
T["E4_est_B1000_s"] = 1000 * (T["E4_fastreml_fit_s"] + T["E4_cells_s"] * 1.5)
# E5
fe = {"mu": f["mu"], "beta": f["beta"], "gamma": f["gamma"], "x_center": 0.0, "z_center": 0.0}
t = time.time(); dat = E5Data(bf, panel, ev, fe); T["E5_setup_s"] = time.time() - t
t = time.time(); F, _ = parametric_floor(dat, P=10, rho=0.3, nu_t=0.02); T["E5_floor_P10_s"] = time.time() - t
t = time.time(); core = e5_core(dat, ev.b_aa, ev.b_t, F, B=100, B_cv=50); T["E5_core_B100_Bcv50_s"] = time.time() - t
T["E5_est_full_s(P=100,B=1000,Bcv=200,spike 4x50xB200)"] = (T["E5_floor_P10_s"] * 10 + T["E5_core_B100_Bcv50_s"] * (1000 + 200) / 150
                                                       + 4 * 50 * T["E5_core_B100_Bcv50_s"] * 200 / 150)
T["E5_n_lineages"] = len(dat.lineages); T["E5_n_contrasts"] = core["n_contrasts"]
write_json(T, a.out)
print(json.dumps(T, indent=1))
