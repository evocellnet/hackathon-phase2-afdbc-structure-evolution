"""§4.5 E5 spike-in power on the tuning set (D8 lineage precision rule, min_t_snr = 8 by default).

(A) REAL tune: branch engine on split == 'tune' (--aln), E4 fixed effects (run_E4 point, small B), E5 lineages
    (>= --min-hogs HOGs, D_e >= 0.05, D8 SNR rule). If fewer than --min-lineages eligible lineages, the real-data
    power is NOT computed and the record is flagged 'insufficient_lineages'. Otherwise: parametric floor F (P),
    c_v/vx_v (B_cv), spike-in power for f in --fs with --reps seeded replicates (B per replicate).
(B) SYNTHETIC world sized like the production tune split (--syn-hogs, default 1400 HOGs, 8-30 species) on the REAL
    panel tree, N0 truth (b_3Di = r_h b_AA), observed with the tuned noise (kappa, rho, nu, tip SD; provisional
    where the SE-calibration fragments are absent); then the same E5 pipeline and the same spike-in procedure.
The yaml values (e5.power_*) are taken from (A) when it is sufficient, else from (B) (recorded in `power_source`).
Writes ROOT/results/tuning/e5_power.json.
"""
import argparse
import time

import numpy as np

import tuning_common as tc
from strucevo.stats.branch import fit_branches
from strucevo.stats.common import stage_rng
from strucevo.stats.endpoints_branch import run_E4, edge_view
from strucevo.stats.endpoints_lineage import E5Data, parametric_floor, e5_core, spikein_power
from strucevo.stats.synth import make_world

ap = tc.add_common_args(argparse.ArgumentParser())
ap.add_argument("--min-t-snr", type=float, default=8.0)
ap.add_argument("--min-hogs", type=int, default=50)
ap.add_argument("--min-lineages", type=int, default=5)
ap.add_argument("--fs", type=float, nargs="+", default=[1.5, 1.25, 2.0, 1.0])
ap.add_argument("--reps", type=int, default=50)
ap.add_argument("--B", type=int, default=200, help="bootstrap replicates per spike-in replicate (plan 200)")
ap.add_argument("--P", type=int, default=100, help="parametric floor replicates (plan 100)")
ap.add_argument("--B-cv", type=int, default=200)
ap.add_argument("--e4-B", type=int, default=20, help="E4 bootstrap reps for the fixed effects (point only matters)")
ap.add_argument("--syn-hogs", type=int, default=1400)
ap.add_argument("--skip-synthetic", action="store_true")
a = ap.parse_args()
t0 = time.time()
root = a.root
tuned = tc.load_tuned(root, a.aln)
tau = tc.tau_from_fragment(root)
rho = tuned["summary"]["rho_tdi"]


def e5_power_block(bf, panel, label, stage_off):
    t = time.time()
    e4 = run_E4(bf, panel, tau_misfit=tau, B=a.e4_B, sensitivities=False, tip_var=tuned["tip_sd"] ** 2)
    fe = {"mu": e4["mu"], "gamma": e4["gamma"], "x_center": e4["x_center"], "z_center": e4["z_center"],
          "beta": e4["estimate"]}
    ev = edge_view(bf, panel, tau_misfit=tau)
    dat = E5Data(bf, panel, ev, fe, None, min_hogs=a.min_hogs, min_hogs_contrast=a.min_hogs, min_t_snr=a.min_t_snr)
    dat_nosnr = E5Data(bf, panel, ev, fe, None, min_hogs=a.min_hogs, min_hogs_contrast=a.min_hogs, min_t_snr=None)
    out = {"label": label, "n_hogs": int(dat.H), "n_lineages": int(len(dat.lineages)),
           "n_lineages_without_D8": int(len(dat_nosnr.lineages)), "n_contrasts_kept": int(len(dat.c_keep)),
           "e4_beta": e4["estimate"], "e4_flags": e4["flags"]}
    if len(dat.lineages) < a.min_lineages:
        out["status"] = "insufficient_lineages"
        out["message"] = (f"{len(dat.lineages)} eligible lineages (need >= {a.min_lineages}; rule >= {a.min_hogs} HOGs "
                          f"per lineage, D_e >= 0.05, D8 SNR >= {a.min_t_snr}); spike-in power not computed")
        out["runtime_s"] = time.time() - t
        return out
    F, Fv = parametric_floor(dat, rho=rho, P=a.P, nu_t=tuned["nu_tdi"], rng=stage_rng("bootstrap_E5", 1, stage_off))
    core = e5_core(dat, ev.b_aa, ev.b_t, F, B=a.B, B_cv=a.B_cv, rng=stage_rng("bootstrap_E5", 0, stage_off))
    out.update({"floor_F": F, "sigma2_obs": core["sigma2_obs"], "sigma2_L": core["sigma2_L"],
                "sigma_L_upper99": core["sigma_L_upper99"], "beta_contrast": core["beta"],
                "beta_ci98": core["beta_ci98"], "beta_ci98_hog_only": core["beta_ci98_hog_only"],
                "n_contrasts": core["n_contrasts"]})
    power = {}
    for i, f in enumerate(a.fs):
        ts = time.time()
        pw, ac, det = spikein_power(dat, F, core["cv"], f=f, reps=a.reps, B=a.B, stage_extra=100 * stage_off + i)
        power[f"f={f}"] = {"detect_rate": pw, "accept_rate": ac, "reps": a.reps, "runtime_s": time.time() - ts}
        print(f"  [{label}] f={f}: detect {pw:.2f} accept {ac:.2f}", flush=True)
    out["power"] = power
    out["status"] = "ok"
    out["runtime_s"] = time.time() - t
    return out


rec = tc.header(root, "e5_power", a)
rec["tau_misfit_used"] = tau
rec["noise"] = {"summary": tuned["summary"], "provisional": tuned["provisional"], "sources": tuned["sources"]}
bf, panel, tab, info = tc.fit_tune_branches(root, a.aln, tuned, tau_misfit=tau, metrics=["aa_ml", "tdi_ml"])
real = e5_power_block(bf, panel, "real_tune", 0)
real["fit_info"] = info
rec["real_tune"] = real
print("real tune:", real["status"], real["n_lineages"], "lineages", flush=True)
if not a.skip_synthetic:
    t = time.time()
    w = make_world(panel, a.syn_hogs, rng=tc.tuning_rng("e5_world", 0), min_sp=8, max_sp=30)
    tabs = tc.observe_world_tuned(w, tuned, tc.tuning_rng("e5_world", 1))
    bfs = fit_branches(tabs, panel, metrics=["aa_ml", "tdi_ml"], nu={"tdi_ml": tuned["nu_tdi"]}, se_method="model")
    syn = e5_power_block(bfs, panel, "synthetic_prod_tune_size", 1)
    syn["world"] = {"n_hogs": a.syn_hogs, "n_pairs": int(len(tabs["hog"])), "panel": "real panel tree",
                    "truth": "N0 (b_3Di = r_h b_AA), lineage AA-rate SD 0.3, HOG-edge SD 0.2 (synth.make_world defaults)",
                    "build_fit_s": time.time() - t - syn["runtime_s"]}
    rec["synthetic"] = syn
src = rec["real_tune"] if rec["real_tune"]["status"] == "ok" else rec.get("synthetic")
rec["power_source"] = (None if src is None or src.get("status") != "ok" else src["label"])
if rec["power_source"]:
    p = src["power"]
    g = lambda f, k="detect_rate": p.get(f"f={f}", {}).get(k)
    rec["yaml"] = {"power_f1.5": g(1.5), "power_f1.25": g(1.25), "power_f2": g(2.0), "null_rate_f1": g(1.0),
                   "p_accept_f1.5": g(1.5, "accept_rate"), "min_t_snr": a.min_t_snr}
flags = []
if rec["real_tune"]["status"] != "ok":
    flags.append("insufficient_lineages_real_tune")
if tuned["provisional"]:
    flags.append("provisional_noise:" + ",".join(tuned["provisional"]))
rec["flags"] = flags
rec["runtime_s"] = time.time() - t0
tc.dump(rec, tc.out_path(root, "e5_power"))
print(rec.get("yaml"), rec["flags"], "runtime", rec["runtime_s"])
