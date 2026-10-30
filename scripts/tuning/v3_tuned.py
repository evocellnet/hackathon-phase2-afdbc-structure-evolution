"""V3 with tuned noise on REAL tune HOG designs (ANALYSIS_PLAN §7 V3, §3.3; extends scripts/verify/v3_curvature.py).

Part 1 (E1 curvature). Designs = the real split == 'tune' HOGs (--aln): their species sets S_h, their primary pairs
(including which pairs failed -> weight 0), their alignment core lengths (n_core) and their AA distance profile:
true D_AA = X b_nnls (non-negative species-tree fit of the observed aa_ml), r_h = the HOG's fitted 3Di/AA ratio of
sums. Truth: D_3Di = r_h D_AA^c (c = 1: linear; power runs c in --power-c). Observation (as synth.observe_gauss):
x = D_AA + kappa_aa(bin) SE_LG e1, y = D_3Di + u_i + u_j + kappa_tdi(bin) SE_3Di e2, corr(e1, e2) = rho(bin),
u = nu + tip_sd N(0,1) per protein, aa_ml_seq = D_AA + independent noise; then the full E1 estimator (prepare_curve_data,
two-pass, SIMEX, pigeonhole bootstrap, curve_decision). Each dataset draws --hogs designs with replacement
(default 0 = every tune HOG once; use --hogs-scale 4 for a confirm-sized dataset from the production tune split).
Reports: false-concave rate (pass <= 0.10), minimal detectable |c-1| at 80% power over --power-c.

Part 2 (E4 |beta|). Worlds on the REAL panel with the real tune species sets (resampled to --beta-hogs), AA lengths
as synth.make_world, b_3Di = r_h b_AA exp(age_beta (log tau - mean)), tuned noise, branch engine, run_E4 (B =
--beta-B, sensitivities off). Reports: null SE of beta_hat, analytic minimal detectable |beta_hat| at 80% power
(two-sided 98%: (z_0.99 + z_0.80) SE), and empirical detection rates for generating age_beta in --beta-grid
(with the mean fitted beta, since the cell aggregation attenuates the generating value).
Writes ROOT/results/tuning/v3_tuned.json.
"""
import argparse
import time
import warnings

import numpy as np
from scipy.stats import norm

import tuning_common as tc
from strucevo.stats.branch import fit_branches, nnls_fit
from strucevo.stats.common import Z98
from strucevo.stats.endpoints_branch import run_E4
from strucevo.stats.endpoints_curve import prepare_curve_data, fit_curve_full, curve_decision
from strucevo.stats.synth import World, se_of, true_pair_distances, _panel_edge_lengths_for_hog

warnings.simplefilter("ignore")
ap = tc.add_common_args(argparse.ArgumentParser())
ap.add_argument("--n-datasets", type=int, default=100)
ap.add_argument("--hogs", type=int, default=0, help="HOG designs per dataset (0 = all tune HOGs once)")
ap.add_argument("--hogs-scale", type=float, default=None, help="designs per dataset = scale x n tune HOGs")
ap.add_argument("--B", type=int, default=200, help="E1 bootstrap reps (plan 200)")
ap.add_argument("--Bs", type=int, default=50, help="SIMEX B_s (plan 50)")
ap.add_argument("--Bs-boot", type=int, default=10, help="SIMEX B_s inside bootstrap (plan 10)")
ap.add_argument("--power-c", type=float, nargs="*", default=[0.90, 0.92, 0.94, 0.96])
ap.add_argument("--power-datasets", type=int, default=20)
ap.add_argument("--beta-hogs", type=int, default=0, help="HOGs per E4 world (0 = n tune HOGs)")
ap.add_argument("--beta-datasets", type=int, default=10, help="null E4 worlds (SE) and worlds per beta")
ap.add_argument("--beta-grid", type=float, nargs="*", default=[-0.05, -0.1])
ap.add_argument("--beta-B", type=int, default=100)
ap.add_argument("--tip-sd", type=float, default=None, help="override the tuned tip SD (diagnostics)")
ap.add_argument("--scalar-noise", action="store_true",
                help="diagnostic: simulate with the scalar kappa/rho summaries the estimator uses (not per bin)")
ap.add_argument("--out-name", default="v3_tuned", help="fragment name (diagnostic runs: use another name)")
ap.add_argument("--skip-e1", action="store_true")
ap.add_argument("--skip-e4", action="store_true")
a = ap.parse_args()
t0 = time.time()
root = a.root
tuned = tc.load_tuned(root, a.aln)
kaa, ktd, rh_ = tuned["kappa_aa"], tuned["kappa_tdi"], tuned["rho_tdi"]
if a.tip_sd is not None:
    tuned["tip_sd"] = tuned["summary"]["tip_sd"] = a.tip_sd
if a.scalar_noise:
    from tuning_common import BinLookup
    for k, sk in (("kappa_aa", "kappa_aa"), ("kappa_tdi", "kappa_tdi"), ("rho_tdi", "rho_tdi")):
        tuned[k] = BinLookup(None, tuned["summary"][sk], 1.0)
nu, tip_sd = tuned["nu_tdi"], tuned["tip_sd"]
bf, panel, tab, info = tc.fit_tune_branches(root, a.aln, tuned, metrics=["aa_ml", "tdi_ml"], misfit=False)
ka, kt = bf.metrics.index("aa_ml"), bf.metrics.index("tdi_ml")

# ------------------------------------------------------------------ real designs
designs = []
ncore = tab["n_core"].astype(float)
for hf in bf.hogs:
    if hf.fit["under"][ka] or hf.fit["under"][kt]:
        continue
    bn = nnls_fit(hf.X, hf.d[ka], hf.fit["w"][ka], hf.ok[ka])
    if not np.all(np.isfinite(bn)):
        continue
    Daa = hf.X @ bn
    sa, st = np.nansum(hf.fit["b"][ka]), np.nansum(hf.fit["b"][kt])
    r = st / sa if sa > 0 and st > 0 else np.nan
    L = ncore[hf.pair_index]
    L = np.where(np.isfinite(L) & (L >= 20), L, np.nanmedian(L))
    designs.append({"hog": hf.hog, "S": hf.ht.species, "ia": hf.ia, "ib": hf.ib, "Daa": Daa, "r": r,
                    "ok": hf.ok[ka] & hf.ok[kt], "L": L})
rs = np.array([d["r"] for d in designs])
r_med = float(np.nanmedian(rs))
for d in designs:
    d["r"] = float(np.clip(d["r"], 0.05, 3.0)) if np.isfinite(d["r"]) else r_med
nH = len(designs)
n_per = a.hogs if a.hogs > 0 else (int(round(a.hogs_scale * nH)) if a.hogs_scale else nH)
print(f"{nH} real designs; {n_per} HOGs per dataset; r_h median {r_med:.3f}", flush=True)


def simulate(rng, c_true=1.0):
    pick = np.arange(nH) if (n_per == nH and not a.hogs_scale and a.hogs <= 0) else rng.integers(0, nH, n_per)
    cols = {k: [] for k in ("hog", "species_a", "species_b", "status", "aa_ml", "aa_ml_se", "tdi_ml", "tdi_ml_se",
                            "aa_ml_seq", "len_a")}
    for j, k in enumerate(pick):
        d = designs[k]
        n = len(d["ia"])
        u = nu + tip_sd * rng.standard_normal(len(d["S"]))
        D3 = d["r"] * np.maximum(d["Daa"], 0.0) ** c_true
        y_true = D3 + u[d["ia"]] + u[d["ib"]]
        sx = se_of("LG", d["Daa"], d["L"]); sy = se_of("Q3Di", y_true, d["L"])
        kx = kaa(d["Daa"]); ky = ktd(d["Daa"]); rr = rh_(d["Daa"])
        z1 = rng.standard_normal(n); z2 = rr * z1 + np.sqrt(1 - rr ** 2) * rng.standard_normal(n)
        cols["hog"].append(np.full(n, f"{d['hog']}#{j}")); cols["species_a"].append(d["S"][d["ia"]])
        cols["species_b"].append(d["S"][d["ib"]]); cols["status"].append(np.where(d["ok"], "ok", "fail_cov"))
        cols["aa_ml"].append(d["Daa"] + kx * sx * z1); cols["aa_ml_se"].append(sx)
        cols["tdi_ml"].append(y_true + ky * sy * z2); cols["tdi_ml_se"].append(sy)
        cols["aa_ml_seq"].append(d["Daa"] + kx * sx * rng.standard_normal(n)); cols["len_a"].append(d["L"])
    t = {k: np.concatenate(v) for k, v in cols.items()}
    n = len(t["hog"])
    t.update({"split": np.full(n, "tune"), "aa_sat": np.zeros(n, bool), "tdi_sat": np.zeros(n, bool),
              "one_to_one": np.ones(n, bool)})
    return t


kx_s, ky_s, rho_s = tuned["summary"]["kappa_aa"], tuned["summary"]["kappa_tdi"], tuned["summary"]["rho_tdi"]


def one(i, c_true=1.0, tag=0):
    rng = tc.tuning_rng("v3_tuned", tag, i)
    t = simulate(rng, c_true)
    data, inf, idx, sp = prepare_curve_data(t, "E1")
    res = fit_curve_full(data, 2 * nu, kx_s, ky_s, rho_s, B=a.B, B_s=a.Bs, B_s_boot=a.Bs_boot, B_naive=a.B,
                         seed_stage="v3_sim", sp=sp, seed_extra=(500 + tag, i))
    cls, p, ex = curve_decision(res["c"], res["se_boot"])
    return {"c_naive": res["c_naive"], "c_simex": res["c"], "se": res["se_boot"], "class": cls,
            "n_pairs": inf["n_pairs_fit"], "n_hogs": inf["n_hogs"]}


rec = tc.header(root, "v3_tuned", a)
rec["noise"] = {"summary": tuned["summary"], "provisional": tuned["provisional"], "sources": tuned["sources"],
                "kappa_rho_in_estimator": "scalar summaries (median over bins); simulation uses the per-bin values"}
rec["designs"] = {"n_real_hogs": nH, "hogs_per_dataset": n_per, "r_h_median": r_med,
                  "D_aa_quantiles": np.quantile(np.concatenate([d["Daa"] for d in designs]), [0.05, 0.5, 0.95]).tolist()}
if not a.skip_e1:
    t1 = time.time()
    runs = [one(i) for i in range(a.n_datasets)]
    t_e1 = time.time() - t1
    cls = [r["class"] for r in runs]
    cs = np.array([r["c_simex"] for r in runs]); cn = np.array([r["c_naive"] for r in runs])
    se = np.array([r["se"] for r in runs])
    e1 = {"n_datasets": a.n_datasets, "false_concave_rate": float(np.mean([c == "concave" for c in cls])),
          "class_counts": {k: cls.count(k) for k in set(cls)}, "mean_c_naive": float(cn.mean()),
          "mean_c_simex": float(cs.mean()), "sd_c_simex": float(cs.std(ddof=1)) if len(cs) > 1 else None,
          "mean_se_boot": float(se.mean()), "coverage98_of_1": float(np.mean(np.abs(cs - 1) <= Z98 * se)),
          "mean_n_pairs": float(np.mean([r["n_pairs"] for r in runs])),
          "mean_n_hogs": float(np.mean([r["n_hogs"] for r in runs])), "runtime_s": t_e1,
          "per_dataset_s": t_e1 / max(1, a.n_datasets)}
    power = {}
    for j, c in enumerate(a.power_c):
        rr = [one(i, c_true=c, tag=10 + j) for i in range(a.power_datasets)]
        power[str(c)] = {"concave_rate": float(np.mean([r["class"] == "concave" for r in rr])),
                         "mean_c_simex": float(np.mean([r["c_simex"] for r in rr])),
                         "mean_c_naive": float(np.mean([r["c_naive"] for r in rr])),
                         "mean_se": float(np.mean([r["se"] for r in rr]))}
        print(f"  c={c}: concave rate {power[str(c)]['concave_rate']:.2f}", flush=True)
    e1["power"] = power
    det = [abs(1 - float(c)) for c, v in power.items() if v["concave_rate"] >= 0.80]
    e1["min_detectable_abs_c_minus_1"] = min(det) if det else None
    e1["min_detectable_note"] = (None if det else f"no c in {a.power_c} reached 80% power")
    # analytic companion: |c-1| detectable at 80% with one-sided 99% bound: (z_.99 + z_.80) * mean SE
    e1["analytic_mde_abs_c_minus_1"] = float((Z98 + norm.ppf(0.8)) * se.mean())
    e1["pass"] = e1["false_concave_rate"] <= 0.10
    rec["E1"] = e1
    print("E1:", {k: e1[k] for k in ("false_concave_rate", "mean_c_simex", "mean_se_boot", "min_detectable_abs_c_minus_1",
                                     "analytic_mde_abs_c_minus_1")}, flush=True)

# ------------------------------------------------------------------ E4 |beta|
if not a.skip_e4:
    t1 = time.time()
    Sets = [d["S"] for d in designs]
    nb = a.beta_hogs if a.beta_hogs > 0 else nH
    cand_root_age = None

    def world(rng, age_beta):
        N = panel.n_nodes
        lmult = np.exp(0.3 * rng.standard_normal(N) - 0.045)
        pick = rng.integers(0, len(Sets), nb) if nb != len(Sets) else np.arange(len(Sets))
        hogs, species, baa, rh = [], {}, {}, {}
        for j, k in enumerate(pick):
            S = np.sort(Sets[k])
            ht0 = __import__("strucevo.stats.trees", fromlist=["prune"]).prune(panel, S)
            rate = rng.uniform(0.5, 1.9) / (2 * panel.age[ht0.root])
            _, b = _panel_edge_lengths_for_hog(panel, S, rate, lmult, 0.2, rng)
            hid = f"W{j:05d}"
            hogs.append(hid); species[hid] = S; baa[hid] = b; rh[hid] = 0.35 * np.exp(0.3 * rng.standard_normal())
        w = World(panel=panel, hogs=hogs, species=species, b_aa=baa, r_h=rh, lineage_mult=lmult,
                  spiked=np.ones(N), age_beta=age_beta)
        for hid in hogs:
            w.pairs_true[hid] = true_pair_distances(w, hid)
        return w

    def e4_once(age_beta, i, tag):
        rng = tc.tuning_rng("v3_beta", tag, i)
        w = world(rng, age_beta)
        t = tc.observe_world_tuned(w, tuned, rng)
        bfw = fit_branches(t, panel, metrics=["aa_ml", "tdi_ml"], nu={"tdi_ml": nu}, se_method="model", misfit=False)
        r = run_E4(bfw, panel, B=a.beta_B, sensitivities=False)
        return r["estimate"], r["se_boot"], r["p_value"], r["n_units"]

    nulls = [e4_once(0.0, i, 0) for i in range(a.beta_datasets)]
    se0 = float(np.mean([s for _, s, _, _ in nulls]))
    e4 = {"hogs_per_world": nb, "B": a.beta_B, "n_null_worlds": a.beta_datasets,
          "null_mean_beta": float(np.mean([b for b, *_ in nulls])), "null_mean_se": se0,
          "null_sd_beta": float(np.std([b for b, *_ in nulls], ddof=1)) if len(nulls) > 1 else None,
          "null_false_positive_rate_98": float(np.mean([p < 0.02 for *_, p, _ in nulls])),
          "mean_n_cells": float(np.mean([n for *_, n in nulls])),
          "analytic_mde_abs_beta_hat": float((Z98 + norm.ppf(0.8)) * se0)}
    grid = {}
    for j, bt in enumerate(a.beta_grid):
        rr = [e4_once(bt, i, 10 + j) for i in range(a.beta_datasets)]
        grid[str(bt)] = {"detect_rate_98": float(np.mean([p < 0.02 for *_, p, _ in rr])),
                         "mean_beta_hat": float(np.mean([b for b, *_ in rr])),
                         "mean_se": float(np.mean([s for _, s, _, _ in rr]))}
        print(f"  age_beta={bt}: detect {grid[str(bt)]['detect_rate_98']:.2f} beta_hat {grid[str(bt)]['mean_beta_hat']:.4f}", flush=True)
    e4["grid_generating_beta"] = grid
    det = [abs(v["mean_beta_hat"]) for v in grid.values() if v["detect_rate_98"] >= 0.80]
    e4["min_detectable_abs_beta_hat_empirical"] = min(det) if det else None
    e4["runtime_s"] = time.time() - t1
    rec["E4"] = e4
    print("E4:", {k: e4[k] for k in ("null_mean_se", "analytic_mde_abs_beta_hat", "min_detectable_abs_beta_hat_empirical")}, flush=True)

e1r = rec.get("E1", {})
mdc = e1r.get("min_detectable_abs_c_minus_1")
rec["min_detectable_c_basis"] = ("empirical (smallest |c-1| in --power-c with >= 80% concave calls)" if mdc is not None
                                 else "analytic (z_.99 + z_.80) x mean SE_boot (no grid value reached 80%)")
rec["min_detectable_beta_basis"] = ("empirical mean |beta_hat| at the smallest grid beta with >= 80% detection"
                                    if rec.get("E4", {}).get("min_detectable_abs_beta_hat_empirical") is not None
                                    else "analytic (z_.99 + z_.80) x null SE_boot of beta_hat")
rec["yaml"] = {"false_concave_rate": e1r.get("false_concave_rate"),
               "min_detectable_c": mdc if mdc is not None else e1r.get("analytic_mde_abs_c_minus_1"),
               "min_detectable_beta_e4": (rec.get("E4", {}).get("min_detectable_abs_beta_hat_empirical")
                                          or rec.get("E4", {}).get("analytic_mde_abs_beta_hat"))}
flags = []
if "E1" in rec:
    e1 = rec["E1"]
    bias = e1["mean_c_simex"] - 1
    if e1["false_concave_rate"] > 0.10:
        flags.append("V3_FAIL_false_concave_gt_0.10")
    if e1["n_datasets"] > 1 and abs(bias) > 2 * e1["mean_se_boot"] / max(1, e1["n_datasets"]) ** 0.5:
        flags.append(f"E1_simex_bias_under_linear_truth:{bias:+.3f}")
if "E4" in rec:
    e4 = rec["E4"]
    if e4["n_null_worlds"] > 1 and abs(e4["null_mean_beta"]) > 2 * e4["null_mean_se"] / e4["n_null_worlds"] ** 0.5:
        flags.append(f"E4_null_bias:{e4['null_mean_beta']:+.3f}(se {e4['null_mean_se']:.3f})")
    if e4["null_false_positive_rate_98"] > 0.10:
        flags.append(f"E4_null_false_positive_rate:{e4['null_false_positive_rate_98']:.2f}")
rec["flags"] = flags
print("flags:", flags)
rec["fit_info"] = info
rec["runtime_s"] = time.time() - t0
tc.dump(rec, tc.out_path(root, a.out_name))
print(rec["yaml"], "runtime", rec["runtime_s"])
