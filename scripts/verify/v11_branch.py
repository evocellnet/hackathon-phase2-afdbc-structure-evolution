"""V11: independent E3 and E4 point estimates from the confirm branch table (§2.5, §2.6, §4.3, §4.4, D32)."""
import sys, warnings
import numpy as np
from v11_common import *

br = load_branches(PRIMARY, SPLIT)
fit = pl.read_parquet(ROOT / f"data/interim/branches/{PRIMARY}_fit.parquet")
assert set(br["split"].unique()) == {"confirm"} and set(fit["split"].unique()) == {"confirm"}
hc = hogcols()
TAU = TUN["tau_misfit"]

# ---- HOG inclusion: no exclusion reason in the AA or 3Di fit, misfit_ratio <= tau (§4.3 "excluded for misfit or failures")
f2 = fit.filter(pl.col("metric").is_in(["aa_ml", "tdi_ml"]))
bad = f2.filter(pl.col("excluded") | (pl.col("misfit_ratio").fill_nan(None) > TAU).fill_null(False))["hog"].unique().to_list()
reasons = {}
for h, r in f2.filter(pl.col("hog").is_in(bad)).group_by("hog").agg(pl.col("reason").str.join(";")).iter_rows():
    for x in set(t for t in r.split(";") if t):
        reasons[x] = reasons.get(x, 0) + 1
# misfit from hogs_endpoint_cols as cross-check
mis_hc = set(hc.filter(pl.col("excl_misfit"))["hog"].to_list())
mis_fit = set(fit.filter((pl.col("metric") == "aa_ml") & (pl.col("misfit_ratio").fill_nan(None) > TAU).fill_null(False))["hog"].to_list())
keep = sorted(set(fit["hog"].unique().to_list()) - set(bad))
print("HOGs kept", len(keep), "excl reasons", reasons, "misfit sets equal", mis_hc == mis_fit)

E = br.filter(pl.col("hog").is_in(keep))
A = E.filter(pl.col("metric") == "aa_ml").sort("hog", "edge_id")
T = E.filter(pl.col("metric") == "tdi_ml").sort("hog", "edge_id")
assert (A["hog"] == T["hog"]).all() and (A["edge_id"] == T["edge_id"]).all()
bA = A["b"].to_numpy(); bT = T["b"].to_numpy()
assert np.isfinite(bA).all() and np.isfinite(bT).all()
# tip correction check: b = b_raw - nu on terminal 3Di edges, nu = eta_tdi/2
nu = TUN["eta"][PRIMARY]["tdi_ml"] / 2
term = T["is_terminal"].to_numpy()
tipchk = float(np.max(np.abs(T["b"].to_numpy() - (T["b_raw"].to_numpy() - nu * term))))
tipchkA = float(np.max(np.abs(A["b"].to_numpy() - A["b_raw"].to_numpy())))
res = {"n_hogs_own": len(keep), "n_edges_own": A.height, "tip_correction_maxabs": tipchk, "aa_b_eq_braw_maxabs": tipchkA,
       "misfit_sets_equal_hogcols": mis_hc == mis_fit}

# ================= E3
p3 = endpoint("E3")
k = bA.sum() / bT.sum()
res["E3"] = {"k_own": k, "k_project": p3["estimate"], "absdiff": abs(k - p3["estimate"]),
             "n_units_project": p3["n_units"], "n_hogs_project": p3["n_hogs"], "n_excluded_project": p3["n_excluded"],
             "n_excluded_own": reasons}
inte = ~A["is_terminal"].to_numpy()
ki = bA[inte].sum() / bT[inte].sum()
res["E3"]["internal_only_own"] = ki; res["E3"]["internal_only_project"] = p3["sensitivity"]["internal_edges_only"]["estimate"]
# per-HOG k_h and sigma2_k (delta method). v_h with full HOG-total (co)variances from the fit table (var_sum, cov_sum_aa)
hogA = A["hog"].to_numpy()
SA = pl.DataFrame({"hog": hogA, "a": bA, "t": bT, "va": A["b_se"].to_numpy() ** 2, "vt": T["b_se"].to_numpy() ** 2,
                   "c": T["cov_aa"].to_numpy()}).group_by("hog").agg(pl.all().sum()).sort("hog")
fa = fit.filter(pl.col("metric") == "aa_ml").select("hog", pl.col("var_sum").alias("VA"))
ft = fit.filter(pl.col("metric") == "tdi_ml").select("hog", pl.col("var_sum").alias("VT"), pl.col("cov_sum_aa").alias("CT"))
SA = SA.join(fa, on="hog").join(ft, on="hog")
a_, t_ = SA["a"].to_numpy(), SA["t"].to_numpy()
pos = (a_ > 0) & (t_ > 0)
lk = np.log(a_[pos] / t_[pos])
for nm, va, vt, cv in (("diag", SA["va"], SA["vt"], SA["c"]), ("fulltotals", SA["VA"], SA["VT"], SA["CT"])):
    va, vt, cv = va.to_numpy()[pos], vt.to_numpy()[pos], cv.to_numpy()[pos]
    vh = va / a_[pos] ** 2 + vt / t_[pos] ** 2 - 2 * cv / (a_[pos] * t_[pos])
    for dd in (0, 1):
        res["E3"][f"sigma2_k_own_{nm}_ddof{dd}"] = float(np.var(lk, ddof=dd) - np.mean(vh))
res["E3"]["n_hogs_kh_positive"] = int(pos.sum())
res["E3"]["sigma2_k_project"] = p3["sigma2_k"]
# CATH contrast (C24 class from hogs_endpoint_cols)
cl = dict(zip(hc["hog"].to_list(), hc["cath_class"].to_list()))
cls = np.array([cl.get(h, "none") for h in hogA])
ka = bA[cls == "alpha"].sum() / bT[cls == "alpha"].sum(); kb = bA[cls == "beta"].sum() / bT[cls == "beta"].sum()
cc = p3["cath_contrast"]
res["E3"]["cath"] = {"k_alpha_own": ka, "k_beta_own": kb, "delta_own": float(np.log(ka) - np.log(kb)),
                     "delta_project": cc["delta"], "absdiff": abs(np.log(ka) - np.log(kb) - cc["delta"]),
                     "n_alpha_own": len(set(hogA[cls == "alpha"])), "n_beta_own": len(set(hogA[cls == "beta"])),
                     "n_alpha_project": cc["n_alpha"], "n_beta_project": cc["n_beta"]}
print(res["E3"], flush=True)

# ================= E4 cells
p4 = endpoint("E4")
BINS = np.array([0, 30, 100, 300, 1000, np.inf])
tau = A["tau"].to_numpy(); tc = A["t_child"].to_numpy(); tp = A["t_parent"].to_numpy()
res["tau_eq_midpoint_maxabs"] = float(np.nanmax(np.abs(tau - (tc + tp) / 2)))
res["tau_nan"] = int(np.isnan(tau).sum())
gbin = np.searchsorted(BINS, tau, side="right") - 1
r_bar = bT.sum() / bA.sum()                 # D32: global r-bar
tip_var = p4["tip_var"]                     # per-protein tip-noise variance used by the project (check below)
vt_edge = T["b_se"].to_numpy() ** 2 + tip_var * term


def cells(rule="info", internal_only=False):
    rows = []
    sel_all = np.ones(len(bA), bool) if not internal_only else ~term
    df = pl.DataFrame({"hog": hogA, "g": gbin, "a": bA, "t": bT, "tau": tau, "v": vt_edge, "sel": sel_all}).filter(pl.col("sel"))
    for (h,), sub in df.group_by(["hog"], maintain_order=True):
        cs = {}
        for g, a, t, ta, v in sub.select("g", "a", "t", "tau", "v").iter_rows():
            c = cs.setdefault(g, [0.0, 0.0, [], 0.0, 0])
            c[0] += a; c[1] += t; c[2].append(ta); c[3] += v; c[4] += 1
        def weak(c):
            if c[0] < 0.05:
                return True
            return rule == "info" and r_bar * c[0] < 4 * np.sqrt(c[3])
        while len(cs) > 1:
            gs = sorted(cs)
            w = [g for g in gs if weak(cs[g])]
            if not w:
                break
            g = w[0]                                   # youngest weak cell first
            tgt = gs[gs.index(g) + 1] if g != gs[-1] else gs[-2]   # next older; oldest -> next younger
            c, d = cs.pop(g), cs[tgt]
            d[0] += c[0]; d[1] += c[1]; d[2] += c[2]; d[3] += c[3]; d[4] += c[4]
        for g, c in cs.items():
            rows.append((h, g, c[0], c[1], float(np.mean(c[2])), c[4]))
    return pl.DataFrame(rows, schema=["hog", "g", "SA", "ST", "tau_mean", "n_edges"], orient="row")


def mixed(C):
    import statsmodels.formula.api as smf
    SA_, ST_ = C["SA"].to_numpy(), C["ST"].to_numpy()
    cens = ST_ <= 0
    ST_ = np.where(cens, 0.005, ST_)
    y = np.log(ST_ / SA_); x = np.log(C["tau_mean"].to_numpy()); z = np.log(SA_ / C["n_edges"].to_numpy())
    xc, zc = x.mean(), z.mean()
    D = pd.DataFrame({"y": y, "x": x - xc, "z": z - zc, "hog": C["hog"].to_list()})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        md = smf.mixedlm("y ~ x + z", D, groups=D["hog"], re_formula="~x")
        best = None
        for meth in ("lbfgs", "bfgs", "powell", "nm"):
            try:
                r = md.fit(reml=True, method=meth, maxiter=5000)
                if best is None or r.llf > best.llf:
                    best = r
            except Exception as e:
                pass
    return {"beta": float(best.fe_params["x"]), "mu": float(best.fe_params["Intercept"]), "gamma": float(best.fe_params["z"]),
            "llf": float(best.llf), "converged": bool(best.converged), "x_center": xc, "z_center": zc,
            "n_cells": len(y), "censored_frac": float(cens.mean()), "SA_nonpos": int((SA_ <= 0).sum())}


import pandas as pd
C = cells("info")
m = mixed(C)
res["E4"] = {"own": m, "beta_project": p4["estimate"], "absdiff": abs(m["beta"] - p4["estimate"]),
             "n_cells_project": p4["n_units"], "censored_frac_project": p4["censored_frac"],
             "mu_project": p4["mu"], "gamma_project": p4["gamma"], "x_center_project": p4["x_center"], "z_center_project": p4["z_center"]}
print(res["E4"], flush=True)
Cp = cells("plan")
mp = mixed(Cp)
res["E4"]["cells_plan_rule_own"] = mp; res["E4"]["cells_plan_rule_project"] = p4["sensitivity"]["cells_plan_rule"]
print("plan rule", mp, p4["sensitivity"]["cells_plan_rule"]["estimate"], flush=True)
wjson(res, OUT / "V11_branch.json")
