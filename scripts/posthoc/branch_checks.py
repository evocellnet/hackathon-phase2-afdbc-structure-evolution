#!/usr/bin/env python3
"""Post-hoc verification (independent of scripts/review): branch-level checks for claims 3, 4 and the C35 sensitivity.
- Refits the production branch engine in memory (no files written to data/interim) at several tip corrections nu.
- Own code: edge-type x age-bin k (ratio of sums), terminal excess per edge, internal-only k with HOG bootstrap,
  sigma_k estimators (plan-style moment, DerSimonian-Laird, median-v, precise-HOG subset).
- Production estimator (run_E4) for the E4 point slope beta at each nu, internal-only beta, and C35 nu +/- 0.005 with B.
  CONFIRM_RUN=1 python scripts/posthoc/branch_checks.py  -> results/posthoc/branch_checks.json
"""
import sys, json, time, tempfile
from pathlib import Path
import numpy as np, polars as pl
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from strucevo.stats.driver import EndpointDriver
from strucevo.stats.branch import fit_branches
from strucevo.stats.endpoints_branch import run_E4, edge_view

OUT = ROOT / "results" / "posthoc"; OUT.mkdir(parents=True, exist_ok=True)
BINS = [0, 30, 100, 300, 1000, np.inf]
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]", *a, flush=True)

d = EndpointDriver(split="confirm", run_alt_alignments=False, outdir=tempfile.mkdtemp(), log=lambda *a: None)
d.load()
d.tau = d.tun.tau_misfit()
tab = d.tabs["tdiaa"]
eta = d.tun.eta("tdiaa")
tipvar = d.tun.tip_sd("tdiaa") ** 2
hs = {h: d.hog_species[h] for h in np.unique(tab["hog"].to_numpy()) if h in d.hog_species}

def fit(nu_t):
    return fit_branches(tab, d.panel, metrics=("aa_ml", "tdi_ml", "lddt_D"),
                        nu={"tdi_ml": nu_t, "lddt_D": eta["lddt_D"] / 2}, aln="tdiaa", hog_species=hs,
                        tau_misfit=d.tau, se_method="model", misfit=True)

def own_edges(bf):
    """my own edge table: non-excluded HOGs (any reason in aa or tdi, or misfit > tau), AA and 3Di side by side."""
    E = bf.edge_frame(d.panel); F = bf.fit_frame()
    F = F.filter(pl.col("metric").is_in(["aa_ml", "tdi_ml"]))
    bad = F.filter((pl.col("reason").fill_null("") != "") | (pl.col("misfit_ratio").fill_nan(0) > d.tau))["hog"].unique()
    a = E.filter(pl.col("metric") == "aa_ml").select("hog", "edge_id", "is_terminal", "t_child", "t_parent",
                                                    pl.col("b").alias("b_aa"), pl.col("b_raw").alias("braw_aa"))
    t = E.filter(pl.col("metric") == "tdi_ml").select("hog", "edge_id", pl.col("b").alias("b_t"), pl.col("b_raw").alias("braw_t"))
    j = a.join(t, on=["hog", "edge_id"], how="inner").filter(~pl.col("hog").is_in(bad))
    j = j.with_columns(((pl.col("t_child") + pl.col("t_parent")) / 2).alias("mid"))
    return j, F

def k_table(j, k_ref):
    rows = []
    mid = j["mid"].to_numpy(); term = j["is_terminal"].to_numpy()
    ba = j["b_aa"].to_numpy(); bt = j["b_t"].to_numpy(); brt = j["braw_t"].to_numpy()
    for lo, hi in zip(BINS[:-1], BINS[1:]):
        for et, m0 in (("terminal", term), ("internal", ~term)):
            m = m0 & (mid >= lo) & (mid < hi)
            if m.sum() < 20:
                rows.append({"bin": f"{lo}-{hi}", "type": et, "n": int(m.sum())}); continue
            rows.append({"bin": f"{lo}-{hi}", "type": et, "n": int(m.sum()), "k": float(ba[m].sum() / bt[m].sum()),
                         "excess_per_edge_vs_kint": float(np.mean(bt[m] - ba[m] / k_ref)),
                         "nu_share_of_raw_3di": float((brt[m] - bt[m]).sum() / brt[m].sum())})
    return rows

def hog_boot_k(j, B=1000, seed=7):
    h, inv = np.unique(j["hog"].to_numpy(), return_inverse=True)
    SA = np.bincount(inv, j["b_aa"].to_numpy(), len(h)); ST = np.bincount(inv, j["b_t"].to_numpy(), len(h))
    rng = np.random.default_rng(seed); W = rng.multinomial(len(h), np.full(len(h), 1 / len(h)), size=B)
    ks = (W @ SA) / (W @ ST)
    return float(SA.sum() / ST.sum()), [float(x) for x in np.percentile(ks, [1, 99])], float(np.std(np.log(ks)))

def sigma_k(j, F):
    """per-HOG k_h over all edges; v_h by delta method from the production HOG var/cov sums."""
    h, inv = np.unique(j["hog"].to_numpy(), return_inverse=True)
    SA = np.bincount(inv, j["b_aa"].to_numpy(), len(h)); ST = np.bincount(inv, j["b_t"].to_numpy(), len(h))
    fa = F.filter(pl.col("metric") == "aa_ml"); ft = F.filter(pl.col("metric") == "tdi_ml")
    va = dict(zip(fa["hog"], fa["var_sum"])); vt = dict(zip(ft["hog"], ft["var_sum"])); cv = dict(zip(ft["hog"], ft["cov_sum_aa"]))
    VA = np.array([va[x] for x in h]); VT = np.array([vt[x] for x in h]); CV = np.array([cv[x] for x in h])
    ok = (SA > 0) & (ST > 0)
    lk = np.log(SA[ok] / ST[ok]); v = (VA / SA ** 2 + VT / ST ** 2 - 2 * CV / (SA * ST))[ok]
    n = len(lk)
    out = {"n_hogs": int(n), "v_h_mean": float(v.mean()), "v_h_median": float(np.median(v)),
           "v_h_q99": float(np.quantile(v, .99)), "var_logk": float(lk.var(ddof=1))}
    out["sigma2_plan_style_mean_v"] = float(lk.var(ddof=1) - v.mean())
    out["sigma2_median_v"] = float(lk.var(ddof=1) - np.median(v))
    # robust spread: (MAD-based SD)^2 - median v
    mad = 1.4826 * np.median(np.abs(lk - np.median(lk)))
    out["sigma_robust_mad_minus_medv"] = float(np.sqrt(max(0, mad ** 2 - np.median(v))))
    # DerSimonian-Laird
    vv = np.maximum(v, 1e-8); w = 1 / vv; mu = np.sum(w * lk) / w.sum(); Qs = np.sum(w * (lk - mu) ** 2)
    tau2 = max(0.0, (Qs - (n - 1)) / (w.sum() - np.sum(w ** 2) / w.sum()))
    out["sigma_DL"] = float(np.sqrt(tau2))
    for thr in (0.02, 0.05):
        m = v < thr
        s2 = lk[m].var(ddof=1) - v[m].mean()
        out[f"precise_v_lt_{thr}"] = {"n": int(m.sum()), "sigma_k": float(np.sqrt(max(s2, 0))),
                                     "k_h_q10_q90": [float(x) for x in np.exp(np.quantile(lk[m], [.1, .9]))]}
    return out

res = {"eta_tdi": eta["tdi_ml"], "tip_sd": float(np.sqrt(tipvar)), "tau_misfit": d.tau}
nus = [0.0, 0.013, eta["tdi_ml"] / 2 - 0.005, eta["tdi_ml"] / 2, eta["tdi_ml"] / 2 + 0.005, 0.035, 0.039, 0.040, 0.045, 0.052]
res["nu_grid"] = {}
fits = {}
for nu in nus:
    bf = fit(nu); fits[nu] = bf
    j, F = own_edges(bf)
    k_all = float(j["b_aa"].sum() / j["b_t"].sum())
    ji = j.filter(~pl.col("is_terminal")); jt = j.filter(pl.col("is_terminal"))
    k_int = float(ji["b_aa"].sum() / ji["b_t"].sum()); k_term = float(jt["b_aa"].sum() / jt["b_t"].sum())
    r4 = run_E4(bf, d.panel, d.tau, B=2, sensitivities=False, tip_var=tipvar)
    rec = {"k_all": k_all, "k_internal": k_int, "k_terminal": k_term, "E4_beta_point": r4["estimate"],
           "E4_n_cells": r4["n_units"], "E4_flags": r4["flags"], "n_hogs": int(j["hog"].n_unique())}
    res["nu_grid"][str(nu)] = rec
    log("nu", nu, json.dumps(rec, default=float))
    if abs(nu - eta["tdi_ml"] / 2) < 1e-9:
        res["production_nu"] = nu
        res["k_table_prod_nu"] = k_table(j, k_int)
        res["k_internal_boot"] = hog_boot_k(ji)
        res["k_all_boot"] = hog_boot_k(j)
        res["sigma_k_all_edges"] = sigma_k(j, F)
        # strata
        st = dict(zip(d.hogs["hog"].to_list(), d.hogs["stratum"].to_list()))
        js = j.with_columns(pl.col("hog").replace_strict(st, default=None).alias("stratum"))
        res["k_by_stratum"] = {str(s): {"terminal": float(g.filter(pl.col("is_terminal"))["b_aa"].sum() / g.filter(pl.col("is_terminal"))["b_t"].sum()),
                                        "internal": float(g.filter(~pl.col("is_terminal"))["b_aa"].sum() / g.filter(~pl.col("is_terminal"))["b_t"].sum()),
                                        "all": float(g["b_aa"].sum() / g["b_t"].sum()), "n_hogs": g["hog"].n_unique()}
                               for (s,), g in js.group_by("stratum")}
        # internal-only per-HOG k spread (whole-HOG v_h -> lower bound on precision)
        hh, inv = np.unique(ji["hog"].to_numpy(), return_inverse=True)
        SA = np.bincount(inv, ji["b_aa"].to_numpy(), len(hh)); ST = np.bincount(inv, ji["b_t"].to_numpy(), len(hh))
        ok = (SA > 0.05) & (ST > 0.01)
        res["internal_kh_q10_50_90"] = [float(x) for x in np.quantile(SA[ok] / ST[ok], [.1, .5, .9])]
        log("k_table", json.dumps(res["k_table_prod_nu"], default=float))
        log("sigma", json.dumps(res["sigma_k_all_edges"], default=float))
    (OUT / "branch_checks.json").write_text(json.dumps(res, indent=1, default=float))

# internal-only E4 point (nu-independent by construction) via the production estimator, using an internal view:
# run_E4 has no internal flag for the point; reproduce through its sensitivity block with small B_sens
bfp = fits[eta["tdi_ml"] / 2]
r = run_E4(bfp, d.panel, d.tau, B=20, sensitivities=True, B_sens=20, tip_var=tipvar)
res["E4_prod_nu_reproduced"] = {"beta": r["estimate"], "internal_edges_only": r["sensitivity"].get("internal_edges_only")}
log("E4 prod", json.dumps(res["E4_prod_nu_reproduced"], default=float))
(OUT / "branch_checks.json").write_text(json.dumps(res, indent=1, default=float))
# C35: nu +/- 0.005 with bootstrap B=200
res["C35"] = {}
for nu in (eta["tdi_ml"] / 2 - 0.005, eta["tdi_ml"] / 2 + 0.005):
    bf = fits.get(nu) or fit(nu)
    r = run_E4(bf, d.panel, d.tau, B=200, sensitivities=False, tip_var=tipvar)
    res["C35"][f"{nu:.4f}"] = {"estimate": r["estimate"], "ci98": r["ci98"], "p": r["p_value"], "se": r["se_boot"], "B": 200}
    log("C35", nu, json.dumps(res["C35"][f"{nu:.4f}"], default=float))
    (OUT / "branch_checks.json").write_text(json.dumps(res, indent=1, default=float))
log("done")
