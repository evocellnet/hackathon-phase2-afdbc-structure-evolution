#!/usr/bin/env python3
"""Post-hoc: between-HOG spread of log k_h on INTERNAL edges only, with internal-only sampling variances from the
production branch fits (Var of sum of internal b = 1' V_b 1 per metric, V_b = A diag(1/w) A'). The AA/3Di covariance is
bracketed (HC1-score cross-covariance as in production, and 0). Estimators: DL, median-v moment, precise-HOG subset.
  CONFIRM_RUN=1 python scripts/posthoc/sigma_k_internal.py -> results/posthoc/sigma_k_internal.json"""
import sys, json, tempfile
from pathlib import Path
import numpy as np, polars as pl
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from strucevo.stats.driver import EndpointDriver

d = EndpointDriver(split="confirm", run_alt_alignments=False, outdir=tempfile.mkdtemp(), log=lambda *a: None)
d.load(); d.tau = d.tun.tau_misfit()
bf = d.fit(d.tabs["tdiaa"], "tdiaa")
F = bf.fit_frame().filter(pl.col("metric").is_in(["aa_ml", "tdi_ml"]))
bad = set(F.filter((pl.col("reason").fill_null("") != "") | (pl.col("misfit_ratio").fill_nan(0) > d.tau))["hog"].to_list())
ka, kt = bf.metrics.index("aa_ml"), bf.metrics.index("tdi_ml")
rows = []
for hf in bf.hogs:
    if hf.hog in bad or hf.fit["under"][ka] or hf.fit["under"][kt]:
        continue
    ht = hf.ht; term = np.asarray(ht.is_terminal); m = ~term
    if m.sum() < 2:
        continue
    cvec = np.zeros(ht.n_params); np.add.at(cvec, ht.param_of_edge[m], ht.share[m])
    fit = hf.fit; out = {}
    for k, name in ((ka, "aa"), (kt, "t")):
        out["S" + name] = float(cvec @ fit["b"][k]); out["V" + name] = float(cvec @ fit["V"][k] @ cvec)
    G1 = fit["A"][ka] * fit["r"][ka][None] * np.sqrt(fit["hc1"][ka]); G2 = fit["A"][kt] * fit["r"][kt][None] * np.sqrt(fit["hc1"][kt])
    out["C"] = float(cvec @ (G1 @ G2.T) @ cvec)
    rows.append(out | {"hog": hf.hog, "n_int": int(m.sum())})
df = pl.DataFrame(rows)
SA, ST, VA, VT, CC = (df[c].to_numpy() for c in ("Saa", "St", "Vaa", "Vt", "C"))
ok = (SA > 0) & (ST > 0)
res = {"n_hogs_internal": int(ok.sum()), "k_internal_ratio_of_sums_subset": float(SA[ok].sum() / ST[ok].sum())}
for rho in ("hc1_cross", 0.0):
    lk = np.log(SA[ok] / ST[ok])
    v = VA[ok] / SA[ok] ** 2 + VT[ok] / ST[ok] ** 2 - 2 * (CC[ok] if rho == "hc1_cross" else 0.0) / (SA[ok] * ST[ok])
    n = len(lk); w = 1 / np.maximum(v, 1e-8); mu = np.sum(w * lk) / w.sum(); Q = np.sum(w * (lk - mu) ** 2)
    tau2 = max(0.0, (Q - (n - 1)) / (w.sum() - np.sum(w ** 2) / w.sum()))
    r = {"sigma_DL": float(np.sqrt(tau2)), "v_median": float(np.median(v)), "v_mean": float(v.mean()),
         "sigma_median_v": float(np.sqrt(max(0, lk.var(ddof=1) - np.median(v)))),
         "sigma_plan_mean_v": float(np.sqrt(max(0, lk.var(ddof=1) - v.mean()))), "sigma2_plan_mean_v": float(lk.var(ddof=1) - v.mean())}
    for thr in (0.01, 0.02, 0.05):
        mm = v < thr
        r[f"precise_v_lt_{thr}"] = {"n": int(mm.sum()), "sigma": float(np.sqrt(max(0, lk[mm].var(ddof=1) - v[mm].mean()))),
                                   "k_h_q10_q90": np.exp(np.quantile(lk[mm], [.1, .9])).tolist()}
    res[f"rho_{rho}"] = r
(ROOT / "results/posthoc/sigma_k_internal.json").write_text(json.dumps(res, indent=1))
print(json.dumps(res, indent=1))
