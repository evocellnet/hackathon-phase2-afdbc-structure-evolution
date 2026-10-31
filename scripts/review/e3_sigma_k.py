#!/usr/bin/env python3
"""Statistical review (read-only): why sigma_k = 0 in E3? Distribution of per-HOG v_h, robust / precision-weighted
between-HOG variance of log k_h (DerSimonian-Laird, Paule-Mandel-like REML-free moment, trimmed), all edges vs
internal edges only, and with an extra terminal 3Di excess removed.
  CONFIRM_RUN=1 python scripts/review/e3_sigma_k.py > logs/review/e3_sigma_k.json
"""
import sys, json
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts" / "review"))
from e4_tip_diagnostics import load_bf, with_nu
from strucevo.stats.driver import Tuning
from strucevo.stats.endpoints_branch import edge_view

def stats(ev):
    H = len(ev.hog_ids)
    SA = np.bincount(ev.hog_code, ev.b_aa, H); ST = np.bincount(ev.hog_code, ev.b_t, H)
    ok = (SA > 0) & (ST > 0)
    lk = np.log(SA[ok] / ST[ok])
    v = (ev.hog_var_aa / SA ** 2 + ev.hog_var_t / ST ** 2 - 2 * ev.hog_cov / (SA * ST))[ok]
    o = {"n": int(ok.sum()), "n_nonpos_ST": int((ST <= 0).sum()), "var_logk": float(np.var(lk, ddof=1)),
         "mean_v": float(np.mean(v)), "median_v": float(np.median(v)), "v_q": np.percentile(v, [10, 50, 90, 99]).tolist(),
         "plan_sigma2": float(np.var(lk, ddof=1) - np.mean(v))}
    # DerSimonian-Laird
    w = 1 / np.maximum(v, 1e-8); mu = np.sum(w * lk) / w.sum()
    Qs = np.sum(w * (lk - mu) ** 2); df = len(lk) - 1
    tau2 = max(0.0, (Qs - df) / (w.sum() - np.sum(w ** 2) / w.sum()))
    o["DL_tau2"] = tau2; o["DL_tau"] = float(np.sqrt(tau2)); o["DL_Q_over_df"] = float(Qs / df)
    # precise subset: v_h < 0.02 (SE of log k < 0.14)
    for thr in (0.01, 0.02, 0.05):
        m = v < thr
        if m.sum() > 20:
            o[f"precise_v<{thr}"] = {"n": int(m.sum()), "var_logk": float(np.var(lk[m], ddof=1)), "mean_v": float(v[m].mean()),
                                     "sigma_k": float(np.sqrt(max(0, np.var(lk[m], ddof=1) - v[m].mean()))),
                                     "q10_q90_k": np.exp(np.percentile(lk[m], [10, 90])).tolist()}
    # robust: MAD-based variance minus median v
    mad = 1.4826 * np.median(np.abs(lk - np.median(lk)))
    o["robust_sigma2_mad_minus_median_v"] = float(mad ** 2 - np.median(v))
    o["k_h_q10_50_90"] = np.exp(np.percentile(lk, [10, 50, 90])).tolist()
    return o

def main():
    tun = Tuning(ROOT); tau = tun.tau_misfit(); nu0 = tun.eta("tdiaa")["tdi_ml"] / 2
    bf, panel = load_bf()
    out = {"all_edges": stats(edge_view(bf, panel, tau_misfit=tau)),
           "internal_only": stats(edge_view(bf, panel, tau_misfit=tau, internal_only=True)),
           "all_edges_nu0.040": stats(edge_view(with_nu(bf, 0.040, nu0), panel, tau_misfit=tau))}
    print(json.dumps(out, indent=1, default=float))

if __name__ == "__main__":
    main()
