#!/usr/bin/env python3
"""Statistical review (read-only): E1 concavity diagnostics + 'is E4 the time-domain image of E1?' simulation.
  CONFIRM_RUN=1 python scripts/review/e1_concavity_diagnostics.py > logs/review/e1_diag.json
Part A (E1, naive fits, eta fixed or profiled; no SIMEX): replicate c_naive; Q profile over eta for c free
  and c = 1 (offset + linear); x-window fits; binned residual ratio vs a per-HOG linear fit.
Part B: synthetic 3Di = eta + a_h * aa_ml^c (noise-free, c in {1, c_hat}) on the real confirm pair design,
  branch engine refit with the production tip correction, then E3 k (terminal / internal) and E4 point beta.
  Optionally, an additive per-protein excess delta on every non-identical pair (y += 2*delta).
Writes nothing under results/ or data/.
"""
import os, sys, json, tempfile
from pathlib import Path
import numpy as np, polars as pl, torch
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from strucevo.stats.driver import EndpointDriver
from strucevo.stats.endpoints_curve import prepare_curve_data
from strucevo.stats.pairwise import CurveEngine, two_pass
from strucevo.stats.branch import fit_branches
from strucevo.stats.endpoints_branch import edge_view, build_cells, cell_design, fit_mixed, BIN_EDGES

def main():
    out = {}
    d = EndpointDriver(split="confirm", run_alt_alignments=False, outdir=tempfile.mkdtemp(), log=lambda *a: None)
    d.load()
    tab = d.tabs["tdiaa"]
    eta = d.tun.eta("tdiaa")["tdi_ml"]
    hc = pl.read_parquet(ROOT / "data/interim/hogs_endpoint_cols_confirm.parquet")
    excl = set(hc.filter(pl.col("excl_tdi_tdiaa"))["hog"].to_list())

    def fit(x_range=(0.05, 2.0), eta_v=eta, c_fixed=None):
        data, info, idx, sp = prepare_curve_data(tab, "E1", excluded_hogs=excl, x_range=x_range)
        data.eta = eta_v
        eng = CurveEngine(data)
        f = two_pass(eng)
        res = {"c": f["c"], "Q": f["Qmin"], "s2": f["s2"], "n": info["n_pairs_fit"], "H": info["n_hogs"]}
        if c_fixed is not None:
            res["Q_at_c"] = float(eng.Q(torch.tensor([c_fixed]), w=f["w"])[0, 0])
        return res, eng, f, data, idx

    base, eng, f, data, idx = fit(c_fixed=1.0)
    out["naive_eta_fixed"] = base
    w = f["w"]
    # Q profile with weights held at base pass-2 weights
    prof = []
    for ev in [0.0, 0.03, 0.0524, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12]:
        eng.eta = ev
        c, b, q = eng.fit_single(w=w)
        q1 = float(eng.Q(torch.tensor([1.0]), w=w)[0, 0])
        prof.append({"eta": ev, "c_hat": c, "Q_min": q, "Q_c1": q1})
    eng.eta = eta
    out["Q_profile_fixed_w"] = prof
    # windows
    for lo, hi in [(0.05, 0.5), (0.05, 1.0), (0.3, 2.0), (0.5, 2.0), (1.0, 2.0)]:
        try:
            r, *_ = fit(x_range=(lo, hi))
            out[f"window_{lo}_{hi}"] = r
        except Exception as e:
            out[f"window_{lo}_{hi}"] = {"error": repr(e)}
    # binned ratio of observed (y - eta) to per-HOG linear fit a_h^lin x
    x = data.x; y = data.y; g = data.g; el = data.w0 > 0
    wn = w.cpu().numpy()
    a1 = eng.a_h(1.0, w=w).cpu().numpy()
    pred = a1[g] * x
    bins = [0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0]
    rb = []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = el & (x >= lo) & (x < hi)
        rb.append({"x_lo": lo, "x_hi": hi, "n": int(m.sum()),
                   "ratio_obs_over_linear": float(np.sum(y[m] - eta) / np.sum(pred[m])),
                   "mean_y_minus_eta_over_x": float(np.sum(y[m] - eta) / np.sum(x[m]))})
    out["binned_vs_linear"] = rb
    # same with offset model: per-HOG linear with free global eta from profile
    print(json.dumps(out, indent=1, default=float), flush=True)

    # ---------------- Part B: synthetic worlds on the real design
    chat = 0.7226
    tabB = tab
    hog_all = tabB["hog"].to_numpy().astype(str)
    xa = tabB["aa_ml"].to_numpy().astype(float)
    pos = {h: i for i, h in enumerate(data.hog_ids)}
    keep_h = np.array([h in pos for h in hog_all])
    gi = np.array([pos.get(h, 0) for h in hog_all])
    resB = {}
    tv = d.tun.tip_sd("tdiaa") ** 2
    d.tau = d.tun.tau_misfit()
    for name, c_w, delta in [("c1", 1.0, 0.0), ("c0.72", chat, 0.0), ("c1_delta0.014", 1.0, 0.014)]:
        eng.eta = eta
        ah = eng.a_h(c_w, w=w).cpu().numpy()
        ysim = eta + ah[gi] * np.power(np.clip(xa, 1e-6, None), c_w) + 2 * delta
        ysim = np.where(np.isfinite(xa) & keep_h, ysim, np.nan)
        t2 = tabB.with_columns(pl.Series("tdi_ml", ysim))
        t2 = t2.filter(pl.Series(keep_h))
        bf = d.fit(t2, "tdiaa")
        ev = edge_view(bf, d.panel, tau_misfit=d.tau)
        T = ev.is_terminal
        cells = build_cells(ev, d.panel.age, rule="info", info_T=4.0, tip_var=tv)
        yy, xx, zz, _, _ = cell_design(cells)
        fm = fit_mixed(yy, xx, zz, cells["hog"])
        evi = edge_view(bf, d.panel, tau_misfit=d.tau, internal_only=True)
        ci = build_cells(evi, d.panel.age, rule="info", info_T=4.0, tip_var=tv)
        y2, x2, z2, _, _ = cell_design(ci)
        fi = fit_mixed(y2, x2, z2, ci["hog"])
        tau_e = 0.5 * (d.panel.age[ev.child] + d.panel.age[ev.parent])
        gb = np.clip(np.searchsorted(BIN_EDGES, tau_e, side="right") - 1, 0, 4)
        kb = {f"bin{j}_{lab}": float(ev.b_aa[(gb == j) & mm].sum() / ev.b_t[(gb == j) & mm].sum())
              for j in range(5) for lab, mm in (("term", T), ("int", ~T)) if ((gb == j) & mm).sum() > 0}
        resB[name] = {"beta_all": fm["beta"], "beta_internal": fi["beta"], "n_cells": len(yy),
                      "k_all": float(ev.b_aa.sum() / ev.b_t.sum()), "k_term": float(ev.b_aa[T].sum() / ev.b_t[T].sum()),
                      "k_int": float(ev.b_aa[~T].sum() / ev.b_t[~T].sum()), "k_by_bin": kb,
                      "cens": float(cells["censored"].mean())}
        print(json.dumps({name: resB[name]}, indent=1, default=float), flush=True)
    out["synthetic"] = resB

if __name__ == "__main__":
    main()
