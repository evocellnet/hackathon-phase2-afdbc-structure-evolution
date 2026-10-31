#!/usr/bin/env python3
"""Statistical review (read-only): E3/E4 diagnostics on the saved confirm branch table (tdiaa).
Varies the terminal tip correction nu, splits by edge type / stratum / age bin, and refits the E4 point
estimate (MixedLM point only; no bootstrap). Writes nothing under results/endpoints.
  CONFIRM_RUN=1 python scripts/review/e4_tip_diagnostics.py > logs/review_e4_tip.log
"""
import os, sys, json
from pathlib import Path
import numpy as np, polars as pl
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from strucevo.stats.common import guard_split
from strucevo.stats.driver import load_panel, Tuning
from strucevo.stats.endpoints_branch import edge_view, build_cells, cell_design, fit_mixed, BIN_EDGES

class BF:  # minimal stand-in for BranchFit
    pass

def load_bf(aln="tdiaa"):
    E = pl.read_parquet(ROOT / f"data/interim/branches/{aln}.parquet")
    guard_split(E["split"].to_numpy())
    F = pl.read_parquet(ROOT / f"data/interim/branches/{aln}_fit.parquet")
    flags = []
    panel = load_panel(ROOT, flags)
    e = {c: E[c].to_numpy() for c in E.columns}
    e["panel_edge"] = np.array([None if v is None else v for v in E["panel_edge"].to_list()], dtype=object)
    e["_child_idx"] = np.array([panel.index[str(v)] if str(v) in panel.index else panel.index[int(v)] for v in e["panel_child"]])
    e["_parent_idx"] = np.array([panel.index[str(v)] if str(v) in panel.index else panel.index[int(v)] for v in e["panel_parent"]])
    bf = BF(); bf.edges = e; bf.fits = {c: F[c].to_numpy() for c in F.columns}
    return bf, panel

def with_nu(bf, nu_new, nu_old, metric="tdi_ml"):
    e = dict(bf.edges)
    m = (e["metric"] == metric) & e["is_terminal"]
    b = e["b"].astype(float).copy()
    b[m] = b[m] + nu_old - nu_new
    e["b"] = b
    out = BF(); out.edges = e; out.fits = bf.fits
    return out

def e4_point(bf, panel, tau, tip_var, internal_only=False, terminal_scale=None):
    ev = edge_view(bf, panel, tau_misfit=tau, internal_only=internal_only)
    c = build_cells(ev, panel.age, rule="info", info_T=4.0, tip_var=tip_var)
    y, x, z, _, _ = cell_design(c)
    f = fit_mixed(y, x, z, c["hog"])
    return f["beta"], len(y), float(c["censored"].mean()), ev, c

def ksum(ev, mask=None):
    m = np.ones(len(ev.b_aa), bool) if mask is None else mask
    return ev.b_aa[m].sum() / ev.b_t[m].sum()

def main():
    tun = Tuning(ROOT)
    tau = tun.tau_misfit()
    tv = tun.tip_sd("tdiaa") ** 2
    eta = tun.eta("tdiaa")["tdi_ml"]; nu0 = eta / 2
    bf, panel = load_bf()
    out = {"nu0": nu0, "tip_var": tv}
    ev = edge_view(bf, panel, tau_misfit=tau)
    tau_e = 0.5 * (panel.age[ev.child] + panel.age[ev.parent])
    g = np.clip(np.searchsorted(BIN_EDGES, tau_e, side="right") - 1, 0, 4)
    T = ev.is_terminal
    # --- E3 decomposition
    out["k_all"] = ksum(ev); out["k_terminal"] = ksum(ev, T); out["k_internal"] = ksum(ev, ~T)
    out["share_AA_terminal"] = ev.b_aa[T].sum() / ev.b_aa.sum()
    out["share_3Di_terminal"] = ev.b_t[T].sum() / ev.b_t.sum()
    # raw (uncorrected) terminal k
    btraw = ev.b_t + np.where(T, nu0, 0)
    out["k_terminal_nu0"] = ev.b_aa[T].sum() / btraw[T].sum()
    # nu that equalises terminal k with internal k
    need = ev.b_t[T].sum() + nu0 * T.sum() - ev.b_aa[T].sum() / out["k_internal"]
    out["nu_equalising"] = need / T.sum()
    out["n_terminal_edges"] = int(T.sum())
    out["mean_terminal_b_aa"] = float(ev.b_aa[T].mean()); out["median_terminal_b_aa"] = float(np.median(ev.b_aa[T]))
    out["mean_internal_b_aa"] = float(ev.b_aa[~T].mean())
    # by age bin x edge type
    tab = []
    for j in range(5):
        for lab, mm in (("term", T), ("int", ~T)):
            m = (g == j) & mm
            if m.sum() == 0: continue
            tab.append({"bin": j, "type": lab, "n": int(m.sum()), "S_aa": float(ev.b_aa[m].sum()),
                        "S_t": float(ev.b_t[m].sum()), "k": float(ev.b_aa[m].sum() / ev.b_t[m].sum()),
                        "frac_3Di_from_nu": float(nu0 * m.sum() / (ev.b_t[m].sum() + nu0 * m.sum())) if lab == "term" else 0.0})
    out["by_bin_type"] = tab
    # --- E4 point fits
    res = {}
    b, n, cf, _, c = e4_point(bf, panel, tau, tv); res["primary_nu0"] = {"beta": b, "n_cells": n, "cens": cf}
    for mult in (0.0, 0.5, 1.5, 2.0, 3.0):
        b, n, cf, _, _ = e4_point(with_nu(bf, nu0 * mult, nu0), panel, tau, tv)
        res[f"nu_x{mult}"] = {"nu": nu0 * mult, "beta": b, "n_cells": n, "cens": cf}
    b, n, cf, _, _ = e4_point(with_nu(bf, out["nu_equalising"], nu0), panel, tau, tv)
    res["nu_equalising"] = {"nu": out["nu_equalising"], "beta": b, "n_cells": n, "cens": cf}
    b, n, cf, _, _ = e4_point(bf, panel, tau, tv, internal_only=True); res["internal_only"] = {"beta": b, "n_cells": n, "cens": cf}
    # terminal-share of each cell vs residual: fraction of cell 3Di that is terminal
    out["e4"] = res
    # E4 by stratum: hog stratum from hogs.parquet
    hogs = pl.read_parquet(ROOT / "data/interim/hogs.parquet").filter(pl.col("split") == "confirm")
    scol = [c_ for c_ in hogs.columns if "strat" in c_.lower()]
    out["stratum_col"] = scol
    if scol:
        sm = dict(zip(hogs["hog"].to_list(), hogs[scol[0]].to_list()))
        st = {}
        for s in sorted(set(sm.values())):
            sub = [h for h, v in sm.items() if v == s]
            try:
                evs = edge_view(bf, panel, tau_misfit=tau, hog_subset=sub)
                cs = build_cells(evs, panel.age, rule="info", info_T=4.0, tip_var=tv)
                y, x, z, _, _ = cell_design(cs)
                f = fit_mixed(y, x, z, cs["hog"])
                Ts = evs.is_terminal
                st[str(s)] = {"n_hogs": len(evs.hog_ids), "n_cells": len(y), "beta": f["beta"],
                              "k": ksum(evs), "k_term": ksum(evs, Ts), "k_int": ksum(evs, ~Ts),
                              "terminal_share_AA": float(evs.b_aa[Ts].sum() / evs.b_aa.sum())}
            except Exception as ex:
                st[str(s)] = {"error": repr(ex)}
        out["by_stratum"] = st
    print(json.dumps(out, indent=1, default=float))

if __name__ == "__main__":
    main()
