#!/usr/bin/env python3
"""X2 (ANALYSIS_PLAN §6): path and root-to-tip distances vs divergence age per metric; time-curvature fits.

(a) Path distances. Pairs: split == 'tune', primary 1:1 orthologs, status ok, not saturated, HOG not excluded
    for the metric family (no aa_ml_seq window: age, not distance, is the axis). Divergence time
    T = 2 × age(MRCA(species_a, species_b)) from data/interim/timetree/{panel_tree.nwk, node_ages.tsv} (median).
    Age bins (§0) are assigned on the MRCA age t = T/2 (§0's branch-midpoint rule applies to branches; for a
    pair path the node age is the natural analogue). Pairs whose MRCA is undated are counted and dropped.
(b) Root-to-tip. strucevo.stats.branch.fit_branches on all tune primary pairs of --aln (default tdiaa), additive
    metrics only (aa_ml, tdi_ml, lddt_D when present); edge values tip-corrected with ν = η/2 (tuning.yaml, 0 if
    null); HOG×metric fits flagged excluded (underdetermined / fail_frac / misfit) are dropped. Per tip,
    RTT = Σ edge values on the HOG-root→tip path (merged bifurcating roots apportioned by time span, D1); unit =
    per-HOG mean RTT, x = age of the HOG root (MRCA of S_h).
Fits (strucevo.models.fit_model): d = a t^c (power), d = a(1 − e^{−bt}) (saturating), d = a + bt (lin),
    d = bt (lin0). Pair-level least squares is computed exactly from per-unique-age means with weights n (the
    within-age SS is common to all models). Windows: all ages, and ages < 1000 Mya (bins 1-4); bin 5 (≥ 1000
    Mya) is summarised separately. Floor-bearing metrics (tdi_ml, lddt_D, 1 − lddt_global) are fitted as
    d − η. Unadjusted 95% CIs: HOG-cluster bootstrap (B replicates, multinomial HOG weights).
Licence: JSON holds binned summaries and fit parameters only; a bin's mean age is reported only when the bin
    spans >= 3 distinct dated nodes (else null).

Run:  python scripts/exploratory/x2_age.py [--root R] [--aln tdiaa|seqonly|tmalign] [--B 200] [--min-n 20]
Out:  ROOT/results/exploratory/X2_<aln>.json ; ROOT/figures/exploratory/X2_<aln>_{path,rtt}.png
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl

from xcommon import (METRICS, OI, AGE_EDGES, AGE_LABELS, common_args, load_tuning, load_tune_pairs, primary_1to1,
                     load_hog_flags, excluded_hogs, metric_values, ci95, write_json, out_dirs, mpl_style, log, rng,
                     multinomial, load_panel, mrca_age, age_bin)
from strucevo.models import fit_model
from strucevo.stats.branch import fit_branches, edge_values_from_params

PATH_METRICS = ["aa_ml", "aa_ml_seq"] + list(METRICS)
EXTRA = {"aa_ml": dict(col="aa_ml", label="AA ML distance (core, LG+Γ4)", unit="subst./site", sat="aa_sat",
                       excl="aa", eta_key=None),
         "aa_ml_seq": dict(col="aa_ml_seq", label="AA ML distance (seqonly, aa_ml_seq)", unit="subst./site",
                           sat=None, excl=None, eta_key=None)}
MODELS = ("power", "sat", "lin", "lin0")
MIN_UNIQUE_AGES = 3


def mlabel(m):
    s = EXTRA.get(m) or METRICS[m]
    return s["label"] + (f" [{s['unit']}]" if s["unit"] else "")


def path_values(df, m, tuning):
    if m in EXTRA:
        s = EXTRA[m]
        y = df[s["col"]].to_numpy().astype(float)
        sat = df[s["sat"]].to_numpy().astype(bool) if s["sat"] else np.zeros(len(y), bool)
        return y, sat, s["excl"], 0.0
    y, _, sat = metric_values(df, m, tuning)      # metric's own saturation flag (none for p/TM/RMSD/LDDT)
    s = METRICS[m]
    eta = tuning["eta"][s["eta_key"]] if s["eta_key"] else 0.0
    return y, sat, s["excl"], eta


def age_bins_summary(t, y, hog, min_n, x_is_T):
    """Per §0 bin (on node age t): n, n_hogs, n_nodes, distance quantiles, rate quantiles d/T."""
    b = age_bin(t)
    rows = []
    for k, lab in enumerate(AGE_LABELS):
        m = b == k
        n = int(m.sum())
        r = dict(bin=lab, n=n, n_hogs=int(len(np.unique(hog[m]))) if n else 0,
                 n_nodes=int(len(np.unique(t[m]))) if n else 0)
        if n >= min_n:
            q = np.quantile(y[m], [0.25, 0.5, 0.75])
            xx = (2 * t[m]) if x_is_T else t[m]
            rate = y[m] / np.where(xx > 0, xx, np.nan)
            r.update(y_q25=q[0], y_median=q[1], y_q75=q[2], rate_per_My_median=float(np.nanmedian(rate)))
            r["x_mean"] = float(np.mean(xx)) if r["n_nodes"] >= 3 else None
            r["_x_plot"] = float(np.mean(xx))
        rows.append(r)
    return rows


def grouped(x, y, unit_idx, wunit=None):
    """Means of y per unique x (weights = summed unit weights). Returns ux, ybar, n."""
    ux, gi = np.unique(x, return_inverse=True)
    w = np.ones(len(x)) if wunit is None else wunit[unit_idx]
    n = np.bincount(gi, weights=w, minlength=len(ux))
    s = np.bincount(gi, weights=w * y, minlength=len(ux))
    ok = n > 0
    return ux[ok], s[ok] / n[ok], n[ok]


def fit_window(x, y, unit, n_units, B, r, label):
    """Fit MODELS to (x, y) aggregated by unique x; HOG/unit-cluster bootstrap of the key parameters."""
    out = {"window": label, "n": int(len(x)), "n_unique_ages": int(len(np.unique(x)))}
    if out["n_unique_ages"] < MIN_UNIQUE_AGES or len(x) < 10:
        out["skipped"] = f"< {MIN_UNIQUE_AGES} distinct ages (or < 10 points): curvature not identifiable"
        return out, None
    ux, yb, n = grouped(x, y, unit)
    fits = {m: fit_model(m, ux, yb, weights=n) for m in MODELS}
    best = min(fits.values(), key=lambda f: f.aicc).aicc
    out["fits"] = {m: {**f.params, "aicc_grouped": f.aicc, "dAICc": f.aicc - best, "success": f.success}
                   for m, f in fits.items()}
    W = multinomial(r, n_units, B)
    reps = {"power_c": [], "sat_b": [], "lin_b": [], "lin0_b": []}
    for b in range(B):
        ux_b, yb_b, n_b = grouped(x, y, unit, W[b])
        if len(ux_b) < MIN_UNIQUE_AGES:
            continue
        reps["power_c"].append(fit_model("power", ux_b, yb_b, n_b).params["c"])
        reps["sat_b"].append(fit_model("sat", ux_b, yb_b, n_b).params["b"])
        reps["lin_b"].append(fit_model("lin", ux_b, yb_b, n_b).params["b"])
        reps["lin0_b"].append(fit_model("lin0", ux_b, yb_b, n_b).params["b"])
    out["ci95"] = {"power_c": ci95(fits["power"].params["c"], reps["power_c"]),
                   "sat_b_per_My": ci95(fits["sat"].params["b"], reps["sat_b"]),
                   "lin_slope_per_My": ci95(fits["lin"].params["b"], reps["lin_b"]),
                   "lin0_slope_per_My": ci95(fits["lin0"].params["b"], reps["lin0_b"])}
    return out, fits


def analyse(x_age, y, hog, t_node, B, min_n, seed, x_is_T):
    """x_age: fit axis (T for path, t_root for RTT); t_node: node age used for §0 bins."""
    uh, unit = np.unique(hog, return_inverse=True)
    rec = {"bins": age_bins_summary(t_node, y, hog, min_n, x_is_T), "windows": {}}
    fitobjs = {}
    for wi, (wl, m) in enumerate([("all_ages", np.ones(len(y), bool)), ("lt1000Mya_bins1-4", t_node < 1000.0)]):
        sel = m & (x_age > 0)
        uh2, unit2 = np.unique(unit[sel], return_inverse=True)
        w, fo = fit_window(x_age[sel], y[sel], unit2, len(uh2), B, rng("X2", *seed, wi), wl)
        rec["windows"][wl] = w
        fitobjs[wl] = fo
    b5 = t_node >= 1000.0
    rec["bin5_separate"] = {"n": int(b5.sum()), "n_hogs": int(len(np.unique(hog[b5]))),
                            "y_median": float(np.median(y[b5])) if b5.any() else None,
                            "note": "bin 5 (>= 1000 Mya) shown separately in figures (distinct marker); "
                                    "'lt1000Mya' window excludes it"}
    return rec, fitobjs


def rtt_table(bf, panel, metrics, nu):
    """Rows per HOG x metric x node: level 'root' (HOG root -> all tips) and 'internal' (every internal HOG node,
    root included -> its descendant tips). node_to_tip(e) = mean over tips k below e of the summed edge values
    of the edges strictly below e on the path to k."""
    rows = {"hog": [], "metric": [], "level": [], "t_node": [], "rtt_mean": [], "rtt_sd": [], "n_tips": []}
    n_excl = {m: 0 for m in metrics}
    for hf in bf.hogs:
        ht = hf.ht
        M = ht.member.astype(float)                                  # (n_sp, E)
        Mb = ht.member
        # sub[e, f]: edge f lies strictly below the child node of edge e
        sub = (Mb[:, None, :] <= Mb[:, :, None]).all(0) & ~np.eye(Mb.shape[1], dtype=bool)
        internal = np.flatnonzero(~ht.is_terminal)
        for k, m in enumerate(metrics):
            if hf.excluded.get(m):
                n_excl[m] += 1
                continue
            ev = edge_values_from_params(hf, hf.fit["b"][k], nu.get(m, 0.0))
            nodes = [("root", float(panel.age[ht.root]), np.ones(M.shape[0], bool), M @ ev)]
            for e in internal:
                tips = Mb[:, e]
                nodes.append(("internal", float(panel.age[ht.edge_child[e]]), tips, M @ (ev * sub[e])))
            for lvl, age, tips, d in nodes:
                for lv in ([lvl, "internal"] if lvl == "root" else [lvl]):
                    rows["hog"].append(hf.hog); rows["metric"].append(m); rows["level"].append(lv)
                    rows["t_node"].append(age); rows["rtt_mean"].append(float(np.mean(d[tips])))
                    rows["rtt_sd"].append(float(np.std(d[tips]))); rows["n_tips"].append(int(tips.sum()))
    return pl.DataFrame(rows), n_excl


def plot(recs, fitobjs, kind, xlab, path, title):
    plt = mpl_style()
    ms = [m for m in recs if "bins" in recs[m]]
    nc = min(5, max(1, len(ms)))
    nr = int(np.ceil(len(ms) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(3.6 * nc, 3.2 * nr), constrained_layout=True, squeeze=False)
    for ax, m in zip(axs.flat, ms):
        r = recs[m]
        for bi, b in enumerate(r["bins"]):
            if "_x_plot" not in b:
                continue
            is5 = bi == 4
            ax.errorbar([b["_x_plot"]], [b["y_median"]], yerr=[[b["y_median"] - b["y_q25"]], [b["y_q75"] - b["y_median"]]],
                        fmt="s" if is5 else "o", color=OI[6] if is5 else OI[0], ms=6, capsize=3,
                        mfc="white" if is5 else OI[0], label=("bin 5 (≥1000 Mya)" if is5 else
                                                              ("bin median ± IQR (§0 bins)" if bi == min(
                                                                  i for i, q in enumerate(r["bins"]) if "_x_plot" in q) else None)))
        xmax = max([b["_x_plot"] for b in r["bins"] if "_x_plot" in b] + [1.0]) * 1.05
        xs = np.linspace(1e-3, xmax, 300)
        fo = fitobjs.get(m, {})
        for wl, ls in [("all_ages", "-"), ("lt1000Mya_bins1-4", "--")]:
            f = fo.get(wl)
            if not f:
                continue
            for mi, mn in enumerate(("power", "sat", "lin")):
                if wl != "all_ages" and mn != "power":
                    continue
                p = f[mn].params
                lab = {"power": f"a·t^c (c={p.get('c', np.nan):.2f})", "sat": "a(1−e^{−bt})",
                       "lin": "a + b·t"}[mn] + ("" if wl == "all_ages" else " [<1000 Mya]")
                ax.plot(xs, f[mn].predict(xs) + r.get("eta_subtracted", 0.0), color=OI[1 + mi], ls=ls, lw=1.6, label=lab)
        ax.set_title(f"{m}", fontsize=9)
        ax.set_xlabel(xlab, fontsize=8)
        ax.set_ylabel(mlabel(m), fontsize=7.5)
        ax.set_xlim(0, xmax)
        ax.legend(fontsize=6, loc="upper left")
    for ax in list(axs.flat)[len(ms):]:
        ax.axis("off")
    import textwrap
    fig.suptitle(textwrap.fill(title, int(18 * 3.6 * nc)), fontsize=9.5)
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = common_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter))
    a = ap.parse_args()
    t0 = time.time()
    root = Path(a.root)
    res, figd = out_dirs(root)
    tuning = load_tuning(root, a.aln)
    panel = load_panel(root)
    allp = load_tune_pairs(root, a.aln)
    hflags = load_hog_flags(root)
    df = primary_1to1(allp)
    sa, sb = df["species_a"].to_numpy(), df["species_b"].to_numpy()
    t = mrca_age(panel, sa, sb)
    hog = df["hog"].to_numpy().astype(str)
    status_ok = df["status"].to_numpy() == "ok"
    out = {"analysis": "X2", "aln": a.aln, "split": "tune", "root": str(root), "B": a.B, "min_n_bin": a.min_n,
           "tuning": tuning, "age_bins_Mya": AGE_EDGES[:-1] + ["inf"], "age_bin_labels": AGE_LABELS,
           "ages": "panel node_ages.tsv median; T = 2 x age(MRCA)", "ci": "unadjusted 95% (HOG-cluster bootstrap)",
           "seed": "SeedSequence(20261030, spawn_key=(1002, part, metric_idx, window))",
           "n_primary_1to1": df.height, "n_pairs_mrca_undated": int(np.isnan(t).sum())}
    log(f"X2 {a.aln}: {df.height} primary 1:1 pairs; undated MRCA {out['n_pairs_mrca_undated']}")
    # ---------------- (a) path distances
    precs, pfits = {}, {}
    for j, m in enumerate(PATH_METRICS):
        y, sat, fam, eta = path_values(df, m, tuning)
        excl = excluded_hogs(hflags, fam, a.aln) if fam else set()
        ok = status_ok & ~sat & np.isfinite(y) & np.isfinite(t) & ~np.isin(hog, list(excl))
        if ok.sum() < 10:
            precs[m] = {"skipped": "no eligible values (e.g. lddt_D before floors)"}
            log(f"  path {m}: skipped")
            continue
        r, fo = analyse(2 * t[ok], y[ok] - eta, hog[ok], t[ok], a.B, a.min_n, (0, j), True)
        r.update(n_eligible=int(ok.sum()), n_hogs=int(len(np.unique(hog[ok]))), eta_subtracted=eta,
                 n_hogs_excluded_flag=len(excl))
        precs[m], pfits[m] = r, fo
        w = r["windows"]["all_ages"]
        log(f"  path {m:12s} n={ok.sum()} c_all={w.get('ci95', {}).get('power_c')} "
            f"c_lt1000={r['windows']['lt1000Mya_bins1-4'].get('ci95', {}).get('power_c', {}).get('estimate')} "
            f"{w.get('skipped', '')}")
    # ---------------- (b) root-to-tip from branch fits
    prim = allp.filter(pl.col("pair_set") == "primary")
    bmetrics = ["aa_ml", "tdi_ml"] + (["lddt_D"] if np.isfinite(prim["lddt_D"].to_numpy().astype(float)).any() else [])
    nu = {"tdi_ml": tuning["eta"]["tdi_ml"] / 2, "lddt_D": tuning["eta"]["lddt_D"] / 2}
    sel = {h: list(s) for h, s in zip(hflags["hog"].to_list(), hflags["selected"].to_list())} \
        if "selected" in hflags.columns else {}
    hs = {}
    for h, g in prim.group_by("hog"):
        sp = set(np.r_[g["species_a"].to_numpy(), g["species_b"].to_numpy()].tolist())
        S = sel.get(h[0])
        hs[h[0]] = np.asarray(sorted(S) if S and sp <= set(S) else sorted(sp), np.int64)
    t1 = time.time()
    bf = fit_branches(prim, panel, metrics=bmetrics, nu=nu, aln=a.aln, hog_species=hs,
                      tau_misfit=tuning["tau_misfit"], misfit=tuning["tau_misfit"] is not None, display_nnls=False)
    out["rtt_branch_fit"] = {"metrics": bmetrics, "nu": nu, "n_hogs": len(bf.hogs), "n_pairs": prim.height,
                             "runtime_s": time.time() - t1, "tau_misfit": tuning["tau_misfit"],
                             "kappa": "1 (Fisher SE; tuning kappa null)" if tuning["kappa"]["tdi_ml"] is None else "tuned",
                             "merged_root": "bifurcating HOG roots apportioned by panel time span (D1)"}
    rt, n_excl = rtt_table(bf, panel, bmetrics, nu)
    out["rtt_branch_fit"]["n_hog_metric_excluded"] = n_excl
    lvl_recs = {}
    for li, lvl in enumerate(["root", "internal"]):
        rrecs, rfits = {}, {}
        for j, m in enumerate(bmetrics):
            s_ = rt.filter((pl.col("metric") == m) & (pl.col("level") == lvl) & pl.col("t_node").is_finite())
            if s_.height < 3:
                rrecs[m] = {"skipped": "too few HOG nodes"}
                continue
            tr = s_["t_node"].to_numpy()
            r, fo = analyse(tr, s_["rtt_mean"].to_numpy(), s_["hog"].to_numpy().astype(str), tr, a.B,
                            max(3, a.min_n // 10), (1 + li, j), False)
            r.update(n_units=s_.height, n_hogs=int(s_["hog"].n_unique()),
                     rtt_sd_within_node_median=float(np.median(s_["rtt_sd"].to_numpy())),
                     n_nodes_undated=int(rt.filter((pl.col("metric") == m) & (pl.col("level") == lvl)
                                                   & ~pl.col("t_node").is_finite()).height))
            rrecs[m], rfits[m] = r, fo
            log(f"  {lvl:8s} {m:8s} units={s_.height} c_all="
                f"{r['windows']['all_ages'].get('ci95', {}).get('power_c', {}).get('estimate')} "
                f"c_lt1000={r['windows']['lt1000Mya_bins1-4'].get('ci95', {}).get('power_c', {}).get('estimate')} "
                f"{r['windows']['all_ages'].get('skipped', '')}")
        lvl_recs[lvl] = (rrecs, rfits)
    plot(precs, pfits, "path", "divergence time T = 2 × age(MRCA) [My]", figd / f"X2_{a.aln}_path.png",
         f"X2: pair (path) distance vs divergence time  (tune, primary 1:1 orthologs, {a.aln}; "
         "points = §0 age-bin median ± IQR at the bin's mean T; lines = fits on all ages (solid) / <1000 Mya (dashed))")
    plot(*lvl_recs["root"], "rtt", "age of HOG root [My]", figd / f"X2_{a.aln}_rtt.png",
         f"X2: per-HOG mean root-to-tip distance (branch engine, {a.aln}) vs HOG root age "
         "(points = §0 bin median ± IQR over HOGs; fits need >= 3 distinct root ages)")
    plot(*lvl_recs["internal"], "node", "age of HOG node [My]", figd / f"X2_{a.aln}_node_to_tip.png",
         f"X2: node-to-tip distance (every internal HOG node, branch engine, {a.aln}) vs node age "
         "(points = §0 bin median ± IQR over HOG nodes; bootstrap clusters = HOGs)")
    strip = lambda d: {m: {k: ([{kk: vv for kk, vv in b.items() if not kk.startswith("_")} for b in v]
                               if k == "bins" else v) for k, v in r.items()} for m, r in d.items()}
    out["path"] = strip(precs)
    out["rtt"] = strip(lvl_recs["root"][0])
    out["node_to_tip"] = strip(lvl_recs["internal"][0])
    out["keys"] = {"bins.x_mean": "mean fit-axis value in the bin (path: T = 2 x MRCA age; rtt/node: node age), "
                   "null unless the bin spans >= 3 dated nodes (licence)",
                   "node_to_tip": "unit = HOG x internal node (root included); bins on node age; "
                   "bootstrap resamples HOGs"}
    out["runtime_s"] = time.time() - t0
    write_json(out, res / f"X2_{a.aln}.json")
    log(f"X2 done in {out['runtime_s']:.1f}s")


if __name__ == "__main__":
    main()
