#!/usr/bin/env python3
"""X1 (ANALYSIS_PLAN §6): curves of structural / sequence distances vs AA ML distance, with naive power-law ĉ.

Pairs: split == 'tune', pair_set == 'primary', one_to_one (1:1 orthologs), E1 filters (§4.1) applied by
strucevo.stats.endpoints_curve.prepare_curve_data: aa_ml_seq in [0.05, 2]; HOG not excluded for the metric
family (hogs.parquet excl_{aa|tdi|lddt}_{aln}); failed (status != ok) or saturated pairs get weight 0;
HOG kept if >= 10 eligible pairs and max(x)/min(x) >= 4.

x in {aa_ml (core ML; aa_sat pairs weight 0), aa_ml_seq (seqonly selection distance)}.
y in {1 − tm_short, rmsd_core [Å], aa_p, tdi_p, 1 − lddt_core, 1 − lddt_global (raw 1 − S), tdi_ml, lddt_D}.

Model (same profiling as E1, strucevo.stats.pairwise): y = η + a_h x^c, a_h profiled per HOG, one global c;
two-pass weights w = 1/(SE_y² + s²).
  - SE_y: tdi_ml, lddt_D = Fisher SE (× kappa if tuned); aa_p, tdi_p = binomial sqrt(p(1−p)/n_core), floor
    1/n_core; tm, rmsd, lddt_core, lddt_global have no per-pair SE -> SE_y := 1, i.e. equal weights (OLS).
  - η: tdi_ml, lddt_D, lddt_global from tuning.yaml eta.<aln> (0 if null, recorded); η = 0 for the metrics
    without a floor estimate in the plan (tm, rmsd, aa_p, tdi_p, lddt_core) -- a non-zero AF-vs-AF floor
    in those biases the naive ĉ downwards (documented, not corrected).
  - No SIMEX: naive ĉ only. Unadjusted 95% CI = HOG-cluster bootstrap (multinomial HOG weights, B
    replicates; pass-2 weights held fixed, D5), normal (± 1.96 SE_boot) and percentile.

Run:  python scripts/exploratory/x1_curves.py [--root R] [--aln tdiaa|seqonly|tmalign] [--B 200] [--min-n 20]
Out:  ROOT/results/exploratory/X1_<aln>.json ; ROOT/figures/exploratory/X1_<aln>_<xvar>.png
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import polars as pl
import torch
torch.set_num_threads(2)

from xcommon import (METRICS, XVARS, OI, common_args, load_tuning, load_tune_pairs, primary_1to1, load_hog_flags,
                     excluded_hogs, metric_values, binned_summary, ci95, write_json, out_dirs, mpl_style, log, rng,
                     multinomial)
from strucevo.stats.endpoints_curve import prepare_curve_data
from strucevo.stats.pairwise import CurveEngine, two_pass, fit_batch, reference_fit, weighted_median, fit_eta_free

X_EDGES = np.round(np.arange(0.0, 2.51, 0.1), 2)


def build_tab(df, xv, m, tuning):
    xs = XVARS[xv]
    y, se, sat = metric_values(df, m, tuning)
    n = df.height
    x = df[xs["col"]].to_numpy().astype(float)
    xse = df[xs["se"]].to_numpy().astype(float) if xs["se"] else np.ones(n)
    xsat = df[xs["sat"]].to_numpy().astype(bool) if xs["sat"] else np.zeros(n, bool)
    return pl.DataFrame({
        "split": df["split"], "hog": df["hog"], "one_to_one": df["one_to_one"], "status": df["status"],
        "species_a": df["species_a"], "species_b": df["species_b"], "aa_ml_seq": df["aa_ml_seq"],
        "aa_ml": x, "aa_ml_se": xse, "aa_sat": xsat,
        "y": y, "y_se": se if se is not None else np.ones(n), "y_sat": sat})


def fit_one(df, xv, m, tuning, hflags, aln, B, min_n, seed_key):
    spec = METRICS[m]
    y_all, _, _ = metric_values(df, m, tuning)
    if not np.isfinite(y_all).any():
        return {"skipped": f"{spec['col']} is all-NaN in this table (e.g. lddt_D before floors are computed)"}
    eta = tuning["eta"][spec["eta_key"]] if spec["eta_key"] else 0.0
    tab = build_tab(df, xv, m, tuning)
    excl = excluded_hogs(hflags, spec["excl"], aln)
    data, info, idx, sp = prepare_curve_data(tab, endpoint="E1", excluded_hogs=excl, ycols=("y", "y_se", "y_sat"))
    rec = {"metric": m, "x": xv, "eta": eta, "eta_source": ("tuning.yaml eta.%s.%s" % (aln, spec["eta_key"]))
           if spec["eta_key"] else "fixed 0 (no floor in plan)", "weights": (
               "1/(SE^2+s^2), SE=" + spec["se"]) if spec["se"] else "equal (no per-pair SE)", **info,
           "n_hogs_excluded_flag": len(excl)}
    if data.H < 2 or info["n_pairs_fit"] < 20:
        rec["skipped"] = "too few HOGs/pairs after filters"
        return rec
    data.eta = eta
    eng = CurveEngine(data)
    f = two_pass(eng)
    c = f["c"]
    w2 = f["w"].cpu().numpy()
    c_ref = reference_fit(data.g, data.x, data.y, w2, eta)
    # HOG-cluster bootstrap
    W = multinomial(rng("X1", *seed_key), data.H, B)
    wb = f["w"][None] * torch.as_tensor(W[:, data.g], dtype=torch.float64, device=eng.dev)
    reps, bnd = fit_batch(eng, w=wb, window=(c - 0.3, c + 0.3))
    a_h = f["a_h"]
    npairs_h = np.bincount(data.g, weights=data.w0, minlength=data.H)
    a_disp = weighted_median(a_h, npairs_h)
    el = data.w0 > 0
    # sensitivity (point only): eta free in [-median(y), q05(y)] -- shows how much a floor drives c
    ye = data.y[el]
    try:
        ef = fit_eta_free(eng, f["w"], eta_bounds=(-float(np.median(ye)) - eta, float(np.quantile(ye, 0.05)) - eta))
    except Exception as e:  # pragma: no cover
        ef = {"error": repr(e)}
    rec["c_eta_free_point"] = ef
    rec.update({
        "c_naive": ci95(c, reps), "grid_boundary": bool(f["boundary"]), "boot_boundary_frac": float(np.mean(bnd)),
        "s2": f["s2"], "c_reference_recompute": c_ref, "c_reference_absdiff": abs(c_ref - c),
        "a_h_summary": {"median_npair_weighted": a_disp, "q25": float(np.nanquantile(a_h, 0.25)),
                        "q75": float(np.nanquantile(a_h, 0.75))},
        "binned": binned_summary(data.x[el], data.y[el], X_EDGES, min_n),
    })
    rec["_plot"] = (data.x[el], data.y[el], eta, a_disp, c)
    return rec


def plot(recs, xv, aln, path, n_tot):
    plt = mpl_style()
    ms = [m for m in METRICS if m in recs]
    nc = 4
    nr = int(np.ceil(len(ms) / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(4.0 * nc, 3.3 * nr), constrained_layout=True, squeeze=False)
    for ax, m in zip(axs.flat, ms):
        r = recs[m]
        ylab = METRICS[m]["label"] + (f" [{METRICS[m]['unit']}]" if METRICS[m]["unit"] else "")
        ax.set_ylabel(ylab, fontsize=8)
        ax.set_xlabel(XVARS[xv]["label"], fontsize=8)
        if "_plot" not in r:
            import textwrap
            ax.text(0.5, 0.5, "not available:\n" + textwrap.fill(r.get("skipped", ""), 40), ha="center",
                    va="center", transform=ax.transAxes, fontsize=8)
            continue
        x, y, eta, a, c = r["_plot"]
        ax.scatter(x, y, s=3, color=OI[7], alpha=0.25, lw=0, rasterized=True, label="pairs (eligible)")
        b = [q for q in r["binned"] if "y_median" in q]
        if b:
            xm = [q["x_median"] for q in b]
            ax.fill_between(xm, [q["y_q25"] for q in b], [q["y_q75"] for q in b], color=OI[0], alpha=0.2, lw=0,
                            label="IQR (0.1 bins)")
            ax.plot(xm, [q["y_median"] for q in b], color=OI[0], lw=2, marker="o", ms=3, label="bin median")
        xs = np.linspace(max(1e-3, np.min(x)), np.max(x), 200)
        cc = r["c_naive"]
        ci = cc["ci95_pct"] or [np.nan, np.nan]
        ax.plot(xs, eta + a * xs ** c, color=OI[1], lw=2,
                label=f"η + ã·x^ĉ, ĉ={c:.2f} [{ci[0]:.2f}, {ci[1]:.2f}]")
        ax.set_ylim(min(0, float(np.min(y))), float(np.quantile(y, 0.995)) * 1.1)
        ax.set_title(f"{m}  (n={r['n_pairs_fit']}, HOGs={r['n_hogs']})", fontsize=9)
        ax.legend(fontsize=6.5, loc="upper left")
    for ax in list(axs.flat)[len(ms):]:
        ax.axis("off")
    fig.suptitle(f"X1: distance curves vs {xv}  (tune split, primary 1:1 orthologs, {aln}; "
                 f"ã = pair-weighted median of per-HOG a_h; 95% CI = HOG bootstrap percentile)", fontsize=10)
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = common_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter))
    ap.add_argument("--metrics", nargs="*", default=list(METRICS))
    a = ap.parse_args()
    t0 = time.time()
    root = Path(a.root)
    res, figd = out_dirs(root)
    tuning = load_tuning(root, a.aln)
    df = load_tune_pairs(root, a.aln)
    n_in = df.height
    df = primary_1to1(df)
    hflags = load_hog_flags(root)
    log(f"X1 {a.aln}: {n_in} tune rows -> {df.height} primary 1:1")
    out = {"analysis": "X1", "aln": a.aln, "split": "tune", "root": str(root), "B": a.B, "min_n_bin": a.min_n,
           "n_rows_tune": n_in, "n_primary_1to1": df.height, "tuning": {k: v for k, v in tuning.items()},
           "bins_x": X_EDGES.tolist(), "ci": "unadjusted 95% (HOG-cluster bootstrap)",
           "seed": "SeedSequence(20261030, spawn_key=(1001, xvar_idx, metric_idx))",
           "device": "cuda" if torch.cuda.is_available() else "cpu", "results": {}}
    for i, xv in enumerate(XVARS):
        recs = {}
        for j, m in enumerate(a.metrics):
            t1 = time.time()
            r = fit_one(df, xv, m, tuning, hflags, a.aln, a.B, a.min_n, (i, j))
            r["runtime_s"] = time.time() - t1
            recs[m] = r
            cc = r.get("c_naive", {})
            log(f"  {xv:9s} {m:12s} c={cc.get('estimate', float('nan')):.3f} CI={cc.get('ci95_pct')} "
                f"etafree={r.get('c_eta_free_point', {}).get('c')} eta={r.get('c_eta_free_point', {}).get('eta')} "
                f"n={r.get('n_pairs_fit')} H={r.get('n_hogs')} {r.get('skipped', '')} ({r['runtime_s']:.1f}s)")
        plot(recs, xv, a.aln, figd / f"X1_{a.aln}_{xv}.png", df.height)
        out["results"][xv] = {m: {k: v for k, v in r.items() if not k.startswith("_")} for m, r in recs.items()}
    out["runtime_s"] = time.time() - t0
    write_json(out, res / f"X1_{a.aln}.json")
    log(f"X1 done in {out['runtime_s']:.1f}s")


if __name__ == "__main__":
    main()
