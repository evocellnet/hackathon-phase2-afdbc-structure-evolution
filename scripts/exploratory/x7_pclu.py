#!/usr/bin/env python3
"""X7 (ANALYSIS_PLAN §6): P-clu (AFDB-cluster member-vs-representative) curves vs ortholog curves.

The ortholog curve is computed with the SAME code as scripts/pclu/analyse_pclu.py (imported, not copied):
  x = fident binned in 0.02 bins (fbin = floor(fident/0.02), clipped 0..49), per-bin n, median and IQR of TM and
  LDDT (`binned`), bins with n < min_n dropped (P-clu used 100; --min-n here), and the same curve families
  (`FAM`/`fit_table`) fitted to n-weighted bin medians of 1 − TM and 1 − LDDT against p = 1 − fident,
  dP = −ln(fident) and dK (Kimura; only for windows with fident >= 0.16), windows as in `run_fits`.
Ortholog pairs: split == 'tune', primary 1:1 orthologs, status ok (no aa_ml_seq window, as P-clu has none).
Variable mapping (P-clu is Foldseek file 8; documented mismatches):
  fident -> 1 − aa_p (AA identity over core columns of the --aln alignment; primary)    [P-clu: all aligned cols]
         -> aa_pid_seq (seqonly identity over all aligned columns; sensitivity 'x=aa_pid_seq')
  tm     -> tm_short (TM normalised by the shorter protein = query)    [P-clu: alntmscore, normalised by aln length]
  lddt   -> lddt_core (mean per-residue LDDT over aligned core residues)  [P-clu: Foldseek lddt over aligned residues]
         (lddt_global = Σ/len_a reported as a sensitivity)
P-clu strata overlaid: overall, cluster-balanced (<=5/cluster), Euk-Euk pairs, cluster LCA = one superkingdom,
  representative pLDDT >= 90.  Differences per bin = ortholog − P-clu median (same fbin).

Run:  python scripts/exploratory/x7_pclu.py [--root R] [--aln tdiaa] [--min-n 20] [--pclu-dir R/results/pclu]
Out:  ROOT/results/exploratory/X7_<aln>.json ; ROOT/figures/exploratory/X7_<aln>_<xvar>.png
"""
from __future__ import annotations

import argparse
import importlib.util
import time
from pathlib import Path

import numpy as np
import polars as pl

from xcommon import (R_DEFAULT, OI, common_args, load_tune_pairs, primary_1to1, write_json, out_dirs, mpl_style, log)


def load_pclu_module():
    spec = importlib.util.spec_from_file_location("analyse_pclu", R_DEFAULT / "scripts" / "pclu" / "analyse_pclu.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)          # top level only defines constants/functions (and mkdirs existing dirs)
    return mod


def strata(pdir: Path):
    rd = lambda f: pl.read_csv(pdir / f, separator="\t")
    s = {"overall": rd("binned_curve_overall.tsv"),
         "cluster_balanced_max5": rd("binned_curve_cluster_balanced_max5.tsv"),
         "repPlddt_ge90": rd("binned_curve_repPlddt_ge90.tsv")}
    ps = rd("binned_curve_by_pair_superkingdoms.tsv")
    s["pair_Euk-Euk"] = ps.filter(pl.col("level") == "Euk-Euk").drop(["level", "stratum_var"])
    la = rd("binned_curve_by_lca_age_class.tsv")
    s["lca_one_superkingdom"] = la.filter(pl.col("level") == "4 one superkingdom").drop(["level", "stratum_var"])
    return s


def ortholog_frame(df, xv, lddt_col, BIN):
    f = (1.0 - df["aa_p"].to_numpy().astype(float)) if xv == "x=1-aa_p(core)" else df["aa_pid_seq"].to_numpy().astype(float)
    p = pl.DataFrame({"fident": f, "tm": df["tm_short"].to_numpy().astype(float),
                      "lddt": df[lddt_col].to_numpy().astype(float)})
    p = p.filter(pl.col("fident").is_finite())
    if p["tm"].is_nan().all():
        p = p.with_columns(pl.lit(None, dtype=pl.Float64).alias("tm"))
    else:
        p = p.with_columns(pl.col("tm").fill_nan(None))
    return p.with_columns(pl.col("lddt").fill_nan(None),
                          (-pl.col("fident").clip(1e-6, 1).log()).alias("dP"),
                          (pl.col("fident") / BIN).floor().clip(0, 49).cast(pl.Int32).alias("fbin"))


def fits_bin_medians(mod, b, windows):
    """Same as analyse_pclu.run_fits, bin-median part only."""
    rows = []
    for wname, (lo, hi) in windows.items():
        for ycol, ylab in [("tm", "1-TM"), ("lddt", "1-LDDT")]:
            bb = b.filter((pl.col("fident_center") >= lo) & (pl.col("fident_center") < hi)).drop_nulls(f"{ycol}_med")
            if bb.height < 4:
                continue
            fm = bb["fident_med"].to_numpy()
            y = 1 - bb[f"{ycol}_med"].to_numpy()
            n = bb["n"].to_numpy().astype(float)
            for xname, x in [("p=1-fident", 1 - fm), ("dP=-ln(fident)", -np.log(np.clip(fm, 1e-6, 1))),
                             ("dK(Kimura)", -np.log(np.clip(1 - (1 - fm) - 0.2 * (1 - fm) ** 2, 1e-12, None)))]:
                if xname.startswith("dK") and lo < 0.16:
                    continue
                ok = np.isfinite(x) & np.isfinite(y)
                rows += mod.fit_table(x[ok], y[ok], n[ok], dict(data="bin_medians_weighted_n", window=wname,
                                                                y=ylab, x=xname))
    return rows


def diffs(orth, pc, name):
    j = orth.join(pc, on="fbin", how="inner", suffix="_pclu")
    out = []
    for r in j.sort("fbin").iter_rows(named=True):
        out.append(dict(fbin=r["fbin"], fident_center=r["fident_center"], n_orth=r["n"], n_pclu=r["n_pclu"],
                        tm_med_orth=r["tm_med"], tm_med_pclu=r["tm_med_pclu"],
                        d_tm=(None if r["tm_med"] is None else r["tm_med"] - r["tm_med_pclu"]),
                        lddt_med_orth=r["lddt_med"], lddt_med_pclu=r["lddt_med_pclu"],
                        d_lddt=(None if r["lddt_med"] is None else r["lddt_med"] - r["lddt_med_pclu"])))
    summ = {}
    for k in ("d_tm", "d_lddt"):
        v = np.array([o[k] for o in out if o[k] is not None], float)
        w = np.array([o["n_orth"] for o in out if o[k] is not None], float)
        summ[k] = dict(n_bins=int(v.size), n_weighted_mean=float(np.sum(v * w) / np.sum(w)) if v.size else None,
                       min=float(v.min()) if v.size else None, max=float(v.max()) if v.size else None)
    return {"stratum": name, "per_bin": out, "summary": summ}


def plot(orth, pcs, dd, xv, aln, lddt_col, path):
    plt = mpl_style()
    fig, axs = plt.subplots(2, 2, figsize=(11.5, 7.6), constrained_layout=True, sharex=True,
                            gridspec_kw={"height_ratios": [2, 1]})
    cols = {"overall": OI[6], "cluster_balanced_max5": OI[7], "pair_Euk-Euk": OI[2], "lca_one_superkingdom": OI[3],
            "repPlddt_ge90": OI[4]}
    for j, (yc, lab) in enumerate([("tm", "median TM-score"), ("lddt", "median LDDT")]):
        ax = axs[0, j]
        for name, b in pcs.items():
            ax.plot(b["fident_center"], b[f"{yc}_med"], color=cols[name], lw=1.4,
                    ls="-" if name == "overall" else "--", label=f"P-clu {name}")
        o = orth.drop_nulls(f"{yc}_med")
        if o.height:
            ax.fill_between(o["fident_center"], o[f"{yc}_q25"], o[f"{yc}_q75"], color=OI[0], alpha=0.2, lw=0,
                            label="orthologs IQR")
            ax.plot(o["fident_center"], o[f"{yc}_med"], color=OI[0], lw=2.4, marker="o", ms=3.5,
                    label=f"orthologs (tune, {aln})")
        else:
            ax.text(0.5, 0.5, f"no ortholog {yc} for {aln}", transform=ax.transAxes, ha="center")
        ax.set_ylabel(lab + (" (orthologs: tm_short; P-clu: alntmscore)" if yc == "tm"
                             else f" (orthologs: {lddt_col}; P-clu: Foldseek lddt)"), fontsize=8)
        ax.set_ylim(0.2, 1.0)
        ax.legend(fontsize=7, loc="lower right")
        ax2 = axs[1, j]
        for d in dd:
            v = [(r["fident_center"], r[f"d_{yc}"]) for r in d["per_bin"] if r[f"d_{yc}"] is not None]
            if v:
                ax2.plot(*zip(*v), color=cols[d["stratum"]], lw=1.4, marker=".",
                         ls="-" if d["stratum"] == "overall" else "--", label=d["stratum"])
        ax2.axhline(0, color="k", lw=0.8)
        ax2.set_ylabel(f"Δ median {yc.upper()} (orthologs − P-clu)", fontsize=8)
        ax2.set_xlabel(f"sequence identity fident (0.02 bins; orthologs: {xv[2:]})", fontsize=8)
        ax2.legend(fontsize=6.5)
    axs[1, 0].set_xlim(0, 1)
    fig.suptitle(f"X7: AFDB-cluster (P-clu) vs ortholog median curves  [{aln}, {xv}]", fontsize=10)
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = common_args(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter))
    ap.add_argument("--pclu-dir", default=str(R_DEFAULT / "results" / "pclu"))
    a = ap.parse_args()
    t0 = time.time()
    root = Path(a.root)
    res, figd = out_dirs(root)
    mod = load_pclu_module()
    pdir = Path(a.pclu_dir)
    pcs = strata(pdir)
    pfits = pl.read_csv(pdir / "curve_fits_overall.tsv", separator="\t").filter(
        pl.col("data") == "bin_medians_weighted_n")
    df = primary_1to1(load_tune_pairs(root, a.aln)).filter(pl.col("status") == "ok")
    log(f"X7 {a.aln}: {df.height} primary 1:1 status-ok tune pairs")
    windows = {"fident[0.10,1.0]": (0.10, 1.01), "fident[0.16,1.0] (Kimura-valid)": (0.16, 1.01),
               "fident[0.30,1.0]": (0.30, 1.01)}
    out = {"analysis": "X7", "aln": a.aln, "split": "tune", "root": str(root), "min_n_bin": a.min_n,
           "pclu_min_n_bin": mod.MIN_N_BIN, "bin_width": mod.BIN, "pclu_dir": str(pdir),
           "n_pairs": df.height, "mapping": {
               "fident": "1 - aa_p (core columns, this alignment); sensitivity aa_pid_seq (seqonly, all columns)",
               "tm": "tm_short (qTM, normalised by shorter protein) vs P-clu alntmscore (normalised by aln length)",
               "lddt": "lddt_core (primary) / lddt_global (sensitivity) vs P-clu Foldseek lddt (aligned residues)"},
           "variants": {}}
    for xv in ["x=1-aa_p(core)", "x=aa_pid_seq"]:
        for lddt_col in ["lddt_core", "lddt_global"]:
            key = f"{xv}|lddt={lddt_col}"
            p = ortholog_frame(df, xv, lddt_col, mod.BIN)
            orth = mod.binned(p, min_n=a.min_n)
            dd = [diffs(orth, b, n) for n, b in pcs.items()]
            frows = fits_bin_medians(mod, orth, windows)
            side = []
            for r in frows:
                q = pfits.filter((pl.col("window") == r["window"]) & (pl.col("y") == r["y"]) &
                                 (pl.col("x") == r["x"]) & (pl.col("family") == r["family"]))
                side.append(dict(window=r["window"], y=r["y"], x=r["x"], family=r["family"],
                                 params_orth=r["params"], se_orth=r["se"], dAIC_orth=r.get("dAIC"), n_bins_orth=r["n_fit"],
                                 params_pclu=q["params"][0] if q.height else None,
                                 dAIC_pclu=float(q["dAIC"][0]) if q.height else None,
                                 n_bins_pclu=int(q["n_fit"][0]) if q.height else None))
            best = {}
            for r in side:
                k2 = f"{r['window']}|{r['y']}|{r['x']}"
                if r["dAIC_orth"] == 0:
                    best.setdefault(k2, {})["orth"] = r["family"]
                if r["dAIC_pclu"] == 0:
                    best.setdefault(k2, {})["pclu"] = r["family"]
            out["variants"][key] = {"ortholog_binned": orth.to_dicts(), "differences": dd,
                                    "fits_side_by_side": side, "best_family": best,
                                    "fident_range_orth": [float(p["fident"].min()), float(p["fident"].max())]}
            ov = dd[0]["summary"]
            log(f"  {key:40s} bins={orth.height} Δtm(overall)={ov['d_tm']['n_weighted_mean']} "
                f"Δlddt={ov['d_lddt']['n_weighted_mean']}")
            if lddt_col == "lddt_core":
                plot(orth, pcs, dd, xv, a.aln, lddt_col, figd / f"X7_{a.aln}_{xv.split('=')[1].replace('(core)', '_core').replace('1-', '1m')}.png")
    out["runtime_s"] = time.time() - t0
    write_json(out, res / f"X7_{a.aln}.json")
    log(f"X7 done in {out['runtime_s']:.1f}s")


if __name__ == "__main__":
    main()
