#!/usr/bin/env python3
"""P-clu analysis: structure-vs-sequence relationship inside AFDB clusters (v3 file 8).

Run:  python scripts/pclu/build_pairs.py   (stage 1, ~70 s; makes data/processed/pclu/pairs.parquet)
      python scripts/pclu/analyse_pclu.py  (stage 2: tables -> results/pclu/, figures -> figures/pclu/)
Seed: 20261030 (all sampling).  Environment: conda_envs/main (polars, scipy, matplotlib, numpy).

File 8 columns (no header; verified from the file + file 1): query=cluster REPRESENTATIVE, target=cluster MEMBER,
fident alnlen mismatch gapopen qstart qend tstart tend evalue bits lddt alntmscore.
Divergence quantities used (x axes), with p = 1 - fident:
  p      : p-distance
  dP     : Poisson-corrected distance, -ln(fident)
  dK     : Kimura protein distance, -ln(1 - p - 0.2 p^2)  (undefined for fident < ~0.146)
Responses: 1 - TM (alntmscore) and 1 - LDDT.
"""
import json, time, warnings
from pathlib import Path
import numpy as np
import polars as pl
from scipy.optimize import curve_fit
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
SEED = 20261030
ROOT = Path(__file__).resolve().parents[2]
PAIRS = ROOT / "data/processed/pclu/pairs.parquet"
RES = ROOT / "results/pclu"; FIG = ROOT / "figures/pclu"
RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
BIN = 0.02
MIN_N_BIN = 100          # minimum pairs for a bin median to be reported
MIN_N_BIN_STRAT = 300
OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#999999"]


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ---------------------------------------------------------------- loading
def load():
    df = pl.read_parquet(PAIRS)
    checks = {}
    s = df.filter(pl.col("is_self"))
    checks["self_pairs_n"] = s.height
    checks["self_fident_min"] = float(s["fident"].min())
    checks["self_tm_min"] = float(s["tm"].min())
    checks["self_tm_mean"] = float(s["tm"].mean())
    checks["self_lddt_min"] = float(s["lddt"].drop_nulls().min())
    checks["self_lddt_frac_ge_0.999"] = float((s["lddt"] >= 0.999).mean())
    checks["self_alnlen_ne_repLen"] = int((s["alnlen"] != s["repLen"]).sum())
    p = df.filter(~pl.col("is_self"))
    checks["pairs_n"] = p.height
    # fident bookkeeping: fident <= 1 - mismatch/alnlen must hold
    checks["frac_fident_gt_1_minus_mm_over_aln(+1e-3)"] = float(
        (p["fident"] > 1 - p["mismatch"] / p["alnlen"] + 1e-3).mean())
    checks["n_pairs_fident_eq_1"] = int((p["fident"] >= 0.9995).sum())
    checks["n_clusters"] = int(p["query"].n_unique())
    p = p.with_columns(
        (1 - pl.col("fident")).alias("p"),
        (-pl.col("fident").clip(1e-6, 1).log()).alias("dP"),
        (1 - pl.col("tm")).alias("y_tm"),
        (1 - pl.col("lddt")).alias("y_lddt"),
    )
    pp = p["p"].to_numpy()
    arg = 1 - pp - 0.2 * pp ** 2
    dK = np.where(arg > 1e-3, -np.log(np.clip(arg, 1e-3, None)), np.nan)
    p = p.with_columns(pl.Series("dK", dK))
    checks["frac_pairs_kimura_undefined(fident<0.146)"] = float(np.isnan(dK).mean())
    # ---- strata labels
    dm, dr = p["dom_member"].to_numpy().astype(str), p["dom_rep"].to_numpy().astype(str)
    short = {"Bacteria": "Bac", "Archaea": "Arc", "Eukaryota": "Euk", "Viruses": "Vir", "Other": "Unc"}
    pd_ = np.array(["-".join(sorted([short[a], short[b]])) for a, b in zip(dm[:0], dr[:0])])  # placeholder
    pair_dom = np.empty(len(dm), dtype=object)
    sm = np.vectorize(short.get)(dm); sr = np.vectorize(short.get)(dr)
    lo = np.where(sm <= sr, sm, sr); hi = np.where(sm <= sr, sr, sm)
    pair_dom = np.char.add(np.char.add(lo.astype(str), "-"), hi.astype(str))
    known = np.isin(sm, ["Bac", "Arc", "Euk"]) & np.isin(sr, ["Bac", "Arc", "Euk"])
    pair_dom = np.where(known, pair_dom, "other(unclassified/virus)")
    sd = np.where(~known, "unknown", np.where(sm == sr, "same_superkingdom", "different_superkingdom"))
    lca_rank = p["lca_rank"].to_numpy().astype(str)
    age = np.select(
        [np.isin(lca_rank, ["species", "genus"]), np.isin(lca_rank, ["family", "order", "class"]),
         np.isin(lca_rank, ["phylum", "kingdom"]), lca_rank == "domain", lca_rank == "root_or_cellular"],
        ["1 species/genus", "2 family-class", "3 phylum-kingdom", "4 one superkingdom", "5 cross-superkingdom (root/cellular)"],
        "other")
    def binlab(x, edges, fmt):
        idx = np.digitize(x, edges)
        labs = [f"<{edges[0]}"] + [f"{edges[i]}-{edges[i+1]}" for i in range(len(edges) - 1)] + [f">={edges[-1]}"]
        return np.array([f"{i} {l}" for i, l in enumerate(labs)])[idx]
    p = p.with_columns(
        pl.Series("pair_superkingdoms", pair_dom), pl.Series("same_superkingdom", sd),
        pl.Series("lca_age_class", age),
        pl.Series("dark", np.where(p["isDark"].to_numpy(), "dark", "annotated")),
        pl.Series("repPlddt_bin", binlab(p["repPlddt"].to_numpy(), [60, 70, 80, 90], "")),
        pl.Series("repLen_bin", binlab(p["repLen"].to_numpy(), [100, 200, 300, 500, 800], "")),
        pl.Series("nMem_bin", binlab(p["nMem"].to_numpy(), [5, 20, 100, 1000], "")),
        pl.Series("qcov_bin", binlab(np.clip(p["qcov"].to_numpy(), 0, 1.5), [0.5, 0.7, 0.9, 0.95], "")),
        pl.Series("evalue_bin", binlab(np.log10(np.clip(p["evalue"].to_numpy(), 1e-300, None)), [-30, -10, -5, -3, -2, 0], "")),
        (pl.col("fident") / BIN).floor().clip(0, 49).cast(pl.Int32).alias("fbin"),
    )
    return p, checks


# ---------------------------------------------------------------- binned curves
def binned(p, by=None, ycols=("tm", "lddt"), min_n=MIN_N_BIN):
    g = ["fbin"] + ([by] if by else [])
    aggs = [pl.len().alias("n")]
    for y in ycols:
        aggs += [pl.col(y).median().alias(f"{y}_med"), pl.col(y).quantile(0.25).alias(f"{y}_q25"),
                 pl.col(y).quantile(0.75).alias(f"{y}_q75")]
    aggs += [pl.col("fident").median().alias("fident_med"), pl.col("dP").median().alias("dP_med")]
    out = p.group_by(g).agg(aggs).filter(pl.col("n") >= min_n).sort(g)
    return out.with_columns((pl.col("fbin") * BIN + BIN / 2).alias("fident_center"))


# ---------------------------------------------------------------- curve families
FAM = {
    "linear":        (lambda x, a, b: a + b * x,                              2, [0.1, 0.5]),
    "linear_origin": (lambda x, b: b * x,                                     1, [0.5]),
    "saturating":    (lambda x, A, b: A * (1 - np.exp(-b * x)),                2, [0.8, 1.0]),
    "saturating_offset": (lambda x, y0, A, b: y0 + A * (1 - np.exp(-b * x)),   3, [0.05, 0.8, 1.0]),
    "power":         (lambda x, a, c: a * np.power(np.clip(x, 1e-9, None), c), 2, [0.5, 0.7]),
    "exponential_CL": (lambda x, a, b: a * np.exp(b * x),                      2, [0.1, 2.0]),  # Chothia-Lesk form
    # plateau followed by linear rise (hockey stick): y0 + s*max(0, x - xc)
    "broken_stick":  (lambda x, y0, s_, xc: y0 + s_ * np.clip(x - xc, 0, None), 3, [0.17, 0.3, 0.8]),
    "hill_sigmoid":  (lambda x, A, K, h: A * np.power(np.clip(x, 1e-9, None), h) / (K ** h + np.power(np.clip(x, 1e-9, None), h)), 3, [0.9, 1.5, 2.0]),
}


def fit_family(name, x, y, w=None, start=None):
    f, k, p0 = FAM[name]
    best = None
    starts = [start] if start is not None else [p0] + [[v * s for v in p0] for s in (0.3, 3.0)]
    for st in starts:
        try:
            sig = None if w is None else 1 / np.sqrt(w)
            popt, pcov = curve_fit(f, x, y, p0=st, sigma=sig, maxfev=20000)
            r = y - f(x, *popt)
            rss = float(np.sum((r ** 2) * (1 if w is None else w)))
            if best is None or rss < best[1]:
                best = (popt, rss, pcov)
        except Exception:
            continue
    if best is None:
        return None
    popt, rss, pcov = best
    n = len(x)
    rss_n = rss / (1 if w is None else w.mean())  # weights normalised to mean 1 (effective N = bins)
    aic = n * np.log(rss_n / n) + 2 * k
    bic = n * np.log(rss_n / n) + k * np.log(n)
    se = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else np.full(k, np.nan)
    return dict(family=name, k=k, params=";".join(f"{v:.5g}" for v in popt), se=";".join(f"{v:.3g}" for v in se),
                rss=rss_n, rmse=float(np.sqrt(rss_n / n)), aic=float(aic), bic=float(bic), n_fit=n)


def fit_table(x, y, w, label, fams=FAM, starts=None):
    rows = []
    for fn in fams:
        r = fit_family(fn, x, y, w, None if starts is None else starts.get(fn))
        if r:
            r.update(label); rows.append(r)
    if rows:
        m = min(r["aic"] for r in rows)
        for r in rows:
            r["dAIC"] = r["aic"] - m
    return rows


def run_fits(p, bins_all, sample):
    """Fits of 1-TM / 1-LDDT vs p, dP, dK on (a) bin medians weighted by n; (b) 1M-row random sample."""
    rows = []; starts = {}
    windows = {"fident[0.10,1.0]": (0.10, 1.01), "fident[0.16,1.0] (Kimura-valid)": (0.16, 1.01),
               "fident[0.10,0.50]": (0.10, 0.50), "fident[0.30,1.0]": (0.30, 1.01)}
    for wname, (lo, hi) in windows.items():
        for ycol, ylab in [("tm", "1-TM"), ("lddt", "1-LDDT")]:
            b = bins_all.filter((pl.col("fident_center") >= lo) & (pl.col("fident_center") < hi))
            fm = b["fident_med"].to_numpy()
            y = 1 - b[f"{ycol}_med"].to_numpy(); n = b["n"].to_numpy().astype(float)
            for xname, x in [("p=1-fident", 1 - fm), ("dP=-ln(fident)", -np.log(fm)),
                             ("dK(Kimura)", -np.log(1 - (1 - fm) - 0.2 * (1 - fm) ** 2))]:
                ok = np.isfinite(x) & np.isfinite(y)
                if xname.startswith("dK") and lo < 0.16:
                    continue
                tb = fit_table(x[ok], y[ok], n[ok] , dict(data="bin_medians_weighted_n", window=wname, y=ylab, x=xname))
                rows += tb
                starts[(wname, ylab, xname)] = {r["family"]: [float(v) for v in r["params"].split(";")] for r in tb}
            if wname not in SAMPLE_WINDOWS:
                continue
            # sample
            s = sample.filter((pl.col("fident") >= lo) & (pl.col("fident") < hi)).drop_nulls(["lddt"])
            ys = s["y_tm"].to_numpy() if ycol == "tm" else s["y_lddt"].to_numpy()
            for xname, xs in [("p=1-fident", s["p"].to_numpy()), ("dP=-ln(fident)", s["dP"].to_numpy()),
                              ("dK(Kimura)", s["dK"].to_numpy())]:
                ok = np.isfinite(xs) & np.isfinite(ys)
                if xname.startswith("dK") and lo < 0.16:
                    ok &= (s["fident"].to_numpy() >= 0.16)  # restrict to Kimura-valid rows
                if xname.startswith("dK") and lo < 0.16:
                    pass
                rows += fit_table(xs[ok], ys[ok], None, dict(data="random_1M_rows", window=wname, y=ylab, x=xname),
                                  starts=starts.get((wname, ylab, xname)))
    return pl.DataFrame(rows)


# ---------------------------------------------------------------- stratification
STRATA = ["same_superkingdom", "pair_superkingdoms", "pair_mrca", "lca_rank", "lca_age_class", "dark",
          "repPlddt_bin", "repLen_bin", "nMem_bin", "qcov_bin", "evalue_bin"]
SAMPLE_WINDOWS = ["fident[0.10,1.0]", "fident[0.16,1.0] (Kimura-valid)"]
FWIN = [(0.10, 0.20), (0.20, 0.30), (0.30, 0.40), (0.40, 0.60), (0.60, 1.01)]


def stratum_summary(p):
    out = []
    for sv in STRATA:
        agg = [pl.len().alias("n"), pl.col("fident").median().alias("fident_med"),
               pl.col("tm").median().alias("tm_med"), pl.col("lddt").median().alias("lddt_med"),
               pl.col("tm").quantile(0.25).alias("tm_q25"), pl.col("tm").quantile(0.75).alias("tm_q75"),
               (pl.col("tm") < 0.5).mean().alias("frac_tm_lt_0.5"),
               (pl.col("evalue") > 0.01).mean().alias("frac_E_gt_0.01"),
               (pl.col("qcov") < 0.9).mean().alias("frac_qcov_lt_0.9")]
        for lo, hi in FWIN:
            m = (pl.col("fident") >= lo) & (pl.col("fident") < hi)
            agg += [pl.col("tm").filter(m).median().alias(f"tm_med@f{lo:.1f}-{min(hi,1.0):.1f}"),
                    pl.col("lddt").filter(m).median().alias(f"lddt_med@f{lo:.1f}-{min(hi,1.0):.1f}"),
                    m.sum().alias(f"n@f{lo:.1f}-{min(hi,1.0):.1f}")]
        t = p.group_by(sv).agg(agg).sort(sv).rename({sv: "level"}).with_columns(pl.lit(sv).alias("stratum_var"))
        out.append(t)
    return pl.concat(out, how="diagonal").select(["stratum_var", "level"] + [c for c in out[0].columns if c not in ("level", "stratum_var")])


def stratum_fits(p):
    """Per stratum level: best family (AIC on weighted bin medians) for 1-TM and 1-LDDT vs dP on fident>=0.10, + linear slope."""
    rows = []
    for sv in STRATA:
        b = binned(p.filter(pl.col("fident") >= 0.10), by=sv, min_n=MIN_N_BIN_STRAT)
        for lev in b[sv].unique().sort().to_list():
            bl = b.filter(pl.col(sv) == lev)
            if bl.height < 8:
                continue
            fm = bl["fident_med"].to_numpy(); n = bl["n"].to_numpy().astype(float)
            for ycol, ylab in [("tm", "1-TM"), ("lddt", "1-LDDT")]:
                y = 1 - bl[f"{ycol}_med"].to_numpy()
                x = -np.log(fm)
                ft = fit_table(x, y, n, dict(stratum_var=sv, level=lev, y=ylab, x="dP", n_bins=bl.height),
                               fams=["linear", "saturating", "power", "exponential_CL"])
                if not ft:
                    continue
                lin = [r for r in ft if r["family"] == "linear"]
                best = min(ft, key=lambda r: r["aic"])
                rows.append(dict(stratum_var=sv, level=lev, y=ylab, x="dP", n_bins=bl.height,
                                 best_family=best["family"], best_params=best["params"],
                                 linear_params=(lin[0]["params"] if lin else ""),
                                 **{f"dAIC_{r['family']}": r["dAIC"] for r in ft}))
    return pl.DataFrame(rows)


# ---------------------------------------------------------------- truncation
def truncation(p):
    f = lambda e: float(p.select(e.mean())[0, 0])
    nominal = (pl.col("qcov") >= 0.9) & (pl.col("evalue") < 0.01)
    d = {
        "n_pairs_non_self": p.height,
        "frac_E_gt_0.01": f(pl.col("evalue") > 0.01), "frac_E_gt_1e-3": f(pl.col("evalue") > 1e-3),
        "frac_E_in_[1e-3,1e-2]": f((pl.col("evalue") >= 1e-3) & (pl.col("evalue") <= 1e-2)),
        "frac_E_gt_1": f(pl.col("evalue") > 1), "frac_E_gt_10": f(pl.col("evalue") > 10),
        "frac_qcov_lt_0.9": f(pl.col("qcov") < 0.9), "frac_qcov_lt_0.5": f(pl.col("qcov") < 0.5),
        "frac_qcov_in_[0.9,0.95)": f((pl.col("qcov") >= 0.9) & (pl.col("qcov") < 0.95)),
        "frac_pass_nominal_qcov>=0.9_&_E<0.01": f(nominal),
        "frac_fident_in_[0.45,0.55)": f((pl.col("fident") >= 0.45) & (pl.col("fident") < 0.55)),
        "frac_fident_ge_0.5": f(pl.col("fident") >= 0.5), "frac_fident_ge_0.7": f(pl.col("fident") >= 0.7),
        "frac_fident_ge_0.9": f(pl.col("fident") >= 0.9), "frac_fident_lt_0.1": f(pl.col("fident") < 0.1),
        "frac_fident_lt_0.15": f(pl.col("fident") < 0.15),
        "frac_tm_lt_0.3": f(pl.col("tm") < 0.3), "frac_tm_lt_0.5": f(pl.col("tm") < 0.5),
        "frac_lddt_lt_0.5": f(pl.col("lddt") < 0.5),
        "fident_quantiles_1_5_25_50_75_95_99": [float(v) for v in np.quantile(p["fident"].to_numpy(), [.01, .05, .25, .5, .75, .95, .99])],
        "tm_quantiles_1_5_25_50_75_95_99": [float(v) for v in np.quantile(p["tm"].to_numpy(), [.01, .05, .25, .5, .75, .95, .99])],
    }
    # fine histogram of fident around 0.5 (0.01 bins) to look for a step at the AFDB50 threshold
    h = p.with_columns((pl.col("fident") / 0.01).floor().cast(pl.Int32).alias("b")).group_by("b").len().sort("b")
    hv = dict(zip(h["b"].to_list(), h["len"].to_list()))
    d["n_per_0.01_fident_bin_0.40-0.60"] = {f"{b/100:.2f}": hv.get(b, 0) for b in range(40, 60)}
    by = (p.group_by("fbin").agg(pl.len().alias("n"), (pl.col("qcov") < 0.9).mean().alias("frac_qcov_lt_0.9"),
                                  (pl.col("evalue") > 0.01).mean().alias("frac_E_gt_0.01"),
                                  (pl.col("tm") < 0.3).mean().alias("frac_tm_lt_0.3"),
                                  (pl.col("tm") < 0.5).mean().alias("frac_tm_lt_0.5"),
                                  pl.col("qcov").median().alias("qcov_med"),
                                  pl.col("evalue").median().alias("evalue_med"),
                                  pl.col("tm").median().alias("tm_med"),
                                  pl.col("tm").quantile(0.05).alias("tm_q05"))
          .sort("fbin").with_columns(((pl.col("fbin") * BIN)).alias("fident_lo")))
    return d, by


# ---------------------------------------------------------------- figures
def style():
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 100, "savefig.dpi": 300})


def fig_hexbin(p, bins_all):
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.2), constrained_layout=True)
    x = p["fident"].to_numpy()
    for ax, ycol, lab in zip(axs, ["tm", "lddt"], ["Foldseek alignment TM-score (alntmscore)", "Foldseek LDDT"]):
        y = p[ycol].to_numpy()
        ok = np.isfinite(y)
        hb = ax.hexbin(x[ok], y[ok], gridsize=(90, 60), extent=(0, 1, 0, 1), bins="log", mincnt=1, cmap="viridis", linewidths=0)
        c = bins_all
        ax.fill_between(c["fident_center"], c[f"{ycol}_q25"], c[f"{ycol}_q75"], color="white", alpha=0.35, lw=0)
        ax.plot(c["fident_center"], c[f"{ycol}_med"], color="#D55E00", lw=2.2, label="median (0.02 bins)")
        ax.plot(c["fident_center"], c[f"{ycol}_q25"], color="#D55E00", lw=1, ls="--", label="IQR")
        ax.plot(c["fident_center"], c[f"{ycol}_q75"], color="#D55E00", lw=1, ls="--")
        ax.set_xlabel("Sequence identity of member vs representative (fident, Foldseek alignment)")
        ax.set_ylabel(lab); ax.set_xlim(0, 1); ax.set_ylim(0, 1.0)
        cb = fig.colorbar(hb, ax=ax); cb.set_label("pairs per hexagon (log)")
        ax.legend(loc="lower right", frameon=True)
    fig.suptitle(f"AFDB clusters v3 (file 8): {p.height/1e6:.1f}M member-representative pairs (self-pairs excluded)")
    fig.savefig(FIG / "pclu_hexbin_fident_vs_tm_lddt.png"); plt.close(fig)


def fig_fits(bins_all, fits):
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    cols = {"linear": OI[0], "saturating": OI[1], "power": OI[2], "exponential_CL": OI[3], "hill_sigmoid": OI[4], "broken_stick": OI[5]}
    w = "fident[0.10,1.0]"
    b = bins_all.filter(pl.col("fident_center") >= 0.10)
    fm = b["fident_med"].to_numpy(); n = b["n"].to_numpy()
    for i, (ycol, ylab) in enumerate([("tm", "1 - TM"), ("lddt", "1 - LDDT")]):
        yv = 1 - b[f"{ycol}_med"].to_numpy()
        iqr_lo = 1 - b[f"{ycol}_q75"].to_numpy(); iqr_hi = 1 - b[f"{ycol}_q25"].to_numpy()
        for j, (xn, xlab, xv) in enumerate([("p=1-fident", "p = 1 - fident", 1 - fm),
                                            ("dP=-ln(fident)", "Poisson distance  -ln(fident)", -np.log(fm)),
                                            ("dK(Kimura)", "Kimura distance (fident >= 0.16)", None)]):
            ax = axs[i, j]
            if xv is None:
                msk = fm >= 0.16
                xv_ = -np.log(1 - (1 - fm[msk]) - 0.2 * (1 - fm[msk]) ** 2); yy = yv[msk]; lo_ = iqr_lo[msk]; hi_ = iqr_hi[msk]; ww = "fident[0.16,1.0] (Kimura-valid)"
            else:
                xv_, yy, lo_, hi_, ww = xv, yv, iqr_lo, iqr_hi, w
            ax.fill_between(xv_, lo_, hi_, color="#999999", alpha=0.25, lw=0, label="IQR of pairs")
            ax.scatter(xv_, yy, s=np.clip(8 + 40 * np.sqrt(n[:len(xv_)] if xv is not None else n[msk]) / np.sqrt(n.max()), 6, 50),
                       color="k", zorder=3, label="bin median (0.02 fident bins)")
            sub = fits.filter((pl.col("data") == "bin_medians_weighted_n") & (pl.col("window") == ww) &
                              (pl.col("y") == ("1-TM" if ycol == "tm" else "1-LDDT")) & (pl.col("x") == xn))
            xs = np.linspace(0, xv_.max() * 1.02, 200)
            for r in sub.iter_rows(named=True):
                if r["family"] not in cols:
                    continue
                prm = [float(v) for v in r["params"].split(";")]
                ax.plot(xs, FAM[r["family"]][0](xs, *prm), color=cols[r["family"]], lw=1.8,
                        label=f"{r['family']} (dAIC={r['dAIC']:.0f})")
            ax.set_xlabel(xlab); ax.set_ylabel(ylab); ax.set_ylim(0, 1.0); ax.set_xlim(0, xs.max())
            if i == 0 and j == 0:
                ax.legend(fontsize=8, loc="lower right")
            else:
                ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("Structural divergence vs sequence divergence in AFDB clusters: curve families fitted to n-weighted bin medians")
    fig.savefig(FIG / "pclu_curve_fits.png"); plt.close(fig)


def fig_strat(p):
    specs = [("same_superkingdom", "Pair superkingdoms"), ("lca_age_class", "Cluster LCA rank (age class)"),
             ("dark", "Cluster dark flag"), ("repPlddt_bin", "Representative pLDDT"), ("repLen_bin", "Representative length (aa)")]
    fig, axs = plt.subplots(2, len(specs), figsize=(4.2 * len(specs), 7.6), constrained_layout=True, sharex=True, sharey="row")
    for j, (sv, title) in enumerate(specs):
        b = binned(p, by=sv, min_n=MIN_N_BIN_STRAT)
        levels = b[sv].unique().sort().to_list()
        levels = [l for l in levels if not l.startswith("unknown")]
        cmap = plt.get_cmap("viridis")
        for i, (ycol, ylab) in enumerate([("tm", "median TM-score"), ("lddt", "median LDDT")]):
            ax = axs[i, j]
            for k, lev in enumerate(levels):
                bl = b.filter(pl.col(sv) == lev)
                if sv in ("lca_age_class", "repPlddt_bin", "repLen_bin"):
                    col = cmap(k / max(len(levels) - 1, 1) * 0.9)
                else:
                    col = OI[k % len(OI)]
                ax.plot(bl["fident_center"], bl[f"{ycol}_med"], color=col, lw=1.8, label=f"{lev}")
                ax.fill_between(bl["fident_center"], bl[f"{ycol}_q25"], bl[f"{ycol}_q75"], color=col, alpha=0.08, lw=0)
            ax.set_ylim(0.2, 1.0); ax.set_xlim(0, 1)
            if i == 1:
                ax.set_xlabel("fident (member vs representative)")
            if j == 0:
                ax.set_ylabel(ylab)
            if i == 0:
                ax.set_title(title, fontsize=10)
                ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("Stratified median curves (line = median, band = IQR; bins of 0.02 fident, >=300 pairs/bin/stratum)")
    fig.savefig(FIG / "pclu_stratified_median_curves.png"); plt.close(fig)


def fig_strat2(p):
    """pair relationship (MRCA of the two species) and cluster-size / filter strata."""
    specs = [("pair_superkingdoms", "Pair superkingdom combination"), ("pair_mrca", "Species-pair shared lineage rank"),
             ("qcov_bin", "Query (rep) alignment coverage"), ("evalue_bin", "log10 E-value of file-8 alignment"),
             ("nMem_bin", "Cluster size (members)")]
    fig, axs = plt.subplots(2, len(specs), figsize=(4.2 * len(specs), 7.6), constrained_layout=True, sharex=True, sharey="row")
    for j, (sv, title) in enumerate(specs):
        b = binned(p, by=sv, min_n=MIN_N_BIN_STRAT)
        levels = [l for l in b[sv].unique().sort().to_list()]
        if sv == "pair_superkingdoms":
            levels = [l for l in levels if not l.startswith("other")]
        cmap = plt.get_cmap("viridis")
        for i, (ycol, ylab) in enumerate([("tm", "median TM-score"), ("lddt", "median LDDT")]):
            ax = axs[i, j]
            for k, lev in enumerate(levels):
                bl = b.filter(pl.col(sv) == lev)
                col = OI[k % len(OI)] if sv == "pair_superkingdoms" else cmap(k / max(len(levels) - 1, 1) * 0.9)
                ax.plot(bl["fident_center"], bl[f"{ycol}_med"], color=col, lw=1.6, label=lev)
            ax.set_ylim(0.0, 1.0); ax.set_xlim(0, 1)
            if i == 1:
                ax.set_xlabel("fident")
            if j == 0:
                ax.set_ylabel(ylab)
            if i == 0:
                ax.set_title(title, fontsize=10); ax.legend(fontsize=6.5, loc="lower right")
    fig.suptitle("Further strata: species-pair relationship, alignment filters and cluster size (medians)")
    fig.savefig(FIG / "pclu_stratified_median_curves_2.png"); plt.close(fig)


def fig_trunc(by, p):
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6), constrained_layout=True)
    c = by["fident_lo"].to_numpy() + BIN / 2
    ax = axs[0]
    ax.bar(c, by["n"].to_numpy() / 1e6, width=BIN * 0.95, color=OI[0])
    ax.set_xlabel("fident"); ax.set_ylabel("pairs per 0.02 bin (millions)"); ax.set_title("Pair density: truncated at high identity")
    ax.axvline(0.5, color="k", ls=":", lw=1); ax.text(0.51, ax.get_ylim()[1] * 0.9, "AFDB50 (50% id)", fontsize=8)
    ax = axs[1]
    ax.plot(c, by["frac_qcov_lt_0.9"], color=OI[1], label="query coverage < 0.9")
    ax.plot(c, by["frac_E_gt_0.01"], color=OI[2], label="E-value > 0.01")
    ax.plot(c, by["frac_tm_lt_0.3"], color=OI[3], label="TM < 0.3")
    ax.plot(c, by["frac_tm_lt_0.5"], color=OI[4], label="TM < 0.5")
    ax.set_xlabel("fident"); ax.set_ylabel("fraction of pairs"); ax.set_title("Pairs violating nominal clustering criteria / TM floors"); ax.legend(fontsize=8)
    ax = axs[2]
    ax.plot(c, by["tm_med"], color="k", label="median TM")
    ax.plot(c, by["tm_q05"], color="k", ls="--", label="5th percentile TM")
    ax.plot(c, by["qcov_med"], color=OI[1], label="median query coverage")
    ax.set_xlabel("fident"); ax.set_ylim(0, 1.05); ax.set_title("Low-identity tail: TM / coverage vs fident"); ax.legend(fontsize=8)
    fig.savefig(FIG / "pclu_selection_truncation.png"); plt.close(fig)


def fig_sens(sens):
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.6), constrained_layout=True, sharex=True)
    for ax, ycol, lab in zip(axs, ["tm", "lddt"], ["median TM-score", "median LDDT"]):
        for k, (name, b) in enumerate(sens.items()):
            ax.plot(b["fident_center"], b[f"{ycol}_med"], color=OI[k], lw=1.8, label=f"{name} (n={int(b['n'].sum())/1e6:.1f}M)")
        ax.set_xlabel("fident"); ax.set_ylabel(lab); ax.set_xlim(0, 1); ax.set_ylim(0.2, 1.0); ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("Sensitivity of the median curves to cluster size weighting and to filtering")
    fig.savefig(FIG / "pclu_sensitivity_curves.png"); plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    t0 = time.time()
    p, checks = load()
    log("loaded", p.height)
    rho = {}
    smp = p.sample(1_000_000, seed=SEED)
    for y in ["tm", "lddt"]:
        s = smp.drop_nulls(y)
        rho[f"spearman_fident_vs_{y}_1Msample"] = float(spearmanr(s["fident"].to_numpy(), s[y].to_numpy())[0])
        rho[f"pearson_fident_vs_{y}_1Msample"] = float(np.corrcoef(s["fident"].to_numpy(), s[y].to_numpy())[0, 1])
    rho["spearman_tm_vs_lddt_1Msample"] = float(spearmanr(smp["tm"].to_numpy(), smp["lddt"].to_numpy(), nan_policy="omit")[0])
    # --- overall binned curves
    bins_all = binned(p)
    bins_all.write_csv(RES / "binned_curve_overall.tsv", separator="\t", float_precision=5)
    log("binned overall", bins_all.height)
    # cluster-balanced sensitivity: at most 5 random pairs per cluster
    bal = (p.with_columns(pl.int_range(pl.len()).shuffle(seed=SEED).alias("rnd"))
             .with_columns(pl.col("rnd").rank("ordinal").over("query").alias("rk")).filter(pl.col("rk") <= 5))
    binned(bal).write_csv(RES / "binned_curve_cluster_balanced_max5.tsv", separator="\t", float_precision=5)
    checks["n_pairs_cluster_balanced_max5"] = bal.height
    nom = p.filter((pl.col("qcov") >= 0.9) & (pl.col("evalue") < 0.01))
    checks["n_pairs_nominal_pass"] = nom.height
    bn = binned(nom); bn.write_csv(RES / "binned_curve_nominal_pass_qcov0.9_E0.01.tsv", separator="\t", float_precision=5)
    run_fits(nom, bn, nom.head(1000).clone()).filter(pl.col("data") == "bin_medians_weighted_n").write_csv(
        RES / "curve_fits_nominal_pass.tsv", separator="\t", float_precision=6)
    # LDDT-based sensitivity and per-pLDDT-class high-confidence subset (rep pLDDT>=90)
    hc = p.filter(pl.col("repPlddt") >= 90)
    bh = binned(hc); bh.write_csv(RES / "binned_curve_repPlddt_ge90.tsv", separator="\t", float_precision=5)
    run_fits(hc, bh, hc.head(1000).clone()).filter(pl.col("data") == "bin_medians_weighted_n").write_csv(
        RES / "curve_fits_repPlddt_ge90.tsv", separator="\t", float_precision=6)
    sens = {"all": bins_all, "cluster-balanced (<=5 pairs/cluster)": binned(bal),
            "qcov>=0.9 & E<0.01": bn, "rep pLDDT>=90": bh}
    # --- fits
    log("fitting")
    fits = run_fits(p, bins_all, smp)
    fits.write_csv(RES / "curve_fits_overall.tsv", separator="\t", float_precision=6)
    # fits on the cluster-balanced set (bin medians) for sensitivity
    fits_bal = run_fits(bal, binned(bal), bal.head(1000).clone()).filter(pl.col("data") == "bin_medians_weighted_n")
    fits_bal.write_csv(RES / "curve_fits_cluster_balanced_max5.tsv", separator="\t", float_precision=6)
    # --- stratification
    log("strata")
    stratum_summary(p).write_csv(RES / "stratum_summary.tsv", separator="\t", float_precision=4)
    stratum_fits(p).write_csv(RES / "stratum_curve_fits.tsv", separator="\t", float_precision=4)
    for sv in STRATA:
        binned(p, by=sv, min_n=MIN_N_BIN_STRAT).rename({sv: "level"}).with_columns(pl.lit(sv).alias("stratum_var")).write_csv(
            RES / f"binned_curve_by_{sv}.tsv", separator="\t", float_precision=5)
    # --- truncation
    log("truncation")
    d, by = truncation(p)
    by.write_csv(RES / "truncation_by_fident_bin.tsv", separator="\t", float_precision=5)
    (RES / "truncation_summary.json").write_text(json.dumps(d, indent=1))
    (RES / "sanity_checks.json").write_text(json.dumps({**checks, **rho}, indent=1))
    # --- figures
    log("figures")
    style()
    fig_hexbin(p, bins_all)
    fig_fits(bins_all, fits)
    fig_strat(p); fig_strat2(p)
    fig_trunc(by, p)
    fig_sens(sens)
    log("done in", f"{time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
