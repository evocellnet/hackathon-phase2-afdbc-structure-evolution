#!/usr/bin/env python3
"""Main-text figures A-E for the technical report (exploratory; ALL data = tune + confirm; needs CONFIRM_RUN=1).

  A  Goal-1 shape: 3Di-ML and lddt_D vs AA-ML (1:1 orthologs, E1/E2 filters), binned medians + IQR, confirmed
     E1/E2 curve eta + A x^c (A = pair-weighted median of per-HOG a_h at c fixed), linear reference (c = 1), floor eta.
  B  X7: ortholog pairs vs AFDB in-cluster pairs (P-clu), TM / LDDT vs identity (0.02 bins) + identity density.
  C  X2: path distances vs TimeTree MRCA node age, node-level medians across HOGs; AA/3Di ratio vs age.
  D  E3: per-HOG k_h distribution; pooled k by CATH class, age stratum, internal vs all edges (HOG bootstrap).
  E  E4/E5 caveat: pooled k by edge age bin, terminal vs internal edges.

Outputs: figures/main/fig{A..E}_*.png (300 dpi) + figures/main/fig{A..E}_*.tsv (plotted aggregates only; TimeTree
ages appear only as bin-level summaries, a bin's mean age only if it spans >= 3 nodes) + figures/main/key_numbers.json.
Run:  CONFIRM_RUN=1 python scripts/exploratory/figures_main.py [--only A,B,...] [--B 1000]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("POLARS_MAX_THREADS", "8")
import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from strucevo.stats.common import load_pairs, confirm_allowed  # noqa: E402
from strucevo.stats.endpoints_curve import prepare_curve_data  # noqa: E402
from strucevo.stats.pairwise import weighted_median  # noqa: E402
from xcommon import load_panel  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

OUT = ROOT / "figures" / "main"
OUT.mkdir(parents=True, exist_ok=True)
KEY = {}
SEED = 20261031

# palette (dataviz reference, light mode; <= 3 categorical slots per chart) and ink tokens
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7, "axes.edgecolor": INK2,
    "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "axes.axisbelow": True, "legend.frameon": False, "lines.linewidth": 2,
    "savefig.dpi": 300, "savefig.bbox": "tight", "figure.facecolor": "white"})


def logx(ax):
    from matplotlib.ticker import FixedLocator, NullFormatter
    ax.xaxis.set_major_locator(FixedLocator([3, 10, 30, 100, 300, 1000]))
    ax.set_xticklabels(["3", "10", "30", "100", "300", "1000"])
    ax.xaxis.set_minor_formatter(NullFormatter())


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def panel_label(ax, s):
    ax.text(-0.13, 1.04, s, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def tuning():
    import yaml
    return yaml.safe_load((ROOT / "config" / "tuning.yaml").read_text())


def ep(name):
    return json.loads((ROOT / "results" / "endpoints" / f"{name}.json").read_text())


def hog_cols():
    return pl.concat([pl.read_parquet(ROOT / "data" / "interim" / f"hogs_endpoint_cols_{s}.parquet")
                      for s in ("tune", "confirm")], how="diagonal_relaxed")


def all_pairs(aln, columns=None):
    return pl.concat([load_pairs(aln, s, columns=columns, pair_set="primary") for s in ("tune", "confirm")])


def binned(x, y, edges, min_n=30, extra=None):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    idx = np.digitize(x, edges) - 1
    rows = []
    for k in range(len(edges) - 1):
        m = idx == k
        n = int(m.sum())
        if n < min_n:
            continue
        q = np.quantile(y[m], [0.25, 0.5, 0.75])
        rows.append(dict(bin_lo=edges[k], bin_hi=edges[k + 1], n=n, x_median=float(np.median(x[m])),
                         y_q25=q[0], y_median=q[1], y_q75=q[2]))
    return pl.DataFrame(rows)


# ================================================================================================ Fig A
def per_hog_a(data, c, s2):
    """a_h by WLS of (y - eta) on x^c within HOG (w = 1/(sy^2 + s2)); pair-weighted median over HOGs."""
    el = data.w0 > 0
    g, x, y, sy = data.g[el], data.x[el], data.y[el] - data.eta, data.sy[el]
    w = 1.0 / (sy ** 2 + s2)
    xc = x ** c
    num = np.bincount(g, w * xc * y, data.H)
    den = np.bincount(g, w * xc * xc, data.H)
    a = num / np.where(den > 0, den, np.nan)
    npair = np.bincount(g, minlength=data.H).astype(float)
    return weighted_median(a, npair), a


def figA():
    T = tuning()
    E1, E2 = ep("E1"), ep("E2")
    chat = {("tdi_ml", "tdiaa"): E1["estimate"], ("tdi_ml", "seqonly"): E1["sensitivity"]["aln_seqonly"]["estimate"],
            ("tdi_ml", "tmalign"): E1["sensitivity"]["aln_tmalign"]["estimate"],
            ("lddt_D", "tdiaa"): E2["estimate"], ("lddt_D", "seqonly"): E2["sensitivity"]["aln_seqonly"]["estimate"],
            ("lddt_D", "tmalign"): E2["sensitivity"]["aln_tmalign"]["estimate"]}
    s2 = {"tdi_ml": E1["s2"], "lddt_D": E2["s2"]}
    fam = {"tdi_ml": "tdi", "lddt_D": "lddt"}
    hc = hog_cols()
    edges = np.round(np.arange(0.0, 2.51, 0.1), 2)
    res, rows = {}, []
    cols = ["split", "hog", "one_to_one", "status", "species_a", "species_b", "aa_ml_seq", "aa_ml", "aa_ml_se",
            "aa_sat", "tdi_ml", "tdi_ml_se", "tdi_sat", "lddt_D", "lddt_D_se", "lddt_sat"]
    for aln in ("tdiaa", "seqonly", "tmalign"):
        df = all_pairs(aln, cols)
        log(f"A {aln}: {df.height} primary pairs")
        for m in ("tdi_ml", "lddt_D"):
            ex = set(hc.filter(pl.col(f"excl_{fam[m]}_{aln}").fill_null(False))["hog"].to_list())
            data, info, idx, _ = prepare_curve_data(df, ycols=(m, f"{m}_se", m.split("_")[0] + "_sat"
                                                               if m == "tdi_ml" else "lddt_sat"), excluded_hogs=ex)
            data.eta = float(T["eta"][aln][m])
            c = chat[(m, aln)]
            A, _ = per_hog_a(data, c, s2[m])
            A1, _ = per_hog_a(data, 1.0, s2[m])
            el = data.w0 > 0
            b = binned(data.x[el], data.y[el], edges)
            res[(m, aln)] = dict(b=b, eta=data.eta, c=c, A=A, A1=A1, n=int(el.sum()), H=data.H,
                                 xmax=float(np.quantile(data.x[el], 0.995)))
            for r in b.iter_rows(named=True):
                rows.append(dict(metric=m, aln=aln, **r))
            KEY[f"A_{m}_{aln}"] = dict(c_confirmed=c, A=A, A_linear=A1, eta=data.eta, n_pairs=int(el.sum()),
                                       n_hogs=data.H)
            log(f"   {m}: n={el.sum()} H={data.H} c={c:.3f} A={A:.3f} A1={A1:.3f}")
    pl.DataFrame(rows).write_csv(OUT / "figA_shape_binned.tsv", separator="\t")

    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.0), constrained_layout=True)
    ylab = {"tdi_ml": "3Di ML distance (subst./site)", "lddt_D": "LDDT distance $D$ (floor-scaled)"}
    xl = "AA ML distance (subst./site)"
    col = {"tdiaa": C1, "seqonly": C2, "tmalign": C3}
    nm = {"tdiaa": "3Di+AA alignment", "seqonly": "sequence-only alignment", "tmalign": "TM-align"}
    for i, m in enumerate(("tdi_ml", "lddt_D")):
        ax = axs[i, 0]
        r = res[(m, "tdiaa")]
        b = r["b"]
        xs = np.linspace(0.001, 2.05, 300)
        ax.fill_between(b["x_median"], b["y_q25"], b["y_q75"], color=C1, alpha=0.18, lw=0, label="IQR of pairs")
        ax.plot(b["x_median"], b["y_median"], "o", color=C1, ms=4, mec="white", mew=0.6, label="bin median",
                zorder=4)
        ax.plot(xs, r["eta"] + r["A"] * xs ** r["c"], color=INK, lw=1.8,
                label=f"fit $\\eta + A x^{{\\hat c}}$, $\\hat c$ = {r['c']:.2f}")
        ax.plot(xs, r["eta"] + r["A1"] * xs, color=MUTED, lw=1.4, ls="--", label="linear reference ($c$ = 1)")
        ax.axhline(r["eta"], color=C2, lw=1.2, ls=":")
        ax.text(2.02, r["eta"], f"floor $\\eta$ = {r['eta']:.3f}", color=INK2, fontsize=7, ha="right", va="bottom")
        ax.set_xlim(0, 2.08)
        ax.set_ylim(0, None)
        ax.set_xlabel(xl); ax.set_ylabel(ylab[m])
        ax.set_title(f"{'3Di' if m == 'tdi_ml' else 'LDDT'} vs AA (3Di+AA alignment)", loc="left")
        ax.text(0.98, 0.16, f"{r['n']:,} pairs\n{r['H']:,} HOGs", transform=ax.transAxes, ha="right", fontsize=6.5, color=INK2)
        ax.legend(loc="upper left")
        panel_label(ax, "ab"[i * 1] if i == 0 else "c")
        ax2 = axs[i, 1]
        for aln in ("tdiaa", "seqonly", "tmalign"):
            r = res[(m, aln)]
            b = r["b"]
            ax2.plot(b["x_median"], b["y_median"], "o", color=col[aln], ms=3.5, mec="white", mew=0.5, zorder=4)
            ax2.plot(xs, r["eta"] + r["A"] * xs ** r["c"], color=col[aln], lw=1.6,
                     label=f"{nm[aln]}: $\\hat c$ = {r['c']:.2f}")
        ax2.set_xlim(0, 2.08); ax2.set_ylim(0, None)
        ax2.set_xlabel(xl); ax2.set_ylabel(ylab[m])
        ax2.set_title("Alignment sensitivity", loc="left")
        ax2.legend(loc="upper left")
        panel_label(ax2, "b" if i == 0 else "d")
    fig.savefig(OUT / "figA_shape.png")
    plt.close(fig)


# ================================================================================================ Fig B
def figB():
    pdir = ROOT / "results" / "pclu"
    po = pl.read_csv(pdir / "binned_curve_overall.tsv", separator="\t")
    pe = pl.read_csv(pdir / "binned_curve_by_pair_superkingdoms.tsv", separator="\t").filter(pl.col("level") == "Euk-Euk")
    df = all_pairs("tdiaa", ["split", "hog", "one_to_one", "status", "aa_p", "tm_short", "lddt_core"])
    df = df.filter((pl.col("one_to_one").cast(pl.Int64) == 1) & (pl.col("status") == "ok"))
    f = 1.0 - df["aa_p"].to_numpy()
    tm, ld = df["tm_short"].to_numpy(), df["lddt_core"].to_numpy()
    edges = np.round(np.arange(0, 1.0001, 0.02), 2)
    bt, bl = binned(f, tm, edges, 200), binned(f, ld, edges, 200)
    orth = bt.rename({"y_q25": "tm_q25", "y_median": "tm_med", "y_q75": "tm_q75"}).join(
        bl.select("bin_lo", pl.col("y_q25").alias("lddt_q25"), pl.col("y_median").alias("lddt_med"),
                  pl.col("y_q75").alias("lddt_q75")), on="bin_lo")
    ok = np.isfinite(f)
    hist_o, _ = np.histogram(f[ok], bins=edges)
    out = orth.with_columns(pl.lit("orthologs_1to1").alias("set"), ((pl.col("bin_lo") + pl.col("bin_hi")) / 2).alias("fident_center"))
    pclu_rows = []
    for nm_, d in (("pclu_overall", po), ("pclu_euk_euk", pe)):
        dd = d.filter(pl.col("n") >= 100).select("fident_center", "n", "tm_q25", "tm_med", "tm_q75", "lddt_q25",
                                                 "lddt_med", "lddt_q75").with_columns(pl.lit(nm_).alias("set"))
        pclu_rows.append(dd)
    tsv = pl.concat([out.select("set", "fident_center", "n", "tm_q25", "tm_med", "tm_q75", "lddt_q25", "lddt_med",
                                "lddt_q75")] + [p.select("set", "fident_center", "n", "tm_q25", "tm_med", "tm_q75",
                                                         "lddt_q25", "lddt_med", "lddt_q75") for p in pclu_rows])
    tsv.write_csv(OUT / "figB_selection_bias.tsv", separator="\t")

    # key numbers: medians at fixed identity
    def at(d, fc, col):
        r = d.filter((pl.col("fident_center") - fc).abs() < 1e-6)
        return float(r[col][0]) if r.height else None
    for fc in (0.21, 0.31, 0.51, 0.71):
        KEY[f"B_at_fident_{fc}"] = {s: {c: at(tsv.filter(pl.col("set") == s), fc, c) for c in ("tm_med", "lddt_med", "n")}
                                    for s in ("orthologs_1to1", "pclu_overall", "pclu_euk_euk")}
    KEY["B_orth_fident_quantiles"] = dict(zip(["q05", "q25", "q50", "q75", "q95"],
                                              np.quantile(f[ok], [.05, .25, .5, .75, .95]).tolist()))
    KEY["B_orth_frac_tm_lt_0.5"] = float(np.mean(tm[np.isfinite(tm)] < 0.5))
    KEY["B_orth_n"] = int(ok.sum())
    KEY["B_pclu_n"] = int(po["n"].sum())

    fig = plt.figure(figsize=(7.2, 5.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.3])
    sets = [("orthologs_1to1", C1, "1:1 orthologs (this study)"),
            ("pclu_overall", C2, "AFDB cluster pairs (all)"),
            ("pclu_euk_euk", C3, "AFDB cluster pairs (euk–euk)")]
    for j, (yc, ylabel) in enumerate((("tm", "TM-score"), ("lddt", "LDDT"))):
        ax = fig.add_subplot(gs[0, j])
        ax.axvspan(0, 0.10, color=GRID, alpha=0.7, lw=0)
        ax.text(0.05, 0.03, "cluster\nalignments\nunreliable", ha="center", va="bottom", fontsize=6, color=INK2)
        for s, c, lab in sets:
            d = tsv.filter(pl.col("set") == s)
            if s != "pclu_euk_euk":
                ax.fill_between(d["fident_center"], d[f"{yc}_q25"], d[f"{yc}_q75"], color=c, alpha=0.15, lw=0)
            ax.plot(d["fident_center"], d[f"{yc}_med"], color=c, lw=1.8 if s != "pclu_euk_euk" else 1.3,
                    ls="-" if s != "pclu_euk_euk" else "--", label=lab)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
        ax.set_xlabel("Sequence identity of aligned residues")
        ax.set_ylabel(f"{ylabel} (median, band = IQR)")
        ax.set_title(f"{ylabel} vs identity", loc="left")
        panel_label(ax, "ab"[j])
        if j == 0:
            ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.04), fontsize=6.5)
    ax = fig.add_subplot(gs[1, :])
    ctr = (edges[:-1] + edges[1:]) / 2
    dens_o = hist_o / hist_o.sum() / 0.02
    pp = po.sort("fident_center")
    dens_p = pp["n"].to_numpy() / pp["n"].sum() / 0.02
    ax.plot(ctr, dens_o, color=C1, lw=1.8, label="1:1 orthologs")
    ax.plot(pp["fident_center"], dens_p, color=C2, lw=1.8, label="AFDB cluster pairs")
    ax.set_xlim(0, 1); ax.set_ylim(0, None)
    ax.set_xlabel("Sequence identity of aligned residues"); ax.set_ylabel("Density of pairs")
    ax.set_title("Where the pairs are: identity distribution", loc="left")
    ax.legend(loc="upper right")
    panel_label(ax, "c")
    fig.savefig(OUT / "figB_selection_bias.png")
    plt.close(fig)


# ================================================================================================ Fig C
def figC():
    T = tuning()
    eta = T["eta"]["tdiaa"]
    hc = hog_cols()
    df = all_pairs("tdiaa", ["split", "hog", "one_to_one", "status", "species_a", "species_b", "aa_ml", "aa_sat",
                             "tdi_ml", "tdi_sat", "lddt_D", "lddt_sat"])
    df = df.filter((pl.col("one_to_one").cast(pl.Int64) == 1) & (pl.col("status") == "ok"))
    panel = load_panel(ROOT)
    from xcommon import mrca_age
    t = mrca_age(panel, df["species_a"].to_numpy(), df["species_b"].to_numpy())
    df = df.with_columns(pl.Series("age", t))
    log(f"C: {df.height} pairs; undated {np.isnan(t).sum()}")
    mets = {"aa_ml": ("aa_sat", "aa", 0.0), "tdi_ml": ("tdi_sat", "tdi", eta["tdi_ml"]),
            "lddt_D": ("lddt_sat", "lddt", eta["lddt_D"])}
    node = None
    for m, (sat, fam, e) in mets.items():
        ex = hc.filter(pl.col(f"excl_{fam}_tdiaa").fill_null(False))["hog"].to_list()
        d = df.filter(~pl.col(sat) & pl.col(m).is_finite() & ~pl.col("hog").is_in(ex) & pl.col("age").is_finite())
        d = d.with_columns((pl.col(m) - e).alias("y"))
        hn = d.group_by("age", "hog").agg(pl.col("y").median())
        nd = hn.group_by("age").agg(pl.col("y").median().alias(m), pl.len().alias(f"nh_{m}"))
        node = nd if node is None else node.join(nd, on="age", how="full", coalesce=True)
    # ratio per node: per-HOG means of AA and floor-corrected 3Di on the same pairs, ratio of medians across HOGs
    exa = set(hc.filter(pl.col("excl_aa_tdiaa").fill_null(False) | pl.col("excl_tdi_tdiaa").fill_null(False))["hog"].to_list())
    d = df.filter(~pl.col("aa_sat") & ~pl.col("tdi_sat") & pl.col("aa_ml").is_finite() & pl.col("tdi_ml").is_finite()
                  & ~pl.col("hog").is_in(list(exa)) & pl.col("age").is_finite())
    hn = d.group_by("age", "hog").agg(pl.col("aa_ml").mean().alias("a"), (pl.col("tdi_ml") - eta["tdi_ml"]).mean().alias("s"))
    rn = hn.group_by("age").agg((pl.col("a").sum() / pl.col("s").sum()).alias("ratio_aa_tdi"), pl.len().alias("nh_ratio"))
    node = node.join(rn, on="age", how="left").filter(pl.col("nh_aa_ml") >= 5).sort("age")
    # private node-level table (TimeTree ages) -> scratch only
    node.write_csv(Path(os.environ.get("SCRATCH_DIR", "/tmp")) / "figC_nodes_private.tsv", separator="\t")
    # committed aggregate: log-age bins
    bed = [0, 10, 30, 100, 300, 1000, 2000]
    rows = []
    a = node["age"].to_numpy()
    for lo, hi in zip(bed[:-1], bed[1:]):
        mm = (a >= lo) & (a < hi)
        sub = node.filter(pl.Series(mm))
        if sub.height == 0:
            continue
        r = dict(age_bin_lo=lo, age_bin_hi=hi, n_nodes=sub.height,
                 mean_node_age=float(sub["age"].mean()) if sub.height >= 3 else None)
        for m in ("aa_ml", "tdi_ml", "lddt_D", "ratio_aa_tdi"):
            r[f"{m}_median_of_nodes"] = float(sub[m].median())
        r["n_hogs_aa_total"] = int(sub["nh_aa_ml"].sum())
        rows.append(r)
    agg = pl.DataFrame(rows)
    agg.write_csv(OUT / "figC_time_binned.tsv", separator="\t")
    KEY["C_age_bins"] = agg.to_dicts()
    deep = node.filter(pl.col("age") >= 1000)
    shallow = node.filter(pl.col("age") < 1000)
    KEY["C_n_nodes"] = node.height
    KEY["C_ratio_median_nodes_lt1Ga"] = float(shallow["ratio_aa_tdi"].median())
    KEY["C_ratio_median_nodes_ge1Ga"] = float(deep["ratio_aa_tdi"].median())
    KEY["C_n_nodes_ge1Ga"] = deep.height
    # log-log slopes (node level, < 1 Ga, nodes >= 10 Mya)
    sl = {}
    for m in ("aa_ml", "tdi_ml", "lddt_D"):
        s = shallow.filter((pl.col("age") >= 10) & (pl.col(m) > 0))
        sl[m] = float(np.polyfit(np.log(s["age"].to_numpy()), np.log(s[m].to_numpy()), 1)[0])
    KEY["C_loglog_slope_10_1000Mya"] = sl

    fig, ax = plt.subplots(1, 1, figsize=(4.8, 3.6), constrained_layout=True)
    lab = {"aa_ml": "AA distance", "tdi_ml": "3Di distance − $\\eta$", "lddt_D": "LDDT distance − $\\eta$"}
    cc = {"aa_ml": C1, "tdi_ml": C2, "lddt_D": C3}
    for m in ("aa_ml", "tdi_ml", "lddt_D"):
        for sub, filled in ((shallow, True), (deep, False)):
            s = sub.filter(pl.col(m) > 0)
            ax.scatter(s["age"], s[m], s=10 + 2 * np.sqrt(s[f"nh_{m}"].to_numpy()),
                       facecolor=cc[m] if filled else "white", edgecolor=cc[m] if not filled else "white",
                       lw=0.8 if not filled else 0.3, alpha=0.9, zorder=3)
        gm = agg.filter(pl.col("age_bin_hi") <= 1000)
        xs = [np.sqrt(max(lo, 3) * hi) for lo, hi in zip(gm["age_bin_lo"], gm["age_bin_hi"])]
        ax.plot(xs, gm[f"{m}_median_of_nodes"], color=cc[m], lw=1.5, zorder=4,
                label=f"{lab[m]} (slope {sl[m]:.2f})")
    ax.axvspan(1000, 2000, color=GRID, alpha=0.6, lw=0)
    ax.text(1400, 0.012, "≥ 1 Ga\n(open)", ha="center", fontsize=6.5, color=INK2)
    xr = np.array([10, 1000])
    ax.plot(xr, 0.004 * xr / 10, color=MUTED, ls="--", lw=1, label="slope 1 (clock)")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(3, 2000); logx(ax)
    ax.set_xlim(3, 2000)
    ax.set_xlabel("Age of last common ancestor of the pair (Mya, TimeTree)")
    ax.set_ylabel("Ortholog distance (median over families)")
    ax.set_title("Distance vs divergence time", loc="left")
    ax.legend(loc="lower right", fontsize=6.5)
    fig.savefig(OUT / "figC_time.png")
    plt.close(fig)
    return
    ax = None
    E3 = ep("E3")
    for sub, filled in ((shallow, True), (deep, False)):
        ax.scatter(sub["age"], sub["ratio_aa_tdi"], s=10 + 2 * np.sqrt(sub["nh_ratio"].to_numpy()),
                   facecolor=C4 if filled else "white", edgecolor="white" if filled else C4, lw=0.8 if not filled else 0.3,
                   zorder=3)
    xs = [np.sqrt(max(lo, 3) * hi) for lo, hi in zip(agg["age_bin_lo"], agg["age_bin_hi"])]
    ax.plot(xs, agg["ratio_aa_tdi_median_of_nodes"], color=C4, lw=1.5, label="log-age bin median")
    ax.axhline(E3["estimate"], color=INK, lw=1, ls="--")
    ax.text(3.3, E3["estimate"], f" E3 branch estimate k = {E3['estimate']:.2f}", va="bottom", fontsize=6.5, color=INK2)
    ax.axvspan(1000, 2000, color=GRID, alpha=0.6, lw=0)
    ax.set_xscale("log"); ax.set_xlim(3, 2000); ax.set_ylim(0, None); logx(ax)
    ax.set_xlabel("Age of last common ancestor of the pair (Mya, TimeTree)")
    ax.set_ylabel("AA / (3Di − $\\eta$) path-distance ratio")
    ax.set_title("Sequence-to-structure ratio vs time", loc="left")
    ax.legend(loc="upper left", fontsize=6.5)
    panel_label(ax, "b")
    fig.savefig(OUT / "figC_time.png")
    plt.close(fig)


# ================================================================================================ branches (D, E)
def branch_edges():
    tau_m = float(tuning()["tau_misfit"])
    out = []
    for sub, split in (("", "confirm"), ("tune/", "tune")):
        E = pl.read_parquet(ROOT / "data" / "interim" / "branches" / f"{sub}tdiaa.parquet")
        F = pl.read_parquet(ROOT / "data" / "interim" / "branches" / f"{sub}tdiaa_fit.parquet")
        F = F.filter(pl.col("metric").is_in(["aa_ml", "tdi_ml"]))
        bad = F.filter((pl.col("reason").fill_null("") != "") | pl.col("excluded").fill_null(False)
                       | (pl.col("misfit_ratio").is_finite() & (pl.col("misfit_ratio") > tau_m)))["hog"].unique()
        ea = E.filter(pl.col("metric") == "aa_ml").select("hog", "edge_id", "is_terminal", "tau", pl.col("b").alias("b_aa"))
        et = E.filter(pl.col("metric") == "tdi_ml").select("hog", "edge_id", pl.col("b").alias("b_t"))
        e = ea.join(et, on=["hog", "edge_id"]).filter(~pl.col("hog").is_in(bad.to_list())).with_columns(pl.lit(split).alias("split"))
        out.append(e)
        log(f"branches {split}: {e['hog'].n_unique()} HOGs kept, {len(bad)} excluded")
    return pl.concat(out)


def boot_ratio(codes, SA, ST, B, rng, groups=None):
    """Ratio-of-sums k = sum SA / sum ST with HOG bootstrap (SA, ST per HOG; groups: per HOG label)."""
    H = len(SA)
    W = rng.multinomial(H, np.full(H, 1.0 / H), size=B).astype(float)
    res = {}
    labs = np.unique(groups) if groups is not None else ["all"]
    for g in labs:
        m = (groups == g) if groups is not None else np.ones(H, bool)
        k = SA[m].sum() / ST[m].sum()
        kb = (W[:, m] @ SA[m]) / (W[:, m] @ ST[m])
        res[g] = dict(k=float(k), lo=float(np.quantile(kb, 0.025)), hi=float(np.quantile(kb, 0.975)), n_hogs=int(m.sum()))
    return res


def figD(B):
    rng = np.random.default_rng(SEED)
    e = branch_edges()
    hogs = pl.read_parquet(ROOT / "data" / "interim" / "hogs.parquet").select("hog", "stratum")
    hc = hog_cols().select("hog", "cath_class")
    per = e.group_by("hog", "split").agg(pl.col("b_aa").sum().alias("SA"), pl.col("b_t").sum().alias("ST"),
                                         pl.col("b_aa").filter(~pl.col("is_terminal")).sum().alias("SA_int"),
                                         pl.col("b_t").filter(~pl.col("is_terminal")).sum().alias("ST_int"))
    per = per.join(hogs, on="hog", how="left").join(hc, on="hog", how="left").sort("hog")
    SA, ST = per["SA"].to_numpy(), per["ST"].to_numpy()
    conf = per["split"].to_numpy() == "confirm"
    KEY["D_check_confirm_k_recomputed"] = float(SA[conf].sum() / ST[conf].sum())
    KEY["D_E3_confirm_k"] = ep("E3")["estimate"]
    allk = boot_ratio(None, SA, ST, B, rng)["all"]
    intk = boot_ratio(None, per["SA_int"].to_numpy(), per["ST_int"].to_numpy(), B, rng)["all"]
    cath = boot_ratio(None, SA, ST, B, rng, per["cath_class"].fill_null("none").to_numpy())
    strat = boot_ratio(None, SA, ST, B, rng, per["stratum"].cast(pl.Utf8).fill_null("?").to_numpy())
    okh = (SA > 0) & (ST > 0)
    kh = SA[okh] / ST[okh]
    KEY["D_k_all"] = allk; KEY["D_k_internal"] = intk; KEY["D_k_cath"] = cath; KEY["D_k_stratum"] = strat
    KEY["D_kh_quantiles"] = dict(zip(["q10", "q25", "q50", "q75", "q90"], np.quantile(kh, [.1, .25, .5, .75, .9]).tolist()))
    KEY["D_kh_n"] = int(okh.sum()); KEY["D_kh_n_nonpos"] = int((~okh).sum()); KEY["D_n_hogs"] = len(SA)
    rows = [dict(panel="all_edges", group="all", **allk), dict(panel="internal_edges", group="internal", **intk)]
    rows += [dict(panel="cath", group=g, **v) for g, v in cath.items()]
    rows += [dict(panel="stratum", group=g, **v) for g, v in strat.items()]
    hb = np.arange(-1.0, 1.81, 0.05)
    hcnt, _ = np.histogram(np.log10(kh), hb)
    for lo, hi, c in zip(hb[:-1], hb[1:], hcnt):
        rows.append(dict(panel="kh_hist_log10", group=f"[{lo:.2f},{hi:.2f})", k=None, lo=None, hi=None, n_hogs=int(c)))
    pl.DataFrame(rows).write_csv(OUT / "figD_ratio_k.tsv", separator="\t")

    E3 = ep("E3")
    s3 = E3["sensitivity"]
    fig = plt.figure(figsize=(7.2, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.2, 1, 1], height_ratios=[2.2, 1])
    ax = fig.add_subplot(gs[:, 0])
    ax.hist(np.log10(kh), bins=hb, color=C1, edgecolor="white", lw=0.4)
    for v, ls, lab in ((np.log10(allk["k"]), "-", f"pooled k = {allk['k']:.2f}"),
                       (np.log10(np.median(kh)), "--", f"median per-family k = {np.median(kh):.2f}")):
        ax.axvline(v, color=INK, ls=ls, lw=1.2, label=lab)
    ax.set_xticks(np.log10([0.1, 0.3, 1, 3, 10, 30])); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30"])
    ax.set_xlim(-1, 1.8)
    ax.set_ylim(0, hcnt.max() * 1.3)
    ax.set_xlabel("k = AA / 3Di branch length (per family, log scale)"); ax.set_ylabel("Number of families")
    ax.set_title(f"Per-family ratio ({len(kh):,} families)", loc="left")
    ax.legend(loc="upper right", fontsize=6.5)
    panel_label(ax, "a")

    def dots(ax, items, title, lab):
        y = np.arange(len(items))[::-1]
        for yi, (name, v, e3v) in zip(y, items):
            ax.plot([v["lo"], v["hi"]], [yi, yi], color=C1, lw=2, solid_capstyle="round")
            ax.plot(v["k"], yi, "o", color=C1, ms=6, mec="white", mew=0.8, zorder=4)
            if e3v is not None:
                ax.plot(e3v, yi, "D", mfc="white", mec=INK, ms=4.5, mew=0.9, zorder=5)
            ax.text(1.02, yi, f"n={v['n_hogs']:,}", transform=ax.get_yaxis_transform(), fontsize=6, color=INK2, va="center")
        ax.set_yticks(y); ax.set_yticklabels([n for n, _, _ in items])
        ax.axvline(allk["k"], color=MUTED, lw=1, ls="--")
        ax.set_xlabel("pooled k (ratio of summed branch lengths)")
        ax.set_title(title, loc="left"); ax.grid(axis="y", visible=False)
        panel_label(ax, lab)

    cl = [("α", "alpha", s3.get("k_alpha")), ("β", "beta", None), ("α/β", "alphabeta", None), ("mixed", "mixed", None),
          ("no CATH", "none", None)]
    ccd = E3.get("cath_contrast") or {}
    e3c = {"alpha": ccd.get("k_alpha"), "beta": ccd.get("k_beta")}
    dots(fig.add_subplot(gs[0, 1]), [(n, cath[k], e3c.get(k)) for n, k, _ in cl if k in cath], "By CATH class", "b")
    sl = {"1": "<100", "2": "100–300", "3": "300–1000", "4": "≥1000"}
    hogroot = pl.read_parquet(ROOT / "data" / "interim" / "hogs.parquet")
    dots(fig.add_subplot(gs[0, 2]), [(f"{sl[s]} Mya", strat[s], s3.get(f"k_stratum_{s}", {}).get("estimate"))
                                     for s in ("1", "2", "3", "4") if s in strat], "By age of sampled clade", "c")
    dots(fig.add_subplot(gs[1, 1:]), [("all branches", allk, E3["estimate"]),
                                      ("internal branches only", intk, s3["internal_edges_only"]["estimate"])],
         "All vs internal branches", "d")
    fig.legend(handles=[Line2D([], [], marker="o", color=C1, ms=5, lw=2, label="all data (tuning + confirmation sets), 95% family-bootstrap CI; n = families"),
                        Line2D([], [], marker="D", ls="", mfc="white", mec=INK, ms=4.5, label="pre-registered value (confirmation set only)")],
               loc="lower center", bbox_to_anchor=(0.6, -0.06), ncol=2, fontsize=6.5)
    KEY["D_stratum_def"] = "hogs.parquet stratum 1..4 (root-age strata; see ANALYSIS_PLAN §0)"
    fig.savefig(OUT / "figD_ratio_k.png")
    plt.close(fig)
    return e


# ================================================================================================ Fig E
def figE(B, e=None):
    rng = np.random.default_rng(SEED + 1)
    e = branch_edges() if e is None else e
    bed = [0, 10, 30, 100, 200, 300, 600, 1000, 2000]
    tau = e["tau"].to_numpy()
    e = e.with_columns(pl.Series("abin", np.digitize(tau, bed) - 1))
    hogs = np.unique(e["hog"].to_numpy())
    hidx = {h: i for i, h in enumerate(hogs)}
    H = len(hogs)
    code = np.array([hidx[h] for h in e["hog"].to_numpy()])
    W = rng.multinomial(H, np.full(H, 1.0 / H), size=B).astype(float)
    term = e["is_terminal"].to_numpy()
    ab = e["abin"].to_numpy()
    baa, bt = e["b_aa"].to_numpy(), e["b_t"].to_numpy()
    rows = []
    for kind, m0 in (("all", np.ones(len(term), bool)), ("terminal", term), ("internal", ~term)):
        for k in range(len(bed) - 1):
            m = m0 & (ab == k)
            nh = len(np.unique(code[m]))
            if nh < 30:
                continue
            SA = np.bincount(code[m], baa[m], H); ST = np.bincount(code[m], bt[m], H)
            kk = SA.sum() / ST.sum()
            kb = (W @ SA) / (W @ ST)
            rows.append(dict(edges=kind, age_bin_lo=bed[k], age_bin_hi=bed[k + 1], n_edges=int(m.sum()), n_hogs=nh,
                             k=kk, lo=float(np.quantile(kb, .025)), hi=float(np.quantile(kb, .975)),
                             frac_terminal=float(term[(ab == k)].mean()) if kind == "all" else None,
                             mean_SA_per_edge=float(baa[m].mean()), mean_ST_per_edge=float(bt[m].mean())))
    r = pl.DataFrame(rows)
    r.write_csv(OUT / "figE_terminal_internal.tsv", separator="\t")
    KEY["E_bins"] = r.to_dicts()
    E4 = ep("E4")
    KEY["E_E4_beta"] = E4["estimate"]; KEY["E_E4_internal_only_beta"] = E4["sensitivity"]["internal_edges_only"]["estimate"]
    KEY["E_E4_internal_only_p"] = E4["sensitivity"]["internal_edges_only"]["p_unadjusted"]
    E5 = ep("E5")
    KEY["E_E5_decision"] = E5["decision"]; KEY["E_E5_internal_decision"] = E5["sensitivity"]["internal_edges_only"]["decision"]

    fig, axs = plt.subplots(2, 1, figsize=(5.2, 4.8), constrained_layout=True, sharex=True,
                            gridspec_kw=dict(height_ratios=[3, 1]))
    ax = axs[0]
    off = {"all": 1.0, "terminal": 0.94, "internal": 1.06}
    cc = {"all": MUTED, "terminal": C2, "internal": C1}
    lab = {"all": "all edges", "terminal": "terminal edges (tips)", "internal": "internal edges"}
    for kind in ("all", "terminal", "internal"):
        d = r.filter(pl.col("edges") == kind)
        x = np.sqrt(np.maximum(d["age_bin_lo"].to_numpy(), 3) * d["age_bin_hi"].to_numpy()) * off[kind]
        ax.plot(x, d["k"], marker="o", color=cc[kind], ms=5, mec="white", mew=0.7, lw=1.6 if kind != "all" else 1.2,
                ls="-" if kind != "all" else ":", label=lab[kind], zorder=4 if kind != "all" else 3)
        ax.vlines(x, d["lo"], d["hi"], color=cc[kind], lw=1.6, alpha=0.8)
    ax.set_xscale("log"); ax.set_xlim(3, 2000); logx(ax)
    ax.axhline(1, color=MUTED, lw=0.8, ls="--")
    ax.set_ylabel("k = Σ AA / Σ 3Di branch length")
    ax.set_title("Sequence-to-structure ratio by branch age and branch type", loc="left")
    ax.legend(loc="upper left")
    panel_label(ax, "a")
    ax = axs[1]
    d = r.filter(pl.col("edges") == "all")
    x = np.sqrt(np.maximum(d["age_bin_lo"].to_numpy(), 3) * d["age_bin_hi"].to_numpy())
    ax.bar(x, d["frac_terminal"], width=x * 0.45, color=C2, edgecolor="white")
    ax.set_ylim(0, 1); ax.set_ylabel("Share terminal")
    ax.set_xlabel("Branch midpoint age (Mya, TimeTree; log scale)")
    ax.grid(axis="x", visible=False)
    panel_label(ax, "b")
    fig.savefig(OUT / "figE_terminal_internal.png")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="A,B,C,D,E")
    ap.add_argument("--B", type=int, default=1000)
    a = ap.parse_args()
    assert confirm_allowed(), "set CONFIRM_RUN=1 (exploratory on all data after the confirmation run)"
    kp = OUT / "key_numbers.json"
    if kp.exists():
        KEY.update(json.loads(kp.read_text()))
    todo = a.only.split(",")
    e = None
    for f in todo:
        t0 = time.time()
        if f == "A": figA()
        if f == "B": figB()
        if f == "C": figC()
        if f == "D": e = figD(a.B)
        if f == "E": figE(a.B, e)
        log(f"fig {f} done {time.time() - t0:.0f}s")
        kp.write_text(json.dumps(KEY, indent=1, default=float))


if __name__ == "__main__":
    main()
