#!/usr/bin/env python3
"""floors.json (ANALYSIS_PLAN §1.4) + figures/floors/floor_distributions.png from data/interim/pairs_floors/{aln}."""
import glob, json, os
import numpy as np, polars as pl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
P = os.path.join(ROOT, "data/interim/pairs_floors")
ALNS = ["tdiaa", "tmalign", "seqonly"]
BINS = [(0, 100), (100, 200), (200, 300), (300, 500), (500, 800), (800, 10**6)]
NAMES = ["<100", "100-200", "200-300", "300-500", "500-800", ">=800"]


def load(aln):
    df = pl.concat([pl.read_parquet(f) for f in sorted(glob.glob(f"{P}/{aln}/shard_*.parquet"))])
    return df


def binof(L):
    return np.array([next(i for i, (lo, hi) in enumerate(BINS) if lo <= v < hi) for v in L])


def qs(x):
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return dict(n=0, median=None, q05=None, q95=None, q99=None)
    return dict(n=int(len(x)), median=float(np.median(x)), q05=float(np.quantile(x, .05)),
                q95=float(np.quantile(x, .95)), q99=float(np.quantile(x, .99)))


def main():
    out = {"bins": NAMES, "bin_edges_len_a": [[lo, hi if hi < 10**6 else None] for lo, hi in BINS],
           "notes": {
               "S0": "floor_random pairs, all pairs with a Foldseek/seqonly row (status != fail_align); binned by len_a (shorter protein). "
                     "S0_ok restricts to status == ok (n_core>=50, cov_short>=0.5), the population on which lddt_D is evaluated.",
               "eta": "floor_ident pairs with status == ok; nu = eta/2; eta_aa = 0 by construction (reported measured value).",
               "TM": "tm_short = qtmscore normalised by len_a; undefined for seqonly (fixed sequence alignment).",
               "not_computed_here": "rho_m(bin), kappa_m(bin) and the lddt_D SE model are produced by the §2.3 SE calibration."}}
    data = {}
    for aln in ALNS:
        r = load(aln)
        rr = r.filter(pl.col("pair_set") == "floor_random")
        ii = r.filter(pl.col("pair_set") == "floor_ident")
        data[aln] = (rr, ii)
        res = {"n_random": rr.height, "n_ident": ii.height,
               "random_status": dict(rr["status"].value_counts().iter_rows()),
               "ident_status": dict(ii["status"].value_counts().iter_rows()),
               "S0": {}, "S0_ok": {}, "eta": {}}
        for tag, sub in (("S0", rr.filter(pl.col("status") != "fail_align")), ("S0_ok", rr.filter(pl.col("status") == "ok"))):
            b = binof(sub["len_a"].to_numpy())
            for k, nm in enumerate(NAMES):
                s = sub.filter(pl.Series(b == k))
                res[tag][nm] = {"LDDT_global": qs(s["lddt_global"].to_numpy()), "TM": qs(s["tm_short"].to_numpy())}
        iok = ii.filter(pl.col("status") == "ok")
        metrics = {"tdi_ml": iok["tdi_ml"].to_numpy(), "aa_ml": iok["aa_ml"].to_numpy(),
                   "lddt_raw": 1 - iok["lddt_global"].to_numpy(),
                   "tm_short": iok["tm_short"].to_numpy(), "lddt_global": iok["lddt_global"].to_numpy(),
                   "aa_ml_seq": iok["aa_ml_seq"].to_numpy()}
        for m, v in metrics.items():
            d = qs(v)
            res["eta"][m] = {**d, "mean": float(np.nanmean(v)) if len(v) else None}
        b = binof(iok["len_a"].to_numpy())
        res["eta_by_bin"] = {nm: {m: qs(metrics[m][b == k]) for m in ("tdi_ml", "lddt_raw", "tm_short")} for k, nm in enumerate(NAMES)}
        res["eta_tdi_ml"] = res["eta"]["tdi_ml"]["median"]; res["nu_tdi_ml"] = res["eta_tdi_ml"] / 2 if res["eta_tdi_ml"] is not None else None
        res["eta_lddt_raw"] = res["eta"]["lddt_raw"]["median"]; res["nu_lddt_raw"] = res["eta_lddt_raw"] / 2 if res["eta_lddt_raw"] is not None else None
        res["eta_aa_ml"] = res["eta"]["aa_ml"]["median"]
        res["frac_ident_exact"] = {"tm_ge_0.999": float(np.mean(metrics["tm_short"] >= 0.999)) if len(iok) else None,
                                   "lddt_ge_0.999": float(np.mean(metrics["lddt_global"] >= 0.999)) if len(iok) else None,
                                   "tdi_ml_eq_0": float(np.mean(metrics["tdi_ml"] == 0)) if len(iok) else None}
        out[aln] = res
    os.makedirs(os.path.join(ROOT, "results/floors"), exist_ok=True)
    json.dump(out, open(os.path.join(ROOT, "results/floors/floors.json"), "w"), indent=1)

    # ---------------- figure ----------------
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": .6})
    C = {"tdiaa": "#1f6fb4", "tmalign": "#d9822b"}
    fig, ax = plt.subplots(2, 3, figsize=(12, 7))
    for col, (metric, lab) in enumerate((("lddt_global", "LDDT_global"), ("tm_short", "TM (norm. by shorter)"))):
        a = ax[0, col]
        for j, aln in enumerate(("tdiaa", "tmalign")):
            rr = data[aln][0].filter(pl.col("status") != "fail_align")
            b = binof(rr["len_a"].to_numpy()); v = rr[metric].to_numpy()
            groups = [v[b == k][np.isfinite(v[b == k])] for k in range(6)]
            pos = np.arange(6) + (j - .5) * .36
            bp = a.boxplot(groups, positions=pos, widths=.3, whis=(5, 95), showfliers=False, patch_artist=True,
                           medianprops=dict(color="k", lw=1), boxprops=dict(lw=.6), whiskerprops=dict(lw=.6), capprops=dict(lw=.6))
            for p in bp["boxes"]:
                p.set_facecolor(C[aln]); p.set_alpha(.55)
            a.plot([], [], color=C[aln], lw=6, alpha=.55, label=aln)
            q = [np.quantile(g, .95) for g in groups]
            a.plot(pos, q, "v", color=C[aln], ms=4)
        a.set_xticks(range(6)); a.set_xticklabels(NAMES); a.set_xlabel("length of shorter protein (aa)")
        a.set_ylabel(lab); a.set_title(f"floor_random: {lab}\n(box 25-75, whiskers 5-95, triangle = q95)")
        a.legend(frameon=False)
    a = ax[0, 2]
    for aln in ("tdiaa", "tmalign", "seqonly"):
        rr = data[aln][0]; rr = rr.filter(pl.col("status") != "fail_align")
        a.hist(rr["lddt_global"].to_numpy(), bins=60, range=(0, 1), histtype="step", lw=1.2, density=True,
               color={"tdiaa": C["tdiaa"], "tmalign": C["tmalign"], "seqonly": "#555"}[aln], label=aln)
    a.set_xlabel("LDDT_global, all floor_random pairs"); a.set_ylabel("density"); a.legend(frameon=False)
    a.set_title("floor_random LDDT_global by alignment")
    for col, (metric, lab, rng) in enumerate((("tm_short", "TM (floor_ident)", (0, 1)), ("lddt_raw", "1 - LDDT_global (floor_ident)", (0, .5)),
                                               ("tdi_ml", "tdi_ml distance (floor_ident)", (0, .5)))):
        a = ax[1, col]
        for aln in ALNS:
            ii = data[aln][1].filter(pl.col("status") == "ok")
            v = ii["lddt_global"].to_numpy() if metric == "lddt_raw" else ii[metric].to_numpy()
            if metric == "lddt_raw":
                v = 1 - v
            if np.isfinite(v).sum() == 0:
                continue
            col_ = {"tdiaa": C["tdiaa"], "tmalign": C["tmalign"], "seqonly": "#555"}[aln]
            a.hist(np.clip(v[np.isfinite(v)], rng[0], rng[1]), bins=60, range=rng, histtype="step", lw=1.2, color=col_,
                   label=f"{aln} (median {np.nanmedian(v):.4f})")
        a.set_yscale("log"); a.set_xlabel(lab); a.set_ylabel("pairs"); a.legend(frameon=False, fontsize=7)
        a.set_title(lab + " (clipped at range end)")
    fig.tight_layout()
    os.makedirs(os.path.join(ROOT, "figures/floors"), exist_ok=True)
    fig.savefig(os.path.join(ROOT, "figures/floors/floor_distributions.png"), dpi=150)
    print(json.dumps({k: {x: out[k][x] for x in ("eta_tdi_ml", "nu_tdi_ml", "eta_lddt_raw", "eta_aa_ml", "n_random", "n_ident", "random_status", "ident_status")} for k in ALNS}, indent=1))


if __name__ == "__main__":
    main()
