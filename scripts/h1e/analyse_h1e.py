#!/usr/bin/env python3
"""X4 (ANALYSIS_PLAN §6; DESIGN H1e, §11.10) -- EXPLORATORY, unadjusted 95% CIs.

Do AlphaFold models distort structural divergence between orthologs?  Same 1:1 ortholog pairs, residue-matched
(inputs from scripts/h1e/prep_x4.py: AF models trimmed to the residues observed in the PDB chain, identical
sequences, identical pLDDT>=70 core masks), compared PDB-PDB vs AF-AF.

  python scripts/h1e/analyse_h1e.py [--B 2000]

Distances d (larger = more divergent): d_TM = 1 - TM(short), d_LDDT = 1 - global LDDT, d_LDDTcore = 1 - core LDDT,
d_3Di = 3Di ML distance (Q.3Di.AF+G4, tuned alpha 1.60; saturated rows dropped), RMSD_core (A).
Delta = d_AF - d_PDB per ortholog pair (negative = AF pair looks MORE similar than the experimental pair).
d_AA = AA ML distance (LG+G4, alpha 1.64) of the matched AF pair on the same alignment (identical sequences, so the
AF and PDB pairs differ only through the alignment).
Inference: cluster bootstrap over top-level HOGs (all pairs sharing a protein share a HOG), B resamples, percentile
95% CI; seed SeedSequence(20261030, spawn_key=(1010,)) [stage 'x4_bootstrap', outside config/seeds.yaml -> C-item].
Chain choices: post (released > 2018-04-30, after the AF2 training cutoff; §11.10 primary), pre (<= cutoff, inside
AF2 training), any (best chain of any date). Floor: same protein, matched AF model vs its own PDB chain.
Outputs: results/h1e/x4_*.tsv, results/exploratory/x4.json, figures/h1e/*.png (300 dpi).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
from scipy import stats  # noqa: E402

PROJ = Path(__file__).resolve().parents[2]
WD = PROJ / "data/interim/h1e_x4"
RES = PROJ / "results/h1e"; FIG = PROJ / "figures/h1e"
SEED = 20261030; STAGE = 1010
MISMATCH_MAX = 0.05          # protein excluded if > 5% of observed PDB residues differ from the AF sequence
METRICS = {  # name: (expression on a pairs frame, axis label)
    "d_TM": (1 - pl.col("tm_short"), "1 − TM-score (norm. by shorter)"),
    "d_LDDT": (1 - pl.col("lddt_global"), "1 − global LDDT"),
    "d_LDDTcore": (1 - pl.col("lddt_core"), "1 − core LDDT"),
    "d_3Di": (pl.when(pl.col("tdi_sat")).then(None).otherwise(pl.col("tdi_ml")), "3Di ML distance (subst./site)"),
    "RMSD_core": (pl.col("rmsd_core"), "core RMSD (Å)"),
}
ALNS = ("tdiaa", "tmalign")
CHOICES = ("post", "pre", "any")
C_POST, C_PRE, C_GREY = "#2a78d6", "#eb6834", "#8a8984"


def acc_of(name: str) -> str:
    return name.split("_")[1]


def load(aln: str) -> pl.DataFrame:
    d = pl.concat([pl.read_parquet(f) for f in sorted((WD / "pairs" / aln).glob("shard_*.parquet"))], how="diagonal_relaxed")
    # structure pairs were run once; re-expand to every pair set they belong to (prep_x4.py pair_meta.tsv)
    meta = pl.read_csv(WD / "pair_meta.tsv", separator="\t", infer_schema_length=0).select("pair_set", "g1", "g2")
    d = d.drop("pair_set").with_columns(pl.min_horizontal("gene_a", "gene_b").alias("g1"),
                                        pl.max_horizontal("gene_a", "gene_b").alias("g2")).join(meta, on=["g1", "g2"])
    d = d.with_columns(pl.col("gene_a").map_elements(acc_of, return_dtype=pl.Utf8).alias("acc1"),
                       pl.col("gene_b").map_elements(acc_of, return_dtype=pl.Utf8).alias("acc2"))
    d = d.with_columns(pl.min_horizontal("acc1", "acc2").alias("k1"), pl.max_horizontal("acc1", "acc2").alias("k2"))
    return d.with_columns(*[e.alias(m) for m, (e, _) in METRICS.items()])


def boot_idx(hogs: np.ndarray, B: int, rng):
    u, inv = np.unique(hogs, return_inverse=True)
    members = [np.where(inv == i)[0] for i in range(len(u))]
    for _ in range(B):
        pick = rng.integers(0, len(u), len(u))
        yield np.concatenate([members[i] for i in pick])


def ols_slope(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    xm = x - x.mean()
    return float((xm * (y - y.mean())).sum() / (xm ** 2).sum())


def ci(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if v.size else [None, None]


def paired_stats(j: pl.DataFrame, m: str, B: int, rng) -> dict:
    j = j.filter(pl.col(f"{m}_af").is_not_null() & pl.col(f"{m}_pdb").is_not_null()
                 & pl.col(f"{m}_af").is_not_nan() & pl.col(f"{m}_pdb").is_not_nan())
    a = j[f"{m}_af"].to_numpy(); p = j[f"{m}_pdb"].to_numpy(); x = j["d_aa"].to_numpy(); h = j["hog"].to_numpy()
    D = a - p
    def est(ix):
        return (D[ix].mean(), np.median(D[ix]), a[ix].mean() / p[ix].mean(), ols_slope(x[ix], D[ix]),
                ols_slope(x[ix], a[ix]), ols_slope(x[ix], p[ix]))
    pt = est(np.arange(len(D)))
    bs = np.array([est(ix) for ix in boot_idx(h, B, rng)])
    names = ["mean_delta", "median_delta", "ratio_mean_af_over_pdb", "slope_delta_vs_daa", "slope_af_vs_daa", "slope_pdb_vs_daa"]
    out = {"n_pairs": int(len(D)), "n_hogs": int(len(np.unique(h))), "median_d_af": float(np.median(a)),
           "median_d_pdb": float(np.median(p)), "frac_af_smaller": float(np.mean(D < 0)),
           "wilcoxon_p_unclustered": float(stats.wilcoxon(D).pvalue) if np.any(D != 0) else None,
           "spearman_delta_daa": float(stats.spearmanr(x, D).statistic)}
    for i, nm in enumerate(names):
        out[nm] = float(pt[i]); out[nm + "_ci95"] = ci(bs[:, i])
    out["slope_ratio_af_over_pdb"] = float(pt[4] / pt[5]) if pt[5] else None
    out["slope_ratio_af_over_pdb_ci95"] = ci(bs[:, 4] / bs[:, 5])
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--B", type=int, default=2000)
    a = ap.parse_args()
    RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    (PROJ / "results/exploratory").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(STAGE,)))

    ms = pl.read_csv(WD / "match_summary.tsv", separator="\t", infer_schema_length=0).filter(pl.col("status") == "ok")
    ms = ms.with_columns((pl.col("n_mismatch").cast(pl.Float64) / pl.col("n_pdb_bb").cast(pl.Float64)).alias("mis"))
    bad_tags = set(ms.filter(pl.col("mis") > MISMATCH_MAX)["tag"].to_list())
    chains = pl.read_csv(WD / "chains_selected.tsv", separator="\t", infer_schema_length=0)
    out = {"label": "EXPLORATORY (X4, unadjusted 95% CIs)", "B": a.B, "seed": [SEED, STAGE],
           "mismatch_max": MISMATCH_MAX, "n_chain_tags_excluded_mismatch": len(bad_tags),
           "match": {"n_chain_tags": ms.height, "median_frac_af_residues_kept": float(ms["frac_keep_of_af"].cast(pl.Float64).median()),
                     "median_frac_core_kept": float(ms["frac_core_keep"].cast(pl.Float64).median())},
           "alpha": {"aa": 1.64, "tdi": 1.60}, "results": {}, "floor": {}, "post_vs_pre_same_pairs": {}, "counts": {}}

    rows, floor_rows, pp_rows = [], [], []
    J = {}
    for aln in ALNS:
        d = load(aln)
        d = d.with_columns(pl.concat_str([pl.col("gene_a"), pl.col("gene_b")]).alias("_g"))
        tag_ok = lambda df: df.filter(~pl.col("gene_a").str.slice(2).is_in(list(bad_tags)) &  # noqa: E731
                                      ~pl.col("gene_b").str.slice(2).is_in(list(bad_tags)))
        if aln == "tdiaa":
            daa = d.filter(pl.col("pair_set").str.starts_with("x4m_af_") & (pl.col("status") == "ok")).select(
                pl.col("pair_set").str.replace("x4m_af_", "").alias("choice"), "k1", "k2", pl.col("aa_ml").alias("d_aa"))
        for c in CHOICES:
            af = tag_ok(d.filter(pl.col("pair_set") == f"x4m_af_{c}"))
            pdb = tag_ok(d.filter(pl.col("pair_set") == f"x4m_pdb_{c}"))
            n_all = af.height
            cols = ["k1", "k2", "status", "aa_ml", "n_core"] + list(METRICS)
            j = af.select("k1", "k2", "hog", "split", "species_a", "species_b", *[pl.col(x).alias(x + "_af") for x in cols[2:]]).join(
                pdb.select("k1", "k2", *[pl.col(x).alias(x + "_pdb") for x in cols[2:]]), on=["k1", "k2"])
            j = j.join(daa.filter(pl.col("choice") == c).drop("choice"), on=["k1", "k2"], how="left")
            j = j.filter((pl.col("status_af") == "ok") & (pl.col("status_pdb") == "ok") & pl.col("d_aa").is_not_null())
            J[(aln, c)] = j
            out["counts"][f"{aln}/{c}"] = {"pairs_in_set": n_all, "pairs_both_ok": j.height,
                                           "by_split": dict(zip(*j["split"].value_counts().to_dict(as_series=False).values())),
                                           "aa_ml_gap_af_minus_pdb_median": float((j["aa_ml_af"] - j["aa_ml_pdb"]).median())}
            for m in METRICS:
                for sub, jj in (("all", j), ("tune_none", j.filter(pl.col("split") != "confirm"))):
                    if jj.height < 10:
                        continue
                    s = paired_stats(jj, m, a.B, rng)
                    out["results"][f"{aln}/{c}/{sub}/{m}"] = s
                    rows.append({"aln": aln, "choice": c, "subset": sub, "metric": m} | {
                        k: (json.dumps([round(v, 5) if v is not None else None for v in vv]) if isinstance(vv, list) else vv)
                        for k, vv in s.items()})
            # floor: same protein, matched AF vs own PDB chain
            fl = tag_ok(d.filter((pl.col("pair_set") == f"x4m_floor_{c}") & (pl.col("status") == "ok")))
            for m in METRICS:
                v = fl[m].drop_nulls().drop_nans().to_numpy()
                if v.size == 0:
                    continue
                bs = [np.median(v[rng.integers(0, v.size, v.size)]) for _ in range(a.B)]
                f = {"n_proteins": int(v.size), "median": float(np.median(v)), "median_ci95": ci(bs),
                     "mean": float(v.mean()), "q90": float(np.quantile(v, 0.9))}
                out["floor"][f"{aln}/{c}/{m}"] = f
                floor_rows.append({"aln": aln, "choice": c, "metric": m} | {k: (json.dumps(vv) if isinstance(vv, list) else vv) for k, vv in f.items()})
            if c == "post":
                J[(aln, "floor_post")] = fl
            if c == "pre":
                J[(aln, "floor_pre")] = fl
        # same ortholog pairs with both a post and a pre chain on both sides: Delta_post - Delta_pre
        jp, jr = J[(aln, "post")], J[(aln, "pre")]
        both = jp.join(jr, on=["k1", "k2"], suffix="_pre")
        for m in METRICS:
            bb = both.filter(pl.all_horizontal([pl.col(x).is_not_null() & pl.col(x).is_not_nan() for x in
                                                (f"{m}_af", f"{m}_pdb", f"{m}_af_pre", f"{m}_pdb_pre")]))
            if bb.height < 5:
                continue
            dd = ((bb[f"{m}_af"] - bb[f"{m}_pdb"]) - (bb[f"{m}_af_pre"] - bb[f"{m}_pdb_pre"])).to_numpy()
            h = bb["hog"].to_numpy()
            bs = [dd[ix].mean() for ix in boot_idx(h, a.B, rng)]
            r = {"n_pairs": int(bb.height), "mean_delta_post_minus_pre": float(dd.mean()), "ci95": ci(bs),
                 "mean_delta_post": float((bb[f"{m}_af"] - bb[f"{m}_pdb"]).mean()),
                 "mean_delta_pre": float((bb[f"{m}_af_pre"] - bb[f"{m}_pdb_pre"]).mean())}
            out["post_vs_pre_same_pairs"][f"{aln}/{m}"] = r
            pp_rows.append({"aln": aln, "metric": m} | {k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r.items()})

    # unmatched reference (node-1 style: full AF models vs raw PDB chains), tdiaa, post
    d = load("tdiaa")
    u_af = d.filter((pl.col("pair_set") == "x4u_af") & (pl.col("status") == "ok"))
    u_pdb = d.filter((pl.col("pair_set") == "x4u_pdb_post") & (pl.col("status") == "ok"))
    uj = u_af.select("k1", "k2", *[pl.col(m).alias(m + "_af") for m in METRICS]).join(
        u_pdb.select("k1", "k2", *[pl.col(m).alias(m + "_pdb") for m in METRICS]), on=["k1", "k2"])
    out["unmatched_reference_tdiaa_post"] = {m: {"n": uj.height, "median_delta": float((uj[m + "_af"] - uj[m + "_pdb"]).drop_nans().median())}
                                             for m in METRICS}

    pl.DataFrame(rows).write_csv(RES / "x4_paired_summary.tsv", separator="\t")
    pl.DataFrame(floor_rows).write_csv(RES / "x4_floor.tsv", separator="\t")
    pl.DataFrame(pp_rows).write_csv(RES / "x4_post_vs_pre_same_pairs.tsv", separator="\t")
    J[("tdiaa", "post")].write_csv(WD / "x4_pairs_tdiaa_post.tsv", separator="\t")  # per-pair rows stay in data/interim
    jp = J[("tdiaa", "post")].with_columns(pl.min_horizontal(pl.col("species_a").cast(pl.Int64), pl.col("species_b").cast(pl.Int64)).alias("s1"),
                                           pl.max_horizontal(pl.col("species_a").cast(pl.Int64), pl.col("species_b").cast(pl.Int64)).alias("s2"))
    out["species_pairs_post"] = {f"{r[0]}-{r[1]}": r[2] for r in jp.group_by("s1", "s2").len().sort("len", descending=True).iter_rows()}
    out["d_aa_quantiles"] = {c: [float(q) for q in np.quantile(J[("tdiaa", c)]["d_aa"].to_numpy(), [0.1, 0.25, 0.5, 0.75, 0.9, 1.0])]
                             for c in CHOICES}
    json.dump(out, open(PROJ / "results/exploratory/x4.json", "w"), indent=1)
    json.dump(out, open(RES / "x4_results.json", "w"), indent=1)
    figures(J, out)
    print(json.dumps({k: v for k, v in out.items() if k in ("counts", "match", "n_chain_tags_excluded_mismatch")}, indent=1))


# ------------------------------------------------------------------------------------------------ figures
def style():
    plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.edgecolor": "#52514e", "axes.linewidth": 0.6,
                         "xtick.color": "#52514e", "ytick.color": "#52514e", "legend.frameon": False,
                         "savefig.dpi": 300, "font.family": "DejaVu Sans"})


def scatter_fig(J, aln, fn):
    style()
    fig, axs = plt.subplots(1, 5, figsize=(13, 2.9))
    for ax, (m, (_, lab)) in zip(axs, METRICS.items()):
        for c, col, nm in (("pre", C_PRE, "pre-2018-04-30 chains (in AF2 training)"),
                           ("post", C_POST, "post-2018-04-30 chains (after AF2 cutoff)")):
            j = J[(aln, c)].filter(pl.col(f"{m}_af").is_not_null() & pl.col(f"{m}_pdb").is_not_null())
            ax.scatter(j[f"{m}_pdb"], j[f"{m}_af"], s=9, color=col, alpha=0.65, lw=0.3, ec="white",
                       label=f"{nm} (n={j.height})", zorder=3 if c == "post" else 2)
        lo, hi = ax.get_xlim(); lo2, hi2 = ax.get_ylim(); L = (min(lo, lo2, 0), max(hi, hi2))
        ax.plot(L, L, color="#52514e", lw=0.8, ls="--", zorder=1, label="identity (d_AF = d_PDB)")
        ax.set_xlim(L); ax.set_ylim(L); ax.set_aspect("equal")
        ax.set_xlabel(f"PDB–PDB: {lab}"); ax.set_ylabel(f"AF–AF: {lab}")
        ax.grid(color="#e6e5e0", lw=0.4)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.07), markerscale=1.6)
    fig.suptitle(f"Ortholog divergence, AF models vs experimental structures (residue-matched, {aln} alignment; points below the line: AF more similar)",
                 y=1.13, fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / fn, bbox_inches="tight"); plt.close(fig)


def delta_fig(J, out, aln, fn):
    style()
    fig, axs = plt.subplots(1, 5, figsize=(13, 2.9))
    for ax, (m, (_, lab)) in zip(axs, METRICS.items()):
        ax.axhline(0, color="#52514e", lw=0.8, ls="--", zorder=1)
        ys, txt = [], []
        for c, col, nm in (("pre", C_PRE, "pre-cutoff chains"), ("post", C_POST, "post-cutoff chains")):
            j = J[(aln, c)].filter(pl.col(f"{m}_af").is_not_null() & pl.col(f"{m}_pdb").is_not_null())
            x = j["d_aa"].to_numpy(); y = (j[f"{m}_af"] - j[f"{m}_pdb"]).to_numpy()
            ax.scatter(x, y, s=9, color=col, alpha=0.6, lw=0.3, ec="white", zorder=2,
                       label=f"{nm} (n={j.height})" if m == "d_TM" else None)
            ys.append(y)
            s = out["results"].get(f"{aln}/{c}/all/{m}")
            if s:
                xs = np.linspace(0, np.nanmax(x), 50)
                b = s["slope_delta_vs_daa"]; y0 = np.mean(y) - b * np.mean(x)
                lo, hi = s["slope_delta_vs_daa_ci95"]
                ax.plot(xs, y0 + b * xs, color=col, lw=2, zorder=4)
                txt.append(f"{c}: slope {b:+.3g}\n  [{lo:+.3g}, {hi:+.3g}]")
        yy = np.concatenate(ys); ql, qh = np.nanquantile(yy, [0.01, 0.99]); pad = 0.15 * (qh - ql)
        ax.set_ylim(ql - pad, max(qh, 0) + pad)
        ax.text(0.98, 0.03, "\n".join(txt), transform=ax.transAxes, ha="right", va="bottom", fontsize=6, color="#52514e",
                bbox=dict(fc="white", ec="none", alpha=0.8))
        ax.set_xlabel("AA ML distance (subst./site)"); ax.set_ylabel(f"Δ = AF − PDB\n{lab}")
        ax.grid(color="#e6e5e0", lw=0.4)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.06), markerscale=1.6)
    fig.suptitle(f"Paired AF − PDB difference vs sequence divergence ({aln}; OLS slope, 95% CI from HOG-cluster bootstrap; "
                 "y-axis clipped to the 1–99% range; exploratory)", y=1.12, fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / fn, bbox_inches="tight"); plt.close(fig)


def floor_fig(J, out, fn):
    style()
    fig, axs = plt.subplots(1, 5, figsize=(13, 2.6))
    for ax, (m, (_, lab)) in zip(axs, METRICS.items()):
        data, labs = [], []
        for c in ("post", "pre"):
            fl = J[("tdiaa", f"floor_{c}")][m].drop_nulls().drop_nans().to_numpy()
            data.append(fl); labs.append(f"floor\n{c}\nn={fl.size}")
        jd = J[("tdiaa", "post")]
        for col, nm in ((f"{m}_pdb", "orth.\nPDB–PDB"), (f"{m}_af", "orth.\nAF–AF")):
            v = jd[col].drop_nulls().drop_nans().to_numpy(); data.append(v); labs.append(f"{nm}\nn={v.size}")
        bp = ax.boxplot(data, widths=0.55, showfliers=False, patch_artist=True, medianprops=dict(color="#0b0b0b", lw=1.2))
        for p, col in zip(bp["boxes"], (C_POST, C_PRE, C_GREY, C_GREY)):
            p.set_facecolor(col); p.set_alpha(0.45); p.set_edgecolor("#52514e")
        for i, v in enumerate(data):
            ax.scatter(np.full(v.size, i + 1) + np.random.default_rng(i).uniform(-0.18, 0.18, v.size), v, s=3,
                       color="#52514e", alpha=0.35, lw=0, zorder=3)
        ax.set_xticks(range(1, 5)); ax.set_xticklabels(labs, fontsize=6); ax.set_ylabel(lab)
        ax.grid(axis="y", color="#e6e5e0", lw=0.4)
    fig.suptitle("Same-protein AF-vs-PDB floor (matched residues, tdiaa) compared with post-cutoff ortholog divergences", y=1.03, fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / fn, bbox_inches="tight"); plt.close(fig)


def figures(J, out):
    scatter_fig(J, "tdiaa", "x4_paired_scatter_tdiaa.png")
    scatter_fig(J, "tmalign", "x4_paired_scatter_tmalign.png")
    delta_fig(J, out, "tdiaa", "x4_delta_vs_daa_tdiaa.png")
    delta_fig(J, out, "tmalign", "x4_delta_vs_daa_tmalign.png")
    floor_fig(J, out, "x4_floor_tdiaa.png")


if __name__ == "__main__":
    main()
