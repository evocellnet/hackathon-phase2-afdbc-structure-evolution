#!/usr/bin/env python3
"""Goal 2-ii item 2: why do HOGs split across AFDB v3 structural clusters?
Pair-level logistic regression on within-HOG pairs (<= 8 per HOG, both members eligible) from HOGs whose eligible
members span > 1 cluster:  P(different cluster) ~ TED CATH-architecture string differs + |#TED domains diff|
+ |log length ratio| + mean pLDDT + AA distance (aa_ml_seq, seqonly LG+G alpha=1; capped at 5, saturation flag)
+ relation (ortholog / out-paralog / in-paralog) [+ 'no CATH label on either' indicator].
SEs clustered by HOG (sandwich). Dominance: standardized coefficients, drop-one deviance, single-predictor AUC.
Unadjusted 95% CIs (exploratory).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import statsmodels.api as sm

sys.path.insert(0, str(Path(__file__).resolve().parent))
from g2_common import OUT_F, OUT_I, OUT_R, ROOT  # noqa: E402

C1, C2, INK, INK2 = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"
p = pl.read_parquet(OUT_I / "pairs_sample_dist.parquet")
ann = pl.read_parquet(ROOT / "data/processed/annotations/protein_annot.parquet",
                      columns=["acc", "ted_arch_cath", "ted_n_domains", "ted_n_cath", "len_v4", "mean_plddt"])


def side(s):
    return ann.rename({c: f"{c}_{s}" for c in ann.columns if c != "acc"}).rename({"acc": f"gene_{s}"})


p = p.join(side("a"), on="gene_a", how="left").join(side("b"), on="gene_b", how="left")
p = p.with_columns(
    clu_diff=(pl.col("clu_a") != pl.col("clu_b")),
    arch_diff=pl.col("ted_arch_cath_a").fill_null("") != pl.col("ted_arch_cath_b").fill_null(""),
    # topology-level (C.A.T) architecture: removes H-vs-T label-depth differences of the same domain
    archT_diff=pl.col("ted_arch_cath_a").fill_null("").str.split("|").list.eval(
        pl.element().str.split(".").list.head(3).list.join(".")).list.join("|")
        != pl.col("ted_arch_cath_b").fill_null("").str.split("|").list.eval(
        pl.element().str.split(".").list.head(3).list.join(".")).list.join("|"),
    both_labelled=(pl.col("ted_n_cath_a") > 0) & (pl.col("ted_n_cath_b") > 0),
    no_label_either=(pl.col("ted_n_cath_a") == 0) | (pl.col("ted_n_cath_b") == 0),
    ndom_diff=(pl.col("ted_n_domains_a").cast(pl.Int64) - pl.col("ted_n_domains_b").cast(pl.Int64)).abs(),
    abs_log_len=(pl.col("len_v4_a").log() - pl.col("len_v4_b").log()).abs(),
    plddt=(pl.col("mean_plddt_a") + pl.col("mean_plddt_b")) / 2,
    d_aa=pl.when(pl.col("aa_sat") | pl.col("aa_ml_seq").is_nan() | pl.col("aa_ml_seq").is_null()).then(5.0)
          .otherwise(pl.col("aa_ml_seq").clip(0, 5)),
    sat=(pl.col("aa_sat") | pl.col("aa_ml_seq").is_nan()).fill_null(True) | (pl.col("aa_ml_seq") >= 5).fill_null(True),
)
p.write_parquet(OUT_I / "pairs_sample_feat.parquet")
hs = pl.read_parquet(OUT_I / "hog_cluster_stats.parquet").select("hog", "split", "n_elig")
d = p.join(hs, on="hog").filter(pl.col("split") & pl.col("relation").is_in(["ortholog", "out_paralog", "in_paralog"]))
d = d.with_columns(out_para=(pl.col("relation") == "out_paralog").cast(pl.Float64),
                   in_para=(pl.col("relation") == "in_paralog").cast(pl.Float64),
                   y=pl.col("clu_diff").cast(pl.Float64))
X_COLS = ["arch_diff", "no_label_either", "ndom_diff", "abs_log_len", "plddt", "d_aa", "sat", "out_para", "in_para"]
d = d.with_columns([pl.col(c).cast(pl.Float64) for c in X_COLS]).drop_nulls(X_COLS + ["y"])
out = {"n_pairs": d.height, "n_hogs": d["hog"].n_unique(), "frac_diff_cluster": float(d["y"].mean()),
       "relation_counts": dict(zip(*d["relation"].value_counts().to_dict(as_series=False).values())),
       "frac_arch_diff": float(d["arch_diff"].mean()), "frac_sat": float(d["sat"].mean()),
       "median_d_aa": float(d["d_aa"].median())}
groups = d["hog"].cast(pl.Categorical).to_physical().to_numpy()
y = d["y"].to_numpy()


def fit(cols, data=d, grp=groups, yy=y):
    X = sm.add_constant(data.select(cols).to_numpy())
    return sm.Logit(yy, X).fit(disp=0, maxiter=200, cov_type="cluster", cov_kwds={"groups": grp})


def auc(score, yy):
    from scipy.stats import rankdata
    r = rankdata(score)
    n1 = yy.sum(); n0 = len(yy) - n1
    return float((r[yy == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


full = fit(X_COLS)
ci = full.conf_int()
sd = {c: float(d[c].std()) for c in X_COLS}
coef = []
for i, c in enumerate(["const"] + X_COLS):
    row = dict(term=c, coef=float(full.params[i]), se_hogcl=float(full.bse[i]), lo95=float(ci[i][0]),
               hi95=float(ci[i][1]), OR=float(np.exp(full.params[i])), p=float(full.pvalues[i]))
    if c != "const":
        binary = set(np.unique(d[c].to_numpy())) <= {0.0, 1.0}
        row["per"] = "0->1" if binary else "1 SD"
        row["sd"] = sd[c]
        row["std_coef"] = row["coef"] * (1 if binary else sd[c])
        # drop-one deviance (unweighted model deviance difference; descriptive)
        red = fit([x for x in X_COLS if x != c])
        row["drop1_dDev"] = float((-2 * red.llf) - (-2 * full.llf))
        row["auc_alone"] = auc(d[c].to_numpy(), y)
    coef.append(row)
ctab = pl.DataFrame(coef)
ctab.write_csv(OUT_R / "split_model_coefs.tsv", separator="\t")
null_dev = -2 * sm.Logit(y, np.ones((len(y), 1))).fit(disp=0).llf
out["mcfadden_R2_full"] = float(1 - (-2 * full.llf) / null_dev)
out["auc_full"] = auc(full.predict(), y)
# group shares of explained deviance (architecture block vs length vs pLDDT vs AA vs relation)
blocks = {"architecture": ["arch_diff", "no_label_either", "ndom_diff"], "length": ["abs_log_len"],
          "pLDDT": ["plddt"], "AA divergence": ["d_aa", "sat"], "relation": ["out_para", "in_para"]}
bl = {}
for b, cols in blocks.items():
    red = fit([x for x in X_COLS if x not in cols])
    alone = sm.Logit(y, sm.add_constant(d.select(cols).to_numpy())).fit(disp=0)
    bl[b] = dict(drop_block_dDev=float((-2 * red.llf) - (-2 * full.llf)),
                 R2_alone=float(1 - (-2 * alone.llf) / null_dev),
                 auc_alone=auc(alone.predict(), y))
out["blocks"] = bl
# sensitivity: orthologs only; and both-labelled only
for nm, flt in [("orthologs_only", pl.col("relation") == "ortholog"), ("both_labelled", pl.col("both_labelled"))]:
    dd = d.filter(flt)
    cols = [c for c in X_COLS if dd[c].n_unique() > 1]
    m = fit(cols, dd, dd["hog"].cast(pl.Categorical).to_physical().to_numpy(), dd["y"].to_numpy())
    cc = m.conf_int()
    out[f"sens_{nm}"] = {"n": dd.height, **{c: dict(coef=float(m.params[i + 1]), lo95=float(cc[i + 1][0]),
                                                      hi95=float(cc[i + 1][1])) for i, c in enumerate(cols)}}
# marginal descriptive: fraction different cluster by arch_diff x relation; by length-ratio bins
desc = (d.with_columns(lenbin=pl.col("abs_log_len").cut([0.05, 0.1, 0.2, 0.4], labels=["<0.05", "0.05-0.1", "0.1-0.2", "0.2-0.4", ">0.4"]).cast(pl.Utf8))
        .group_by("arch_diff", "lenbin").agg(n=pl.len(), frac_diff_cluster=pl.col("y").mean()).sort("arch_diff", "lenbin"))
desc.write_csv(OUT_R / "split_desc_arch_len.tsv", separator="\t")
json.dump(out, open(OUT_R / "split_model_summary.json", "w"), indent=1)
print(ctab); print(json.dumps(out, indent=1)); print(desc)

# figure: standardized log-odds with 95% CI + fraction by length bin x arch
fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6), dpi=300, gridspec_kw={"width_ratios": [1.2, 1]})
c2 = ctab.filter(pl.col("term") != "const")
lab = {"arch_diff": "CATH architecture differs", "no_label_either": "no CATH label (either)",
       "ndom_diff": "|Δ # TED domains| (per SD)", "abs_log_len": "|log length ratio| (per SD)",
       "plddt": "mean pLDDT (per SD)", "d_aa": "AA distance (per SD)", "sat": "AA distance saturated",
       "out_para": "out-paralog vs ortholog", "in_para": "in-paralog vs ortholog"}
yy = np.arange(c2.height)[::-1]
scale = np.array([1 if r["per"] == "0->1" else r["sd"] for r in c2.iter_rows(named=True)])
ax[0].errorbar(c2["std_coef"], yy, xerr=[c2["std_coef"] - c2["lo95"] * scale, c2["hi95"] * scale - c2["std_coef"]],
               fmt="o", color=C1, ms=4, capsize=2, lw=1.5)
ax[0].axvline(0, color=INK2, lw=0.8)
ax[0].set_yticks(yy, [lab[t] for t in c2["term"]])
ax[0].set_xlabel("log-odds of different cluster (standardized)")
ax[0].set_title(f"a  Within split HOGs: {out['n_pairs']:,} pairs, {out['n_hogs']:,} HOGs", loc="left", fontsize=9)
for k, col, nm in [(False, C1, "same CATH architecture"), (True, C2, "different CATH architecture")]:
    q = desc.filter(pl.col("arch_diff") == float(k))
    order = ["<0.05", "0.05-0.1", "0.1-0.2", "0.2-0.4", ">0.4"]
    q = q.with_columns(o=pl.col("lenbin").replace_strict({v: i for i, v in enumerate(order)})).sort("o")
    ax[1].plot(q["o"], q["frac_diff_cluster"], "o-", color=col, ms=4, lw=1.5, label=nm)
ax[1].set_xticks(range(5), ["<0.05", "0.05-0.1", "0.1-0.2", "0.2-0.4", ">0.4"])
ax[1].set_xlabel("|log length ratio|"); ax[1].set_ylabel("fraction in different AFDB cluster"); ax[1].set_ylim(0, 1)
ax[1].legend(frameon=False, fontsize=7)
ax[1].set_title("b  Pairs in different clusters", loc="left", fontsize=9)
fig.text(0.01, -0.03, "Logistic regression, HOG-clustered sandwich 95% CIs (unadjusted). Continuous terms per SD.",
         fontsize=6.5, color=INK2)
fig.tight_layout(); fig.savefig(OUT_F / "hog_split_model.png", bbox_inches="tight")
