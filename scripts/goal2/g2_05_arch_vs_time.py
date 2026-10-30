#!/usr/bin/env python3
"""Goal 2-ii item 4: domain-architecture (TED CATH string) change and AFDB-cluster split in TUNE-split 1:1 ortholog
pairs vs species-pair age and vs AA distance; orthologs vs paralogs at matched AA distance (tune E6 pairs and the
all-HOG pair sample). Unadjusted 95% CIs: HOG-cluster bootstrap (B=500, seed 20261030) for binned fractions;
HOG-clustered sandwich SEs for logistic models. Results tables contain no ages, only age-bin labels.
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
from g2_common import OUT_F, OUT_I, OUT_R, SEED  # noqa: E402

C1, C2, C3, INK, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2})
B = 500
rng = np.random.default_rng(SEED)
out = {}


def boot_frac(df, ycol, bcol, labels):
    """HOG-cluster bootstrap of the per-bin mean of ycol."""
    rows = []
    for lab in labels:
        d = df.filter(pl.col(bcol) == lab)
        if d.height == 0:
            continue
        g = d.group_by("hog").agg(s=pl.col(ycol).sum(), n=pl.len())
        s, n = g["s"].to_numpy().astype(float), g["n"].to_numpy().astype(float)
        est = s.sum() / n.sum()
        idx = rng.integers(0, len(s), size=(B, len(s)))
        bs = s[idx].sum(1) / n[idx].sum(1)
        rows.append(dict(bin=lab, n_pairs=int(n.sum()), n_hogs=len(s), frac=est,
                         lo95=float(np.quantile(bs, .025)), hi95=float(np.quantile(bs, .975))))
    return pl.DataFrame(rows)


t = pl.read_parquet(OUT_I / "tune_pairs_arch.parquet")
o = t.filter((pl.col("pair_set") == "primary") & (pl.col("relation") == "ortholog") & (pl.col("one_to_one") == 1))
out["tune_primary_1to1_pairs"] = o.height
out["tune_primary_1to1_hogs"] = o["hog"].n_unique()
AGEB = [0, 100, 300, 500, 700, 1000, 2000]
AGEL = ["<100", "100-300", "300-500", "500-700", "700-1000", ">1000"]
DB = [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 100]
DL = ["0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1", "1-1.5", "1.5-2", ">2"]


def binc(col, br, lab):
    e = pl.lit(None, dtype=pl.Utf8)
    for lo, hi, l in zip(br[:-1], br[1:], lab):
        e = pl.when((pl.col(col) >= lo) & (pl.col(col) < hi)).then(pl.lit(l)).otherwise(e)
    return e


o = o.with_columns(age_bin=binc("age_mya", AGEB, AGEL), d_bin=binc("aa_ml_seq", DB, DL),
                   arch_diff_i=pl.col("arch_diff").cast(pl.Int8), clu_diff_i=pl.col("clu_diff").cast(pl.Int8),
                   archT_diff_i=pl.col("archT_diff").cast(pl.Int8), ndom_chg_i=(pl.col("ndom_diff") > 0).cast(pl.Int8))
ol = o.filter(pl.col("both_labelled"))
oc = o.filter(pl.col("both_in_cluster"))
out["frac_both_labelled"] = float(o["both_labelled"].mean())
tabs = {}
tabs["arch_by_age"] = boot_frac(ol, "arch_diff_i", "age_bin", AGEL)
tabs["arch_by_age_alllabelstates"] = boot_frac(o, "arch_diff_i", "age_bin", AGEL)
tabs["arch_by_daa"] = boot_frac(ol, "arch_diff_i", "d_bin", DL)
tabs["archT_by_age"] = boot_frac(ol, "archT_diff_i", "age_bin", AGEL)
tabs["archT_by_daa"] = boot_frac(ol, "archT_diff_i", "d_bin", DL)
tabs["ndomchg_by_age"] = boot_frac(o, "ndom_chg_i", "age_bin", AGEL)
tabs["ndomchg_by_daa"] = boot_frac(o, "ndom_chg_i", "d_bin", DL)
tabs["clusplit_by_age"] = boot_frac(oc, "clu_diff_i", "age_bin", AGEL)
tabs["clusplit_by_daa"] = boot_frac(oc, "clu_diff_i", "d_bin", DL)
for k, v in tabs.items():
    v.write_csv(OUT_R / f"orth1to1_tune_{k}.tsv", separator="\t")
    out[k] = v.to_dicts()


def logit_cluster(df, y, xcols):
    d = df.drop_nulls([y] + xcols + ["hog"])
    X = sm.add_constant(d.select(xcols).to_numpy().astype(float))
    g = d["hog"].cast(pl.Categorical).to_physical().to_numpy()
    m = sm.Logit(d[y].to_numpy().astype(float), X).fit(
        disp=0, maxiter=200, cov_type="cluster", cov_kwds={"groups": g})
    ci = m.conf_int()
    return {nm: dict(coef=float(m.params[i]), se=float(m.bse[i]), lo95=float(ci[i][0]), hi95=float(ci[i][1]),
                     p=float(m.pvalues[i])) for i, nm in enumerate(["const"] + xcols)}, d.height


# age vs AA distance as predictors of architecture change (both labelled, finite non-saturated d)
ol2 = ol.filter(pl.col("age_mya").is_not_nan() & (pl.col("aa_ml_seq") < 5)).with_columns(
    age_100=pl.col("age_mya") / 100, d_aa=pl.col("aa_ml_seq"))
r1, n1 = logit_cluster(ol2, "arch_diff_i", ["age_100"])
r2, n2 = logit_cluster(ol2, "arch_diff_i", ["d_aa"])
r3, n3 = logit_cluster(ol2, "arch_diff_i", ["age_100", "d_aa"])
out["logit_arch_orth_age_only"] = r1; out["logit_arch_orth_daa_only"] = r2; out["logit_arch_orth_age_daa"] = r3
out["logit_arch_orth_n"] = n3
out["logit_archT_orth_age_daa"], _ = logit_cluster(ol2, "archT_diff_i", ["age_100", "d_aa"])

# ---- orthologs vs paralogs at matched AA distance
# (a) tune E6: e6_para (in-paralogs, recent duplications) vs e6_orth (orthologs of the copies)
e6 = t.filter(pl.col("pair_set").is_in(["e6_para", "e6_orth"]) & pl.col("both_labelled") & (pl.col("aa_ml_seq") < 5))
e6 = e6.with_columns(para=(pl.col("pair_set") == "e6_para").cast(pl.Int8), d_aa=pl.col("aa_ml_seq"),
                     d_aa2=pl.col("aa_ml_seq") ** 2, arch_diff_i=pl.col("arch_diff").cast(pl.Int8),
                     archT_diff_i=pl.col("archT_diff").cast(pl.Int8), d_bin=binc("aa_ml_seq", DB, DL))
out["logit_archT_e6_para_vs_orth"], _ = logit_cluster(e6, "archT_diff_i", ["para", "d_aa", "d_aa2"])
# near-identical noise floor (AA distance < 0.02) for 1:1 orthologs
nf = o.filter(pl.col("both_labelled") & (pl.col("aa_ml_seq") < 0.02))
out["noise_floor_dlt0.02"] = dict(n=nf.height, n_hogs=nf["hog"].n_unique(), arch=float(nf["arch_diff"].mean()),
                                  archT=float(nf["archT_diff"].mean()), ndom=float((nf["ndom_diff"] > 0).mean()),
                                  clu=float(nf["clu_diff"].mean()))
re6, ne6 = logit_cluster(e6, "arch_diff_i", ["para", "d_aa", "d_aa2"])
out["logit_arch_e6_para_vs_orth"] = re6; out["logit_arch_e6_n"] = ne6
e6tab = pl.concat([boot_frac(e6.filter(pl.col("para") == k), "arch_diff_i", "d_bin", DL).with_columns(
    relation=pl.lit(nm)) for k, nm in [(0, "ortholog (E6)"), (1, "in-paralog (E6)")]])
e6tab.write_csv(OUT_R / "e6_tune_arch_by_daa.tsv", separator="\t")

# (b) all-HOG pair sample (not split-restricted): ortholog vs paralog
ps = pl.read_parquet(OUT_I / "pairs_sample_feat.parquet")
ps = ps.filter(pl.col("both_labelled") & (pl.col("aa_ml_seq") < 5) & pl.col("relation").is_in(
    ["ortholog", "in_paralog", "out_paralog"])).with_columns(
    para=(pl.col("relation") != "ortholog").cast(pl.Int8), d_aa=pl.col("aa_ml_seq"), d_aa2=pl.col("aa_ml_seq") ** 2,
    arch_diff_i=pl.col("arch_diff").cast(pl.Int8), archT_diff_i=pl.col("archT_diff").cast(pl.Int8),
    d_bin=binc("aa_ml_seq", DB, DL))
out["logit_archT_sample_para_vs_orth"], _ = logit_cluster(ps, "archT_diff_i", ["para", "d_aa", "d_aa2"])
rps, nps = logit_cluster(ps, "arch_diff_i", ["para", "d_aa", "d_aa2"])
out["logit_arch_sample_para_vs_orth"] = rps; out["logit_arch_sample_n"] = nps
pstab = pl.concat([boot_frac(ps.filter(pl.col("relation") == r), "archT_diff_i", "d_bin", DL).with_columns(
    relation=pl.lit(r)) for r in ["ortholog", "out_paralog", "in_paralog"]])
pstab.write_csv(OUT_R / "sample_arch_by_daa_relation.tsv", separator="\t")
json.dump(out, open(OUT_R / "arch_vs_time_summary.json", "w"), indent=1)

# ---- figure
fig, ax = plt.subplots(1, 3, figsize=(11, 3.4), dpi=300)
for k, col, lab in [("arch_by_age", C1, "CATH architecture string differs"),
                    ("archT_by_age", C3, "topology-level (C.A.T) architecture differs"),
                    ("clusplit_by_age", C2, "different AFDB v3 cluster")]:
    d = tabs[k]
    x = np.arange(d.height)
    ax[0].errorbar(x, d["frac"], yerr=[d["frac"] - d["lo95"], d["hi95"] - d["frac"]], color=col, marker="o",
                   ms=4, lw=1.5, capsize=2, label=lab)
ax[0].set_xticks(np.arange(len(AGEL)), AGEL, rotation=30)
ax[0].set_xlabel("species-pair divergence time (Mya, TimeTree)")
ax[0].set_ylabel("fraction of 1:1 ortholog pairs")
ax[0].set_title("a  1:1 orthologs (tune split) vs time", loc="left", fontsize=9)
ax[0].legend(frameon=False, fontsize=7); ax[0].set_ylim(0, 1)
for k, col, lab in [("arch_by_daa", C1, "CATH string differs"), ("archT_by_daa", C3, "C.A.T architecture differs"),
                    ("clusplit_by_daa", C2, "different cluster")]:
    d = tabs[k]
    x = np.arange(d.height)
    ax[1].errorbar(x, d["frac"], yerr=[d["frac"] - d["lo95"], d["hi95"] - d["frac"]], color=col, marker="o",
                   ms=4, lw=1.5, capsize=2, label=lab)
ax[1].set_xticks(np.arange(len(DL) - 1), DL[:-1], rotation=30)  # primary pairs are capped at d_AA <= 2
ax[1].set_xlabel("AA distance (LG+Γ ML, subst./site)")
ax[1].set_title("b  1:1 orthologs vs AA distance", loc="left", fontsize=9)
ax[1].legend(frameon=False, fontsize=7); ax[1].set_ylim(0, 1)
for r, col in [("ortholog", C1), ("out_paralog", C2), ("in_paralog", C3)]:
    d = pstab.filter(pl.col("relation") == r)
    x = np.array([DL.index(b) for b in d["bin"]])
    ax[2].errorbar(x, d["frac"], yerr=[d["frac"] - d["lo95"], d["hi95"] - d["frac"]], color=col, marker="o",
                   ms=4, lw=1.5, capsize=2, label=r.replace("_", "-"))
ax[2].set_xticks(np.arange(len(DL)), DL, rotation=30)
ax[2].set_xlabel("AA distance (LG+Γ ML, subst./site)")
ax[2].set_title("c  C.A.T architecture differs: orth. vs paralogs", loc="left", fontsize=9)
ax[2].legend(frameon=False, fontsize=7); ax[2].set_ylim(0, 1)
fig.text(0.01, -0.02, "Pairs with ≥1 CATH-labelled TED domain on both proteins; error bars: 95% HOG-cluster bootstrap CI "
         "(unadjusted). c: within-HOG pair sample, all HOGs ≤600 genes.", fontsize=6.5, color=INK2)
fig.tight_layout()
fig.savefig(OUT_F / "arch_change_vs_time_and_daa.png", bbox_inches="tight")
print(json.dumps({k: v for k, v in out.items() if k.startswith("logit") or k.startswith("frac") or "n" in k[:5]}, indent=1))
for k, v in tabs.items():
    print(k); print(v)
print(pstab)
