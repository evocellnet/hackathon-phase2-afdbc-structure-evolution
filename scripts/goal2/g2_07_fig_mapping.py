#!/usr/bin/env python3
"""Goal 2-ii figures for items 1 and 3 (HOG <-> cluster mapping; cluster LCA vs HOG root)."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from g2_common import OUT_F, OUT_I, OUT_R  # noqa: E402

C1, C2, C3, C4, INK, INK2 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2})
ph = pl.read_parquet(OUT_I / "hog_cluster_stats.parquet")
cl = pl.read_parquet(OUT_I / "cluster_hog_stats.parquet")
sz = pl.read_csv(OUT_R / "split_by_hog_size.tsv", separator="\t")

fig, ax = plt.subplots(1, 3, figsize=(11, 3.3), dpi=300)
bins = np.unique(np.round(np.logspace(0, np.log10(ph["n_clu"].max() + 1), 30)).astype(int))
ax[0].hist(ph["n_clu"], bins=bins, color=C1, edgecolor="white", lw=0.5)
ax[0].set_xscale("log"); ax[0].set_yscale("log")
ax[0].set_xlabel("# AFDB v3 clusters per HOG"); ax[0].set_ylabel("# HOGs")
ax[0].set_title(f"a  {ph.height:,} HOGs (≥2 eligible members)", loc="left", fontsize=9)
x = np.arange(sz.height)
ax[1].plot(x, sz["frac_split"], "o-", color=C1, ms=4, lw=1.5, label="HOG spans >1 cluster")
ax[1].plot(x, sz["frac_largest"], "s-", color=C2, ms=4, lw=1.5, label="mean fraction in largest cluster")
ax[1].set_xticks(x, sz["bin"], rotation=30); ax[1].set_ylim(0, 1.02)
ax[1].set_xlabel("HOG size (eligible members)"); ax[1].legend(frameon=False, fontsize=7)
ax[1].set_title("b  Splitting grows with HOG size", loc="left", fontsize=9)
cb = np.arange(1, min(cl["n_hog"].max(), 30) + 2) - 0.5
ax[2].hist(cl["n_hog"].clip(upper_bound=30), bins=cb, color=C3, edgecolor="white", lw=0.5)
ax[2].set_yscale("log"); ax[2].set_xlabel("# HOGs per cluster (panel members; 30 = ≥30)")
ax[2].set_ylabel("# clusters")
ax[2].set_title(f"c  {cl.height:,} clusters", loc="left", fontsize=9)
fig.tight_layout(); fig.savefig(OUT_F / "hog_cluster_mapping.png", bbox_inches="tight")

# LCA consistency: categories for clusters and HOGs; heatmap of HOG root age vs cluster LCA age (per cluster)
order = ["consistent", "consistent_censored_pre_euk", "cluster_younger", "cluster_older_within_euk",
         "cluster_older_pre_euk", "incomparable"]
labels = {"consistent": "same node", "consistent_censored_pre_euk": "pre-euk. LCA, HOG at euk. root (censored)",
          "cluster_younger": "cluster LCA younger", "cluster_older_within_euk": "cluster LCA older (within Eukaryota)",
          "cluster_older_pre_euk": "cluster LCA pre-eukaryotic, HOG younger", "incomparable": "incomparable"}
cols = [C1, "#86b6ef", C3, C4, C2, "#c3c2b7"]
fig, ax = plt.subplots(1, 2, figsize=(11, 3.6), dpi=300, gridspec_kw={"width_ratios": [1, 1.1]})
for i, (f, nm) in enumerate([("lca_consistency_clusters.tsv", "clusters\n(vs dominant HOG)"),
                             ("lca_consistency_hogs.tsv", "HOGs\n(vs largest cluster)")]):
    t = pl.read_csv(OUT_R / f, separator="\t")
    left = 0
    for c, col in zip(order, cols):
        v = t.filter(pl.col("cons") == c)["frac"]
        v = float(v[0]) if len(v) else 0
        ax[0].barh(i, v, left=left, color=col, edgecolor="white", lw=1, label=labels[c] if i == 0 else None)
        if v > 0.06:
            ax[0].text(left + v / 2, i, f"{v:.0%}", ha="center", va="center", fontsize=7, color="white")
        left += v
ax[0].set_yticks([0, 1], ["clusters\n(vs dominant HOG)", "HOGs\n(vs largest cluster)"])
ax[0].set_xlim(0, 1); ax[0].set_xlabel("fraction")
ax[0].legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
ax[0].set_title("a  AFDB cluster LCA vs Ancestra HOG root", loc="left", fontsize=9)
a = pl.read_parquet(OUT_I / "cluster_lca_vs_hog_ages.parquet")
AB = [-1, 50, 200, 500, 800, 1100, 1400, 1e9]
AL = ["<50", "50-200", "200-500", "500-800", "800-1100", "1100-1400", "≥1400\n(euk root)"]
def b(v):
    return np.digitize(np.nan_to_num(v, nan=-1), AB[1:-1]) if True else None
ha = a["hog_age"].to_numpy(); ca = a["full_lca_age"].to_numpy()
ok = ~np.isnan(ha) & ~np.isnan(ca)
hb = np.digitize(ha[ok], AB[1:-1]); cbn = np.where(np.isinf(ca[ok]), len(AL), np.digitize(ca[ok], AB[1:-1]))
M = np.zeros((len(AL), len(AL) + 1))
np.add.at(M, (hb, cbn), 1)
Mr = M / M.sum(1, keepdims=True)
im = ax[1].imshow(Mr, cmap="Blues", vmin=0, vmax=1, aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        if Mr[i, j] >= 0.05:
            ax[1].text(j, i, f"{Mr[i, j]:.2f}", ha="center", va="center", fontsize=6,
                       color="white" if Mr[i, j] > 0.5 else INK)
ax[1].set_xticks(range(len(AL) + 1), [l.replace("\n", " ") for l in AL] + ["pre-eukaryotic"], fontsize=6.5, rotation=35, ha="right")
ax[1].set_yticks(range(len(AL)), [f"{l.splitlines()[0]} (n={int(M[i].sum()):,})" for i, l in enumerate(AL)], fontsize=6.5)
ax[1].set_xlabel("AFDB cluster LCA, projected on the panel (Mya)"); ax[1].set_ylabel("dominant HOG root age (Mya)")
ax[1].set_title("b  Row-normalised, per cluster", loc="left", fontsize=9)
plt.colorbar(im, ax=ax[1], fraction=0.04, label="row fraction")
fig.tight_layout(); fig.savefig(OUT_F / "cluster_lca_vs_hog_root.png", bbox_inches="tight")
print("ok", int(ok.sum()), "clusters in heatmap")
