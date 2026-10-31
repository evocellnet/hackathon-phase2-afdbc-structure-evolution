"""X6 node 2 figures (figures/casestudy/):
  {sf}_tree_rho.{pdf,png}   ML tree (outgroup-rooted for display; x = IQ-TREE AA length, MAFFT columns) with branches
                            coloured by log2 rho (3Di/AA relative ratio); grey = not comparable (UFBoot<95 or
                            b_AA<0.02), dashed black = saturated (b_AA>2); stars = robust flags; named clades labelled.
  rho_distributions.{pdf,png}  rho vs b_AA for comparable branches (internal / terminal), comparable-length 95 % window,
                            named clades with bootstrap 95 % CIs.
  method_comparison.{pdf,png}  IQ-TREE (MAFFT columns) vs branch engine (pairwise tdiaa) b_AA, b_3Di and rho, and rho of
                            the named clades across alignments / methods.
"""
import json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.collections import LineCollection
from ete3 import Tree

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
OUT = f"{R}/results/casestudy"
FIG = f"{R}/figures/casestudy"
CMAP = LinearSegmentedColormap.from_list("div", ["#1f5fb4", "#8fb3e3", "#e4e3df", "#eb9a7a", "#c2361b"])
NORM = Normalize(-2, 2)
DOM_COL = {"Eukaryota": "#2a78d6", "Asgard": "#e87ba4", "Archaea": "#eda100", "Bacteria": "#1baf7a"}
plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False})
sel_all = pd.read_csv(f"{OUT}/sample_members.tsv", sep="\t")


def short(nm):
    return nm.split(" stem")[0].replace(" (F1=1.00)", "")


def tree_fig(sf):
    sel = sel_all[sel_all.superfamily == sf].set_index("seqid")
    leaves = sorted(sel.index); bit = {l: 1 << i for i, l in enumerate(leaves)}
    B = pd.read_csv(f"{OUT}/branch_ratios_{sf}_ml.tsv", sep="\t")
    B["mm"] = [int(h, 16) for h in B.clade_hex]
    info = B.set_index("mm")
    OG = json.load(open(f"{OUT}/branch_ratio_summary_{sf}_ml.json"))["outgroup_leaf_for_orientation"]
    t = Tree(f"{CS}/iqfix/{sf}_linsi_ml_aa_full.treefile", format=1)
    t.set_outgroup(t & OG); t.ladderize()
    ALL = (1 << len(leaves)) - 1
    for n in t.traverse("postorder"):
        n.m = bit[n.name] if n.is_leaf() else sum(c.m for c in n.children)
    lv = t.get_leaves(); y = {l: i for i, l in enumerate(lv)}
    x = {}
    for n in t.traverse("preorder"):
        x[n] = 0.0 if n.is_root() else x[n.up] + n.dist
    for n in t.traverse("postorder"):
        if not n.is_leaf():
            y[n] = np.mean([y[c] for c in n.children])
    fig, ax = plt.subplots(figsize=(8.5, 15))
    segs, cols, lws, dsegs, stars, labs = [], [], [], [], [], []
    for n in t.traverse():
        if n.is_root():
            continue
        m = n.m if not (n.m & bit[OG]) else ALL ^ n.m
        r = info.loc[m] if m in info.index else None
        seg = [(x[n.up], y[n]), (x[n], y[n])]
        if r is None:
            col, lw = "#bbbbbb", 0.4
        else:
            fl = r.iq_linsi_flag
            comp = bool(r.supported) and 0.02 <= r.iq_linsi_b_aa <= 2
            if fl == "saturated":
                dsegs.append(seg); col = None
            elif comp:
                col, lw = CMAP(NORM(np.log2(r.iq_linsi_rho))), 1.3
            else:
                col, lw = "#c8c8c8", 0.4
            if bool(r.get("robust", False)):
                stars.append(((x[n.up] + x[n]) / 2, y[n], "#c2361b" if fl == "high" else "#1f5fb4"))
            if isinstance(r.clade_name, str) and r.clade_name and r.n_clade > 1:
                labs.append(((x[n.up] + x[n]) / 2, y[n], short(r.clade_name), r.iq_linsi_rho))
        if col is not None:
            segs.append(seg); cols.append(col); lws.append(lw)
        if not n.is_leaf():
            ys = [y[c] for c in n.children]
            segs.append([(x[n], min(ys)), (x[n], max(ys))]); cols.append("#888888"); lws.append(0.3)
    ys = [y[c] for c in t.children]
    segs.append([(0, min(ys)), (0, max(ys))]); cols.append("#888888"); lws.append(0.3)
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=lws))
    ax.add_collection(LineCollection(dsegs, colors="k", linewidths=1.0, linestyles="dashed"))
    for sx, sy, c in stars:
        ax.plot(sx, sy, marker="*", ms=7, color=c, mec="k", mew=0.3, zorder=5)
    for lx, ly, s, rho in labs:
        ax.text(lx, ly - 0.8, f"{s} ({rho:.2f})", fontsize=5, ha="center", va="bottom", color="#222222")
    xmax = max(x.values())
    for i, l in enumerate(lv):
        ax.add_patch(plt.Rectangle((xmax * 1.02, i - 0.5), xmax * 0.03, 1, color=DOM_COL[sel.at[l.name, "domain"]], lw=0))
    subs = [sel.at[l.name, "subfamily"] for l in lv]
    i = 0
    while i < len(subs):
        j = i
        while j + 1 < len(subs) and subs[j + 1] == subs[i]:
            j += 1
        if j - i + 1 >= 3 or sel.at[lv[i].name, "domain"] != "Eukaryota":
            ax.text(xmax * 1.07, (i + j) / 2, f"{subs[i]} ({j - i + 1})", fontsize=4.2, va="center")
        i = j + 1
    ax.set_ylim(len(lv), -1); ax.set_xlim(-xmax * 0.01, xmax * 1.35)
    ax.set_yticks([]); ax.spines["left"].set_visible(False)
    ax.set_xlabel("AA substitutions/site (IQ-TREE LG+G4, fixed ML topology, MAFFT columns)")
    ax.set_title(f"{sf}: relative 3Di/AA ratio per branch (rho; colour = log2 rho; outgroup-rooted for display)", fontsize=8)
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=NORM)
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.01, shrink=0.3)
    cb.set_label("log2 rho (3Di/AA relative to tree)")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    h = [Line2D([], [], color="#c8c8c8", lw=0.6, label="UFBoot<95 or b_AA<0.02"),
         Line2D([], [], color="k", ls="--", lw=1, label="b_AA>2 (saturated, excluded)"),
         Line2D([], [], marker="*", color="#c2361b", ls="", mec="k", label="robust high flag"),
         Line2D([], [], marker="*", color="#1f5fb4", ls="", mec="k", label="robust low flag")] + \
        [Patch(color=c, label=k) for k, c in DOM_COL.items()]
    ax.legend(handles=h, loc="lower left", frameon=False, fontsize=6)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIG}/{sf}_tree_rho.{ext}", dpi=220)
    plt.close(fig)


def dist_fig():
    fig, axs = plt.subplots(2, 2, figsize=(9, 7.5))
    for i, sf in enumerate(("actin", "tubulin")):
        B = pd.read_csv(f"{OUT}/branch_ratios_{sf}_ml.tsv", sep="\t")
        comp = B.supported & B.iq_linsi_b_aa.between(0.02, 2)
        # display clip: b_3Di at IQ-TREE's lower bound (1e-6) gives rho ~1e-5; drawn at the axis floor
        for c_ in ("iq_linsi_rho", "iq_linsi_rho_lo", "iq_linsi_rho_hi", "iq_linsi_win_lo", "iq_linsi_win_hi"):
            B[c_] = B[c_].clip(0.012, 80)
        for j, term in enumerate((False, True)):
            ax = axs[i, j]
            c = comp & (B.terminal == term)
            d = B[c].sort_values("iq_linsi_b_aa")
            ax.scatter(d.iq_linsi_b_aa, d.iq_linsi_rho, s=6, color="#8c8c8c", alpha=0.6, lw=0, label="comparable branches")
            ax.plot(d.iq_linsi_b_aa, d.iq_linsi_win_lo.rolling(9, center=True, min_periods=1).median(), color="#555555", lw=1)
            ax.plot(d.iq_linsi_b_aa, d.iq_linsi_win_hi.rolling(9, center=True, min_periods=1).median(), color="#555555", lw=1,
                    label="95 % range, 50 nearest-length branches")
            for fl, colr in (("high", "#c2361b"), ("low", "#1f5fb4")):
                f = d[d.iq_linsi_flag == fl]
                ax.scatter(f.iq_linsi_b_aa, f.iq_linsi_rho, s=10, color=colr, lw=0, label=f"flag {fl}")
            if "robust" in B:
                rb = d[d.robust]
                ax.scatter(rb.iq_linsi_b_aa, rb.iq_linsi_rho, s=40, facecolor="none", edgecolor="k", lw=0.7, label="robust")
            if not term:
                nm = B[(B.clade_name.fillna("") != "") & ~B.terminal]
                for _, r in nm.iterrows():
                    ax.errorbar(r.iq_linsi_b_aa, r.iq_linsi_rho,
                                yerr=[[r.iq_linsi_rho - r.iq_linsi_rho_lo], [r.iq_linsi_rho_hi - r.iq_linsi_rho]],
                                fmt="D", ms=3.5, color="#eda100" if r.iq_linsi_b_aa <= 2 else "k", ecolor="#555555", lw=0.6, zorder=4)
                    if r.n_clade >= 3 and (r.iq_linsi_b_aa > 0.15 or abs(np.log(r.iq_linsi_rho)) > 0.6):
                        ax.annotate(short(r.clade_name), (r.iq_linsi_b_aa, r.iq_linsi_rho), fontsize=5, xytext=(3, 2),
                                    textcoords="offset points")
            ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(0.01, 100)
            ax.axvline(2, color="k", ls=":", lw=0.6)
            ax.set_xlabel("b_AA (subs/site, IQ-TREE, MAFFT columns)"); ax.set_ylabel("rho = (b_3Di/b_AA) / tree ratio")
            ax.set_title(f"{sf}: {'terminal' if term else 'internal (UFBoot>=95)'} branches (n={int(c.sum())})")
            if i == 0 and j == 0:
                ax.legend(fontsize=5.5, frameon=False, loc="lower left")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIG}/rho_distributions.{ext}", dpi=200)
    plt.close(fig)


def method_fig():
    """Report version (Figure 11): 2 x 3 panels, Spearman rank correlation, plain-language labels."""
    from scipy.stats import spearmanr
    fig, axs = plt.subplots(2, 3, figsize=(10, 6.6))
    names = {"b_AA": "AA branch length", "b_3Di": "3Di branch length", "rho": "3Di/AA ratio relative to tree (\u03c1)"}
    for i, sf in enumerate(("actin", "tubulin")):
        B = pd.read_csv(f"{OUT}/branch_ratios_{sf}_ml.tsv", sep="\t")
        c = B.supported & ~B.terminal & B.eng_tdiaa_b_aa.notna()
        for j, (a, b, lab) in enumerate((("iq_linsi_b_aa", "eng_tdiaa_b_aa", "b_AA"), ("iq_linsi_b_3di", "eng_tdiaa_b_3di", "b_3Di"),
                                         ("iq_linsi_rho", "eng_tdiaa_rho", "rho"))):
            ax = axs[i, j]
            d = B[c & (B[a] > 0) & (B[b] > 0)]
            ax.scatter(d[a], d[b], s=8, color="#2a78d6", alpha=0.6, lw=0)
            lim = [min(d[a].min(), d[b].min()), max(d[a].max(), d[b].max())]
            ax.plot(lim, lim, color="#888888", lw=0.6)
            ax.set_xscale("log"); ax.set_yscale("log")
            rs = spearmanr(d[a], d[b]).correlation
            print(sf, lab, "n", len(d), "spearman", round(rs, 3))
            ax.set_title(f"{sf.capitalize()}: {names[lab]}\nSpearman {rs:.2f}, n = {len(d)} internal branches", fontsize=8)
            ax.set_xlabel("IQ-TREE on fixed tree (MAFFT alignment)", fontsize=7)
            ax.set_ylabel("Pairwise branch engine (3Di+AA alignment)", fontsize=7)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(f"{FIG}/method_comparison.{ext}", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    for sf in ("actin", "tubulin"):
        tree_fig(sf)
    dist_fig()
    method_fig()
