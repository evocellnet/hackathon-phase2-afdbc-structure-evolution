"""X6 figures: (a) ML trees (AA, LG+G4, midpoint-rooted for display only) with subfamily and domain strips;
(b) per-branch b_3Di vs b_AA on the fixed topology."""
import json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from ete3 import Tree

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
OUT = f"{R}/results/casestudy"
FIG = f"{R}/figures/casestudy"
import os; os.makedirs(FIG, exist_ok=True)
SRC = "linsi"
sel = pd.read_csv(f"{OUT}/sample_members.tsv", sep="\t").set_index("seqid")

DOM_COL = {"Eukaryota": "#4C72B0", "Asgard": "#C44E52", "Archaea": "#DD8452", "Bacteria": "#55A868"}
PAL = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
       "#17becf", "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173", "#3182bd", "#e6550d", "#31a354",
       "#756bb1", "#636363", "#9c9ede", "#cedb9c", "#e7ba52", "#ad494a", "#a55194", "#6baed6", "#fd8d3c"]


def coarse(sub):
    return sub


def draw_tree(ax, t, sf):
    t.set_outgroup(t.get_midpoint_outgroup())
    t.ladderize()
    leaves = t.get_leaves()
    y = {}
    for i, l in enumerate(leaves):
        y[l] = i
    x = {}
    for n in t.traverse("preorder"):
        x[n] = 0.0 if n.is_root() else x[n.up] + n.dist
    for n in t.traverse("postorder"):
        if not n.is_leaf():
            y[n] = np.mean([y[c] for c in n.children])
    for n in t.traverse():
        if n.is_root():
            continue
        ax.plot([x[n.up], x[n]], [y[n], y[n]], color="k", lw=0.4)
        if not n.is_leaf():
            ys = [y[c] for c in n.children]
            ax.plot([x[n], x[n]], [min(ys), max(ys)], color="k", lw=0.4)
            try:
                sup = float(str(n.name).split("/")[-1])
            except ValueError:
                sup = np.nan
            if sup >= 95:
                ax.plot(x[n], y[n], "o", ms=1.4, color="#333333")
    if not t.is_leaf():
        ys = [y[c] for c in t.children]
        ax.plot([0, 0], [min(ys), max(ys)], color="k", lw=0.4)
    xmax = max(x.values())
    subs = [sel.at[l.name, "subfamily"] for l in leaves]
    doms = [sel.at[l.name, "domain"] for l in leaves]
    usub = sorted(set(subs))
    scol = {s: PAL[i % len(PAL)] for i, s in enumerate(usub)}
    x1, x2 = xmax * 1.03, xmax * 1.08
    for i, (s, d) in enumerate(zip(subs, doms)):
        ax.add_patch(plt.Rectangle((x1, i - 0.5), xmax * 0.04, 1, color=DOM_COL[d], lw=0))
        ax.add_patch(plt.Rectangle((x2, i - 0.5), xmax * 0.04, 1, color=scol[s], lw=0))
    # run labels
    i = 0
    while i < len(subs):
        j = i
        while j + 1 < len(subs) and subs[j + 1] == subs[i]:
            j += 1
        if j - i + 1 >= 3 or doms[i] != "Eukaryota":
            ax.text(xmax * 1.14, (i + j) / 2, f"{subs[i]} ({j - i + 1})", fontsize=4.2, va="center")
        i = j + 1
    ax.set_ylim(-1, len(leaves))
    ax.set_xlim(-xmax * 0.01, xmax * 1.55)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.set_xlabel("substitutions/site (LG+G4)")
    ax.set_title(f"{sf} superfamily: AA ML tree ({len(leaves)} seqs; dots: UFBoot>=95; midpoint root, display only)",
                 fontsize=8)
    h = [Patch(color=c, label=k) for k, c in DOM_COL.items()]
    ax.legend(handles=h, loc="lower left", fontsize=6, frameon=False, title="domain strip", title_fontsize=6)


for sf in ["actin", "tubulin"]:
    t = Tree(f"{CS}/iqtree/{sf}_{SRC}_aa_ml.treefile", format=1)
    fig, ax = plt.subplots(figsize=(8, 16))
    draw_tree(ax, t, sf)
    fig.tight_layout()
    fig.savefig(f"{FIG}/{sf}_ml_tree.pdf"); fig.savefig(f"{FIG}/{sf}_ml_tree.png", dpi=200)
    plt.close(fig)

# ---- branch scatter -------------------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(10, 9.2))
CCOL = {"within_eukaryote": "#4C72B0", "within_prok_or_mixed": "#55A868", "mixed_split": "#999999",
        "EUK_STEM": "#C44E52", "EUK+ASGARD_STEM": "#8172B2"}
for ax, (src, sf) in zip(axes.ravel(), [(a, b_) for a in ("linsi", "foldmason") for b_ in ("actin", "tubulin")]):
    SRC_ = src
    if not os.path.exists(f"{OUT}/branches_{sf}_{SRC_}_binary.tsv"):
        ax.text(0.5, 0.5, f"{sf} / FoldMason: not run\n(FoldMason left 67 prokaryotic sequences,\n50 FtsA + 15 MreB + 2 salactin,\nin separate columns: all-gap after the occupancy filter)",
                ha="center", va="center", fontsize=8, transform=ax.transAxes); ax.set_axis_off()
        continue
    b = pd.read_csv(f"{OUT}/branches_{sf}_{SRC_}_binary.tsv", sep="\t").fillna({"clade_label": "", "focal": ""})
    s = json.load(open(f"{OUT}/branch_summary_{sf}_{SRC_}_binary.json"))
    term = b[b.terminal]
    ax.scatter(term.b_aa, term.b_3di, s=3, c="#cccccc", label="terminal", lw=0)
    inn = b[~b.terminal]
    ax.scatter(inn[~inn.supported].b_aa, inn[~inn.supported].b_3di, s=4, facecolors="none",
               edgecolors="#bbbbbb", lw=0.4, label="internal, UFBoot<95")
    for c, col in [("within_eukaryote", "#4C72B0"), ("within_prok_or_mixed", "#55A868"), ("mixed_split", "#999999")]:
        x = inn[inn.supported & (inn["class"] == c) & (inn.focal == "")]
        if len(x):
            ax.scatter(x.b_aa, x.b_3di, s=6, c=col, lw=0, label=f"{c} (UFBoot>=95)", zorder=3)
    foc = inn[inn.focal != ""].drop_duplicates(["b_aa", "b_3di"])
    for k, (_, r) in enumerate(foc.iterrows()):
        ax.scatter([r.b_aa], [r.b_3di], s=45, c=PAL[(k + 3) % len(PAL)], edgecolors="k", lw=0.6, zorder=6,
                   label=f"{r.focal} (UFBoot {r.ufboot:.0f}, rho {r.rho:.2f})")
    lab = b[(b.clade_label != "") & ~b.terminal & b.supported]
    for _, r in lab.iterrows():
        ax.annotate(r.clade_label.replace("stem:", ""), (r.b_aa, r.b_3di), fontsize=4.5, alpha=0.8)
    lim = np.array([1e-4, 5])
    ax.plot(lim, lim * s["scale_sum3di_over_sumaa"], "k--", lw=0.6, label="tree-wide b3Di/bAA")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(1e-3, 5); ax.set_ylim(1e-3, 10)
    ax.set_xlabel("b_AA (LG+G4, fixed topology)"); ax.set_ylabel("b_3Di (Q.3Di.AF+G4, same columns)")
    ax.set_title(f"{sf}: {'MAFFT L-INS-i' if SRC_ == 'linsi' else 'FoldMason'} columns, fixed AA-ML topology", fontsize=8)
    ax.legend(fontsize=5, frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(f"{FIG}/branch_b3di_vs_baa.pdf"); fig.savefig(f"{FIG}/branch_b3di_vs_baa.png", dpi=200)
print("figures written")
