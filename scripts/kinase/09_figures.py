"""X7 figures.
 (a) figures/kinase/pkl_tree_rho.{png,pdf}: AA ML tree (rooted on the ELK/atypical clade, ladderized; for display),
     each branch coloured by log2 rho (diverging blue-grey-red, clipped at +-2) from the primary (MAFFT L-INS-i)
     fixed-topology run; strips: tip category and eukaryotic ePK group; dots = UFBoot >= 95.
 (b) figures/kinase/rho_by_group.{png,pdf}: rho distributions by branch class (terminal / supported internal) and,
     for eukaryotic ePK groups, terminal-branch rho per group, with the group-stem rho marked.
 (c) figures/kinase/rho_alignment_sensitivity.{png,pdf}: rho of supported internal branches, L-INS-i vs FAMSA and
     vs FoldMason (same topology).
"""
import json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch
from ete3 import Tree

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
K = f"{R}/data/interim/kinase"
OUT = f"{R}/results/kinase"
FIG = f"{R}/figures/kinase"
sel = pd.read_csv(f"{OUT}/sample_members.tsv", sep="\t").set_index("seqid")
INK, INK2 = "#0b0b0b", "#52514e"
CMAP = LinearSegmentedColormap.from_list("div", ["#184f95", "#6da7ec", "#c9c8c4", "#ec8f8e", "#b52d2d"])
NORM = Normalize(-2, 2)
CATC = {"euk_ePK": "#2a78d6", "prok_ePK": "#eb6834", "euk_ELK": "#1baf7a", "prok_ELK": "#eda100"}
GRPS = ["AGC", "CAMK", "CK1", "CMGC", "STE", "TK", "TKL", "Other", "TKlike_nonholozoan"]
GRPC = dict(zip(GRPS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#e34948", "#4a3aa7", "#8a8984",
                       "#c9c8c4"]))


def cat(k):
    c = sel.at[k, "klass"]
    if c.startswith("ePK_"):
        return "euk_ePK"
    if c in ("prok_ePK", "Asgard_ePK"):
        return "prok_ePK"
    return "euk_ELK" if c.endswith("_euk") else "prok_ELK"


br = pd.read_csv(f"{OUT}/branches_linsi.tsv", sep="\t")
t = Tree(f"{K}/iqtree/pkl_linsi_aa_ml.treefile", format=1)
leaves = t.get_leaf_names()
ALL = frozenset(leaves)
CAT = {k: cat(k) for k in leaves}
t.set_outgroup(t & next(k for k in leaves if CAT[k] == "euk_ePK"))
t.set_outgroup(t.get_common_ancestor([k for k in leaves if CAT[k].endswith("ELK")]))
t.ladderize()
aaf = Tree(f"{K}/iqtree/pkl_linsi_aa_fixed.treefile", format=1)
tdf = Tree(f"{K}/iqtree/pkl_linsi_3di_fixed.treefile", format=1)
bip = lambda s: frozenset(s) if leaves[0] not in s else ALL - frozenset(s)
La = {bip(n.get_leaf_names()): n.dist for n in aaf.traverse() if not n.is_root()}
Lt = {bip(n.get_leaf_names()): n.dist for n in tdf.traverse() if not n.is_root()}
scale = sum(Lt.values()) / sum(La.values())

fig, ax = plt.subplots(figsize=(9, 16))
lv = t.get_leaves()
y = {l: i for i, l in enumerate(lv)}
x = {}
for n in t.traverse("preorder"):
    x[n] = 0.0 if n.is_root() else x[n.up] + La[bip(n.get_leaf_names())]
for n in t.traverse("postorder"):
    if not n.is_leaf():
        y[n] = np.mean([y[c] for c in n.children])
for n in t.traverse():
    if n.is_root():
        continue
    b = bip(n.get_leaf_names())
    r = np.log2(max(Lt[b] / La[b] / scale, 1e-3)) if La[b] > 1e-6 else 0
    ax.plot([x[n.up], x[n]], [y[n], y[n]], color=CMAP(NORM(np.clip(r, -2, 2))), lw=1.1 if not n.is_leaf() else 0.7,
            solid_capstyle="butt")
    if not n.is_leaf():
        ys = [y[c] for c in n.children]
        ax.plot([x[n], x[n]], [min(ys), max(ys)], color="#c9c8c4", lw=0.4)
        try:
            if float(str(n.name).split("/")[-1]) >= 95:
                ax.plot(x[n], y[n], "o", ms=1.6, color=INK2)
        except ValueError:
            pass
ys = [y[c] for c in t.children]
ax.plot([0, 0], [min(ys), max(ys)], color="#c9c8c4", lw=0.4)
xm = max(x.values())
for l in lv:
    k = l.name
    ax.add_patch(plt.Rectangle((xm * 1.02, y[l] - 0.5), xm * 0.03, 1, color=CATC[CAT[k]], lw=0))
    g = sel.at[k, "klass"][4:] if CAT[k] == "euk_ePK" else ""
    if g:
        ax.add_patch(plt.Rectangle((xm * 1.06, y[l] - 0.5), xm * 0.03, 1, color=GRPC[g], lw=0))
    else:
        fam = sel.at[k, "klass"].rsplit("_", 1)[0] if CAT[k].endswith("ELK") else ""
        if fam and (l == lv[0] or sel.at[lv[y[l] - 1].name, "klass"].rsplit("_", 1)[0] != fam):
            ax.text(xm * 1.10, y[l], fam.replace("Kdo", "Kdo/Bud32"), fontsize=5, va="center", color=INK2)
# group labels at the median y of each group block
for g in GRPS:
    yy = [y[l] for l in lv if CAT[l.name] == "euk_ePK" and sel.at[l.name, "klass"][4:] == g]
    if yy and g != "TKlike_nonholozoan":
        ax.text(xm * 1.10, np.median(yy), f"{g} ({len(yy)})", fontsize=6, va="center", color=INK)
ax.set_ylim(-1, len(lv)); ax.invert_yaxis(); ax.set_yticks([])
ax.set_xlim(-0.02 * xm, xm * 1.35)
ax.set_xlabel("AA substitutions / site (LG+G4, fixed topology)", color=INK2)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
sm = plt.cm.ScalarMappable(cmap=CMAP, norm=NORM)
cb = fig.colorbar(sm, ax=ax, fraction=0.015, pad=0.04, shrink=0.3, location="bottom")
cb.set_label("log2 rho  (b_3Di / b_AA, relative to tree-wide ratio)", color=INK2, fontsize=8)
ax.legend(handles=[Patch(color=v, label=k) for k, v in CATC.items()], loc="lower left", fontsize=7, frameon=False,
          title="tip category (inner strip)", title_fontsize=7)
ax.set_title("PKL kinase domains (n = 500): branch colour = 3Di vs AA change; dots UFBoot >= 95", fontsize=9,
             color=INK, loc="left")
fig.tight_layout()
fig.savefig(f"{FIG}/pkl_tree_rho.png", dpi=200); fig.savefig(f"{FIG}/pkl_tree_rho.pdf")
plt.close(fig)

# (b) distributions
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[1, 1.3]))
ax = axs[0]
order = [("terminal_euk_ePK", "terminal\neuk ePK"), ("terminal_prok_ePK", "terminal\nprok ePK"),
         ("terminal_euk_ELK", "terminal\neuk ELK"), ("terminal_prok_ELK", "terminal\nprok ELK"),
         ("within_euk_ePK", "internal\neuk ePK"), ("within_prok_ePK", "internal\nprok ePK"),
         ("within_euk_ELK", "internal\neuk ELK"), ("within_prok_ELK", "internal\nprok ELK"), ("mixed", "internal\nmixed")]
sub = br[br.supported & (br.b_aa >= 0.02)]
data, labs, cols = [], [], []
for c, l in order:
    v = np.log2(sub[sub.bclass == c].rho).clip(lower=-4)
    if len(v):
        data.append(v); labs.append(f"{l}\nn={len(v)}")
        cols.append(CATC.get(c.replace("terminal_", "").replace("within_", ""), "#8a8984"))
bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.2))
for p_, c in zip(bp["boxes"], cols):
    p_.set_facecolor(c); p_.set_alpha(0.55); p_.set_edgecolor(c)
for i, v in enumerate(data):
    ax.scatter(np.random.default_rng(1).normal(i + 1, 0.06, len(v)), v, s=3, color=INK2, alpha=0.35, lw=0)
ax.axhline(0, color="#c9c8c4", lw=0.8, zorder=0)
ax.set_xticks(range(1, len(labs) + 1)); ax.set_xticklabels(labs, fontsize=6.5)
ax.set_ylabel("log2 rho (clipped at -4; b_3Di ~ 0)", color=INK2)
ax.set_title("Branch classes (supported, b_AA >= 0.02)", fontsize=9, loc="left"); ax.set_ylim(-4.3, 3.2)
ax = axs[1]
st = pd.read_csv(f"{OUT}/group_stems_linsi.tsv", sep="\t")
term = br[br.terminal & (br.bclass == "terminal_euk_ePK")].copy()
term["grp"] = [sel.at[k, "klass"][4:] for k in term.tip]
data, labs = [], []
for i, g in enumerate(GRPS):
    v = np.log2(term[term.grp == g].rho)
    data.append(v); labs.append(g.replace("TKlike_nonholozoan", "TK-like\n(non-holo)") + f"\nn={len(v)}")
bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False, medianprops=dict(color=INK, lw=1.2))
for p_, g in zip(bp["boxes"], GRPS):
    p_.set_facecolor(GRPC[g]); p_.set_alpha(0.5); p_.set_edgecolor(GRPC[g])
for i, g in enumerate(GRPS):
    for nm, mk in ((f"{g} MRCA", "D"), (f"{g} (largest pure clade)", "^")):
        r = st[st.stem == nm]
        if len(r) and pd.notna(r.rho.iloc[0]) and r.purity.iloc[0] >= 0.8:
            ax.scatter(i + 1.3, np.log2(r.rho.iloc[0]), marker=mk, s=28, color=INK, zorder=5)
r = st[st.stem == "TK holozoan MRCA"]
if len(r) and pd.notna(r.rho.iloc[0]):
    ax.scatter(GRPS.index("TK") + 1.3, np.log2(r.rho.iloc[0]), marker="*", s=70, color="#e34948", zorder=6,
               edgecolor=INK, lw=0.5)
ax.axhline(0, color="#c9c8c4", lw=0.8, zorder=0)
bge = sub[sub.bclass == "within_euk_ePK"].rho
ax.axhspan(np.log2(bge.quantile(.25)), np.log2(bge.quantile(.75)), color="#c9c8c4", alpha=0.3, zorder=0)
ax.set_xticks(range(1, len(labs) + 1)); ax.set_xticklabels(labs, fontsize=6.5)
ax.set_ylabel("log2 rho", color=INK2)
ax.set_ylim(-4.3, 3.2)
ax.set_title("Euk ePK groups: terminal branches (boxes); stems: diamond = MRCA (purity >= 0.8),\n"
             "triangle = largest pure clade, star = holozoan TK MRCA; grey band = IQR of supported internal euk-ePK branches", fontsize=8, loc="left")
for a in axs:
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(f"{FIG}/rho_by_group.png", dpi=200); fig.savefig(f"{FIG}/rho_by_group.pdf")
plt.close(fig)

# (c) alignment sensitivity
srcs = [s for s in ("famsa", "foldmason") if __import__("os").path.exists(f"{OUT}/branches_{s}.tsv")]
if srcs:
    fig, axs = plt.subplots(1, len(srcs), figsize=(4.6 * len(srcs), 4.4), squeeze=False)
    a0 = br.set_index(["clade_size", "b_aa"])  # placeholder to keep pandas happy
    for ax, s in zip(axs[0], srcs):
        o = pd.read_csv(f"{OUT}/branches_{s}.tsv", sep="\t")
        assert len(o) == len(br)
        m = br.supported & ~br.terminal & (br.b_aa >= 0.02)
        xv, yv = np.log2(br.rho[m]).clip(-5, 5), np.log2(o.rho[m]).clip(-5, 5)
        c = [CATC.get(k.replace("within_", ""), "#8a8984") for k in br.bclass[m]]
        ax.scatter(xv, yv, s=9, c=c, alpha=0.8, lw=0)
        lim = [min(xv.min(), yv.min()) - 0.2, max(xv.max(), yv.max()) + 0.2]
        ax.plot(lim, lim, color="#c9c8c4", lw=0.8)
        rs = pd.Series(xv.values).corr(pd.Series(yv.values), method="spearman")
        ax.set_xlabel("log2 rho, MAFFT L-INS-i columns (clipped to +-5)", color=INK2)
        ax.set_ylabel(f"log2 rho, {s} columns", color=INK2)
        ax.set_title(f"supported internal branches (n={m.sum()}), Spearman {rs:.2f}", fontsize=8, loc="left")
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
    axs[0][0].legend(handles=[Patch(color=v, label=k) for k, v in CATC.items()] +
                     [Patch(color="#8a8984", label="mixed")], fontsize=6.5, frameon=False)
    fig.tight_layout()
    fig.savefig(f"{FIG}/rho_alignment_sensitivity.png", dpi=200); fig.savefig(f"{FIG}/rho_alignment_sensitivity.pdf")
print("figures written")
