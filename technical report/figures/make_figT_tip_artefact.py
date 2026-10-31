#!/usr/bin/env python3
"""Schematic of the terminal-branch ("tip") artefact (report Figure 6). Illustrative only: k = 4.6 (internal-branch
estimate) and an extra 3Di offset of 0.014 per terminal branch (the value that reproduces the observed pattern in the
post-hoc simulations, section 4.3). No project data are plotted."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C1, C2, INK, INK2, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e0"
K, DELTA = 4.6, 0.014
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": INK2})
fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw=dict(width_ratios=[1, 1.15]))

# (a) tree: root -> internal node -> two young tips; a third, older tip from the root
a.set_xlim(0, 10); a.set_ylim(-1.2, 10); a.axis("off")
def seg(x0, y0, x1, y1, **k): a.plot([x0, x1], [y0, y1], solid_capstyle="butt", **k)
seg(0.5, 5, 0.5, 8.0, color=INK, lw=2); seg(0.5, 2.0, 0.5, 5, color=INK, lw=2)
seg(0.5, 8.0, 5.0, 8.0, color=INK, lw=2)            # internal branch
seg(5.0, 6.8, 5.0, 9.2, color=INK, lw=2)
for y in (9.2, 6.8):                                   # young terminal branches: short
    seg(5.0, y, 7.0, y, color=INK, lw=2)
    seg(7.0, y, 8.3, y, color=C2, lw=6)
seg(0.5, 2.0, 7.0, 2.0, color=INK, lw=2)              # old terminal branch: long
seg(7.0, 2.0, 8.3, 2.0, color=C2, lw=6)
a.text(2.75, 8.25, "internal branch", ha="center", va="bottom", color=INK2)
a.text(6.0, 9.45, "young tip", ha="center", va="bottom", color=INK2)
a.text(3.75, 2.25, "old tip", ha="center", va="bottom", color=INK2)
for y in (9.2, 6.8, 2.0):
    a.text(8.5, y, "+ δ", va="center", color=C2, fontsize=10, weight="bold")
a.text(0.2, -1.1, "Black: real change; 3Di length = AA length / k.\nOrange: the same extra 3Di length δ on every tip\n"
       "(one AlphaFold model per tip, beyond the fixed\ncorrection), whatever the tip's age.", fontsize=8.5, color=INK2, va="bottom")
a.set_title("a  Every terminal branch carries the same extra 3Di length", loc="left", fontsize=10)

# (b) apparent k on a branch vs its AA length
x = np.logspace(-3, 0.3, 300)
b.axhline(K, color=C1, lw=2, label=f"internal branches: k = {K}")
b.plot(x, x / (x / K + DELTA), color=C2, lw=2, label=f"terminal branches: AA / (AA/k + δ), δ = {DELTA}")
b.set_xscale("log"); b.set_ylim(0, 5.2); b.set_xlim(1e-3, 2)
b.set_xlabel("AA length of the branch (substitutions/site; short = young)")
b.set_ylabel("Apparent ratio k = AA / 3Di")
b.grid(color=GRID, lw=0.6)
b.annotate("young tips: δ dominates,\nk looks small", xy=(0.01, 0.01 / (0.01 / K + DELTA)), xytext=(0.0013, 2.6),
           arrowprops=dict(arrowstyle="->", color=INK2), fontsize=8.5, color=INK2)
b.annotate("old tips: δ is negligible", xy=(0.8, 0.8 / (0.8 / K + DELTA)), xytext=(0.0015, 3.95),
           arrowprops=dict(arrowstyle="->", color=INK2), fontsize=8.5, color=INK2)
b.legend(loc="lower right", frameon=False, fontsize=8.5)
b.set_title("b  So k appears to rise with branch age", loc="left", fontsize=10)
fig.text(0.5, -0.04, "Most young branches are terminal, so an age trend fitted over all branches turns negative "
         "(3Di/AA falls with age) even if the true ratio is constant.", ha="center", fontsize=8.5, color=INK2)
fig.tight_layout()
fig.savefig("figT_tip_artefact.png", dpi=250, bbox_inches="tight")
