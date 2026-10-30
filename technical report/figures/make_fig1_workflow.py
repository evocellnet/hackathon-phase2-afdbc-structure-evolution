"""Fig 1: workflow schematic (data -> orthology -> pairs -> metrics -> engines -> endpoints).

All counts are taken from project files (results/hogs/hogs_build_summary.json,
logs/filters/*.tsv, agents/panel_report.md, results/ancestra/euk180/manifest.json).
Run: python "technical report/figures/make_fig1_workflow.py"
"""
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path(__file__).resolve().parent

INK = "#0b0b0b"
INK2 = "#52514e"
EDGE = "#8a8984"
C_DATA = "#dbe8f8"    # light tint of slot-1 blue
C_ORTH = "#d7f0e5"    # light tint of aqua
C_PAIR = "#fde7da"    # light tint of orange
C_ENG = "#e6e2f6"     # light tint of violet
C_END = "#f3f2ee"     # neutral
C_GUARD = "#fff4d6"   # light tint of yellow

fig, ax = plt.subplots(figsize=(13, 7.4))
ax.set_xlim(0, 13)
ax.set_ylim(0, 7.4)
ax.axis("off")


def box(x, y, w, h, title, body, fc):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                       fc=fc, ec=EDGE, lw=1.0)
    ax.add_patch(p)
    ax.text(x + 0.12, y + h - 0.14, title, ha="left", va="top", fontsize=10,
            fontweight="bold", color=INK)
    ax.text(x + 0.12, y + h - 0.47, body, ha="left", va="top", fontsize=8.6,
            color=INK2, linespacing=1.35)
    return (x, y, w, h)


def arrow(x0, y0, x1, y1, rad=0.0):
    a = FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=13,
                        lw=1.3, color=INK2, connectionstyle=f"arc3,rad={rad}")
    ax.add_patch(a)


def elbow(x0, y0, x1, y1, ymid):
    ax.plot([x0, x0, x1], [y0, ymid, ymid], color=INK2, lw=1.3, solid_capstyle="round")
    arrow(x1, ymid, x1, y1)


# ---- row 1: data and orthology -------------------------------------------
y1, h1 = 5.25, 1.9
b1 = box(0.15, y1, 2.4, h1, "1  Species panel",
         "180 eukaryotes (UniProt\n2026_03 reference\nproteomes; 3.46 M proteins)\n"
         "TimeTree ages ->\ndated species tree", C_DATA)
b2 = box(2.75, y1, 2.4, h1, "2  Structures",
         "AlphaFold DB v4 models\n(3.29 M; sequence identical\nto UniProt)\n"
         "per-residue pLDDT;\n3Di letters (Foldseek)", C_DATA)
b3 = box(5.35, y1, 2.4, h1, "3  Orthology",
         "Ancestra, sequence only\n109,598 HOGs (families)\n"
         "ortholog / paralog label\nfor every gene pair", C_ORTH)
b4 = box(7.95, y1, 2.5, h1, "4  Family sampling",
         "filters (pLDDT >= 70, no\nfragments) -> 12,336 HOGs\n"
         "1 gene/species, AA distance\n<= 2 -> 8,211 HOGs\n"
         "7,000 sampled by age", C_ORTH)
b5 = box(10.65, y1, 2.2, h1, "5  Split (once)",
         "random 20 / 80 within\nage strata\n1,400 tuning HOGs\n5,600 confirmation HOGs", C_GUARD)
for (a, b) in [(b1, b2), (b2, b3), (b3, b4), (b4, b5)]:
    arrow(a[0] + a[2] + 0.02, y1 + h1 / 2, b[0] - 0.02, y1 + h1 / 2)

# ---- row 2: pairs and metrics ---------------------------------------------
y2, h2 = 2.75, 1.9
b6 = box(0.15, y2, 3.0, h2, "6  Protein pairs",
         "all pairs within each family:\n1.09 M ortholog pairs\n"
         "+ 0.31 M paralog-control pairs\n(duplicated copies + outgroups)", C_PAIR)
b7 = box(3.4, y2, 3.0, h2, "7  Three alignments",
         "primary: Foldseek 3Di+AA\nsensitivity: sequence-only;\nFoldseek TM-align\n"
         "core = residues with\npLDDT >= 70 in both", C_PAIR)
b8 = box(6.65, y2, 3.2, h2, "8  Divergence per pair",
         "AA: ML distance (LG+G)\n3Di: ML distance (Q.3Di.AF+G)\n"
         "3D: LDDT (floor-corrected),\nTM-score, RMSD\n+ a standard error for each", C_PAIR)
b9 = box(10.1, y2, 2.75, h2, "9  Noise floors",
         "100,000 random unrelated\npairs (chance similarity)\n"
         "17,768 identical-sequence\npairs (AlphaFold noise)", C_PAIR)
elbow(b5[0] + b5[2] / 2, y1 - 0.02, b6[0] + b6[2] / 2, y2 + h2 + 0.02, (y1 + y2 + h2) / 2)
for (a, b) in [(b6, b7), (b7, b8)]:
    arrow(a[0] + a[2] + 0.02, y2 + h2 / 2, b[0] - 0.02, y2 + h2 / 2)
arrow(b9[0] - 0.02, y2 + h2 / 2, b8[0] + b8[2] + 0.02, y2 + h2 / 2)

# ---- row 3: engines and endpoints ------------------------------------------
y3, h3 = 0.15, 1.95
b10 = box(0.15, y3, 3.7, h3, "10a  Pairwise engine",
          "structure dist. = floor + a x AA dist.^c\n"
          "c = 1: linear;  c < 1: saturating\n"
          "corrected for noise on both axes\n-> E1 (3Di), E2 (LDDT), E6 (paralogs)", C_ENG)
b12 = box(4.1, y3, 2.55, h3, "11  Endpoints E1-E6",
          "decision rules fixed\nin advance; Holm\ncorrection over five\n"
          "tests -> interpretation", C_END)
b11 = box(6.9, y3, 3.5, h3, "10b  Branch engine",
          "fit each family onto the dated\nspecies tree -> AA and 3Di\n"
          "branch lengths on every branch\n-> E3 (ratio), E4 (age), E5 (lineages)", C_ENG)
b13 = box(10.65, y3, 2.2, h3, "Safeguards",
          "tune on 20% only, then\none frozen run on 80%\nknown-answer\nsimulations\ndeviation log", C_GUARD)
ymid3 = (y2 + y3 + h3) / 2
ax.plot([b8[0] + 1.6, b8[0] + 1.6], [y2 - 0.02, ymid3], color=INK2, lw=1.3)
ax.plot([b10[0] + b10[2] / 2, b11[0] + b11[2] / 2], [ymid3, ymid3], color=INK2, lw=1.3)
arrow(b10[0] + b10[2] / 2, ymid3, b10[0] + b10[2] / 2, y3 + h3 + 0.02)
arrow(b11[0] + b11[2] / 2, ymid3, b11[0] + b11[2] / 2, y3 + h3 + 0.02)
arrow(b10[0] + b10[2] + 0.02, y3 + h3 / 2, b12[0] - 0.02, y3 + h3 / 2)
arrow(b11[0] - 0.02, y3 + h3 / 2, b12[0] + b12[2] + 0.02, y3 + h3 / 2)

fig.savefig(OUT / "fig1_workflow.png", dpi=200, bbox_inches="tight", facecolor="white")
fig.savefig(OUT / "fig1_workflow.svg", bbox_inches="tight", facecolor="white")
print("wrote", OUT / "fig1_workflow.png")
