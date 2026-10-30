"""Fig 2: species panel overview.

(a) species per group; (b) ages of the internal nodes of the dated, collapsed panel tree;
(c) ages of all dated species pairs; (d) HOGs per root-age stratum (eligible vs sampled, tune/confirm).

Reads TimeTree-derived tables from data/interim/timetree/ (gitignored, never committed) and
plots only binned counts; no per-pair or per-node age values are written to any file.
Run from the project root: python "technical report/figures/make_fig2_panel.py"
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
GREY = "#c9c8c2"

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
                     "axes.spines.right": False})

panel = pd.read_csv(ROOT / "data/processed/panel/species_panel.tsv", sep="\t")
nodes = pd.read_csv(ROOT / "data/interim/timetree/node_ages.tsv", sep="\t")
pairs = pd.read_csv(ROOT / "data/interim/timetree/pairwise_ages.tsv", sep="\t")
samp = pd.read_csv(ROOT / "results/hogs/hogs_sampling.tsv", sep="\t")

labels = {"Metazoa": "Metazoa (animals)", "Fungi": "Fungi", "Land_plants": "Land plants",
          "Alveolata": "Alveolata", "Discoba": "Discoba", "Stramenopiles": "Stramenopiles",
          "Chlorophyta": "Chlorophyta", "Amoebozoa": "Amoebozoa",
          "Holozoa_other": "Other holozoans", "Rhodophyta": "Rhodophyta (red algae)",
          "Metamonada": "Metamonada", "Charophytes": "Charophytes", "Rhizaria": "Rhizaria",
          "Haptista": "Haptista", "Cryptista": "Cryptista"}

fig, axs = plt.subplots(2, 2, figsize=(11, 7.6))

# (a) species per group
ax = axs[0, 0]
vc = panel["group"].value_counts().sort_values()
ax.barh([labels.get(g, g) for g in vc.index], vc.values, color=BLUE, height=0.7)
for i, v in enumerate(vc.values):
    ax.text(v + 0.8, i, str(v), va="center", fontsize=8, color=INK2)
ax.set_xlabel("species")
ax.set_title(f"a  Species panel ({len(panel)} eukaryotes)", loc="left", fontweight="bold",
             color=INK)
ax.grid(axis="x", color=GRID, lw=0.6)
ax.set_axisbelow(True)

bins = np.logspace(np.log10(0.3), np.log10(2000), 23)


def loghist(ax, values, color, title, ylabel):
    values = values[np.isfinite(values) & (values > 0)]
    ax.hist(values, bins=bins, color=color, edgecolor="white", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("divergence time (million years ago, log scale)")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", fontweight="bold", color=INK)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    return len(values)


# (b) internal nodes of the collapsed panel tree
ax = axs[0, 1]
inner = nodes.loc[~nodes["is_leaf"].astype(bool), "median"].astype(float).to_numpy()
n_dated = loghist(ax, inner, AQUA, "b  Splits in the dated species tree", "internal nodes")
ax.text(0.02, 0.95, f"{n_dated} dated internal nodes\n({len(inner) - n_dated} undated)",
        transform=ax.transAxes, va="top", fontsize=8, color=INK2)

# (c) all dated species pairs
ax = axs[1, 0]
pa = pairs["age_mya"].astype(float).to_numpy()
n_pairs = loghist(ax, pa, ORANGE, "c  All species pairs", "species pairs (log scale)")
ax.set_yscale("log")
ax.text(0.02, 0.95, f"{n_pairs:,} dated pairs of {len(pairs):,}\n"
        "most pairs share a few deep splits:\nanalyse by branch, not by pair",
        transform=ax.transAxes, va="top", fontsize=8, color=INK2)

# (d) HOGs per root-age stratum
ax = axs[1, 1]
el = samp[samp.step == "final_quota"].set_index("stratum")
sp = samp[samp.step == "split"].set_index("stratum")
strata = [1, 2, 3, 4]
names = ["< 100", "100-300", "300-1000", ">= 1000"]
x = np.arange(4)
w = 0.38
elig = el.loc[strata, "eligible"].to_numpy()
tune = sp.loc[strata, "tune"].to_numpy()
conf = sp.loc[strata, "confirm"].to_numpy()
ax.bar(x - w / 2 - 0.01, elig, w, color=GREY, label="eligible")
ax.bar(x + w / 2 + 0.01, conf, w, color=BLUE, label="sampled: confirmation (80%)")
ax.bar(x + w / 2 + 0.01, tune, w, bottom=conf + 15, color=ORANGE, label="sampled: tuning (20%)")
for i in range(4):
    ax.text(x[i] - w / 2, elig[i] + 60, f"{int(elig[i]):,}", ha="center", fontsize=8, color=INK2)
    ax.text(x[i] + w / 2, conf[i] + tune[i] + 75, f"{int(conf[i] + tune[i]):,}", ha="center",
            fontsize=8, color=INK2)
ax.set_xticks(x)
ax.set_xticklabels(names)
ax.set_xlabel("age of the family's sampled clade (million years ago)")
ax.set_ylabel("gene families (HOGs)")
ax.set_title("d  Families per age stratum", loc="left", fontweight="bold", color=INK)
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.grid(axis="y", color=GRID, lw=0.6)
ax.set_axisbelow(True)

fig.tight_layout()
fig.savefig(OUT / "fig2_panel.png", dpi=200, facecolor="white")
print("wrote", OUT / "fig2_panel.png")
