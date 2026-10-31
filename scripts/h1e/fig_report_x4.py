#!/usr/bin/env python3
"""Report figures for the AF-vs-PDB comparison (X4), redrawn from the same paired data as analyse_h1e.py
(no new statistics; same pair construction and mismatch filter).
  figures/h1e/x4_report_main.png  : 3Di ML distance and 1 - global LDDT (report Figure 7)
  figures/h1e/x4_report_all.png   : all five measures (report Figure S1)
Run: python scripts/h1e/fig_report_x4.py
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import polars as pl
import analyse_h1e as A

LAB = {"d_TM": "1 − TM-score", "d_LDDT": "1 − global LDDT", "d_LDDTcore": "1 − core LDDT",
       "d_3Di": "3Di distance (substitutions/site)", "RMSD_core": "core RMSD (Å)"}

TIT = {"d_TM": "1 − TM-score", "d_LDDT": "1 − global LDDT", "d_LDDTcore": "1 − core LDDT",
       "d_3Di": "3Di distance", "RMSD_core": "Core RMSD"}


def build():
    ms = pl.read_csv(A.WD / "match_summary.tsv", separator="\t", infer_schema_length=0).filter(pl.col("status") == "ok")
    ms = ms.with_columns((pl.col("n_mismatch").cast(pl.Float64) / pl.col("n_pdb_bb").cast(pl.Float64)).alias("mis"))
    bad = list(set(ms.filter(pl.col("mis") > A.MISMATCH_MAX)["tag"].to_list()))
    d = A.load("tdiaa")
    tag_ok = lambda df: df.filter(~pl.col("gene_a").str.slice(2).is_in(bad) & ~pl.col("gene_b").str.slice(2).is_in(bad))  # noqa
    daa = d.filter(pl.col("pair_set").str.starts_with("x4m_af_") & (pl.col("status") == "ok")).select(
        pl.col("pair_set").str.replace("x4m_af_", "").alias("choice"), "k1", "k2", pl.col("aa_ml").alias("d_aa"))
    J = {}
    for c in ("post", "pre"):
        af = tag_ok(d.filter(pl.col("pair_set") == f"x4m_af_{c}"))
        pdb = tag_ok(d.filter(pl.col("pair_set") == f"x4m_pdb_{c}"))
        cols = ["status", "aa_ml", "n_core"] + list(A.METRICS)
        j = af.select("k1", "k2", "hog", "split", *[pl.col(x).alias(x + "_af") for x in cols]).join(
            pdb.select("k1", "k2", *[pl.col(x).alias(x + "_pdb") for x in cols]), on=["k1", "k2"])
        j = j.join(daa.filter(pl.col("choice") == c).drop("choice"), on=["k1", "k2"], how="left")
        J[c] = j.filter((pl.col("status_af") == "ok") & (pl.col("status_pdb") == "ok") & pl.col("d_aa").is_not_null())
    return J


def draw(J, metrics, fn, size, fs):
    A.style()
    plt.rcParams.update({"font.size": fs, "axes.labelsize": fs, "xtick.labelsize": fs - 1, "ytick.labelsize": fs - 1})
    fig, axs = plt.subplots(1, len(metrics), figsize=size)
    for ax, m in zip(axs, metrics):
        for c, col, nm in (("pre", A.C_PRE, "Chains released before AlphaFold2's training cutoff"),
                           ("post", A.C_POST, "Chains released after the cutoff, primary set")):
            j = J[c].filter(pl.col(f"{m}_af").is_not_null() & pl.col(f"{m}_pdb").is_not_null())
            ax.scatter(j[f"{m}_pdb"], j[f"{m}_af"], s=14 if len(metrics) < 3 else 9, color=col, alpha=0.65, lw=0.3,
                       ec="white", label=f"{nm}" + (f" (n = {j.height})" if len(metrics) < 3 else ""), zorder=3 if c == "post" else 2)
        lo, hi = ax.get_xlim(); lo2, hi2 = ax.get_ylim(); L = (min(lo, lo2, 0), max(hi, hi2))
        ax.plot(L, L, color="#52514e", lw=0.9, ls="--", zorder=1, label="Equal divergence (AF = crystal)")
        ax.set_xlim(L); ax.set_ylim(L); ax.set_aspect("equal")
        ax.set_xlabel(f"Crystal structures: {LAB[m]}"); ax.set_ylabel(f"AlphaFold models: {LAB[m]}")
        ax.set_title(TIT[m], loc="left")
        ax.grid(color="#e6e5e0", lw=0.4)
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=1 if len(metrics) < 3 else 3, bbox_to_anchor=(0.5, 1.12 if len(metrics) < 3 else 1.08),
               markerscale=1.4)
    fig.tight_layout(); fig.savefig(A.FIG / fn, bbox_inches="tight", dpi=300); plt.close(fig)


if __name__ == "__main__":
    J = build()
    print({c: J[c].height for c in J})
    draw(J, ["d_3Di", "d_LDDT"], "x4_report_main.png", (8.4, 4.4), 10)
    draw(J, list(A.METRICS), "x4_report_all.png", (15, 3.6), 8)
