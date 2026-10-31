#!/usr/bin/env python3
"""Map the literature deep-node ages (agents/deep_ages_report.md; data/interim/timetree/deep_ages_alt.tsv, curated on
the PRE-D30 collapsed panel tree in data/interim/timetree/old_ci_collapse/) onto the CURRENT panel tree
(data/interim/timetree/panel_tree.nwk / node_ages.tsv, D30 collapse rule) for the E4 deep-age sensitivities
(ANALYSIS_PLAN §4.4 "Strassert 2021 / Betts 2018 deep ages").

Mapping is by identical panel-leaf sets (old node -> new node). Two clades that the old tree had collapsed into the
root polytomy and that exist as nodes in the current tree are added from the curator's reference table
(deep_ages_report.md, "Collapsed-out deep nodes"):
  Opisthokonta  = Fungi + Holozoa (+ Fonticula)   S21 n242 1359 (1051-1679); B18 1125-1707   [exact]
  Archaeplastida = Rhodophyta + Viridiplantae     S21 n191 1925 (1799-2063); B18 1118-1667   [approximate: the
                   panel has no glaucophyte; if Glaucophyta branch first in S21 the Rh+Vir MRCA is younger]
Only DEEP nodes (TimeTree median >= 500 My, the curator's own cut: Embryophyta "TT < 500, reference only") are
overridden. Betts 2018 gives intervals only: point = geometric midpoint sqrt(lo*hi) (the median of the plan's
lognormal age model with that 95% interval); lo/hi = the interval.

Output: data/interim/timetree/deep_ages_panel.tsv (ages stay under data/interim/timetree; trap #7).
  node, clade, timetree_median, s21_median, s21_lo, s21_hi, b18_median, b18_lo, b18_hi, mapping, old_node
"""
import argparse
import math
from pathlib import Path

import polars as pl

PROJECT = Path(__file__).resolve().parents[2]
EXTRA = [  # (clade, union of old nodes, s21 (med, lo, hi), b18 (lo, hi), mapping)
    ("Opisthokonta (Fungi + Fonticula + Holozoa)", ("N026", "N044"), (1359.0, 1051.0, 1679.0), (1125.0, 1707.0), "exact"),
    ("Archaeplastida (Rhodophyta + Viridiplantae)", ("N006", "N008"), (1925.0, 1799.0, 2063.0), (1118.0, 1667.0),
     "approximate"),
]
DEEP_MIN = 500.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=str(PROJECT))
    a = ap.parse_args()
    tt = Path(a.root) / "data" / "interim" / "timetree"
    old = pl.read_csv(tt / "old_ci_collapse" / "node_ages.tsv", separator="\t", infer_schema_length=0)
    new = pl.read_csv(tt / "node_ages.tsv", separator="\t", infer_schema_length=0)
    alt = pl.read_csv(tt / "deep_ages_alt.tsv", separator="\t", infer_schema_length=0)
    L = lambda s: frozenset((s or "").split(","))
    old_leaves = {r["node"]: L(r["leaves"]) for r in old.iter_rows(named=True)}
    new_of = {L(r["leaves"]): r["node"] for r in new.iter_rows(named=True) if r["is_leaf"] == "False"}
    tt_med = {r["node"]: r["median"] for r in new.iter_rows(named=True)}
    f = lambda v: float(v) if v not in (None, "", "NA") else None
    rows = []

    def add(node, clade, s21, b18, mapping, old_node):
        tm = f(tt_med.get(node))
        if tm is None or tm < DEEP_MIN:
            print(f"  skip {node} {clade[:40]}: TimeTree median {tm} < {DEEP_MIN} (not a deep node)")
            return
        bmed = math.sqrt(b18[0] * b18[1]) if b18 and None not in b18 else None
        rows.append({"node": node, "clade": clade, "timetree_median": tm,
                     "s21_median": s21[0] if s21 else None, "s21_lo": s21[1] if s21 else None, "s21_hi": s21[2] if s21 else None,
                     "b18_median": bmed, "b18_lo": b18[0] if bmed else None, "b18_hi": b18[1] if bmed else None,
                     "mapping": mapping, "old_node": old_node})

    for r in alt.iter_rows(named=True):
        node = new_of.get(old_leaves[r["panel_node_id"]])
        if node is None:
            raise SystemExit(f"no current node with the leaf set of old {r['panel_node_id']}")
        s21 = (f(r["strassert2021_median"]), f(r["strassert2021_lo"]), f(r["strassert2021_hi"]))
        b18 = (f(r["betts2018_lo"]), f(r["betts2018_hi"]))
        s21 = s21 if s21[0] is not None else None
        b18 = b18 if b18[0] is not None else None
        if s21 is None and b18 is None:
            continue
        mt = r["source_detail"].split("]")[0].strip('"[') if r["source_detail"] else ""
        add(node, r["clade"], s21, b18, mt, r["panel_node_id"])
    for clade, olds, s21, b18, mt in EXTRA:
        ls = frozenset().union(*[old_leaves[o] for o in olds])
        node = new_of.get(ls)
        if node is None:      # Opisthokonta: current node also contains Fonticula (Holomycota) -> MRCA = crown
            cand = [r["node"] for r in new.iter_rows(named=True) if r["is_leaf"] == "False" and ls <= L(r["leaves"])]
            node = min(cand, key=lambda k: len(L(new.filter(pl.col("node") == k)["leaves"][0])))
            extra = L(new.filter(pl.col("node") == node)["leaves"][0]) - ls
            assert len(extra) <= 1, (clade, extra)
        add(node, clade, s21, b18, mt, "+".join(olds))
    df = pl.DataFrame(rows)
    out = tt / "deep_ages_panel.tsv"
    df.write_csv(out, separator="\t", null_value="NA")
    print(df)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
