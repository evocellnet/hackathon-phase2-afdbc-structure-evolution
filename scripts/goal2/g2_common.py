"""Shared helpers for Goal 2-ii (cluster audit): NCBI taxonomy, Ancestra guide-tree node leaf sets,
panel-tree (TimeTree) node ages. Ages are read from data/interim/timetree only and never written under results/."""
from __future__ import annotations

from pathlib import Path

import polars as pl
from ete3 import Tree

ROOT = Path(__file__).resolve().parents[2]
SEED = 20261030
TAXDUMP = ROOT / "data/interim/ancestra/taxdump"
GUIDE = ROOT / "data/interim/ancestra_euk180/guide_tree.nwk"
NODE_AGES = ROOT / "data/interim/timetree/node_ages.tsv"
OUT_I = ROOT / "data/interim/goal2"
OUT_R = ROOT / "results/goal2"
OUT_F = ROOT / "figures/goal2"
EUK = 2759


def load_taxonomy():
    par, rank = {}, {}
    with open(TAXDUMP / "nodes.dmp") as fh:
        for line in fh:
            f = line.split("\t|\t")
            t = int(f[0]); par[t] = int(f[1]); rank[t] = f[2]
    with open(TAXDUMP / "merged.dmp") as fh:
        for line in fh:
            f = line.split("\t|\t")
            o, n = int(f[0]), int(f[1].split("\t")[0])
            par.setdefault(o, par.get(n, 1))
    names = {}
    with open(TAXDUMP / "names.dmp") as fh:
        for line in fh:
            f = line.split("\t|\t")
            if f[3].startswith("scientific name"):
                names[int(f[0])] = f[1]
    return par, rank, names


def lineage(t, par):
    out = [t]
    while t != 1 and t in par:
        t = par[t]
        out.append(t)
    return out  # leaf -> root


def guide_tree():
    """Return (leafsets: node name -> frozenset(leaf names), leaves list, tree)."""
    t = Tree(open(GUIDE).read().strip(), format=1, quoted_node_names=True)
    leafsets = {}
    for n in t.traverse():
        nm = n.name.strip("'")
        leafsets[nm] = frozenset(x.name for x in n.get_leaves())
    return leafsets, [x.name for x in t.get_leaves()], t


def leaf_taxid(leaf: str) -> int:
    return int(leaf.rsplit("_", 1)[1])


class PanelAges:
    """MRCA age of a set of panel leaves on the dated panel tree (node_ages.tsv)."""

    def __init__(self):
        na = pl.read_csv(NODE_AGES, separator="\t", infer_schema_length=0)
        self.nodes = []
        for r in na.iter_rows(named=True):
            if r["is_leaf"] == "True":
                continue
            lv = frozenset(r["leaves"].split(","))
            age = float(r["median"]) if r["median"] not in (None, "") else float("nan")
            self.nodes.append((len(lv), lv, r["node"], age))
        self.nodes.sort(key=lambda x: x[0])
        self.cache = {}

    def mrca(self, leaves: frozenset):
        if leaves in self.cache:
            return self.cache[leaves]
        if len(leaves) == 1:
            res = (next(iter(leaves)), 0.0)
        else:
            res = (None, float("nan"))
            for n, lv, name, age in self.nodes:
                if n >= len(leaves) and leaves <= lv:
                    res = (name, age)
                    break
        self.cache[leaves] = res
        return res
