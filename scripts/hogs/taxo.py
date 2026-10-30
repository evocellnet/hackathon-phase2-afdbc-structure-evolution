"""NCBI taxonomy helpers (data/raw/ncbi_taxonomy/new_taxdump)."""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TAX = ROOT / "data/raw/ncbi_taxonomy/new_taxdump"


class Taxonomy:
    def __init__(self, taxdir=TAX):
        self.parent, self.rank, self.name = {}, {}, {}
        with open(taxdir / "nodes.dmp") as fh:
            for line in fh:
                f = line.split("\t|\t")
                t = int(f[0]); self.parent[t] = int(f[1]); self.rank[t] = f[2]
        with open(taxdir / "names.dmp") as fh:
            for line in fh:
                f = line.split("\t|\t")
                if f[3].startswith("scientific name"):
                    self.name[int(f[0])] = f[1]
        self.merged = {}
        with open(taxdir / "merged.dmp") as fh:
            for line in fh:
                f = line.split("\t|\t")
                self.merged[int(f[0])] = int(f[1].split("\t")[0])

    def norm(self, t: int) -> int:
        return self.merged.get(t, t)

    def lineage(self, t: int) -> list[int]:
        """[t, parent, ..., 1]"""
        t = self.norm(t); out = [t]
        while t in self.parent and self.parent[t] != t:
            t = self.parent[t]; out.append(t)
        return out
