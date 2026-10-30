#!/usr/bin/env python3
"""Resolution-aware comparison of duplication timing (Ancestra vs Ensembl Compara) for human paralog pairs.

Ensembl and our guide tree sample different species, so the same duplication gets different node names.
Both are converted to an interval on the human NCBI lineage (depth = index in the root->human lineage):
  Ensembl level L  -> duplication in (depth(L_prev), depth(L)], L_prev = next older level among all Ensembl
                      taxonomy levels observed on the human lineage in the REST results;
  Ancestra level P -> duplication in (depth(parent(P)), depth(P)] (P = node below duplication branch).
The two are *compatible* if the intervals overlap.
Usage: python scripts/ancestra/ensembl_level_compat.py RUN_DIR results/ancestra/ensembl_check
"""
import os, re, sys
import pandas as pd
from ancestra.newick import read_newick
from ancestra.speciestree import load_taxdump

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
run, d = sys.argv[1:3]
df = pd.read_csv(os.path.join(d, "ensembl_vs_ancestra_pairs.tsv"), sep="\t")
parent, rank, name, merged = load_taxdump(os.path.join(ROOT, "data", "interim", "ancestra", "taxdump"))
lin = [9606]
while lin[-1] != 1:
    lin.append(parent[lin[-1]])
lin = lin[::-1]
depth_of_tid = {t: i for i, t in enumerate(lin)}
name2depth = {name.get(t, str(t)): i for i, t in enumerate(lin)}
name2depth["Homo sapiens"] = depth_of_tid[9606]
tree = read_newick(os.path.join(run, "guide_tree.nwk"))
nodes = {n.name: n for n in tree.preorder()}
def gdepth(nm):
    if nm == "UP000005640_9606":
        return depth_of_tid[9606]
    return name2depth.get(re.sub(r"_\d+$", "", nm))
hum = df[(df.species == "homo_sapiens")]
ens_levels = sorted({l for l in df.ens_level.unique() if l in name2depth}, key=lambda l: name2depth[l])
prev = {l: (name2depth[ens_levels[i - 1]] if i else -1) for i, l in enumerate(ens_levels)}
rows = []
for r in hum[hum.anc_relation == "in_paralog"].itertuples():
    if r.ens_level not in name2depth or r.anc_level not in nodes:
        continue
    e1, e2 = prev[r.ens_level], name2depth[r.ens_level]
    p = nodes[r.anc_level]
    a2 = gdepth(r.anc_level)
    a1 = gdepth(p.parent.name) if p.parent is not None else -1
    if a2 is None or a1 is None:
        continue
    compat = (e1 < a2) and (a1 < e2)
    rel = "compatible" if compat else ("ancestra_older" if a2 <= e1 else "ancestra_younger")
    rows.append((r.ens_level, r.anc_level, rel))
res = pd.DataFrame(rows, columns=["ens_level", "anc_level", "verdict"])
out = [f"Ensembl levels observed on human lineage: {ens_levels}",
       f"human-human pairs paralogous in both: {len(res)}",
       res.verdict.value_counts(normalize=True).round(3).to_string(),
       "", pd.crosstab(res.ens_level, res.verdict).loc[[l for l in ens_levels if l in set(res.ens_level)]].to_string()]
txt = "\n".join(out)
open(os.path.join(d, "ensembl_level_compat.txt"), "w").write(txt + "\n")
print(txt)
