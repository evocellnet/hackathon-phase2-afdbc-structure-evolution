#!/usr/bin/env python3
"""Goal 2-ii step 1 (run in conda_envs/ancestra). Parse the euk180 Ancestra OrthoXML once and write
  data/interim/goal2/hog_members.tsv : acc, species, hog, hog_kind, hog_level, hog_size
  data/interim/goal2/hog_pairs_sample.tsv : within-HOG labelled pairs (hog_pair_labels.label_pairs semantics)
     for HOGs of size 2..MAXSIZE; per HOG at most CAP pairs (Bernoulli thinning, seed 20261030 + crc32(hog)).
Usage: python g2_01_hog_extract.py [--cap 300] [--maxsize 600]
"""
import argparse, os, sys, zlib, random, pickle
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts/ancestra"))
from ancestra.orthoxml import Group, read_orthoxml
from ancestra.newick import TreeIndex, read_newick
from hog_pair_labels import label_pairs, parent_map, COLS

ap = argparse.ArgumentParser()
ap.add_argument("--cap", type=int, default=300)
ap.add_argument("--maxsize", type=int, default=600)
a = ap.parse_args()
X = ROOT / "data/interim/ancestra_euk180/ancestra/hogs.orthoxml"
OUT = ROOT / "data/interim/goal2"; OUT.mkdir(parents=True, exist_ok=True)
ox = read_orthoxml(str(X))
print("genes", len(ox.genes), "top groups", len(ox.groups), flush=True)
sp = {k: v.species for k, v in ox.genes.items()}
n = 0
with open(OUT / "hog_members.tsv", "w") as fh:
    fh.write("acc\tspecies\thog\thog_kind\thog_level\thog_size\n")
    for top in ox.groups:
        if not isinstance(top, Group):
            continue
        g = list(top.iter_genes())
        n += 1
        for x in g:
            fh.write(f"{x}\t{sp[x]}\t{top.id}\t{top.kind}\t{top.level}\t{len(g)}\n")
print("top-level HOGs written", n, flush=True)

tree = ox.species_tree if ox.species_tree is not None else \
    read_newick(str(ROOT / "data/interim/ancestra_euk180/guide_tree.nwk"))
tix = TreeIndex(tree); par = parent_map(tree)
# thinning probabilities per HOG
prob, rng = {}, {}
keep = []
for top in ox.groups:
    if not isinstance(top, Group):
        continue
    s = sum(1 for _ in top.iter_genes())
    if 2 <= s <= a.maxsize:
        npairs = s * (s - 1) / 2
        prob[top.id] = min(1.0, a.cap / npairs)
        rng[top.id] = random.Random(20261030 + zlib.crc32(top.id.encode()))
        keep.append(top)
ox.groups = keep
fh = open(OUT / "hog_pairs_sample.tsv", "w"); fh.write("\t".join(COLS) + "\n")
def emit(*r):
    h = r[-1]
    if rng[h].random() < prob[h]:
        fh.write("\t".join(map(str, r)) + "\n")
st = label_pairs(ox, tix, par, between_species_only=False, emit=emit)
fh.close()
print(dict(st), flush=True)
