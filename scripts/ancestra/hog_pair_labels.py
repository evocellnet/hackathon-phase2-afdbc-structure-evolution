#!/usr/bin/env python3
"""Label every pair of proteins in an Ancestra ``hogs.orthoxml`` as ortholog / paralog, with the
species-tree node of the event that separates them (speciation or duplication).

Semantics (Ancestra docs, ancestra/orthoxml.py; identical to the QfO / pyham flattening):
  * The *deepest common group* (DCG) of two genes decides the relation.
  * DCG = orthologGroup and different species  -> ``ortholog``. The group's TaxRange is the
    speciation node (== MRCA of the two species; checked, column ``lca_ok``).
  * DCG = paralogGroup                          -> ``paralog``. The paralogGroup level P is the
    species-tree node *directly below the branch on which the duplication happened*; the duplication
    is therefore on the branch parent(P) -> P (``dup_branch_from`` -> ``dup_level``). If P is a leaf,
    it is a species-specific (terminal-branch) duplication. ``enclosing_og_level`` is the level of
    the nearest enclosing orthologGroup (the duplication is younger than that speciation).
    Paralogs are further split into ``in_paralog`` (same species) and ``out_paralog`` (different
    species: the duplication predates their speciation).
  * DCG = orthologGroup and same species         -> ``inconsistent`` (QfO drops these; rare).
  * Different top-level HOGs                      -> not listed (unrelated in this file). If the
    run's family table is given (``--families-from-run OUTDIR``), same-family pairs that ended in
    different root HOGs can be emitted as ``root_split`` (= duplication older than the root of the
    guide tree, or a merge split; Ancestra does not distinguish the two).

Extra columns useful for divergence studies:
  ``one_to_one``  ortholog pair whose two sides of the DCG each hold exactly one gene of the
                  pair's species (the classic 1:1 ortholog)
  ``dcg_id``      OrthoXML id of the DCG (pairs from the same duplication node share it)
  ``hog``         top-level HOG id.

Usage:
  python scripts/ancestra/hog_pair_labels.py OUTDIR/hogs.orthoxml pairs_labelled.tsv.gz \
      [--between-species-only] [--families-from-run OUTDIR] [--max-hog-size N]
Requires the ``ancestra`` package (conda_envs/ancestra).
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import os
import sys
from collections import Counter

from ancestra.newick import TreeIndex, read_newick
from ancestra.orthoxml import Group, read_orthoxml

COLS = ["gene_a", "gene_b", "species_a", "species_b", "relation", "event_level", "species_lca",
        "lca_ok", "dup_branch_from", "enclosing_og_level", "one_to_one", "dcg_id", "hog"]


def parent_map(tree):
    par = {}
    for n in tree.preorder():
        for c in n.children:
            par[c.name] = n.name
    return par


def label_pairs(ox, tix, par, between_species_only=False, max_hog_size=None, emit=None):
    sp = {k: v.species for k, v in ox.genes.items()}
    stats = Counter()
    for top in ox.groups:
        if not isinstance(top, Group):
            continue
        hog = top.id or ""
        if max_hog_size and sum(1 for _ in top.iter_genes()) > max_hog_size:
            stats["hogs_skipped_size"] += 1
            continue
        # iterative post-order, carrying the level of the enclosing orthologGroup
        stack = [(top, 0, None)]
        memo = {}
        while stack:
            g, st, enc = stack.pop()
            if st == 0:
                stack.append((g, 1, enc))
                sub_enc = g.level if g.kind == "ortholog" else enc
                stack.extend((c, 0, sub_enc) for c in g.children if isinstance(c, Group))
                continue
            kids = [memo.pop(id(c)) if isinstance(c, Group) else [c] for c in g.children]
            if g.kind == "ortholog":
                # per child: species counts (for 1:1 test)
                cnt = [Counter(sp[x] for x in k) for k in kids]
                for i, j in itertools.combinations(range(len(kids)), 2):
                    for a in kids[i]:
                        sa = sp[a]
                        for b in kids[j]:
                            sb = sp[b]
                            if sa == sb:
                                if between_species_only:
                                    continue
                                rel, ok = "inconsistent", ""
                                lca = sa
                            else:
                                rel = "ortholog"
                                lca = tix.mrca((sa, sb))
                                ok = int(lca == g.level)
                            o2o = int(rel == "ortholog" and cnt[i][sa] == 1 and cnt[j][sb] == 1
                                      and cnt[i][sb] == 0 and cnt[j][sa] == 0)
                            stats[rel] += 1
                            stats["one_to_one"] += o2o
                            emit(a, b, sa, sb, rel, g.level, lca, ok, "", enc or "", o2o, g.id or "", hog)
            else:
                lvl = g.level
                frm = par.get(lvl, "") if lvl else ""
                for i, j in itertools.combinations(range(len(kids)), 2):
                    for a in kids[i]:
                        sa = sp[a]
                        for b in kids[j]:
                            sb = sp[b]
                            if sa == sb:
                                if between_species_only:
                                    continue
                                rel, lca = "in_paralog", sa
                            else:
                                rel, lca = "out_paralog", tix.mrca((sa, sb))
                            stats[rel] += 1
                            emit(a, b, sa, sb, rel, lvl or "", lca, "", frm, enc or "", 0, g.id or "", hog)
            memo[id(g)] = [x for k in kids for x in k]
    return stats


def root_split_pairs(ox, run_dir, tix, emit, between_species_only, max_family_size):
    """Same S0/S2 family but different root HOG (or singleton)."""
    import numpy as np
    from ancestra.io import read_gene_table
    gt = read_gene_table(os.path.join(run_dir, "stages", "input", "genes.tsv"))
    f = os.path.join(run_dir, "stages", "scores", "families_final.npy")
    if not os.path.exists(f):
        f = os.path.join(run_dir, "stages", "families", "families.npy")
    fam = np.load(f)
    hog_of = {}
    for top in ox.groups:
        if isinstance(top, Group):
            for x in top.iter_genes():
                hog_of[x] = top.id
    by_fam = {}
    for gi, fi in enumerate(fam):
        by_fam.setdefault(int(fi), []).append(gt.names[gi])
    sp = {k: v.species for k, v in ox.genes.items()}
    n = 0
    for fi, mem in by_fam.items():
        if len(mem) < 2 or (max_family_size and len(mem) > max_family_size):
            continue
        for a, b in itertools.combinations(mem, 2):
            ha, hb = hog_of.get(a), hog_of.get(b)
            if ha is not None and ha == hb:
                continue
            sa, sb = sp[a], sp[b]
            if between_species_only and sa == sb:
                continue
            emit(a, b, sa, sb, "root_split", "", tix.mrca((sa, sb)), "", "", "", 0, f"fam{fi}",
                 f"{ha or 'singleton'}|{hb or 'singleton'}")
            n += 1
    return n


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("orthoxml")
    ap.add_argument("out", help="TSV output (.gz ok)")
    ap.add_argument("--tree", help="guide tree (default: taxonomy block of the OrthoXML)")
    ap.add_argument("--between-species-only", action="store_true")
    ap.add_argument("--families-from-run", help="Ancestra OUTDIR: also emit same-family root_split pairs")
    ap.add_argument("--max-hog-size", type=int, default=None, help="skip top-level HOGs larger than this")
    a = ap.parse_args(argv)
    ox = read_orthoxml(a.orthoxml)
    tree = read_newick(a.tree) if a.tree else ox.species_tree
    if tree is None:
        sys.exit("no species tree: pass --tree")
    tix, par = TreeIndex(tree), parent_map(tree)
    op = gzip.open if a.out.endswith(".gz") else open
    with op(a.out, "wt") as fh:
        fh.write("\t".join(COLS) + "\n")
        buf = []

        def emit(a, b, sa, sb, *rest):
            if b < a:   # canonical order gene_a < gene_b
                a, b, sa, sb = b, a, sb, sa
            row = (a, b, sa, sb) + rest
            buf.append("\t".join(str(x) for x in row))
            if len(buf) >= 200_000:
                fh.write("\n".join(buf) + "\n")
                buf.clear()
        stats = label_pairs(ox, tix, par, a.between_species_only, a.max_hog_size, emit)
        if a.families_from_run:
            stats["root_split"] = root_split_pairs(ox, a.families_from_run, tix, emit,
                                                   a.between_species_only, a.max_hog_size)
        if buf:
            fh.write("\n".join(buf) + "\n")
    print(dict(stats), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
