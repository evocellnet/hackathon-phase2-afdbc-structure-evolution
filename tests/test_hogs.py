"""Tests for the panel tree and HOG construction (scripts/hogs/).

Run:  python -m pytest tests/test_hogs.py -q
Needs: data/interim/timetree/{panel_tree.nwk,node_ages.tsv} (scripts/hogs/panel_tree.py) and the euk20 test build
       data/interim/hogs_test_euk20 (scripts/hogs/build_hogs.py on results/ancestra/runs/euk20_accurate).
"""
from __future__ import annotations

import glob
import itertools
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts/hogs"))
sys.path.insert(0, str(ROOT / "tools/ancestra-mvp-ai"))
sys.path.insert(0, str(ROOT / "scripts/ancestra"))
import seqdist  # noqa: E402
import treeutil as T  # noqa: E402

TT = ROOT / "data/interim/timetree"
TEST = ROOT / "data/interim/hogs_test_euk20"
OXML = ROOT / "results/ancestra/runs/euk20_accurate/ancestra/hogs.orthoxml"


# ------------------------------------------------------------------------------------------- panel tree
@pytest.fixture(scope="module")
def panel():
    tr = T.parse_newick((TT / "panel_tree.nwk").read_text())
    na = pd.read_csv(TT / "node_ages.tsv", sep="\t").set_index("node")
    return tr, na


def test_panel_tree_leaves(panel):
    tr, na = panel
    sp = pd.read_csv(ROOT / "data/processed/panel/species_panel.tsv", sep="\t")
    assert sorted(l.name for l in tr.leaves()) == sorted(sp.leaf)


def test_panel_tree_monotone(panel):
    tr, na = panel
    for n in tr.preorder():
        if n.parent is None or n.is_leaf():
            continue
        a, b = na.at[n.name, "median"], na.at[n.parent.name, "median"]
        if np.isfinite(a) and np.isfinite(b):
            assert a < b, (n.name, a, b)
            assert b - a >= 5.0 - 1e-9, (n.name, a, b)       # collapse rule: no internode < 5 My


def test_human_mouse_age(panel):
    tr, na = panel
    nodes, ls = T.index_tree(tr)
    m = T.mrca(nodes, ls, {"UP000005640_9606", "UP000000589_10090"})
    assert abs(na.at[m.name, "median"] - 87) <= 10


def test_newick_roundtrip():
    s = "((A:1,B:2)'x y':3,(C:1,D:1)E:4)root;"
    r = T.parse_newick(s)
    assert T.to_newick(r) == "((A:1,B:2)'x y':3,(C:1,D:1)E:4)root;"
    p = T.prune_to(r, {"A", "C", "D"})
    assert T.to_newick(p, internal_names=False) == "(A:4,(C:1,D:1):4);"


# ------------------------------------------------------------------------------------------- seqdist
def test_seqdist_matches_pairs_pipeline():
    files = sorted(glob.glob(str(ROOT / "data/interim/pairs_demo/pairs/seqonly/shard_*.parquet")))
    if not files:
        pytest.skip("no pairs demo")
    df = pd.concat([pd.read_parquet(f, columns=["gene_a", "gene_b", "aa_ml_seq", "aa_pid_seq", "pair_set"])
                    for f in files[:2]])
    pt = pd.read_parquet(ROOT / "data/processed/panel/protein_table.parquet", columns=["acc", "seq_match"])
    ok = set(pt.acc[pt.seq_match])
    df = df[(df.pair_set != "v5_self") & df.gene_a.isin(ok) & df.gene_b.isin(ok)].sample(300, random_state=1)
    seqdist.SEQ.update(seqdist.read_fasta_dir(ROOT / "data/interim/ancestra/euk20", set(df.gene_a) | set(df.gene_b)))
    ml, pid, _ = seqdist.seq_distances(list(zip(df.gene_b, df.gene_a)), threads=1)   # reversed on purpose
    assert np.allclose(ml, df.aa_ml_seq, atol=1e-9, equal_nan=True)
    assert np.allclose(pid, df.aa_pid_seq, atol=1e-12, equal_nan=True)


# ------------------------------------------------------------------------------------------- labels
def test_hogindex_equals_hog_pair_labels():
    from ancestra.newick import TreeIndex
    from ancestra.orthoxml import Group, read_orthoxml
    from build_hogs import HogIndex
    from hog_pair_labels import label_pairs, parent_map
    ox = read_orthoxml(str(OXML))
    tops = [g for g in ox.groups if isinstance(g, Group)]
    random.seed(3)
    pick = random.sample([g for g in tops if 10 <= sum(1 for _ in g.iter_genes()) <= 400], 150)
    tix, par = TreeIndex(ox.species_tree), parent_map(ox.species_tree)
    sp_of = {k: v.species for k, v in ox.genes.items()}

    class O:
        pass
    o = O(); o.genes = ox.genes; o.groups = pick
    ref = {}

    def emit(a, b, sa, sb, rel, lvl, lca, ok, frm, enc, o2o, dcg, hog):
        ref[(hog, frozenset((a, b)))] = (rel, lvl or "", frm, int(o2o), dcg)
    label_pairs(o, tix, par, False, None, emit)
    n = 0
    for g in pick:
        hi = HogIndex(g, sp_of, par)
        for a, b in itertools.combinations(list(g.iter_genes()), 2):
            assert hi.label(a, b) == ref[(g.id, frozenset((a, b)))], (g.id, a, b)
            n += 1
    assert n > 10000


# ------------------------------------------------------------------------------------------- build outputs
@pytest.fixture(scope="module")
def build():
    if not (TEST / "hogs.parquet").exists():
        pytest.skip("run the euk20 test build first")
    h = pd.read_parquet(TEST / "hogs.parquet")
    p = pd.read_csv(TEST / "pairlists/primary.tsv", sep="\t")
    e = pd.read_csv(TEST / "pairlists/e6.tsv", sep="\t")
    t = pd.read_csv(TEST / "pairlists/e6_triplets.tsv", sep="\t")
    return h, p, e, t


def test_contract_columns(build):
    h, *_ = build
    need = ["hog", "n_species_1to1", "presample", "selected", "sel_clade", "root_age_mya", "stratum", "in_sample",
            "split", "cath_class", "mean_plddt", "misfit_ratio", "excl_misfit"]
    assert all(c in h.columns for c in need)
    assert h.hog.is_unique
    s = h[h.in_sample]
    assert s.selected.map(len).between(8, 30).all()
    assert (s.presample.map(len) <= 60).all()
    assert s.stratum.isin([1, 2, 3, 4]).all()
    assert set(s.split) <= {"tune", "confirm"} and (h.loc[~h.in_sample, "split"] == "none").all()
    assert (h.n_species_1to1 >= 8).all()
    assert all(set(x) <= set(y) for x, y in zip(s.selected, s.presample))


def test_clade_closed_and_dmax(build):
    h, *_ = build
    tr = T.parse_newick((TT / "panel_tree.nwk").read_text())
    nodes, ls = T.index_tree(tr)
    byname = {n.name: n for n in nodes}
    leaf_of_taxid = {int(l.rsplit("_", 1)[1]): l for l in ls[0]}
    s = h[h.in_sample].sample(min(60, int(h.in_sample.sum())), random_state=0)
    pairs = []
    for r in s.itertuples():
        C = byname[r.sel_clade]
        assert {leaf_of_taxid[t] for t in r.selected} <= ls[C.id]
        pairs += list(itertools.combinations(r.selected_genes, 2))
    seqdist.SEQ.update(seqdist.read_fasta_dir(ROOT / "data/interim/ancestra/euk20", {g for p in pairs for g in p}))
    ml, pid, sat = seqdist.seq_distances(pairs, threads=4)
    assert np.all(np.isfinite(ml)) and not sat.any() and np.all(ml <= 2.0)


def test_primary_pairs(build):
    h, p, e, t = build
    s = h[h.in_sample]
    assert len(p) == int((s.selected.map(len) * (s.selected.map(len) - 1) // 2).sum())
    key = p.apply(lambda r: frozenset((r.gene_a, r.gene_b)), axis=1)
    assert key.is_unique
    assert (p.pair_set == "primary").all() and (p.species_a != p.species_b).all()
    assert set(p.hog) == set(s.hog)


def test_e6_triplets(build):
    h, p, e, t = build
    from ancestra.orthoxml import read_orthoxml
    ox = read_orthoxml(str(OXML))
    leaves = set(ox.species_tree.leaf_names())
    sp_of = {k: v.species for k, v in ox.genes.items()}
    assert not t.event_level.isin(leaves).any()                       # no leaf-level in-paralogs
    for r in t.itertuples():
        assert sp_of[r.P1] == sp_of[r.P2]
        o1 = set(str(r.O1).split(",")); o2 = set(str(r.O2).split(",")) if isinstance(r.O2, str) else set()
        assert len(o1) >= 1
        assert not ({sp_of[x] for x in o1} & {sp_of[x] for x in o2})   # O1 ∩ O2 = ∅
    assert set(e.pair_set) <= {"e6_para", "e6_orth", "e6_ldo"}
    k = e.apply(lambda r: frozenset((r.gene_a, r.gene_b)), axis=1)
    assert k.is_unique
    kp = set(map(frozenset, zip(p.gene_a, p.gene_b)))
    assert not any(x in kp for x in k)
    para = e[e.pair_set == "e6_para"]
    assert (para.relation == "in_paralog").all()
    assert (e[e.pair_set != "e6_para"].relation == "ortholog").all()


def test_split_fraction(build):
    h, *_ = build
    s = h[h.in_sample]
    for k, g in s.groupby("stratum"):
        assert abs((g.split == "tune").sum() - round(0.2 * len(g))) == 0
    sp = pd.read_csv(TEST / "split.tsv", sep="\t")
    assert dict(zip(sp.hog, sp.split)) == dict(zip(s.hog, s.split))
