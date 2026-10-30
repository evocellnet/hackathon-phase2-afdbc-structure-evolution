#!/usr/bin/env python3
"""Goal 2-ii item 4 (part A): TUNE-split pairs only (pair tables are scanned with a split=='tune' predicate; confirm
rows are never materialised). For 1:1 ortholog primary pairs: is the TED CATH-architecture string different, and are
the two proteins in different AFDB v3 clusters, vs TimeTree age of the species pair (panel-tree MRCA median) and vs
AA distance (aa_ml_seq, seqonly). E6 tune pairs give paralog (e6_para) vs ortholog (e6_orth) at matched AA distance.
Output: data/interim/goal2/tune_pairs_arch.parquet (contains ages: interim only).
"""
from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from g2_common import OUT_I, ROOT, PanelAges, guide_tree  # noqa: E402

cols = ["pair_set", "hog", "split", "gene_a", "gene_b", "species_a", "species_b", "relation", "one_to_one",
        "aa_ml_seq", "aa_pid_seq", "aa_ml", "status"]
t = (pl.scan_parquet(str(ROOT / "data/interim/pairs/seqonly/shard_*.parquet"))
     .filter(pl.col("split") == "tune").select(cols).collect())
assert (t["split"] == "tune").all()
print(t.group_by("pair_set", "relation").len())
_, leaves, _ = guide_tree()
tx2leaf = {}
for l in leaves:
    tx2leaf.setdefault(int(l.rsplit("_", 1)[1]), []).append(l)
assert all(len(v) == 1 for v in tx2leaf.values())
pa = PanelAges()
sp = t.select("species_a", "species_b").unique()
ages = [pa.mrca(frozenset([tx2leaf[a][0], tx2leaf[b][0]]))[1] for a, b in sp.iter_rows()]
sp = sp.with_columns(age_mya=pl.Series(ages))
ann = pl.read_parquet(ROOT / "data/processed/annotations/protein_annot.parquet",
                      columns=["acc", "ted_arch_cath", "ted_n_domains", "ted_n_cath", "afdb_cluster_rep",
                               "afdb_in_cluster", "afdb_cluflag", "len_v4", "mean_plddt"])


def side(s):
    return ann.rename({c: f"{c}_{s}" for c in ann.columns if c != "acc"}).rename({"acc": f"gene_{s}"})


t = (t.join(sp, on=["species_a", "species_b"], how="left").join(side("a"), on="gene_a", how="left")
     .join(side("b"), on="gene_b", how="left"))
t = t.with_columns(
    both_labelled=(pl.col("ted_n_cath_a") > 0) & (pl.col("ted_n_cath_b") > 0),
    arch_diff=pl.col("ted_arch_cath_a").fill_null("") != pl.col("ted_arch_cath_b").fill_null(""),
    # topology-level (C.A.T) architecture: removes H-vs-T label-depth differences of the same domain
    archT_diff=pl.col("ted_arch_cath_a").fill_null("").str.split("|").list.eval(
        pl.element().str.split(".").list.head(3).list.join(".")).list.join("|")
        != pl.col("ted_arch_cath_b").fill_null("").str.split("|").list.eval(
        pl.element().str.split(".").list.head(3).list.join(".")).list.join("|"),
    ndom_diff=(pl.col("ted_n_domains_a").cast(pl.Int64) - pl.col("ted_n_domains_b").cast(pl.Int64)).abs(),
    both_in_cluster=pl.col("afdb_in_cluster_a") & pl.col("afdb_in_cluster_b")
                    & (pl.col("afdb_cluflag_a") != 3) & (pl.col("afdb_cluflag_b") != 3),
    clu_diff=pl.col("afdb_cluster_rep_a") != pl.col("afdb_cluster_rep_b"),
)
t.write_parquet(OUT_I / "tune_pairs_arch.parquet")
print(t.select("age_mya", "aa_ml_seq", "arch_diff", "both_labelled", "clu_diff").describe())
