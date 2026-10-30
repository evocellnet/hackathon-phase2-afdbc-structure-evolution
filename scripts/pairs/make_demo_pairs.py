#!/usr/bin/env python3
"""Demo / benchmark pair list for the pairs pipeline, drawn from the euk20 Ancestra pilot.

Selection (seeded, np.random.SeedSequence(20261030) stage 'pairs_demo' -> spawn index 99, demo only):
  * 45,000 random ortholog pairs (any species pair)            pair_set=primary
  * 5,000 human-mouse 1:1 orthologs (V5 identity check)        pair_set=primary
  * 5,000 random paralog pairs (in_ + out_paralog)             pair_set=e6_para (exercises split halves)
  * 1,000 self pairs (gene vs itself)                          pair_set=v5_self, relation=self
Only proteins with a single-fragment AFDB v4 model (F1, no F2) are used, so every pair has structures.

Writes data/interim/pairs_demo/pairs.tsv and accessions.tsv (acc, fc_key, entry).
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import polars as pl

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--labels", default=os.path.join(ROOT, "results/ancestra/runs/euk20_accurate/pair_labels.tsv.gz"))
    ap.add_argument("--out", default=os.path.join(ROOT, "data/interim/pairs_demo"))
    ap.add_argument("--n-orth", type=int, default=45000)
    ap.add_argument("--n-hm", type=int, default=5000)
    ap.add_argument("--n-para", type=int, default=5000)
    ap.add_argument("--n-self", type=int, default=1000)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(np.random.SeedSequence(20261030).spawn(100)[99])

    lab = pl.read_csv(a.labels, separator="\t", infer_schema_length=0)
    accs = pl.concat([lab["gene_a"], lab["gene_b"]]).unique()
    lk = (pl.scan_parquet(os.path.join(ROOT, "data/interim/panel/afdb_v4_lookup.parquet"))
          .filter(pl.col("acc").is_in(accs.to_list())).collect())
    multi = set(lk.filter(pl.col("frag") > 1)["acc"].to_list())
    lk = lk.filter((pl.col("frag") == 1) & ~pl.col("acc").is_in(list(multi))).unique("acc")
    have = set(lk["acc"].to_list())
    print(f"label accessions {len(accs)}; with single-fragment v4 model {len(have)}")
    lab = lab.filter(pl.col("gene_a").is_in(list(have)) & pl.col("gene_b").is_in(list(have)))

    def take(df, n):
        idx = rng.choice(df.height, size=min(n, df.height), replace=False)
        return df[np.sort(idx)]

    hm = lab.filter((pl.col("one_to_one") == "1") & (
        ((pl.col("species_a") == "UP000005640_9606") & (pl.col("species_b") == "UP000000589_10090")) |
        ((pl.col("species_b") == "UP000005640_9606") & (pl.col("species_a") == "UP000000589_10090"))))
    hm = take(hm, a.n_hm)
    orth = take(lab.filter(pl.col("relation") == "ortholog").join(hm, on=["gene_a", "gene_b"], how="anti"), a.n_orth)
    para = take(lab.filter(pl.col("relation").is_in(["in_paralog", "out_paralog"])), a.n_para)
    pairs = pl.concat([orth.with_columns(pl.lit("primary").alias("pair_set")),
                       hm.with_columns(pl.lit("primary").alias("pair_set")),
                       para.with_columns(pl.lit("e6_para").alias("pair_set"))])
    genes = np.array(sorted(set(pairs["gene_a"].to_list()) | set(pairs["gene_b"].to_list())))
    sp_of = dict(zip(pairs["gene_a"].to_list() + pairs["gene_b"].to_list(),
                     pairs["species_a"].to_list() + pairs["species_b"].to_list()))
    selfg = rng.choice(genes, size=a.n_self, replace=False)
    selfp = pl.DataFrame({"gene_a": selfg, "gene_b": selfg,
                          "species_a": [sp_of[g] for g in selfg], "species_b": [sp_of[g] for g in selfg]}).with_columns(
        pl.lit("self").alias("relation"), pl.lit("v5_self").alias("pair_set"), pl.lit("NA").alias("hog"))
    pairs = pl.concat([pairs, selfp], how="diagonal")
    # taxid = suffix of the proteome id
    for c in ("species_a", "species_b"):
        pairs = pairs.with_columns(pl.col(c).str.split("_").list.last().cast(pl.Int64))
    pairs = pairs.with_columns(pl.lit("none").alias("split"))
    keep = ["gene_a", "gene_b", "pair_set", "hog", "split", "species_a", "species_b", "relation", "one_to_one",
            "dcg_id", "event_level", "dup_branch_from"]
    pairs = pairs.select(keep)
    pairs.write_csv(os.path.join(a.out, "pairs.tsv"), separator="\t")
    used = set(pairs["gene_a"].to_list()) | set(pairs["gene_b"].to_list())
    lk.filter(pl.col("acc").is_in(list(used))).select("acc", "fc_key", "entry").sort("fc_key").write_csv(
        os.path.join(a.out, "accessions.tsv"), separator="\t")
    print(pairs.group_by("pair_set", "relation").len().sort("pair_set"))
    print("unique proteins", len(used))


if __name__ == "__main__":
    main()
