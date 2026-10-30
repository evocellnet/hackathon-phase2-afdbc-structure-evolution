#!/usr/bin/env python3
"""Per-split pair lists + structure accession list for run_pairs.py (orchestration step c0).

  python scripts/pipeline/prep_pair_lists.py --pairlists data/interim/pairlists --outdir data/interim/pairruns

Writes <outdir>/pairs_tune.tsv, pairs_confirm.tsv (primary + e6_*; columns as produced by build_hogs) and
<outdir>/accessions.tsv (acc, fc_key, entry: AFDB v4 FoldComp keys of every gene in any list), plus counts.json.
Genes without an AFDB v4 entry are listed in missing_structures.tsv (they become fail_align rows, never dropped).
"""
from __future__ import annotations

import argparse
import json
import os

import polars as pl

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairlists", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--lists", nargs="+", default=["primary.tsv", "e6.tsv"])
    ap.add_argument("--lookup", default=os.path.join(ROOT, "data/interim/panel/afdb_v4_lookup.parquet"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    dfs = [pl.read_csv(os.path.join(a.pairlists, f), separator="\t", infer_schema_length=0) for f in a.lists]
    allp = pl.concat(dfs, how="diagonal")
    counts = {"total": allp.height}
    for sp in ("tune", "confirm"):
        d = allp.filter(pl.col("split") == sp)
        d.write_csv(os.path.join(a.outdir, f"pairs_{sp}.tsv"), separator="\t")
        counts[sp] = {"pairs": d.height, "by_set": dict(d.group_by("pair_set").len().iter_rows()),
                      "hogs": d["hog"].n_unique()}
    other = allp.filter(~pl.col("split").is_in(["tune", "confirm"])).height
    assert other == 0, f"{other} pairs with split not in tune/confirm"
    genes = pl.concat([allp["gene_a"], allp["gene_b"]]).unique()
    lk = (pl.scan_parquet(a.lookup).filter(pl.col("acc").is_in(genes.to_list()))
          .select("acc", "fc_key", "entry").collect().unique("acc").sort("fc_key"))
    lk.write_csv(os.path.join(a.outdir, "accessions.tsv"), separator="\t")
    miss = genes.filter(~genes.is_in(lk["acc"].to_list()))
    pl.DataFrame({"acc": miss}).write_csv(os.path.join(a.outdir, "missing_structures.tsv"), separator="\t")
    counts.update({"genes": genes.len(), "with_structure": lk.height, "missing_structure": miss.len()})
    json.dump(counts, open(os.path.join(a.outdir, "counts.json"), "w"), indent=1)
    print(json.dumps(counts))


if __name__ == "__main__":
    main()
