#!/usr/bin/env python3
"""Fetch the AFDB v4 FoldComp entries listed in a TSV (columns acc, fc_key, entry) by HTTP range requests
into a local FoldComp DB. Thin wrapper around scripts/panel/fetch_foldcomp.py (same resumable shards/merge).

  python scripts/pairs/fetch_structures.py data/interim/pairs_demo/accessions.tsv data/interim/pairs_demo/fc
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import sys
import time

import polars as pl

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "scripts/panel"))
from fetch_foldcomp import INDEX, merge, worker  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("accessions")
    ap.add_argument("out", help="output FoldComp DB path (shards go to <out>_shards)")
    ap.add_argument("--procs", type=int, default=8)
    ap.add_argument("--threads", type=int, default=64)
    a = ap.parse_args()
    acc = pl.read_csv(a.accessions, separator="\t", schema_overrides={"fc_key": pl.Int64})
    idx = (pl.scan_csv(INDEX, separator="\t", has_header=False, new_columns=["fc_key", "off", "len"],
                       schema_overrides={"fc_key": pl.Int64, "off": pl.Int64, "len": pl.Int64})
           .join(acc.lazy(), on="fc_key", how="inner").collect().sort("off"))
    assert idx.height == acc.height, (idx.height, acc.height)
    rows = list(idx.select("fc_key", "off", "len", "entry").iter_rows())
    sdir = a.out + "_shards"
    os.makedirs(sdir, exist_ok=True)
    shards = [[] for _ in range(a.procs)]
    for r in rows:
        shards[r[0] % a.procs].append(r)
    t = time.time()
    with mp.Pool(a.procs) as pool:
        for wid, ok, fail in pool.imap_unordered(worker, [(i, s, sdir, a.threads) for i, s in enumerate(shards)]):
            print(f"[w{wid}] fetched {ok}, failed {fail}", flush=True)
    print(f"fetch {len(rows)} entries ({idx['len'].sum()/1e9:.2f} GB) in {time.time()-t:.0f}s", flush=True)
    merge(sdir, os.path.abspath(a.out))


if __name__ == "__main__":
    main()
