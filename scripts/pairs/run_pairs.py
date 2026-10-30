#!/usr/bin/env python3
"""Per-pair divergence table (ANALYSIS_PLAN §1.2) for alignments seqonly (first), tmalign (primary), tdiaa.

  python scripts/pairs/run_pairs.py --pairs pairs.tsv --foldcomp <FoldComp DB> --store <store dir> \
      --out data/interim/pairs [--alns seqonly tmalign tdiaa] [--n-shards 200] [--alpha-aa 1.0 --alpha-tdi 1.0]

pairs.tsv: tab-separated with header; required gene_a, gene_b (UniProt accessions); optional metadata columns
pair_set, hog, split, species_a, species_b, relation, one_to_one, dcg_id, event_level, dup_branch_from.
The structure store is built once from the FoldComp DB (skipped if <store>/BUILD.json exists).
Shards by md5(hog) (pair_id for hog = NA); reruns skip shards whose .done sentinel matches. Rows out == rows in
is asserted per shard and overall; pairs without structures or without a Foldseek row get status fail_align.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))
# one BLAS thread per forked worker (set before numpy is imported); Foldseek gets --threads explicitly
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import polars as pl  # noqa: E402

from strucevo.pairs.pipeline import ALIGNMENTS, Runner, prepare_pairs  # noqa: E402
from strucevo.pairs.store import Store, build_store  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--foldcomp", help="FoldComp DB of all involved structures (needed only to build the store)")
    ap.add_argument("--store", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--alns", nargs="+", default=list(ALIGNMENTS), choices=ALIGNMENTS)
    ap.add_argument("--n-shards", type=int, default=200)
    ap.add_argument("--shards", type=int, nargs="*", help="run only these shard ids")
    ap.add_argument("--shard-offset", type=int, default=0,
                    help="shard id = offset + md5 %% n_shards; lets separate runs (e.g. tune / confirm split) share --out")
    ap.add_argument("--threads", type=int, default=48)
    ap.add_argument("--alpha-aa", type=float, default=1.0, help="placeholder until §Tuning")
    ap.add_argument("--alpha-tdi", type=float, default=1.0, help="placeholder until §Tuning")
    ap.add_argument("--tmalign-fast", type=int, default=1)
    ap.add_argument("--device", default=None)
    ap.add_argument("--dump-columns", metavar="DIR", default=None,
                    help="also save per-aligned-column states of every pair of this run to DIR/{aln}/cols_shard_k.npz "
                         "(ANALYSIS_PLAN §2.3 column bootstrap; default off -> outputs unchanged). Use with a separate "
                         "--out (e.g. data/interim/pairs_bootsub), never the production pairs dir")
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(a.out, exist_ok=True)
    if not os.path.exists(os.path.join(a.store, "BUILD.json")):
        if not a.foldcomp:
            sys.exit("store missing and no --foldcomp given")
        info = build_store(a.foldcomp, a.store, a.threads)
        print("store built:", json.dumps(info), flush=True)
    store = Store(a.store)
    pairs = pl.read_csv(a.pairs, separator="\t", infer_schema_length=0)
    n_in = pairs.height
    pairs = prepare_pairs(pairs, store, a.n_shards)
    if a.shard_offset:
        pairs = pairs.with_columns(pl.col("shard") + a.shard_offset)
    n_nostruct = int(((pairs["qidx"] < 0) | (pairs["tidx"] < 0)).sum())
    print(f"{n_in} pairs, {n_nostruct} without structure in store, {pairs['shard'].n_unique()} shards", flush=True)
    del store
    with Runner(a.store, a.out, a.threads, a.alpha_aa, a.alpha_tdi, a.device, a.tmalign_fast,
                log=lambda s: print(s, flush=True), dump_dir=a.dump_columns) as R:
        infos = R.run(pairs, a.alns, a.shards)
    summ = {}
    for i in infos:
        s = summ.setdefault(i["aln"], {"rows_in": 0, "rows_out": 0, "status": {}, "timings": {}})
        s["rows_in"] += i["rows_in"]; s["rows_out"] += i["rows_out"]
        for k, v in i["status"].items():
            s["status"][k] = s["status"].get(k, 0) + v
        if not i["skipped"]:
            for k, v in i["timings"].items():
                s["timings"][k] = round(s["timings"].get(k, 0) + v, 2)
    for aln, s in summ.items():
        if a.shards is None:
            assert s["rows_in"] == s["rows_out"] == n_in, (aln, s, n_in)
    summ["wall_s"] = round(time.time() - t0, 1)
    print(json.dumps(summ, indent=1))
    json.dump(summ, open(os.path.join(a.out, f"run_summary_{int(time.time())}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
