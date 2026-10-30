#!/usr/bin/env python3
"""V8 (parsers) and V9 (in-house lDDT vs Foldseek lddtfull) on n pairs (ANALYSIS_PLAN §7).

V8  * Foldseek column order: convertalis is re-run with the pipeline's FS_FIELDS plus qaln,taln,qseq,tseq;
      query/target names map to store entries whose length equals qlen/tlen; qseq/tseq equal the store sequence;
      qstart <= qend <= qlen etc.; the number of lddtfull values equals the number of cigar M columns.
    * cigar-derived alignment (strucevo.pairs.align.pairs_to_aln_strings) reproduces qaln/taln exactly.
V9  * strucevo.lddt per-residue vs lddtfull (3-decimal) and global LDDT (sum/len_a): |Δ| <= 0.01 on every pair.
Both alignments tmalign and tdiaa are checked. Pairs: a seeded random sample of rows from a pair table
(parquet glob or TSV with gene_a, gene_b; oriented as in the pipeline).

  python scripts/verify/v8_v9_foldseek.py --store S --pairs 'data/interim/pairs_demo/pairs/tmalign/*.parquet' \
      --n 1000 --json results/verify/V8_V9.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

from strucevo.lddt import per_residue_lddt  # noqa: E402
from strucevo.pairs.align import FS_FIELDS, cigar_to_pairs, pairs_to_aln_strings, run_foldseek  # noqa: E402
from strucevo.pairs.pipeline import prepare_pairs  # noqa: E402
from strucevo.pairs.store import Store  # noqa: E402

EXTRA = ["qaln", "taln", "qseq", "tseq"]


def check(store, P, method, threads, tmpd):
    tsv, _ = run_foldseek(store.db, store.key[P["qidx"].to_numpy()], store.key[P["tidx"].to_numpy()],
                          os.path.join(tmpd, method), method, threads=threads, extra_fields=EXTRA)
    fs = pl.read_csv(tsv, separator="\t", has_header=False, new_columns=FS_FIELDS + EXTRA, quote_char=None,
                     infer_schema_length=0)
    r = {"n_pairs_in": P.height, "n_rows_out": fs.height}
    bad = {"name": 0, "len": 0, "seq": 0, "bounds": 0, "lddtfull_len": 0, "qaln": 0, "taln": 0}
    d_res, d_glob, d_aln = [], [], []
    for row in fs.iter_rows(named=True):
        qi_, ti_ = store.idx_of_name.get(row["query"]), store.idx_of_name.get(row["target"])
        if qi_ is None or ti_ is None:
            bad["name"] += 1
            continue
        qs, ts = store.seq(qi_), store.seq(ti_)
        if int(row["qlen"]) != len(qs) or int(row["tlen"]) != len(ts):
            bad["len"] += 1
        if row["qseq"] != qs or row["tseq"] != ts:
            bad["seq"] += 1
        q0, q1, t0, t1 = (int(row[c]) for c in ("qstart", "qend", "tstart", "tend"))
        if not (1 <= q0 <= q1 <= len(qs) and 1 <= t0 <= t1 <= len(ts)):
            bad["bounds"] += 1
        qa, ta = pairs_to_aln_strings(qs, ts, q0, t0, row["cigar"])
        bad["qaln"] += qa != row["qaln"]
        bad["taln"] += ta != row["taln"]
        qi, ti = cigar_to_pairs(q0, t0, row["cigar"])
        lf = np.array(row["lddtfull"].split(","), float) if row["lddtfull"] else np.zeros(0)
        if len(lf) != len(qi):
            bad["lddtfull_len"] += 1
            continue
        xq = store.get(qi_)[0]; xt = store.get(ti_)[0]
        pr = per_residue_lddt(xq, xt, qi, ti)
        if len(pr):
            d_res.append(np.abs(pr - lf).max())
            d_aln.append(abs(pr.mean() - float(row["lddt"])))
        d_glob.append(abs(pr.sum() / len(qs) - lf.sum() / len(qs)))
    d_res, d_glob, d_aln = map(np.array, (d_res, d_glob, d_aln))
    r["V8_bad_counts"] = bad
    r["V8_pass"] = bool(r["n_rows_out"] == r["n_pairs_in"] and sum(bad.values()) == 0)
    r["V9_global_absdiff"] = {"median": float(np.median(d_glob)), "q99": float(np.quantile(d_glob, 0.99)),
                              "max": float(d_glob.max()), "n": int(len(d_glob))}
    r["V9_per_residue_maxabsdiff"] = {"median": float(np.median(d_res)), "max": float(d_res.max())}
    r["V9_lddt_aln_vs_foldseek_lddt"] = {"median": float(np.median(d_aln)), "max": float(d_aln.max())}
    r["V9_pass"] = bool(d_glob.max() <= 0.01 and len(d_glob) == P.height)
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--store", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=48)
    ap.add_argument("--json", required=True)
    a = ap.parse_args()
    store = Store(a.store)
    P = pl.read_parquet(a.pairs) if a.pairs.endswith(".parquet") else pl.read_csv(a.pairs, separator="\t",
                                                                                 infer_schema_length=0)
    P = P.select("gene_a", "gene_b").unique().sort("gene_a", "gene_b")
    rng = np.random.default_rng(20261030)
    P = P[np.sort(rng.choice(P.height, size=min(a.n, P.height), replace=False))]
    P = prepare_pairs(P, store, 1).filter((pl.col("qidx") >= 0) & (pl.col("tidx") >= 0))
    out = {}
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(a.json))) as tmpd:
        for m in ("tmalign", "tdiaa"):
            out[m] = check(store, P, m, a.threads, tmpd)
            print(m, json.dumps(out[m], indent=1), flush=True)
    out["V8_pass"] = all(out[m]["V8_pass"] for m in ("tmalign", "tdiaa"))
    out["V9_pass"] = all(out[m]["V9_pass"] for m in ("tmalign", "tdiaa"))
    os.makedirs(os.path.dirname(os.path.abspath(a.json)), exist_ok=True)
    json.dump(out, open(a.json, "w"), indent=1)
    print("V8", "PASS" if out["V8_pass"] else "FAIL", "| V9", "PASS" if out["V9_pass"] else "FAIL")


if __name__ == "__main__":
    main()
