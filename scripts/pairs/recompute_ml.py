#!/usr/bin/env python3
"""Recompute ML distances from stored counts (no realignment), e.g. after §Tuning fixes alpha.

  python scripts/pairs/recompute_ml.py --out data/interim/pairs --aln tmalign --alpha-aa A --alpha-tdi B \
      [--s0 results/floors/floors.json] [--inplace | --suffix _a2]

Rewrites aa_*/tdi_* (core + split halves) and, for seqonly, aa_ml_seq/aa_pid_seq; in the other tables aa_ml_seq
is re-joined from the (recomputed) seqonly shards, so run seqonly first. With --s0, fills lddt_D/lddt_sat
from the S0 table (floors.json key 'S0_lddt_global': {bin_lo, median, q95}).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

from strucevo.pairs.ml import lddt_D, ml_from_counts  # noqa: E402
from strucevo.pairs.pipeline import SCHEMA_ORDER, md5_file  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--alpha-aa", type=float, required=True)
    ap.add_argument("--alpha-tdi", type=float, required=True)
    ap.add_argument("--s0")
    ap.add_argument("--suffix", default="", help="write shard_XXXX{suffix}.parquet instead of overwriting")
    a = ap.parse_args()
    alpha = {"aa": a.alpha_aa, "tdi": a.alpha_tdi}
    s0 = json.load(open(a.s0))["S0_lddt_global"] if a.s0 else None
    for pq in sorted(glob.glob(os.path.join(a.out, a.aln, "shard_[0-9][0-9][0-9][0-9].parquet"))):
        k = os.path.basename(pq)[6:10]
        df = pl.read_parquet(pq)
        z = np.load(os.path.join(a.out, a.aln, f"counts_shard_{k}.npz"))
        assert list(z["pair_id"]) == df["pair_id"].to_list(), f"pair_id order mismatch in shard {k}"
        halves = df["pair_set"].str.starts_with("e6_").to_numpy()
        n = df.height
        cols = {}
        for key in ("aa", "tdi"):
            r = ml_from_counts(z[key], key, alpha[key])
            cols.update({f"{key}_p": r["p"], f"{key}_ml": r["ml"], f"{key}_ml_se": r["se"], f"{key}_sat": r["sat"]})
            for h in ("hA", "hB"):
                rh = ml_from_counts(z[f"{key}_{h}"][halves], key, alpha[key])
                v = np.full(n, np.nan); s = np.full(n, np.nan)
                v[halves] = rh["ml"]; s[halves] = rh["se"]
                cols[f"{key}_ml_{h}"] = v; cols[f"{key}_ml_{h}_se"] = s
        if a.aln == "seqonly":
            r = ml_from_counts(z["aa_all"], "aa", alpha["aa"])
            cols["aa_ml_seq"] = r["ml"]; cols["aa_pid_seq"] = 1.0 - r["p"]
        if s0 is not None:
            cols.update(lddt_D(df["lddt_global"].to_numpy(), df["len_a"].to_numpy(), s0))
        df = df.with_columns([pl.Series(c, v) for c, v in cols.items()])
        if a.aln != "seqonly":
            sp = os.path.join(a.out, "seqonly", f"shard_{k}{a.suffix}.parquet")
            sq = pl.read_parquet(sp, columns=["pair_id", "aa_ml_seq", "aa_pid_seq"])
            df = df.drop("aa_ml_seq", "aa_pid_seq").join(sq, on="pair_id", how="left")
        df = df.select(SCHEMA_ORDER)
        outp = pq.replace(".parquet", f"{a.suffix}.parquet")
        df.write_parquet(outp + ".tmp"); os.replace(outp + ".tmp", outp)
        done = pq.replace(".parquet", ".done")
        if not a.suffix and os.path.exists(done):
            info = json.load(open(done)); info["parquet_md5"] = md5_file(outp); info["alpha"] = alpha
            json.dump(info, open(done, "w"), indent=1)
        print(f"shard {k}: {n} rows recomputed (alpha {alpha})", flush=True)


if __name__ == "__main__":
    main()
