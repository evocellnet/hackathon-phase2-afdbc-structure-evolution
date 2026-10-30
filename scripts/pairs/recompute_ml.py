#!/usr/bin/env python3
"""Recompute ML distances from stored counts (no realignment), e.g. after §Tuning fixes alpha.

  python scripts/pairs/recompute_ml.py --out data/interim/pairs --aln tdiaa --alpha-aa A --alpha-tdi B \
      [--s0 config/tuning.yaml] [--se-model config/tuning.yaml] [--out-dir DIR] [--recompute-selection]

Rewrites aa_*/tdi_* (core + split halves) at the given alpha. The selection variables aa_ml_seq / aa_pid_seq are
KEPT as stored (fixed alpha = 1.0, C26) unless --recompute-selection is given (then recomputed from the seqonly
aa_all counts at --alpha-selection, default 1.0).
--s0 fills lddt_D / lddt_sat; accepted layouts: config/tuning.yaml (s0.{aln}), results/tuning/floors_tuned.json
(s0.{aln}), results/floors/floors.json (C28 table built: S0_ok, fallback to S0 for bins with < 200 ok pairs),
or a legacy JSON with 'S0_lddt_global'. --se-model fills lddt_D_se (tuning.yaml or se_calibration.json,
lddt_se_model.{aln}); requires --s0.
Output: in place (default; the .done sentinel's parquet_md5 / alpha are updated) or --out-dir DIR, which writes
DIR/{aln}/shard_k.parquet + a copy of counts_shard_k.npz + an updated .done (never extra files in the source dir).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

from strucevo.pairs.ml import lddt_D, ml_from_counts  # noqa: E402
from strucevo.pairs.pipeline import SCHEMA_ORDER, md5_file  # noqa: E402

C28_MIN_N = 200


def _load_any(path):
    if path.endswith((".yaml", ".yml")):
        import yaml
        return yaml.safe_load(open(path))
    return json.load(open(path))


def read_s0(path: str, aln: str) -> dict:
    """S0 table {bin_lo, median, q95[, source]} for `aln` from any supported layout (see module doc)."""
    d = _load_any(path)
    if isinstance(d.get("s0"), dict) and isinstance(d["s0"].get(aln), dict) and d["s0"][aln].get("bin_lo") is not None:
        return d["s0"][aln]
    if "S0_lddt_global" in d:
        return d["S0_lddt_global"]
    if aln in d and "S0_ok" in d[aln]:      # results/floors/floors.json -> C28 table
        f = d[aln]
        lo, med, q95, src = [], [], [], []
        for (b0, _), name in zip(d["bin_edges_len_a"], d["bins"]):
            ok, un = f["S0_ok"][name]["LDDT_global"], f["S0"][name]["LDDT_global"]
            use, s = (ok, "S0_ok") if ok["n"] >= C28_MIN_N else (un, "S0")
            lo.append(b0); med.append(use["median"]); q95.append(use["q95"]); src.append(s)
        return {"bin_lo": lo, "median": med, "q95": q95, "source": src}
    raise KeyError(f"{path}: no S0 table for {aln}")


def read_se_model(path: str, aln: str) -> dict:
    d = _load_any(path)
    m = (d.get("lddt_se_model") or {}).get(aln)
    if not m or m.get("sd_S") is None:
        raise KeyError(f"{path}: no lddt_se_model for {aln}")
    return m


def recompute_shard(pq: str, aln: str, alpha: dict, dst_dir: str | None = None, s0=None, se_model=None,
                    recompute_selection: bool = False, alpha_selection: float = 1.0, device=None,
                    boot_sub: set | None = None) -> dict:
    """Recompute one shard. Writes in place (dst_dir None) or to dst_dir/{aln}/. Returns a small summary."""
    src_dir = os.path.dirname(pq)
    k = os.path.basename(pq)[6:10]
    df = pl.read_parquet(pq)
    npz = os.path.join(src_dir, f"counts_shard_{k}.npz")
    z = np.load(npz)
    assert list(z["pair_id"]) == df["pair_id"].to_list(), f"pair_id order mismatch in shard {k}"
    halves = df["pair_set"].str.starts_with("e6_").to_numpy()
    n = df.height
    cols = {}
    for key in ("aa", "tdi"):
        r = ml_from_counts(z[key], key, alpha[key], device)
        cols.update({f"{key}_p": r["p"], f"{key}_ml": r["ml"], f"{key}_ml_se": r["se"], f"{key}_sat": r["sat"]})
        for h in ("hA", "hB"):
            rh = ml_from_counts(z[f"{key}_{h}"][halves], key, alpha[key], device)
            v = np.full(n, np.nan); s = np.full(n, np.nan)
            v[halves] = rh["ml"]; s[halves] = rh["se"]
            cols[f"{key}_ml_{h}"] = v; cols[f"{key}_ml_{h}_se"] = s
    if recompute_selection:   # C26: selection variable at fixed alpha (default 1.0), never the tuned alpha
        zs = z if aln == "seqonly" else np.load(os.path.join(os.path.dirname(src_dir), "seqonly",
                                                             f"counts_shard_{k}.npz"))
        pos = {p: i for i, p in enumerate(zs["pair_id"])}
        ii = np.array([pos[p] for p in df["pair_id"].to_list()], dtype=np.int64)
        r = ml_from_counts(zs["aa_all"][ii], "aa", alpha_selection, device)
        cols["aa_ml_seq"] = r["ml"]; cols["aa_pid_seq"] = 1.0 - r["p"]
    if s0 is not None:
        cols.update(lddt_D(df["lddt_global"].to_numpy(), df["len_a"].to_numpy(), s0, se_model))
    if boot_sub is not None:   # §1.2 boot_sub flag = pair is in the §2.3 column-bootstrap subsample
        cols["boot_sub"] = np.array([p in boot_sub for p in df["pair_id"].to_list()], dtype=bool)
    df = df.with_columns([pl.Series(c, v) for c, v in cols.items()]).select(SCHEMA_ORDER)
    if dst_dir is None:
        out_dir = src_dir
    else:
        out_dir = os.path.join(dst_dir, aln)
        if os.path.realpath(out_dir) == os.path.realpath(src_dir):
            raise ValueError("--out-dir resolves to the source directory; omit it for in-place")
        os.makedirs(out_dir, exist_ok=True)
        shutil.copyfile(npz, os.path.join(out_dir, f"counts_shard_{k}.npz.tmp"))
        os.replace(os.path.join(out_dir, f"counts_shard_{k}.npz.tmp"), os.path.join(out_dir, f"counts_shard_{k}.npz"))
    outp = os.path.join(out_dir, f"shard_{k}.parquet")
    df.write_parquet(outp + ".tmp")
    os.replace(outp + ".tmp", outp)
    done_src = os.path.join(src_dir, f"shard_{k}.done")
    info = json.load(open(done_src)) if os.path.exists(done_src) else {"shard": int(k), "aln": aln, "rows_in": n}
    info["parquet_md5"] = md5_file(outp)
    info["counts_md5"] = md5_file(os.path.join(out_dir, f"counts_shard_{k}.npz"))
    info["rows_out"] = n
    info["alpha"] = alpha
    info["recomputed"] = {"alpha": alpha, "selection": "recomputed" if recompute_selection else "kept (C26, alpha=1)",
                          "alpha_selection": alpha_selection if recompute_selection else 1.0,
                          "s0": s0 is not None, "lddt_D_se": se_model is not None,
                          "boot_sub": boot_sub is not None}
    with open(os.path.join(out_dir, f"shard_{k}.done.tmp"), "w") as fh:
        json.dump(info, fh, indent=1)
    os.replace(os.path.join(out_dir, f"shard_{k}.done.tmp"), os.path.join(out_dir, f"shard_{k}.done"))
    return {"shard": k, "rows": n}


def shard_files(out: str, aln: str) -> list[str]:
    return sorted(glob.glob(os.path.join(out, aln, "shard_[0-9][0-9][0-9][0-9].parquet")))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="pairs dir containing {aln}/shard_k.parquet (source)")
    ap.add_argument("--aln", required=True)
    ap.add_argument("--alpha-aa", type=float, required=True)
    ap.add_argument("--alpha-tdi", type=float, required=True)
    ap.add_argument("--s0", help="tuning.yaml / floors_tuned.json / floors.json (C28) with the S0 table")
    ap.add_argument("--se-model", help="tuning.yaml / se_calibration.json with lddt_se_model (fills lddt_D_se)")
    ap.add_argument("--out-dir", help="write DIR/{aln}/ instead of in place (variants)")
    ap.add_argument("--recompute-selection", action="store_true",
                    help="recompute aa_ml_seq/aa_pid_seq (default: keep as stored, alpha = 1, C26)")
    ap.add_argument("--alpha-selection", type=float, default=1.0)
    a = ap.parse_args()
    alpha = {"aa": a.alpha_aa, "tdi": a.alpha_tdi}
    s0 = read_s0(a.s0, a.aln) if a.s0 else None
    if a.se_model and s0 is None:
        sys.exit("--se-model needs --s0")
    se = read_se_model(a.se_model, a.aln) if a.se_model else None
    for pq in shard_files(a.out, a.aln):
        r = recompute_shard(pq, a.aln, alpha, a.out_dir, s0, se, a.recompute_selection, a.alpha_selection)
        print(f"shard {r['shard']}: {r['rows']} rows recomputed (alpha {alpha})", flush=True)


if __name__ == "__main__":
    main()
