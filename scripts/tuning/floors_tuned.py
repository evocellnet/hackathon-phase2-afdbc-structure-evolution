#!/usr/bin/env python3
"""ANALYSIS_PLAN §1.4 / §3.2: eta / nu and S0 at the tuned alpha (floors are split-independent).

- S0 per alignment (C28): S0_ok (floor_random pairs with status ok) per length bin, fallback to the
  unconditional S0 for bins with < 200 ok random pairs (results/floors/floors.json, LDDT_global).
- eta_tdi_ml: median tdi_ml of floor_ident pairs (status ok) recomputed from data/interim/pairs_floors/{aln}/
  counts_shard_*.npz at the tuned alpha_3Di (nu = eta/2; eta_aa = 0, measured value reported).
- eta_lddtD: median lddt_D of floor_ident ok pairs with the C28 S0 table (also with unconditional S0: sensitivity).
- eta_lddt_raw: median 1 - lddt_global (alpha-independent; as in floors.json).
Floors dir: ROOT/data/interim/pairs_floors if present, else the project's (production floors are allowed).

  python scripts/tuning/floors_tuned.py [--root R] [--alpha-aa A --alpha-tdi B]
Output ROOT/results/tuning/floors_tuned.json. Top-level keys = config/tuning.yaml nesting:
  eta.{aln}.{tdi_ml, lddt_D, lddt_raw};  s0.{aln}.{bin_lo, median, q95, source}  (+ s0.{aln}.n);  details.
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

from se_common import (ALNS, LEN_BIN_LO, PROJECT, add_root_arg, bin_index, np, out_json, pl, provenance, read_alpha,
                       s0_table, write_json)

from strucevo.pairs.ml import lddt_D, ml_from_counts


def unconditional_s0(root, aln):
    t = s0_table(root, aln, min_n=10 ** 12)   # forces the unconditional table
    return t


def main():
    ap = add_root_arg(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter))
    ap.add_argument("--alns", nargs="+", default=list(ALNS))
    ap.add_argument("--alpha-aa", type=float)
    ap.add_argument("--alpha-tdi", type=float)
    ap.add_argument("--floors-dir")
    a = ap.parse_args()
    alpha, src = read_alpha(a.root, a.alpha_aa, a.alpha_tdi)
    fdir = Path(a.floors_dir) if a.floors_dir else Path(a.root) / "data" / "interim" / "pairs_floors"
    if not fdir.exists():
        fdir = PROJECT / "data" / "interim" / "pairs_floors"
    res = {"step": "floors_tuned", "alpha_used": alpha, "alpha_source": src, "floors_dir": str(fdir),
           "eta": {}, "s0": {}, "details": {}}
    for aln in a.alns:
        dfs, ca, ct = [], [], []
        for f in sorted(glob.glob(str(fdir / aln / "shard_[0-9][0-9][0-9][0-9].parquet"))):
            df = pl.read_parquet(f, columns=["pair_id", "pair_set", "status", "len_a", "lddt_global", "tdi_ml"])
            z = np.load(f.replace("shard_", "counts_shard_").replace(".parquet", ".npz"))
            assert list(z["pair_id"]) == df["pair_id"].to_list()
            m = ((df["pair_set"] == "floor_ident") & (df["status"] == "ok")).to_numpy()
            dfs.append(df.filter(pl.Series(m))); ca.append(z["aa"][m]); ct.append(z["tdi"][m])
        df = pl.concat(dfs)
        ra = ml_from_counts(np.concatenate(ca), "aa", alpha["aa"])
        rt = ml_from_counts(np.concatenate(ct), "tdi", alpha["tdi"])
        r1 = ml_from_counts(np.concatenate(ct), "tdi", 1.0)
        S = df["lddt_global"].to_numpy(); L = df["len_a"].to_numpy()
        s0 = s0_table(a.root, aln)
        s0u = unconditional_s0(a.root, aln)
        D = lddt_D(S, L, s0)["lddt_D"]
        Du = lddt_D(S, L, s0u)["lddt_D"]
        b = bin_index(L, LEN_BIN_LO)
        by_bin = {"bin_lo": list(LEN_BIN_LO),
                  "tdi_ml": [float(np.nanmedian(rt["ml"][b == i])) if (b == i).any() else None for i in range(len(LEN_BIN_LO))],
                  "lddt_D": [float(np.nanmedian(D[b == i])) if (b == i).any() else None for i in range(len(LEN_BIN_LO))],
                  "n": [int((b == i).sum()) for i in range(len(LEN_BIN_LO))]}
        eta = {"tdi_ml": float(np.nanmedian(rt["ml"])), "lddt_D": float(np.nanmedian(D)),
               "lddt_raw": float(np.nanmedian(1 - S))}
        res["eta"][aln] = eta
        res["s0"][aln] = {k: s0[k] for k in ("bin_lo", "median", "q95", "source", "n")}
        res["details"][aln] = {
            "n_floor_ident_ok": df.height, "nu": {k: v / 2 for k, v in eta.items()},
            "eta_aa_ml_measured": float(np.nanmedian(ra["ml"])), "eta_aa_ml_mean": float(np.nanmean(ra["ml"])),
            "eta_tdi_ml_alpha1_recomputed": float(np.nanmedian(r1["ml"])),
            "eta_tdi_ml_alpha1_table": float(df["tdi_ml"].median()),
            "eta_lddtD_unconditional_S0": float(np.nanmedian(Du)),
            "frac_lddtD_nan": float(np.isnan(D).mean()), "eta_by_len_bin": by_bin,
            "s0_unconditional": {k: s0u[k] for k in ("bin_lo", "median", "q95", "n")}}
        print(aln, eta, "aa", res["details"][aln]["eta_aa_ml_measured"], "tdi@1", res["details"][aln]["eta_tdi_ml_alpha1_recomputed"],
              "vs table", res["details"][aln]["eta_tdi_ml_alpha1_table"], flush=True)
    res["provenance"] = provenance()
    print(write_json(res, out_json(a.root, "floors_tuned")))


if __name__ == "__main__":
    main()
