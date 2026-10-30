#!/usr/bin/env python3
"""ANALYSIS_PLAN §3.2 / D22: tune the discrete-Gamma shape alpha for LG (aa) and Q.3Di.AF (tdi).

Pooled profile likelihood (per-pair d free; strucevo.distances.profile_alpha, validated in
profile_alpha_validation.json) over the core count matrices of split == 'tune', pair_set == 'primary',
status == 'ok' pairs with aa_ml_seq <= 2 (aa_ml_seq at fixed alpha = 1, C26). Coarse grid then Brent
refinement on log(alpha) (precision << 0.05) and a profile 95% CI. Primary alignment = tdiaa (D22), but
alpha is taken from seqonly (DEVIATIONS C29); tdiaa / tmalign values are reported for information only.

  python scripts/tuning/tune_alpha.py [--root R] [--alns tdiaa seqonly tmalign]
Output: ROOT/results/tuning/alpha.json
Fragment keys -> config/tuning.yaml (top level, same nesting): alpha.aa, alpha.tdi (rounded to 0.01), split_md5.
alpha_selection_aa stays 1.0 (C26). Details per alignment under by_alignment. If the likelihood is maximal at the
grid end (alpha_at_grid_boundary), alpha = 1000 means "no rate heterogeneity" and must be flagged to the PI.
"""
from __future__ import annotations

import argparse
import time

from se_common import (ALNS, add_root_arg, fit_alpha, load_pairs, load_tune_counts, models, np, out_json,
                           pl, provenance, split_md5, write_json)


def main():
    ap = add_root_arg(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter))
    ap.add_argument("--alns", nargs="+", default=list(ALNS))
    ap.add_argument("--max-aa-ml-seq", type=float, default=2.0)
    a = ap.parse_args()
    t0 = time.time()
    res = {"step": "tune_alpha", "root": a.root, "split": "tune", "split_md5": split_md5(a.root),
           "selection": {"pair_set": "primary", "status": "ok", "aa_ml_seq_max": a.max_aa_ml_seq,
                         "aa_ml_seq_alpha": 1.0},
           "method": "pooled profile likelihood over pairs (per-pair d profiled, d_max=20, ncat=4), coarse grid "
                     "+ bounded Brent on log(alpha); CI = profile likelihood drop 1.92",
           "by_alignment": {}}
    M = models()
    for aln in a.alns:
        df = load_pairs(aln, "tune", columns=["pair_id", "split", "pair_set", "status", "aa_ml_seq", "n_core"],
                        root=a.root, pair_set="primary")
        n_all = df.height
        df = df.filter((pl.col("status") == "ok") & (pl.col("aa_ml_seq") <= a.max_aa_ml_seq))
        C = load_tune_counts(a.root, aln, df["pair_id"].to_list())
        out = {"n_primary_tune": n_all, "n_used": df.height}
        for key in ("aa", "tdi"):
            r = fit_alpha(C[key], M[key])
            out[key] = r
            print(f"{aln} {key}: alpha={r['alpha']:.3f} CI95={r['ci95'][0]:.3f}-{r['ci95'][1]:.3f} "
                  f"n={df.height} ({r['seconds']} s)", flush=True)
        res["by_alignment"][aln] = out
    # DEVIATIONS C29: alpha (AA and 3Di) is taken from the seqonly alignment (3Di columns not optimised), applied to all
    res["alpha_source_alignment"] = "seqonly"
    p = res["by_alignment"].get("seqonly")
    # top-level keys = config/tuning.yaml nesting (collector merges them); 2 decimals (precision target +-0.05)
    res["alpha"] = {"aa": round(p["aa"]["alpha"], 2) if p else None, "tdi": round(p["tdi"]["alpha"], 2) if p else None}
    res["alpha_at_grid_boundary"] = {k: bool(p[k]["at_grid_boundary"]) for k in ("aa", "tdi")} if p else None
    res["provenance"] = provenance()
    res["seconds"] = round(time.time() - t0, 1)
    print(write_json(res, out_json(a.root, "alpha")))


if __name__ == "__main__":
    main()
