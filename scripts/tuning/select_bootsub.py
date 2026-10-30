#!/usr/bin/env python3
"""ANALYSIS_PLAN §2.3 step 1: choose boot_sub = up to N (20k) tune primary pairs, stratified over 10 aa_ml bins.

Population: split == 'tune', pair_set == 'primary', status == 'ok' in the primary alignment (tdiaa, D22).
aa_ml is recomputed from the tdiaa core counts at the tuned alpha_aa (results/tuning/alpha.json or tuning.yaml).
Allocation: equal per bin (N/10), water-filled (bins with fewer pairs give their share to the others);
sampling without replacement with strucevo.stats.common.stage_rng('se_calibration').
Writes ROOT/data/interim/pairs_bootsub/bootsub_pairs.tsv (run_pairs.py input; same metadata columns) and
ROOT/data/interim/pairs_bootsub/bootsub_selection.json. Then run (separate --out, never data/interim/pairs):
  python scripts/pairs/run_pairs.py --pairs ROOT/data/interim/pairs_bootsub/bootsub_pairs.tsv --store <store> \
      --out ROOT/data/interim/pairs_bootsub --dump-columns ROOT/data/interim/pairs_bootsub/columns --n-shards 8
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from se_common import add_root_arg, load_pairs, load_tune_counts, np, pl, read_alpha, split_md5, write_json

from strucevo.pairs.ml import ml_from_counts
from strucevo.stats.common import stage_rng

AA_BIN_EDGES = (0.0, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 20.0001)
META = ["gene_a", "gene_b", "pair_set", "hog", "split", "species_a", "species_b", "relation", "one_to_one", "dcg_id",
        "event_level", "dup_branch_from"]


def allocate(avail, total):
    k = len(avail)
    alloc = np.zeros(k, int)
    left = total
    open_ = np.array(avail) > 0
    while left > 0 and open_.any():
        share = max(left // open_.sum(), 1)
        for i in np.where(open_)[0]:
            add = min(share, avail[i] - alloc[i], left)
            alloc[i] += add
            left -= add
            if alloc[i] >= avail[i]:
                open_[i] = False
            if left == 0:
                break
    return alloc


def main():
    ap = add_root_arg(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter))
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--aln", default="tdiaa")
    ap.add_argument("--alpha-aa", type=float)
    a = ap.parse_args()
    alpha, src = ({"aa": a.alpha_aa}, "cli") if a.alpha_aa is not None else read_alpha(a.root)
    df = load_pairs(a.aln, "tune", root=a.root, pair_set="primary").filter(pl.col("status") == "ok").sort("pair_id")
    C = load_tune_counts(a.root, a.aln, df["pair_id"].to_list(), keys=("aa",))
    aa = ml_from_counts(C["aa"], "aa", alpha["aa"])["ml"]
    b = np.clip(np.searchsorted(AA_BIN_EDGES, aa, side="right") - 1, 0, len(AA_BIN_EDGES) - 2)
    avail = np.bincount(b, minlength=len(AA_BIN_EDGES) - 1)
    alloc = allocate(avail, min(a.n, len(b)))
    rng = stage_rng("se_calibration")
    pick = np.concatenate([rng.choice(np.where(b == i)[0], alloc[i], replace=False) for i in range(len(alloc))])
    pick.sort()
    sel = df[pick].select(META).with_columns(pl.col("one_to_one").cast(pl.Utf8))
    out = Path(a.root) / "data" / "interim" / "pairs_bootsub"
    out.mkdir(parents=True, exist_ok=True)
    sel.write_csv(out / "bootsub_pairs.tsv", separator="\t")
    info = {"n_population": df.height, "n_selected": sel.height, "aln": a.aln, "alpha_aa": alpha["aa"],
            "alpha_source": src, "bin_edges": AA_BIN_EDGES, "available": avail, "allocated": alloc,
            "seed_stage": "se_calibration", "split_md5": split_md5(a.root)}
    write_json(info, out / "bootsub_selection.json")
    print(json.dumps({k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in info.items()}))


if __name__ == "__main__":
    main()
