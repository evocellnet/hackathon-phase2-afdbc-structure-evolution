#!/usr/bin/env python3
"""Query TimeTree 6 pairwise divergence times for all species pairs of the eukaryotic panel.

Output (NOT committed; TimeTree licence forbids redistribution):
  data/interim/timetree/pairwise_ages.tsv : leaf_a, leaf_b, taxid_a, taxid_b (species-level NCBI), tt_name_a,
      tt_name_b, n_studies, median_age (TimeTree 'precomputed' median of studies), ci_low, ci_high,
      adjusted_age (TimeTree estimated/synthesis time), age_mya (= adjusted_age if > 0 else median_age),
      ncbi_mrca, ncbi_mrca_name, ncbi_mrca_rank
Summary (aggregate counts only, committed): results/panel/timetree_age_summary.json
"""
from __future__ import annotations

import itertools
import json
import math
import os
import sys

import numpy as np
import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, Taxonomy, TimeTree  # noqa: E402


def main():
    panel = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data/processed/panel/species_panel.tsv")
    p = pl.read_csv(panel, separator="\t")
    tx = Taxonomy()
    tt = TimeTree()
    sp = dict(zip(p["leaf"], p["species_taxid"]))
    pairs = list(itertools.combinations(sorted(sp), 2))
    res = tt.pairwise_many([(sp[a], sp[b]) for a, b in pairs])
    rows = []
    for a, b in pairs:
        ta, tb = sp[a], sp[b]
        r = res.get((min(ta, tb), max(ta, tb))) or {}
        adj, med = r.get("adjusted_age"), r.get("precomputed_age")
        age = adj if adj and adj > 0 else (med if med and med > 0 else None)
        m = tx.mrca(ta, tb)
        nm_a, nm_b = (r.get("name_a"), r.get("name_b")) if ta <= tb else (r.get("name_b"), r.get("name_a"))
        rows.append({"leaf_a": a, "leaf_b": b, "taxid_a": ta, "taxid_b": tb, "tt_name_a": nm_a, "tt_name_b": nm_b,
                     "n_studies": r.get("all_total"), "median_age": med, "ci_low": r.get("ci_low"),
                     "ci_high": r.get("ci_high"), "adjusted_age": adj, "age_mya": age,
                     "age_source": "adjusted" if adj and adj > 0 else ("median" if age else None),
                     "ncbi_mrca": m, "ncbi_mrca_name": tx.name.get(m), "ncbi_mrca_rank": tx.rank.get(m)})
    df = pl.DataFrame(rows)
    out = os.path.join(ROOT, "data/interim/timetree/pairwise_ages.tsv")
    df.write_csv(out, separator="\t")
    ok = df.filter(pl.col("age_mya").is_not_null())
    ages = ok["age_mya"].to_numpy()
    la = np.log10(ages)
    edges = np.arange(0, 3.31 + 1e-9, 0.3)
    h_pairs = np.histogram(np.clip(la, 0, 3.3), bins=edges)[0]
    # distinct divergence events: distinct (NCBI MRCA, rounded age) combinations and distinct ages
    nodes = ok.group_by("ncbi_mrca", "ncbi_mrca_name").agg(pl.col("age_mya").median().alias("age"), pl.len().alias("n_pairs"))
    distinct_ages = np.unique(np.round(ages, 3))
    h_nodes = np.histogram(np.clip(np.log10(distinct_ages), 0, 3.3), bins=edges)[0]
    # nearest-neighbour age of each species
    nn = (pl.concat([ok.select(pl.col("leaf_a").alias("leaf"), "age_mya"), ok.select(pl.col("leaf_b").alias("leaf"), "age_mya")])
          .group_by("leaf").agg(pl.col("age_mya").min().alias("nn_age")))
    h_nn = np.histogram(np.clip(np.log10(nn["nn_age"].to_numpy()), 0, 3.3), bins=edges)[0]

    def ks_loguniform(x):
        x = np.sort(x)
        lo, hi = x.min(), x.max()
        u = (x - lo) / (hi - lo)
        n = len(x)
        return float(max(np.max(np.arange(1, n + 1) / n - u), np.max(u - np.arange(n) / n)))
    lab = [f"{10**edges[i]:.0f}-{10**edges[i+1]:.0f}" for i in range(len(edges) - 1)]
    summ = {
        "panel": os.path.relpath(panel, ROOT), "n_species": len(sp), "n_pairs": len(pairs),
        "n_pairs_with_age": ok.height, "n_pairs_no_answer": df.height - ok.height,
        "pairs_no_answer": df.filter(pl.col("age_mya").is_null()).select("leaf_a", "leaf_b", "ncbi_mrca_name").rows()[:200],
        "n_pairs_age_from_median_only": int((ok["age_source"] == "median").sum()),
        "age_min": float(ages.min()), "age_max": float(ages.max()), "age_median": float(np.median(ages)),
        "log_bins_mya": lab, "hist_pairs": h_pairs.tolist(), "hist_distinct_ages": h_nodes.tolist(),
        "hist_nearest_neighbour_age": h_nn.tolist(),
        "n_distinct_ages": int(len(distinct_ages)), "n_distinct_ncbi_mrca_nodes": nodes.height,
        "ks_distance_loguniform_distinct_ages": ks_loguniform(np.log10(distinct_ages)),
        "ks_distance_loguniform_pairs": ks_loguniform(la),
        "n_pairs_lt_10mya": int((ages < 10).sum()), "n_pairs_10_100mya": int(((ages >= 10) & (ages < 100)).sum()),
        "n_pairs_100_1000mya": int(((ages >= 100) & (ages < 1000)).sum()), "n_pairs_ge_1000mya": int((ages >= 1000).sum()),
        "timetree_api_version": "20260820", "note": "aggregate counts only; per-pair TimeTree ages are not redistributed",
    }
    os.makedirs(os.path.join(ROOT, "results/panel"), exist_ok=True)
    with open(os.path.join(ROOT, "results/panel/timetree_age_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k != "pairs_no_answer"}, indent=1))


if __name__ == "__main__":
    main()
