"""§Tuning rows 'HOGs per stratum; total primary + E6 pairs' and 'pairs/s' (ANALYSIS_PLAN §3.2).

Sample: ROOT/data/interim/hogs_build_summary.json (final.per_stratum_sampled, primary_pairs, e6_pairs) plus the
tune/confirm split per stratum from ROOT/data/interim/split.tsv.
Throughput: every ROOT/data/interim/pairs/{aln}/*.done sentinel (rows_out, timings.total_s). Shards run one after
another, each with all threads, so pairs/s = sum(rows_out) / sum(total_s) at the thread count of the run
(--threads-used; the .done files do not record it; run_pairs.py default = 48). Scaled to 48 threads linearly
(an optimistic bound when the run used fewer threads). Writes ROOT/results/tuning/sample_throughput.json.
"""
import argparse
import json
from pathlib import Path

import numpy as np

import tuning_common as tc

ap = argparse.ArgumentParser()
ap.add_argument("--root", default=str(tc.PROJECT))
ap.add_argument("--threads-used", type=int, default=48, help="--threads of the run_pairs.py run that wrote the .done files")
ap.add_argument("--alns", nargs="+", default=["tmalign", "tdiaa", "seqonly"])
a = ap.parse_args()
root = Path(a.root)
rec = tc.header(root, "sample_throughput", a)
summ = json.loads((root / "data" / "interim" / "hogs_build_summary.json").read_text())
fin = summ.get("final", {})
import pandas as pd
sp = pd.read_csv(root / "data" / "interim" / "split.tsv", sep="\t")
by = sp.groupby(["stratum", "split"]).size().unstack(fill_value=0)
rec["sample"] = {
    "hogs_per_stratum": {str(k): int(v) for k, v in fin.get("per_stratum_sampled", {}).items()},
    "hogs_per_stratum_eligible": {str(k): int(v) for k, v in fin.get("per_stratum_eligible", {}).items()},
    "split_per_stratum": {str(int(s)): {c: int(by.loc[s, c]) for c in by.columns} for s in by.index},
    "n_hogs_sampled": fin.get("n_sampled"), "n_primary_pairs": fin.get("primary_pairs"),
    "n_e6_pairs": fin.get("e6_pairs"), "n_total_pairs": fin.get("total_pairs"), "e6_triplets": fin.get("e6_triplets"),
    "split_md5_in_summary": summ.get("split_md5"), "split_md5_now": tc.split_md5(root),
}
thr = {}
for aln in a.alns:
    files = sorted((root / "data" / "interim" / "pairs" / aln).glob("*.done"))
    rows, secs, n = 0, 0.0, 0
    for f in files:
        try:
            d = json.loads(f.read_text())
        except Exception:
            continue
        if d.get("skipped"):
            continue
        ts = d.get("timings", {}).get("total_s")
        if ts is None:
            continue
        rows += int(d.get("rows_out", 0)); secs += float(ts); n += 1
    if n == 0:
        thr[aln] = {"n_shards": 0, "pairs_per_s": None}
        continue
    pps = rows / secs if secs > 0 else None
    thr[aln] = {"n_shards": n, "rows": rows, "seconds": secs, "threads_used": a.threads_used,
                "pairs_per_s": pps, "pairs_per_s_48threads": (pps * 48 / a.threads_used) if pps else None,
                "scaling": "linear in threads (optimistic)" if a.threads_used != 48 else "measured at 48"}
rec["throughput"] = thr
npairs = rec["sample"]["n_total_pairs"]
rec["projection_full_run_h_48threads"] = {aln: (npairs / v["pairs_per_s_48threads"] / 3600 if v.get("pairs_per_s_48threads") and npairs else None)
                                          for aln, v in thr.items()}
tc.dump(rec, tc.out_path(root, "sample_throughput"))
print(json.dumps({"sample": rec["sample"], "throughput": thr}, indent=1, default=str))
