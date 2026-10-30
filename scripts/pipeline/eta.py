#!/usr/bin/env python3
"""ETA / progress for the post-Ancestra pipeline (scripts/run_after_ancestra.sh).

  python scripts/pipeline/eta.py predict --pairs N --alns seqonly tdiaa --threads 48 [--measured OUT --offset-range 0 999]
  python scripts/pipeline/eta.py progress --out data/interim/pairs --alns seqonly tdiaa tmalign --offset-range 1000 1999 --expected N
  python scripts/pipeline/eta.py rates --out data/interim/pairs [--offset-range a b]      # measured pairs/s per aln

Rates: measured from the .done sentinels (rows_in / timings.total_s of shards that were actually computed) when
available; otherwise the demo benchmark at 48 threads (agents/pairs_pipeline_report.md: seqonly 1,400, tdiaa 1,140,
tmalign 322 pairs/s), scaled linearly by threads/48 (an optimistic bound under node load).
"""
from __future__ import annotations

import argparse
import glob
import json
import os

BENCH48 = {"seqonly": 1400.0, "tdiaa": 1140.0, "tmalign": 322.0}


def done_infos(out, aln, lo=0, hi=10 ** 9):
    res = []
    for f in glob.glob(os.path.join(out, aln, "shard_*.done")):
        k = int(os.path.basename(f)[6:-5])
        if lo <= k <= hi:
            try:
                res.append(json.load(open(f)))
            except Exception:
                pass
    return res


def measured_rate(out, aln, lo=0, hi=10 ** 9):
    rows = secs = 0.0
    for i in done_infos(out, aln, lo, hi):
        if i.get("skipped"):
            continue
        t = i.get("timings", {}).get("total_s")
        if t:
            rows += i["rows_in"]; secs += t
    return (rows / secs, rows, secs) if secs > 0 and rows >= 500 else (None, rows, secs)


def rate_for(aln, threads, measured=None, lo=0, hi=10 ** 9):
    if measured:
        r, _, _ = measured_rate(measured, aln, lo, hi)
        if r:
            return r, "measured"
    return BENCH48[aln] * threads / 48.0, f"benchmark@48 scaled to {threads} thr"


def fmt(s):
    s = int(s)
    return f"{s // 3600}h{(s % 3600) // 60:02d}m" if s >= 3600 else f"{s // 60}m{s % 60:02d}s"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("predict")
    p.add_argument("--pairs", type=int, required=True)
    p.add_argument("--alns", nargs="+", required=True)
    p.add_argument("--threads", type=int, default=48)
    p.add_argument("--measured")
    p.add_argument("--offset-range", type=int, nargs=2, default=[0, 10 ** 9])
    p.add_argument("--label", default="")
    g = sub.add_parser("progress")
    g.add_argument("--out", required=True)
    g.add_argument("--alns", nargs="+", required=True)
    g.add_argument("--offset-range", type=int, nargs=2, default=[0, 10 ** 9])
    g.add_argument("--expected", type=int, required=True)
    r = sub.add_parser("rates")
    r.add_argument("--out", required=True)
    r.add_argument("--offset-range", type=int, nargs=2, default=[0, 10 ** 9])
    a = ap.parse_args()
    if a.cmd == "predict":
        tot = 0.0
        parts = []
        for aln in a.alns:
            rt, src = rate_for(aln, a.threads, a.measured, *a.offset_range)
            t = a.pairs / rt
            tot += t
            parts.append(f"{aln} {rt:,.0f} pairs/s ({src}) -> {fmt(t)}")
        print(f"[ETA] {a.label} {a.pairs:,} pairs: " + "; ".join(parts) + f"; total ~{fmt(tot)}", flush=True)
    elif a.cmd == "progress":
        for aln in a.alns:
            inf = done_infos(a.out, aln, *a.offset_range)
            rows = sum(i["rows_in"] for i in inf)
            rt, _, _ = measured_rate(a.out, aln, *a.offset_range)
            rem = max(0, a.expected - rows)
            eta = fmt(rem / rt) if rt else "n/a"
            print(f"[progress] {aln}: {rows:,}/{a.expected:,} pairs in {len(inf)} shards; "
                  f"rate {rt and round(rt)} pairs/s; remaining ~{eta}", flush=True)
    else:
        for aln in BENCH48:
            rt, rows, secs = measured_rate(a.out, aln, *a.offset_range)
            print(json.dumps({"aln": aln, "pairs_per_s": rt, "rows": rows, "seconds": round(secs, 1)}))


if __name__ == "__main__":
    main()
