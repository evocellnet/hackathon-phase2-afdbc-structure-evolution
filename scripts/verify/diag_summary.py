"""Summarise a bias-diagnostic jsonl: mean (SE) of each estimate per config."""
import json
import sys
from collections import defaultdict

import numpy as np

rows = defaultdict(list)
for line in open(sys.argv[1]):
    r = json.loads(line)
    rows[r["config"]].append(r)
keys = sys.argv[2].split(",") if len(sys.argv) > 2 else ["c_naive", "c_simex", "c_rational", "c_cubic"]
for name, rr in rows.items():
    s = [f"{name:28s} n={len(rr):3d}"]
    for k in keys:
        v = np.array([r.get(k, np.nan) for r in rr], float)
        v = v[np.isfinite(v)]
        if v.size:
            s.append(f"{k}={v.mean():+.4f}({v.std(ddof=1) / np.sqrt(v.size) if v.size > 1 else 0:.4f}|sd {v.std(ddof=1) if v.size > 1 else 0:.3f})")
    print("  ".join(s))
