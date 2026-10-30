#!/usr/bin/env python3
"""Project clock. START.txt holds one START_EPOCH line per node (appended at each node start).
Prints total elapsed and per-node elapsed/remaining. Check at every milestone."""
import re, time, pathlib
root = pathlib.Path(__file__).resolve().parents[1]
txt = (root / "START.txt").read_text()
starts = [int(x) for x in re.findall(r"START_EPOCH: (\d+)", txt)]
now = time.time()
print(f"Wall-clock since project start: {(now - starts[0]) / 3600:.2f} h")
node_budget = [11.5] + [47.5] * (len(starts) - 1)   # safe limits per node (h)
used = 0.0
for i, s in enumerate(starts):
    end = starts[i + 1] if i + 1 < len(starts) else now
    used += (end - s) / 3600
cur = (now - starts[-1]) / 3600
print(f"Node {len(starts)}: elapsed {cur:.2f} h, remaining before safe limit {node_budget[len(starts)-1] - cur:.2f} h")
print(f"Compute-time used across nodes: {used:.2f} h of 48 h project budget (+12 h node-1 bonus)")
