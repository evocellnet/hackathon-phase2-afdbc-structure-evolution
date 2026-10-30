#!/usr/bin/env python3
"""Print elapsed time since START.txt (project clock). Check at every milestone."""
import re, time, pathlib
root = pathlib.Path(__file__).resolve().parents[1]
txt = (root / "START.txt").read_text()
starts = [int(x) for x in re.findall(r"START_EPOCH: (\d+)", txt)]
now = time.time()
el = (now - starts[0]) / 3600
print(f"Total elapsed since project start: {el:.2f} h (budget 48 h of compute + node switches)")
if len(starts) > 1:
    print(f"Elapsed on current node: {(now - starts[-1]) / 3600:.2f} h")
else:
    print(f"Node-1 remaining before vacate (11.5 h safe limit): {11.5 - el:.2f} h")
