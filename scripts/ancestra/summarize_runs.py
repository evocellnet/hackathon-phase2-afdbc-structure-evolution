#!/usr/bin/env python3
"""Runtime/memory table for Ancestra runs made with run_ancestra.py (reads each run's
ancestra/manifest.json and the /usr/bin/time block of ancestra.log).
Usage: python scripts/ancestra/summarize_runs.py results/ancestra/runs/* > results/ancestra/runtime_table.tsv"""
import json, os, re, sys
cols = ["run", "preset", "threads", "n_species", "n_genes", "n_families", "max_family", "sum_g2",
        "wall_s", "cpu_s", "peak_rss_gb", "families_s", "scores_s", "merge_s", "merge_tree_s", "output_s",
        "n_root_hogs", "n_ortholog_pairs"]
print("\t".join(cols))
for d in sys.argv[1:]:
    man = os.path.join(d, "ancestra", "manifest.json")
    if not os.path.exists(man):
        continue
    m = json.load(open(man))
    log = open(os.path.join(d, "ancestra.log")).read() if os.path.exists(os.path.join(d, "ancestra.log")) else ""
    g = lambda pat: re.search(pat, log)
    el = g(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\): ([\d:.]+)")
    wall = None
    if el:
        p = [float(x) for x in el.group(1).split(":")]
        wall = sum(v * 60 ** i for i, v in enumerate(reversed(p)))
    cpu = sum(float(x) for x in re.findall(r"(?:User|System) time \(seconds\): ([\d.]+)", log)) or None
    rss = g(r"Maximum resident set size \(kbytes\): (\d+)")
    st = m["stages"]
    fs = st["families"]["stats"]
    row = [os.path.basename(d.rstrip("/")), m["config"]["preset"], m["config"]["threads"], m["n_species"], m["n_genes"],
           fs.get("n_families"), fs.get("max_family"), fs.get("sum_g2"), round(wall or m["total_wall_s"]),
           round(cpu) if cpu else "", round(int(rss.group(1)) / 1e6, 2) if rss else "",
           st["families"]["wall_s"], st["scores"]["wall_s"], st["merge"]["wall_s"],
           st["merge"]["stats"].get("tree_s"), st["output"]["wall_s"],
           st["merge"]["stats"].get("n_root_groups"), st["output"]["stats"].get("n_pairs")]
    print("\t".join(str(x) for x in row))
