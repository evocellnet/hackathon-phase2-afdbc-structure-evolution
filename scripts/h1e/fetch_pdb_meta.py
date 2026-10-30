#!/usr/bin/env python
"""Bulk per-entry metadata (method, resolution, initial release date) for all PDB entries
in SIFTS pdb_chain_uniprot, via RCSB Data API GraphQL (the wwPDB entries.idx only has the
*deposition* date, not the release date). Output: data/raw/pdb_meta/pdb_entry_meta.tsv.gz"""
import gzip, sys, time, concurrent.futures as cf, requests
ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution/data/raw/pdb_meta"
ids = set()
with gzip.open(f"{ROOT}/pdb_chain_uniprot.tsv.gz", "rt") as f:
    next(f); next(f)
    for l in f: ids.add(l.split("\t", 1)[0].upper())
ids = sorted(ids); B = 400
Q = ('query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id exptl{method} '
     'rcsb_entry_info{resolution_combined experimental_method} '
     'rcsb_accession_info{initial_release_date deposit_date}}}')
def run(chunk):
    for t in range(6):
        try:
            r = requests.post("https://data.rcsb.org/graphql", json={"query": Q, "variables": {"ids": chunk}}, timeout=120)
            r.raise_for_status(); d = r.json()
            if d.get("data") and d["data"].get("entries") is not None: return d["data"]["entries"]
        except Exception as e: err = e
        time.sleep(2 * (t + 1))
    raise RuntimeError(f"failed chunk {chunk[:2]} {err}")
rows = []
with cf.ThreadPoolExecutor(6) as ex:
    for i, ents in enumerate(ex.map(run, [ids[i:i+B] for i in range(0, len(ids), B)])):
        for e in ents:
            if e is None: continue
            res = (e["rcsb_entry_info"] or {}).get("resolution_combined") or []
            ai = e["rcsb_accession_info"] or {}
            rows.append((e["rcsb_id"].lower(), "|".join(x["method"] for x in (e["exptl"] or [])),
                         (e["rcsb_entry_info"] or {}).get("experimental_method") or "",
                         f"{min(res):.3f}" if res else "", (ai.get("initial_release_date") or "")[:10],
                         (ai.get("deposit_date") or "")[:10]))
        if i % 50 == 0: print(i, len(rows), file=sys.stderr, flush=True)
with gzip.open(f"{ROOT}/pdb_entry_meta.tsv.gz", "wt") as o:
    o.write("pdb_id\texptl_method\trcsb_method\tresolution\tinitial_release_date\tdeposit_date\n")
    for r in rows: o.write("\t".join(r) + "\n")
print("entries requested", len(ids), "returned", len(rows))
